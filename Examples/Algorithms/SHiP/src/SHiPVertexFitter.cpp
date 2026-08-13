#include "ActsExamples/SHiP/SHiPVertexFitter.hpp"

#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/NumericalTrackLinearizer.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <iostream>

namespace ActsExamples {

SHiPVertexFitter::SHiPVertexFitter(const Config& cfg, std::shared_ptr<const Acts::MagneticFieldProvider> bField)
  : m_cfg(cfg), m_bField(std::move(bField)) {}

std::vector<Acts::Vertex> SHiPVertexFitter::fit(
    const std::vector<ActsExamples::ConstTrackContainer::ConstTrackProxy>& proxies,
    const Acts::GeometryContext& geoCtx,
    const Acts::TrackingGeometry& trackingGeometry) const
{
    std::vector<Acts::Vertex> out;
    if (proxies.size() < 2) return out;

    Acts::MagneticFieldContext magCtx;

    // 1) compute PCA seed (analytical two-line closest point for 2 tracks, general PCA otherwise)
    Acts::Vector3 seedPos = Acts::Vector3::Zero();
    if (proxies.size() == 2) {
        std::vector<Acts::BoundTrackParameters> bparams;
        bparams.reserve(2);
        for (size_t i = 0; i < 2; ++i) {
            const auto& proxy = proxies[i];
            Acts::BoundVector paramsVec = proxy.parameters();
            Acts::BoundMatrix covMat = proxy.covariance();
            auto surfacePtr = proxy.hasReferenceSurface() ? proxy.referenceSurface().getSharedPtr() : std::shared_ptr<const Acts::Surface>();
            if (!surfacePtr) {
                surfacePtr = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Transform3(Acts::Translation3((m_cfg.seedXMin + m_cfg.seedXMax) / 2.0, 0.0, 0.0)));
            }
            bparams.emplace_back(surfacePtr, paramsVec, std::optional<Acts::BoundMatrix>(covMat), Acts::ParticleHypothesis::pion());
        }
        Acts::Vector3 p1 = bparams[0].position(geoCtx);
        Acts::Vector3 u1 = bparams[0].momentum().normalized();
        Acts::Vector3 p2 = bparams[1].position(geoCtx);
        Acts::Vector3 u2 = bparams[1].momentum().normalized();

        Acts::Vector3 w0 = p1 - p2;
        double a = u1.dot(u1);
        double b = u1.dot(u2);
        double c = u2.dot(u2);
        double d = u1.dot(w0);
        double e = u2.dot(w0);

        double denom = a * c - b * b;
        double t = 0.0, s = 0.0;
        if (std::abs(denom) < 1e-12) {
            t = (u1.dot(p2 - p1)) / a;
            s = 0.0;
        } else {
            t = (b * e - c * d) / denom;
            s = (a * e - b * d) / denom;
        }
        Acts::Vector3 c1 = p1 + u1 * t;
        Acts::Vector3 c2 = p2 + u2 * s;
        seedPos = Acts::Vector3(0.5 * (c1.x() + c2.x()), 0.5 * (c1.y() + c2.y()), 0.5 * (c1.z() + c2.z()));
        // clamp
        seedPos.x() = std::clamp(seedPos.x(), m_cfg.seedXMin, m_cfg.seedXMax);
    } else {
        // general PCA seeder
        Acts::SquareMatrix3 A = Acts::SquareMatrix3::Zero();
        Acts::Vector3 bvec = Acts::Vector3::Zero();
        for (const auto& proxy : proxies) {
            Acts::BoundVector paramsVec = proxy.parameters();
            auto surfacePtr = proxy.hasReferenceSurface() ? proxy.referenceSurface().getSharedPtr() : std::shared_ptr<const Acts::Surface>();
            Acts::BoundTrackParameters btp(surfacePtr, paramsVec, std::optional<Acts::BoundMatrix>(proxy.covariance()), Acts::ParticleHypothesis::pion());
            Acts::Vector3 p = btp.position(geoCtx);
            Acts::Vector3 n = btp.momentum().normalized();
            Acts::SquareMatrix3 projection = Acts::SquareMatrix3::Identity() - (n * n.transpose());
            A += projection;
            bvec += projection * p;
        }
        seedPos = A.colPivHouseholderQr().solve(bvec);
        if (!std::isfinite(seedPos.x())) seedPos = Acts::Vector3((m_cfg.seedXMin + m_cfg.seedXMax) / 2.0, 0.0, 0.0);
        seedPos.x() = std::clamp(seedPos.x(), m_cfg.seedXMin, m_cfg.seedXMax);
    }

    // Build perigee surface at seed
    auto perigeeSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Transform3(Acts::Translation3(seedPos)));

    // Build extracted perigee parameters per proxy
    std::vector<Acts::BoundTrackParameters> extractedParams;
    extractedParams.reserve(proxies.size());
    for (const auto& proxy : proxies) {
        Acts::Vector3 globalPos = proxy.hasReferenceSurface() ? proxy.referenceSurface().center(geoCtx) : seedPos;
        Acts::Vector3 globalMom = proxy.momentum();
        double charge = proxy.charge();
        double delta_x = seedPos.x() - globalPos.x();
        Acts::Vector3 dir = globalMom.normalized();
        Acts::Vector3 extrapolatedPos = globalPos + dir * (delta_x / (dir.x() + 1e-12));
        auto localRes = perigeeSurface->globalToLocal(geoCtx, extrapolatedPos, dir);
        if (!localRes.ok()) continue;
        Acts::BoundVector perigeeParams = Acts::BoundVector::Zero();
        perigeeParams[Acts::eBoundLoc0] = localRes.value()[Acts::eBoundLoc0];
        perigeeParams[Acts::eBoundLoc1] = localRes.value()[Acts::eBoundLoc1];
        perigeeParams[Acts::eBoundPhi] = std::atan2(globalMom.y(), globalMom.x());
        perigeeParams[Acts::eBoundTheta] = std::acos(globalMom.z() / std::max(1e-12, globalMom.norm()));
        perigeeParams[Acts::eBoundQOverP] = charge / std::max(1e-12, globalMom.norm());
        perigeeParams[Acts::eBoundTime] = 0.0;
        Acts::BoundMatrix perigeeCov = Acts::BoundMatrix::Identity();
        Acts::BoundMatrix originalCov = proxy.covariance();
        perigeeCov(Acts::eBoundLoc0, Acts::eBoundLoc0) = originalCov(Acts::eBoundLoc0, Acts::eBoundLoc0) * m_cfg.perigeeInflation;
        perigeeCov(Acts::eBoundLoc1, Acts::eBoundLoc1) = originalCov(Acts::eBoundLoc1, Acts::eBoundLoc1) * m_cfg.perigeeInflation;
        perigeeCov(Acts::eBoundPhi,   Acts::eBoundPhi)   = originalCov(Acts::eBoundPhi,   Acts::eBoundPhi);
        perigeeCov(Acts::eBoundTheta, Acts::eBoundTheta) = originalCov(Acts::eBoundTheta, Acts::eBoundTheta);
        perigeeCov(Acts::eBoundQOverP,Acts::eBoundQOverP)= originalCov(Acts::eBoundQOverP,Acts::eBoundQOverP);
        perigeeCov(Acts::eBoundTime,  Acts::eBoundTime)  = m_cfg.timeInflation;
        extractedParams.emplace_back(perigeeSurface, perigeeParams, std::optional<Acts::BoundMatrix>(perigeeCov), Acts::ParticleHypothesis::pion());
    }

    if (extractedParams.empty()) return out;

    // Prepare vertex fitter
    Acts::Vertex seedVertex(seedPos);
    Acts::SquareMatrix4 seedCov = Acts::SquareMatrix4::Identity() * 1000000.0;
    seedVertex.setFullCovariance(seedCov);
    Acts::VertexingOptions vtxOptions(geoCtx, magCtx, seedVertex);

    using Stepper = Acts::EigenStepper<>;
    using Navigator = Acts::VoidNavigator;
    using Propagator = Acts::Propagator<Stepper, Navigator>;
    Stepper::Config stepperConfig;
    auto propagator = std::make_shared<Propagator>(Stepper(stepperConfig));

    Acts::NumericalTrackLinearizer::Config ltConfig(m_bField, propagator);
    Acts::NumericalTrackLinearizer linearizer(ltConfig, Acts::getDefaultLogger("NumLin", Acts::Logging::INFO));

    Acts::FullBilloirVertexFitter::Config fCfg;
    fCfg.maxIterations = m_cfg.maxIterations;
    // extractor callable: convert InputTrack back to BoundTrackParameters
    struct LocalTrackExtractor {
        Acts::BoundTrackParameters operator()(const Acts::InputTrack& it) const {
            return *(it.as<Acts::BoundTrackParameters>());
        }
    };
    static const LocalTrackExtractor lTrackExtractor;
    fCfg.extractParameters.connect(lTrackExtractor);
    fCfg.trackLinearizer.connect<&Acts::NumericalTrackLinearizer::linearizeTrack>(&linearizer);

    Acts::FullBilloirVertexFitter fitter(fCfg);

    // Build input tracks
    std::vector<Acts::InputTrack> inputTracks;
    inputTracks.reserve(extractedParams.size());
    for (size_t i = 0; i < extractedParams.size(); ++i) inputTracks.emplace_back(&extractedParams[i]);

    auto fieldCache = m_bField->makeCache(magCtx);
    auto result = fitter.fit(inputTracks, vtxOptions, fieldCache);
    if (result.ok()) {
        out.push_back(result.value());
    } else {
        std::cout << "SHiPVertexFitter: vertex fit failed: " << result.error().message() << std::endl;
    }
    return out;
}

} // namespace ActsExamples
