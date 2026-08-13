#include <algorithm>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>

#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <cstdint>
#include <unordered_map>

#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/VoidFitterComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Utilities/AnnealingUtility.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/NumericalTrackLinearizer.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

#include "ActsExamples/DetectorCommons/Detector.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/EventData/Vertex.hpp"
#include "ActsExamples/SHiP/HGCBuilder.hpp"
#include "ActsExamples/SHiP/HGCDetector.hpp"
#include "ActsExamples/SHiP/RecoTrack.hpp"
#include "ActsExamples/SHiP/RecoVertex.hpp"
#include "ActsExamples/SHiP/SHiPFieldProvider.hpp"
#include "ActsExamples/SHiP/SHiPMeasurementProvider.hpp"
#include "ActsExamples/SHiP/StrawtubeBuilder.hpp"
#include "ActsExamples/SHiP/StrawtubeDetector.hpp"
#include "ActsExamples/SHiP/DeterministicAnnealingFitter.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"


struct OwnedInput {
    int originalIndex = -1;
    Acts::BoundTrackParameters params;
    OwnedInput(int idx, const Acts::BoundTrackParameters& p)
        : originalIndex(idx), params(p) {}
};


static std::unordered_map<uintptr_t, std::vector<int>> g_vertexMatchedIndices;

static std::vector<Acts::BoundTrackParameters>* g_lastHeapParams = nullptr;

namespace ActsExamples {

namespace py = pybind11;


struct GlobalTrackExtractor {
    Acts::BoundTrackParameters operator()(const Acts::InputTrack& it) const {
        return *(it.as<Acts::BoundTrackParameters>());
    }
};
static const GlobalTrackExtractor gTrackExtractor;
// =========================================================================

} // namespace ActsExamples


namespace ActsPython {
    void addSHiP(pybind11::module& mex) {

    using namespace ActsExamples;
    auto m = pybind11::module_::import("acts");


    struct VertexTrackData {
        int trackIndex = -1;
        py::object position = py::none();
        py::object momentum = py::none();
        py::object covariance = py::none();
        std::string error = "";
    };

    py::class_<VertexTrackData>(m, "VertexTrack")
        .def(py::init<>())
        .def_readwrite("trackIndex", &VertexTrackData::trackIndex)
        .def_readwrite("position", &VertexTrackData::position)
        .def_readwrite("momentum", &VertexTrackData::momentum)
        .def_readwrite("covariance", &VertexTrackData::covariance)
        .def_readwrite("error", &VertexTrackData::error);


    m.def("processMeasurements", [](const std::vector<std::vector<float>>& hits, 
                                     std::shared_ptr<const Acts::TrackingGeometry> tg) {
        ActsExamples::SHiPMeasurementProvider::Config cfg;
        cfg.trackingGeometry = tg;
        ActsExamples::SHiPMeasurementProvider provider(cfg);
        return provider.process(hits, Acts::GeometryContext::dangerouslyDefaultConstruct());

    });

    py::class_<SHiPFieldProvider, Acts::MagneticFieldProvider, std::shared_ptr<SHiPFieldProvider>>(m, "SHiPFieldProvider")
        .def(py::init<const std::string&, double>());

    m.def("createShipFieldProvider", [](const std::string& filename, double scale) {
        auto provider = std::shared_ptr<SHiPFieldProvider>(new SHiPFieldProvider(filename, scale));
        return std::static_pointer_cast<Acts::MagneticFieldProvider>(provider);
    });

    m.def("getMeasurementGeoId", [](const ActsExamples::MeasurementContainer& measurements, size_t index) {
        if (index >= measurements.size()) {
            throw std::out_of_range("Measurement index out of range");
        }
        return measurements.getMeasurement(index).geometryId();
    });

    m.def("extrapolateTrack", [](const Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator>& prop,
                                  const Acts::BoundTrackParameters& start,
                                  const Acts::Surface& target,
                                  const Acts::GeometryContext& geoCtx,
                                  const Acts::MagneticFieldContext& magCtx) -> py::object {
        using Propagator = Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator>;
        using PropagatorOptions = Propagator::Options<>;
        PropagatorOptions options(geoCtx, magCtx);
        options.pathLimit = 100000.0; // 100m limit
        auto result = prop.propagate(start, target, options);
        if (result.ok()) {
            return py::cast(result.value().endParameters);
        } else {
            return py::none();
        }
    }, py::arg("propagator"), py::arg("start"), py::arg("target"), py::arg("geoCtx"), py::arg("magCtx"));

    // Change the first argument from const ActsExamples::RecoTrack& to uintptr_t
m.def("extrapolateTrackToZ", [](uintptr_t track_ptr_addr, double targetZ_ship) -> py::tuple {

        // Cast the raw memory address integer back into a proper C++ pointer
        const auto* track_ptr = reinterpret_cast<const ActsExamples::RecoTrack*>(track_ptr_addr);

        if (!track_ptr) {
            return py::make_tuple(false, py::none(), py::make_tuple(0.0, 0.0, 0.0));
        }

        // Dereference the pointer to safely use the original logic
        const ActsExamples::RecoTrack& track = *track_ptr;


        double px = track.px();
        double py = track.py();
        double pz = track.pz();

        if (std::abs(pz) < 1e-10) {
            return py::make_tuple(false, py::none(), py::make_tuple(0.0, 0.0, 0.0));
        }

        double track_x = track.x();
        double track_y = track.y();
        double track_z = track.z();

        double lambda = (targetZ_ship - track_z) / pz;

        double final_x = track_x + lambda * px;
        double final_y = track_y + lambda * py;
        double final_z = targetZ_ship;

        return py::make_tuple(
            true,
            py::make_tuple(final_x, final_y, final_z),
            py::make_tuple(px, py, pz)
        );
    }, py::arg("track_ptr_addr"), py::arg("targetZ_ship"));



    py::class_<ActsExamples::RecoTrack, std::shared_ptr<ActsExamples::RecoTrack>>(mex, "RecoTrack")
        .def(py::init<>());

    mex.def("pushRecoTrack", [](long vectorAddress, const Acts::GeometryContext& gctx, size_t track_idx, pybind11::object containerObj) {
        auto* container = reinterpret_cast<std::vector<ActsExamples::RecoTrack>*>(vectorAddress);
        if (!container) throw std::runtime_error("CRITICAL: Null pointer!");
        if (containerObj.is_none()) return;

        try {

            const auto& master_container = containerObj.cast<const ActsExamples::ConstTrackContainer&>();

            auto track = master_container.getTrack(track_idx);

            bool has_surface = (&track.referenceSurface() != nullptr);
            unsigned int nMeasurements = track.nMeasurements();

            if (has_surface && nMeasurements > 0) {
                auto params = track.parameters();
                Acts::BoundMatrix cov = track.covariance();
                const auto& surface = track.referenceSurface();

                float chi2        = track.chi2();
                unsigned int ndof = track.nDoF();

                std::vector<Double_t> residuals;
                std::vector<Double_t> pulls;
                residuals.reserve(nMeasurements);
                pulls.reserve(nMeasurements);

                bool states_parsed = false;

                try {
                    const auto& multitrajectory = master_container.trackStateContainer();
                    auto tipIndex = track.tipIndex();

                    unsigned int true_states_found = 0;

                    multitrajectory.visitBackwards(tipIndex, [&](const auto& state) {
                        if (!state.hasUncalibratedSourceLink() || !state.hasCalibrated() || !state.hasSmoothed()) {
                            return true; // Keep iterating backwards
                        }

                        double meas_loc0 = state.template calibrated<1>()(0);
                        double meas_err  = state.template calibratedCovariance<1>()(0, 0);

                        double smoothed_loc0 = state.smoothed()(Acts::eBoundLoc0);
                        double smoothed_err  = state.smoothedCovariance()(Acts::eBoundLoc0, Acts::eBoundLoc0);

                        double res_val = (meas_loc0 - smoothed_loc0);
                        double res_cov = (meas_err - smoothed_err);

                        residuals.push_back(res_val * 0.1); // Convert mm to cm for FairShip validation

                        if (res_cov > 1e-9) {
                            pulls.push_back(res_val / std::sqrt(res_cov));
                        } else {
                            pulls.push_back(meas_err > 0.0 ? (res_val / std::sqrt(meas_err)) : 0.0);
                        }

                        true_states_found++;
                        return true;
                    });

                    if (true_states_found > 0) {
                        states_parsed = true;
                        std::reverse(residuals.begin(), residuals.end());
                        std::reverse(pulls.begin(), pulls.end());
                    }
                }
                catch (...) {
                    states_parsed = false;
                }

                if (!states_parsed) {
                    residuals.clear();
                    pulls.clear();
                    for (unsigned int i = 0; i < nMeasurements; ++i) {
                        double res_val = params(Acts::eBoundLoc0);
                        residuals.push_back(res_val * 0.1);
                        pulls.push_back(res_val / 0.012);
                    }
                }

                container->emplace_back(
                    params,
                    cov,
                    surface,
                    chi2,
                    ndof,
                    gctx,
                    residuals,
                    pulls,
                    residuals.size()
                );
            }
        }
        catch (const std::exception& e) {
            std::cout << "WARNING: pushRecoTrack failed or data corrupt: " << e.what() << std::endl;
            return;
        }
    }, pybind11::arg("vectorAddress"), pybind11::arg("gctx"), pybind11::arg("track_idx"), pybind11::arg("containerObj"));

    mex.def("makeBoundTrackParameters", [](std::shared_ptr<const Acts::Surface> surface,
                                           const Acts::BoundVector& params,
                                           const Acts::BoundMatrix& cov) {
        if (!surface) {
            throw std::runtime_error("makeBoundTrackParameters: Surface pointer is null!");
        }
        return Acts::BoundTrackParameters(surface, params, cov, Acts::ParticleHypothesis::pion());
    }, pybind11::arg("surface"), pybind11::arg("params"), pybind11::arg("cov"));


    m.def("getSurface", [](std::shared_ptr<const Acts::TrackingGeometry> geometry, 
                           Acts::GeometryIdentifier geoId) {
        const auto* surface = geometry->findSurface(geoId);
        return surface ? surface->getSharedPtr() : std::shared_ptr<const Acts::Surface>();
    });


    m.def("makePassThroughCalibrator", []() -> std::shared_ptr<ActsExamples::MeasurementCalibrator> {
        return std::make_shared<ActsExamples::PassThroughCalibrator>();
    });

    m.def("getMagneticFieldAt", [](std::shared_ptr<const Acts::MagneticFieldProvider> bField, double x, double y, double z) {
        Acts::Vector3 pos(x, y, z);
        Acts::MagneticFieldContext magCtx;
        auto cache = bField->makeCache(magCtx);
        auto result = bField->getField(pos, cache);
        if (!result.ok()) {
            throw std::runtime_error("Field lookup failed: " + result.error().message());
        }
        return result.value();
    }, py::arg("bField"), py::arg("x"), py::arg("y"), py::arg("z"));


            m.def("fitVertex", [](py::object proxiesObj,
                       std::shared_ptr<const Acts::MagneticFieldProvider> bField,
                       const Acts::GeometryContext& geoCtx,
                       [[maybe_unused]] const Acts::TrackingGeometry& trackingGeometry) {

        Acts::MagneticFieldContext magCtx;


        // Convert Python sequence of proxies into a C++ vector of proxies
        std::vector<ActsExamples::ConstTrackContainer::ConstTrackProxy> proxies;
        std::vector<int> inputTrackIndices;
        inputTrackIndices.reserve(proxies.size());
        try {
            py::sequence seq = proxiesObj;
            proxies.reserve(seq.size());
            for (auto item : seq) {
                proxies.push_back(item.cast<ActsExamples::ConstTrackContainer::ConstTrackProxy>());
            }
        } catch (const std::exception& e) {
            std::cout << "ERROR: fitVertex failed to cast proxies: " << e.what() << std::endl;
            return ActsExamples::VertexContainer();
        }

        if (proxies.size() < 2) {
            std::cout << "WARNING: Vertex fit requires at least 2 tracks." << std::endl;
            return ActsExamples::VertexContainer();
        }


        Acts::Vector3 seedPos = Acts::Vector3::Zero();

        if (proxies.size() == 2) {
            // Build BoundTrackParameters for both proxies to get global positions
            // and momenta evaluated in the geometry context.
            std::vector<Acts::BoundTrackParameters> bparams;
            bparams.reserve(2);
            for (size_t i = 0; i < 2; ++i) {
                const auto& proxy = proxies[i];
                Acts::BoundVector paramsVec = proxy.parameters();
                Acts::BoundMatrix covMat = proxy.covariance();
                auto surfacePtr = proxy.hasReferenceSurface() ? proxy.referenceSurface().getSharedPtr() : std::shared_ptr<const Acts::Surface>();
                if (!surfacePtr) {
                    // fallback perigee at mid-vessel
                    surfacePtr = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3(45000.0, 0.0, 0.0));
                }
                bparams.emplace_back(surfacePtr, paramsVec, std::optional<Acts::BoundMatrix>(covMat), Acts::ParticleHypothesis::pion());
            }

            Acts::Vector3 p1 = bparams[0].position(geoCtx);
            Acts::Vector3 u1 = bparams[0].momentum().normalized();
            Acts::Vector3 p2 = bparams[1].position(geoCtx);
            Acts::Vector3 u2 = bparams[1].momentum().normalized();

            Acts::Vector3 w0 = p1 - p2;
            double a = u1.dot(u1);
            double bdot = u1.dot(u2);
            double c = u2.dot(u2);
            double d = u1.dot(w0);
            double e = u2.dot(w0);

            double denom = a * c - bdot * bdot;

            double t = 0.0, s = 0.0;
            if (std::abs(denom) < 1e-12) {
                // Nearly parallel: project midpoint between p1 and projection of p2 onto line1
                t = (u1.dot(p2 - p1)) / a;
                s = 0.0;
            } else {
                t = (bdot * e - c * d) / denom;
                s = (a * e - bdot * d) / denom;
            }

            Acts::Vector3 c1 = p1 + u1 * t;
            Acts::Vector3 c2 = p2 + u2 * s;
            seedPos = Acts::Vector3(0.5 * (c1.x() + c2.x()), 0.5 * (c1.y() + c2.y()), 0.5 * (c1.z() + c2.z()));
            // Clamp along vessel longitudinal X bounds (20m to 90m)
            if (seedPos.x() < 20000.0) seedPos.x() = 20000.0;
            if (seedPos.x() > 90000.0) seedPos.x() = 90000.0;
        } else {

            Acts::SquareMatrix3 A = Acts::SquareMatrix3::Zero();
            Acts::Vector3 b = Acts::Vector3::Zero();

            for (const auto& proxy : proxies) {

                try {
                    inputTrackIndices.push_back(static_cast<int>(proxy.index()));
                    std::cout<<"Proxy index: "<<proxy.index()<<std::endl;
                } catch (...) {
                    inputTrackIndices.push_back(-1);
                }
                Acts::BoundVector paramsVec = proxy.parameters();
                Acts::BoundMatrix covMat = proxy.covariance();
                auto surfacePtr = proxy.hasReferenceSurface() ? proxy.referenceSurface().getSharedPtr() : std::shared_ptr<const Acts::Surface>();
                if (!surfacePtr) {
                    surfacePtr = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3(45000.0, 0.0, 0.0));
                }
                Acts::BoundTrackParameters btp(surfacePtr, paramsVec, std::optional<Acts::BoundMatrix>(covMat), Acts::ParticleHypothesis::pion());

                Acts::Vector3 p = btp.position(geoCtx);
                Acts::Vector3 n = btp.momentum().normalized();

                Acts::SquareMatrix3 projection = Acts::SquareMatrix3::Identity() - (n * n.transpose());
                A += projection;
                b += projection * p;
            }

            seedPos = A.colPivHouseholderQr().solve(b);

            if (std::isnan(seedPos.x()) || seedPos.x() < 20000.0 || seedPos.x() > 90000.0) {
                seedPos = Acts::Vector3(45000.0, 0.0, 0.0);
            }
        }

        //std::cout << "[VERTEX SEEDER] Analytic two-line PCA Seed: X = " << seedPos.x() << " mm" << std::endl;


        std::vector<Acts::BoundTrackParameters> extractedParams;
        extractedParams.reserve(proxies.size());

        // Perigee surface at seed
        auto perigeeSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(seedPos);

        using StepperDirect = Acts::EigenStepper<>;
        using NavigatorDirect = Acts::DirectNavigator;
        using PropagatorDirect = Acts::Propagator<StepperDirect, NavigatorDirect>;

        StepperDirect stepperDirect(bField);
        NavigatorDirect directNav;
        auto propagatorDirect = std::make_shared<PropagatorDirect>(std::move(stepperDirect), std::move(directNav), Acts::getDefaultLogger("PropagatorDirect", Acts::Logging::INFO));

        for (size_t idx = 0; idx < proxies.size(); ++idx) {
            const auto& proxy = proxies[idx];

            // Start from the proxy's bound parameters
            Acts::BoundVector params = proxy.parameters();
            params[Acts::eBoundTime] = 0.0; // Clear time coordinate
            Acts::BoundMatrix cov = proxy.covariance();
            cov(Acts::eBoundTime, Acts::eBoundTime) = 1e6; // De-weight time completely

            std::shared_ptr<const Acts::Surface> startSurface;
            if (proxy.hasReferenceSurface()) {
                startSurface = proxy.referenceSurface().getSharedPtr();
            } else {

                std::cout << "[VTX DEBUG] Proxy missing reference surface; using proxy/seed position as start surface fallback for track " << idx << std::endl;
                Acts::Vector3 proxyPos = seedPos;

                try {
                    Acts::BoundVector pv = proxy.parameters();
                    // Reconstruct a global position along the track direction from the proxy parameters
                    double phi = pv[Acts::eBoundPhi];
                    double theta = pv[Acts::eBoundTheta];
                    // approximate short displacement along momentum to create a starting point
                    Acts::Vector3 dir(std::cos(phi)*std::sin(theta), std::sin(phi)*std::sin(theta), std::cos(theta));
                    proxyPos = seedPos + dir * 10.0; // 10 mm offset from seed as a guess
                } catch (...) {
                    proxyPos = seedPos;
                }

                startSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(proxyPos);
            }
            Acts::BoundTrackParameters startBP(startSurface, params, std::optional<Acts::BoundMatrix>(cov), Acts::ParticleHypothesis::pion());

            // Propagate from startBP to the perigee surface using the direct propagator
            using PropOptions = PropagatorDirect::Options<>;
            PropOptions popts(geoCtx, magCtx);
            popts.pathLimit = 1e6; // allow long propagation
            // Provide navigator hints: start surface and an ordered list containing the perigee target
            popts.navigation.startSurface = startSurface.get();
            popts.navigation.externalSurfaces.clear();
            popts.navigation.appendExternalSurface(*perigeeSurface);


            try {
                auto startPos = startBP.position(geoCtx);
                auto startMom = startBP.momentum();
                Acts::Vector3 toTarget = seedPos - startPos;
                double dp = startMom.dot(toTarget);
                popts.direction = (dp >= 0.0) ? Acts::Direction::Forward() : Acts::Direction::Backward();
                //std::cout << "[VTX DEBUG] Chosen propagation direction for track " << idx << " dp=" << dp << " -> " << (dp >= 0.0 ? "Forward" : "Backward") << std::endl;
            } catch (...) {

                std::cout << "[VTX DEBUG] Failed to compute propagation direction heuristic for track " << idx << ", using default" << std::endl;
            }

            try {
                auto pres = propagatorDirect->propagate(startBP, *perigeeSurface, popts);
                if (pres.ok() && pres.value().endParameters.has_value()) {
                    // Use propagated parameters at perigee
                    extractedParams.push_back(pres.value().endParameters.value());
                    continue;
                } else {
                    std::cout << "[VTX DEBUG] Direct propagation to perigee failed for track " << idx << ": " << (pres.ok() ? "no endParameters" : pres.error().message()) << std::endl;
                    // Try the opposite propagation direction as a fallback
                    popts.direction = Acts::Direction::Backward();
                    std::cout << "[VTX DEBUG] Retrying direct propagation in backward direction for track " << idx << std::endl;
                    auto pres2 = propagatorDirect->propagate(startBP, *perigeeSurface, popts);
                    if (pres2.ok() && pres2.value().endParameters.has_value()) {
                        extractedParams.push_back(pres2.value().endParameters.value());
                        continue;
                    } else {
                        std::cout << "[VTX DEBUG] Backward propagation also failed for track " << idx << ": " << (pres2.ok() ? "no endParameters" : pres2.error().message()) << std::endl;
                    }
                }
            } catch (const std::exception& e) {
                std::cout << "[VTX DEBUG] Exception during direct propagation for track " << idx << ": " << e.what() << std::endl;
            }

            // Fallback: build approximate perigee parameters analytically (straight-line)
            Acts::Vector3 globalPos = proxy.hasReferenceSurface() ? proxy.referenceSurface().center(geoCtx) : seedPos;
            Acts::Vector3 globalMom = proxy.momentum();
            double charge = proxy.charge();

            double delta_x = seedPos.x() - globalPos.x();
            Acts::Vector3 extrapolatedPos = globalPos + globalMom.normalized() * (delta_x / globalMom.normalized().x());

            auto localResult = perigeeSurface->globalToLocal(geoCtx, extrapolatedPos, globalMom.normalized());
            if (localResult.ok()) {
                Acts::BoundVector perigeeParams = Acts::BoundVector::Zero();
                perigeeParams[Acts::eBoundLoc0]   = localResult.value()[Acts::eBoundLoc0];
                perigeeParams[Acts::eBoundLoc1]   = localResult.value()[Acts::eBoundLoc1];
                perigeeParams[Acts::eBoundPhi]    = std::atan2(globalMom.y(), globalMom.x());
                perigeeParams[Acts::eBoundTheta]  = std::acos(globalMom.z() / globalMom.norm());
                perigeeParams[Acts::eBoundQOverP] = charge / globalMom.norm();
                perigeeParams[Acts::eBoundTime]   = 0.0;

                Acts::BoundMatrix perigeeCov = Acts::BoundMatrix::Identity();
                Acts::BoundMatrix originalCov = proxy.covariance();
                perigeeCov(Acts::eBoundLoc0, Acts::eBoundLoc0) = originalCov(Acts::eBoundLoc0, Acts::eBoundLoc0) * 10.0;
                perigeeCov(Acts::eBoundLoc1, Acts::eBoundLoc1) = originalCov(Acts::eBoundLoc1, Acts::eBoundLoc1) * 10.0;
                perigeeCov(Acts::eBoundPhi,   Acts::eBoundPhi)   = originalCov(Acts::eBoundPhi,   Acts::eBoundPhi);
                perigeeCov(Acts::eBoundTheta, Acts::eBoundTheta) = originalCov(Acts::eBoundTheta, Acts::eBoundTheta);
                perigeeCov(Acts::eBoundQOverP,Acts::eBoundQOverP)= originalCov(Acts::eBoundQOverP,Acts::eBoundQOverP);
                perigeeCov(Acts::eBoundTime,  Acts::eBoundTime)  = 1e6;

                extractedParams.emplace_back(perigeeSurface, perigeeParams, std::optional<Acts::BoundMatrix>(perigeeCov), Acts::ParticleHypothesis::pion());
            } else {
                // As last resort, use the startBP unchanged
                extractedParams.push_back(startBP);
            }
        }


        Acts::Vertex seedVertex(seedPos);
        Acts::SquareMatrix4 seedCov = Acts::SquareMatrix4::Identity() * 1000000.0;
        seedVertex.setFullCovariance(seedCov);


        Acts::VertexingOptions vtxOptions(geoCtx, magCtx, seedVertex);


        std::shared_ptr<const Acts::MagneticFieldProvider> usedB = bField;
        try {
            auto cacheCheck = bField->makeCache(magCtx);
            auto fieldRes = bField->getField(Acts::Vector3(seedPos.x(), seedPos.y(), seedPos.z()), cacheCheck);
            if (fieldRes.ok()) {
                double bmag = fieldRes.value().norm();
                if (bmag < 1e-6) {
                    std::cout << "[VTX DEBUG] Local B-field magnitude " << bmag << " T is negligible; using NullBField for linearization." << std::endl;
                    usedB = std::make_shared<Acts::NullBField>();
                }
            }
        } catch (...) {
            std::cout << "[VTX DEBUG] Failed to query B-field; proceeding with provided field for linearization." << std::endl;
        }

        using Stepper = Acts::EigenStepper<>;
        using Navigator = Acts::VoidNavigator; // use VoidNavigator for debugging/stability
        using Propagator = Acts::Propagator<Stepper, Navigator>;


        // VoidNavigator (geometry-free) so we can isolate navigation issues.
        Stepper stepper(usedB);
        Navigator voidNav;
        auto propagator = std::make_shared<Propagator>(std::move(stepper), std::move(voidNav), Acts::getDefaultLogger("Propagator", Acts::Logging::INFO));


        Acts::HelicalTrackLinearizer::Config ltConfig;
        ltConfig.bField = usedB;
        ltConfig.propagator = propagator;

        ltConfig.targetTolerance = 1e-4;

        Acts::HelicalTrackLinearizer linearizer(ltConfig);

        Acts::FullBilloirVertexFitter::Config fCfg;
        fCfg.maxIterations = 200; // allow more iterations for convergence

        fCfg.extractParameters.connect(gTrackExtractor);
        fCfg.trackLinearizer.connect<&Acts::HelicalTrackLinearizer::linearizeTrack>(&linearizer);


        auto vlogger = Acts::getDefaultLogger("Billoir", Acts::Logging::INFO);
        Acts::FullBilloirVertexFitter fitter(fCfg, std::move(vlogger));


        std::vector<Acts::InputTrack> inputTracks;
        auto* heapParams = new std::vector<Acts::BoundTrackParameters>(std::move(extractedParams));
        // Record globally so Python can fetch this address if needed
        g_lastHeapParams = heapParams;
        inputTracks.reserve(heapParams->size());
        for (size_t i = 0; i < heapParams->size(); ++i) {
            inputTracks.emplace_back(&((*heapParams)[i]));
        }

        auto fieldCache = bField->makeCache(magCtx);


        auto tstart = std::chrono::steady_clock::now();
        auto result = fitter.fit(inputTracks, vtxOptions, fieldCache);
        auto tend = std::chrono::steady_clock::now();
        (void)tstart; (void)tend; // silence unused warnings if any

        ActsExamples::VertexContainer vertexContainer;
        if (result.ok()) {
           auto v = result.value();
           try {
               py::object pyVtx = py::cast(v);
               pyVtx.attr("_input_track_indices") = py::cast(inputTrackIndices);

               try {
                   if (heapParams) {

                       py::capsule cap(reinterpret_cast<void*>(heapParams), [](void* p){ delete static_cast<std::vector<Acts::BoundTrackParameters>*>(p); });
                       py::module_::import("acts").attr("_last_extracted_params_capsule") = cap;
                       py::module_::import("acts").attr("_last_extracted_params_addr") = reinterpret_cast<uintptr_t>(heapParams);
                       g_lastHeapParams = heapParams;

                       pyVtx.attr("_extracted_params_addr") = reinterpret_cast<uintptr_t>(heapParams);
                   }
               } catch (...) {}
           } catch (...) {
               // best-effort attach; ignore failures
           }

           // Diagnostic: Fit quality (chi2, ndof)
           auto fq = v.fitQuality();
           //std::cout << "[VTX DEBUG] FitQuality: chi2=" << fq.first << " ndof=" << fq.second << std::endl;

           // Diagnostic: seed vs final
           //std::cout << "[VTX DEBUG] seedPos(mm): (" << seedPos.x() << ", " << seedPos.y() << ", " << seedPos.z() << ")" << std::endl;
           auto vpos = v.position();
           //std::cout << "[VTX DEBUG] fittedPos(mm): (" << vpos.x() << ", " << vpos.y() << ", " << vpos.z() << ")" << std::endl;


           // Per-track fitted params
           try {
               size_t nt = v.tracks().size();
             //  std::cout << "[VTX DEBUG] vtx.tracks size = " << nt << std::endl;
               for (size_t i = 0; i < nt; ++i) {
                   const auto& tav = v.tracks()[i];
                   auto bp = tav.fittedParams;
                   auto ppos = bp.position(geoCtx);
                   auto pmom = bp.momentum();
               //    std::cout << "[VTX DEBUG] track " << i << " fitted pos(mm): (" << ppos.x() << ", " << ppos.y() << ", " << ppos.z() << ") ";
               //    std::cout << "mom: (" << pmom.x() << ", " << pmom.y() << ", " << pmom.z() << ")" << std::endl;
                   auto ocov = bp.covariance();
                   if (ocov.has_value()) {
                       const auto& bcov = ocov.value();
                //       std::cout << "[VTX DEBUG] track " << i << " cov00=" << bcov(0,0) << std::endl;
                   }
                   // Extrapolate original (extracted) track parameters to the fitted vertex and
                   // print the distance between the extrapolated position and the fitted vertex.
                   try {

                       if (heapParams && i < heapParams->size()) {
                           const auto& start = (*heapParams)[i];
                           // Create a perigee surface at the fitted vertex to receive the extrapolation
                           auto targetSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(vpos);
                           using PropagatorT = Propagator;
                           using PropOptions = PropagatorT::Options<>;
                           PropOptions popts(geoCtx, magCtx);
                           popts.pathLimit = 1e6; // large path limit
                           auto pres = propagator->propagate(start, *targetSurface, popts);
                           if (pres.ok()) {
                               if (pres.value().endParameters.has_value()) {
                                   auto endPar = pres.value().endParameters.value();
                                   auto endPos = endPar.position(geoCtx);
                                   double dx = endPos.x() - vpos.x();
                                   double dy = endPos.y() - vpos.y();
                                   double dz = endPos.z() - vpos.z();
                                   double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                  //                 std::cout << "[VTX DEBUG] extrapolated distance to vertex for track " << i << " = " << dist << " mm" << std::endl;
                               } else {
                  //                 std::cout << "[VTX DEBUG] extrapolation returned no endParameters for track " << i << std::endl;
                               }
                           } else {
                               //std::cout << "[VTX DEBUG] extrapolation failed for track " << i << " (" << pres.error().message() << ")" << std::endl;
                           }
                       }
                   } catch (const std::exception& e) {
                       std::cout << "[VTX DEBUG] extrapolation exception for track " << i << ": " << e.what() << std::endl;
                   }
               }
           } catch (const std::exception& e) {
               std::cout << "[VTX DEBUG] Failed to print per-track fitted params: " << e.what() << std::endl;
           }

           vertexContainer.push_back(std::move(v));
        } else {
           std::cout << "DEBUG: Vertex Fit failed! Error: " << result.error().message()
                     << " (category: " << result.error().category().name() << ")"
                     << " at Seed X: " << seedPos.x() << std::endl;
        }
        return vertexContainer;
    }, py::arg("proxies"), py::arg("bField"), py::arg("geoCtx"), py::arg("trackingGeometry"));




    py::class_<Acts::Vertex>(m, "Vertex")
        .def(py::init<>())
        .def("position", [](const Acts::Vertex& vtx) -> Acts::Vector3 {
            return vtx.position();
        })
        .def("fullPosition", [](const Acts::Vertex& vtx) -> Acts::Vector4 {
            return vtx.fullPosition();
        })
        .def("covariance", [](const Acts::Vertex& vtx) -> Acts::SquareMatrix3 {
            return vtx.covariance();
        })
                .def("fullCovariance", [](const Acts::Vertex& vtx) -> Acts::SquareMatrix4 {
            return vtx.fullCovariance();
        })

        .def("tracks", [](const Acts::Vertex& vtx) {
            py::list out;
                        py::object SimpleNS = py::module_::import("types").attr("SimpleNamespace");
            // try to retrieve matched indices attached by pushRecoVertex
            py::list matched;
            try {
                py::object pyVtx = py::cast(vtx);
                matched = pyVtx.attr("_matched_track_indices");
            } catch (...) {
                matched = py::list();
            }

                       const auto& tracks = vtx.tracks();
            for (size_t i = 0; i < tracks.size(); ++i) {
                const auto& tv = tracks[i];
                py::object item = SimpleNS();
                // momentum as tuple (fitted / smoothed at vertex)
                auto mom = tv.fittedParams.momentum();
                item.attr("momentum") = py::make_tuple(mom.x(), mom.y(), mom.z());
                // attempt to expose fitted position (may throw if unavailable)
                Acts::GeometryContext defGeo = Acts::GeometryContext::dangerouslyDefaultConstruct();
                try {
                    auto pos = tv.fittedParams.position(defGeo);
                    item.attr("position") = py::make_tuple(pos.x(), pos.y(), pos.z());
                } catch (...) {
                    item.attr("position") = py::none();
                }

                // expose original (pre-vertex) parameters when available
                try {
                    const Acts::BoundTrackParameters* origBp = nullptr;
                    int originalIndex = -1;
                    // Try OwnedInput wrapper first
                    try {
                        const OwnedInput* oi = tv.originalParams.as<OwnedInput>();
                        if (oi) {
                            origBp = &oi->params;
                            originalIndex = oi->originalIndex;
                        }
                    } catch (...) {}
                    // Fallback to direct BoundTrackParameters stored in originalParams
                    if (!origBp) {
                        try {
                            const Acts::BoundTrackParameters* bp = tv.originalParams.as<Acts::BoundTrackParameters>();
                            if (bp) origBp = bp;
                        } catch (...) {}
                    }

                    if (origBp) {
                        try {
                            auto op = origBp->position(defGeo);
                            auto om = origBp->momentum();
                            item.attr("originalPosition") = py::make_tuple(op.x(), op.y(), op.z());
                            item.attr("originalMomentum") = py::make_tuple(om.x(), om.y(), om.z());
                            item.attr("originalIndex") = originalIndex;
                        } catch (...) {
                            item.attr("originalPosition") = py::none();
                            item.attr("originalMomentum") = py::none();
                            item.attr("originalIndex") = -1;
                        }
                    } else {
                        item.attr("originalPosition") = py::none();
                        item.attr("originalMomentum") = py::none();
                        item.attr("originalIndex") = -1;
                    }
                } catch (...) {
                    item.attr("originalPosition") = py::none();
                    item.attr("originalMomentum") = py::none();
                    item.attr("originalIndex") = -1;
                }

                // try to compute and expose a 6x6 position+momentum covariance derived from the fitted parameter covariance
                try {
                    auto ocov = tv.fittedParams.covariance();
                    if (ocov.has_value()) {
                        const auto& pCov = ocov.value();
                        int nParams = static_cast<int>(pCov.rows());
                        // read parameter vector
                        auto params = tv.fittedParams.parameters();
                        // Build numerical Jacobian J6: rows [x,y,z,px,py,pz] x cols [params]
                        Eigen::MatrixXd J6(6, nParams);
                        Acts::GeometryContext geoCtx = Acts::GeometryContext::dangerouslyDefaultConstruct();
                        // base values
                        auto baseMom = tv.fittedParams.momentum();
                        auto basePos = tv.fittedParams.position(geoCtx);
                        for (int k = 0; k < nParams; ++k) {
                            double pv = params.coeff(k);
                            double eps = std::max(1e-8, std::abs(pv) * 1e-6);
                            auto newParams = params;
                            newParams.coeffRef(k) = pv + eps;
                            std::shared_ptr<const Acts::Surface> surfPtr;
                            try {
                                // attempt to obtain reference surface if available
                                try { surfPtr = tv.fittedParams.referenceSurface().getSharedPtr(); } catch (...) { surfPtr.reset(); }
                            } catch (...) { surfPtr.reset(); }
                            Acts::BoundTrackParameters temp(surfPtr, newParams, std::optional<Acts::BoundMatrix>(), tv.fittedParams.particleHypothesis());
                            // perturbed values
                            auto pm = temp.momentum();
                            auto pos = temp.position(geoCtx);
                            J6(0, k) = (pos.x() - basePos.x()) / eps;
                            J6(1, k) = (pos.y() - basePos.y()) / eps;
                            J6(2, k) = (pos.z() - basePos.z()) / eps;
                            J6(3, k) = (pm.x() - baseMom.x()) / eps;
                            J6(4, k) = (pm.y() - baseMom.y()) / eps;
                            J6(5, k) = (pm.z() - baseMom.z()) / eps;
                        }

                        Eigen::MatrixXd Pcov(nParams, nParams);
                        for (int ii = 0; ii < nParams; ++ii) for (int jj = 0; jj < nParams; ++jj) Pcov(ii, jj) = pCov(ii, jj);

                        Eigen::MatrixXd cov6 = J6 * Pcov * J6.transpose();

                        // expose as Eigen 6x6 (numpy) to Python
                        item.attr("covariance") = py::cast(cov6);
                    } else {
                        item.attr("covariance") = py::none();
                    }
                } catch (...) {
                    item.attr("covariance") = py::none();
                }

                int tidx = -1;
                try {
                    if (i < static_cast<size_t>(py::len(matched))) tidx = matched[i].cast<int>();
                } catch (...) { tidx = -1; }


                if (tidx == -1) {
                    try {
                        uintptr_t addr = reinterpret_cast<uintptr_t>(std::addressof(vtx));
                        auto it = g_vertexMatchedIndices.find(addr);
                        if (it != g_vertexMatchedIndices.end()) {
                            const auto& cached = it->second;
                            if (i < cached.size()) tidx = cached[i];
                        }
                    } catch (...) {}
                }

                item.attr("trackIndex") = tidx;
                out.append(item);
            }
            return out;
        });

    py::class_<ActsExamples::RecoVertex>(m, "RecoVertex")
        .def(py::init<const Acts::Vertex&>())
        .def_property_readonly("x", &ActsExamples::RecoVertex::x)
        .def_property_readonly("y", &ActsExamples::RecoVertex::y)
        .def_property_readonly("z", &ActsExamples::RecoVertex::z)
        .def_property_readonly("chi2", &ActsExamples::RecoVertex::chi2)
        .def_property_readonly("ndof", &ActsExamples::RecoVertex::nDoF)
        .def("trackIds", &ActsExamples::RecoVertex::trackIds)
        .def("trackPx", &ActsExamples::RecoVertex::trackPx)
        .def("trackPy", &ActsExamples::RecoVertex::trackPy)
        .def("trackPz", &ActsExamples::RecoVertex::trackPz)
        .def("trackX", &ActsExamples::RecoVertex::trackX)
        .def("trackY", &ActsExamples::RecoVertex::trackY)
        .def("trackZ", &ActsExamples::RecoVertex::trackZ);

    m.def("pushRecoVertex", [](long vectorAddr,
                               const Acts::Vertex& vtx,
                               const ActsExamples::ConstTrackContainer& outputTracks,
                               long inputParamsAddr) {

       auto* vertexVector = reinterpret_cast<std::vector<ActsExamples::RecoVertex>*>(vectorAddr);
       if (!vertexVector) {
           throw std::runtime_error("CRITICAL: Null vertex vector pointer provided!");
       }

       ActsExamples::RecoVertex recoVtx(vtx);
       recoVtx.clearTrackIds();


       std::vector<int> inputIndices;
       try {
           py::object pyVtx = py::cast(vtx);
       if (py::hasattr(pyVtx, "_input_track_indices")) {
               inputIndices = py::cast<std::vector<int>>(pyVtx.attr("_input_track_indices"));
           }
       } catch (...) {
           inputIndices.clear();
       }


       std::vector<Acts::BoundTrackParameters>* inputParamsPtr = nullptr;
       if (inputParamsAddr != 0) {
           inputParamsPtr = reinterpret_cast<std::vector<Acts::BoundTrackParameters>*>(inputParamsAddr);
       }

       std::vector<int> matchedIndices;
       matchedIndices.reserve(vtx.tracks().size());

       for (size_t trk_idx = 0; trk_idx < vtx.tracks().size(); ++trk_idx) {

           Acts::Vector3 vtxTrackMom = vtx.tracks()[trk_idx].fittedParams.momentum();

           int matched_index = -1;


           if (inputParamsPtr) {
               const auto& tav = vtx.tracks()[trk_idx];
               const Acts::BoundTrackParameters* origPtr = nullptr;
               try {
                   const OwnedInput* oi = tav.originalParams.as<OwnedInput>();
                   if (oi) origPtr = &oi->params;
               } catch (...) {}
               if (!origPtr) {
                   try {
                       const Acts::BoundTrackParameters* bp = tav.originalParams.as<Acts::BoundTrackParameters>();
                       if (bp) origPtr = bp;
                   } catch (...) {}
               }
               if (origPtr) {
                   for (size_t j = 0; j < inputParamsPtr->size(); ++j) {
                       if (&((*inputParamsPtr)[j]) == origPtr) {
                           matched_index = static_cast<int>(j);
                           break;
                       }
                   }
               }
           }

           matchedIndices.push_back(matched_index);

           if (matched_index != -1) {
               recoVtx.addTrackId(matched_index);
           } else {

           }
       }


       try {
           uintptr_t addr = reinterpret_cast<uintptr_t>(std::addressof(vtx));
           g_vertexMatchedIndices[addr] = matchedIndices;
       } catch (...) {}


       try {
           py::object pyVtx = py::cast(vtx);
           pyVtx.attr("_matched_track_indices") = py::cast(matchedIndices);
       } catch (...) {

       }

       vertexVector->push_back(std::move(recoVtx));

    }, py::arg("vectorAddr"), py::arg("vtx"), py::arg("outputTracks"), py::arg("inputParamsAddr") = 0);

    m.def("getTrackResiduals", [](const ActsExamples::ConstTrackProxy& track) {
        std::vector<double> residuals;
        for (const auto& state : track.trackStatesReversed()) {
            if (state.hasUncalibratedSourceLink()) {
                double pred = state.predicted()[Acts::eBoundLoc0];
                double meas = state.calibrated<1>()[0];
                residuals.push_back(meas - pred);
            }
        }
        return residuals;
    });




    m.def("getTrackResiduals", [](const ActsExamples::ConstTrackProxy& track) {
        std::vector<double> residuals;
        for (const auto& state : track.trackStatesReversed()) {
            if (state.hasUncalibratedSourceLink()) {
                double pred = state.predicted()[Acts::eBoundLoc0];
                double meas = state.calibrated<1>()[0];
                residuals.push_back(meas - pred);
            }
        }
        return residuals;
    });

    m.def("fitTrackDAF", [](
        const ActsExamples::MeasurementContainer& measurements,
        const std::vector<unsigned int>& indices,
        py::object initialParamsObj,
        ActsExamples::TrackContainer& outputTracks,
        std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
        std::shared_ptr<const Acts::MagneticFieldProvider> bField,
        int daf_max_iter /*=6*/,
        py::object daf_anneal /*=py::none()*/, 
        double daf_cutoff /*=9.0*/, 
        double daf_min_variance /*=1e-4*/) {
        
        const auto& initialParams = initialParamsObj.cast<const Acts::BoundTrackParameters&>();

        // Build annealing schedule
        std::vector<double> annealSchedule;
        try {
            if (!daf_anneal.is_none()) {
                annealSchedule = daf_anneal.cast<std::vector<double>>();
            }
        } catch(...) {}
        if (annealSchedule.empty()) {
            double bStart = 100.0, bFinal = 0.1;
            unsigned int nSteps = 10;
            for (unsigned int i = 0; i < nSteps; ++i) {
                annealSchedule.push_back(bStart * pow(bFinal / bStart, double(i)/(nSteps-1)));
            }
        }

        // Create DAF fitter and run
        ActsExamples::DeterministicAnnealingFitter::Config cfg;
        cfg.annealingSchedule = annealSchedule;
        cfg.maxIterations = daf_max_iter;
        cfg.gateThreshold = daf_cutoff;
        cfg.minBaseVariance = daf_min_variance;
        cfg.convergenceTolerance = 1e-3;
        cfg.priorVarianceScale = 0.8;

        ActsExamples::DeterministicAnnealingFitter fitter(cfg, Acts::Logging::INFO);
        auto result = fitter.fit(measurements, indices, initialParams, tGeometry, bField, outputTracks);

        // Export diagnostics to Python module (for debugging)
        try {
            auto acts_mod = py::module_::import("acts");
            acts_mod.attr("_last_fit_diagnostics") = result.diagnostics;
        } catch(...) {}

        return result.success;

    }, py::arg("measurements"), py::arg("indices"), py::arg("initialParams"), py::arg("outputTracks"), py::arg("trackingGeometry"), py::arg("magneticField"), py::arg("daf_max_iter") = 6, py::arg("daf_anneal") = py::none(), py::arg("daf_cutoff") = 9.0, py::arg("daf_min_variance") = 1e-4);

    //If used with strawHits requires drift to be set to the correct side//
    m.def("fitTrack", [](
        const ActsExamples::MeasurementContainer& measurements,
        const std::vector<unsigned int>& indices,
        py::object initialParamsObj,
        ActsExamples::TrackContainer& outputTracks,
        std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
        std::shared_ptr<const Acts::MagneticFieldProvider> bField) {
   

        Acts::GeometryContext geoCtx = Acts::GeometryContext::dangerouslyDefaultConstruct();
        Acts::MagneticFieldContext magCtx;
        Acts::CalibrationContext calibCtx;
        const auto& initialParams = initialParamsObj.cast<const Acts::BoundTrackParameters&>();

 
        std::vector<Acts::SourceLink> concreteSourceLinks;
        concreteSourceLinks.reserve(indices.size());
        for (auto idx : indices) {
            const auto& meas = measurements.getMeasurement(idx);
            auto* surface = tGeometry->findSurface(meas.geometryId());
            if (surface) {
                concreteSourceLinks.push_back(Acts::SourceLink(ActsExamples::IndexSourceLink{meas.geometryId(), idx}));
            }
        }

        Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions;
        
        auto accessor = [tg = tGeometry.get()](const Acts::SourceLink& sl) -> const Acts::Surface* {
            auto geoId = sl.template get<ActsExamples::IndexSourceLink>().geometryId();
            const auto* surf = tg->findSurface(geoId);
            if (surf) {
                const_cast<Acts::Surface*>(surf)->assignGeometryId(geoId);
            }
            return surf;
        };
        
        auto calibrator = [&measurements](const Acts::GeometryContext&, 
                                         const Acts::CalibrationContext&, 
                                         const Acts::SourceLink& sl, 
                                         Acts::TrackStateProxy<Acts::VectorMultiTrajectory, 6, false> ts) {
            const auto& islink = sl.template get<ActsExamples::IndexSourceLink>();
            const auto& meas = measurements.getMeasurement(islink.index());
            auto measDim = meas.size();
            if (measDim == 1) {
                ts.allocateCalibrated(1);
                ts.template calibrated<1>() = meas.parameters();
                ts.template calibratedCovariance<1>() = meas.covariance();
            } else if (measDim == 2) {
                ts.allocateCalibrated(2);
                ts.template calibrated<2>() = meas.parameters();
                ts.template calibratedCovariance<2>() = meas.covariance();
            }
            Acts::SourceLink slCopy = sl;
            ts.setUncalibratedSourceLink(std::move(slCopy));

        };

        extensions.surfaceAccessor.connect(accessor);
        extensions.calibrator.connect(calibrator);

        Acts::GainMatrixUpdater updater;
        Acts::GainMatrixSmoother smoother;
   
        extensions.updater.template connect<
            &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(&updater);
   
        extensions.smoother.template connect<
            &Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(&smoother);

        auto outlierFinder = [](Acts::TrackStateProxy<Acts::VectorMultiTrajectory, 6, true>) -> bool {
            return false;
        };
        extensions.outlierFinder.connect(outlierFinder);

        Acts::Navigator::Config navCfg{tGeometry};
        navCfg.resolvePassive = true;
        navCfg.resolveMaterial = true;
        navCfg.resolveSensitive = true;
    
        Acts::Navigator navigator(navCfg, Acts::getDefaultLogger("Navigator", Acts::Logging::INFO));
        auto logger = Acts::getDefaultLogger("KalmanFitter", Acts::Logging::INFO);

        Acts::DirectNavigator navigator2;

        using Stepper = Acts::EigenStepper<>;
        using Propagator = Acts::Propagator<Stepper, Acts::DirectNavigator>;
        using Fitter = Acts::KalmanFitter<Propagator, Acts::VectorMultiTrajectory>;

        Stepper stepper(bField);
        Propagator propagator(std::move(stepper), std::move(navigator2));
        Fitter fitter(std::move(propagator), std::move(logger));

        Acts::PropagatorPlainOptions pOptions(geoCtx, magCtx);
       
        Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory> options(
            geoCtx, magCtx, std::ref(calibCtx), extensions,
            pOptions, 
            &initialParams.referenceSurface());
   
        std::vector<const Acts::Surface*> surfaceSequence;
        for (auto idx : indices) {
            const auto& meas = measurements.getMeasurement(idx);
            auto* surf = tGeometry->findSurface(meas.geometryId());
            if (surf) {
                surfaceSequence.push_back(surf);
            }
        }

        options.multipleScattering = true;
        options.energyLoss = true;
        options.referenceSurfaceStrategy = Acts::TrackExtrapolationStrategy::first;

        auto result = fitter.fit(concreteSourceLinks.begin(), concreteSourceLinks.end(), initialParams, options, surfaceSequence, outputTracks);
        if (result.ok()) {
            const auto& track = result.value();
//            std::cout << "=== FIT SUCCESSFUL ===" << std::endl;
//            std::cout << "Measurements added: " << track.nMeasurements() << std::endl;
//            std::cout << "Holes: " << track.nHoles() << std::endl;
//            std::cout << "Chi2: " << track.chi2() << std::endl;
//            std::cout << "Final Q/P: " << track.parameters()[Acts::eBoundQOverP] << std::endl;
        } else {
//            std::cout << "Fit failed: " << result.error().message() << std::endl;
        }
          return result.ok();
    }, py::arg("measurements"), py::arg("indices"), py::arg("initialParams"),py::arg("outputTracks"), py::arg("trackingGeometry"), py::arg("magneticField"));


    m.def("makeIndexSourceLink", [](Acts::GeometryIdentifier geoId, std::size_t index) {
        return ActsExamples::IndexSourceLink{geoId, static_cast<unsigned int>(index)};
    });

    m.def("createSourceLinks", [](const ActsExamples::MeasurementContainer& measurements,
                                  const std::vector<unsigned int>& indices) {
        std::vector<ActsExamples::IndexSourceLink> sourceLinks;
        sourceLinks.reserve(indices.size());
        for (auto idx : indices) {
            const auto& meas = measurements.getMeasurement(idx);
            sourceLinks.push_back(ActsExamples::IndexSourceLink{meas.geometryId(), idx});
        }
        return sourceLinks;
    });


   
    mex.def("getTrackParameters", [](const ActsExamples::TrackContainer& container) {
        std::vector<Acts::BoundTrackParameters> params;
        for (const auto& track : container) {
            if (track.hasReferenceSurface()) {
                params.emplace_back(
                    track.referenceSurface().getSharedPtr(),
                    track.parameters(),
                    track.covariance(),
                    Acts::ParticleHypothesis::muon()
                );
            }
        }
        return params;
    });
   
    m.def("createTargetSurface", [](double z) -> std::shared_ptr<Acts::Surface>{
        auto transform = Acts::Transform3(Acts::Translation3(z, 0.0, 0.0));
        return Acts::Surface::makeShared<Acts::PerigeeSurface>(transform);
    }, py::arg("z"));

    m.def("createTrackParameters", [](double gx, double gy, double gz,
                                      double px, double py, double pz,
                                      double charge,
                                      std::shared_ptr<const Acts::Surface> surface,
                                      const std::vector<double>& covVec,
                                      const Acts::GeometryContext& gctx) {
        Acts::Vector3 globalPos(gx, gy, gz);
        Acts::Vector3 mom(px, py, pz);

        auto localPosRes = surface->globalToLocal(gctx, globalPos, mom);
        double loc0 = localPosRes.ok() ? localPosRes.value()[0] : 0.0;
        double loc1 = localPosRes.ok() ? localPosRes.value()[1] : 0.0;

        double phi = Acts::VectorHelpers::phi(mom);
        double theta = Acts::VectorHelpers::theta(mom);

        Acts::BoundVector params = Acts::BoundVector::Zero();
        params[Acts::eBoundLoc0] = loc0;
        params[Acts::eBoundLoc1] = loc1;
        params[Acts::eBoundPhi] = phi;
        params[Acts::eBoundTheta] = theta;
        params[Acts::eBoundQOverP] = charge / (mom.norm() + 1e-9);
        params[Acts::eBoundTime] = 0.0;

        Acts::BoundMatrix cov = Acts::BoundMatrix::Identity() * 0.1;
        if (covVec.size() == 36) {
            Eigen::Map<const Acts::BoundMatrix> covMap(covVec.data());
            cov = covMap;
        }

        return Acts::BoundTrackParameters(surface, params, cov, Acts::ParticleHypothesis::muon());
    }, py::arg("x"), py::arg("y"), py::arg("z"),
       py::arg("px"), py::arg("py"), py::arg("pz"),
       py::arg("charge"), py::arg("surface"), py::arg("cov"), py::arg("geo_ctx"));

    m.def("createPlaneSurface", [](Acts::Vector3 center, Acts::Vector3 normal) -> std::shared_ptr<Acts::Surface> {
        Acts::Transform3 transform{Acts::Translation3{center}};
        if (normal.norm() > 1e-6) {
            Acts::Vector3 n = normal.normalized();
            auto rotation = Acts::RotationMatrix3(
                Eigen::Quaternion<double>::FromTwoVectors(Acts::Vector3::UnitZ(), n)
            );
            transform.rotate(rotation);
        }
        return Acts::Surface::makeShared<Acts::PlaneSurface>(transform, nullptr);
    }, py::arg("center"), py::arg("normal"));

    {
        using Builder = ActsExamples::StrawtubeBuilder;
        auto b = py::class_<Builder, ActsExamples::Detector, std::shared_ptr<Builder>>(mex, "StrawtubeBuilder")
                     .def(py::init<const Builder::Config&>(), py::arg("config"))
                     .def("layers", &Builder::layers);
        auto c = py::class_<Builder::Config>(b, "Config").def(py::init<>());
        ACTS_PYTHON_STRUCT(c, fileName, logLevel, layerLogLevel);
    }

    {
        using Detector = ActsExamples::StrawtubeDetector;
        auto d = py::class_<Detector, ActsExamples::Detector, std::shared_ptr<Detector>>(mex, "StrawtubeDetector")
                     .def(py::init<const Detector::Config&>(), py::arg("config"));
        auto c = py::class_<Detector::Config>(d, "Config").def(py::init<>());
        ACTS_PYTHON_STRUCT(c, fileName, logLevel);
    }

    py::class_<HGCBuilder::Config>(mex, "HGCBuilderConfig")
        .def(py::init<>())
        .def_readwrite("fileName", &HGCBuilder::Config::fileName)
        .def_readwrite("logLevel", &HGCBuilder::Config::logLevel);
   
    py::class_<HGCBuilder, ActsExamples::Detector, std::shared_ptr<HGCBuilder>>(mex, "HGCBuilder")
        .def(py::init<const HGCBuilder::Config&>())
        .def("layers", &HGCBuilder::layers);

    auto hgcDetector = py::class_<HGCDetector, ActsExamples::Detector, std::shared_ptr<HGCDetector>>(mex, "HGCDetector")
        .def(py::init<const HGCDetector::Config&>())
        .def("trackingGeometry", &HGCDetector::trackingGeometry);
   
    py::class_<HGCDetector::Config>(hgcDetector, "Config")
        .def(py::init<>())
        .def_readwrite("fileName", &HGCDetector::Config::fileName)
        .def_readwrite("logLevel", &HGCDetector::Config::logLevel);

    m.def("dumpGeometry", [](std::shared_ptr<const Acts::TrackingGeometry> geometry) {
        std::cout << "=== ACTS GEOMETRY DUMP ===" << std::endl;
        geometry->visitSurfaces([](const Acts::Surface* surface) {
            auto geoId = surface->geometryId();
            auto center = surface->center(Acts::GeometryContext::dangerouslyDefaultConstruct());
            std::cout << "Volume: " << geoId.volume()
                      << " | Layer: " << geoId.layer()
                      << " | Sensitive: " << geoId.sensitive()
                      << " | Center: (" << center.x() << ", " << center.y() << ", " << center.z() << ")"
                      << std::endl;
        });
        std::cout << "==========================" << std::endl;
    });

    // Expose accessor for last heap-stored extracted parameters pointer
    m.def("get_last_extracted_params_addr", []() -> uintptr_t {
        return reinterpret_cast<uintptr_t>(g_lastHeapParams);
    });


}
}
