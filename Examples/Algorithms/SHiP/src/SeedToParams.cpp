#include "ActsExamples/SHiP/SeedToParams.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include <cmath>

namespace ActsExamples {

SeedToParams::SeedToParams(Config config, Acts::Logging::Level level)
    : IAlgorithm("SeedToParams", level), m_cfg(std::move(config)) {
  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_outputParams.initialize(m_cfg.outputParams);
}

ProcessCode SeedToParams::execute(const AlgorithmContext& ctx) const {
  const auto& protoTracks = m_inputProtoTracks(ctx);
  const auto& measurements = m_inputMeasurements(ctx);

  auto fieldCache = m_cfg.magneticField->makeCache(ctx.magFieldContext);

  std::vector<Acts::BoundTrackParameters> initialParams;
  initialParams.reserve(protoTracks.size());

  for (const auto& indices : protoTracks) {
    if (indices.size() < 4) continue;

    // Use a vector of indices to avoid taking addresses of temporary proxies
    std::vector<std::size_t> idxs = {indices[0], indices[1], indices[2], indices[3]};

    // Setup Matrix Solver
    Eigen::Matrix4d M;
    Eigen::Vector4d B;
    auto surf1 = m_cfg.trackingGeometry->findSurface(measurements.getMeasurement(idxs[0]).geometryId());
    double zRef = surf1->center(ctx.geoContext).z();

    for (int k = 0; k < 4; ++k) {
      auto m = measurements.getMeasurement(idxs[k]);
      auto surf = m_cfg.trackingGeometry->findSurface(m.geometryId());
      double zk = surf->center(ctx.geoContext).z() - zRef;
      Acts::Vector3 uAxis = surf->transform(ctx.geoContext).linear().col(0);
      M.row(k) << uAxis.x(), zk * uAxis.x(), uAxis.y(), zk * uAxis.y();
      B(k) = m.parameters()[Acts::eBoundLoc0];
    }

    Eigen::Vector4d sol = M.colPivHouseholderQr().solve(B);
    Acts::Vector3 dir = Acts::Vector3(sol[1], sol[3], 1.0).normalized();

    // Momentum Estimation
    auto surf4 = m_cfg.trackingGeometry->findSurface(measurements.getMeasurement(idxs[3]).geometryId());
    Acts::Vector3 midPoint(sol[0], sol[2], zRef + (surf4->center(ctx.geoContext).z() - zRef) / 2.0);

    auto fieldResult = m_cfg.magneticField->getField(midPoint, fieldCache);
    double B_val = (fieldResult.ok()) ? fieldResult.value().norm() : 0.0;

    double qOverP = 1.0 / (m_cfg.defaultMomentum * Acts::UnitConstants::GeV);
    if (B_val > 0.01 * Acts::UnitConstants::T) {
      double dz = surf4->center(ctx.geoContext).z() - zRef;
      double x_pred = sol[0] + sol[1] * dz;
      double sagitta = std::abs(measurements.getMeasurement(idxs[3]).parameters()[Acts::eBoundLoc0] - x_pred);
      double p_est = (0.3 * B_val * std::pow(dz, 2)) / (8.0 * std::max(sagitta, 1e-4));
      qOverP = 1.0 / std::clamp(p_est, 0.5 * Acts::UnitConstants::GeV, 100.0 * Acts::UnitConstants::GeV);
    }

    // Construct Parameters
    Acts::BoundVector params;
    params << measurements.getMeasurement(idxs[3]).parameters()[Acts::eBoundLoc0], 0.0,
              std::atan2(dir.y(), dir.x()), std::acos(dir.z()), qOverP, 0.0;

    Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Identity();
    cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = 0.05 * 0.05;
    cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = 0.05 * 0.05;
    cov(Acts::eBoundPhi, Acts::eBoundPhi) = 0.1;
    cov(Acts::eBoundTheta, Acts::eBoundTheta) = 0.1;
    cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = 1e-2;

    initialParams.emplace_back(surf4->getSharedPtr(), params,
                               std::optional<Acts::BoundSquareMatrix>(cov),
                               Acts::ParticleHypothesis::muon());
  }

  m_outputParams(ctx, std::move(initialParams));
  return ProcessCode::SUCCESS;
}

} // namespace ActsExamples
