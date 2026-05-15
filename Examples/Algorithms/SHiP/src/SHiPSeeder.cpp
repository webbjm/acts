#include "ActsExamples/SHiP/SHiPSeeder.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include <algorithm>
#include <map>
#include <set>
#include <vector>

namespace ActsExamples {

SHiPSeeder::SHiPSeeder(Config config, Acts::Logging::Level level)
    : IAlgorithm("SHiPSeeder", level), m_cfg(std::move(config)) {
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
}

ProcessCode SHiPSeeder::execute(const AlgorithmContext& ctx) const {
  const auto& measurements = m_inputMeasurements(ctx);

  // Group hit indices by GeometryIdentifier
  std::map<Acts::GeometryIdentifier, std::vector<std::size_t>> layerHitIndices;
  for (const auto& sl : measurements.orderedIndices()) {
    layerHitIndices[sl.geometryId()].push_back(sl.index());
  }

  // Identify unique active layers (already sorted by GeometryIdentifier hierarchy)
  std::vector<Acts::GeometryIdentifier> activeLayers;
  for (auto const& [id, indices] : layerHitIndices) {
    activeLayers.push_back(id);
  }

  // Setup Propagator for confirmation through magnetic field
  using Stepper = Acts::EigenStepper<>;
  using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
  Stepper stepper(m_cfg.magneticField);
  Propagator propagator(stepper, Acts::Navigator({m_cfg.trackingGeometry}));
  using Options = typename Propagator::Options<>;
  Options pOptions(ctx.geoContext, ctx.magFieldContext);

  ProtoTrackContainer protoTracks;
  std::set<Acts::GeometryIdentifier> lockedLayers;

  // Sliding Window: Layers i to i+3 form the seed; Layer i+4 confirms
  for (size_t i = 0; i + 4 < activeLayers.size(); ++i) {
    const auto& id1 = activeLayers[i];
    if (lockedLayers.count(id1)) continue; // Prune layer if already seeded

    const auto& ids = std::vector<Acts::GeometryIdentifier>{
        id1, activeLayers[i + 1], activeLayers[i + 2], activeLayers[i + 3]};
    const auto& idConf = activeLayers[i + 4];

    // Seed hits must be within the same physical Volume (Station)
    if (ids[0].volume() != ids[3].volume()) continue;

    bool foundAnyForLayer = false;
    for (const auto idx1 : layerHitIndices[ids[0]]) {
      for (const auto idx2 : layerHitIndices[ids[1]]) {
        for (const auto idx3 : layerHitIndices[ids[2]]) {
          for (const auto idx4 : layerHitIndices[ids[3]]) {
            
            // 4-Hit 3D Line Solver
            Eigen::Matrix4d M;
            Eigen::Vector4d B;
            std::vector<std::size_t> indices = {idx1, idx2, idx3, idx4};
            auto surfRef = m_cfg.trackingGeometry->findSurface(ids[0]);
            double zRef = surfRef->center(ctx.geoContext).z();

            for (int k = 0; k < 4; ++k) {
              auto surf = m_cfg.trackingGeometry->findSurface(ids[k]);
              double zk = surf->center(ctx.geoContext).z() - zRef;
              Acts::Vector3 uAxis = surf->transform(ctx.geoContext).linear().col(0);
              M.row(k) << uAxis.x(), zk * uAxis.x(), uAxis.y(), zk * uAxis.y();
              B(k) = measurements.getMeasurement(indices[k]).parameters()[Acts::eBoundLoc0];
            }

            Eigen::Vector4d sol = M.colPivHouseholderQr().solve(B);
            Acts::Vector3 dir = Acts::Vector3(sol[1], sol[3], 1.0).normalized();

            // Next-Layer Propagation Confirmation
            auto surf4 = m_cfg.trackingGeometry->findSurface(ids[3]);
            Acts::BoundVector bParams;
            bParams << measurements.getMeasurement(idx4).parameters()[Acts::eBoundLoc0], 0, 
                       std::atan2(dir.y(), dir.x()), std::acos(dir.z()),
                       1.0 / (m_cfg.defaultMomentum * Acts::UnitConstants::GeV), 0;
            Acts::BoundTrackParameters startPar(surf4->getSharedPtr(), bParams, std::nullopt, 
                                                 Acts::ParticleHypothesis::muon());

            auto surfConf = m_cfg.trackingGeometry->findSurface(idConf);
            auto result = propagator.propagate(startPar, *surfConf, pOptions);

            if (result.ok() && result.value().endParameters.has_value()) {
              double predX = result.value().endParameters->parameters()[Acts::eBoundLoc0];
              for (const auto idxConf : layerHitIndices[idConf]) {
                if (std::abs(measurements.getMeasurement(idxConf).parameters()[Acts::eBoundLoc0] - predX) < m_cfg.projectionWindow) {
                  protoTracks.push_back({static_cast<unsigned int>(idx1), 
                                         static_cast<unsigned int>(idx2), 
                                         static_cast<unsigned int>(idx3), 
                                         static_cast<unsigned int>(idx4), 
                                         static_cast<unsigned int>(idxConf)});
                  foundAnyForLayer = true;
                  goto next_combination; 
                }
              }
            }
            next_combination:;
          }
        }
      }
    }
    if (foundAnyForLayer) lockedLayers.insert(id1);
  }

  m_outputProtoTracks(ctx, std::move(protoTracks));
  return ProcessCode::SUCCESS;
}

} // namespace ActsExamples
