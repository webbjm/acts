#include "ActsExamples/SHiP/StrawtubeDetector.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Definitions/Units.hpp"

#include <limits>

namespace ActsExamples {

StrawtubeDetector::StrawtubeDetector(const Config& cfg)
    : Detector(Acts::getDefaultLogger("StrawtubeDetector", cfg.logLevel)), m_cfg(cfg) {
     // Build the layers using the StrawtubeBuilder
     StrawtubeBuilder::Config builderCfg;
     builderCfg.fileName = m_cfg.fileName;
     builderCfg.logLevel = m_cfg.logLevel;
     StrawtubeBuilder builder(builderCfg);
     const auto& layers = builder.layers();

     // Determine geometric Z-range dynamically for the "Directory" (LayerArray)
     double minZ = std::numeric_limits<double>::max();
     double maxZ = -std::numeric_limits<double>::max();
   
     Acts::GeometryContext context; // Needed for the surface lookup if alignment exists
     for (const auto& layer : layers) {
         double z = layer->surfaceArray()->transform().translation().z();
         minZ = std::min(minZ, z);
         maxZ = std::max(maxZ, z);
     }

     // Create LayerArray (Tight binning for fast lookup)
     Acts::LayerArrayCreator::Config lacConfig;
     Acts::LayerArrayCreator layArrCreator(lacConfig,
         Acts::getDefaultLogger("LayerCreator", Acts::Logging::VERBOSE));

     auto layArr = layArrCreator.layerArray(context, layers,
                                            minZ - 50.0 * Acts::UnitConstants::mm,
                                            maxZ + 50.0 * Acts::UnitConstants::mm,
                                            Acts::BinningType::arbitrary,
                                            Acts::AxisDirection::AxisX);

     // Define Volume Bounds: From Target (Z=0) to past final layer (Z<100m)
     // This allows tracking particles from the target and through the trackers.
     double volStartZ = 0.0 * Acts::UnitConstants::m;
     double volEndZ = 100.0 * Acts::UnitConstants::m;
     double volHalfZ = (volEndZ - volStartZ) / 2.0;
     double volCenterZ = (volEndZ + volStartZ) / 2.0;

     auto bounds = std::make_shared<Acts::CuboidVolumeBounds>(
         500.0 * Acts::UnitConstants::cm, // X half-width
         500.0 * Acts::UnitConstants::cm, // Y half-width
         volHalfZ);                       // Z half-length

     // Build TrackingVolume
     auto trackVolume = std::make_shared<Acts::TrackingVolume>(
         Acts::Transform3(Acts::Translation3(0, 0, volCenterZ)),
         bounds,
         nullptr,
         std::move(layArr),
         nullptr,
         Acts::MutableTrackingVolumeVector{},
         "StrawVolume");

     m_trackingGeometry = std::make_shared<Acts::TrackingGeometry>(trackVolume);
}
} // namespace ActsExamples
