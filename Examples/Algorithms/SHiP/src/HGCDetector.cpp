#include "ActsExamples/SHiP/HGCDetector.hpp"
#include "ActsExamples/SHiP/HGCBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Definitions/Units.hpp"
#include <limits>
#include <algorithm>

namespace ActsExamples {

HGCDetector::HGCDetector(const Config& cfg)
    : Detector(Acts::getDefaultLogger("HGCDetector", cfg.logLevel)), m_cfg(cfg) {

    // 1. Build the layers and sensitive elements using the HGCBuilder
    HGCBuilder::Config builderCfg;
    builderCfg.fileName = m_cfg.fileName;
    builderCfg.logLevel = m_cfg.logLevel;
    HGCBuilder builder(builderCfg);

    const auto& layers = builder.layers();
    m_detElements = builder.detectorElements();
    
    // 2. Determine geometric range dynamically
    double minX = std::numeric_limits<double>::max();
    double maxX = -std::numeric_limits<double>::max();

    Acts::GeometryContext context;
    for (const auto& layer : layers) {
        double x = layer->surfaceRepresentation().center(context).x();
        minX = std::min(minX, x);
        maxX = std::max(maxX, x);
    }

    // 3. Create LayerArray
    Acts::LayerArrayCreator::Config lacConfig;
    Acts::LayerArrayCreator layArrCreator(lacConfig, 
        Acts::getDefaultLogger("LayerCreator", Acts::Logging::VERBOSE));
    
    auto layArr = layArrCreator.layerArray(context, layers, 
                                           minX - 100.0 * Acts::UnitConstants::mm, 
                                           maxX + 100.0 * Acts::UnitConstants::mm, 
                                           Acts::BinningType::arbitrary, 
                                           Acts::AxisDirection::AxisX);

    // 4. Build TrackingVolume
    auto bounds = std::make_shared<Acts::CuboidVolumeBounds>(
        100.0 * Acts::UnitConstants::m,
        100.0 * Acts::UnitConstants::m,
        100.0 * Acts::UnitConstants::m);
        
    auto trackVolume = std::make_shared<Acts::TrackingVolume>(
        Acts::Transform3::Identity(),
        bounds, 
        nullptr, 
        std::move(layArr), 
        nullptr, 
        Acts::MutableTrackingVolumeVector{}, 
        "HGCVolume");
    
    // 5. Construct the final TrackingGeometry
    m_trackingGeometry = std::make_shared<Acts::TrackingGeometry>(trackVolume);
    
    ACTS_INFO("HGCDetector built with " << layers.size() << " active/passive layers.");
}

} // namespace ActsExamples
