#include "ActsExamples/SHiP/StrawtubeDetector.hpp"
#include "ActsExamples/SHiP/StrawtubeBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Definitions/Units.hpp"

#include <limits>

namespace ActsExamples {

StrawtubeDetector::StrawtubeDetector(const Config& cfg)
    : Detector(Acts::getDefaultLogger("StrawtubeDetector", cfg.logLevel)), m_cfg(cfg) {
    
    // Build layers using the builder
    StrawtubeBuilder::Config builderCfg;
    builderCfg.fileName = m_cfg.fileName;
    builderCfg.logLevel = m_cfg.logLevel;
    
    StrawtubeBuilder builder(builderCfg);
    const auto& layers = builder.layers();
    m_detElements = builder.detectorElements();

    // Determine bounds dynamically (Use a small buffer)
    double minX = std::numeric_limits<double>::max();
    double maxX = -std::numeric_limits<double>::max();
    for (const auto& layer : layers) {
        double x = layer->surfaceArray()->transform().translation().x();
        minX = std::min(minX, x);
        maxX = std::max(maxX, x);
    }

    // Create LayerArray
    Acts::LayerArrayCreator::Config lacConfig;
    Acts::LayerArrayCreator layArrCreator(lacConfig, 
        Acts::getDefaultLogger("LayerCreator", Acts::Logging::VERBOSE));
    
    auto layArr = layArrCreator.layerArray(Acts::GeometryContext(), layers,
                                           minX - 50.0, maxX + 50.0,
                                           Acts::BinningType::arbitrary,
                                           Acts::AxisDirection::AxisX);

    // Volume Definition: Centered on the stations to keep Navigation fast
    double volCenter = (minX + maxX) / 2.0;
    double volHalfX = (maxX - minX) / 2.0 + 1000.0; // 1m buffer
    
    auto bounds = std::make_shared<Acts::CuboidVolumeBounds>(volHalfX, 5000.0, 5000.0);
    Acts::Transform3 volumeTransform(Acts::Translation3(volCenter, 0.0, 0.0));

    // Build and Index Volume
    auto trackVolume = std::make_shared<Acts::TrackingVolume>(
        volumeTransform, bounds, nullptr, std::move(layArr), 
        nullptr, Acts::MutableTrackingVolumeVector{}, "StrawVolume");

    trackVolume->assignGeometryId(Acts::GeometryIdentifier().withVolume(1));

    // Build final geometry
    m_trackingGeometry = std::make_shared<Acts::TrackingGeometry>(trackVolume);
    
    ACTS_INFO("StrawtubeDetector built with " << layers.size() << " layers.");
}
}
