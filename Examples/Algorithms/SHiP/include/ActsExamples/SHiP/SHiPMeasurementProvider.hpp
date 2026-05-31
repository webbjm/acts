#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include <vector>
#include <string>
#include <memory>

namespace ActsExamples {

class SHiPMeasurementProvider {
public:
    struct Config {
        // Detector-specific configuration (mm)
        double strawRes = 0.12;
        double strawRadius = 10.0;
        double scifiPitch = 1.0;
        double siliconPitch = 0.075;
        double scifiRes = 0.05;
        double siliconRes = 0.015;
        double minBound = -125.0;
        double noiseFloor = 0.0;

        std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    };

    SHiPMeasurementProvider(Config config);

    // Standalone process method
    MeasurementContainer process(const std::vector<std::vector<float>>& inputHits, 
                                 const Acts::GeometryContext& geoCtx) const;

    const Config& config() const { return m_cfg; }

private:
    Config m_cfg;
    Acts::GeometryIdentifier getGeoId(const std::vector<float>& iHit) const;
};

} // namespace ActsExamples
