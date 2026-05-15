 #pragma once


// ACTS Core Definitions
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

// Framework & Event Data
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include <vector>
#include <string>
#include <memory>
namespace ActsExamples {

class SHiPMeasurementProvider : public IAlgorithm {
public:
    struct Config {
        std::vector<std::vector<float>>* inputHitArray = nullptr;
        std::string outputMeasurements = "measurements";
        std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

        // Detector-specific configuration (mm)
        double strawRes = 0.12;
        double scifiPitch = 1.0;
        double siliconPitch = 0.075;
        double scifiRes = 0.05;
        double siliconRes = 0.015;
        double minBound = -125.0; 
        double noiseFloor = 0;//0.01;
    };

    SHiPMeasurementProvider(Config config, Acts::Logging::Level level);
    ProcessCode execute(const AlgorithmContext& ctx) const override;
    const Config& config() const { return m_cfg; }

private:
    Config m_cfg;
    Acts::GeometryIdentifier getGeoId(const std::vector<float>& iHit) const;
    WriteDataHandle<MeasurementContainer> m_outputMeasurements{this, "outputMeasurements"};
};

} // namespace ActsExamples
