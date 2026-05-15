#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <memory>
#include <string>
#include <vector>

namespace ActsExamples {

/**
 * @class SeedToParams
 * @brief Converts ProtoTracks into BoundTrackParameters.
 *
 * Uses a 4-hit 3D line solver to derive track direction. Dynamically estimates
 * q/p based on local magnetic field strength (B-field sagitta) or defaults
 * to a high-momentum straight-line assumption in field-free silicon regions.
 */
class SeedToParams : public IAlgorithm {
 public:
  struct Config {
    std::string inputProtoTracks;
    std::string inputMeasurements;
    std::string outputParams;
    
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    std::shared_ptr<Acts::MagneticFieldProvider> magneticField;

    double defaultMomentum = 3.0; // GeV
  };

  SeedToParams(Config config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const override;

 private:
  Config m_cfg;

  ReadDataHandle<ProtoTrackContainer> m_inputProtoTracks{this, "InputProtoTracks"};
  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this, "InputMeasurements"};
  WriteDataHandle<std::vector<Acts::BoundTrackParameters>> m_outputParams{this, "OutputParams"};
};

}  // namespace ActsExamples
