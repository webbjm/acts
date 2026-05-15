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
 * @class SHiPSeeder
 * @brief Idiomatic ACTS seeder for SHiP Silicon (X/Y) and SciFi (U/V) detectors.
 *
 * Implements a geometry-aware sliding window that:
 * 1. Groups hits by GeometryIdentifier (Volume/Layer).
 * 2. Fits a 3D line to 4-hit candidates within a single physical volume.
 * 3. Confirms candidates by propagating to the very next active layer.
 */
class SHiPSeeder : public IAlgorithm {
 public:
  struct Config {
    /// Input measurements collection
    std::string inputMeasurements;
    /// Output proto tracks collection
    std::string outputProtoTracks;
    /// Tracking geometry
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Magnetic field
    std::shared_ptr<Acts::MagneticFieldProvider> magneticField;

    /// Seeding window for matching projected hits in the confirmation layer (mm)
    double projectionWindow = 5.0;
    /// Default momentum for initial propagation guess (GeV)
    double defaultMomentum = 3.0;
  };

  SHiPSeeder(Config config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const override;

 private:
  Config m_cfg;

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this, "InputMeasurements"};
  WriteDataHandle<ProtoTrackContainer> m_outputProtoTracks{this, "OutputProtoTracks"};
};

}  // namespace ActsExamples
