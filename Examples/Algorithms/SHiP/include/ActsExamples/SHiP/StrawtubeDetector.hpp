#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/DetectorCommons/Detector.hpp"
#include "ActsExamples/SHiP/SHiPDetectorElement.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>
#include <vector>

namespace ActsExamples {

class StrawtubeBuilder; 

class StrawtubeDetector : public Detector {
 public:
  struct Config {
    std::string fileName = "";
    Acts::Logging::Level logLevel = Acts::Logging::WARNING;
  };

  using DetectorElementPtr = std::shared_ptr<const SHiPDetectorElement>;

  explicit StrawtubeDetector(const Config& cfg);

  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry() const override {
      return m_trackingGeometry;
  }

  const std::vector<DetectorElementPtr>& detectorElements() const {
      return m_detElements;
  }

 private:
  Config m_cfg;
  std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeometry;
  std::vector<DetectorElementPtr> m_detElements;
};

}  // namespace ActsExamples
