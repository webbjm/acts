#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DetectorCommons/Detector.hpp"
#include "ActsExamples/SHiP/StrawtubeBuilder.hpp"

#include <memory>
#include <string>

namespace ActsExamples {

class StrawtubeDetector : public Detector {
 public:
  struct Config {
    std::string fileName = "";
    Acts::Logging::Level logLevel = Acts::Logging::WARNING;
  };
  
  explicit StrawtubeDetector(const Config& cfg);

  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry() const override {
      return m_trackingGeometry;
  }


 private:
  Config m_cfg;
  std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeometry;
};

}  // namespace ActsExamples
