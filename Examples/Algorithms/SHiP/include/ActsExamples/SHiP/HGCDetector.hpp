#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/DetectorCommons/Detector.hpp"
#include "ActsExamples/SHiP/SHiPDetectorElement.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <vector>

namespace ActsExamples {

class HGCBuilder; 

class HGCDetector : public Detector {
 public:
  struct Config {
    std::string fileName = "";
    Acts::Logging::Level logLevel = Acts::Logging::WARNING;
  };

  explicit HGCDetector(const Config& cfg);

  /// Return the built tracking geometry
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry() const {
    return m_trackingGeometry;
  }

 private:
  Config m_cfg;
  std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeometry;
  
  /// Keeps detector elements alive for the lifetime of the geometry
  std::vector<std::shared_ptr<const SHiPDetectorElement>> m_detElements;
};

}  // namespace ActsExamples
