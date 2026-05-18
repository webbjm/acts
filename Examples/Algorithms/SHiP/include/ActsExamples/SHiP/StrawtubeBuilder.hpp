#pragma once

#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DetectorCommons/Detector.hpp"

#include <memory>
#include <string>
#include <vector>

namespace ActsExamples {

class StrawtubeBuilder : public Detector {
 public:
  struct Config {
    Acts::Logging::Level logLevel = Acts::Logging::WARNING;
    Acts::Logging::Level layerLogLevel = Acts::Logging::WARNING;
    std::string fileName = "";
  };
  
  explicit StrawtubeBuilder(const Config& cfg);
  const std::vector<Acts::LayerPtr>& layers() const { return m_layers; }

 private:
  Config m_cfg;
  std::vector<Acts::LayerPtr> m_layers;
};

} // namespace ActsExamples

