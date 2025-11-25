#pragma once

#include "Acts/Plugins/Root/TGeoLayerBuilder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DetectorCommons/Detector.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <cstddef>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ActsExamples {

class MTCBuilder : public Detector  {
 public:
  struct Config {
    Acts::Logging::Level logLevel = Acts::Logging::WARNING;
    Acts::Logging::Level surfaceLogLevel = Acts::Logging::WARNING;
    Acts::Logging::Level layerLogLevel = Acts::Logging::WARNING;
    Acts::Logging::Level volumeLogLevel = Acts::Logging::WARNING;
    std::string fileName;
  };
  explicit MTCBuilder(const Config& cfg);

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
