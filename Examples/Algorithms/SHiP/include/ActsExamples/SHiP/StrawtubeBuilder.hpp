#pragma once

#include "Acts/Geometry/Layer.hpp" 
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DetectorCommons/Detector.hpp"
#include "ActsExamples/SHiP/SHiPDetectorElement.hpp"

#include <memory>
#include <string>
#include <vector>

// Forward declarations
class TGeoVolume;
class TGeoNode;

namespace ActsExamples {

struct StrawNode {
    TGeoNode* node;
    Acts::Transform3 transform;
    double x;
};

class StrawtubeBuilder : public Detector {
 public:
  struct Config {
    Acts::Logging::Level logLevel = Acts::Logging::WARNING;
    Acts::Logging::Level layerLogLevel = Acts::Logging::WARNING;
    std::string fileName = "";
  };

  explicit StrawtubeBuilder(const Config& cfg);

  const std::vector<Acts::LayerPtr>& layers() const { return m_layers; }
  
  using DetectorElementPtr = std::shared_ptr<const SHiPDetectorElement>;
  const std::vector<DetectorElementPtr>& detectorElements() const {
      return m_detElements;
  }

  void findStraws(TGeoVolume* vol, const Acts::Transform3& parentTrans, std::vector<StrawNode>& straws);

 private:
  Config m_cfg;
  std::vector<Acts::LayerPtr> m_layers;
  std::vector<DetectorElementPtr> m_detElements;
};

} // namespace ActsExamples
