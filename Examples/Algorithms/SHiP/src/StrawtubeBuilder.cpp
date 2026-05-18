#include "ActsExamples/SHiP/StrawtubeBuilder.hpp"
#include "ActsExamples/SHiP/HelperFunctions.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Plugins/Root/TGeoSurfaceConverter.hpp"
#include "Acts/Plugins/Root/TGeoPrimitivesHelper.hpp"
#include "Acts/Geometry/ProtoLayerHelper.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Definitions/Units.hpp"
#include "TGeoManager.h"
#include "TGeoTube.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace ActsExamples {

struct StrawNode {
    TGeoNode* node;
    Acts::Transform3 transform;
    double z;
};

void findStraws(TGeoVolume* vol, const Acts::Transform3& parentTrans, std::vector<StrawNode>& straws) {
    for (int i = 0; i < vol->GetNdaughters(); ++i) {
        TGeoNode* node = vol->GetNode(i);
        const Double_t* r = node->GetMatrix()->GetRotationMatrix();
        const Double_t* t = node->GetMatrix()->GetTranslation();
        Acts::RotationMatrix3 rot;
        rot << r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8];
        Acts::Vector3 trans(t[0] * Acts::UnitConstants::cm, t[1] * Acts::UnitConstants::cm, t[2] * Acts::UnitConstants::cm);
        Acts::Transform3 nodeTrans = parentTrans * Acts::Transform3(Acts::Translation3(trans) * rot);
        
        Acts::Transform3 globalTrans = parentTrans * nodeTrans;
   
        if (std::string(node->GetName()).find("straw") != std::string::npos) {
            // Apply the reconstruction rotation helper here
            Acts::Transform3 recoTrans = ActsExamples::ApplyRecoRotation(globalTrans);
            straws.push_back({node, recoTrans, recoTrans.translation().z()});
        } else {
            findStraws(node->GetVolume(), globalTrans, straws);
        }

    }
}

StrawtubeBuilder::StrawtubeBuilder(const Config& cfg)
    : Detector(Acts::getDefaultLogger("StrawtubeBuilder", cfg.logLevel)), m_cfg(cfg) {
    
    if (gGeoManager == nullptr && !m_cfg.fileName.empty()) TGeoManager::Import(m_cfg.fileName.c_str());

    //  Define materials for Gas and Straw Wall
    Acts::Material strawTube = Acts::Material::fromMassDensity(
        28.5548 * Acts::UnitConstants::cm, 56.387 * Acts::UnitConstants::cm,
        12.8772, 6.45628, 6.45628 * Acts::UnitConstants::g / Acts::UnitConstants::cm3);
  
    Acts::Material strawGas = Acts::Material::fromMassDensity(
        6119.96 * Acts::UnitConstants::cm, 34106.3 * Acts::UnitConstants::cm,
        37.194, 16.84, 0.00336 * Acts::UnitConstants::g / Acts::UnitConstants::cm3);
  
    //  Combine layers to get effective slab (gas + wall)
    // Note: Assuming a dual straw wall thickness (in and out)
    Acts::MaterialSlab matGas(strawGas, 0.9964 * Acts::UnitConstants::cm);
    Acts::MaterialSlab matTube(strawTube, 0.0018 * Acts::UnitConstants::cm);
  
    Acts::MaterialSlab totalStraw = Acts::MaterialSlab::combineLayers(matTube, matGas);
    totalStraw = Acts::MaterialSlab::combineLayers(totalStraw, matTube);
   
    auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(totalStraw);

    std::vector<StrawNode> allStraws;
    findStraws(gGeoManager->GetTopVolume(), Acts::Transform3::Identity(), allStraws);
    std::sort(allStraws.begin(), allStraws.end(), [](auto& a, auto& b) { return a.z < b.z; });

    std::vector<std::vector<StrawNode>> clusters;
    for (const auto& straw : allStraws) {
        if (clusters.empty() || std::abs(straw.z - clusters.back().front().z) > 1.0 * Acts::UnitConstants::cm) {
            clusters.push_back({straw});
        } else {
            clusters.back().push_back(straw);
        }
    }

    Acts::LayerCreator::Config lcConfig;
    lcConfig.surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>();
    auto layerCreator = std::make_shared<const Acts::LayerCreator>(lcConfig, Acts::getDefaultLogger("LayerCreator", m_cfg.layerLogLevel));
    
    Acts::GeometryContext context;
    for (auto& cluster : clusters) {
       std::vector<std::shared_ptr<const Acts::Surface>> surfaces;
       for (auto& s : cluster) {
           auto* tube = dynamic_cast<const TGeoTube*>(s.node->GetVolume()->GetShape());
   
           auto surface = Acts::Surface::makeShared<Acts::StrawSurface>(
               s.transform,
               tube->GetRmax() * Acts::UnitConstants::cm,
               tube->GetDZ() * Acts::UnitConstants::cm);
   
           surface->assignSurfaceMaterial(surfaceMaterial);
           surfaces.push_back(std::move(surface));
       }
       Acts::ProtoLayer pl(context, surfaces);
       m_layers.push_back(layerCreator->planeLayer(context, surfaces, 0, 0, Acts::AxisDirection::AxisX, pl));
   }
}
} // namespace ActsExamples
