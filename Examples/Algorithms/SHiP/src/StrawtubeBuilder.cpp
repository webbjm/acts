#include "ActsExamples/SHiP/StrawtubeBuilder.hpp"
#include "ActsExamples/SHiP/HelperFunctions.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Plugins/Root/TGeoMaterialConverter.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "TGeoManager.h"
#include "TGeoTube.h"
#include "TGeoNode.h"
#include <algorithm>
#include "Acts/Geometry/ProtoLayer.hpp"     
#include "Acts/Geometry/Extent.hpp"  
#include "Acts/Utilities/AxisDefinitions.hpp"

namespace ActsExamples {

StrawtubeBuilder::StrawtubeBuilder(const Config& cfg)
 : Detector(Acts::getDefaultLogger("StrawtubeBuilder", cfg.logLevel)), m_cfg(cfg) {

 if (gGeoManager == nullptr && !m_cfg.fileName.empty()) {
     TGeoManager::Import(m_cfg.fileName.c_str());
 }

 Acts::Material strawGas = Acts::Material::fromMassDensity(
     30000.0 * Acts::UnitConstants::mm, // X0
     50000.0 * Acts::UnitConstants::mm, // L0
     30.0, 15.0, 
     0.0018 * Acts::UnitConstants::g / Acts::UnitConstants::cm3);
 
 Acts::Material strawWall = Acts::Material::fromMassDensity(
     300.0 * Acts::UnitConstants::mm,
     400.0 * Acts::UnitConstants::mm,
     27.0, 13.0,
     1.4 * Acts::UnitConstants::g / Acts::UnitConstants::cm3);
 
 Acts::MaterialSlab matGas(strawGas, 10.0 * Acts::UnitConstants::mm); 
 Acts::MaterialSlab matTube(strawWall, 0.0036 * Acts::UnitConstants::cm);
 
 Acts::MaterialSlab totalStraw = Acts::MaterialSlab::combineLayers(matTube, matGas);
 totalStraw = Acts::MaterialSlab::combineLayers(totalStraw, matTube);
 
 auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(totalStraw);

 // Sort straws
 std::vector<StrawNode> allStraws;
 findStraws(gGeoManager->GetTopVolume(), Acts::Transform3::Identity(), allStraws);
 std::sort(allStraws.begin(), allStraws.end(), [](const auto& a, const auto& b) { return a.x < b.x; });

 // Cluster by 1mm window
 std::vector<std::vector<StrawNode>> clusters;
 for (const auto& straw : allStraws) {
     if (clusters.empty() || std::abs(straw.x - clusters.back().front().x) > 1.0 * Acts::UnitConstants::mm) {
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
        double stationX = 0;
        std::vector<std::shared_ptr<const Acts::Surface>> surfaces;

        for (auto& s : cluster) {
            auto* tube = dynamic_cast<const TGeoTube*>(s.node->GetVolume()->GetShape());
            if (!tube) continue;

            auto bounds = std::make_shared<const Acts::LineBounds>(
                tube->GetRmax() * Acts::UnitConstants::cm,
                tube->GetDZ() * Acts::UnitConstants::cm);

            auto detElement = std::make_shared<SHiPDetectorElement>(
                std::make_shared<const Acts::Transform3>(s.transform), bounds, 2.0, surfaceMaterial);

            m_detElements.push_back(detElement);
            surfaces.push_back(detElement->surface().getSharedPtr());
            
            stationX += s.transform.translation().x();
        }
        stationX /= cluster.size();

        std::sort(surfaces.begin(), surfaces.end(), [](const auto& a, const auto& b) {
            return a->center(Acts::GeometryContext()).y() > b->center(Acts::GeometryContext()).y();
        });


        Acts::ProtoLayer pl(context, surfaces);

  
    pl.extent.range(Acts::AxisDirection::AxisX).setMin(stationX - 1.0);
    pl.extent.range(Acts::AxisDirection::AxisX).setMax(stationX + 1.0);
   
    pl.extent.range(Acts::AxisDirection::AxisY).setMin(-5000.0);
    pl.extent.range(Acts::AxisDirection::AxisY).setMax(5000.0);
   
    pl.extent.range(Acts::AxisDirection::AxisZ).setMin(-5000.0);
    pl.extent.range(Acts::AxisDirection::AxisZ).setMax(5000.0);

        Acts::Transform3 layerTransf = Acts::Transform3(Acts::Translation3(stationX, 0, 0));
        auto layer = layerCreator->planeLayer(context, surfaces, 300, 1, Acts::AxisDirection::AxisX, pl, layerTransf);
        for (const auto& surf : surfaces) {
                        const_cast<Acts::Surface*>(surf.get())->associateLayer(*layer);
                    }
        m_layers.push_back(layer);
   }
}


 void StrawtubeBuilder::findStraws(TGeoVolume* vol, const Acts::Transform3& parentTrans, std::vector<StrawNode>& straws) {
     for (int i = 0; i < vol->GetNdaughters(); ++i) {
         TGeoNode* node = vol->GetNode(i);
         const Double_t* r = node->GetMatrix()->GetRotationMatrix();
         const Double_t* t = node->GetMatrix()->GetTranslation();
         Acts::RotationMatrix3 rot;
         rot << r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8];
         Acts::Vector3 trans(t[0] * Acts::UnitConstants::cm, t[1] * Acts::UnitConstants::cm, t[2] * Acts::UnitConstants::cm);

         Acts::Transform3 nodeGlobalTrans = parentTrans * Acts::Transform3(Acts::Translation3(trans) * rot);

         if (std::string(node->GetName()).find("straw") != std::string::npos) {
             Acts::Transform3 recoTrans = ActsExamples::ApplyRecoRotation(nodeGlobalTrans);
             straws.push_back({node, recoTrans, recoTrans.translation().x()});
         } else {
             findStraws(node->GetVolume(), nodeGlobalTrans, straws);
         }
     }
}

} // namespace ActsExamples
