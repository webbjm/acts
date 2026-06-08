#include "ActsExamples/SHiP/HGCBuilder.hpp"
#include "ActsExamples/SHiP/HelperFunctions.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "TGeoNode.h"
#include <algorithm>                
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Geometry/Extent.hpp"
#include <iostream>
namespace ActsExamples {

using namespace Acts::UnitLiterals;

HGCBuilder::HGCBuilder(const Config& cfg)
 : Detector(Acts::getDefaultLogger("HGCBuilder", cfg.logLevel)), m_cfg(cfg) {
    
 if (gGeoManager == nullptr && !m_cfg.fileName.empty()) {
     TGeoManager::Import(m_cfg.fileName.c_str());
 }

 // 1. Setup Materials (Values from original user script)
 Acts::Material silicon = Acts::Material::fromMassDensity(9.37_cm, 46.52_cm, 28.08, 14, 2.329_g / 1_cm3);
 auto siMat = std::make_shared<Acts::HomogeneousSurfaceMaterial>(Acts::MaterialSlab(silicon, 300_um));
    
 Acts::Material tungsten = Acts::Material::fromMassDensity(0.369_cm, 10.35_cm, 177.65, 71.72, 18.78_g / 1_cm3);
 auto wMat = std::make_shared<Acts::HomogeneousSurfaceMaterial>(Acts::MaterialSlab(tungsten, 3.5_mm));

 // 2. Recursive extraction and sorting
 std::vector<HGCNode> sensorNodes, absorberNodes;
 findHGCModules(gGeoManager->GetTopVolume(), Acts::Transform3::Identity(), sensorNodes, absorberNodes);
 
 auto sortByX = [](const auto& a, const auto& b) { return a.x < b.x; };
 std::sort(sensorNodes.begin(), sensorNodes.end(), sortByX);
 std::sort(absorberNodes.begin(), absorberNodes.end(), sortByX);

 Acts::LayerCreator::Config lcConfig;
 lcConfig.surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>();
 auto layerCreator = std::make_shared<const Acts::LayerCreator>(lcConfig, Acts::getDefaultLogger("LayerCreator", m_cfg.layerLogLevel));

 Acts::GeometryContext context;

 // 3. Build Active Layers (8 sensors per layer)
 for (size_t i = 0; i < sensorNodes.size(); i += 8) {
     std::vector<std::shared_ptr<const Acts::Surface>> surfaces;
     double avgX = 0;
     for (size_t j = 0; j < 8 && (i + j) < sensorNodes.size(); ++j) {
         auto& s = sensorNodes[i + j];
         auto* bbox = dynamic_cast<const TGeoBBox*>(s.node->GetVolume()->GetShape());
         auto bounds = std::make_shared<const Acts::RectangleBounds>(bbox->GetDX() * 10, bbox->GetDY() * 10);
         
         auto detElement = std::make_shared<SHiPDetectorElement>(
             std::make_shared<const Acts::Transform3>(s.transform), bounds, 300_um, siMat);

         m_detElements.push_back(detElement);
         surfaces.push_back(detElement->surface().getSharedPtr());
         avgX += s.x;
     }
     avgX /= 8.0;

     Acts::ProtoLayer pl(context, surfaces);
     pl.extent.range(Acts::AxisDirection::AxisX).setMin(avgX - 1.0);
     pl.extent.range(Acts::AxisDirection::AxisX).setMax(avgX + 1.0);
     pl.extent.range(Acts::AxisDirection::AxisY).setMin(-5000.0);
     pl.extent.range(Acts::AxisDirection::AxisY).setMax(5000.0);
     pl.extent.range(Acts::AxisDirection::AxisZ).setMin(-5000.0);
     pl.extent.range(Acts::AxisDirection::AxisZ).setMax(5000.0);

     auto layer = layerCreator->planeLayer(context, surfaces, 0, 0, Acts::AxisDirection::AxisX, pl, Acts::Transform3(Acts::Translation3(avgX, 0, 0)));
     for (const auto& surf : surfaces) const_cast<Acts::Surface*>(surf.get())->associateLayer(*layer);
     m_layers.push_back(layer);
 }

 // 4. Build Passive Tungsten Layers
 for (auto& anode : absorberNodes) {
     auto* bbox = dynamic_cast<const TGeoBBox*>(anode.node->GetVolume()->GetShape());
     auto bounds = std::make_shared<const Acts::RectangleBounds>(bbox->GetDX() * 10, bbox->GetDY() * 10);
     
     auto pLayer = Acts::PlaneLayer::create(anode.transform, bounds, nullptr, 3.5_mm);
     pLayer->surfaceRepresentation().assignSurfaceMaterial(wMat);
     m_layers.push_back(pLayer);
 }

 // 5. Final global sort for Navigator stability
 std::sort(m_layers.begin(), m_layers.end(), [](const auto& a, const auto& b) {
     return a->surfaceRepresentation().center(Acts::GeometryContext()).x() < 
            b->surfaceRepresentation().center(Acts::GeometryContext()).x();
 });
}


void HGCBuilder::findHGCModules(TGeoVolume* vol, const Acts::Transform3& parentTrans,
                                std::vector<HGCNode>& sensors,
                                std::vector<HGCNode>& absorbers) {

    // We only need the top-level volume to start the search
    // Using a TGeoIterator is safer and cleaner than manual recursion
    TGeoIterator it(vol);
    TGeoNode* node;

    while ((node = it.Next())) {
        std::string name = node->GetName();
        if (name.find("SensorVolume") == std::string::npos &&
            name.find("TargetVol") == std::string::npos) {
            continue;
        }

        // 1. Get the absolute global matrix for this node
        // GetCurrentMatrix() in an iterator returns the concatenated global transform
//        const TGeoHMatrix* globalMatrix = it.GetCurrentMatrix();
        const TGeoHMatrix* globalMatrix = static_cast<const TGeoHMatrix*>(it.GetCurrentMatrix());
        const Double_t* t = globalMatrix->GetTranslation();
        const Double_t* r = globalMatrix->GetRotationMatrix();

        // 2. Map ROOT(Z, Y, X) to ACTS(X, Y, Z)
        // Adjust these mappings based on your specific SHiP -> ACTS coordinate requirement
        Acts::Vector3 trans(t[2] * 10.0, t[1] * 10.0, -t[0] * 10.0);

        Acts::RotationMatrix3 rot;
        // Map the rotation correctly
        rot << r[0], r[1], r[2],
               r[3], r[4], r[5],
               r[6], r[7], r[8];

        // 3. Create the final Transform
        Acts::Transform3 trafo = Acts::Transform3::Identity();
        trafo.linear() = rot; // Add fixed rotation here if SHiP/ACTS axes don't match
        trafo.translation() = trans;

        if (name.find("SensorVolume") != std::string::npos) {
            sensors.push_back({node, trafo, trans.x()});
        } else {
            absorbers.push_back({node, trafo, trans.x()});
        }
    }
}


/*
void HGCBuilder::findHGCModules(TGeoVolume* vol, const Acts::Transform3& parentTrans, 
                               std::vector<HGCNode>& sensors, 
                               std::vector<HGCNode>& absorbers) {
    for (int i = 0; i < vol->GetNdaughters(); ++i) {
        TGeoNode* node = vol->GetNode(i);
        const Double_t* r = node->GetMatrix()->GetRotationMatrix();
        const Double_t* t = node->GetMatrix()->GetTranslation();
        Acts::RotationMatrix3 rot;
        rot << r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8];
        Acts::Vector3 trans(t[0] * 10, t[1] * 10, t[2] * 10);
        Acts::Transform3 nodeGlobalTrans = parentTrans * Acts::Transform3(Acts::Translation3(trans) * rot);

        std::string name = node->GetName();
        if (name.find("SensorVolume") != std::string::npos) {
            Acts::Transform3 recoTrans = ActsExamples::ApplyRecoRotation(nodeGlobalTrans);
                      std::cout << "DEBUG: Node: " << name
                           << " | Original Translation: " << trans.transpose()
                           << " | Final recoTrans Translation: " << recoTrans.translation().transpose()
                           << std::endl;
            sensors.push_back({node, recoTrans, recoTrans.translation().x()});
        } else if (name.find("TargetVol") != std::string::npos) {
            Acts::Transform3 recoTrans = ActsExamples::ApplyRecoRotation(nodeGlobalTrans);
            absorbers.push_back({node, recoTrans, recoTrans.translation().x()});
        } else {
            findHGCModules(node->GetVolume(), nodeGlobalTrans, sensors, absorbers);
        }
    }
}

*/
} // namespace ActsExamples

