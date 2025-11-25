#include "ActsExamples/TGeoDetector/MTCBuilder.hpp"
#include "ActsExamples/TGeoDetector/ShipDetectorElement.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Root/TGeoParser.hpp"
#include "Acts/Plugins/Root/TGeoPrimitivesHelper.hpp"
#include "Acts/Plugins/Root/TGeoSurfaceConverter.hpp"
#include "Acts/Plugins/Root/TGeoDetectorElement.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/ITrackingVolumeHelper.hpp"
#include "Acts/Geometry/ProtoLayerHelper.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Plugins/Root/TGeoMaterialConverter.hpp"
#include "Acts/Plugins/Root/TGeoSurfaceConverter.hpp"


#include "TGeoShape.h"
#include "TGeoTube.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "TGeoManager.h"
#include "TGeoMatrix.h"


namespace ActsExamples {


namespace {    
/// @brief struct to load the global geometry
  std::shared_ptr<const Acts::TrackingGeometry> buildMTCDetector(
      const MTCBuilder::Config& config, const Acts::GeometryContext& context,
    std::vector<std::shared_ptr<const Acts::DetectorElementBase>>&
        detElementStore,
      const Acts::Logger& logger) {

  using namespace Acts::UnitLiterals;
  TGeoManager::Import(config.fileName.c_str());

  if (gGeoManager != nullptr) {
    std::string volumeName = "*";//"*";
    Acts::TGeoParser::Options tgpOptions;
    tgpOptions.volumeNames = {volumeName};
    tgpOptions.targetNames = {"FiberVol*"};


    Acts::TGeoParser::Options tgpOptions_passive;
    tgpOptions_passive.volumeNames = {volumeName};
    tgpOptions_passive.targetNames = {"MTC_iron*"};

    Acts::TGeoParser::State tgpState_passive;
    tgpState_passive.volume = gGeoManager->GetTopVolume();
    Acts::TGeoParser::select(tgpState_passive, tgpOptions_passive);



    Acts::TGeoParser::State tgpState;
    tgpState.volume = gGeoManager->GetTopVolume();
    Acts::TGeoParser::select(tgpState, tgpOptions);//, gMatrix);

    Acts::ProtoLayerHelper::Config plhConfig;
    Acts::ProtoLayerHelper plHelper(
      plhConfig, Acts::getDefaultLogger("ProtoLayerHelper", Acts::Logging::VERBOSE));

    Acts::LayerCreator::Config lcConfig;
    lcConfig.surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>();
    auto layerCreator = std::make_shared<const Acts::LayerCreator>(
        lcConfig, logger.clone("LayerCreator", config.layerLogLevel));


    std::vector<double> positions;
    std::vector<double> positionsW;

    Acts::Material silicon = Acts::Material::fromMassDensity(
      9.370_cm, 46.52_cm, 28.0855, 14, 2.329_g / 1_cm3);
    Acts::MaterialSlab matProp(silicon, 225._um);
    auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);

    Acts::Material iron = Acts::Material::fromMassDensity(
      1.757_cm, 16.77_cm, 55.845, 26, 7.874_g / 1_cm3);
    Acts::MaterialSlab matPropW(iron, 5._mm);
    auto passiveMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matPropW);



    std::vector<Acts::LayerPtr> layers;
    std::vector<Acts::LayerPtr> passiveLayers;
    int n = 0;
      using LayerSurfaceVector = std::vector<std::shared_ptr<const Acts::Surface>>;
      LayerSurfaceVector layerSurfaces;
      Acts::SurfaceArrayCreator sac;


    for (auto& snode : tgpState.selectedNodes) {
      const auto& shape = (snode.node->GetVolume()->GetShape());
      const auto& transform = *snode.transform;
      n++;

      Acts::RotationMatrix3 rotations = Acts::RotationMatrix3::Identity();
      const Double_t* rotation = transform.GetRotationMatrix();
      const Double_t* translation = transform.GetTranslation();
      Acts::Translation3 trans(translation[0]*10, translation[1]*10, translation[2]*10);
      positions.push_back(translation[2]*10); //Vector of z positions //Changed from 2 -> 0 index
      rotations.col(0) = Acts::Vector3(rotation[0],rotation[3],rotation[6]);
      rotations.col(1) = Acts::Vector3(rotation[1],rotation[4],rotation[7]);
      rotations.col(2) = Acts::Vector3(rotation[2],rotation[5],rotation[8]);
      //Rotate about y axis by ppi/2 and then x axis by -pi/2, combined matrix below
      Acts::RotationMatrix3 rotateFrame;
      rotateFrame.col(0) = Acts::Vector3(0,0,-1);
      rotateFrame.col(1) = Acts::Vector3(0,1,0);
      rotateFrame.col(2) = Acts::Vector3(1,0,0);

      auto* tube = dynamic_cast<const TGeoTube*>(snode.node->GetVolume()->GetShape());

      Acts::RotationMatrix3 combinedRot =  rotateFrame * rotations;
      //Acts::Transform3 trafo(combinedRot * trans);
      Acts::Vector3 t(1_cm*(translation[2]+3.5*tube->GetRmax()),1_cm*translation[1],-1_cm*translation[0]);
      auto trafo = Acts::TGeoPrimitivesHelper::makeTransform(combinedRot.col(0),combinedRot.col(1),combinedRot.col(2),t);

      Acts::Translation3 lTrans(0,0,(translation[2]-3.5*tube->GetRmax())*1_cm);
      Acts::Transform3 lTransf(rotateFrame * lTrans);


      const auto ppBounds =
      std::make_shared<const Acts::LineBounds>(tube->GetRmax() * 1_cm, tube->GetDZ() * 1_cm);//(bounds[0], bounds[1]);
      double thickness = 2* tube->GetRmax() * 1_cm;
      std::shared_ptr<ShipDetectorElement> detElement = nullptr;
      detElement = std::make_shared<ShipDetectorElement>(
          std::make_shared<const Acts::Transform3>(trafo), ppBounds, thickness, 
          surfaceMaterial); //No idea how thick the fibres are

      auto surface = detElement->surface().getSharedPtr();
      detElementStore.push_back(std::move(detElement));

      layerSurfaces.push_back(surface);
      
      if (n%13137 == 0){
          //Configure protolayer
          std::cout<<"Creating layer"<<std::endl; 
          std::cout<<layerSurfaces.size()<<std::endl;
          Acts::ProtoLayer pl(context, layerSurfaces);
          
          layers.push_back(layerCreator->planeLayer(context, layerSurfaces, 1, 1, Acts::AxisDirection::AxisX,pl,lTransf));//,pl)); //Three lines changed Z to X

          layerSurfaces.clear();
      }
    }


    for (auto& snode : tgpState_passive.selectedNodes) {
      const auto& transform = *snode.transform;

      Acts::RotationMatrix3 rotations = Acts::RotationMatrix3::Identity();
      const Double_t* rotation = transform.GetRotationMatrix();
      const Double_t* translation = transform.GetTranslation();

      Acts::RotationMatrix3 rotateFrame;
      rotateFrame.col(0) = Acts::Vector3(0,0,-1);
      rotateFrame.col(1) = Acts::Vector3(0,1,0);
      rotateFrame.col(2) = Acts::Vector3(1,0,0);
      positionsW.push_back(translation[2]*10); //Vector of z positions

      Acts::RotationMatrix3 combinedRot =  rotateFrame * rotations;
      Acts::Vector3 t(1_cm*translation[2],1_cm*translation[1],-1_cm*translation[0]);
      auto trafo = Acts::TGeoPrimitivesHelper::makeTransform(combinedRot.col(0),combinedRot.col(1),combinedRot.col(2),t);

      auto* tube = dynamic_cast<const TGeoBBox*>(snode.node->GetVolume()->GetShape());

      const auto pBounds =
      std::make_shared<const Acts::RectangleBounds>(tube->GetDX() * 1_cm, tube->GetDY() * 1_cm);//(bounds[0], bounds[1]);

      Acts::MutableLayerPtr nLayer = Acts::PlaneLayer::create(trafo, pBounds, nullptr, 5._mm);//, nullptr, Acts::passive));
      
      std::shared_ptr<const Acts::ISurfaceMaterial> passiveMat = passiveMaterial;//nullptr; 
      nLayer->surfaceRepresentation().assignSurfaceMaterial(passiveMat);
      passiveLayers.push_back(nLayer);
    }


//    std::cout<<"After snode loop"<<std::endl;
    std::cout<<"Length of surface array: "<<layers.size()<<std::endl;
    std::cout<<"Length of passive surface array: "<<passiveLayers.size()<<std::endl;

  Acts::Translation3 transVol(0,0, (positionsW.front() + positions.back()) * 0.5);
  //Acts::RotationMatrix3 rotations = Acts::RotationMatrix3::Identity();
      Acts::RotationMatrix3 rotations;
      rotations.col(0) = Acts::Vector3(0,0,-1);
      rotations.col(1) = Acts::Vector3(0,1,0);
      rotations.col(2) = Acts::Vector3(1,0,0);

//  Acts::RotationMatrix3 rotations2 = Acts::RotationMatrix3::Identity();
  Acts::Transform3 trafoVol(rotations * transVol);

  auto length = positions.back() - positionsW.front();

  std::shared_ptr<Acts::VolumeBounds> boundsVol = nullptr;
    boundsVol = std::make_shared<Acts::CuboidVolumeBounds>(
        // 400._mm , 400._mm, 5000._mm);
         400._mm , 400._mm, length + 100._mm);


  Acts::LayerArrayCreator::Config lacConfig;
  Acts::LayerArrayCreator layArrCreator(
      lacConfig,
      Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::VERBOSE));
  Acts::LayerVector layVec;
  for (unsigned int i = 0; i < layers.size(); i++) {
    layVec.push_back(layers[i]);
  }
  for (unsigned int i = 0; i < passiveLayers.size(); i++) {
    layVec.push_back(passiveLayers[i]);
  }
//    std::cout<<"After layer creator loop"<<std::endl;
//    std::cout<<"Length of layer array: "<<layVec.size()<<std::endl;

  // Create the layer array
  Acts::GeometryContext genGctx{context};
  
  std::unique_ptr<const Acts::LayerArray> layArr(layArrCreator.layerArray(
      genGctx, layVec, positionsW.front() - 200._mm, positions.back() + 200._mm,
      //genGctx, layVec, 0._mm,  32000._mm,
      Acts::BinningType::arbitrary, Acts::AxisDirection::AxisX)); //Changed Z to X
  
  // Build the tracking volume
  auto trackVolume = std::make_shared<Acts::TrackingVolume>(
      trafoVol, boundsVol, nullptr, std::move(layArr), nullptr,
      Acts::MutableTrackingVolumeVector{}, "SND");
  // Build and return tracking geometry
  return std::make_shared<Acts::TrackingGeometry>(trackVolume);//, std::move(mdecorator));

    }
  }
} //namespace

MTCBuilder::MTCBuilder(const Config& cfg)
    : Detector(Acts::getDefaultLogger("TGeoDetector", cfg.logLevel)),
      m_cfg(cfg) {

  m_nominalGeometryContext = Acts::GeometryContext();

  m_trackingGeometry =
      buildMTCDetector(m_cfg, m_nominalGeometryContext, m_detectorStore,
                         logger());
}
}  // namespace ActsExamples
