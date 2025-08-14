#include "ActsExamples/TGeoDetector/SiTargetBuilder.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Root/TGeoParser.hpp"
#include "Acts/Plugins/Root/TGeoSurfaceConverter.hpp"
#include "Acts/Plugins/Root/TGeoDetectorElement.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/ITrackingVolumeHelper.hpp"
#include "Acts/Geometry/ProtoLayerHelper.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Plugins/Root/TGeoPrimitivesHelper.hpp"
#include "Acts/Plugins/Root/TGeoMaterialConverter.hpp"
#include "Acts/Plugins/Root/TGeoSurfaceConverter.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "TGeoManager.h"


namespace ActsExamples {


namespace {    

std::shared_ptr<const Acts::TrackingGeometry> buildSiTargetDetector(
    const SiTargetBuilder::Config& config, const Acts::GeometryContext& context,
    std::vector<std::shared_ptr<const Acts::DetectorElementBase>>&
        detElementStore,
    const Acts::Logger& logger) {


    using namespace Acts::UnitLiterals;
    TGeoManager::Import(config.fileName.c_str());

  if (gGeoManager != nullptr) {

    std::string volumeName = "*";
    std::string axes = "XYZ";
    double scale = 10.;

    Acts::TGeoParser::Options tgpOptions;
    tgpOptions.volumeNames = {volumeName};
    tgpOptions.targetNames = {"SensorVolume"};

    Acts::TGeoParser::Options tgpOptions_passive;
    tgpOptions_passive.volumeNames = {volumeName};
    tgpOptions_passive.targetNames = {"volTargetWall"};

    Acts::TGeoParser::State tgpState;
    tgpState.volume = gGeoManager->GetTopVolume();
    Acts::TGeoParser::select(tgpState, tgpOptions);

    Acts::TGeoParser::State tgpState_passive;
    tgpState_passive.volume = gGeoManager->GetTopVolume();
    Acts::TGeoParser::select(tgpState_passive, tgpOptions_passive);

    Acts::ProtoLayerHelper::Config plhConfig;
    Acts::ProtoLayerHelper plHelper(
      plhConfig, Acts::getDefaultLogger("ProtoLayerHelper", Acts::Logging::VERBOSE));

//    Acts::LayerArrayCreator::Config lacConfig;
//  auto layerArrayCreator = std::make_shared<const Acts::LayerArrayCreator>(
//      lacConfig, logger.clone("LayerArrayCreator", config.layerLogLevel));

  Acts::LayerCreator::Config lcConfig;
  lcConfig.surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>();
  auto layerCreator = std::make_shared<const Acts::LayerCreator>(
      lcConfig, logger.clone("LayerCreator", config.layerLogLevel));

    std::vector<double> positions;
    std::vector<double> positionsW;
//    ObjVisualization3D objVis;
//    Acts::MaterialSlab matProp(makeSilicon(), 300_um);
    Acts::Material silicon = Acts::Material::fromMassDensity(
      9.370_cm, 46.52_cm, 28.0855, 14, 2.329_g / 1_cm3);
    Acts::MaterialSlab matProp(silicon, 300_um);
    auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);

    Acts::Material tungsten = Acts::Material::fromMassDensity(
      0.35_cm, 9.9_cm, 183.84, 74, 19.254_g / 1_cm3);
    Acts::MaterialSlab matPropW(tungsten, 1.75_mm);
    auto passiveMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matPropW);


      /*
      auto tgeoMaterial = tgpState.selectedNodes[0].GetMedium()->GetMaterial();
      // Convert the material
      Acts::TGeoMaterialConverter::Options materialOptions;
      materialOptions.unitLengthScalor = 10; //cm to mm
      auto materialSlab = TGeoMaterialConverter::materialSlab(
        *tgeoMaterial, 300_um, 1_mm,
        materialOptions);
    auto surfaceMaterial =
        std::make_shared<HomogeneousSurfaceMaterial>(materialSlab);
      */

    int n=0;
//    std::vector<std::shared_ptr<const Surface>> m_surfaces;
//      std::size_t nLayers = 800;//positions.size();
      std::vector<Acts::LayerPtr> layers;//nLayers);
      using LayerSurfaceVector = std::vector<std::shared_ptr<const Acts::Surface>>;
      LayerSurfaceVector layerSurfaces;
    for (auto& snode : tgpState.selectedNodes) {
      const auto& shape = *(snode.node->GetVolume()->GetShape());
      const auto& transform = *snode.transform;
      n++;
      auto identifier = Acts::TGeoDetectorElement::Identifier();
      auto tgElement = config.elementFactory(
            identifier, *snode.node, *snode.transform, "XYZ",
            10, surfaceMaterial);
//              scale, nullptr);


      //snode.GetMedium()->GetMaterial();
//      surface->assignSurfaceMaterial(surfaceMaterial);

  const auto pBounds =
      std::make_shared<const Acts::RectangleBounds>(100._mm,200._mm);//(bounds[0], bounds[1]);

      Acts::RotationMatrix3 rotations = Acts::RotationMatrix3::Identity();
      const Double_t* rotation = transform.GetRotationMatrix();
      const Double_t* translation = transform.GetTranslation();
      Acts::Translation3 trans(translation[0]*10, translation[1]*10, translation[2]*10);
      positions.push_back(translation[2]*10); //Vector of z positions
      std::cout<< translation[2]*10<<std::endl;
      //      Acts::Translation3 trans(translation[0]*10, translation[1]*10, translation[2]*10);
//      positions.push_back(translation[2]*10); //Vector of z positions
      Acts::Transform3 trafo(rotations * trans);



//      auto surface = tgElement->surface().getSharedPtr();
      layerSurfaces.push_back(tgElement->surface().getSharedPtr());
      detElementStore.push_back(std::move(tgElement));
      //detectorStore.push_back(tgElement->surface().getSharedPtr());



          std::cout<<layerSurfaces.size()<<std::endl;

      if (n%8 == 0){
          //Configure protolayer
          auto protolayer = plHelper.protoLayers(context, Acts::unpackSmartPointers(layerSurfaces),Acts::ProtoLayerHelper::SortingConfig(Acts::AxisDirection::AxisX, 20.)); 
                  
          std::cout<<protolayer.size()<<std::endl;

          for (auto& pLayer : protolayer) {
             layerSurfaces.clear();

          for (const auto& lsurface : pLayer.surfaces()) {
            layerSurfaces.push_back(lsurface->getSharedPtr());
                }
            }

//          layers.push_back(Acts::PlaneLayer::create(trafo, pBounds, layerSurfaces, 300._um));
          Acts::ProtoLayer pl(context, layerSurfaces);
          pl.envelope[Acts::AxisDirection::AxisZ] = {-0.4,
                                                 0.4};
          if(n%2==1){
             layers.push_back(layerCreator->planeLayer(context, layerSurfaces, 4, 2, Acts::AxisDirection::AxisZ));//,pl));
          } else {
             layers.push_back(layerCreator->planeLayer(context, layerSurfaces, 2, 4, Acts::AxisDirection::AxisZ));//,pl));
          }

          layerSurfaces.clear();
      }
    }
//Creation of passive layers vector

    std::vector<Acts::LayerPtr> passiveLayers;

    for (auto& snode : tgpState_passive.selectedNodes) {
      const auto& transform = *snode.transform;
      auto identifier = Acts::TGeoDetectorElement::Identifier();
      auto tgElement = config.elementFactory(
            identifier, *snode.node, *snode.transform, "XYZ",
            10, surfaceMaterial);

  const auto pBounds =
      std::make_shared<const Acts::RectangleBounds>(400._mm,400._mm);

      Acts::RotationMatrix3 rotations = Acts::RotationMatrix3::Identity();
      const Double_t* rotation = transform.GetRotationMatrix();
      const Double_t* translation = transform.GetTranslation();
      Acts::Translation3 trans(translation[0]*10, translation[1]*10, translation[2]*10);
//      Acts::Translation3 trans(translation[0]*10, translation[1]*10, translation[2]*10+ 3100._cm);
//      positions.push_back(translation[2]); //Vector of z positions
      positionsW.push_back(translation[2]*10 ); //Vector of z positions
      std::cout<<translation[2] * 10 <<std::endl;
      Acts::Transform3 trafo(rotations * trans);

//      std::cout<<translation[2]<<std::endl;
      auto surface = tgElement->surface().getSharedPtr();
      detElementStore.push_back(std::move(tgElement));

      std::unique_ptr<Acts::SurfaceArray> surArray(new Acts::SurfaceArray(surface));
      passiveLayers.push_back(Acts::PlaneLayer::create(trafo, pBounds, nullptr, 1.75_mm, nullptr, Acts::passive));

      auto mutableSurface = const_cast<Acts::Surface*>(surface.get());
      mutableSurface->associateLayer(*passiveLayers[-1]);
    }


//    std::cout<<"After snode loop"<<std::endl;
    std::cout<<"Length of surface array: "<<positions.size()<<std::endl;

  Acts::Translation3 transVol(0, 0,
                               (positionsW.front() + positions.back()) * 0.5);
      Acts::RotationMatrix3 rotations = Acts::RotationMatrix3::Identity();
  Acts::Transform3 trafoVol(rotations * transVol);

  auto length = positions.back() - positionsW.front();

  std::shared_ptr<Acts::VolumeBounds> boundsVol = nullptr;
    boundsVol = std::make_shared<Acts::CuboidVolumeBounds>(
        400._mm , 400._mm , length);


  Acts::LayerArrayCreator::Config lacConfig;
  Acts::LayerArrayCreator layArrCreator(
      lacConfig,
      Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::INFO));
  Acts::LayerVector layVec;
  for (unsigned int i = 0; i < layers.size(); i++) {
    layVec.push_back(layers[i]);
  }
  for (unsigned int i = 0; i < passiveLayers.size(); i++) {
//      std::cout<<"W layer"<<std::endl;
    layVec.push_back(passiveLayers[i]);
  }
    std::cout<<"After layer creator loop"<<std::endl;
    std::cout<<"Length of layer array: "<<layVec.size()<<std::endl;

  // Create the layer array
  Acts::GeometryContext genGctx{context};
  
  std::unique_ptr<const Acts::LayerArray> layArr(layArrCreator.layerArray(
      genGctx, layVec, positionsW.front() - 2._mm, positions.back() + 2._mm,
      Acts::BinningType::arbitrary, Acts::AxisDirection::AxisZ));
  
  // Build the tracking volume
  auto trackVolume = std::make_shared<Acts::TrackingVolume>(
      trafoVol, boundsVol, nullptr, std::move(layArr), nullptr,
      Acts::MutableTrackingVolumeVector{}, "SND");
  // Build and return tracking geometry
  return std::make_shared<Acts::TrackingGeometry>(trackVolume);//, std::move(mdecorator));

    }
  }
} //namespace

SiTargetBuilder::SiTargetBuilder(const Config& cfg)
    : Detector(Acts::getDefaultLogger("TGeoDetector", cfg.logLevel)),
      m_cfg(cfg) {

  m_nominalGeometryContext = Acts::GeometryContext();

  m_trackingGeometry =
      buildSiTargetDetector(m_cfg, m_nominalGeometryContext, m_detectorStore,
                        logger());
}
}  // namespace ActsExamples
