#include "ActsExamples/TGeoDetector/StrawtubeBuilder.hpp"
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
  std::shared_ptr<const Acts::TrackingGeometry> buildStrawDetector(
      const StrawtubeBuilder::Config& config, const Acts::GeometryContext& context,
    std::vector<std::shared_ptr<const Acts::DetectorElementBase>>&
        detElementStore,
      const Acts::Logger& logger) {


  using namespace Acts::UnitLiterals;
  TGeoManager::Import(config.fileName.c_str());

  if (gGeoManager != nullptr) {
    std::string volumeName = "*";//"*";
    Acts::TGeoParser::Options tgpOptions;
    tgpOptions.volumeNames = {volumeName};
    //tgpOptions.targetNames = {"straw_100*","straw_101*","straw_110*", "straw_111*","straw_120*","straw_121*","straw_130*","straw_131*"};
    //All the straw nodes seem to be named straw explictly...Makes it difficult to read individual layers
    tgpOptions.targetNames = {"straw"};

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

    //  Gas mixture Ar/CO2 90/10 @ 2Bar
    //  Straw inner diameter 1.9928 cm
    //  Mylar straw material thickness 0.0036 cm
    //
    Acts::Material strawTube = Acts::Material::fromMassDensity(
      28.5548_cm, 56.387_cm, 12.8772, 6.45628, 6.45628_g / 1_cm3);
    Acts::Material strawGas = Acts::Material::fromMassDensity(
      6119.96_cm, 34106.3_cm, 37.194, 16.84, 0.00336_g / 1_cm3); //Using effective quantities for the gas
    //Should also add W wire, thickness 0.003 cm//                                                             
    Acts::MaterialSlab matGas(strawGas, 0.9964_cm); //Assuming this should be the radius
    Acts::MaterialSlab matTube(strawTube, 0.0018_cm); //straw thickness of one side, i.e (r2-r1) 
    Acts::MaterialSlab combinedMaterial = Acts::MaterialSlab::combineLayers(matGas, matTube); //Combined thickness and averaged material constants
    auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(combinedMaterial);
    
    //Straws are arranged in 300/301, 315/316, 315/316, 300/301 pattern for each station.
    int n=0;
    std::vector<Acts::LayerPtr> layers;
    std::vector<std::shared_ptr<const Acts::Surface>> layerSurfaces;
    std::vector<double> positions;

    //Matrix to rotate frame by pi/2 about y-axis
    Acts::RotationMatrix3 rotateFrame;
    rotateFrame.col(0) = Acts::Vector3(0,0,-1);
    rotateFrame.col(1) = Acts::Vector3(0,1,0);
    rotateFrame.col(2) = Acts::Vector3(1,0,0);

    for (auto& snode : tgpState.selectedNodes) {
      n++;
      const auto& transform = *snode.transform;

      const Double_t* rotation = transform.GetRotationMatrix();
      const Double_t* translation = transform.GetTranslation();
      Acts::Translation3 transScaled(translation[0]*1_cm, translation[1]*1_cm, translation[2]*1_cm);
      Acts::RotationMatrix3 rotations;
      rotations.col(0) = Acts::Vector3(rotation[0],rotation[3],rotation[6]);
      rotations.col(1) = Acts::Vector3(rotation[1],rotation[4],rotation[7]);
      rotations.col(2) = Acts::Vector3(rotation[2],rotation[5],rotation[8]);

      Acts::RotationMatrix3 combinedRot = rotateFrame * rotations ;

      Acts::Vector3 t(1_cm*translation[2],1_cm*translation[1],1_cm*translation[0]);

      auto trafo = Acts::TGeoPrimitivesHelper::makeTransform(combinedRot.col(0),combinedRot.col(1),combinedRot.col(2),t);

      positions.push_back(translation[2] * 1_cm); //Vector of z positions to determine volume bounds

      auto* tube = dynamic_cast<const TGeoTube*>(snode.node->GetVolume()->GetShape());
      
      const auto ppBounds =
      std::make_shared<const Acts::LineBounds>(tube->GetRmax() * 1_cm, tube->GetDZ() * 1_cm);//(bounds[0], bounds[1]);
      std::shared_ptr<ShipDetectorElement> detElement = nullptr;
      detElement = std::make_shared<ShipDetectorElement>(
          std::make_shared<const Acts::Transform3>(trafo), ppBounds, 2._cm,
          surfaceMaterial);

      auto surface = detElement->surface().getSharedPtr();
      detElementStore.push_back(std::move(detElement));
      layerSurfaces.push_back(surface);
      
      //Define transformation for layers
      Acts::Vector3 lTrans(translation[2]*1_cm,0,0);
      auto layerTransf = Acts::TGeoPrimitivesHelper::makeTransform(rotateFrame.col(0),rotateFrame.col(1),rotateFrame.col(2),lTrans);

      //Not the most elegant implementation, but cannot use multiple parser states due to naming of nodes
      if (n == 601 || n ==1232 || n==1863 || n==2464 || n==3065 || n== 3696 || n==4327 || n==4928 || n==5529 || n==6160 || n==6791 || n==7392 || n==7993 ||n==8624 ||n==9255 ||n==9856 ){

          Acts::ProtoLayer pl(context, layerSurfaces);
          //Create plane layer comprised of sensitive straw surfaces
          layers.push_back(layerCreator->planeLayer(context, layerSurfaces, 0, 0, Acts::AxisDirection::AxisX,pl,layerTransf));
          layerSurfaces.clear();
      }
    }

    std::sort(positions.begin(), positions.end());
    Acts::Translation3 transVol(0, 0,
                               (positions.front() + positions.back()) * 0.5);
    Acts::RotationMatrix3 rotations = Acts::RotationMatrix3::Identity();
    Acts::Transform3 trafoVol(rotateFrame * transVol); 

    auto length = positions.back() - positions.front();
    // The volume bounds is set to be larger than cubic with planes
    std::shared_ptr<Acts::VolumeBounds> boundsVol = nullptr;
    boundsVol = std::make_shared<Acts::CuboidVolumeBounds>(
        605._cm , 605._cm , length + 100._cm);


    Acts::LayerArrayCreator::Config lacConfig;
    Acts::LayerArrayCreator layArrCreator(
      lacConfig,
      Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::VERBOSE));

    Acts::LayerVector layVec;
    for (unsigned int i = 0; i < layers.size(); i++) {
      layVec.push_back(layers[i]);
    }
      ACTS_INFO("Length of layer array: " << layVec.size());

    // Create the layer array
    Acts::GeometryContext genGctx{context};
    
    std::unique_ptr<const Acts::LayerArray> layArr(layArrCreator.layerArray(
        genGctx, layVec, positions.front() - 50._cm, positions.back() + 50._cm,
        Acts::BinningType::arbitrary, Acts::AxisDirection::AxisX));
    
    // Build the tracking volume
    auto trackVolume = std::make_shared<Acts::TrackingVolume>(
        trafoVol, boundsVol, nullptr, std::move(layArr), nullptr,
        Acts::MutableTrackingVolumeVector{}, "SST");

    // Build and return tracking geometry
    return std::make_shared<Acts::TrackingGeometry>(trackVolume);

    }
  }
} //namespace

StrawtubeBuilder::StrawtubeBuilder(const Config& cfg)
    : Detector(Acts::getDefaultLogger("TGeoDetector", cfg.logLevel)),
      m_cfg(cfg) {

  m_nominalGeometryContext = Acts::GeometryContext();

  m_trackingGeometry =
      buildStrawDetector(m_cfg, m_nominalGeometryContext, m_detectorStore,
                         logger());
}
}  // namespace ActsExamples
