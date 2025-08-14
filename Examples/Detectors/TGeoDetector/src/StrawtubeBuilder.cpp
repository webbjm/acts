#include "ActsExamples/TGeoDetector/StrawtubeBuilder.hpp"

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


namespace ActsExamples {


namespace {    
/// @brief struct to load the global geometry
  std::shared_ptr<const Acts::TrackingGeometry> buildStrawDetector(
      const StrawtubeBuilder::Config& config, const Acts::GeometryContext& context,
      const Acts::Logger& logger) {


  using namespace Acts::UnitLiterals;
  TGeoManager::Import(config.fileName.c_str());

  if (gGeoManager != nullptr) {
    std::string volumeName = "*";
    Acts::TGeoParser::Options tgpOptions;
    tgpOptions.volumeNames = {volumeName};
    tgpOptions.targetNames = {"straw*"};

    Acts::TGeoParser::State tgpState;
    tgpState.volume = gGeoManager->GetTopVolume();
    Acts::TGeoParser::select(tgpState, tgpOptions);

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
    
    //There should be 300 wires per layer, but only 299 in geofile, bug in detector implementation?
    int n=0;
    std::vector<Acts::LayerPtr> layers;
    std::vector<std::shared_ptr<const Acts::Surface>> layerSurfaces;
    std::vector<double> positions;

    for (auto& snode : tgpState.selectedNodes) {
      n++;
      const auto& shape = (snode.node->GetVolume()->GetShape());
      const auto& transform = *snode.transform;

      const Double_t* rotation = transform.GetRotationMatrix();
      const Double_t* translation = transform.GetTranslation();
      Acts::Translation3 transScaled(translation[0]*1_cm, translation[1]*1_cm, translation[2]*1_cm);
      Acts::RotationMatrix3 rotations;
      rotations.col(0) = Acts::Vector3(rotation[0],rotation[3],rotation[6]);
      rotations.col(1) = Acts::Vector3(rotation[1],rotation[4],rotation[7]);
      rotations.col(2) = Acts::Vector3(rotation[2],rotation[5],rotation[8]);
      Acts::Transform3 trafo(transScaled * rotations);
      positions.push_back(translation[2]*1_cm); //Vector of z positions to determine volume bounds

      auto tube = dynamic_cast<const TGeoTube*>(shape);
      //transformation matrix, radius, half length, scale from cm -> mm
      auto aStraw = Acts::Surface::makeShared<Acts::StrawSurface>(trafo, tube->GetRmax() * 1_cm, tube->GetDZ() * 1_cm); 
      aStraw->assignSurfaceMaterial(surfaceMaterial);
      layerSurfaces.push_back(aStraw->getSharedPtr());
      //Not the most elegant implementation, consider using multiple parser states
      if (n%598 == 0){
          //Create plane layer comprised of sensitive straw surfaces
          layers.push_back(layerCreator->planeLayer(context, layerSurfaces, 1, 598, Acts::AxisDirection::AxisZ));
          layerSurfaces.clear();
      }
    }

  // The volume transform
  Acts::Translation3 transVol(0, 0,
                               (positions.front() + positions.back()) * 0.5);
  Acts::RotationMatrix3 rotations = Acts::RotationMatrix3::Identity();
  Acts::Transform3 trafoVol(rotations * transVol);

  auto length = positions.back() - positions.front();
  // The volume bounds is set to be larger than cubic with planes
  std::shared_ptr<Acts::VolumeBounds> boundsVol = nullptr;
    boundsVol = std::make_shared<Acts::CuboidVolumeBounds>(
        1000._cm , 1000._cm , length);


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
      genGctx, layVec, positions.front(), positions.back(),
      Acts::BinningType::arbitrary, Acts::AxisDirection::AxisZ));
  
  // Build the tracking volume
  auto trackVolume = std::make_shared<Acts::TrackingVolume>(
      trafoVol, boundsVol, nullptr, std::move(layArr), nullptr,
      Acts::MutableTrackingVolumeVector{}, "SST");

  // Build and return tracking geometry
  return std::make_shared<Acts::TrackingGeometry>(trackVolume);//, std::move(mdecorator));

    }
  }
} //namespace

StrawtubeBuilder::StrawtubeBuilder(const Config& cfg)
    : Detector(Acts::getDefaultLogger("TGeoDetector", cfg.logLevel)),
      m_cfg(cfg) {

  m_nominalGeometryContext = Acts::GeometryContext();

  m_trackingGeometry =
      buildStrawDetector(m_cfg, m_nominalGeometryContext, 
                         logger());
}
}  // namespace ActsExamples
