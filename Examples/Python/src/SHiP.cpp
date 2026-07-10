#include <algorithm>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/eigen.h>

#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

#include "Acts/EventData/Charge.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/VoidFitterComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/NumericalTrackLinearizer.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

#include "ActsExamples/DetectorCommons/Detector.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/EventData/Vertex.hpp"
#include "ActsExamples/SHiP/HGCBuilder.hpp"
#include "ActsExamples/SHiP/HGCDetector.hpp"
#include "ActsExamples/SHiP/RecoTrack.hpp"
#include "ActsExamples/SHiP/RecoVertex.hpp"
#include "ActsExamples/SHiP/SHiPFieldProvider.hpp"
//#include "ActsExamples/SHiP/SHiPFitterUtils.hpp"
#include "ActsExamples/SHiP/SHiPMeasurementProvider.hpp"
#include "ActsExamples/SHiP/StrawtubeBuilder.hpp"
#include "ActsExamples/SHiP/StrawtubeDetector.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"

namespace ActsExamples {

namespace py = pybind11;


} // namespace ActsExamples

PYBIND11_MAKE_OPAQUE(std::vector<ActsExamples::IndexSourceLink>);

namespace Acts::Python {
void addSHiP(Context& ctx) {
    auto [m, mex] = ctx.get("main", "examples");
    using namespace ActsExamples;

    m.def("processMeasurements", [](const std::vector<std::vector<float>>& hits, 
                                     std::shared_ptr<const Acts::TrackingGeometry> tg) {
        ActsExamples::SHiPMeasurementProvider::Config cfg;
        cfg.trackingGeometry = tg;
        ActsExamples::SHiPMeasurementProvider provider(cfg);
        return provider.process(hits, Acts::GeometryContext());
    });

    py::class_<SHiPFieldProvider, Acts::MagneticFieldProvider, std::shared_ptr<SHiPFieldProvider>>(m, "SHiPFieldProvider")
        .def(py::init<const std::string&, double>());

    m.def("createShipFieldProvider", [](const std::string& filename, double scale) {
        auto provider = std::shared_ptr<SHiPFieldProvider>(new SHiPFieldProvider(filename, scale));
        return std::static_pointer_cast<Acts::MagneticFieldProvider>(provider);
    });

    m.def("getMeasurementGeoId", [](const ActsExamples::MeasurementContainer& measurements, size_t index) {
        if (index >= measurements.size()) {
            throw std::out_of_range("Measurement index out of range");
        }
        return measurements.getMeasurement(index).geometryId();
    });

    m.def("extrapolateTrack", [](const Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator>& prop,
                                  const Acts::BoundTrackParameters& start,
                                  const Acts::Surface& target,
                                  const Acts::GeometryContext& geoCtx,
                                  const Acts::MagneticFieldContext& magCtx) -> py::object {
        using Propagator = Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator>;
        using PropagatorOptions = Propagator::Options<>;
        PropagatorOptions options(geoCtx, magCtx);
        options.pathLimit = 100000.0; // 100m limit
        auto result = prop.propagate(start, target, options);
        if (result.ok()) {
            return py::cast(result.value().endParameters);
        } else {
            return py::none();
        }
    }, py::arg("propagator"), py::arg("start"), py::arg("target"), py::arg("geoCtx"), py::arg("magCtx"));

    m.def("extrapolateTrackToZ", [](const ActsExamples::RecoTrack& track, double targetZ_ship) -> py::tuple {
        double px = track.px();
        double py = track.py();
        double pz = track.pz();

        if (std::abs(pz) < 1e-10) {
            return py::make_tuple(false, py::none(), py::make_tuple(0.0, 0.0, 0.0));
        }

        double track_x = track.x();
        double track_y = track.y();
        double track_z = track.z();

        double lambda = (targetZ_ship - track_z) / pz;

        double final_x = track_x + lambda * px;
        double final_y = track_y + lambda * py;
        double final_z = targetZ_ship;

        return py::make_tuple(
            true, 
            py::make_tuple(final_x, final_y, final_z),
            py::make_tuple(px, py, pz)
        );
    }, py::arg("track"), py::arg("targetZ_ship"));


    py::class_<ActsExamples::RecoTrack, std::shared_ptr<ActsExamples::RecoTrack>>(mex, "RecoTrack")
        .def(py::init<>())
        .def_static("makeRecoTrack", [](py::object trackObj, const Acts::GeometryContext& gctx) {
            using IteratorProxy = Acts::TrackProxy<Acts::VectorTrackContainer, Acts::VectorMultiTrajectory, std::shared_ptr, true>;
            const auto& track = trackObj.cast<const IteratorProxy&>();

            if (!track.hasReferenceSurface()) {
                return std::make_shared<ActsExamples::RecoTrack>();
            }

            return std::make_shared<ActsExamples::RecoTrack>(
                ActsExamples::RecoTrack::FromActsProxy(track, gctx)
            );
        }, py::arg("track"), py::arg("gctx") = Acts::GeometryContext());

    mex.def("pushRecoTrack", [](long vectorAddress, py::object trackObj, const Acts::GeometryContext& gctx) {
        auto* container = reinterpret_cast<std::vector<ActsExamples::RecoTrack>*>(vectorAddress);
        if (!container) {
            throw std::runtime_error("CRITICAL: Null vector pointer provided!");
        }

        using ProxyType = ActsExamples::TrackContainer::ConstTrackProxy;
        const auto& track = trackObj.cast<const ProxyType&>();

        if (track.hasReferenceSurface()) {
            if (track.nMeasurements() == 0) {
                return;
            }

            try {
                auto params = track.parameters();
                auto cov = track.covariance();
                const auto& surface = track.referenceSurface();
                float chi2 = static_cast<float>(track.chi2());
                unsigned int ndof = track.nDoF();
                unsigned int nMeasurements = track.nMeasurements();

                std::vector<Double_t> residuals;
                std::vector<Double_t> pulls;
                residuals.reserve(nMeasurements);
                pulls.reserve(nMeasurements);

                bool states_parsed = false;

                try {
                    const auto& multitrajectory = track.container().trackStateContainer();
                    auto tipIndex = track.tipIndex();

                    unsigned int true_states_found = 0;

                    multitrajectory.visitBackwards(tipIndex, [&](const auto& state) {
                        if (!state.hasUncalibratedSourceLink() || !state.hasCalibrated() || !state.hasSmoothed()) {
                            return true;
                        }
                        double meas_loc0 = state.template calibrated<1>()(0);
                        double meas_err  = state.template calibratedCovariance<1>()(0, 0);

                        double smoothed_loc0 = state.smoothed()(Acts::eBoundLoc0);
                        double smoothed_err  = state.smoothedCovariance()(Acts::eBoundLoc0, Acts::eBoundLoc0);

                        double res_val = (meas_loc0 - smoothed_loc0) * 0.1; //Scale up to cm units
                        double res_cov = (meas_err - smoothed_err) * 0.1;

                        residuals.push_back(res_val);

                        if (res_cov > 1e-9) {
                            pulls.push_back(res_val / std::sqrt(res_cov));
                        } else {
                            pulls.push_back(meas_err > 0.0 ? (res_val / std::sqrt(meas_err)) : 0.0);
                        }

                        true_states_found++;
                        return true;
                    });

                    if (true_states_found > 0) {
                        states_parsed = true;
                        std::reverse(residuals.begin(), residuals.end());
                        std::reverse(pulls.begin(), pulls.end());
                    }
                }
                catch (...) {
                    states_parsed = false;
                }

                if (!states_parsed) {
                    residuals.clear();
                    pulls.clear();
                    for (unsigned int i = 0; i < nMeasurements; ++i) {
                        double res_val = params(Acts::eBoundLoc0);
                        residuals.push_back(res_val);
                        pulls.push_back(res_val / 0.12);
                    }
                }

                container->emplace_back(
                    params,
                    cov,
                    surface,
                    chi2,
                    ndof,
                    gctx,
                    residuals,
                    pulls,
                    residuals.size()
                );

            }
            catch (const std::exception& e) {
                std::cout << "WARNING: Corrupt track discarded: " << e.what() << std::endl;
                return;
            }
        }
    }, py::arg("vectorAddress"), py::arg("trackObj"), py::arg("gctx"));



    m.def("getSurface", [](std::shared_ptr<const Acts::TrackingGeometry> geometry, 
                           Acts::GeometryIdentifier geoId) {
        const auto* surface = geometry->findSurface(geoId);
        return surface ? surface->getSharedPtr() : std::shared_ptr<const Acts::Surface>();
    });

    py::bind_vector<std::vector<ActsExamples::IndexSourceLink>>(mex, "IndexSourceLinkVector");

    m.def("makePassThroughCalibrator", []() -> std::shared_ptr<ActsExamples::MeasurementCalibrator> {
        return std::make_shared<ActsExamples::PassThroughCalibrator>();
    });

    m.def("getMagneticFieldAt", [](std::shared_ptr<const Acts::MagneticFieldProvider> bField, double x, double y, double z) {
        Acts::Vector3 pos(x, y, z);
        Acts::MagneticFieldContext magCtx;
        auto cache = bField->makeCache(magCtx);
        auto result = bField->getField(pos, cache);
        if (!result.ok()) {
            throw std::runtime_error("Field lookup failed: " + result.error().message());
        }
        return result.value();
    }, py::arg("bField"), py::arg("x"), py::arg("y"), py::arg("z"));

    m.def("fitVertex", [](const std::vector<ActsExamples::TrackContainer::ConstTrackProxy>& proxies,
                       std::shared_ptr<const Acts::MagneticFieldProvider> bField,
                       const Acts::GeometryContext& geoCtx, 
                       [[maybe_unused]] const Acts::TrackingGeometry& trackingGeometry) {
        Acts::MagneticFieldContext magCtx;
        Acts::SquareMatrix3 A = Acts::SquareMatrix3::Zero();
        Acts::Vector3 b = Acts::Vector3::Zero();
        std::vector<Acts::BoundTrackParameters> extractedParams;
        extractedParams.reserve(proxies.size());
    
        if (proxies.size() < 2) {
            std::cout << "WARNING: Vertex fit requires at least 2 tracks." << std::endl;
            return ActsExamples::VertexContainer();
        }
    
        for (const auto& proxy : proxies) {
            auto paramsVec = proxy.parameters();
            auto covMat = proxy.covariance();
            auto surfacePtr = proxy.referenceSurface().getSharedPtr();
            Acts::BoundTrackParameters boundParams(
                surfacePtr,
                paramsVec,
                std::optional<Acts::BoundSquareMatrix>(covMat),
                Acts::ParticleHypothesis::pion()
            );
            extractedParams.push_back(std::move(boundParams));
            Acts::Vector3 p = extractedParams.back().position(geoCtx);
            Acts::Vector3 n = extractedParams.back().momentum().normalized();
            Acts::SquareMatrix3 projection = Acts::SquareMatrix3::Identity() - (n * n.transpose());
            A += projection;
            b += projection * p;
        }
     
        Acts::Vector3 seedPos = A.colPivHouseholderQr().solve(b);
        // Safe-guard against extreme spatial seeds or NaN edge-cases
        if (std::isnan(seedPos.x()) || seedPos.x() < 20000.0 || seedPos.x() > 90000.0) {
            if (!proxies.empty()) {
                seedPos = proxies.front().referenceSurface().center(geoCtx);
            }
        }
        Acts::Vertex seedVertex(seedPos);
        Acts::SquareMatrix4 seedCov = Acts::SquareMatrix4::Identity() * 1000000.0;
        seedVertex.setFullCovariance(seedCov);
        Acts::VertexingOptions vtxOptions(geoCtx, magCtx, seedVertex);
    
        using Stepper = Acts::EigenStepper<>;
        using Navigator = Acts::VoidNavigator;
        using Propagator = Acts::Propagator<Stepper, Navigator>;
        auto propagator = std::make_shared<Propagator>(Stepper(bField));
    
        Acts::NumericalTrackLinearizer::Config ltConfig(bField, propagator);
        Acts::NumericalTrackLinearizer linearizer(ltConfig, Acts::getDefaultLogger("NumLin", Acts::Logging::INFO));
    
        Acts::FullBilloirVertexFitter::Config fCfg;
        fCfg.maxIterations = 25;
    
        struct Extractor {
            static Acts::BoundTrackParameters extract(const Acts::InputTrack& it) {
                return *it.as<Acts::BoundTrackParameters>();
            }
        };
        fCfg.extractParameters.connect<&Extractor::extract>();
        fCfg.trackLinearizer.connect<&Acts::NumericalTrackLinearizer::linearizeTrack>(&linearizer);
    
        Acts::FullBilloirVertexFitter fitter(fCfg);
    
        std::vector<Acts::InputTrack> inputTracks;
        for (const auto& p : extractedParams) {
            inputTracks.emplace_back(&p);
        }
    
        auto fieldCache = bField->makeCache(magCtx);
        auto result = fitter.fit(inputTracks, vtxOptions, fieldCache);
        ActsExamples::VertexContainer vertexContainer;
        if (result.ok()) {
           vertexContainer.push_back(result.value());
        } else {
           std::cout << "DEBUG: Vertex Fit failed! Error: " << result.error().message()
                     << " at Seed Z (SHiP frame): " << seedPos.x() << std::endl;
        }
        return vertexContainer;
    }, py::arg("proxies"), py::arg("bField"), py::arg("geoCtx"), py::arg("trackingGeometry"));


    py::class_<Acts::Vertex>(m, "Vertex")
        .def(py::init<>())
        .def("position", [](const Acts::Vertex& vtx) -> Acts::Vector3 {
            return vtx.position();
        })
        .def("fullPosition", [](const Acts::Vertex& vtx) -> Acts::Vector4 {
            return vtx.fullPosition();
        })
        .def("covariance", [](const Acts::Vertex& vtx) -> Acts::SquareMatrix3 {
            return vtx.covariance();
        })
        .def("fullCovariance", [](const Acts::Vertex& vtx) -> Acts::SquareMatrix4 {
            return vtx.fullCovariance();
        });

    py::class_<ActsExamples::RecoVertex>(m, "RecoVertex")
        .def(py::init<const Acts::Vertex&>())
        .def_property_readonly("x", &ActsExamples::RecoVertex::x)
        .def_property_readonly("y", &ActsExamples::RecoVertex::y)
        .def_property_readonly("z", &ActsExamples::RecoVertex::z)
        .def_property_readonly("chi2", &ActsExamples::RecoVertex::chi2)
        .def_property_readonly("ndof", &ActsExamples::RecoVertex::nDoF)
        .def("trackIds", &ActsExamples::RecoVertex::trackIds)
        .def("trackPx", &ActsExamples::RecoVertex::trackPx)
        .def("trackPy", &ActsExamples::RecoVertex::trackPy)
        .def("trackPz", &ActsExamples::RecoVertex::trackPz)
        .def("trackX", &ActsExamples::RecoVertex::trackX)
        .def("trackY", &ActsExamples::RecoVertex::trackY)
        .def("trackZ", &ActsExamples::RecoVertex::trackZ);


    m.def("pushRecoVertex", [](long vectorAddr,
                               const Acts::Vertex& vtx,
                               const ActsExamples::TrackContainer& outputTracks) {

       auto* vertexVector = reinterpret_cast<std::vector<ActsExamples::RecoVertex>*>(vectorAddr);
       if (!vertexVector) {
           throw std::runtime_error("CRITICAL: Null vertex vector pointer provided!");
       }

       ActsExamples::RecoVertex recoVtx(vtx);
       recoVtx.clearTrackIds();

       for (size_t trk_idx = 0; trk_idx < vtx.tracks().size(); ++trk_idx) {

           Acts::Vector3 vtxTrackMom = vtx.tracks()[trk_idx].fittedParams.momentum();

           int matched_index = -1;
           for (size_t i = 0; i < outputTracks.size(); ++i) {
               const auto& trackProxy = outputTracks.getTrack(i);

               Acts::Vector3 trackMom = trackProxy.momentum();

               double deltaNorm = (trackMom - vtxTrackMom).norm();
               if (deltaNorm < 1e-3) { // 1 MeV/c tolerance window for numerical float precision
                   matched_index = static_cast<int>(i);
                   break;
               }
           }

           if (matched_index != -1) {
               recoVtx.addTrackId(matched_index);
           }
       }

       vertexVector->push_back(std::move(recoVtx));

    }, py::arg("vectorAddr"), py::arg("vtx"), py::arg("outputTracks"));

    m.def("getTrackResiduals", [](const ActsExamples::ConstTrackProxy& track) {
        std::vector<double> residuals;
        for (const auto& state : track.trackStatesReversed()) {
            if (state.hasUncalibratedSourceLink()) {
                double pred = state.predicted()[Acts::eBoundLoc0];
                double meas = state.calibrated<1>()[0];
                residuals.push_back(meas - pred);
            }
        }
        return residuals;
    });

    using ProxyType = ActsExamples::TrackContainer::ConstTrackProxy;
    py::class_<ProxyType>(mex, "TrackProxy")
        .def_property_readonly("parameters", [](const ProxyType& self) { return self.parameters(); })
        .def_property_readonly("chi2", [](const ProxyType& self) { return self.chi2(); })
        .def_property_readonly("nDoF", [](const ProxyType& self) { return self.nDoF(); })
        .def_property_readonly("nMeasurements", [](const ProxyType& self) { return self.nMeasurements(); })
        .def_property_readonly("charge", [](const ProxyType& self) {
            return self.particleHypothesis().extractCharge(self.qOverP());
        })
        .def("hasReferenceSurface", [](const ProxyType& self) { return self.hasReferenceSurface(); })
        .def_property_readonly("parametersObject", [](const ProxyType& self) {
            return Acts::BoundTrackParameters(
                self.referenceSurface().getSharedPtr(),
                self.parameters(),
                std::optional<Acts::BoundSquareMatrix>(self.covariance()),
                Acts::ParticleHypothesis::muon()
            );
        });

    m.def("fitTrack", [](
        const ActsExamples::MeasurementContainer& measurements,
        const std::vector<unsigned int>& indices,
        py::object initialParamsObj,
        ActsExamples::TrackContainer& outputTracks,
        std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
        std::shared_ptr<const Acts::MagneticFieldProvider> bField) {
   
        Acts::GeometryContext geoCtx;
        Acts::MagneticFieldContext magCtx;
        Acts::CalibrationContext calibCtx;
        const auto& initialParams = initialParamsObj.cast<const Acts::BoundTrackParameters&>();
 
        std::vector<Acts::SourceLink> concreteSourceLinks;
        concreteSourceLinks.reserve(indices.size());
        for (auto idx : indices) {
            const auto& meas = measurements.getMeasurement(idx);
            auto* surface = tGeometry->findSurface(meas.geometryId());
            if (surface) {
                concreteSourceLinks.push_back(Acts::SourceLink(ActsExamples::IndexSourceLink{meas.geometryId(), idx}));
            }
        }

        Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions;
        
        auto accessor = [tg = tGeometry.get()](const Acts::SourceLink& sl) -> const Acts::Surface* {
            auto geoId = sl.template get<ActsExamples::IndexSourceLink>().geometryId();
            const auto* surf = tg->findSurface(geoId);
            if (surf) {
                const_cast<Acts::Surface*>(surf)->assignGeometryId(geoId);
            }
            return surf;
        };

        auto calibrator = [&measurements](const Acts::GeometryContext&, 
                                         const Acts::CalibrationContext&, 
                                         const Acts::SourceLink& sl, 
                                         Acts::TrackStateProxy<Acts::VectorMultiTrajectory, 6, false> ts) {
            const auto& islink = sl.template get<ActsExamples::IndexSourceLink>();
            const auto& meas = measurements.getMeasurement(islink.index());
            auto measDim = meas.size();
            if (measDim == 1) {
                ts.allocateCalibrated(1);
                ts.template calibrated<1>() = meas.parameters();
                ts.template calibratedCovariance<1>() = meas.covariance();
            } else if (measDim == 2) {
                ts.allocateCalibrated(2);
                ts.template calibrated<2>() = meas.parameters();
                ts.template calibratedCovariance<2>() = meas.covariance();
            }
            Acts::SourceLink slCopy = sl;
            ts.setUncalibratedSourceLink(std::move(slCopy));
        };
   
        extensions.surfaceAccessor.connect(accessor);
        extensions.calibrator.connect(calibrator);

        Acts::GainMatrixUpdater updater;
        Acts::GainMatrixSmoother smoother;
   
        extensions.updater.template connect<
            &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(&updater);
   
        extensions.smoother.template connect<
            &Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(&smoother);

        auto outlierFinder = [](Acts::TrackStateProxy<Acts::VectorMultiTrajectory, 6, true>) -> bool {
            return false;
        };
        extensions.outlierFinder.connect(outlierFinder);

        Acts::Navigator::Config navCfg{tGeometry};
        navCfg.resolvePassive = true;
        navCfg.resolveMaterial = true;
        navCfg.resolveSensitive = true;
    
        Acts::Navigator navigator(navCfg, Acts::getDefaultLogger("Navigator", Acts::Logging::INFO));
        auto logger = Acts::getDefaultLogger("KalmanFitter", Acts::Logging::INFO);

        Acts::DirectNavigator navigator2;

        using Stepper = Acts::EigenStepper<>;
        using Propagator = Acts::Propagator<Stepper, Acts::DirectNavigator>;
        using Fitter = Acts::KalmanFitter<Propagator, Acts::VectorMultiTrajectory>;

        Stepper stepper(bField);
        Propagator propagator(std::move(stepper), std::move(navigator2));
        Fitter fitter(std::move(propagator), std::move(logger));

        Acts::PropagatorPlainOptions pOptions(geoCtx, magCtx);
       
        Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory> options(
            geoCtx, magCtx, std::ref(calibCtx), extensions,
            pOptions, 
            &initialParams.referenceSurface());
   
        std::vector<const Acts::Surface*> surfaceSequence;
        for (auto idx : indices) {
            const auto& meas = measurements.getMeasurement(idx);
            auto* surf = tGeometry->findSurface(meas.geometryId());
            if (surf) {
                surfaceSequence.push_back(surf);
            }
        }

        options.multipleScattering = true;
        options.energyLoss = true;
        options.referenceSurfaceStrategy = Acts::KalmanFitterTargetSurfaceStrategy::first;

        auto result = fitter.fit(concreteSourceLinks.begin(), concreteSourceLinks.end(), initialParams, options, surfaceSequence, outputTracks);
        if (result.ok()) {
            const auto& track = result.value();
//            std::cout << "=== FIT SUCCESSFUL ===" << std::endl;
//            std::cout << "Measurements added: " << track.nMeasurements() << std::endl;
//            std::cout << "Holes: " << track.nHoles() << std::endl;
//            std::cout << "Chi2: " << track.chi2() << std::endl;
//            std::cout << "Final Q/P: " << track.parameters()[Acts::eBoundQOverP] << std::endl;
//        } else {
//            std::cout << "Fit failed: " << result.error().message() << std::endl;
        }
    }, py::arg("measurements"), py::arg("indices"), py::arg("initialParams"), py::arg("outputTracks"), py::arg("trackingGeometry"), py::arg("magneticField"));

    m.def("makeIndexSourceLink", [](Acts::GeometryIdentifier geoId, std::size_t index) {
        return ActsExamples::IndexSourceLink{geoId, static_cast<unsigned int>(index)};
    });

    py::class_<ActsExamples::IndexSourceLink>(mex, "IndexSourceLink")
        .def(py::init<Acts::GeometryIdentifier, std::size_t>(), py::arg("geoId"), py::arg("index"));

    m.def("createSourceLinks", [](const ActsExamples::MeasurementContainer& measurements,
                                  const std::vector<unsigned int>& indices) {
        std::vector<ActsExamples::IndexSourceLink> sourceLinks;
        sourceLinks.reserve(indices.size());
        for (auto idx : indices) {
            const auto& meas = measurements.getMeasurement(idx);
            sourceLinks.push_back(ActsExamples::IndexSourceLink{meas.geometryId(), idx});
        }
        return sourceLinks;
    });

    py::class_<ActsExamples::VariableMeasurementProxy<6, true>>(m, "MeasurementProxy")
        .def("geometryId", &ActsExamples::VariableMeasurementProxy<6, true>::geometryId)
        .def("parameters", [](const ActsExamples::VariableMeasurementProxy<6, true>& self) {
            return self.parameters();
        })
        .def("covariance", [](const ActsExamples::VariableMeasurementProxy<6, true>& self) {
            return self.covariance();
        });
   
    py::class_<ActsExamples::MeasurementContainer>(m, "MeasurementContainer")
        .def("__len__", &ActsExamples::MeasurementContainer::size)
        .def("__getitem__", [](const ActsExamples::MeasurementContainer& c, size_t i) {
            if (i >= c.size()) throw py::index_error();
            auto it = c.begin();
            std::advance(it, i);
            return *it;
        });

    py::class_<ActsExamples::TrackContainer, std::shared_ptr<ActsExamples::TrackContainer>>(mex, "TrackContainer")
        .def("__iter__", [](const ActsExamples::TrackContainer& c) {
            return py::make_iterator(c.begin(), c.end());
        }, py::keep_alive<0, 1>());
   
    mex.def("makeTrackContainer", []() {
        auto tc = std::make_shared<Acts::VectorTrackContainer>();
        auto mtj = std::make_shared<Acts::VectorMultiTrajectory>();
        return std::make_shared<ActsExamples::TrackContainer>(tc, mtj);
    });
   
    mex.def("getTrackParameters", [](const ActsExamples::TrackContainer& container) {
        std::vector<Acts::BoundTrackParameters> params;
        for (const auto& track : container) {
            if (track.hasReferenceSurface()) {
                params.emplace_back(
                    track.referenceSurface().getSharedPtr(),
                    track.parameters(),
                    track.covariance(),
                    Acts::ParticleHypothesis::muon()
                );
            }
        }
        return params;
    });
   
    m.def("createTargetSurface", [](double z) -> std::shared_ptr<Acts::Surface>{
        auto transform = Acts::Transform3(Acts::Translation3(z, 0.0, 0.0));
        return Acts::Surface::makeShared<Acts::PerigeeSurface>(transform);
    }, py::arg("z"));

    py::class_<Acts::BoundTrackParameters, std::shared_ptr<Acts::BoundTrackParameters>>(mex, "BoundTrackParameters")
        .def("position", [](const Acts::BoundTrackParameters& self, const Acts::GeometryContext& gctx) {
            return self.position(gctx);
        })
        .def("momentum", [](const Acts::BoundTrackParameters& self) {
            return self.momentum();
        })
        .def_property_readonly("parameters", [](const Acts::BoundTrackParameters& self) {
            return self.parameters();
        });

    m.def("createTrackParameters", [](double gx, double gy, double gz,
                                      double px, double py, double pz,
                                      double charge,
                                      std::shared_ptr<const Acts::Surface> surface,
                                      const std::vector<double>& covVec,
                                      const Acts::GeometryContext& gctx) {
        Acts::Vector3 globalPos(gx, gy, gz);
        Acts::Vector3 mom(px, py, pz);

        auto localPosRes = surface->globalToLocal(gctx, globalPos, mom);
        double loc0 = localPosRes.ok() ? localPosRes.value()[0] : 0.0;
        double loc1 = localPosRes.ok() ? localPosRes.value()[1] : 0.0;

        double phi = Acts::VectorHelpers::phi(mom);
        double theta = Acts::VectorHelpers::theta(mom);

        Acts::BoundVector params = Acts::BoundVector::Zero();
        params[Acts::eBoundLoc0] = loc0;
        params[Acts::eBoundLoc1] = loc1;
        params[Acts::eBoundPhi] = phi;
        params[Acts::eBoundTheta] = theta;
        params[Acts::eBoundQOverP] = charge / (mom.norm() + 1e-9);
        params[Acts::eBoundTime] = 0.0;

        Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Identity() * 0.1;
        if (covVec.size() == 36) {
            Eigen::Map<const Acts::BoundSquareMatrix> covMap(covVec.data());
            cov = covMap;
        }

        return Acts::BoundTrackParameters(surface, params, cov, Acts::ParticleHypothesis::muon());
    }, py::arg("x"), py::arg("y"), py::arg("z"),
       py::arg("px"), py::arg("py"), py::arg("pz"),
       py::arg("charge"), py::arg("surface"), py::arg("cov"), py::arg("geo_ctx"));

    m.def("createPlaneSurface", [](Acts::Vector3 center, Acts::Vector3 normal) -> std::shared_ptr<Acts::Surface> {
        Acts::Transform3 transform{Acts::Translation3{center}};
        if (normal.norm() > 1e-6) {
            Acts::Vector3 n = normal.normalized();
            auto rotation = Acts::RotationMatrix3(
                Eigen::Quaternion<double>::FromTwoVectors(Acts::Vector3::UnitZ(), n)
            );
            transform.rotate(rotation);
        }
        return Acts::Surface::makeShared<Acts::PlaneSurface>(transform, nullptr);
    }, py::arg("center"), py::arg("normal"));

    {
        using Builder = ActsExamples::StrawtubeBuilder;
        auto b = py::class_<Builder, ActsExamples::Detector, std::shared_ptr<Builder>>(mex, "StrawtubeBuilder")
                     .def(py::init<const Builder::Config&>(), py::arg("config"))
                     .def("layers", &Builder::layers);
        auto c = py::class_<Builder::Config>(b, "Config").def(py::init<>());
        ACTS_PYTHON_STRUCT(c, fileName, logLevel, layerLogLevel);
    }

    {
        using Detector = ActsExamples::StrawtubeDetector;
        auto d = py::class_<Detector, ActsExamples::Detector, std::shared_ptr<Detector>>(mex, "StrawtubeDetector")
                     .def(py::init<const Detector::Config&>(), py::arg("config"));
        auto c = py::class_<Detector::Config>(d, "Config").def(py::init<>());
        ACTS_PYTHON_STRUCT(c, fileName, logLevel);
    }

    py::class_<HGCBuilder::Config>(mex, "HGCBuilderConfig")
        .def(py::init<>())
        .def_readwrite("fileName", &HGCBuilder::Config::fileName)
        .def_readwrite("logLevel", &HGCBuilder::Config::logLevel);
   
    py::class_<HGCBuilder, Detector, std::shared_ptr<HGCBuilder>>(mex, "HGCBuilder")
        .def(py::init<const HGCBuilder::Config&>())
        .def("layers", &HGCBuilder::layers);
   
    auto hgcDetector = py::class_<HGCDetector, Detector, std::shared_ptr<HGCDetector>>(mex, "HGCDetector")
        .def(py::init<const HGCDetector::Config&>())
        .def("trackingGeometry", &HGCDetector::trackingGeometry);
   
    py::class_<HGCDetector::Config>(hgcDetector, "Config")
        .def(py::init<>())
        .def_readwrite("fileName", &HGCDetector::Config::fileName)
        .def_readwrite("logLevel", &HGCDetector::Config::logLevel);

    m.def("dumpGeometry", [](std::shared_ptr<const Acts::TrackingGeometry> geometry) {
        std::cout << "=== ACTS GEOMETRY DUMP ===" << std::endl;
        geometry->visitSurfaces([](const Acts::Surface* surface) {
            auto geoId = surface->geometryId();
            auto center = surface->center(Acts::GeometryContext());
            std::cout << "Volume: " << geoId.volume()
                      << " | Layer: " << geoId.layer()
                      << " | Sensitive: " << geoId.sensitive()
                      << " | Center: (" << center.x() << ", " << center.y() << ", " << center.z() << ")"
                      << std::endl;
        });
        std::cout << "==========================" << std::endl;
    });
}
}
