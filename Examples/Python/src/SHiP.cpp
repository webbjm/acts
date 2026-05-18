#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "ActsExamples/SHiP/SHiPMeasurementProvider.hpp"
#include "ActsExamples/SHiP/StrawtubeDetector.hpp"
#include "ActsExamples/SHiP/StrawtubeBuilder.hpp"
#include "ActsExamples/DetectorCommons/Detector.hpp"
#include "Acts/Plugins/Python/Utilities.hpp" // Required for Context

namespace py = pybind11;

namespace Acts::Python {
void addSHiP(Context& ctx) {
    auto [m, mex] = ctx.get("main", "examples"); // Get the modules from the context
    using namespace ActsExamples;

//    py::enum_<SHiPDetectorType>(m, "SHiPDetectorType")
//        .value("Straw", SHiPDetectorType::Straw)
//        .value("SciFi", SHiPDetectorType::SciFi)
//        .value("Silicon", SHiPDetectorType::Silicon)
//        .export_values();

    using Alg = SHiPMeasurementProvider;
    py::class_<Alg, IAlgorithm, std::shared_ptr<Alg>>(mex, "SHiPMeasurementProvider")
        .def(py::init<Alg::Config, Acts::Logging::Level>(), py::arg("config"), py::arg("level"))
        .def_property_readonly("config", &Alg::config);

    py::class_<Alg::Config>(mex, "SHiPMeasurementProviderConfig")
        .def(py::init<>())
        .def_readwrite("inputHitArray", &Alg::Config::inputHitArray)
        .def_readwrite("outputMeasurements", &Alg::Config::outputMeasurements)
//        .def_readwrite("detectorType", &Alg::Config::detectorType)
        .def_readwrite("trackingGeometry", &Alg::Config::trackingGeometry)
        .def_readwrite("strawRes", &Alg::Config::strawRes)
        .def_readwrite("scifiPitch", &Alg::Config::scifiPitch)
        .def_readwrite("siliconPitch", &Alg::Config::siliconPitch)
        .def_readwrite("scifiRes", &Alg::Config::scifiRes)
        .def_readwrite("siliconRes", &Alg::Config::siliconRes);

     // Bind StrawtubeBuilder
     {
         using Builder = ActsExamples::StrawtubeBuilder;

         auto b = py::class_<Builder, ActsExamples::Detector, std::shared_ptr<Builder>>(mex, "StrawtubeBuilder")
                      .def(py::init<const Builder::Config&>(), py::arg("config"))
                      .def("layers", &Builder::layers);

         auto c = py::class_<Builder::Config>(b, "Config").def(py::init<>());
         ACTS_PYTHON_STRUCT(c, fileName, logLevel, layerLogLevel);
     }

     // Bind StrawtubeDetector
     {
         using Detector = ActsExamples::StrawtubeDetector;

         auto d = py::class_<Detector, ActsExamples::Detector, std::shared_ptr<Detector>>(mex, "StrawtubeDetector")
                      .def(py::init<const Detector::Config&>(), py::arg("config"));

         auto c = py::class_<Detector::Config>(d, "Config").def(py::init<>());
         ACTS_PYTHON_STRUCT(c, fileName, logLevel);
     }
}
}
