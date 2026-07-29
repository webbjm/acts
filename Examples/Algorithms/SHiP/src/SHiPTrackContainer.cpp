#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "ShipTrackContainer.hpp"

namespace py = pybind11;
using namespace ActsExamples;

// This macro registers your classes into your compiled Acts Examples Python Bindings module namespace
PYBIND11_MODULE(ShipTrackContainerModule, m) {
    m.doc() = "Custom tracking container extension module for FairShip ACTS upgrades.";

    // 1. Bind ShipTrackProxy
    py::class_<ShipTrackProxy>(m, "ShipTrackProxy")
        .def(py::init<>())
        // Explicitly map properties as read/write to match pushRecoTrack lookups
        .def_property("referenceSurface", 
            [](const ShipTrackProxy& self) { return self.referenceSurface; },
            [](ShipTrackProxy& self, std::shared_ptr<const Acts::Surface> srf) { self.referenceSurface = std::move(srf); })
        
        .def_property("parameters",
            [](const ShipTrackProxy& self) { return self.parameters; },
            [](ShipTrackProxy& self, const Acts::BoundVector& v) { self.parameters = v; })
        
        .def_property("covariance",
            [](const ShipTrackProxy& self) { return self.covariance; },
            [](ShipTrackProxy& self, const Acts::BoundMatrix& m) { self.covariance = m; })
        
        .def_readwrite("nMeasurements", &ShipTrackProxy::nMeasurements)
        .def_readwrite("nHoles", &ShipTrackProxy::nHoles)
        .def_readwrite("chi2", &ShipTrackProxy::chi2)

        // Read-only synthetic properties requested by pushRecoTrack lookups
        .def_property_readonly("hasReferenceSurface", &ShipTrackProxy::hasReferenceSurface)
        .def_property_readonly("nDoF", &ShipTrackProxy::nDoF)
        
        // Return an empty list for trackStatesReversed to bypass residuals lookup loop and trigger fallback
        .def_property_readonly("trackStatesReversed", [](const ShipTrackProxy&) {
            return std::vector<py::object>();
        });

    // 2. Bind ShipTrackContainer
    py::class_<ShipTrackContainer>(m, "ShipTrackContainer")
        .def(py::init<>())
        .def("append", &ShipTrackContainer::append, py::arg("track"))
        .def("clear", &ShipTrackContainer::clear)
        .def("size", &ShipTrackContainer::size)
        .def("empty", &ShipTrackContainer::empty)
        .def("getTrack", py::overload_cast<size_t>(&ShipTrackContainer::getTrack), py::return_value_policy::reference_internal)
        
        // Native Python length hook: len(container)
        .def("__len__", &ShipTrackContainer::size)
        
        // Native Python indexing hook: container[i]
        .def("__getitem__", [](ShipTrackContainer& self, size_t index) {
            if (index >= self.size()) {
                throw py::index_error("Index out of range");
            }
            return self.getTrack(index);
        }, py::return_value_policy::reference_internal)

        // Native Python iteration hook: for track in container:
        .def("__iter__", [](ShipTrackContainer& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());
}

