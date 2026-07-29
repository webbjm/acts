#pragma once

#include <vector>
#include <memory>
#include <stdexcept>

// Include core ACTS track types needed by pushRecoTrack
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/Surface.hpp"

namespace ActsExamples {

/**
 * @brief A light, python-friendly proxy wrapper representing a single fitted track.
 * Safely mimics the property names expected by pushRecoTrack.
 */
struct ShipTrackProxy {
    std::shared_ptr<const Acts::Surface> referenceSurface{nullptr};
    Acts::BoundVector parameters{Acts::BoundVector::Zero()};
    Acts::BoundMatrix covariance{Acts::BoundMatrix::Zero()};
    unsigned int nMeasurements{0};
    unsigned int nHoles{0};
    float chi2{0.0f};

    // Helper functions to mimic python attributes
    bool hasReferenceSurface() const { return referenceSurface != nullptr; }
    unsigned int nDoF() const { return (nMeasurements > 5) ? (nMeasurements - 5) : 0; }
};

/**
 * @brief An iterable, indexable collection of ShipTrackProxy items.
 */
class ShipTrackContainer {
public:
    ShipTrackContainer() = default;
    ~ShipTrackContainer() = default;

    // Vector modifiers
    void append(const ShipTrackProxy& track) { m_tracks.push_back(track); }
    void clear() { m_tracks.clear(); }
    
    // Size and bounds checking
    size_t size() const { return m_tracks.size(); }
    bool empty() const { return m_tracks.empty(); }

    // Index-based accessor
    ShipTrackProxy& getTrack(size_t index) {
        if (index >= m_tracks.size()) {
            throw std::out_of_range("ShipTrackContainer: Index out of range");
        }
        return m_tracks[index];
    }

    const ShipTrackProxy& getTrack(size_t index) const {
        if (index >= m_tracks.size()) {
            throw std::out_of_range("ShipTrackContainer: Index out of range");
        }
        return m_tracks[index];
    }

    // Iteration support for C++ STL (and pybind11)
    auto begin() { return m_tracks.begin(); }
    auto end() { return m_tracks.end(); }
    auto begin() const { return m_tracks.begin(); }
    auto end() const { return m_tracks.end(); }

private:
    std::vector<ShipTrackProxy> m_tracks;
};

} // namespace ActsExamples

