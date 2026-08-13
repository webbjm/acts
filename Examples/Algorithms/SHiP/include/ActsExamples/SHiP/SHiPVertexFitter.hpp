#ifndef ACTSEXAMPLES_SHIP_VERTEXFITTER_HPP
#define ACTSEXAMPLES_SHIP_VERTEXFITTER_HPP

#include <memory>
#include <vector>

#include "ActsExamples/EventData/Track.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Definitions/Algebra.hpp"

namespace ActsExamples {

class SHiPVertexFitter {
  public:
    struct Config {
        size_t maxIterations = 25;
        double perigeeInflation = 10.0; // factor to inflate input perigee spatial errors
        double timeInflation = 1e6;     // inflate time variance to de-weight time
        double seedXMin = 20000.0;
        double seedXMax = 90000.0;
    };

    SHiPVertexFitter(const Config& cfg, std::shared_ptr<const Acts::MagneticFieldProvider> bField);

    // Fit vertices from track proxies; returns a vector of Acts::Vertex (empty if failed)
    std::vector<Acts::Vertex> fit(const std::vector<ActsExamples::ConstTrackContainer::ConstTrackProxy>& proxies,
                                   const Acts::GeometryContext& geoCtx,
                                   const Acts::TrackingGeometry& trackingGeometry) const;

  private:
    Config m_cfg;
    std::shared_ptr<const Acts::MagneticFieldProvider> m_bField;
};

} // namespace ActsExamples

#endif // ACTSEXAMPLES_SHIP_VERTEXFITTER_HPP
