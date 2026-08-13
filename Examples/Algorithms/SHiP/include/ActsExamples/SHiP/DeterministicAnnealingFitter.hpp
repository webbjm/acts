#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"

#include <memory>
#include <vector>
#include <unordered_map>

namespace ActsExamples {

/**
 * @class DeterministicAnnealingFitter
 * @brief Deterministic Annealing Filter (DAF) for resolving left/right ambiguity
 *        in drift-tube measurements using iterative Kalman fitting with probabilistic
 *        weight updates inspired by GenFit.
 *
 * Implements a DAF loop that:
 * - Uses annealing schedule to gradually sharpen L/R probability weights
 * - Applies physics-informed priors based on drift radius
 * - Performs iterative Kalman smoothing between weight updates
 * - Converges when max probability change falls below tolerance
 *
 * GenFit-inspired approach for handling left/right ambiguity in wire chambers.
 */
class DeterministicAnnealingFitter {
 public:
  struct Config {
    /// Annealing temperature schedule (default: geometric 100 -> 0.1, 10 steps)
    std::vector<double> annealingSchedule;
    
    /// Maximum iterations in DAF loop (default: 6)
    int maxIterations = 6;
    
    /// Convergence tolerance for max probability change (default: 1e-3)
    double convergenceTolerance = 1e-3;
    
    /// Gate threshold (chi2 gating in sigma units, default: 9.0)
    double gateThreshold = 9.0;
    
    /// Minimum base variance floor (default: 1e-4)
    double minBaseVariance = 1e-4;
    
    /// Variance scaling factor for small drift radii (GenFit-style, default: 0.8)
    double priorVarianceScale = 0.8;
  };

  struct Result {
    bool success = false;
    std::string diagnostics;
  };

  DeterministicAnnealingFitter(Config config, Acts::Logging::Level logLevel = Acts::Logging::INFO);

  /**
   * Fit track using DAF with left/right ambiguity resolution
   *
   * @param measurements Container of drift-tube measurements
   * @param indices Indices of measurements to use
   * @param initialParams Initial track parameters (seed)
   * @param trackingGeometry Geometry for surface lookup
   * @param magneticField Magnetic field for propagation
   * @param outputTracks Output container for fitted track
   * @return Result with success flag and diagnostic message
   */
  Result fit(
      const MeasurementContainer& measurements,
      const std::vector<unsigned int>& indices,
      const Acts::BoundTrackParameters& initialParams,
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
      TrackContainer& outputTracks) const;

 private:
  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace ActsExamples
