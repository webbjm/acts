#include "ActsExamples/SHiP/DeterministicAnnealingFitter.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/Utilities/AnnealingUtility.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <sstream>
#include <algorithm>
#include <cmath>

namespace ActsExamples::SHiP {

DeterministicAnnealingFitter::DeterministicAnnealingFitter(
    Config config, Acts::Logging::Level logLevel)
    : m_cfg(std::move(config)),
      m_logger(Acts::getDefaultLogger("DeterministicAnnealingFitter", logLevel)) {
  // Set default annealing schedule if not provided
  if (m_cfg.annealingSchedule.empty()) {
    double bStart = 100.0, bFinal = 0.1;
    unsigned int nSteps = 10;
    for (unsigned int i = 0; i < nSteps; ++i) {
      m_cfg.annealingSchedule.push_back(
          bStart * std::pow(bFinal / bStart, double(i) / (nSteps - 1)));
    }
  }
}

DeterministicAnnealingFitter::Result DeterministicAnnealingFitter::fit(
    const MeasurementContainer& measurements,
    const std::vector<unsigned int>& indices,
    const Acts::BoundTrackParameters& initialParams,
    std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> bField,
    TrackContainer& outputTracks) const {
  
  if (m_logger) {
    ACTS_LOG_WITH_LOGGER(*m_logger, Acts::Logging::DEBUG,
                         "DeterministicAnnealingFitter::fit called with " << indices.size() << " measurements");
  }
  
  Result result;
  result.success = false;

  Acts::GeometryContext geoCtx = Acts::GeometryContext::dangerouslyDefaultConstruct();
  Acts::MagneticFieldContext magCtx;
  Acts::CalibrationContext calibCtx;

  // Collect source links and surfaces
  std::vector<Acts::SourceLink> concreteSourceLinks;
  std::vector<const Acts::Surface*> surfaceSequence;
  std::unordered_map<unsigned int, size_t> measIndexToPos;

  for (size_t p = 0; p < indices.size(); ++p) {
    unsigned int idx = indices[p];
    measIndexToPos[idx] = p;
    try {
      const auto& meas = measurements.getMeasurement(idx);
      auto* surface = tGeometry->findSurface(meas.geometryId());
      if (surface) {
        concreteSourceLinks.push_back(
            Acts::SourceLink(IndexSourceLink{meas.geometryId(), idx}));
        surfaceSequence.push_back(surface);
      }
    } catch (...) {
    }
  }

  if (concreteSourceLinks.empty()) {
    result.diagnostics = "[DAF] No valid source links";
    return result;
  }

  // Compute GenFit-style physics-informed priors based on drift radius
  std::vector<double> priorWeight(indices.size(), 0.0);
  {
    double maxR = 0.0;
    for (size_t p = 0; p < indices.size(); ++p) {
      try {
        const auto& meas = measurements.getMeasurement(indices[p]);
        double r = std::abs(meas.parameters()(0));
        if (r > maxR) maxR = r;
      } catch (...) {
      }
    }
    if (maxR <= 0.0) maxR = 1.0;

    for (size_t p = 0; p < indices.size(); ++p) {
      try {
        const auto& meas = measurements.getMeasurement(indices[p]);
        double r = std::abs(meas.parameters()(0));
        double frac = std::max(0.0, 1.0 - r / maxR);
        priorWeight[p] = 0.5 * frac * frac;
      } catch (...) {
        priorWeight[p] = 0.0;
      }
    }
  }

  // Initialize probRight uniformly (will be refined by DAF)
  std::vector<double> probRight(indices.size(), 0.5);

  // Compute base variances
  double minBaseVar = m_cfg.minBaseVariance;
  {
    std::vector<double> measuredVars;
    for (auto idx : indices) {
      try {
        const auto& meas = measurements.getMeasurement(idx);
        measuredVars.push_back(meas.covariance()(0, 0));
      } catch (...) {
        measuredVars.push_back(m_cfg.minBaseVariance);
      }
    }
    if (!measuredVars.empty()) {
      std::sort(measuredVars.begin(), measuredVars.end());
      double median = measuredVars[measuredVars.size() / 2];
      minBaseVar = std::max(m_cfg.minBaseVariance, 0.1 * median);
      minBaseVar = std::max(minBaseVar, 1e-6);
    }
  }

  // Setup annealing utility
  Acts::AnnealingUtility::Config aCfg(m_cfg.gateThreshold, m_cfg.annealingSchedule);
  Acts::AnnealingUtility annealer(aCfg);
  Acts::AnnealingUtility::State aState;

  // DAF main loop
  bool converged = false;
  std::vector<double> lastGoodProbRight;
  bool haveLastGood = false;

  for (int iter = 0; iter < m_cfg.maxIterations; ++iter) {
    // Setup Kalman fitter for this iteration
    Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions_iter;

    auto accessor = [tg = tGeometry.get()](const Acts::SourceLink& sl) -> const Acts::Surface* {
      auto geoId = sl.template get<IndexSourceLink>().geometryId();
      const auto* surf = tg->findSurface(geoId);
      if (surf) const_cast<Acts::Surface*>(surf)->assignGeometryId(geoId);
      return surf;
    };

    // Calibrator: mix left/right measurements per probRight
    auto calibrator_iter = [&measurements, &measIndexToPos, &probRight,
                            &priorWeight, minBaseVar,
                            priorVarScale = m_cfg.priorVarianceScale](
        const Acts::GeometryContext&, const Acts::CalibrationContext&,
        const Acts::SourceLink& sl,
        Acts::TrackStateProxy<Acts::VectorMultiTrajectory, 6, false> ts) {
      const auto& islink = sl.template get<IndexSourceLink>();
      const auto& meas = measurements.getMeasurement(islink.index());
      auto measDim = meas.size();
      if (measDim == 1) {
        double r = meas.parameters()(0);
        double pRight = 0.5;
        auto it = measIndexToPos.find(islink.index());
        size_t pos = (it != measIndexToPos.end()) ? it->second : 0;
        if (it != measIndexToPos.end()) pRight = probRight[it->second];

        double pLeft = 1.0 - pRight;
        double left = -r;
        double right = r;
        double mean = pLeft * left + pRight * right;

        ts.allocateCalibrated(1);
        ts.template calibrated<1>()(0) = mean;

        double baseVar = 1e-6;
        try {
          baseVar = meas.covariance()(0, 0);
        } catch (...) {
        }
        baseVar = std::max(baseVar, minBaseVar);

        // Scale variance down for small drift radii (GenFit-style prior)
        double prior = 0.0;
        if (pos < priorWeight.size()) prior = priorWeight[pos];
        double priorNorm = std::min(1.0, prior / 0.5);
        double scale = 1.0 - priorVarScale * priorNorm;
        scale = std::max(scale, 0.1);
        baseVar *= scale;

        double mixVar = pLeft * (left - mean) * (left - mean) +
                        pRight * (right - mean) * (right - mean);
        ts.template calibratedCovariance<1>()(0, 0) = baseVar + mixVar;
      } else if (measDim == 2) {
        ts.allocateCalibrated(2);
        ts.template calibrated<2>() = meas.parameters();
        ts.template calibratedCovariance<2>() = meas.covariance();
      }
      Acts::SourceLink slCopy = sl;
      ts.setUncalibratedSourceLink(std::move(slCopy));
    };

    extensions_iter.surfaceAccessor.connect(accessor);
    extensions_iter.calibrator.connect(calibrator_iter);

    Acts::GainMatrixUpdater updater;
    Acts::GainMatrixSmoother smoother;
    extensions_iter.updater.template connect<
        &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(&updater);
    extensions_iter.smoother.template connect<
        &Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(&smoother);

    auto outlierFinder = [](Acts::TrackStateProxy<Acts::VectorMultiTrajectory, 6, true>) -> bool {
      return false;
    };
    extensions_iter.outlierFinder.connect(outlierFinder);

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
        geoCtx, magCtx, std::ref(calibCtx), extensions_iter, pOptions,
        &initialParams.referenceSurface());

    options.multipleScattering = true;
    options.energyLoss = true;
    options.referenceSurfaceStrategy = Acts::TrackExtrapolationStrategy::first;

    auto trackBackend = std::make_shared<Acts::VectorTrackContainer>();
    auto trajectoryBackend = std::make_shared<Acts::VectorMultiTrajectory>();
    TrackContainer tempOutput(trackBackend, trajectoryBackend);

    auto result_iter = fitter.fit(concreteSourceLinks.begin(), concreteSourceLinks.end(),
                                   initialParams, options, surfaceSequence, tempOutput);

    if (!result_iter.ok()) {
      if (haveLastGood) {
        probRight = lastGoodProbRight;
        converged = true;
        result.success = true;
        result.diagnostics = "[DAF] Iteration fit failed, using last-good state";
        break;
      } else {
        result.diagnostics = "[DAF] Initial fit failed";
        return result;
      }
    }

    if (tempOutput.size() == 0) {
      if (haveLastGood) {
        probRight = lastGoodProbRight;
        converged = true;
        result.success = true;
        result.diagnostics = "[DAF] Iteration produced empty output, using last-good state";
        break;
      } else {
        result.diagnostics = "[DAF] Initial iteration produced empty output";
        return result;
      }
    }

    // Extract predicted positions from fitted track
    const auto& fittedProxy = tempOutput.getTrack(0);
    std::unordered_map<unsigned int, double> predictedMap;
    for (const auto& state : fittedProxy.trackStatesReversed()) {
      if (!state.hasUncalibratedSourceLink()) continue;
      auto sl = state.getUncalibratedSourceLink();
      unsigned int measIndex = sl.get<IndexSourceLink>().index();
      double pred = 0.0;
      if (state.hasSmoothed())
        pred = state.smoothed()[Acts::eBoundLoc0];
      else if (state.hasPredicted())
        pred = state.predicted()[Acts::eBoundLoc0];
      predictedMap[measIndex] = pred;
    }

    // Update probabilities using annealing
    double maxChange = 0.0;
    const double epsP = 1e-6;
    for (size_t p = 0; p < indices.size(); ++p) {
      unsigned int measIdx = indices[p];
      double r_abs = 0.0;
      double baseVar = 1e-6;
      try {
        const auto& meas = measurements.getMeasurement(measIdx);
        r_abs = meas.parameters()(0);
        baseVar = meas.covariance()(0, 0);
      } catch (...) {
      }
      baseVar = std::max(baseVar, minBaseVar);

      double pred = 0.0;
      auto itp = predictedMap.find(measIdx);
      if (itp != predictedMap.end()) pred = itp->second;

      double pRight_curr = probRight[p];
      pRight_curr = std::max(epsP, std::min(pRight_curr, 1.0 - epsP));
      double pLeft_curr = 1.0 - pRight_curr;

      double left = -r_abs;
      double right = r_abs;
      double mean = pLeft_curr * left + pRight_curr * right;
      double mixVar = pLeft_curr * (left - mean) * (left - mean) +
                      pRight_curr * (right - mean) * (right - mean);
      double totalVar = baseVar + mixVar;
      totalVar = std::max(totalVar, minBaseVar);

      double rLeft = pred - left;
      double rRight = pred - right;
      double gateThresh = m_cfg.gateThreshold * std::sqrt(totalVar);

      double chi2L = (totalVar > 0) ? (rLeft * rLeft / totalVar) : 1e300;
      double chi2R = (totalVar > 0) ? (rRight * rRight / totalVar) : 1e300;

      double wL = 0.0, wR = 0.0;
      if (std::abs(rLeft) <= gateThresh) wL = annealer.getWeight(aState, chi2L);
      if (std::abs(rRight) <= gateThresh) wR = annealer.getWeight(aState, chi2R);

      double newProbR = probRight[p];
      double sumw = wL + wR;
      if (sumw > 0.0) newProbR = wR / sumw;

      maxChange = std::max(maxChange, std::abs(newProbR - probRight[p]));
      probRight[p] = newProbR;
    }

    // Advance annealing
    annealer.anneal(aState);

    lastGoodProbRight = probRight;
    haveLastGood = true;

    // Check convergence
    if (maxChange < m_cfg.convergenceTolerance || aState.equilibriumReached) {
      converged = true;
      result.success = true;
      result.diagnostics = "[DAF] Converged";
      break;
    }
  }

  // Final fit with converged probRight
  if (!converged && haveLastGood) {
    probRight = lastGoodProbRight;
    result.success = true;
    result.diagnostics = "[DAF] Max iterations reached, using last-good state";
  }

  if (result.success) {
    // Run final fit with optimized probRight
    Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions_final;

    auto accessor_f = [tg = tGeometry.get()](const Acts::SourceLink& sl) -> const Acts::Surface* {
      auto geoId = sl.template get<IndexSourceLink>().geometryId();
      const auto* surf = tg->findSurface(geoId);
      if (surf) const_cast<Acts::Surface*>(surf)->assignGeometryId(geoId);
      return surf;
    };

    auto calibrator_final = [&measurements, &measIndexToPos, &probRight, &priorWeight,
                             minBaseVar, priorVarScale = m_cfg.priorVarianceScale](
        const Acts::GeometryContext&, const Acts::CalibrationContext&,
        const Acts::SourceLink& sl,
        Acts::TrackStateProxy<Acts::VectorMultiTrajectory, 6, false> ts) {
      const auto& islink = sl.template get<IndexSourceLink>();
      const auto& meas = measurements.getMeasurement(islink.index());
      auto measDim = meas.size();
      if (measDim == 1) {
        double r = meas.parameters()(0);
        double pRight = 0.5;
        auto it = measIndexToPos.find(islink.index());
        size_t pos = (it != measIndexToPos.end()) ? it->second : 0;
        if (it != measIndexToPos.end()) pRight = probRight[it->second];

        double pLeft = 1.0 - pRight;
        double left = -r;
        double right = r;
        double mean = pLeft * left + pRight * right;

        ts.allocateCalibrated(1);
        ts.template calibrated<1>()(0) = mean;

        double baseVar = 1e-6;
        try {
          baseVar = meas.covariance()(0, 0);
        } catch (...) {
        }
        baseVar = std::max(baseVar, minBaseVar);

        // Apply prior-based variance scaling
        double prior = 0.0;
        if (pos < priorWeight.size()) prior = priorWeight[pos];
        double priorNorm = std::min(1.0, prior / 0.5);
        double scale = 1.0 - priorVarScale * priorNorm;
        scale = std::max(scale, 0.1);
        baseVar *= scale;

        double mixVar = pLeft * (left - mean) * (left - mean) +
                        pRight * (right - mean) * (right - mean);
        ts.template calibratedCovariance<1>()(0, 0) = baseVar + mixVar;
      } else if (measDim == 2) {
        ts.allocateCalibrated(2);
        ts.template calibrated<2>() = meas.parameters();
        ts.template calibratedCovariance<2>() = meas.covariance();
      }
      Acts::SourceLink slCopy = sl;
      ts.setUncalibratedSourceLink(std::move(slCopy));
    };

    extensions_final.surfaceAccessor.connect(accessor_f);
    extensions_final.calibrator.connect(calibrator_final);

    Acts::GainMatrixUpdater updater_f;
    Acts::GainMatrixSmoother smoother_f;
    extensions_final.updater.template connect<
        &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(&updater_f);
    extensions_final.smoother.template connect<
        &Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(&smoother_f);

    auto outlierFinderFinal = [](Acts::TrackStateProxy<Acts::VectorMultiTrajectory, 6, true>) -> bool {
      return false;
    };
    extensions_final.outlierFinder.connect(outlierFinderFinal);

    Acts::Navigator::Config navCfg_f{tGeometry};
    navCfg_f.resolvePassive = true;
    navCfg_f.resolveMaterial = true;
    navCfg_f.resolveSensitive = true;
    Acts::Navigator navigator_f(navCfg_f, Acts::getDefaultLogger("Navigator", Acts::Logging::INFO));
    auto logger_f = Acts::getDefaultLogger("KalmanFitter", Acts::Logging::INFO);
    Acts::DirectNavigator navigator2f;

    using Stepper = Acts::EigenStepper<>;
    using Propagator = Acts::Propagator<Stepper, Acts::DirectNavigator>;
    using Fitter = Acts::KalmanFitter<Propagator, Acts::VectorMultiTrajectory>;

    Stepper stepper_f(bField);
    Propagator propagator_f(std::move(stepper_f), std::move(navigator2f));
    Fitter fitter_f(std::move(propagator_f), std::move(logger_f));

    Acts::PropagatorPlainOptions pOptions_f(geoCtx, magCtx);
    Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory> options_f(
        geoCtx, magCtx, std::ref(calibCtx), extensions_final, pOptions_f,
        &initialParams.referenceSurface());

    options_f.multipleScattering = true;
    options_f.energyLoss = true;
    options_f.referenceSurfaceStrategy = Acts::TrackExtrapolationStrategy::first;

    auto result_final = fitter_f.fit(concreteSourceLinks.begin(), concreteSourceLinks.end(),
                                      initialParams, options_f, surfaceSequence, outputTracks);

    if (!result_final.ok()) {
      result.success = false;
      result.diagnostics += "; final fit failed";
    }
  }

  return result;
}

}  // namespace ActsExamples::SHiP
