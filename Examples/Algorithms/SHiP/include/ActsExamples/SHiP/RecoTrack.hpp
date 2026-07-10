#ifndef ACTSEXAMPLES_SHIP_RECOTRACK_HPP
#define ACTSEXAMPLES_SHIP_RECOTRACK_HPP

#include "TObject.h"
#include <vector>
#include <cmath>
#include <algorithm>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace ActsExamples {

class RecoTrack : public TObject {
 public:
  RecoTrack();

  RecoTrack(const Acts::BoundVector& parameters,
            const Acts::BoundSquareMatrix& covariance,
            const Acts::Surface& surface,
            float chi2,
            unsigned int nDoF,
            const Acts::GeometryContext& gctx,
            const std::vector<Double_t>& residuals,
            const std::vector<Double_t>& pulls,
            unsigned int nMeasurements
            );

  ~RecoTrack() override;

  template <typename TrackProxyType>
  static RecoTrack FromActsProxy(const TrackProxyType& track, const Acts::GeometryContext& gctx) {
      std::vector<Double_t> residuals;
      std::vector<Double_t> pulls;
      unsigned int hit_count = 0;

      for (const auto& state : track.trackStates()) {
          if (!state.hasUncalibratedSourceLink() || !state.hasCalibrated() || !state.hasSmoothed()) {
              continue;
          }

          auto meas_params = state.template calibrated<1>();
          auto meas_cov = state.template calibratedCovariance<1>();
          auto smoothed_params = state.parameters();
          auto smoothed_cov = state.covariance();

          double res_val = meas_params(0) - smoothed_params(0);
          double res_cov = meas_cov(0, 0) - smoothed_cov(0, 0);

          residuals.push_back(res_val);
          if (res_cov > 0) {
              pulls.push_back(res_val / std::sqrt(res_cov));
          } else {
              pulls.push_back(-999.0);
          }
          hit_count++;
      }

      return RecoTrack(
          track.parameters(),
          track.covariance(),
          track.referenceSurface(),
          track.chi2(),
          track.nDoF(),
          gctx,
          residuals,
          pulls,
          hit_count
      );
  }

  // Getters
  Double_t x() const { return m_x; }
  Double_t y() const { return m_y; }
  Double_t z() const { return m_z; }
  Double_t px() const { return m_px; }
  Double_t py() const { return m_py; }
  Double_t pz() const { return m_pz; }

  Double_t err_x() const { return m_err_x; }
  Double_t err_y() const { return m_err_y; }
  Double_t err_z() const { return m_err_z; }
  Double_t err_px() const { return m_err_px; }
  Double_t err_py() const { return m_err_py; }
  Double_t err_pz() const { return m_err_pz; }

  Double_t chi2() const { return m_chi2; }
  Double_t nDoF() const { return m_nDoF; }
  unsigned int nMeasurements() const { return m_nMeasurements; }

  const std::vector<Double_t>& GetResiduals() const { return m_residuals; }
  const std::vector<Double_t>& GetPulls() const { return m_pulls; }
  const std::vector<Double_t>& GetCovarianceElements() const { return m_cov_elements; }



 private:
  Double_t m_x{0.};
  Double_t m_y{0.};
  Double_t m_z{0.};
  Double_t m_px{0.};
  Double_t m_py{0.};
  Double_t m_pz{0.};

  Double_t m_err_x{0.};
  Double_t m_err_y{0.};
  Double_t m_err_z{0.};
  Double_t m_err_px{0.};
  Double_t m_err_py{0.};
  Double_t m_err_pz{0.};

  Double_t m_chi2{0.};
  Double_t m_nDoF{0.};

  std::vector<Double_t> m_residuals;
  unsigned int m_nMeasurements{0};
  std::vector<Double_t> m_pulls;
  std::vector<Double_t> m_cov_elements;

  ClassDef(RecoTrack, 1);
};
} // namespace ActsExamples

#endif // ACTSEXAMPLES_SHIP_RECOTRACK_HPP

