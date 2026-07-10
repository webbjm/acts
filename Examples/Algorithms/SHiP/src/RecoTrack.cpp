// -------------------------------------------------------------------------
// -----                   RecoTrack source file                   -----
#include "ActsExamples/SHiP/RecoTrack.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include <TBuffer.h>
#include <TClass.h>
#include <TROOT.h>
#include <algorithm>
#include <cmath>


ClassImp(ActsExamples::RecoTrack);

namespace ActsExamples {

RecoTrack::RecoTrack() : TObject() {}

RecoTrack::~RecoTrack() = default;

RecoTrack::RecoTrack(const Acts::BoundVector& parameters,
                     const Acts::BoundSquareMatrix& covariance,
                     const Acts::Surface& surface,
                     float chi2,
                     unsigned int nDoF,
                     const Acts::GeometryContext& gctx,
                     const std::vector<Double_t>& residuals,
                     const std::vector<Double_t>& pulls,
                     unsigned int nMeasurements
                     )
    : TObject(),
      m_chi2(chi2),
      m_nDoF(nDoF),
      m_residuals(residuals),
      m_nMeasurements(nMeasurements),
      m_pulls(pulls) {

  const auto freeParams = Acts::transformBoundToFreeParameters(surface, gctx, parameters);

  m_x = -freeParams[Acts::eFreePos2] * 0.1; //scale to cm
  m_y =  freeParams[Acts::eFreePos1] * 0.1;
  m_z =  freeParams[Acts::eFreePos0] * 0.1;

  const Double_t invQop = parameters[Acts::eBoundQOverP];
  const Double_t p = (invQop != 0.0) ? std::abs(1.0 / invQop) : 0.0;

  m_px = -freeParams[Acts::eFreeDir2] * p;
  m_py =  freeParams[Acts::eFreeDir1] * p;
  m_pz =  freeParams[Acts::eFreeDir0] * p;

  const auto jacobian = surface.boundToFreeJacobian(
      gctx,
      freeParams.template segment<3>(Acts::eFreePos0),
      freeParams.template segment<3>(Acts::eFreeDir0)
  );
  const auto freeCov = jacobian * covariance * jacobian.transpose();

  m_err_x = std::sqrt(std::max(0.0, covariance(Acts::eBoundLoc0, Acts::eBoundLoc0))) * 0.1;
  m_err_y = std::sqrt(std::max(0.0, covariance(Acts::eBoundLoc1, Acts::eBoundLoc1))) * 0.1;
  m_err_z = 0.0;

  Eigen::Matrix<double, 6, 6> actsSpatialCov = Eigen::Matrix<double, 6, 6>::Zero();
  actsSpatialCov.block<3,3>(0,0) = freeCov.block<3,3>(Acts::eFreePos0, Acts::eFreePos0) * 0.1;
  actsSpatialCov.block<3,3>(3,3) = freeCov.block<3,3>(Acts::eFreeDir0, Acts::eFreeDir0) * 0.1;
  actsSpatialCov.block<3,3>(0,3) = freeCov.block<3,3>(Acts::eFreePos0, Acts::eFreeDir0) * 0.1;
  actsSpatialCov.block<3,3>(3,0) = freeCov.block<3,3>(Acts::eFreeDir0, Acts::eFreePos0) * 0.1;

  Eigen::Matrix<double, 6, 6> R = Eigen::Matrix<double, 6, 6>::Zero();
  R(0, 2) = -1.0;
  R(1, 1) =  1.0;
  R(2, 0) =  1.0;
  R(3, 5) = -1.0;
  R(4, 4) =  1.0;
  R(5, 3) =  1.0;

  const Eigen::Matrix<double, 6, 6> shipSpatialCov = R * actsSpatialCov * R.transpose();

  m_err_px = std::sqrt(std::max(0.0, shipSpatialCov(3, 3)));
  m_err_py = std::sqrt(std::max(0.0, shipSpatialCov(4, 4)));
  m_err_pz = std::sqrt(std::max(0.0, shipSpatialCov(5, 5)));

  m_cov_elements.clear();
  m_cov_elements.reserve(21);
  for (int row = 0; row < 6; ++row) {
      for (int col = row; col < 6; ++col) {
          m_cov_elements.push_back(shipSpatialCov(row, col));
      }
  }
}

} // namespace ActsExamples

