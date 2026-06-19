// -------------------------------------------------------------------------
// -----                   RecoTrack source file                   -----
#include "ActsExamples/SHiP/RecoTrack.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include <TBuffer.h>
#include <TClass.h>
#include <TROOT.h>

namespace ActsExamples {

ClassImp(RecoTrack);

RecoTrack::RecoTrack() : TObject() {}

RecoTrack::RecoTrack(const ConstTrackProxy& track, const Acts::GeometryContext& gctx)
    : TObject() {

  if (!track.hasReferenceSurface()) {
      return;
  }

  const auto& parameters = track.parameters();
  const auto& surface = track.referenceSurface();

  // Transform bound to free parameters for rotation
  auto freeParams = Acts::transformBoundToFreeParameters(surface, gctx, parameters);

  // SHiP Rotation: x=-ActsZ, y=ActsY, z=ActsX
  m_x = -freeParams[Acts::eFreePos2];
  m_y = freeParams[Acts::eFreePos1];
  m_z = freeParams[Acts::eFreePos0];

  float p = std::abs(1.0f / parameters[Acts::eBoundQOverP]);
  m_px = -freeParams[Acts::eFreeDir2] * p;
  m_py = freeParams[Acts::eFreeDir1] * p;
  m_pz = freeParams[Acts::eFreeDir0] * p;

  // Calculate errors using the Jacobian
  const auto& covariance = track.covariance();
  auto jacobian = surface.boundToFreeJacobian(gctx, freeParams.segment<3>(Acts::eFreePos0), freeParams.segment<3>(Acts::eFreeDir0));
  auto freeCov = jacobian * covariance * jacobian.transpose();
    
  m_err_x = std::sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0));
  m_err_y = std::sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1));
  m_err_z = 0.0; // Defined by surface position
    
  auto qop = freeParams[Acts::eFreeQOverP];
  auto calcMomentumErr = [&](double mom_val, int dir_idx) {
      double dir = freeParams[dir_idx];
      return std::abs(mom_val) * std::sqrt(std::abs(
          freeCov(dir_idx, dir_idx) / (dir * dir) + 
          freeCov(Acts::eFreeQOverP, Acts::eFreeQOverP) / (qop * qop) + 
          2 * freeCov(dir_idx, Acts::eFreeQOverP) / (dir * qop)));
  };

  m_err_px = calcMomentumErr(m_px, Acts::eFreeDir2);
  m_err_py = calcMomentumErr(m_py, Acts::eFreeDir1);
  m_err_pz = calcMomentumErr(m_pz, Acts::eFreeDir0);
  
  // Set Chi2 and nDoF from the track proxy
  m_chi2 = track.chi2();
  m_nDoF = track.nDoF();
} 

RecoTrack::RecoTrack(const Acts::BoundVector& parameters,
                    const Acts::BoundSquareMatrix& covariance,
                    const Acts::Surface& surface,
                    float chi2,
                    unsigned int nDoF,
                    const Acts::GeometryContext& gctx)
    : TObject(), m_chi2(chi2), m_nDoF(nDoF) {

  auto freeParams = Acts::transformBoundToFreeParameters(surface, gctx, parameters);
  
  m_x = -freeParams[Acts::eFreePos2];
  m_y = freeParams[Acts::eFreePos1];
  m_z = freeParams[Acts::eFreePos0];

  float p = std::abs(1.0f / parameters[Acts::eBoundQOverP]);
  m_px = -freeParams[Acts::eFreeDir2] * p;
  m_py = freeParams[Acts::eFreeDir1] * p;
  m_pz = freeParams[Acts::eFreeDir0] * p;

  auto jacobian = surface.boundToFreeJacobian(gctx, freeParams.segment<3>(Acts::eFreePos0), freeParams.segment<3>(Acts::eFreeDir0));
  auto freeCov = jacobian * covariance * jacobian.transpose();

  m_err_x = std::sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0));
  m_err_y = std::sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1));
  m_err_z = 0.0;

  auto qop = freeParams[Acts::eFreeQOverP];
  auto calcMomentumErr = [&](double mom_val, int dir_idx) {
      double dir = freeParams[dir_idx];
      return std::abs(mom_val) * std::sqrt(std::abs(
          freeCov(dir_idx, dir_idx) / (dir * dir) + 
          freeCov(Acts::eFreeQOverP, Acts::eFreeQOverP) / (qop * qop) + 
          2 * freeCov(dir_idx, Acts::eFreeQOverP) / (dir * qop)));
  };

  m_err_px = calcMomentumErr(m_px, Acts::eFreeDir2);
  m_err_py = calcMomentumErr(m_py, Acts::eFreeDir1);
  m_err_pz = calcMomentumErr(m_pz, Acts::eFreeDir0);
}

RecoTrack::~RecoTrack() {}

} // namespace ActsExamples

