#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "Rtypes.h"            // for Double_t, Int_t, Double32_t, etc
#include "TObject.h"           // for TObject

namespace ActsExamples {

class RecoTrack : public TObject {
 public:
  RecoTrack();
  // Constructor from the track object (proxy) to get chi2 and nDoF
  RecoTrack(const ConstTrackProxy& track, const Acts::GeometryContext& gctx = Acts::GeometryContext());

  RecoTrack(const Acts::BoundVector& parameters,
            const Acts::BoundSquareMatrix& covariance,
            const Acts::Surface& surface,
            float chi2,
            unsigned int nDoF,
            const Acts::GeometryContext& gctx);

  /**  Destructor  **/
  ~RecoTrack() override;

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

  ClassDef(RecoTrack, 1);
};
} //namespace ActsExamples                             
