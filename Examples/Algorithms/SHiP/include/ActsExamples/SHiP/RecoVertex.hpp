#pragma once

#include "Rtypes.h"            // for Double_t, Int_t, Double32_t, etc
#include "TObject.h"           // for TObject
#include "Acts/Vertexing/Vertex.hpp" // Required for Acts::Vertex
#include <vector>

namespace ActsExamples {

class RecoVertex : public TObject {
 public:
  RecoVertex();
  RecoVertex(const Acts::Vertex& vtx);
  ~RecoVertex() override;

  Double_t x() const { return m_x; }
  Double_t y() const { return m_y; }
  Double_t z() const { return m_z; }
  Double_t t() const { return m_t; }

  Double_t err_x() const { return m_err_x; }
  Double_t err_y() const { return m_err_y; }
  Double_t err_z() const { return m_err_z; }
  Double_t err_t() const { return m_err_t; }

  Double_t chi2() const { return m_chi2; }
  Double_t nDoF() const { return m_nDoF; }

  const std::vector<Int_t>& trackIds() const { return m_trackIds; }

 private:
  Double_t m_x{0.}, m_y{0.}, m_z{0.}, m_t{0.};
  Double_t m_err_x{0.}, m_err_y{0.}, m_err_z{0.}, m_err_t{0.};
  Double_t m_chi2{0.}, m_nDoF{0.};

  std::vector<Int_t> m_trackIds;

  ClassDef(RecoVertex, 1);
};

} // namespace ActsExamples
