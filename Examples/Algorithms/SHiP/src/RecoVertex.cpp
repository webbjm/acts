#include "ActsExamples/SHiP/RecoVertex.hpp"
#include <cmath>

namespace ActsExamples {

ClassImp(RecoVertex);

RecoVertex::RecoVertex() : TObject() {}

RecoVertex::RecoVertex(const Acts::Vertex& vtx) : TObject() {
    const auto& pos = vtx.fullPosition();
    m_x = pos[2] * 0.1; //convert to cm
    m_y = pos[1] * 0.1;
    m_z = -pos[0] * 0.1;
    m_t = pos[3];

    const auto& cov = vtx.fullCovariance();
    m_err_x = (cov(2, 2) > 0.0) ? std::sqrt(cov(2, 2)) * 0.1: 0.0;
    m_err_y = (cov(1, 1) > 0.0) ? std::sqrt(cov(1, 1)) * 0.1: 0.0;
    m_err_z = (cov(0, 0) > 0.0) ? std::sqrt(cov(0, 0)) * 0.1: 0.0;

    auto quality = vtx.fitQuality();
    m_chi2 = quality.first;
    m_nDoF = quality.second;

    m_trackIds.clear(); // Will be populated in the C++ Pybind layer
}

RecoVertex::~RecoVertex() {}

} // namespace ActsExamples

