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

    m_trackIds.clear(); 
    m_trackPx.clear(); m_trackPy.clear(); m_trackPz.clear();
    m_trackX.clear();  m_trackY.clear();  m_trackZ.clear();

    Acts::Vector3 globalVtxPos = vtx.position();
    for (size_t i = 0; i < vtx.tracks().size(); ++i) {
        const auto& trackAtVtx = vtx.tracks()[i];


        Acts::Vector3 mom = trackAtVtx.fittedParams.momentum();  

        m_trackPx.push_back(-mom.z()); // x_ship = -z_acts
        m_trackPy.push_back(mom.y());  // y_ship = y_acts
        m_trackPz.push_back(mom.x());  // z_ship = x_acts

        m_trackX.push_back(globalVtxPos.z() / 10.0);
        m_trackY.push_back(globalVtxPos.y() / 10.0);
        m_trackZ.push_back(-globalVtxPos.x() / 10.0);
    }
}

RecoVertex::~RecoVertex() {}

} // namespace ActsExamples

