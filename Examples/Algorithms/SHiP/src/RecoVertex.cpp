// -------------------------------------------------------------------------
// -----                   RecoTrack source file                   -----
#include "ActsExamples/SHiP/RecoVertex.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include <TBuffer.h>
#include <TClass.h>
#include <TROOT.h>
#include <cmath>


namespace ActsExamples {

ClassImp(RecoVertex);

RecoVertex::RecoVertex() : TObject() {}

RecoVertex::RecoVertex(const Acts::Vertex& vtx) : TObject() {
    // Position and Time from the 4-vector
    const auto& pos = vtx.fullPosition();
    m_x = pos[0];
    m_y = pos[1];
    m_z = pos[2];
    m_t = pos[3];

    // Covariance (Errors) - must use fullCovariance for 4D
    const auto& cov = vtx.fullCovariance();
    m_err_x = std::sqrt(cov(0, 0));
    m_err_y = std::sqrt(cov(1, 1));
    m_err_z = std::sqrt(cov(2, 2));
    //m_err_t = std::sqrt(cov(3, 3));

    // Extract Chi2 and nDoF from fitQuality()
    auto quality = vtx.fitQuality();
    m_chi2 = quality.first;
    m_nDoF = quality.second;

    // Tracks associated with the vertex (currently empty as indexing is external)
    m_trackIds.clear();
    const auto& vertexTracks = vtx.tracks();
    for (const auto& trkAtVtx : vertexTracks) {
        // Acts::TrackAtVertex usually holds a pointer or proxy to the track
        // Depending on your ACTS version, it may store the 'trackIndex'
        // or a pointer to the original track.
        // Assuming your 'trackIndex' is mapped or available:
      //  m_trackIds.push_back(static_cast<Int_t>(trkAtVtx.trackIndex));
    }

}

RecoVertex::~RecoVertex() {}

} // namespace ActsExamples                            
