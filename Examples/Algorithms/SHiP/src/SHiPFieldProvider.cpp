#include "ActsExamples/SHiP/SHiPFieldProvider.hpp"
#include "TFile.h"
#include "TTree.h"
#include <algorithm>
#include <system_error>

SHiPFieldProvider::SHiPFieldProvider(const std::string& filename, double scale) : m_scale(scale) {
    TFile* f = TFile::Open(filename.c_str());
    TTree* tR = static_cast<TTree*>(f->Get("Range"));
    tR->SetBranchAddress("xMin", &m_b.xMin); tR->SetBranchAddress("xMax", &m_b.xMax); tR->SetBranchAddress("dx", &m_b.dx);
    tR->SetBranchAddress("yMin", &m_b.yMin); tR->SetBranchAddress("yMax", &m_b.yMax); tR->SetBranchAddress("dy", &m_b.dy);
    tR->SetBranchAddress("zMin", &m_b.zMin); tR->SetBranchAddress("zMax", &m_b.zMax); tR->SetBranchAddress("dz", &m_b.dz);
    tR->GetEntry(0);

    m_nx = static_cast<int>((m_b.xMax - m_b.xMin) / m_b.dx) + 1;
    m_ny = static_cast<int>((m_b.yMax - m_b.yMin) / m_b.dy) + 1;
    m_nz = static_cast<int>((m_b.zMax - m_b.zMin) / m_b.dz) + 1;

    TTree* tD = static_cast<TTree*>(f->Get("Data"));
    float bx, by, bz;
    tD->SetBranchAddress("Bx", &bx); tD->SetBranchAddress("By", &by); tD->SetBranchAddress("Bz", &bz);
    for(int i=0; i<tD->GetEntries(); ++i) {
        tD->GetEntry(i);
        m_data.push_back(Acts::Vector3(bx, by, bz));
    }
    f->Close();
}

Acts::MagneticFieldProvider::Cache SHiPFieldProvider::makeCache(const Acts::MagneticFieldContext& /*mctx*/) const {
    return Acts::MagneticFieldProvider::Cache();
}

Acts::Result<Acts::Vector3> SHiPFieldProvider::getField(const Acts::Vector3& pos_acts, MagneticFieldProvider::Cache& /*cache*/) const {
    // Transform ACTS position (mm) to SHiP frame coordinates (cm)
    double pos_acts_cm_x = pos_acts.x() / 10.0;
    double pos_acts_cm_y = pos_acts.y() / 10.0;
    double pos_acts_cm_z = pos_acts.z() / 10.0;

    // Coordinate mapping: Rotate Z-Beam (SHiP) to X-Beam (ACTS) about Y-axis
    // X_Ship = -Z_Acts, Y_Ship = Y_Acts, Z_Ship = X_Acts
    double fileX = pos_acts_cm_z;
    double fileY = pos_acts_cm_y;
    double fileZ = pos_acts_cm_x;

    // Apply global offsets to get map-local coordinates
    double offset_z = 8957.0; // cm
    double mapX = fileX - 0.0;
    double mapY = fileY - 0.0;
    double mapZ = fileZ - offset_z;

    // Clamping for interpolation stability
    double clampedX = std::clamp(static_cast<float>(mapX), m_b.xMin, m_b.xMax);
    double clampedY = std::clamp(static_cast<float>(mapY), m_b.yMin, m_b.yMax);
    double clampedZ = std::clamp(static_cast<float>(mapZ), m_b.zMin, m_b.zMax);

    // Trilinear Interpolation
    auto getIdx = [&](int ix, int iy, int iz) { return (ix * m_ny + iy) * m_nz + iz; };
    int ix = std::clamp(static_cast<int>((clampedX - m_b.xMin)/m_b.dx), 0, m_nx-2);
    int iy = std::clamp(static_cast<int>((clampedY - m_b.yMin)/m_b.dy), 0, m_ny-2);
    int iz = std::clamp(static_cast<int>((clampedZ - m_b.zMin)/m_b.dz), 0, m_nz-2);

    double dx = (clampedX - (m_b.xMin + ix*m_b.dx))/m_b.dx;
    double dy = (clampedY - (m_b.yMin + iy*m_b.dy))/m_b.dy;
    double dz = (clampedZ - (m_b.zMin + iz*m_b.dz))/m_b.dz;

    auto i_val = [&](int i, int j, int k) { return m_data[getIdx(i, j, k)]; };

    Acts::Vector3 B_ship = (
        i_val(ix,iy,iz)*(1-dx)*(1-dy)*(1-dz) + i_val(ix+1,iy,iz)*dx*(1-dy)*(1-dz) +
        i_val(ix,iy+1,iz)*(1-dx)*dy*(1-dz) + i_val(ix,iy,iz+1)*(1-dx)*(1-dy)*dz +
        i_val(ix+1,iy+1,iz)*dx*dy*(1-dz) + i_val(ix+1,iy,iz+1)*dx*(1-dy)*dz +
        i_val(ix,iy+1,iz+1)*(1-dx)*dy*dz + i_val(ix+1,iy+1,iz+1)*dx*dy*dz
    ) * m_scale;

      Acts::Vector3 B_acts(B_ship.z(), B_ship.y(), -B_ship.x());
  
      return Acts::Result<Acts::Vector3>::success(B_acts);
}

Acts::Result<Acts::Vector3> SHiPFieldProvider::getFieldGradient(const Acts::Vector3& /*position*/,
                                                              Acts::ActsMatrix<3, 3>& /*derivative*/,
                                                              MagneticFieldProvider::Cache& /*cache*/) const {
    return Acts::Result<Acts::Vector3>::failure(std::make_error_code(std::errc::not_supported));
}
