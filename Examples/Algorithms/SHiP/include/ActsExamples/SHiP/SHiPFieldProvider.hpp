#pragma once
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include <vector>
#include <string>
#include <memory>
#include <system_error>

class SHiPFieldProvider : public Acts::MagneticFieldProvider {
public:
    struct GridBounds { float xMin, xMax, dx, yMin, yMax, dy, zMin, zMax, dz; };

    SHiPFieldProvider(const std::string& filename, double scale);

    // ACTS required methods
    Acts::MagneticFieldProvider::Cache makeCache(const Acts::MagneticFieldContext& mctx) const override;

    Acts::Result<Acts::Vector3> getField(const Acts::Vector3& position,
                                         MagneticFieldProvider::Cache& cache) const override;

    Acts::Result<Acts::Vector3> getFieldGradient(const Acts::Vector3& position,
                                                 Acts::ActsMatrix<3, 3>& derivative,
                                                 MagneticFieldProvider::Cache& cache) const;
private:
    GridBounds m_b;
    int m_nx, m_ny, m_nz;
    std::vector<Acts::Vector3> m_data;
    double m_scale;
};
