#pragma once
#include "Acts/Definitions/Algebra.hpp"
#include <cmath>

namespace ActsExamples {

/// Standard rotation matrix for SHiP reconstruction frame
/// Transforms TGeo coordinates to the reconstruction frame.
inline Acts::RotationMatrix3 GetRecoRotation() {
    Acts::RotationMatrix3 rot;
    // Rotate by 90 degrees (pi/2) about Y-axis
    // The specific rotation frame requested for reconstruction
    rot << 0, 0, 1,
           0, 1, 0,
          -1, 0, 0;
    return rot;
}

/// Helper to apply this rotation and scaling to TGeo transforms
inline Acts::Transform3 ApplyRecoRotation(const Acts::Transform3& original) {
    Acts::RotationMatrix3 recoRot = GetRecoRotation();
    Acts::Transform3 recoTrafo = Acts::Transform3::Identity();
    recoTrafo.linear() = recoRot * original.linear();
    recoTrafo.translation() = recoRot * original.translation();
    return recoTrafo;
}

} // namespace ActsExamples
