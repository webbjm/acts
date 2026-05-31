#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <memory>
#include <utility>
#include <vector>

namespace Acts {
class Surface;
class PlanarBounds;
class DiscBounds;
class LineBounds;
class ISurfaceMaterial;
}  // namespace Acts

namespace ActsExamples {

/// @class ShipDetectorElement
///
/// This is a lightweight type of detector element.
/// It implements the base class to provide surfaces to the tracking geometry.
class SHiPDetectorElement : public Acts::DetectorElementBase {
 public:
  /// Context type for potential alignment
  struct ContextType {
    unsigned int iov = 0;
  };

  /// Constructor for planar detector element
  SHiPDetectorElement(
      std::shared_ptr<const Acts::Transform3> transform,
      std::shared_ptr<const Acts::PlanarBounds> pBounds, double thickness,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr);

  /// Constructor for disc like detector element
  SHiPDetectorElement(
      std::shared_ptr<const Acts::Transform3> transform,
      std::shared_ptr<const Acts::DiscBounds> dBounds, double thickness,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr);

  /// Constructor for straw detector element
  SHiPDetectorElement(
      std::shared_ptr<const Acts::Transform3> transform,
      std::shared_ptr<const Acts::LineBounds> lBounds, double thickness,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr);

  ~SHiPDetectorElement() override = default;

  /// Return surface associated with this detector element
  const Acts::Surface& surface() const final { return *m_elementSurface; }

  /// Non-const access to the surface associated with this detector element
  Acts::Surface& surface() final { return *m_elementSurface; }

  /// Return the element thickness
  double thickness() const final { return m_elementThickness; }

  /// Return local to global transform associated with this context
  const Acts::Transform3& transform(
      const Acts::GeometryContext& gctx) const final;

  /// Return the nominal local to global transform
  const Acts::Transform3& nominalTransform(
      const Acts::GeometryContext& gctx) const;

  /// Add an aligned transform
  void addAlignedTransform(std::unique_ptr<Acts::Transform3> alignedTransform,
                           unsigned int iov);

 private:
  /// The transform for positioning in 3D space
  std::shared_ptr<const Acts::Transform3> m_elementTransform = nullptr;
  /// The aligned transforms
  std::vector<std::unique_ptr<Acts::Transform3>> m_alignedTransforms = {};
  /// The surface represented by it
  std::shared_ptr<Acts::Surface> m_elementSurface = nullptr;
  /// The element thickness
  double m_elementThickness = 0.;
};

}  // namespace ActsExamples

