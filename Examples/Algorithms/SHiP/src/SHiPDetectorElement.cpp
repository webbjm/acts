#include "ActsExamples/SHiP/SHiPDetectorElement.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"

namespace ActsExamples {

SHiPDetectorElement::SHiPDetectorElement(
    std::shared_ptr<const Acts::Transform3> transform,
    std::shared_ptr<const Acts::PlanarBounds> pBounds, double thickness,
    std::shared_ptr<const Acts::ISurfaceMaterial> material)
    : Acts::SurfacePlacementBase(),
      m_elementTransform(std::move(transform)),
      m_elementSurface(
          Acts::Surface::makeShared<Acts::PlaneSurface>(pBounds, *this)),
      m_elementThickness(thickness) {
  auto mutableSurface = std::const_pointer_cast<Acts::Surface>(m_elementSurface);
  mutableSurface->assignSurfaceMaterial(std::move(material));
}

SHiPDetectorElement::SHiPDetectorElement(
    std::shared_ptr<const Acts::Transform3> transform,
    std::shared_ptr<const Acts::DiscBounds> dBounds, double thickness,
    std::shared_ptr<const Acts::ISurfaceMaterial> material)
    : Acts::SurfacePlacementBase(),
      m_elementTransform(std::move(transform)),
      m_elementSurface(
          Acts::Surface::makeShared<Acts::DiscSurface>(dBounds, *this)),
      m_elementThickness(thickness) {
  auto mutableSurface = std::const_pointer_cast<Acts::Surface>(m_elementSurface);
  mutableSurface->assignSurfaceMaterial(std::move(material));
}

SHiPDetectorElement::SHiPDetectorElement(
    std::shared_ptr<const Acts::Transform3> transform,
    std::shared_ptr<const Acts::LineBounds> lBounds, double thickness,
    std::shared_ptr<const Acts::ISurfaceMaterial> material)
    : Acts::SurfacePlacementBase(),
      m_elementTransform(std::move(transform)),
      m_elementSurface(
          Acts::Surface::makeShared<Acts::StrawSurface>(lBounds, *this)),
      m_elementThickness(thickness) {
  auto mutableSurface = std::const_pointer_cast<Acts::Surface>(m_elementSurface);
  mutableSurface->assignSurfaceMaterial(std::move(material));
}

// Fixed method signature from 'transform' to 'localToGlobalTransform'
const Acts::Transform3& SHiPDetectorElement::localToGlobalTransform(
    const Acts::GeometryContext& gctx) const {

  // Check if the type-erased container actually holds data
  if (gctx.hasValue() && !m_alignedTransforms.empty()) {
    try {
      // Cast the std::any container back to your custom alignment struct
      auto alignContext = std::any_cast<ContextType>(gctx);

      if (alignContext.iov < m_alignedTransforms.size() && m_alignedTransforms[alignContext.iov]) {
        return *m_alignedTransforms[alignContext.iov];
      }
    } catch (const std::bad_any_cast&) {
      // Fallback or log if a different context type was passed accidentally
    }
  }
  return nominalTransform(gctx);
}


const Acts::Transform3& SHiPDetectorElement::nominalTransform(
    const Acts::GeometryContext& /*gctx*/) const {
  return *m_elementTransform;
}



void SHiPDetectorElement::addAlignedTransform(
    std::unique_ptr<Acts::Transform3> alignedTransform, unsigned int iov) {
  if (iov >= m_alignedTransforms.size()) {
    m_alignedTransforms.resize(iov + 1);
  }
  m_alignedTransforms[iov] = std::move(alignedTransform);
}

}  // namespace ActsExamples

