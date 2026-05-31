#include "ActsExamples/SHiP/SHiPMeasurementProvider.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include <map>
#include <cmath>
#include <algorithm>

namespace ActsExamples {

SHiPMeasurementProvider::SHiPMeasurementProvider(Config config)
    : m_cfg(std::move(config)) {}

Acts::GeometryIdentifier SHiPMeasurementProvider::getGeoId(const std::vector<float>& iHit) const {
    int detType = static_cast<int>(iHit[0]);
    int station = static_cast<int>(iHit[1]);
    int layer = static_cast<int>(iHit[2]);
    int view = static_cast<int>(iHit[3]);
    int straw = static_cast<int>(iHit[4]);

    if (detType == 0) { // STRAW
        int acts_layer = 16 * station + 4 * view + 2 * layer - 14;
        return Acts::GeometryIdentifier().withVolume(1).withLayer(acts_layer).withSensitive(straw);
    }

    int vol = (detType == 3) ? 2 : ((detType == 4) ? 3 : 1);
    int layer_id = 6 * static_cast<int>(iHit[0]) + 2 * static_cast<int>(iHit[2]) + 4;
    return Acts::GeometryIdentifier().withVolume(vol).withLayer(layer_id).withSensitive(0);
}

MeasurementContainer SHiPMeasurementProvider::process(
    const std::vector<std::vector<float>>& inputHits, 
    const Acts::GeometryContext& geoCtx) const {

    MeasurementContainer measurements;
    measurements.reserve(inputHits.size());

    std::map<Acts::GeometryIdentifier, std::vector<const std::vector<float>*>> surfaceGroups;
    for (const auto& iHit : inputHits) {
        surfaceGroups[getGeoId(iHit)].push_back(&iHit);
    }

    for (auto& [geoId, hits] : surfaceGroups) {
        const auto* surface = m_cfg.trackingGeometry->findSurface(geoId);
        if (!surface){ std::cout<<" not a surf"<<std::endl; continue; }

        auto surfaceId = surface->geometryId();
        if (geoId.volume() == 1) { // STRAW
            for (auto* h : hits) {
                DigitizedParameters d;
                d.indices = {Acts::eBoundLoc0};
                d.values = {0.0f};
                d.variances = {(m_cfg.strawRadius * m_cfg.strawRadius) / 3.0f};
                createMeasurement(measurements, surfaceId, d);
            }
        } else { // SEGMENTED (Silicon/SciFi): CoG Clustering
            std::sort(hits.begin(), hits.end(), [](const auto* a, const auto* b) { return (*a)[4] < (*b)[4]; });
            std::vector<std::vector<const std::vector<float>*>> clusters;
            for (auto* h : hits) {
                if (clusters.empty() || (*h)[4] > (*clusters.back().back())[4] + 1) clusters.push_back({h});
                else clusters.back().push_back(h);
            }

            bool isSil = (geoId.volume() == 2 || static_cast<int>((*hits[0])[0]) == 1);
            double pitch = isSil ? m_cfg.siliconPitch : m_cfg.scifiPitch;
            double res   = isSil ? m_cfg.siliconRes : m_cfg.scifiRes;

            for (const auto& cluster : clusters) {
                double sumLocX = 0, totalWeight = 0;
                for (auto* h : cluster) {
                    Acts::Vector3 gPos((*h)[8], (*h)[7], -1.0 * (*h)[6]);
                    auto lRes = surface->globalToLocal(geoCtx, gPos, Acts::Vector3(0,0,0));
                    if (!lRes.ok()) continue;

                    int binIdx = std::floor((lRes.value()[0] - m_cfg.minBound) / pitch);
                    double segX = m_cfg.minBound + (binIdx + 0.5) * pitch;
                    double weight = ((*h)[17] > 0) ? (*h)[17] : 1.0;
                    sumLocX += segX * weight;
                    totalWeight += weight;
                }
                if (totalWeight <= 0) continue;

                DigitizedParameters d;
                d.indices = {Acts::eBoundLoc0};
                d.values = {sumLocX / totalWeight};
                d.variances = {((res * res) / static_cast<double>(cluster.size())) + (m_cfg.noiseFloor * m_cfg.noiseFloor)};
                createMeasurement(measurements, surfaceId, d);
            }
        }
    }
    return measurements;
}

} // namespace ActsExamples
