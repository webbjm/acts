#include "ActsExamples/SHiP/SHiPMeasurementProvider.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include <algorithm>
#include <cmath>
#include <map>

namespace ActsExamples {

SHiPMeasurementProvider::SHiPMeasurementProvider(Config config, Acts::Logging::Level level)
    : IAlgorithm("SHiPMeasurementProvider", level), m_cfg(std::move(config)) {
    m_outputMeasurements.initialize(m_cfg.outputMeasurements);
}

Acts::GeometryIdentifier SHiPMeasurementProvider::getGeoId(const std::vector<float>& iHit) const {
    int detType = static_cast<int>(iHit[0]);
    int station = static_cast<int>(iHit[1]);
    
    // detType: 0=Straw, 1=Silicon, 2=SciFi, 3=SND_Silicon, 4=SND_SciFi
    if (detType == 0) { // STRAW
        int sens = (static_cast<int>(iHit[3]) == 0 || static_cast<int>(iHit[3]) == 3) ? 
                   (static_cast<int>(iHit[4]) + station * 300) : 
                   (static_cast<int>(iHit[4]) + station * 315);
        return Acts::GeometryIdentifier().withVolume(1).withLayer(8 * static_cast<int>(iHit[1]) + 2 * static_cast<int>(iHit[3]) - 6).withSensitive(sens);
    } 
    
    // Assign Volume IDs: 1 (Standalone), 2 (SND Silicon), 3 (SND SciFi)
    int vol = 1;
    if (detType == 3) vol = 2;
    else if (detType == 4) vol = 3;
    
    return Acts::GeometryIdentifier().withVolume(vol).withLayer(6 * static_cast<int>(iHit[0]) + 2 * static_cast<int>(iHit[2]) + 4).withSensitive(0);
}

ProcessCode SHiPMeasurementProvider::execute(const AlgorithmContext& ctx) const {
    if (!m_cfg.inputHitArray) return ProcessCode::ABORT;

    MeasurementContainer measurements;
    measurements.reserve(m_cfg.inputHitArray->size());

    std::map<Acts::GeometryIdentifier, std::vector<const std::vector<float>*>> surfaceGroups;
    for (const auto& iHit : *m_cfg.inputHitArray) {
        surfaceGroups[getGeoId(iHit)].push_back(&iHit);
    }

    for (auto& [geoId, hits] : surfaceGroups) {
        const auto* surface = m_cfg.trackingGeometry->findSurface(geoId);
        if (!surface) continue;

        if (geoId.volume() == 1 && static_cast<int>((*hits[0])[0]) == 0) { // STRAW
            for (auto* h : hits) {
                DigitizedParameters d;
                d.indices = {Acts::eBoundLoc0};
                d.values = {(*h)[6]}; //5-?6 
                d.variances = {m_cfg.strawRes * m_cfg.strawRes};
                createMeasurement(measurements, geoId, d);
            }
        } else { // SEGMENTED (Silicon/SciFi): CoG Clustering
            std::sort(hits.begin(), hits.end(), [](auto* a, auto* b) { return (*a)[4] < (*b)[4]; }); //3->4
            std::vector<std::vector<const std::vector<float>*>> clusters;
            for (auto* h : hits) {
                if (clusters.empty() || (*h)[4] > (*clusters.back().back())[4] + 1) clusters.push_back({h}); //3->4
                else clusters.back().push_back(h);
            }

            // Volume 2 is Silicon, Volume 3 is SciFi
            bool isSil = (geoId.volume() == 2 || static_cast<int>((*hits[0])[0]) == 1);
            double pitch = isSil ? m_cfg.siliconPitch : m_cfg.scifiPitch;
            double res   = isSil ? m_cfg.siliconRes : m_cfg.scifiRes;

            for (const auto& cluster : clusters) {
                double sumLocX = 0, totalWeight = 0;
                for (auto* h : cluster) {
                    Acts::Vector3 gPos((*h)[8], (*h)[7], -1 * (*h)[6]); //5,6,7 -> 6,7,8 and rotate
                    auto lRes = surface->globalToLocal(ctx.geoContext, gPos, Acts::Vector3(0,0,0));
                    if (!lRes.ok()) continue;
                    
                    int binIdx = std::floor((lRes.value()[0] - m_cfg.minBound) / pitch);
                    double segX = m_cfg.minBound + (binIdx + 0.5) * pitch;
                    double weight = (*h)[17] > 0 ? (*h)[17] : 1.0; //16-?17
                    sumLocX += segX * weight;
                    totalWeight += weight;
                }
                if (totalWeight <= 0) continue;
                DigitizedParameters d;
                d.indices = {Acts::eBoundLoc0};
                d.values = {sumLocX / totalWeight};
                d.variances = {((res * res) / static_cast<double>(cluster.size())) + (m_cfg.noiseFloor * m_cfg.noiseFloor)};
                createMeasurement(measurements, geoId, d);
            }
        }
    }
    m_outputMeasurements(ctx, std::move(measurements));
    return ProcessCode::SUCCESS;
}
} // namespace ActsExamples
