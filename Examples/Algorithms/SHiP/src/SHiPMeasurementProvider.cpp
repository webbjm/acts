#include "ActsExamples/SHiP/SHiPMeasurementProvider.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include <map>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace ActsExamples {

SHiPMeasurementProvider::SHiPMeasurementProvider(Config config)
    : m_cfg(std::move(config)) {}

Acts::GeometryIdentifier SHiPMeasurementProvider::getGeoId(const std::vector<float>& iHit) const {
    int detType = static_cast<int>(iHit[0]);
    int station = static_cast<int>(iHit[1]);
    int layer = static_cast<int>(iHit[2]);
    int view = static_cast<int>(iHit[3]);
    int channel = static_cast<int>(iHit[4]);

    if (detType == 0) { // STRAW
        int acts_layer = 16 * station + 4 * view + 2 * layer - 14;

        return Acts::GeometryIdentifier().withVolume(1).withLayer(acts_layer).withSensitive(channel);
    }

    int vol = (detType == 3) ? 2 : ((detType == 4) ? 3 : 1);
    int layer_id = 6 * station + 2 * layer + 4;
    return Acts::GeometryIdentifier().withVolume(vol).withLayer(layer_id).withSensitive(view);
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
                std::cout<<"Debug: drift distance "<<(*h)[11]*10<<std::endl;
                double sign = 1.0;
                if (h->size() > 16) {
                    sign = (*h)[16];
                }
                d.values = {(*h)[11] * 10 * sign}; //Drift radius in mm, h[16] is the sign for l/r ambiguities
                d.variances = {m_cfg.strawRes * m_cfg.strawRes};
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
                    double sumStrip = 0, totalWeight = 0;
                    for (auto* h : cluster) {
                        double channel = (*h)[4]; // Channel number
                        double weight = ((*h)[10] > 0) ? (*h)[10] : 1.0;
               
                        sumStrip += channel * weight;
                        totalWeight += weight;
                    }
               
                    double avgStrip = sumStrip / totalWeight;
                    //double localX = (avgStrip - center_strip) * pitch;
                    double center_strip = 649.5; 
                    double localX = (avgStrip - center_strip) * pitch;
//                    std::cout<<avgStrip<<std::endl;
//                    std::cout<<avgStrip-center_strip<<std::endl;
                    std::cout<<"Geometry in ACTS: "<<surfaceId.layer()<<" "<<surfaceId.sensitive()<<std::endl;
                    std::cout<<"Cluster: x-local, size, strip num "<<localX<<" "<<cluster.size()<<" "<<avgStrip<<std::endl;

                    DigitizedParameters d;
                    d.indices = {Acts::eBoundLoc0};
                    d.values = {localX};
                    d.variances = {((res * res) / static_cast<double>(cluster.size())) + (m_cfg.noiseFloor * m_cfg.noiseFloor)};
                    // Fill the Cluster object inside DigitizedParameters
                    d.cluster.sizeLoc0 = cluster.size(); 
                    d.cluster.sizeLoc1 = 1;
                    createMeasurement(measurements, surfaceId, d);

                }
            }
    }
    return measurements;
}

} // namespace ActsExamples
