#ifndef MCPAnalyzer_Analyzer_MCPProbQ_h
#define MCPAnalyzer_Analyzer_MCPProbQ_h

// -*- C++ -*-
//
// MCPProbQ - self-contained dE/dx helper for the MCPAnalyzer.
//
// Our own version of the pixel "probQ" discriminant and a minimal harmonic-mean
// dE/dx, ported from HSCPAnalysis/Ntuplizer/src/HSCPDeDxTool.cc but with no
// dependency on that package (no templates / scale-factor / cluster-cleaning
// machinery). Everything needed for studies 1 & 2 lives here:
//   - probQonTrack / probQonTrackNoL1  (pixel cluster shape probability)
//   - Ih harmonic-mean dE/dx, full / pixel-only / strip-only (MeV/cm)
//   - per-subdetector hit counts (does the track reach the strips?)
//
// Header-only so the analyzer plugin can include it directly.

#include <cmath>
#include <vector>

#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
// Pulls in LocalTrajectoryParameters, SiPixelCluster and SiPixelRecHitQuality:
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"

struct MCPDeDxResult {
  // Pixel probQ (track level), with the gated probXY (probQ<0.8 && probXY>0) as in the standard tool
  float probQonTrack = -1.f;     // == probQ_pixel
  float probXYonTrack = -1.f;    // == probXY_pixel (gated)
  float probQonTrackNoL1 = -1.f;
  float probXYonTrackNoL1 = -1.f;

  int nPixHitsUsed = 0;   // pixel hits entering probQ
  int nonL1PixHits = 0;   // pixel hits excluding BPIX layer 1

  // Per-hit probXY diagnostics (to tell a genuine probXY~0 from a degenerate one).
  int nPixClusters = 0;      // pixel hits with a SiPixelCluster in the DeDxHitInfo
  int nPixNoFillProb = 0;    // CPE did not fill the probability word (hasFilledProb == false)
  int nPixSpecInCPE = 0;     // excluded: on edge / bad pixels / spans two ROCs
  int nPixXYpinnedLo = 0;    // raw CPE probXY <= 0  (shape totally inconsistent)
  int nPixXYpinnedHi = 0;    // raw CPE probXY >= 1
  int nPixXYvalid = 0;       // raw CPE probXY in (0,1) -> a real measurement
  float pixXYrawMin = -1.f;  // smallest raw probXY seen among filled hits

  // Harmonic-mean dE/dx (MeV/cm) and number of measurements
  float ih_full = -1.f;  int nom_full = 0;
  float ih_pixel = -1.f; int nom_pixel = 0;
  float ih_strip = -1.f; int nom_strip = 0;

  // Per-hit calibrated dE/dx (MeV/cm), for building a probQ-style strip discriminant downstream.
  std::vector<float> pixelDedxHits;
  std::vector<float> stripDedxHits;

  // Per-subdetector hit counts (from the DeDxHitInfo cluster list)
  int nBPIX = 0, nFPIX = 0, nTIB = 0, nTID = 0, nTOB = 0, nTEC = 0;
  int nPix = 0, nStrip = 0;

  // Cluster-shape observables (independent of charge-based probQ; survive ADC saturation).
  // Pixel: cluster size = number of pixels above threshold; charge in electrons.
  // Strip: cluster width = number of strips; charge = sum of ADC.
  float pixClustSizeMean = -1.f, pixClustSizeXMean = -1.f, pixClustSizeYMean = -1.f;
  int   pixClustSizeMax = 0;
  float pixClustChargeMean = -1.f;
  // Geometric prediction of the r-phi cluster size from the local track angle (thickness*|tan a_x|/pitchX);
  // the angle-corrected size = measured sizeX - predicted removes the trivial angle dependence.
  float pixSizeXpredMean = -1.f;     // mean predicted r-phi size [pixels]
  float pixSizeXresidualMean = -1.f; // mean (measured - predicted) r-phi size [pixels]
  float stripClustWidthMean = -1.f; int stripClustWidthMax = 0;
  float stripClustChargeMean = -1.f;
};

class MCPProbQ {
public:
  // A per-hit charge probability at or below this floor is treated as a saturated /
  // off-template cluster, i.e. NOT a measurement of "probability ~ 0" but a missing
  // measurement. Such hits are skipped from the combination so a single degenerate hit
  // cannot drive the float product to underflow (track-level probQ == 0). Without this,
  // tiny-positive per-hit values multiply down to exactly 0 even for MIP background.
  static constexpr float kProbHitFloor = 1e-6f;

  // Poisson-series combination of per-hit probabilities into a track-level value
  // (identical to HSCPDeDxTool::combineProbs).
  static float combineProbs(float probOnTrackWMulti, int numRecHits) {
    float logprob = (probOnTrackWMulti > 0) ? std::log(probOnTrackWMulti) : 0;
    float factQ = -logprob;
    float term = 0.f;
    if (numRecHits == 1) {
      term = 1.f;
    } else if (numRecHits > 1) {
      term = 1.f + factQ;
      for (int i = 2; i < numRecHits; ++i) {
        factQ *= -logprob / float(i);
        term += factQ;
      }
    }
    return probOnTrackWMulti * term;
  }

  // Compute everything for one track from its DeDxHitInfo + reco direction/charge.
  static MCPDeDxResult compute(const reco::DeDxHitInfo* dedxHits,
                               float track_px, float track_py, float track_pz, int track_charge,
                               const TrackerGeometry* tkGeometry,
                               const TrackerTopology* tkTopo,
                               const PixelClusterParameterEstimator* pixelCPE) {
    MCPDeDxResult r;
    if (!dedxHits || !tkGeometry || !tkTopo) return r;

    // Charge-to-energy conversion (MeV per unit cluster charge); pixel charge is in
    // electrons, strip charge in ADC (~265 e/ADC). Scale factors are the analysis
    // defaults (strip=1.0, pixel-to-strip=1.035).
    const float factorChargeToE[2] = {3.61e-06f, 3.61e-06f * 265.f};
    const float sfStrip = 1.0f;
    const float sfPixel = 1.035f;

    std::vector<float> dedxFull, dedxPixel, dedxStrip;
    std::vector<float> pixSize, pixSizeX, pixSizeY, pixCharge;
    std::vector<float> stripWidth, stripCharge;

    for (unsigned int h = 0; h < dedxHits->size(); ++h) {
      DetId detid(dedxHits->detId(h));
      const uint32_t subdet = detid.subdetId();
      const bool isPixel = (subdet < 3);

      // ---- per-subdetector hit bookkeeping ----
      if (isPixel) {
        r.nPix++;
        if (subdet == PixelSubdetector::PixelBarrel) r.nBPIX++;
        else if (subdet == PixelSubdetector::PixelEndcap) r.nFPIX++;
      } else {
        r.nStrip++;
        if (subdet == SiStripDetId::TIB) r.nTIB++;
        else if (subdet == SiStripDetId::TID) r.nTID++;
        else if (subdet == SiStripDetId::TOB) r.nTOB++;
        else if (subdet == SiStripDetId::TEC) r.nTEC++;
      }

      // ---- cluster-size / cluster-charge observables ----
      if (isPixel) {
        if (auto const* pc = dedxHits->pixelCluster(h)) {
          pixSize.push_back(float(pc->size()));
          pixSizeX.push_back(float(pc->sizeX()));
          pixSizeY.push_back(float(pc->sizeY()));
          pixCharge.push_back(float(pc->charge()));
          if (pc->size() > r.pixClustSizeMax) r.pixClustSizeMax = pc->size();
        }
      } else {
        if (auto const* sc = dedxHits->stripCluster(h)) {
          const auto& amps = sc->amplitudes();
          int w = amps.size();
          float q = 0.f; for (auto a : amps) q += float(a);
          stripWidth.push_back(float(w));
          stripCharge.push_back(q);
          if (w > r.stripClustWidthMax) r.stripClustWidthMax = w;
        }
      }

      // ---- harmonic-mean dE/dx contribution ----
      const float pathlength = dedxHits->pathlength(h);
      if (pathlength > 0.f) {
        float sf = (isPixel ? factorChargeToE[0] : factorChargeToE[1]) * (isPixel ? sfPixel : sfStrip);
        float dedx = sf * dedxHits->charge(h) / pathlength;
        dedxFull.push_back(dedx);
        if (isPixel) dedxPixel.push_back(dedx);
        else         dedxStrip.push_back(dedx);
      }
    }

    r.ih_full  = harmonicMean(dedxFull);   r.nom_full  = dedxFull.size();
    r.ih_pixel = harmonicMean(dedxPixel);  r.nom_pixel = dedxPixel.size();
    r.ih_strip = harmonicMean(dedxStrip);  r.nom_strip = dedxStrip.size();
    r.pixelDedxHits = dedxPixel;
    r.stripDedxHits = dedxStrip;

    r.pixClustSizeMean   = arithMean(pixSize);
    r.pixClustSizeXMean  = arithMean(pixSizeX);
    r.pixClustSizeYMean  = arithMean(pixSizeY);
    r.pixClustChargeMean = arithMean(pixCharge);
    r.stripClustWidthMean  = arithMean(stripWidth);
    r.stripClustChargeMean = arithMean(stripCharge);

    // ---- pixel probQ ----
    if (pixelCPE) computeProbQ(r, dedxHits, track_px, track_py, track_pz, track_charge,
                               tkGeometry, tkTopo, pixelCPE);
    return r;
  }

private:
  static float arithMean(const std::vector<float>& v) {
    if (v.empty()) return -1.f;
    double s = 0.; for (float x : v) s += x;
    return float(s / double(v.size()));
  }

  static float harmonicMean(const std::vector<float>& v) {
    if (v.empty()) return -1.f;
    double s = 0.;
    for (float x : v) {
      if (x > 0.f) s += 1.0 / (double(x) * double(x));
    }
    if (s <= 0.) return -1.f;
    return float(std::pow(s / double(v.size()), -0.5));
  }

  // Port of HSCPDeDxTool::computeProbQ (pixel hits only).
  static void computeProbQ(MCPDeDxResult& r,
                           const reco::DeDxHitInfo* dedxHits,
                           float track_px, float track_py, float track_pz, int track_charge,
                           const TrackerGeometry* tkGeometry,
                           const TrackerTopology* tkTopo,
                           const PixelClusterParameterEstimator* pixelCPE) {
    const int numLayers = tkGeometry->numberOfLayers(PixelSubdetector::PixelBarrel);

    int numRecHitsQ = 0, numRecHitsXY = 0;
    int numRecHitsQNoL1 = 0, numRecHitsXYNoL1 = 0;
    float probQonTrackWMulti = 1.f, probXYonTrackWMulti = 1.f;
    float probQonTrackWMultiNoL1 = 1.f, probXYonTrackWMultiNoL1 = 1.f;
    unsigned int nonL1PixHits = 0;
    double sumPred = 0., sumResidual = 0.; int nPred = 0;  // angle-corrected sizeX accumulators

    for (unsigned int i = 0; i < dedxHits->size(); ++i) {
      DetId detid(dedxHits->detId(i));
      if (detid.subdetId() >= 3) continue;  // pixels only

      auto const* pixelCluster = dedxHits->pixelCluster(i);
      if (pixelCluster == nullptr) continue;
      r.nPixClusters++;

      const GeomDetUnit* geomDet = dynamic_cast<const GeomDetUnit*>(tkGeometry->idToDetUnit(detid));
      if (geomDet == nullptr) continue;

      LocalVector lv = geomDet->toLocal(GlobalVector(track_px, track_py, track_pz));

      // geometric predicted r-phi cluster size from the local track angle: thickness*|tan a_x|/pitchX
      const auto* pgdu = dynamic_cast<const PixelGeomDetUnit*>(geomDet);
      if (pgdu && std::abs(lv.z()) > 1e-6f) {
        float pitchX = pgdu->specificTopology().pitch().first;          // cm
        float thickness = geomDet->surface().bounds().thickness();      // cm
        float pred = std::abs(thickness * lv.x() / lv.z()) / pitchX;     // pixels
        sumPred += pred; sumResidual += (float(pixelCluster->sizeX()) - pred); nPred++;
      }
      auto qword = std::get<2>(pixelCPE->getParameters(
          *pixelCluster, *geomDet, LocalTrajectoryParameters(dedxHits->pos(i), lv, track_charge)));

      if (!SiPixelRecHitQuality::thePacking.hasFilledProb(qword)) { r.nPixNoFillProb++; continue; }  // CPE failed

      float rawProbQ = SiPixelRecHitQuality::thePacking.probabilityQ(qword);
      float rawProbXY = SiPixelRecHitQuality::thePacking.probabilityXY(qword);
      // per-hit probXY diagnostics on the RAW CPE value (before pinning)
      if (rawProbXY <= 0.f)      r.nPixXYpinnedLo++;
      else if (rawProbXY >= 1.f) r.nPixXYpinnedHi++;
      else {
        r.nPixXYvalid++;
        if (r.pixXYrawMin < 0.f || rawProbXY < r.pixXYrawMin) r.pixXYrawMin = rawProbXY;
      }

      float probQ = rawProbQ, probXY = rawProbXY;
      if (probXY <= 0.f || probXY >= 1.f) probXY = 0.f;

      const bool specInCPE = SiPixelRecHitQuality::thePacking.isOnEdge(qword) ||
                             SiPixelRecHitQuality::thePacking.hasBadPixels(qword) ||
                             SiPixelRecHitQuality::thePacking.spansTwoROCs(qword);
      if (specInCPE) r.nPixSpecInCPE++;

      r.nPixHitsUsed++;
      // A per-hit probQ at/below the floor is a saturated / off-template cluster, not a
      // genuine "probability 0" measurement -> skip it from the combination (never combine
      // a 0). The <0.8 upper gate is the standard MIP-rejection requirement.
      const bool qUsable = (probQ > kProbHitFloor && probQ < 0.8f);
      if (!specInCPE && qUsable) { numRecHitsQ++; probQonTrackWMulti *= probQ; }
      if (!specInCPE && qUsable && probXY > 0.f) { numRecHitsXY++; probXYonTrackWMulti *= probXY; }

      // NoL1 variant: exclude BPIX layer 1 (noisy in 2017/2018)
      const bool isPhase0 = (numLayers == 3);
      const bool isPhase1NoL1 = (numLayers == 4) &&
          ((detid.subdetId() == PixelSubdetector::PixelEndcap) ||
           ((detid.subdetId() == PixelSubdetector::PixelBarrel) && (tkTopo->pxbLayer(detid) != 1)));
      if (isPhase0 || isPhase1NoL1) {
        nonL1PixHits++;
        if (!specInCPE && qUsable) { numRecHitsQNoL1++; probQonTrackWMultiNoL1 *= probQ; }
        if (!specInCPE && qUsable && probXY > 0.f) { numRecHitsXYNoL1++; probXYonTrackWMultiNoL1 *= probXY; }
      }
    }

    r.nonL1PixHits = nonL1PixHits;
    // If no usable charge hit survived, probQ is undefined (-1) rather than 0 -- a fully
    // saturated track is "no measurement", not "definitely signal at probQ==0".
    r.probQonTrack     = (numRecHitsQ     > 0) ? combineProbs(probQonTrackWMulti,     numRecHitsQ)     : -1.f;
    r.probXYonTrack    = (numRecHitsXY    > 0) ? combineProbs(probXYonTrackWMulti,    numRecHitsXY)    : -1.f;
    r.probQonTrackNoL1 = (numRecHitsQNoL1 > 0) ? combineProbs(probQonTrackWMultiNoL1, numRecHitsQNoL1) : -1.f;
    r.probXYonTrackNoL1= (numRecHitsXYNoL1> 0) ? combineProbs(probXYonTrackWMultiNoL1,numRecHitsXYNoL1): -1.f;
    if (nPred > 0) {
      r.pixSizeXpredMean = float(sumPred / nPred);
      r.pixSizeXresidualMean = float(sumResidual / nPred);
    }
  }
};

#endif
