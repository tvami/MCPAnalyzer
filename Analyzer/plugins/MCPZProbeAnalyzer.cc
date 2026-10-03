// -*- C++ -*-
//
// Package:    MCPAnalyzer/Analyzer
// Class:      MCPZProbeAnalyzer
//
// Z->mumu tag-and-probe MIP reference on RECO / RAW-RECO (e.g. the ZMu skim):
// probes are generalTracks matched to a reco::Muon, with dE/dx hits from the
// RECO dedxHitInfo association (all tracks with pT>10, no MiniAOD pT>50 cut).
// The track-tree branch names follow MCPAnalyzer, so the same macros run on both.

#include <memory>
#include <vector>
#include <cmath>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "RecoLocalTracker/Records/interface/TkPixelCPERecord.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"

#include "TTree.h"

#include "MCPProbQ.h"

class MCPZProbeAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MCPZProbeAnalyzer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  static bool passAny(const edm::TriggerResults&, const edm::TriggerNames&, const std::vector<std::string>&);

  const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  const edm::EDGetTokenT<reco::DeDxHitInfoAss> dedxToken_;
  const edm::EDGetTokenT<edm::ValueMap<int>> dedxPrescaleToken_;
  const edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  const edm::EDGetTokenT<trigger::TriggerEvent> trigEventToken_;
  const std::string tagFilter_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const edm::ESGetToken<PixelClusterParameterEstimator, TkPixelCPERecord> cpeToken_;
  const double probePtMin_, massMin_, massMax_;

  // same groups as MCPAnalyzer
  const std::vector<std::string> trigMET_ = {
      "HLT_MET105_IsoTrk50", "HLT_MET120_IsoTrk50", "HLT_PFMET105_IsoTrk50",
      "HLT_PFMET120_PFMHT120_IDTight", "HLT_PFMET120_PFMHT120_IDTight_PFHT60", "HLT_PFMET130_PFMHT130_IDTight",
      "HLT_PFMET140_PFMHT140_IDTight", "HLT_PFMET200_BeamHaloCleaned", "HLT_PFMET250_NotCleaned",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight",
      "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF", "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight",
      "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF", "HLT_PFMETTypeOne140_PFMHT140_IDTight",
      "HLT_PFMETTypeOne200_BeamHaloCleaned", "HLT_CaloMET350_NotCleaned"};
  const std::vector<std::string> trigJet_ = {"HLT_PFJet500", "HLT_AK8PFJet500", "HLT_PFHT1050", "HLT_CaloJet500_NoJetID"};
  const std::vector<std::string> trigMuon_ = {"HLT_IsoMu24", "HLT_Mu50"};

  TTree* tT_ = nullptr;
  TTree* tE_ = nullptr;

  unsigned int run_, lumi_; unsigned long long event_;
  int nPV_, nTag_, passMET_, passJet_, passMuon_, tagTrigObj_;
  double pt_, eta_, phi_; int charge_, highPurity_, nPixHit_, nTkLayers_;
  int hasDeDx_, dedxPrescale_;
  float probQ_, probQNoL1_, sizeXres_, sizeXpred_, clSizeX_, ihPixel_, ihStrip_;
  int nPixUsed_, nPixQFloor_, nBPIX_, nFPIX_, nPixClusters_;
  int muMatched_, muTight_; float muRelIso_;
  float tpMass_, tpTagPt_; int tpOS_;
};

MCPZProbeAnalyzer::MCPZProbeAnalyzer(const edm::ParameterSet& iC)
    : trackToken_(consumes<reco::TrackCollection>(iC.getParameter<edm::InputTag>("tracks"))),
      dedxToken_(consumes<reco::DeDxHitInfoAss>(iC.getParameter<edm::InputTag>("dedxHitInfo"))),
      dedxPrescaleToken_(consumes<edm::ValueMap<int>>(iC.getParameter<edm::InputTag>("dedxHitInfoPrescale"))),
      muonToken_(consumes<reco::MuonCollection>(iC.getParameter<edm::InputTag>("muons"))),
      vertexToken_(consumes<reco::VertexCollection>(iC.getParameter<edm::InputTag>("primaryVertices"))),
      triggerResultsToken_(consumes<edm::TriggerResults>(iC.getParameter<edm::InputTag>("triggerResults"))),
      trigEventToken_(consumes<trigger::TriggerEvent>(iC.getParameter<edm::InputTag>("triggerEvent"))),
      tagFilter_(iC.getParameter<std::string>("tagFilter")),
      topoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>()),
      geomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()),
      cpeToken_(esConsumes<PixelClusterParameterEstimator, TkPixelCPERecord>(
          edm::ESInputTag("", iC.getParameter<std::string>("pixelCPE")))),
      probePtMin_(iC.getParameter<double>("probePtMin")),
      massMin_(iC.getParameter<double>("massMin")),
      massMax_(iC.getParameter<double>("massMax")) {
  usesResource("TFileService");
}

bool MCPZProbeAnalyzer::passAny(const edm::TriggerResults& tr, const edm::TriggerNames& names,
                                const std::vector<std::string>& list) {
  for (unsigned int i = 0; i < tr.size(); ++i) {
    if (!tr.accept(i)) continue;
    const std::string& n = names.triggerName(i);
    for (const auto& t : list)
      if (n.compare(0, t.size() + 2, t + "_v") == 0) return true;
  }
  return false;
}

void MCPZProbeAnalyzer::beginJob() {
  edm::Service<TFileService> fs;
  tT_ = fs->make<TTree>("tracks", "Z probes (generalTracks)");
  tT_->Branch("run", &run_); tT_->Branch("lumi", &lumi_); tT_->Branch("event", &event_);
  tT_->Branch("nPV", &nPV_); tT_->Branch("passMET", &passMET_); tT_->Branch("passJet", &passJet_);
  tT_->Branch("passMuon", &passMuon_);
  tT_->Branch("pt", &pt_); tT_->Branch("eta", &eta_); tT_->Branch("phi", &phi_);
  tT_->Branch("charge", &charge_); tT_->Branch("highPurity", &highPurity_);
  tT_->Branch("nValidPixelHits", &nPixHit_); tT_->Branch("nTrackerLayers", &nTkLayers_);
  tT_->Branch("hasDeDx", &hasDeDx_); tT_->Branch("dedxPrescale", &dedxPrescale_);
  tT_->Branch("probQ_pixel", &probQ_); tT_->Branch("probQ_pixelNoL1", &probQNoL1_);
  tT_->Branch("pixSizeXresidual", &sizeXres_); tT_->Branch("pixSizeXpred", &sizeXpred_);
  tT_->Branch("pixClSizeX", &clSizeX_);
  tT_->Branch("ih_pixel", &ihPixel_); tT_->Branch("ih_strip", &ihStrip_);
  tT_->Branch("nPixHitsUsed", &nPixUsed_); tT_->Branch("nPixQFloor", &nPixQFloor_);
  tT_->Branch("nPixClusters", &nPixClusters_);
  tT_->Branch("nBPIX", &nBPIX_); tT_->Branch("nFPIX", &nFPIX_);
  tT_->Branch("muMatched", &muMatched_); tT_->Branch("muTight", &muTight_); tT_->Branch("muRelIso", &muRelIso_);
  tT_->Branch("tpMass", &tpMass_); tT_->Branch("tpTagPt", &tpTagPt_); tT_->Branch("tpOS", &tpOS_);
  tT_->Branch("tagTrigObj", &tagTrigObj_);

  tE_ = fs->make<TTree>("events", "per event");
  tE_->Branch("run", &run_); tE_->Branch("lumi", &lumi_); tE_->Branch("event", &event_);
  tE_->Branch("nPV", &nPV_); tE_->Branch("nTag", &nTag_);
  tE_->Branch("passMET", &passMET_); tE_->Branch("passJet", &passJet_); tE_->Branch("passMuon", &passMuon_);
}

void MCPZProbeAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  run_ = iEvent.run(); lumi_ = iEvent.luminosityBlock(); event_ = iEvent.id().event();

  // good PVs
  nPV_ = 0;
  const reco::Vertex* pv = nullptr;
  const auto vtxH = iEvent.getHandle(vertexToken_);
  if (vtxH.isValid()) {
    for (const auto& v : *vtxH) {
      if (v.isFake() || v.ndof() <= 4 || std::abs(v.z()) >= 24 || v.position().rho() >= 2) continue;
      if (!pv) pv = &v;
      ++nPV_;
    }
  }

  passMET_ = passJet_ = passMuon_ = 0;
  const auto trH = iEvent.getHandle(triggerResultsToken_);
  if (trH.isValid()) {
    const auto& names = iEvent.triggerNames(*trH);
    passMET_ = passAny(*trH, names, trigMET_);
    passJet_ = passAny(*trH, names, trigJet_);
    passMuon_ = passAny(*trH, names, trigMuon_);
  }

  // HLT IsoMu24 objects, if the filter is in the trigger summary
  std::vector<std::pair<float, float>> hltMu;  // (eta, phi)
  tagTrigObj_ = 0;
  const auto teH = iEvent.getHandle(trigEventToken_);
  if (teH.isValid()) {
    const unsigned int idx = teH->filterIndex(edm::InputTag(tagFilter_, "", "HLT"));
    if (idx < teH->sizeFilters()) {
      tagTrigObj_ = 1;
      for (auto k : teH->filterKeys(idx)) {
        const auto& o = teH->getObjects()[k];
        hltMu.emplace_back(o.eta(), o.phi());
      }
    }
  }

  // tags: tight ID, PF iso R04 < 0.15, pT>26, |eta|<2.4, matched to the HLT object (event bit if no summary)
  const auto muH = iEvent.getHandle(muonToken_);
  std::vector<const reco::Muon*> tags;
  if (muH.isValid() && pv) {
    for (const auto& mu : *muH) {
      if (mu.pt() < 26 || std::abs(mu.eta()) > 2.4 || !muon::isTightMuon(mu, *pv)) continue;
      const auto& iso = mu.pfIsolationR04();
      const float relIso = (iso.sumChargedHadronPt +
                            std::max(0.f, iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5f * iso.sumPUPt)) / mu.pt();
      if (relIso > 0.15f) continue;
      bool match = !tagTrigObj_ && passMuon_;
      for (const auto& o : hltMu)
        if (reco::deltaR(o.first, o.second, mu.eta(), mu.phi()) < 0.1) { match = true; break; }
      if (match) tags.push_back(&mu);
    }
  }
  std::sort(tags.begin(), tags.end(), [](auto a, auto b) { return a->pt() > b->pt(); });
  nTag_ = tags.size();
  tE_->Fill();
  if (tags.empty()) return;

  const auto trkH = iEvent.getHandle(trackToken_);
  if (!trkH.isValid()) return;
  const auto dedxH = iEvent.getHandle(dedxToken_);
  const auto dedxPsH = iEvent.getHandle(dedxPrescaleToken_);
  const TrackerTopology* tTopo = &iSetup.getData(topoToken_);
  const TrackerGeometry* tkGeom = &iSetup.getData(geomToken_);
  const PixelClusterParameterEstimator* cpe = &iSetup.getData(cpeToken_);

  for (unsigned int i = 0; i < trkH->size(); ++i) {
    const reco::Track& trk = (*trkH)[i];
    if (trk.pt() < probePtMin_ || std::abs(trk.eta()) > 2.4) continue;

    // T&P mass with the leading tag that is not this track
    const reco::Muon* tag = nullptr;
    for (const auto* t : tags) {
      if (t->innerTrack().isNonnull() && t->innerTrack().key() == i) continue;
      tag = t; break;
    }
    if (!tag) continue;
    reco::Candidate::PolarLorentzVector pTrk(trk.pt(), trk.eta(), trk.phi(), 0.1057);
    tpMass_ = (tag->polarP4() + pTrk).mass();
    if (tpMass_ < massMin_ || tpMass_ > massMax_) continue;
    tpTagPt_ = tag->pt();
    tpOS_ = (tag->charge() * trk.charge() < 0);

    // probe muon match via the inner track
    muMatched_ = muTight_ = 0; muRelIso_ = -1.f;
    for (const auto& mu : *muH) {
      if (mu.innerTrack().isNull() || mu.innerTrack().key() != i) continue;
      muMatched_ = 1;
      muTight_ = muon::isTightMuon(mu, *pv);
      const auto& iso = mu.pfIsolationR04();
      muRelIso_ = (iso.sumChargedHadronPt +
                   std::max(0.f, iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5f * iso.sumPUPt)) / mu.pt();
      break;
    }

    pt_ = trk.pt(); eta_ = trk.eta(); phi_ = trk.phi(); charge_ = trk.charge();
    highPurity_ = trk.quality(reco::TrackBase::highPurity);
    nPixHit_ = trk.hitPattern().numberOfValidPixelHits();
    nTkLayers_ = trk.hitPattern().trackerLayersWithMeasurement();

    const reco::DeDxHitInfo* hits = nullptr;
    reco::TrackRef ref(trkH, i);
    if (dedxH.isValid()) {
      reco::DeDxHitInfoRef hRef = (*dedxH)[ref];
      if (hRef.isNonnull()) hits = &(*hRef);
    }
    hasDeDx_ = hits != nullptr;
    dedxPrescale_ = (hits && dedxPsH.isValid()) ? (*dedxPsH)[(*dedxH)[ref]] : 0;

    const MCPDeDxResult r = MCPProbQ::compute(hits, trk.px(), trk.py(), trk.pz(), trk.charge(), tkGeom, tTopo, cpe);
    probQ_ = r.probQonTrack; probQNoL1_ = r.probQonTrackNoL1;
    sizeXres_ = r.pixSizeXresidualMean; sizeXpred_ = r.pixSizeXpredMean; clSizeX_ = r.pixClustSizeXMean;
    ihPixel_ = r.ih_pixel; ihStrip_ = r.ih_strip;
    nPixUsed_ = r.nPixHitsUsed; nPixQFloor_ = r.nPixQFloor; nPixClusters_ = r.nPixClusters;
    nBPIX_ = r.nBPIX; nFPIX_ = r.nFPIX;
    tT_->Fill();
  }
}

void MCPZProbeAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("dedxHitInfo", edm::InputTag("dedxHitInfo"));
  desc.add<edm::InputTag>("dedxHitInfoPrescale", edm::InputTag("dedxHitInfo", "prescale"));
  desc.add<edm::InputTag>("muons", edm::InputTag("muons"));
  desc.add<edm::InputTag>("primaryVertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("triggerResults", edm::InputTag("TriggerResults", "", "HLT"));
  desc.add<edm::InputTag>("triggerEvent", edm::InputTag("hltTriggerSummaryAOD", "", "HLT"));
  desc.add<std::string>("tagFilter", "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered");
  desc.add<std::string>("pixelCPE", "PixelCPETemplateReco");
  desc.add<double>("probePtMin", 10.);
  desc.add<double>("massMin", 50.);
  desc.add<double>("massMax", 130.);
  descriptions.add("MCPZProbeAnalyzer", desc);
}

DEFINE_FWK_MODULE(MCPZProbeAnalyzer);
