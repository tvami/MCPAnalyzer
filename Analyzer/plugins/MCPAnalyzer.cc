// -*- C++ -*-
//
// Package:    MCPAnalyzer/Analyzer
// Class:      MCPAnalyzer
//
// Characterize HSCP multi-charged particles (MCP) in MiniAOD:
//   (1) charge via gen-pT / reco-pT (curvature assumes |q|=1)
//   (2) discriminators at high charge: pixel probQ, pixel/strip dE/dx,
//       and whether the track's hits reach the strip tracker.
//
// Self-contained: reads the DeDxHitInfo association attached to pat::IsolatedTrack
// directly (no HSCParticleProducer), and uses our own MCPProbQ helper.

#include <memory>
#include <vector>
#include <cmath>
#include <fstream>
#include <set>
#include <map>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "RecoLocalTracker/Records/interface/TkPixelCPERecord.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"

#include "TTree.h"

#include "MCPProbQ.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"



class MCPAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MCPAnalyzer(const edm::ParameterSet&);
  ~MCPAnalyzer() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  // tokens

  //pfmet token: 
  const edm::EDGetTokenT<pat::METCollection> metToken_;
  const edm::EDGetTokenT<pat::METCollection> puppiMetToken_;


  const edm::EDGetTokenT<std::vector<pat::IsolatedTrack>> trackToken_;
  const edm::EDGetTokenT<reco::DeDxHitInfoAss> dedxToken_;
  const edm::EDGetTokenT<std::vector<reco::GenParticle>> prunedGenToken_;
  const edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const std::string pixelCPEName_;
  const edm::ESGetToken<PixelClusterParameterEstimator, TkPixelCPERecord> cpeToken_;
  const int mcpPdgId_;
  const edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  const edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  const edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigObjToken_;
  const bool saveTrigNames_;
  const bool tpOnly_;

  // HLT groups, matched as "<name>_v" prefixes
  static const std::vector<std::string> kTrigOR_, kTrigMET_, kTrigJet_, kTrigTau_, kTrigMuon_;
  static bool passAny(const edm::TriggerResults&, const edm::TriggerNames&, const std::vector<std::string>&);


  // per-track tree
  TTree* tT_ = nullptr;
  // per-gen-MCP tree
  TTree* tG_ = nullptr;
  // per-event tree
  TTree* tE_ = nullptr;

  // track-tree branches
  float b_pfMET_;
  float b_puppiMET_;
  std::vector<std::string> b_trigNames_;
  std::vector<int> b_trigPass_;
  int b_passTrigger_OR;
  int b_passMET_, b_passJet_, b_passTau_, b_passMuon_;
  int b_nPV_;
  // muon tag-and-probe
  int b_muMatched_, b_muTight_, b_muIsTag_; float b_muRelIso_, b_muPt_;
  float b_tpMass_, b_tpTagPt_; int b_tpOS_;

  unsigned int b_run_, b_lumi_; unsigned long long b_event_;
  double b_pt_, b_eta_, b_phi_, b_ptError_, b_normChi2_, b_validFrac_;
  int b_charge_, b_nPixHit_, b_nTkLayers_, b_highPurity_;
  float b_caloEmEnergy_, b_caloHadEnergy_;  // matched calo-jet EM/HAD energy along the track
  float b_pcCaloFrac_, b_pcHcalFrac_;       // packed-candidate calo fractions
  float b_dedxStripBuiltin_, b_dedxPixelBuiltin_;
  float b_probQpixel_, b_probQpixelNoL1_, b_probXYpixel_;
  std::vector<float> b_pixelDedxHits_, b_stripDedxHits_;
  int b_nPixUsed_, b_nonL1Pix_;
  int b_nPixQFloor_;
  int b_nPixClusters_, b_nPixNoFillProb_, b_nPixSpecInCPE_, b_nPixXYpinnedLo_, b_nPixXYpinnedHi_, b_nPixXYvalid_;
  float b_pixXYrawMin_;
  float b_ihFull_, b_ihPixel_, b_ihStrip_;
  int b_nomFull_, b_nomPixel_, b_nomStrip_;
  int b_nBPIX_, b_nFPIX_, b_nTIB_, b_nTID_, b_nTOB_, b_nTEC_, b_nPix_, b_nStrip_;
  float b_pixClSize_, b_pixClSizeX_, b_pixClSizeY_, b_pixClCharge_; int b_pixClSizeMax_;
  float b_pixSizeXpred_, b_pixSizeXresidual_;  // angle-predicted r-phi size and (measured-predicted)
  float b_strClWidth_, b_strClCharge_; int b_strClWidthMax_;
  int b_hasDeDx_;
  int b_genMatched_; double b_genPt_, b_genEta_; int b_genCharge_, b_genPdgId_; double b_dR_;
  double b_chargeFromCurv_;

  // gen-tree branches
  unsigned int g_run_, g_lumi_; unsigned long long g_event_;
  double g_pt_, g_eta_, g_phi_; int g_charge_, g_pdgId_;
  int g_matched_; double g_recoPt_, g_dR_, g_chargeFromCurv_;
  int g_hasDeDx_; float g_probQpixel_, g_sizeXresidual_;  // of the best-matched track

  //event-tree branches
  unsigned int e_run_, e_lumi_; unsigned long long e_event_;
  float e_pfMET_;
  float e_puppiMET_;
  std::vector<std::string> e_trigNames_;
  std::vector<int> e_trigPass_;
  int e_passTrigger_OR;
  int e_passMET_, e_passJet_, e_passTau_, e_passMuon_;
  int e_nPV_, e_nTag_;

};

MCPAnalyzer::MCPAnalyzer(const edm::ParameterSet& iC)
    : metToken_(consumes<pat::METCollection>(iC.getParameter<edm::InputTag>("slimmedMET"))),
      puppiMetToken_(consumes<pat::METCollection>(iC.getParameter<edm::InputTag>("slimmedPuppiMET"))),

      trackToken_(consumes<std::vector<pat::IsolatedTrack>>(iC.getParameter<edm::InputTag>("isolatedTracks"))),
      dedxToken_(consumes<reco::DeDxHitInfoAss>(iC.getParameter<edm::InputTag>("dedxHitInfo"))),
      prunedGenToken_(consumes<std::vector<reco::GenParticle>>(iC.getParameter<edm::InputTag>("prunedGenParticles"))),
      vertexToken_(consumes<std::vector<reco::Vertex>>(iC.getParameter<edm::InputTag>("primaryVertices"))),
      topoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>()),
      geomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()),
      pixelCPEName_(iC.getParameter<std::string>("pixelCPE")),
      cpeToken_(esConsumes<PixelClusterParameterEstimator, TkPixelCPERecord>(edm::ESInputTag("", pixelCPEName_))),
      mcpPdgId_(iC.getParameter<int>("mcpPdgId")),
      triggerResultsToken_(consumes<edm::TriggerResults>(iC.getParameter<edm::InputTag>("triggerResults"))),
      muonToken_(consumes<std::vector<pat::Muon>>(iC.getParameter<edm::InputTag>("muons"))),
      trigObjToken_(consumes<pat::TriggerObjectStandAloneCollection>(iC.getParameter<edm::InputTag>("triggerObjects"))),
      saveTrigNames_(iC.getParameter<bool>("saveTrigNames")),
      tpOnly_(iC.getParameter<bool>("tpOnly")) {
  usesResource("TFileService");
}

// old OR list (MET + ditau), kept for continuity
const std::vector<std::string> MCPAnalyzer::kTrigOR_ = {
    "HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng", "HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1",
    "HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3", "HLT_DoublePNetTauhPFJet30_Tight_L2NN_eta2p3",
    "HLT_MET105_IsoTrk50", "HLT_MET120_IsoTrk50", "HLT_PFMET105_IsoTrk50",
    "HLT_PFMET120_PFMHT120_IDTight", "HLT_PFMET120_PFMHT120_IDTight_PFHT60", "HLT_PFMET130_PFMHT130_IDTight",
    "HLT_PFMET140_PFMHT140_IDTight", "HLT_PFMET200_BeamHaloCleaned", "HLT_PFMET250_NotCleaned",
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",  // NoMu110_FilterHF dropped: disabled from 2024F
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
    "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF",
    "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF",
    "HLT_PFMETTypeOne140_PFMHT140_IDTight", "HLT_PFMETTypeOne200_BeamHaloCleaned"};
const std::vector<std::string> MCPAnalyzer::kTrigMET_ = {
    "HLT_MET105_IsoTrk50", "HLT_MET120_IsoTrk50", "HLT_PFMET105_IsoTrk50",
    "HLT_PFMET120_PFMHT120_IDTight", "HLT_PFMET120_PFMHT120_IDTight_PFHT60", "HLT_PFMET130_PFMHT130_IDTight",
    "HLT_PFMET140_PFMHT140_IDTight", "HLT_PFMET200_BeamHaloCleaned", "HLT_PFMET250_NotCleaned",
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",  // NoMu110_FilterHF dropped: disabled from 2024F
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
    "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF",
    "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF",
    "HLT_PFMETTypeOne140_PFMHT140_IDTight", "HLT_PFMETTypeOne200_BeamHaloCleaned", "HLT_CaloMET350_NotCleaned"};
const std::vector<std::string> MCPAnalyzer::kTrigJet_ = {
    "HLT_PFJet500", "HLT_AK8PFJet500", "HLT_PFHT1050", "HLT_CaloJet500_NoJetID"};
const std::vector<std::string> MCPAnalyzer::kTrigTau_ = {
    "HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng", "HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1",
    "HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3", "HLT_DoublePNetTauhPFJet30_Tight_L2NN_eta2p3"};
const std::vector<std::string> MCPAnalyzer::kTrigMuon_ = {"HLT_IsoMu24", "HLT_Mu50"};

bool MCPAnalyzer::passAny(const edm::TriggerResults& tr, const edm::TriggerNames& names,
                          const std::vector<std::string>& list) {
  for (unsigned int i = 0; i < tr.size(); ++i) {
    if (!tr.accept(i)) continue;
    const std::string& n = names.triggerName(i);
    for (const auto& t : list)
      if (n.compare(0, t.size() + 2, t + "_v") == 0) return true;
  }
  return false;
}


void MCPAnalyzer::beginJob() {
  edm::Service<TFileService> fs;
  std::cout << "MCPAnalyzer new build loaded\n";
  tT_ = fs->make<TTree>("tracks", "per isolated-track");
  tT_->Branch("run", &b_run_);            tT_->Branch("lumi", &b_lumi_);   tT_->Branch("event", &b_event_);
  tT_->Branch("pt", &b_pt_);              tT_->Branch("eta", &b_eta_);     tT_->Branch("phi", &b_phi_);
  tT_->Branch("ptError", &b_ptError_);    tT_->Branch("normChi2", &b_normChi2_);
  tT_->Branch("validFraction", &b_validFrac_);
  tT_->Branch("charge", &b_charge_);      tT_->Branch("nValidPixelHits", &b_nPixHit_);
  tT_->Branch("nTrackerLayers", &b_nTkLayers_); tT_->Branch("highPurity", &b_highPurity_);
  tT_->Branch("caloEmEnergy", &b_caloEmEnergy_); tT_->Branch("caloHadEnergy", &b_caloHadEnergy_);
  tT_->Branch("pcCaloFrac", &b_pcCaloFrac_); tT_->Branch("pcHcalFrac", &b_pcHcalFrac_);
  tT_->Branch("dedxStrip_builtin", &b_dedxStripBuiltin_); tT_->Branch("dedxPixel_builtin", &b_dedxPixelBuiltin_);
  tT_->Branch("probQ_pixel", &b_probQpixel_); tT_->Branch("probQ_pixelNoL1", &b_probQpixelNoL1_);
  tT_->Branch("probXY_pixel", &b_probXYpixel_);
  tT_->Branch("pixelDedxHits", &b_pixelDedxHits_); tT_->Branch("stripDedxHits", &b_stripDedxHits_);
  tT_->Branch("nPixHitsUsed", &b_nPixUsed_); tT_->Branch("nonL1PixHits", &b_nonL1Pix_);
  tT_->Branch("nPixClusters", &b_nPixClusters_); tT_->Branch("nPixNoFillProb", &b_nPixNoFillProb_);
  tT_->Branch("nPixQFloor", &b_nPixQFloor_);
  tT_->Branch("nPixSpecInCPE", &b_nPixSpecInCPE_); tT_->Branch("nPixXYpinnedLo", &b_nPixXYpinnedLo_);
  tT_->Branch("nPixXYpinnedHi", &b_nPixXYpinnedHi_); tT_->Branch("nPixXYvalid", &b_nPixXYvalid_);
  tT_->Branch("pixXYrawMin", &b_pixXYrawMin_);
  tT_->Branch("ih_full", &b_ihFull_);     tT_->Branch("ih_pixel", &b_ihPixel_); tT_->Branch("ih_strip", &b_ihStrip_);
  tT_->Branch("nom_full", &b_nomFull_);   tT_->Branch("nom_pixel", &b_nomPixel_); tT_->Branch("nom_strip", &b_nomStrip_);
  tT_->Branch("nBPIX", &b_nBPIX_); tT_->Branch("nFPIX", &b_nFPIX_);
  tT_->Branch("nTIB", &b_nTIB_); tT_->Branch("nTID", &b_nTID_);
  tT_->Branch("nTOB", &b_nTOB_); tT_->Branch("nTEC", &b_nTEC_);
  tT_->Branch("nPix", &b_nPix_); tT_->Branch("nStrip", &b_nStrip_);
  tT_->Branch("pixClSize", &b_pixClSize_); tT_->Branch("pixClSizeX", &b_pixClSizeX_);
  tT_->Branch("pixClSizeY", &b_pixClSizeY_); tT_->Branch("pixClSizeMax", &b_pixClSizeMax_);
  tT_->Branch("pixClCharge", &b_pixClCharge_);
  tT_->Branch("pixSizeXpred", &b_pixSizeXpred_); tT_->Branch("pixSizeXresidual", &b_pixSizeXresidual_);
  tT_->Branch("strClWidth", &b_strClWidth_); tT_->Branch("strClWidthMax", &b_strClWidthMax_);
  tT_->Branch("strClCharge", &b_strClCharge_);
  tT_->Branch("hasDeDx", &b_hasDeDx_);
  tT_->Branch("genMatched", &b_genMatched_); tT_->Branch("gen_pt", &b_genPt_); tT_->Branch("gen_eta", &b_genEta_);
  tT_->Branch("gen_charge", &b_genCharge_); tT_->Branch("gen_pdgId", &b_genPdgId_); tT_->Branch("dR", &b_dR_);
  tT_->Branch("chargeFromCurvature", &b_chargeFromCurv_);
  
  tT_->Branch("pfMET", &b_pfMET_);
  tT_->Branch("puppiMET", &b_puppiMET_);
  if (saveTrigNames_) { tT_->Branch("trigNames", &b_trigNames_); tT_->Branch("trigPass", &b_trigPass_); }
  tT_->Branch("HLT_trigPass_OR", &b_passTrigger_OR)->SetTitle("OR of MET + ditau paths");
  tT_->Branch("passMET", &b_passMET_); tT_->Branch("passJet", &b_passJet_);
  tT_->Branch("passTau", &b_passTau_); tT_->Branch("passMuon", &b_passMuon_);
  tT_->Branch("nPV", &b_nPV_);
  tT_->Branch("muMatched", &b_muMatched_); tT_->Branch("muTight", &b_muTight_);
  tT_->Branch("muIsTag", &b_muIsTag_); tT_->Branch("muRelIso", &b_muRelIso_); tT_->Branch("muPt", &b_muPt_);
  tT_->Branch("tpMass", &b_tpMass_); tT_->Branch("tpTagPt", &b_tpTagPt_); tT_->Branch("tpOS", &b_tpOS_);


  tE_ = fs->make<TTree>("events", "per event");
  tE_->Branch("run", &e_run_); tE_->Branch("lumi", &e_lumi_); tE_->Branch("event", &e_event_);
  tE_->Branch("pfMET", &e_pfMET_);
  tE_->Branch("puppiMET", &e_puppiMET_);
  if (saveTrigNames_) { tE_->Branch("trigNames", &e_trigNames_); tE_->Branch("trigPass", &e_trigPass_); }
  tE_->Branch("HLT_trigPass_OR", &e_passTrigger_OR)->SetTitle("OR of MET + ditau paths");
  tE_->Branch("passMET", &e_passMET_); tE_->Branch("passJet", &e_passJet_);
  tE_->Branch("passTau", &e_passTau_); tE_->Branch("passMuon", &e_passMuon_);
  tE_->Branch("nPV", &e_nPV_); tE_->Branch("nTag", &e_nTag_);
  

  tG_ = fs->make<TTree>("gen", "per gen MCP");
  tG_->Branch("run", &g_run_); tG_->Branch("lumi", &g_lumi_); tG_->Branch("event", &g_event_);
  tG_->Branch("gen_pt", &g_pt_); tG_->Branch("gen_eta", &g_eta_); tG_->Branch("gen_phi", &g_phi_);
  tG_->Branch("gen_charge", &g_charge_); tG_->Branch("gen_pdgId", &g_pdgId_);
  tG_->Branch("matched", &g_matched_); tG_->Branch("reco_pt", &g_recoPt_); tG_->Branch("dR", &g_dR_);
  tG_->Branch("chargeFromCurvature", &g_chargeFromCurv_);
  tG_->Branch("hasDeDx", &g_hasDeDx_); tG_->Branch("probQ_pixel", &g_probQpixel_);
  tG_->Branch("pixSizeXresidual", &g_sizeXresidual_);
  tG_->Branch("nPV", &e_nPV_);
  
}

void MCPAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const TrackerTopology* tTopo = &iSetup.getData(topoToken_);
  const TrackerGeometry* tkGeom = &iSetup.getData(geomToken_);
  const PixelClusterParameterEstimator* pixelCPE = &iSetup.getData(cpeToken_);

  edm::Handle<std::vector<pat::IsolatedTrack>> tracks = iEvent.getHandle(trackToken_);
  edm::Handle<reco::DeDxHitInfoAss> dedxAss = iEvent.getHandle(dedxToken_);
  edm::Handle<std::vector<reco::GenParticle>> pruned = iEvent.getHandle(prunedGenToken_);

  const unsigned int run = iEvent.run();
  const unsigned int lumi = iEvent.luminosityBlock();
  const unsigned long long event = iEvent.id().event();

  // collect stable gen MCPs
  std::vector<const reco::GenParticle*> mcps;
  if (pruned.isValid()) {
    for (const auto& g : *pruned) {
      if (std::abs(g.pdgId()) == mcpPdgId_ && g.status() == 1) mcps.push_back(&g);
    }
  }

  // invalid MET -> -1, do not drop the event
  const auto pfMETCollection = iEvent.getHandle(metToken_);
  const auto puppiMETCollection = iEvent.getHandle(puppiMetToken_);
  b_pfMET_ = (pfMETCollection.isValid() && !pfMETCollection->empty()) ? pfMETCollection->front().pt() : -1.f;
  b_puppiMET_ = (puppiMETCollection.isValid() && !puppiMETCollection->empty()) ? puppiMETCollection->front().pt() : -1.f;

  // good PVs
  b_nPV_ = 0;
  const auto vtxH = iEvent.getHandle(vertexToken_);
  const reco::Vertex* pv = nullptr;
  if (vtxH.isValid()) {
    for (const auto& v : *vtxH) {
      if (v.isFake() || v.ndof() <= 4 || std::abs(v.z()) >= 24 || v.position().rho() >= 2) continue;
      if (!pv) pv = &v;
      ++b_nPV_;
    }
  }

  b_trigNames_.clear();
  b_trigPass_.clear();
  b_passTrigger_OR = b_passMET_ = b_passJet_ = b_passTau_ = b_passMuon_ = 0;

  const auto triggerH = iEvent.getHandle(triggerResultsToken_);
  const edm::TriggerNames* trigNamesPtr = nullptr;
  if (triggerH.isValid()) {
    const auto& triggerNames = iEvent.triggerNames(*triggerH);
    trigNamesPtr = &triggerNames;
    if (saveTrigNames_) {
      for (unsigned int i = 0; i < triggerH->size(); ++i) {
        b_trigNames_.push_back(triggerNames.triggerName(i));
        b_trigPass_.push_back(triggerH->accept(i) ? 1 : 0);
      }
    }
    b_passTrigger_OR = passAny(*triggerH, triggerNames, kTrigOR_);
    b_passMET_ = passAny(*triggerH, triggerNames, kTrigMET_);
    b_passJet_ = passAny(*triggerH, triggerNames, kTrigJet_);
    b_passTau_ = passAny(*triggerH, triggerNames, kTrigTau_);
    b_passMuon_ = passAny(*triggerH, triggerNames, kTrigMuon_);
  }

  // HLT muon objects (IsoMu24) for tag matching
  std::vector<const pat::TriggerObjectStandAlone*> hltMu;
  std::vector<pat::TriggerObjectStandAlone> trigObjs;
  const auto trigObjH = iEvent.getHandle(trigObjToken_);
  if (trigObjH.isValid() && trigNamesPtr) {
    trigObjs = *trigObjH;
    for (auto& o : trigObjs) o.unpackPathNames(*trigNamesPtr);
    for (const auto& o : trigObjs)
      if (o.hasPathName("HLT_IsoMu24_v*", true, true)) hltMu.push_back(&o);
  }

  // muons: tag = tight ID + tight PF iso, pT>26, |eta|<2.4, matched to IsoMu24 object
  std::vector<const pat::Muon*> muons, tags;
  const auto muH = iEvent.getHandle(muonToken_);
  if (muH.isValid()) {
    for (const auto& mu : *muH) {
      muons.push_back(&mu);
      if (mu.pt() < 26 || std::abs(mu.eta()) > 2.4) continue;
      if (!mu.passed(reco::Muon::CutBasedIdTight) || !mu.passed(reco::Muon::PFIsoTight)) continue;
      bool trigMatch = false;
      for (const auto* o : hltMu)
        if (reco::deltaR(*o, mu) < 0.1) { trigMatch = true; break; }
      if (trigMatch) tags.push_back(&mu);
    }
  }


  //filling event level tree outside of track loop
  e_run_ = run; e_lumi_ = lumi; e_event_ = event;
  e_pfMET_ = b_pfMET_;
  e_puppiMET_ = b_puppiMET_;
  e_trigNames_ = b_trigNames_;
  e_trigPass_ = b_trigPass_;
  e_passTrigger_OR = b_passTrigger_OR;
  e_passMET_ = b_passMET_; e_passJet_ = b_passJet_; e_passTau_ = b_passTau_; e_passMuon_ = b_passMuon_;
  e_nPV_ = b_nPV_; e_nTag_ = tags.size();
  tE_->Fill();
  
  // for gen-tree: track best reco match per MCP
  std::vector<double> mcpBestDR(mcps.size(), 1e9);
  std::vector<double> mcpBestRecoPt(mcps.size(), -1.);
  std::vector<int> mcpBestHasDeDx(mcps.size(), 0);
  std::vector<float> mcpBestProbQ(mcps.size(), -1.f), mcpBestSizeX(mcps.size(), -99.f);

  // ---- per isolated-track loop ----
  if (tracks.isValid()) {
    for (unsigned int i = 0; i < tracks->size(); ++i) {
      const pat::IsolatedTrack& it = (*tracks)[i];

      b_run_ = run; b_lumi_ = lumi; b_event_ = event;
      // pat::IsolatedTrack kinematics are the track-fit (curvature) values, built
      // assuming |q|=1 -> reco pt is gen pt / |Q_true|, which is the charge probe.
      b_pt_ = it.pt(); b_eta_ = it.eta(); b_phi_ = it.phi();
      b_charge_ = it.charge();
      b_nPixHit_ = it.hitPattern().numberOfValidPixelHits();
      b_nTkLayers_ = it.hitPattern().trackerLayersWithMeasurement();
      b_highPurity_ = it.isHighPurityTrack() ? 1 : 0;
      b_dedxStripBuiltin_ = it.dEdxStrip();
      b_dedxPixelBuiltin_ = it.dEdxPixel();
      b_caloEmEnergy_ = it.matchedCaloJetEmEnergy();
      b_caloHadEnergy_ = it.matchedCaloJetHadEnergy();
      // packed-candidate calo fractions (track-associated calo energy fraction)
      b_pcCaloFrac_ = -1.; b_pcHcalFrac_ = -1.;
      if (it.packedCandRef().isNonnull()) {
        b_pcCaloFrac_ = it.packedCandRef()->caloFraction();
        b_pcHcalFrac_ = it.packedCandRef()->hcalFraction();
      }
      // ptError / normChi2 / validFraction need the underlying track (via packedCandRef)
      b_ptError_ = -1.; b_normChi2_ = -1.; b_validFrac_ = -1.;
      if (it.packedCandRef().isNonnull() && it.packedCandRef()->hasTrackDetails()) {
        const reco::Track& trk = it.packedCandRef()->pseudoTrack();
        b_ptError_ = trk.ptError(); b_normChi2_ = trk.normalizedChi2(); b_validFrac_ = trk.validFraction();
      }

      // muon match (dR<0.02) and T&P mass with the leading other tag
      b_muMatched_ = b_muTight_ = b_muIsTag_ = 0; b_muRelIso_ = -1.f; b_muPt_ = -1.f;
      b_tpMass_ = -1.f; b_tpTagPt_ = -1.f; b_tpOS_ = 0;
      const pat::Muon* muMatch = nullptr; double muBestDR = 0.02;
      for (const auto* mu : muons) {
        const double dr = reco::deltaR(*mu, it);
        if (dr < muBestDR) { muBestDR = dr; muMatch = mu; }
      }
      if (muMatch) {
        b_muMatched_ = 1; b_muPt_ = muMatch->pt();
        b_muTight_ = muMatch->passed(reco::Muon::CutBasedIdTight);
        const auto& iso = muMatch->pfIsolationR04();
        b_muRelIso_ = (iso.sumChargedHadronPt +
                       std::max(0.f, iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5f * iso.sumPUPt)) / muMatch->pt();
        b_muIsTag_ = std::find(tags.begin(), tags.end(), muMatch) != tags.end();
      }
      for (const auto* tag : tags) {  // tags are pT-ordered
        if (tag == muMatch || reco::deltaR(*tag, it) < 0.02) continue;
        b_tpMass_ = (tag->p4() + it.p4()).mass();
        b_tpTagPt_ = tag->pt();
        b_tpOS_ = (tag->charge() * it.charge() < 0);
        break;
      }
      if (tpOnly_ && !(b_tpMass_ > 50 && b_tpMass_ < 130)) continue;  // Z probes only

      // dE/dx hits via the isolatedTracks association
      const reco::DeDxHitInfo* hits = nullptr;
      if (dedxAss.isValid()) {
        edm::Ref<std::vector<pat::IsolatedTrack>> ref(tracks, i);
        reco::DeDxHitInfoRef hitRef = (*dedxAss)[ref];
        if (hitRef.isNonnull()) hits = &(*hitRef);
      }
      b_hasDeDx_ = (hits != nullptr) ? 1 : 0;

      MCPDeDxResult r = MCPProbQ::compute(hits, it.px(), it.py(), it.pz(), it.charge(),
                                          tkGeom, tTopo, pixelCPE);
      b_probQpixel_ = r.probQonTrack; b_probQpixelNoL1_ = r.probQonTrackNoL1; b_probXYpixel_ = r.probXYonTrack;
      b_pixelDedxHits_ = r.pixelDedxHits; b_stripDedxHits_ = r.stripDedxHits;
      b_nPixUsed_ = r.nPixHitsUsed; b_nonL1Pix_ = r.nonL1PixHits;
      b_nPixClusters_ = r.nPixClusters; b_nPixNoFillProb_ = r.nPixNoFillProb; b_nPixQFloor_ = r.nPixQFloor;
      b_nPixSpecInCPE_ = r.nPixSpecInCPE; b_nPixXYpinnedLo_ = r.nPixXYpinnedLo;
      b_nPixXYpinnedHi_ = r.nPixXYpinnedHi; b_nPixXYvalid_ = r.nPixXYvalid;
      b_pixXYrawMin_ = r.pixXYrawMin;
      b_ihFull_ = r.ih_full; b_ihPixel_ = r.ih_pixel; b_ihStrip_ = r.ih_strip;
      b_nomFull_ = r.nom_full; b_nomPixel_ = r.nom_pixel; b_nomStrip_ = r.nom_strip;
      b_nBPIX_ = r.nBPIX; b_nFPIX_ = r.nFPIX; b_nTIB_ = r.nTIB; b_nTID_ = r.nTID;
      b_nTOB_ = r.nTOB; b_nTEC_ = r.nTEC; b_nPix_ = r.nPix; b_nStrip_ = r.nStrip;
      b_pixClSize_ = r.pixClustSizeMean; b_pixClSizeX_ = r.pixClustSizeXMean;
      b_pixClSizeY_ = r.pixClustSizeYMean; b_pixClSizeMax_ = r.pixClustSizeMax;
      b_pixClCharge_ = r.pixClustChargeMean;
      b_pixSizeXpred_ = r.pixSizeXpredMean; b_pixSizeXresidual_ = r.pixSizeXresidualMean;
      b_strClWidth_ = r.stripClustWidthMean; b_strClWidthMax_ = r.stripClustWidthMax;
      b_strClCharge_ = r.stripClustChargeMean;

      // gen match (nearest MCP by dR)
      b_genMatched_ = 0; b_genPt_ = -1; b_genEta_ = 0; b_genCharge_ = 0; b_genPdgId_ = 0; b_dR_ = -1;
      double bestDR = 1e9; int bestIdx = -1;
      for (unsigned int m = 0; m < mcps.size(); ++m) {
        double dr = reco::deltaR(it.eta(), it.phi(), mcps[m]->eta(), mcps[m]->phi());
        if (dr < bestDR) { bestDR = dr; bestIdx = m; }
      }
      if (bestIdx >= 0 && bestDR < 0.1) {
        const reco::GenParticle* g = mcps[bestIdx];
        b_genMatched_ = 1; b_genPt_ = g->pt(); b_genEta_ = g->eta();
        b_genCharge_ = g->charge(); b_genPdgId_ = g->pdgId(); b_dR_ = bestDR;
        b_chargeFromCurv_ = (b_pt_ > 0) ? b_genPt_ / b_pt_ : -1.;
        if (bestDR < mcpBestDR[bestIdx]) {
          mcpBestDR[bestIdx] = bestDR; mcpBestRecoPt[bestIdx] = b_pt_;
          mcpBestHasDeDx[bestIdx] = b_hasDeDx_; mcpBestProbQ[bestIdx] = b_probQpixel_;
          mcpBestSizeX[bestIdx] = b_pixSizeXresidual_;
        }

      } else {
        b_chargeFromCurv_ = -1.;
      }

      tT_->Fill();
    }
  }

  // ---- per gen-MCP loop (efficiency / clean charge study) ----
  for (unsigned int m = 0; m < mcps.size(); ++m) {
    const reco::GenParticle* g = mcps[m];
    g_run_ = run; g_lumi_ = lumi; g_event_ = event;
    g_pt_ = g->pt(); g_eta_ = g->eta(); g_phi_ = g->phi();
    g_charge_ = g->charge(); g_pdgId_ = g->pdgId();
    g_matched_ = (mcpBestRecoPt[m] > 0) ? 1 : 0;
    g_recoPt_ = mcpBestRecoPt[m];
    g_dR_ = (mcpBestDR[m] < 1e9) ? mcpBestDR[m] : -1.;
    g_chargeFromCurv_ = (g_matched_ && g_recoPt_ > 0) ? g_pt_ / g_recoPt_ : -1.;
    g_hasDeDx_ = mcpBestHasDeDx[m]; g_probQpixel_ = mcpBestProbQ[m]; g_sizeXresidual_ = mcpBestSizeX[m];
    tG_->Fill();
  }
}

void MCPAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("slimmedMET", edm::InputTag("slimmedMETs"));
  desc.add<edm::InputTag>("slimmedPuppiMET", edm::InputTag("slimmedMETsPuppi"));
  desc.add<edm::InputTag>("isolatedTracks", edm::InputTag("isolatedTracks"));
  desc.add<edm::InputTag>("dedxHitInfo", edm::InputTag("isolatedTracks"));
  desc.add<edm::InputTag>("prunedGenParticles", edm::InputTag("prunedGenParticles"));
  desc.add<edm::InputTag>("primaryVertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<std::string>("pixelCPE", "PixelCPETemplateReco");
  desc.add<int>("mcpPdgId", 10000200);
  desc.add<edm::InputTag>("triggerResults", edm::InputTag("TriggerResults", "", "HLT"));
  desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
  desc.add<edm::InputTag>("triggerObjects", edm::InputTag("slimmedPatTrigger"));
  desc.add<bool>("saveTrigNames", true);
  desc.add<bool>("tpOnly", false);
  descriptions.add("MCPAnalyzer", desc);
}

DEFINE_FWK_MODULE(MCPAnalyzer);
