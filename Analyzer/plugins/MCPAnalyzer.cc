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

  std::vector<std::string> b_trigNames_;
  std::vector<int> b_trigPass_;
  // per-track tree
  TTree* tT_ = nullptr;
  // per-gen-MCP tree
  TTree* tG_ = nullptr;

  // track-tree branches
  float b_pfMET_;
  float b_puppiMET_;

  unsigned int b_run_, b_lumi_; unsigned long long b_event_;
  double b_pt_, b_eta_, b_phi_, b_ptError_, b_normChi2_, b_validFrac_;
  int b_charge_, b_nPixHit_, b_nTkLayers_, b_highPurity_;
  float b_caloEmEnergy_, b_caloHadEnergy_;  // matched calo-jet EM/HAD energy along the track
  float b_pcCaloFrac_, b_pcHcalFrac_;       // packed-candidate calo fractions
  float b_dedxStripBuiltin_, b_dedxPixelBuiltin_;
  float b_probQpixel_, b_probQpixelNoL1_, b_probXYpixel_;
  std::vector<float> b_pixelDedxHits_, b_stripDedxHits_;
  int b_nPixUsed_, b_nonL1Pix_;
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

  std::vector<int> b_passTrigger_OR; 

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
      triggerResultsToken_(consumes<edm::TriggerResults>(iC.getParameter<edm::InputTag>("TriggerResults"))) {
  usesResource("TFileService");
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
  
  tT_->Branch("trigNames", &b_trigNames_);
  tT_->Branch("trigPass", &b_trigPass_);
  tT_->Branch("HLT_trigPass_OR", &b_passTrigger_OR);
  TBranch*br = tT_->GetBranch("HLT_trigPass_OR");
  br->SetTitle("OR_HLT_non-prescaled_triggers");


  tG_ = fs->make<TTree>("gen", "per gen MCP");
  tG_->Branch("run", &g_run_); tG_->Branch("lumi", &g_lumi_); tG_->Branch("event", &g_event_);
  tG_->Branch("gen_pt", &g_pt_); tG_->Branch("gen_eta", &g_eta_); tG_->Branch("gen_phi", &g_phi_);
  tG_->Branch("gen_charge", &g_charge_); tG_->Branch("gen_pdgId", &g_pdgId_);
  tG_->Branch("matched", &g_matched_); tG_->Branch("reco_pt", &g_recoPt_); tG_->Branch("dR", &g_dR_);
  tG_->Branch("chargeFromCurvature", &g_chargeFromCurv_);
  
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

  edm::Handle<pat::METCollection> pfMETCollection;
  iEvent.getByToken(metToken_, pfMETCollection);
  if (!pfMETCollection.isValid()) {
    edm::LogError("DQMClientExample") << "invalid collection: MET"
                                      << "\n";
    return;
  }
  edm::Handle<pat::METCollection> puppiMETCollection;
  iEvent.getByToken(puppiMetToken_, puppiMETCollection);
  if (!puppiMETCollection.isValid()) {
    edm::LogError("DQMClientExample") << "invalid collection: puppiMET"
                                      << "\n";
    return;
  }
  b_pfMET_ = pfMETCollection->front().pt();
  b_puppiMET_ = puppiMETCollection->front().pt();

  b_trigNames_.clear();
  b_trigPass_.clear();
  b_passTrigger_OR.clear(); 

  const auto triggerH = iEvent.getHandle(triggerResultsToken_);
  if (triggerH.isValid()) {
    const auto& triggerNames = iEvent.triggerNames(*triggerH);
    for (unsigned int i = 0; i < triggerH->size(); ++i) {
      TString name(triggerNames.triggerName(i)); 
      b_trigNames_.push_back(triggerNames.triggerName(i));
      b_trigPass_.push_back(triggerH->accept(i) ? 1 : 0);
    }
  }

    int pass_OR = 0;
    
        std::vector<std::string> OR_trigger_list = {
    "HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng",
    "HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1",
    "HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3",
    "HLT_DoublePNetTauhPFJet30_Tight_L2NN_eta2p3",
    "HLT_MET105_IsoTrk50",
    "HLT_MET120_IsoTrk50",
    "HLT_PFMET105_IsoTrk50",
    "HLT_PFMET120_PFMHT120_IDTight",
    "HLT_PFMET120_PFMHT120_IDTight_PFHT60",
    "HLT_PFMET130_PFMHT130_IDTight",
    "HLT_PFMET140_PFMHT140_IDTight",
    "HLT_PFMET200_BeamHaloCleaned",
    "HLT_PFMET250_NotCleaned",
    "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF",
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF",
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
    "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight",
    "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF",
    "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight",
    "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF",
    "HLT_PFMETTypeOne140_PFMHT140_IDTight",
    "HLT_PFMETTypeOne200_BeamHaloCleaned"
  };


    for (unsigned int i = 0; i< b_trigNames_.size(); ++i){
      for (const auto& trigger : OR_trigger_list){
        if (b_trigNames_[i].find(trigger) != std::string::npos && b_trigPass_[i] == 1){
          pass_OR = 1;
        }
      }
    }

  b_passTrigger_OR.push_back(pass_OR); 
  


  // for gen-tree: track best reco match per MCP
  std::vector<double> mcpBestDR(mcps.size(), 1e9);
  std::vector<double> mcpBestRecoPt(mcps.size(), -1.);

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
      b_nPixClusters_ = r.nPixClusters; b_nPixNoFillProb_ = r.nPixNoFillProb;
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
        if (bestDR < mcpBestDR[bestIdx]) { mcpBestDR[bestIdx] = bestDR; mcpBestRecoPt[bestIdx] = b_pt_; }

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
  desc.add<edm::InputTag>("TriggerResults", edm::InputTag("TriggerResults", "", "HLT"));
  descriptions.add("MCPAnalyzer", desc);
}

DEFINE_FWK_MODULE(MCPAnalyzer);
