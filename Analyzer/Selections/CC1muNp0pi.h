#ifndef __CC1muNp0pi_h__
#define __CC1muNp0pi_h__

#include "SelectionBase.h"

class CC1muNp0pi : virtual SelectionBase {
 public:
  CC1muNp0pi();
  
  bool Selection(AnalysisEvent* Event);
  bool DefineSignal(AnalysisEvent* Event);
  void ComputeRecoObservables(AnalysisEvent* Event);
  void ComputeTrueObservables(AnalysisEvent* Event);
  void DefineBranches();
  void DefineConstants();
  
private:
  bool sig_isNuMu_;
  bool sig_inFV_;
  bool sig_leadProtonMomInRange_;
  bool sig_muonInMomRange_;
  bool sig_noFSMesons_;
  bool sig_mc_no_fs_pi0_;
  bool sig_mc_no_charged_pi_above_threshold_;
  
  bool sel_reco_vertex_in_FV_;
  bool sel_pfp_starts_in_PCV_;
  bool sel_has_muon_candidate_;
  bool sel_topo_cut_passed_;
  bool sel_nu_mu_cc_;
  bool sel_muon_contained_;
  bool sel_muon_passed_mom_cuts_;
  bool sel_no_reco_showers_;
  bool sel_has_p_candidate_;
  bool sel_muon_quality_ok_;
  bool sel_protons_contained_;
  bool sel_passed_proton_pid_cut_;
  bool sel_lead_p_passed_mom_cuts_;
  bool sel_cosmic_ip_cut_passed_;
  
  int lead_p_candidate_idx_;
  int muon_candidate_idx_;

  float delta_pT_;
  float delta_phiT_;
  float delta_alphaT_;
  float delta_pL_;
  float pn_;
  float delta_pTx_;
  float delta_pTy_;
  float theta_mu_p_;

  MyPointer<TVector3> p3mu;
  MyPointer<TVector3> p3p;
  MyPointer<std::vector<TVector3>> p3_p_vec_;
  
  float mc_delta_pT_;
  float mc_delta_phiT_;
  float mc_delta_alphaT_;
  float mc_delta_pL_;
  float mc_pn_;
  float mc_delta_pTx_;
  float mc_delta_pTy_;
  float mc_theta_mu_p_;

  MyPointer<TVector3> mc_p3mu;
  MyPointer<TVector3> mc_p3p;
  MyPointer<std::vector<TVector3>> mc_p3_p_vec_;
};

#endif
