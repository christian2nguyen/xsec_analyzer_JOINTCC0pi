// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/CCXp0piBinScheme.hh"

CCXp0piBinScheme::CCXp0piBinScheme() : BinSchemeBase( "CCXp0piBinScheme" ) {}

void CCXp0piBinScheme::DefineBlocks() {

  /////// Set some standard variables before managing the blocks

  // TTree name for the post-processed ntuples
  ntuple_ttree_name_ = "stv_tree";

  // Run numbers to use when plotting migration matrices
  runs_to_use_ = { 1 };

  // Prefix for the output bin and slice configuration text files
 
  out_config_prefix_ = "ccxp0pi_";
//#define TKI_1D
#ifdef TKI_1D
  out_config_prefix_ = "ccxp0pi_TKI_1D_";
#endif
//#define TKI_2D
#ifdef TKI_2D
  out_config_prefix_ = "ccxp0pi_TKI_2D_";
#endif
//#define MUON_PROTON
#ifdef MUON_PROTON
  out_config_prefix_ = "ccxp0pi_MUON_PROTON_";
#endif

//#define MUON_PROTON_0P
//#define MUON_PROTON_MULTI_P_MOM_COSTHETA
#ifdef MUON_PROTON_MULTI_P_MOM_COSTHETA
  out_config_prefix_ = "ccxp0pi_MUON_PROTON_MULTI_P_";
#endif
#define MUON_PROTON_MULTI_P_MOM_COSTHETA_0PNP
#ifdef MUON_PROTON_MULTI_P_MOM_COSTHETA_0PNP
  out_config_prefix_ = "ccxp0pi_MUON_PROTON_MULTI_0PNP_";
#endif

  // Selection to use with this binning scheme
  selection_name_ = "CC1muXp0pi";

  // TDirectory file name to use when producing the univmake output histograms
  out_tdir_name_ = "muon_2d_bin";

  /////// Define the blocks of bins in both true and reco space


  // block of proton multiplicity
  std::vector< double > proton_multi_edges = {0, 1.0, 2.0, 3.0, 10.0  };

  Block1D* proton_multi_t = new Block1D("CC1muXp0pi_sig_mc_num_proton_in_momentum_range", "proton multiplicity", "proton multiplicity", proton_multi_edges, "CC1muXp0pi_MC_Signal", kSignalTrueBin);
  Block1D* proton_multi_r = new Block1D("CC1muXp0pi_sel_num_proton_candidates", "proton multiplicity", "proton multiplicity", proton_multi_edges, "CC1muXp0pi_Selected", kOrdinaryRecoBin);
  vect_block.emplace_back( proton_multi_t, proton_multi_r);

#ifdef MUON_PROTON_MULTI_P_MOM_COSTHETA



  std::map<double, std::vector< double > > cos_theta_mu_2D_edges = { 
    {0., { -1, 0.275, 0.425, 0.575, 0.725, 0.85, 0.9, 0.95, 1}},
    {1., { -1., -0.85, -0.775, -0.7, -0.625, -0.55, -0.475, -0.4, -0.325, -0.25, -0.175, -0.1, -0.025, 0.05, 0.125, 0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.85, 0.875, 0.9, 0.925, 0.950, 0.975, 1. }},
    {2., { -1, 0.275, 0.425, 0.575, 0.725, 0.85, 0.9, 0.95, 1}},
    {3., { -1, 0.275, 0.425, 0.575, 0.725, 0.85, 0.9, 0.95, 1}},
    {10., {}}
  };

// First block: cos_theta_mu and mom_mu in 1D 
//  std::vector< double > cos_theta_mu_1D_edges = { -1., -0.85, -0.775, -0.7,
//    -0.625, -0.55, -0.475, -0.4, -0.325, -0.25, -0.175, -0.1, -0.025, 0.05,
//    0.125, 0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.85,
//    0.875, 0.9, 0.925, 0.950, 0.975, 1. };
  Block2D* mu_costh_0p_t = new Block2D( "CC1muXp0pi_sig_mc_num_proton_in_momentum_range; ; mc_p4_mu->CosTheta(); ",
    " proton multiplicity; ; 0p muon cos#theta;", "proton multiplicity ;; \\cos\\theta_{\\mu};", cos_theta_mu_2D_edges,
    "CC1muXp0pi_MC_Signal", kSignalTrueBin );
  Block2D* mu_costh_0p_r = new Block2D( "CC1muXp0pi_sel_num_proton_candidates;; reco_p4mu->CosTheta(); ",
    "proton multiplicity; ; 0p muon cos#theta;", "proton multiplicity; ; \\cos\\theta_{\\mu};", cos_theta_mu_2D_edges,
    "(CC1muXp0pi_Selected)", kOrdinaryRecoBin );
  vect_block.emplace_back( mu_costh_0p_t, mu_costh_0p_r );


  std::map< double , std::vector< double > > pmu_2D_edges = {
    {0., { 0.1, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2}},
    {1., { 0.1, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2 }},
    {2., { 0.1, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2}},
    {3., { 0.1, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2}},
    {10., {}}
  };

  // Second block: p_mu in 1D
//  std::vector< double > pmu_1D_edges = { 0.1, 0.175, 0.2, 0.225, 0.25,
//    0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.55, 0.6,
//    0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2 };

  Block2D* mu_rho_0p_t = new Block2D( "CC1muXp0pi_sig_mc_num_proton_in_momentum_range; ;mc_p4_mu->Rho(); GeV/c",
    "proton multiplicity; ;0p p_{#mu}; (GeV)", "proton multiplicity; ;p_{\\mu}; (GeV)", pmu_2D_edges,
    "CC1muXp0pi_MC_Signal", kSignalTrueBin );

  Block2D* mu_rho_0p_r = new Block2D( "CC1muXp0pi_sel_num_proton_candidates;;reco_p4mu->Rho(); GeV/c",
    "proton multiplicity; ;0p p_{#mu}; (GeV)", "proton multiplicity; ;p_{\\mu}; (GeV)", pmu_2D_edges,
    "(CC1muXp0pi_Selected)", kOrdinaryRecoBin );

  vect_block.emplace_back( mu_rho_0p_t, mu_rho_0p_r );


#endif


#ifdef MUON_PROTON_MULTI_P_MOM_COSTHETA_0PNP



  std::map<double, std::vector< double > > cos_theta_mu_2D_edges = { 
    {0., { -1, 0.275, 0.425, 0.575, 0.725, 0.85, 0.9, 0.95, 1}},
    {1., { -1., -0.85, -0.775, -0.7, -0.625, -0.55, -0.475, -0.4, -0.325, -0.25, -0.175, -0.1, -0.025, 0.05, 0.125, 0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.85, 0.875, 0.9, 0.925, 0.950, 0.975, 1. }},
    {10., {}}
  };

// First block: cos_theta_mu and mom_mu in 1D 
//  std::vector< double > cos_theta_mu_1D_edges = { -1., -0.85, -0.775, -0.7,
//    -0.625, -0.55, -0.475, -0.4, -0.325, -0.25, -0.175, -0.1, -0.025, 0.05,
//    0.125, 0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.85,
//    0.875, 0.9, 0.925, 0.950, 0.975, 1. };
  Block2D* mu_costh_0p_t = new Block2D( "CC1muXp0pi_sig_mc_num_proton_in_momentum_range; ; mc_p4_mu->CosTheta(); ",
    " proton multiplicity; ; 0p muon cos#theta;", "proton multiplicity ;; \\cos\\theta_{\\mu};", cos_theta_mu_2D_edges,
    "CC1muXp0pi_MC_Signal", kSignalTrueBin );
  Block2D* mu_costh_0p_r = new Block2D( "CC1muXp0pi_sel_num_proton_candidates;; reco_p4mu->CosTheta(); ",
    "proton multiplicity; ; 0p muon cos#theta;", "proton multiplicity; ; \\cos\\theta_{\\mu};", cos_theta_mu_2D_edges,
    "(CC1muXp0pi_Selected)", kOrdinaryRecoBin );
  vect_block.emplace_back( mu_costh_0p_t, mu_costh_0p_r );


  std::map< double , std::vector< double > > pmu_2D_edges = {
    {0., { 0.1, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2}},
    {1., { 0.1, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2 }},
    {10., {}}
  };

  // Second block: p_mu in 1D
//  std::vector< double > pmu_1D_edges = { 0.1, 0.175, 0.2, 0.225, 0.25,
//    0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.55, 0.6,
//    0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2 };

  Block2D* mu_rho_0p_t = new Block2D( "CC1muXp0pi_sig_mc_num_proton_in_momentum_range; ;mc_p4_mu->Rho(); GeV/c",
    "proton multiplicity; ;0p p_{#mu}; (GeV)", "proton multiplicity; ;p_{\\mu}; (GeV)", pmu_2D_edges,
    "CC1muXp0pi_MC_Signal", kSignalTrueBin );

  Block2D* mu_rho_0p_r = new Block2D( "CC1muXp0pi_sel_num_proton_candidates;;reco_p4mu->Rho(); GeV/c",
    "proton multiplicity; ;0p p_{#mu}; (GeV)", "proton multiplicity; ;p_{\\mu}; (GeV)", pmu_2D_edges,
    "(CC1muXp0pi_Selected)", kOrdinaryRecoBin );

  vect_block.emplace_back( mu_rho_0p_t, mu_rho_0p_r );

#endif



#ifdef MUON_PROTON

// First block: cos_theta_mu and mom_mu in 1D 
  std::vector< double > cos_theta_mu_1D_edges = { -1., -0.85, -0.775, -0.7,
    -0.625, -0.55, -0.475, -0.4, -0.325, -0.25, -0.175, -0.1, -0.025, 0.05,
    0.125, 0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.85,
    0.875, 0.9, 0.925, 0.950, 0.975, 1. };
  Block1D* mu_costh_t = new Block1D( "mc_p4_mu->CosTheta()",
    "Xp muon cos#theta", "\\cos\\theta_{\\mu}", cos_theta_mu_1D_edges,
    "CC1muXp0pi_MC_Signal", kSignalTrueBin );
  Block1D* mu_costh_r = new Block1D( "reco_p4mu->CosTheta()",
    "Xp muon cos#theta", "\\cos\\theta_{\\mu}", cos_theta_mu_1D_edges,
    "(CC1muXp0pi_Selected)", kOrdinaryRecoBin );
  vect_block.emplace_back( mu_costh_t, mu_costh_r );

  // Second block: p_mu in 1D
  std::vector< double > pmu_1D_edges = { 0.1, 0.175, 0.2, 0.225, 0.25,
    0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.55, 0.6,
    0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2 };

  Block1D* mu_rho_t = new Block1D( "mc_p4_mu->Rho()",
    "Xp p_{#mu}; (GeV)", "p_{\\mu}; (GeV)", pmu_1D_edges,
    "CC1muXp0pi_MC_Signal", kSignalTrueBin );

  Block1D* mu_rho_r = new Block1D( "reco_p4mu->Rho()",
    "Xp p_{#mu}; (GeV)", "p_{\\mu}; (GeV)", pmu_1D_edges,
    "(CC1muXp0pi_Selected)", kOrdinaryRecoBin );

  vect_block.emplace_back( mu_rho_t, mu_rho_r );


// First block: cos_theta_mu and mom_mu in 1D 
//  std::vector< double > cos_theta_mu_1D_edges = { -1., -0.85, -0.775, -0.7,
//    -0.625, -0.55, -0.475, -0.4, -0.325, -0.25, -0.175, -0.1, -0.025, 0.05,
//    0.125, 0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.85,
//    0.875, 0.9, 0.925, 0.950, 0.975, 1. };
  Block1D* mu_costh_Np_t = new Block1D( "mc_p4_mu->CosTheta()",
    "Np muon cos#theta", "\\cos\\theta_{\\mu}", cos_theta_mu_1D_edges,
    "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range > 0", kSignalTrueBin );
  Block1D* mu_costh_Np_r = new Block1D( "reco_p4mu->CosTheta()",
    "Np muon cos#theta", "\\cos\\theta_{\\mu}", cos_theta_mu_1D_edges,
    "(CC1muXp0pi_Selected && CC1muXp0pi_sel_num_proton_candidates > 0)", kOrdinaryRecoBin );
  vect_block.emplace_back( mu_costh_Np_t, mu_costh_Np_r );

  // Second block: p_mu in 1D
//  std::vector< double > pmu_1D_edges = { 0.1, 0.175, 0.2, 0.225, 0.25,
//    0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.55, 0.6,
//    0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2 };
  Block1D* mu_rho_Np_t = new Block1D( "mc_p4_mu->Rho()",
    "Np p_{#mu}; (GeV)", "p_{\\mu}; (GeV)", pmu_1D_edges,
    "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range > 0", kSignalTrueBin );
  Block1D* mu_rho_Np_r = new Block1D( "reco_p4mu->Rho()",
    "Np p_{#mu}; (GeV)", "p_{\\mu}; (GeV)", pmu_1D_edges,
    "(CC1muXp0pi_Selected && CC1muXp0pi_sel_num_proton_candidates > 0)", kOrdinaryRecoBin );
  vect_block.emplace_back( mu_rho_Np_t, mu_rho_Np_r );



  std::map< double, std::vector<double> > MUON_2D_BIN_EDGES = {
//
    // No need for an underflow bin: due to the signal definition, all muons
    // with reco momentum below 0.1 GeV/c will be lost
    { 0.1, { -1, -0.55, -0.25, 0., 0.25, 0.45, 0.7, 1.00 }, },
    { 0.24, { -1, -0.55, -0.25, 0., 0.25, 0.45, 0.7, 1.00 } },
    { 0.3,  { -1, -0.4, -0.1, 0.1, 0.35, 0.5, 0.7, 0.85, 1. } },
    { 0.38, { -1, 0, 0.5, 0.65, 0.8, 0.92, 1.00 } },
    { 0.48, { -1, 0.2, 0.5, 0.65, 0.8, 0.875, 0.950, 1.00 } },
    { 0.7, { -1, 0.65, 0.8, 0.875, 0.950, 1.00 } },
    { 0.85, { -1, 0.85, 0.9, 0.950, 1.00 } },

    // Upper edge of the last bin. Due to the signal definition, no overflow
    // bin is needed for muons above 1.2 GeV/c
    { 1.2, {} }

  };

  Block2D *b2dmuon_mom_costheta_true = new Block2D("mc_p4_mu->Rho(); GeV/c; mc_p4_mu->CosTheta(); ", 
      "Xp 2D muon mom; GeV/c;Xp 2D muon cos#theta; ",
      "P_{\\mu}; GeV/c; \\cos\\theta; ",
      MUON_2D_BIN_EDGES, "CC1muXp0pi_MC_Signal", kSignalTrueBin);
  Block2D *b2dmuon_mom_costheta_reco = new Block2D("reco_p4mu->Rho(); GeV/c; reco_p4mu->CosTheta(); ",
      "Xp 2D muon mom; GeV/c;Xp 2D muon cos#theta; ", 
      "P_{\\mu}; GeV/c; \\cos\\theta; ",
      MUON_2D_BIN_EDGES, "CC1muXp0pi_Selected", kOrdinaryRecoBin);
  vect_block.emplace_back(b2dmuon_mom_costheta_true, b2dmuon_mom_costheta_reco);





  // cos_theta_p in 1D
  std::vector< double > cos_theta_p_1D_edges = { -1., -0.9, -0.75, -0.6,
  -0.45, -0.3, -0.15, 0.0, 0.15, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9,
  0.925, 0.95, 0.975, 1.0 };

  Block1D* cos_theta_p_t = new Block1D( "mc_p4p->CosTheta()",
    "Np cos#theta_{P};", "\\cos\\theta_{p};", cos_theta_p_1D_edges,
    "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range > 0", kSignalTrueBin );
  Block1D* cos_theta_p_r = new Block1D( "reco_p4p->CosTheta()",
    "Np cos#theta_{P};", "\\cos\\theta_{p};", cos_theta_p_1D_edges,
    "(CC1muXp0pi_Selected && CC1muXp0pi_sel_num_proton_candidates > 0)", kOrdinaryRecoBin );
  vect_block.emplace_back( cos_theta_p_t, cos_theta_p_r );

// p_p in 1D
  std::vector< double > p_rho_1D_edges = { 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
  0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0 };

  Block1D* p_rho_t = new Block1D( "mc_p4p->Rho()",
    "Np p_{P};", "p_{P};", p_rho_1D_edges,
    "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range > 0", kSignalTrueBin );
  Block1D* p_rho_r = new Block1D( "reco_p4p->Rho()",
    "Np p_{P};", "p_{P};", p_rho_1D_edges,
    "(CC1muXp0pi_Selected && CC1muXp0pi_sel_num_proton_candidates > 0)", kOrdinaryRecoBin );
  vect_block.emplace_back( p_rho_t, p_rho_r );


  // Using floating-point numbers as std::map keys is admittedly evil, but it's
  // safe in this case: all we'll do with this map is iterate over the
  // elements. Keys are proton momentum bin edges, values are proton cosine bin
  // edges.
  std::map< double, std::vector<double> > PROTON_2D_BIN_EDGES = {

    // No need for an underflow bin: due to the signal definition, all leading
    // protons with reco momentum below 0.25 GeV/c will be lost
    { 0.250, { -1, 0., 1.0 } },
    { 0.325, { -1, -0.5, 0, 0.5, 0.8, 1.0 } },
    { 0.4,   { -1, -0.6, -0.2, 0.2, 0.5, 0.65, 0.85, 1.0 } },
    { 0.5, { -1, -0.2, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0 } },
    { 0.6,   { -1, 0.1, 0.37, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 } },
    { 0.7,   { -1, 0.45, 0.65, 0.75, 0.82, 0.9, 1.0 } },

    // Upper edge of the last bin. We don't need an overflow bin because the
    // signal definition excludes any leading protons with momenta above
    // 1 GeV/c
    { 1., {} }

  };

  Block2D *proton_mom_costheta_true = new Block2D("mc_p4p->Rho(); GeV/c; mc_p4p->CosTheta(); ", 
      "Np 2D proton mom; GeV/c; Np 2D proton cos#theta; ",
      "P_{p}; GeV/c; \\cos\\theta; ",
      PROTON_2D_BIN_EDGES, "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range > 0", kSignalTrueBin);
  Block2D *proton_mom_costheta_reco = new Block2D("reco_p4p->Rho(); GeV/c; reco_p4p->CosTheta(); ",
      "Np 2D proton mom; GeV/c; Np 2D proton cos#theta; ", 
      "P_{p}; GeV/c; \\cos\\theta; ",
      PROTON_2D_BIN_EDGES, "CC1muXp0pi_Selected  && CC1muXp0pi_sel_num_proton_candidates > 0", kOrdinaryRecoBin);
  vect_block.emplace_back(proton_mom_costheta_true, proton_mom_costheta_reco);

#endif

  std::string branchexpr, title, textitle, selection;

#ifdef TKI_1D



  std::vector< double > pn_1D_edges = { 0., 0.06, 0.12, 0.18, 0.24, 0.3, 0.45, 0.6, 0.75, 0.9, 1.2 };
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_mc_lead_proton_pn ; GeV/c";
  // the title "title; unit" is used in plot in root style
  title = "leading proton Np p_{n}; GeV/c";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton  \\p_{n}; GeV/c";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range > 0";
  Block1D *b1dt_Np_delta_pn = new Block1D(branchexpr, title, textitle, pn_1D_edges, selection, kSignalTrueBin);
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_lead_proton_pn ; GeV/c";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_sel_num_proton_candidates > 0";
  // only the name of branch and the selection is different from true.
  Block1D *b1dr_Np_delta_pn = new Block1D(branchexpr, title, textitle, pn_1D_edges, selection, kOrdinaryRecoBin);
  vect_block.emplace_back(b1dt_Np_delta_pn, b1dr_Np_delta_pn);



//  std::vector< double > pn_1D_edges = { 0., 0.06, 0.12, 0.18, 0.24, 0.3,
//      0.45, 0.6, 0.75, 0.9, 1.2 };
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_mc_lead_proton_pn ; GeV/c";
  // the title "title; unit" is used in plot in root style
  title = "leading proton 1p  p_{n}; GeV/c";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton  \\p_{n}; GeV/c";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range == 1";

  Block1D *b1dt_1p_delta_pn = new Block1D(branchexpr, title, textitle, pn_1D_edges, selection, kSignalTrueBin);

  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_lead_proton_pn ; GeV/c";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_sel_num_proton_candidates == 1";

  // only the name of branch and the selection is different from true.
  Block1D *b1dr_1p_delta_pn = new Block1D(branchexpr, title, textitle, pn_1D_edges, selection, kOrdinaryRecoBin);

  vect_block.emplace_back(b1dt_1p_delta_pn, b1dr_1p_delta_pn);


  std::vector< double > pT_1D_edges = { 0., 0.06, 0.12, 0.18, 0.24, 0.3, 0.45, 0.6, 0.75, 0.9, 1.2};
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_mc_lead_proton_pT ; GeV/c";
  // the title "title; unit" is used in plot in root style
  title = "leading proton Np p_{T}; GeV/c";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton  \\p_{T}; GeV/c";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range > 0";
  Block1D *b1dt_Np_delta_pT = new Block1D(branchexpr, title, textitle, pT_1D_edges, selection, kSignalTrueBin);
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_lead_proton_pT ; GeV/c";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_sel_num_proton_candidates > 0";
  // only the name of branch and the selection is different from true.
  Block1D *b1dr_Np_delta_pT = new Block1D(branchexpr, title, textitle, pT_1D_edges, selection, kOrdinaryRecoBin);
  vect_block.emplace_back(b1dt_Np_delta_pT, b1dr_Np_delta_pT);



//  std::vector< double > qT_1D_edges = { 0., 0.06, 0.12, 0.18, 0.24, 0.3,
//      0.45, 0.6, 0.75, 0.9, 1.2 };
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_mc_lead_proton_pT ; GeV/c";
  // the title "title; unit" is used in plot in root style
  title = "leading proton 1p p_{T}; GeV/c";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton  \\p_{T}; GeV/c";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range == 1";

  Block1D *b1dt_1p_delta_pT = new Block1D(branchexpr, title, textitle, pT_1D_edges, selection, kSignalTrueBin);

  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_lead_proton_pT ; GeV/c";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_sel_num_proton_candidates == 1";

  // only the name of branch and the selection is different from true.
  Block1D *b1dr_1p_delta_pT = new Block1D(branchexpr, title, textitle, pT_1D_edges, selection, kOrdinaryRecoBin);

  vect_block.emplace_back(b1dt_1p_delta_pT, b1dr_1p_delta_pT);


  std::vector< double > delta_alphaT_1D_edges = {0, 0.19635, 0.392699, 0.589049, 0.785398, 0.981748, 1.1781, 1.37445, 1.5708, 1.76715, 1.9635, 2.15984, 2.35619, 2.55254, 2.74889, 2.94524, TMath::Pi()};
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_mc_lead_proton_delta_alphaT ; GeV/c";
  // the title "title; unit" is used in plot in root style
  title = "leading proton Np #delta #alpah_{T}; GeV/c";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton  \\p_{T}; GeV/c";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range > 0";
  Block1D *b1dt_Np_delta_delta_alphaT = new Block1D(branchexpr, title, textitle, delta_alphaT_1D_edges, selection, kSignalTrueBin);
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_lead_proton_delta_alphaT ; GeV/c";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_sel_num_proton_candidates > 0";
  // only the name of branch and the selection is different from true.
  Block1D *b1dr_Np_delta_delta_alphaT = new Block1D(branchexpr, title, textitle, delta_alphaT_1D_edges, selection, kOrdinaryRecoBin);
  vect_block.emplace_back(b1dt_Np_delta_delta_alphaT, b1dr_Np_delta_delta_alphaT);



//  std::vector< double > qT_1D_edges = { 0., 0.06, 0.12, 0.18, 0.24, 0.3,
//      0.45, 0.6, 0.75, 0.9, 1.2 };
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_mc_lead_proton_delta_alphaT ; GeV/c";
  // the title "title; unit" is used in plot in root style
  title = "leading proton  1p #delta #alpah_{T}; GeV/c";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton  \\p_{n}; GeV/c";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range == 1";

  Block1D *b1dt_1p_delta_delta_alphaT = new Block1D(branchexpr, title, textitle, delta_alphaT_1D_edges, selection, kSignalTrueBin);

  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_lead_proton_delta_alphaT ; GeV/c";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_sel_num_proton_candidates == 1";

  // only the name of branch and the selection is different from true.
  Block1D *b1dr_1p_delta_delta_alphaT = new Block1D(branchexpr, title, textitle, delta_alphaT_1D_edges, selection, kOrdinaryRecoBin);

  vect_block.emplace_back(b1dt_1p_delta_delta_alphaT, b1dr_1p_delta_delta_alphaT);


  std::vector< double > delta_phiT_1D_edges = { 0, 0.0616616, 0.127377, 0.197716, 0.27338, 0.35524, 0.444405, 0.542304, 0.650839, 0.772604, 0.911278, 1.07233, 1.26439, 1.50232, 1.81515, 2.27288, TMath::Pi()};
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_mc_lead_proton_delta_phiT ; GeV/c";
  // the title "title; unit" is used in plot in root style
  title = "leading proton  Np #delta #phi_{T}; GeV/c";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton  \\p_{T}; GeV/c";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range > 0";
  Block1D *b1dt_Np_delta_delta_phiT = new Block1D(branchexpr, title, textitle, delta_phiT_1D_edges, selection, kSignalTrueBin);
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_lead_proton_delta_phiT ; GeV/c";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_sel_num_proton_candidates > 0";
  // only the name of branch and the selection is different from true.
  Block1D *b1dr_Np_delta_delta_phiT = new Block1D(branchexpr, title, textitle, delta_phiT_1D_edges, selection, kOrdinaryRecoBin);
  vect_block.emplace_back(b1dt_Np_delta_delta_phiT, b1dr_Np_delta_delta_phiT);



//  std::vector< double > qT_1D_edges = { 0., 0.06, 0.12, 0.18, 0.24, 0.3,
//      0.45, 0.6, 0.75, 0.9, 1.2 };
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_mc_lead_proton_delta_phiT ; GeV/c";
  // the title "title; unit" is used in plot in root style
  title = "leading proton  1p #delta #phi_{T}; GeV/c";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton  \\p_{n}; GeV/c";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range == 1";

  Block1D *b1dt_1p_delta_delta_phiT = new Block1D(branchexpr, title, textitle, delta_phiT_1D_edges, selection, kSignalTrueBin);

  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_lead_proton_delta_phiT ; GeV/c";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_sel_num_proton_candidates == 1";

  // only the name of branch and the selection is different from true.
  Block1D *b1dr_1p_delta_delta_phiT = new Block1D(branchexpr, title, textitle, delta_phiT_1D_edges, selection, kOrdinaryRecoBin);

  vect_block.emplace_back(b1dt_1p_delta_delta_phiT, b1dr_1p_delta_delta_phiT);

#endif

#ifdef TKI_2D


  // kinetic energy of proton
  //

  std::vector< double > ke_proton_edges = {0.03, 0.0539999, 0.0812775, 0.112872, 0.150411, 0.196664, 0.256935, 0.343614, 0.50 };

  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "(mc_p4p->E() - mc_p4p->M()) ; GeV";
  // the title "title; unit" is used in plot in root style
  title = "leading proton 1p kinetic energy ; GeV/c";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton  KE; GeV/c";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range == 1";

  Block1D *b1dt_1p_proton_KE = new Block1D(branchexpr, title, textitle, ke_proton_edges, selection, kSignalTrueBin);
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "(reco_p4p->E() - reco_p4p->M()) ; GeV";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_sel_num_proton_candidates == 1";
  // only the name of branch and the selection is different from true.
  Block1D *b1dr_1p_proton_KE = new Block1D(branchexpr, title, textitle, ke_proton_edges, selection, kOrdinaryRecoBin);
  vect_block.emplace_back(b1dt_1p_proton_KE, b1dr_1p_proton_KE);

  std::map< double, std::vector< double > > proton_KE_pn_2D_edges ={
    {0.03 , { 0., 0.06, 0.12, 0.18, 0.24, 0.3, 0.45, 0.6, 0.75, 0.9, 1.2 }},
    {0.0539999 , { 0., 0.06, 0.12, 0.18, 0.24, 0.3, 0.45, 0.6, 0.75, 0.9, 1.2 }},
    {0.0812775 , { 0., 0.06, 0.12, 0.18, 0.24, 0.3, 0.45, 0.6, 0.75, 0.9, 1.2 }},
    {0.112872 , { 0., 0.06, 0.12, 0.18, 0.24, 0.3, 0.45, 0.6, 0.75, 0.9, 1.2 }},
    {0.150411 , { 0., 0.06, 0.12, 0.18, 0.24, 0.3, 0.45, 0.6, 0.75, 0.9, 1.2 }},
    {0.196664 , { 0., 0.06, 0.12, 0.18, 0.24, 0.3, 0.45, 0.6, 0.75, 0.9, 1.2 }},
    {0.256935 , { 0., 0.06, 0.12, 0.18, 0.24, 0.3, 0.45, 0.6, 0.75, 0.9, 1.2 }},
    { 0.343614, { 0., 0.06, 0.12, 0.18, 0.24, 0.3, 0.45, 0.6, 0.75, 0.9, 1.2 }},
    { 0.5, { }}
  };
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "(mc_p4p->E() - mc_p4p->M()) ; GeV; CC1muXp0pi_mc_lead_proton_pn ; GeV/c";
  // the title "title; unit" is used in plot in root style
  title = "leading proton 1p 2D kinetic energy ; GeV; leading proton 1p 2D p_{n}; GeV/c";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton  KE; GeV; leading proton  \\p_{n}; GeV/c";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range == 1 ";
  Block2D *b2dt_proton_KE_pn = new Block2D(branchexpr, title, textitle, proton_KE_pn_2D_edges, selection, kSignalTrueBin);
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "(reco_p4p->E() - reco_p4p->M()) ; GeV; CC1muXp0pi_lead_proton_pn ; GeV/c";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_sel_num_proton_candidates == 1";
  // only the name of branch and the selection is different from true.
  Block2D *b2dr_proton_KE_pn = new Block2D(branchexpr, title, textitle, proton_KE_pn_2D_edges, selection, kOrdinaryRecoBin);
  vect_block.emplace_back(b2dt_proton_KE_pn, b2dr_proton_KE_pn);


#endif



  // a example to use sideband bins; now we only use sideband bins for reco. 
  // you can definie many blocks of sideband with different selections
  // std::vector< double > test_sideband = {1,2,3,4};
  
  // Block1D *b1r_sideband = new Block1D("sideband", "sideband", "sideband",test_sideband , "selection", kSidebandRecoBin);

  // vect_sideband.emplace_back(b1r_sideband);

  // CATEGORY is the branchexpr
  // background_index is vector of background categories.

#ifdef MUON_PROTON_0P
  CATEGORY = "CC1muXp0pi_EventCategory";
  background_index = {0, 5, 6, 7, 8, 
    9, 10, 11, 12,
    13, 14, 15, 16,
    17, 18, 19, 20,
    21, 22, 23, 24, 25, 26};
#endif
  CATEGORY = "CC1muXp0pi_EventCategory";
  background_index = {0, 21, 22, 23, 24, 25, 26};

}
