//**********************************************
//HEADER FILE FOR SINGLE TOP ANALYSIS
//*********************************************

/*
This header file contains a class called nano9Ana. 
(1) First the TTreeReader is used to declare the variables
(2) These come from three fReaders, a common one and one each for data and MC
(3) Then functions are declared. Focus on the User Added Functions first.
(4) Then variables are declared.
 */



#ifndef SingleTopAna_h
#define SingleTopAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include "TLorentzVector.h"
#include <fstream>
#include <iostream>
#include "TString.h"
#include <bitset>


class SingleTopAna : public TSelector {
public :
  TTreeReader     fReader;       //reads the common branches
  TTreeReader     fReader_MC;    //reads the MC branches
  TTreeReader     fReader_Data;  //reads the Data branches
  TTree          *fChain = 0;    //!pointer to the analyzed TTree or TChain
  
  
  // Readers to access the data (delete the ones you do not need).
  
  //#####################################################
  // The following Common branches are read with fReader
  //#####################################################
  TTreeReaderValue<UInt_t> run = {fReader, "run"};
  TTreeReaderValue<UInt_t> luminosityBlock = {fReader, "luminosityBlock"};
  TTreeReaderValue<ULong64_t> event = {fReader, "event"};

  //HTXS
  /*
  TTreeReaderValue<Float_t> HTXS_Higgs_pt = {fReader, "HTXS_Higgs_pt"};
  TTreeReaderValue<Float_t> HTXS_Higgs_y = {fReader, "HTXS_Higgs_y"};
  TTreeReaderValue<Int_t> HTXS_stage1_1_cat_pTjet25GeV = {fReader, "HTXS_stage1_1_cat_pTjet25GeV"};
  TTreeReaderValue<Int_t> HTXS_stage1_1_cat_pTjet30GeV = {fReader, "HTXS_stage1_1_cat_pTjet30GeV"};
  TTreeReaderValue<Int_t> HTXS_stage1_1_fine_cat_pTjet25GeV = {fReader, "HTXS_stage1_1_fine_cat_pTjet25GeV"};
  TTreeReaderValue<Int_t> HTXS_stage1_1_fine_cat_pTjet30GeV = {fReader, "HTXS_stage1_1_fine_cat_pTjet30GeV"};
  TTreeReaderValue<Int_t> HTXS_stage1_2_cat_pTjet25GeV = {fReader, "HTXS_stage1_2_cat_pTjet25GeV"};
  TTreeReaderValue<Int_t> HTXS_stage1_2_cat_pTjet30GeV = {fReader, "HTXS_stage1_2_cat_pTjet30GeV"};
  TTreeReaderValue<Int_t> HTXS_stage1_2_fine_cat_pTjet25GeV = {fReader, "HTXS_stage1_2_fine_cat_pTjet25GeV"};
  TTreeReaderValue<Int_t> HTXS_stage1_2_fine_cat_pTjet30GeV = {fReader, "HTXS_stage1_2_fine_cat_pTjet30GeV"};
  TTreeReaderValue<Int_t> HTXS_stage_0 = {fReader, "HTXS_stage_0"};
  TTreeReaderValue<Int_t> HTXS_stage_1_pTjet25 = {fReader, "HTXS_stage_1_pTjet25"};
  TTreeReaderValue<Int_t> HTXS_stage_1_pTjet30 = {fReader, "HTXS_stage_1_pTjet30"};
  TTreeReaderValue<UChar_t> HTXS_njets25 = {fReader, "HTXS_njets25"};
  TTreeReaderValue<UChar_t> HTXS_njets30 = {fReader, "HTXS_njets30"};
  */
  //Boosted Tau
  /*
  TTreeReaderValue<UInt_t> nboostedTau = {fReader, "nboostedTau"};
  TTreeReaderArray<Float_t> boostedTau_chargedIso = {fReader, "boostedTau_chargedIso"};
  TTreeReaderArray<Float_t> boostedTau_eta = {fReader, "boostedTau_eta"};
  TTreeReaderArray<Float_t> boostedTau_leadTkDeltaEta = {fReader, "boostedTau_leadTkDeltaEta"};
  TTreeReaderArray<Float_t> boostedTau_leadTkDeltaPhi = {fReader, "boostedTau_leadTkDeltaPhi"};
  TTreeReaderArray<Float_t> boostedTau_leadTkPtOverTauPt = {fReader, "boostedTau_leadTkPtOverTauPt"};
  TTreeReaderArray<Float_t> boostedTau_mass = {fReader, "boostedTau_mass"};
  TTreeReaderArray<Float_t> boostedTau_neutralIso = {fReader, "boostedTau_neutralIso"};
  TTreeReaderArray<Float_t> boostedTau_phi = {fReader, "boostedTau_phi"};
  TTreeReaderArray<Float_t> boostedTau_photonsOutsideSignalCone = {fReader, "boostedTau_photonsOutsideSignalCone"};
  TTreeReaderArray<Float_t> boostedTau_pt = {fReader, "boostedTau_pt"};
  TTreeReaderArray<Float_t> boostedTau_puCorr = {fReader, "boostedTau_puCorr"};
  TTreeReaderArray<Float_t> boostedTau_rawAntiEle2018 = {fReader, "boostedTau_rawAntiEle2018"};
  TTreeReaderArray<Float_t> boostedTau_rawIso = {fReader, "boostedTau_rawIso"};
  TTreeReaderArray<Float_t> boostedTau_rawIsodR03 = {fReader, "boostedTau_rawIsodR03"};
  TTreeReaderArray<Float_t> boostedTau_rawMVAnewDM2017v2 = {fReader, "boostedTau_rawMVAnewDM2017v2"};
  TTreeReaderArray<Float_t> boostedTau_rawMVAoldDM2017v2 = {fReader, "boostedTau_rawMVAoldDM2017v2"};
  TTreeReaderArray<Float_t> boostedTau_rawMVAoldDMdR032017v2 = {fReader, "boostedTau_rawMVAoldDMdR032017v2"};
  TTreeReaderArray<Int_t> boostedTau_charge = {fReader, "boostedTau_charge"};
  TTreeReaderArray<Int_t> boostedTau_decayMode = {fReader, "boostedTau_decayMode"};
  TTreeReaderArray<Int_t> boostedTau_jetIdx = {fReader, "boostedTau_jetIdx"};
  TTreeReaderArray<Int_t> boostedTau_rawAntiEleCat2018 = {fReader, "boostedTau_rawAntiEleCat2018"};
  TTreeReaderArray<UChar_t> boostedTau_idAntiEle2018 = {fReader, "boostedTau_idAntiEle2018"};
  TTreeReaderArray<UChar_t> boostedTau_idAntiMu = {fReader, "boostedTau_idAntiMu"};
  TTreeReaderArray<UChar_t> boostedTau_idMVAnewDM2017v2 = {fReader, "boostedTau_idMVAnewDM2017v2"};
  TTreeReaderArray<UChar_t> boostedTau_idMVAoldDM2017v2 = {fReader, "boostedTau_idMVAoldDM2017v2"};
  TTreeReaderArray<UChar_t> boostedTau_idMVAoldDMdR032017v2 = {fReader, "boostedTau_idMVAoldDMdR032017v2"};
  */
  //btagWeight
  /*
  TTreeReaderValue<Float_t> btagWeight_CSVV2 = {fReader, "btagWeight_CSVV2"};
  TTreeReaderValue<Float_t> btagWeight_DeepCSVB = {fReader, "btagWeight_DeepCSVB"};
  */
  // CaloMET
  /*
  TTreeReaderValue<Float_t> CaloMET_phi = {fReader, "CaloMET_phi"};
  TTreeReaderValue<Float_t> CaloMET_pt = {fReader, "CaloMET_pt"};
  TTreeReaderValue<Float_t> CaloMET_sumEt = {fReader, "CaloMET_sumEt"};
  TTreeReaderValue<Float_t> ChsMET_phi = {fReader, "ChsMET_phi"};
  TTreeReaderValue<Float_t> ChsMET_pt = {fReader, "ChsMET_pt"};
  TTreeReaderValue<Float_t> ChsMET_sumEt = {fReader, "ChsMET_sumEt"};
  */
  //CorrT1METJet
  /*
  TTreeReaderValue<UInt_t> nCorrT1METJet = {fReader, "nCorrT1METJet"};
  TTreeReaderArray<Float_t> CorrT1METJet_area = {fReader, "CorrT1METJet_area"};
  TTreeReaderArray<Float_t> CorrT1METJet_eta = {fReader, "CorrT1METJet_eta"};
  TTreeReaderArray<Float_t> CorrT1METJet_muonSubtrFactor = {fReader, "CorrT1METJet_muonSubtrFactor"};
  TTreeReaderArray<Float_t> CorrT1METJet_phi = {fReader, "CorrT1METJet_phi"};
  TTreeReaderArray<Float_t> CorrT1METJet_rawPt = {fReader, "CorrT1METJet_rawPt"};
  */
  // DeepMET
  /*
  TTreeReaderValue<Float_t> DeepMETResolutionTune_phi = {fReader, "DeepMETResolutionTune_phi"};
  TTreeReaderValue<Float_t> DeepMETResolutionTune_pt = {fReader, "DeepMETResolutionTune_pt"};
  TTreeReaderValue<Float_t> DeepMETResponseTune_phi = {fReader, "DeepMETResponseTune_phi"};
  TTreeReaderValue<Float_t> DeepMETResponseTune_pt = {fReader, "DeepMETResponseTune_pt"};
  */
  //Electron
  TTreeReaderValue<UInt_t> nElectron = {fReader, "nElectron"};
  TTreeReaderArray<Float_t> Electron_dEscaleDown = {fReader, "Electron_dEscaleDown"};
  TTreeReaderArray<Float_t> Electron_dEscaleUp = {fReader, "Electron_dEscaleUp"};
  TTreeReaderArray<Float_t> Electron_dEsigmaDown = {fReader, "Electron_dEsigmaDown"};
  TTreeReaderArray<Float_t> Electron_dEsigmaUp = {fReader, "Electron_dEsigmaUp"};
  TTreeReaderArray<Float_t> Electron_deltaEtaSC = {fReader, "Electron_deltaEtaSC"};
  TTreeReaderArray<Float_t> Electron_dr03EcalRecHitSumEt = {fReader, "Electron_dr03EcalRecHitSumEt"};
  TTreeReaderArray<Float_t> Electron_dr03HcalDepth1TowerSumEt = {fReader, "Electron_dr03HcalDepth1TowerSumEt"};
  TTreeReaderArray<Float_t> Electron_dr03TkSumPt = {fReader, "Electron_dr03TkSumPt"};
  TTreeReaderArray<Float_t> Electron_dr03TkSumPtHEEP = {fReader, "Electron_dr03TkSumPtHEEP"};
  TTreeReaderArray<Float_t> Electron_dxy = {fReader, "Electron_dxy"};
  TTreeReaderArray<Float_t> Electron_dxyErr = {fReader, "Electron_dxyErr"};
  TTreeReaderArray<Float_t> Electron_dz = {fReader, "Electron_dz"};
  TTreeReaderArray<Float_t> Electron_dzErr = {fReader, "Electron_dzErr"};
  TTreeReaderArray<Float_t> Electron_eCorr = {fReader, "Electron_eCorr"};
  TTreeReaderArray<Float_t> Electron_eInvMinusPInv = {fReader, "Electron_eInvMinusPInv"};
  TTreeReaderArray<Float_t> Electron_energyErr = {fReader, "Electron_energyErr"};
  TTreeReaderArray<Float_t> Electron_eta = {fReader, "Electron_eta"};
  TTreeReaderArray<Float_t> Electron_hoe = {fReader, "Electron_hoe"};
  TTreeReaderArray<Float_t> Electron_ip3d = {fReader, "Electron_ip3d"};
  TTreeReaderArray<Float_t> Electron_jetPtRelv2 = {fReader, "Electron_jetPtRelv2"};
  TTreeReaderArray<Float_t> Electron_jetRelIso = {fReader, "Electron_jetRelIso"};
  TTreeReaderArray<Float_t> Electron_mass = {fReader, "Electron_mass"};
  TTreeReaderArray<Float_t> Electron_miniPFRelIso_all = {fReader, "Electron_miniPFRelIso_all"};
  TTreeReaderArray<Float_t> Electron_miniPFRelIso_chg = {fReader, "Electron_miniPFRelIso_chg"};
  TTreeReaderArray<Float_t> Electron_mvaFall17V2Iso = {fReader, "Electron_mvaFall17V2Iso"};
  TTreeReaderArray<Float_t> Electron_mvaFall17V2noIso = {fReader, "Electron_mvaFall17V2noIso"};
  TTreeReaderArray<Float_t> Electron_pfRelIso03_all = {fReader, "Electron_pfRelIso03_all"};
  TTreeReaderArray<Float_t> Electron_pfRelIso03_chg = {fReader, "Electron_pfRelIso03_chg"};
  TTreeReaderArray<Float_t> Electron_phi = {fReader, "Electron_phi"};
  TTreeReaderArray<Float_t> Electron_pt = {fReader, "Electron_pt"};
  TTreeReaderArray<Float_t> Electron_r9 = {fReader, "Electron_r9"};
  TTreeReaderArray<Float_t> Electron_scEtOverPt = {fReader, "Electron_scEtOverPt"};
  TTreeReaderArray<Float_t> Electron_sieie = {fReader, "Electron_sieie"};
  TTreeReaderArray<Float_t> Electron_sip3d = {fReader, "Electron_sip3d"};
  TTreeReaderArray<Float_t> Electron_mvaTTH = {fReader, "Electron_mvaTTH"};
  TTreeReaderArray<Int_t> Electron_charge = {fReader, "Electron_charge"};
  TTreeReaderArray<Int_t> Electron_cutBased = {fReader, "Electron_cutBased"};
  TTreeReaderArray<Int_t> Electron_jetIdx = {fReader, "Electron_jetIdx"};
  TTreeReaderArray<Int_t> Electron_pdgId = {fReader, "Electron_pdgId"};
  TTreeReaderArray<Int_t> Electron_photonIdx = {fReader, "Electron_photonIdx"};
  TTreeReaderArray<Int_t> Electron_tightCharge = {fReader, "Electron_tightCharge"};
  TTreeReaderArray<Int_t> Electron_vidNestedWPBitmap = {fReader, "Electron_vidNestedWPBitmap"};
  TTreeReaderArray<Int_t> Electron_vidNestedWPBitmapHEEP = {fReader, "Electron_vidNestedWPBitmapHEEP"};
  TTreeReaderArray<Bool_t> Electron_convVeto = {fReader, "Electron_convVeto"};
  TTreeReaderArray<Bool_t> Electron_cutBased_HEEP = {fReader, "Electron_cutBased_HEEP"};
  TTreeReaderArray<Bool_t> Electron_isPFcand = {fReader, "Electron_isPFcand"};
  TTreeReaderArray<UChar_t> Electron_jetNDauCharged = {fReader, "Electron_jetNDauCharged"};
  TTreeReaderArray<UChar_t> Electron_lostHits = {fReader, "Electron_lostHits"};
  TTreeReaderArray<Bool_t> Electron_mvaFall17V2Iso_WP80 = {fReader, "Electron_mvaFall17V2Iso_WP80"};
  TTreeReaderArray<Bool_t> Electron_mvaFall17V2Iso_WP90 = {fReader, "Electron_mvaFall17V2Iso_WP90"};
  TTreeReaderArray<Bool_t> Electron_mvaFall17V2Iso_WPL = {fReader, "Electron_mvaFall17V2Iso_WPL"};
  TTreeReaderArray<Bool_t> Electron_mvaFall17V2noIso_WP80 = {fReader, "Electron_mvaFall17V2noIso_WP80"};
  TTreeReaderArray<Bool_t> Electron_mvaFall17V2noIso_WP90 = {fReader, "Electron_mvaFall17V2noIso_WP90"};
  TTreeReaderArray<Bool_t> Electron_mvaFall17V2noIso_WPL = {fReader, "Electron_mvaFall17V2noIso_WPL"};
  TTreeReaderArray<UChar_t> Electron_seedGain = {fReader, "Electron_seedGain"};

  //FatJet
  /*
  TTreeReaderValue<UInt_t> nFatJet = {fReader, "nFatJet"};
  TTreeReaderArray<Float_t> FatJet_area = {fReader, "FatJet_area"};
  TTreeReaderArray<Float_t> FatJet_btagCSVV2 = {fReader, "FatJet_btagCSVV2"};
  TTreeReaderArray<Float_t> FatJet_btagDDBvLV2 = {fReader, "FatJet_btagDDBvLV2"};
  TTreeReaderArray<Float_t> FatJet_btagDDCvBV2 = {fReader, "FatJet_btagDDCvBV2"};
  TTreeReaderArray<Float_t> FatJet_btagDDCvLV2 = {fReader, "FatJet_btagDDCvLV2"};
  TTreeReaderArray<Float_t> FatJet_btagDeepB = {fReader, "FatJet_btagDeepB"};
  TTreeReaderArray<Float_t> FatJet_btagHbb = {fReader, "FatJet_btagHbb"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_H4qvsQCD = {fReader, "FatJet_deepTagMD_H4qvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_HbbvsQCD = {fReader, "FatJet_deepTagMD_HbbvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_TvsQCD = {fReader, "FatJet_deepTagMD_TvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_WvsQCD = {fReader, "FatJet_deepTagMD_WvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_ZHbbvsQCD = {fReader, "FatJet_deepTagMD_ZHbbvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_ZHccvsQCD = {fReader, "FatJet_deepTagMD_ZHccvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_ZbbvsQCD = {fReader, "FatJet_deepTagMD_ZbbvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_ZvsQCD = {fReader, "FatJet_deepTagMD_ZvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_bbvsLight = {fReader, "FatJet_deepTagMD_bbvsLight"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_ccvsLight = {fReader, "FatJet_deepTagMD_ccvsLight"};
  TTreeReaderArray<Float_t> FatJet_deepTag_H = {fReader, "FatJet_deepTag_H"};
  TTreeReaderArray<Float_t> FatJet_deepTag_QCD = {fReader, "FatJet_deepTag_QCD"};
  TTreeReaderArray<Float_t> FatJet_deepTag_QCDothers = {fReader, "FatJet_deepTag_QCDothers"};
  TTreeReaderArray<Float_t> FatJet_deepTag_TvsQCD = {fReader, "FatJet_deepTag_TvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTag_WvsQCD = {fReader, "FatJet_deepTag_WvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTag_ZvsQCD = {fReader, "FatJet_deepTag_ZvsQCD"};
  TTreeReaderArray<Float_t> FatJet_eta = {fReader, "FatJet_eta"};
  TTreeReaderArray<Float_t> FatJet_mass = {fReader, "FatJet_mass"};
  TTreeReaderArray<Float_t> FatJet_msoftdrop = {fReader, "FatJet_msoftdrop"};
  TTreeReaderArray<Float_t> FatJet_n2b1 = {fReader, "FatJet_n2b1"};
  TTreeReaderArray<Float_t> FatJet_n3b1 = {fReader, "FatJet_n3b1"};
  TTreeReaderArray<Float_t> FatJet_particleNetMD_QCD = {fReader, "FatJet_particleNetMD_QCD"};
  TTreeReaderArray<Float_t> FatJet_particleNetMD_Xbb = {fReader, "FatJet_particleNetMD_Xbb"};
  TTreeReaderArray<Float_t> FatJet_particleNetMD_Xcc = {fReader, "FatJet_particleNetMD_Xcc"};
  TTreeReaderArray<Float_t> FatJet_particleNetMD_Xqq = {fReader, "FatJet_particleNetMD_Xqq"};
  TTreeReaderArray<Float_t> FatJet_particleNet_H4qvsQCD = {fReader, "FatJet_particleNet_H4qvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_HbbvsQCD = {fReader, "FatJet_particleNet_HbbvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_HccvsQCD = {fReader, "FatJet_particleNet_HccvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_QCD = {fReader, "FatJet_particleNet_QCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_TvsQCD = {fReader, "FatJet_particleNet_TvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_WvsQCD = {fReader, "FatJet_particleNet_WvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_ZvsQCD = {fReader, "FatJet_particleNet_ZvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_mass = {fReader, "FatJet_particleNet_mass"};
  TTreeReaderArray<Float_t> FatJet_phi = {fReader, "FatJet_phi"};
  TTreeReaderArray<Float_t> FatJet_pt = {fReader, "FatJet_pt"};
  TTreeReaderArray<Float_t> FatJet_rawFactor = {fReader, "FatJet_rawFactor"};
  TTreeReaderArray<Float_t> FatJet_tau1 = {fReader, "FatJet_tau1"};
  TTreeReaderArray<Float_t> FatJet_tau2 = {fReader, "FatJet_tau2"};
  TTreeReaderArray<Float_t> FatJet_tau3 = {fReader, "FatJet_tau3"};
  TTreeReaderArray<Float_t> FatJet_tau4 = {fReader, "FatJet_tau4"};
  TTreeReaderArray<Float_t> FatJet_lsf3 = {fReader, "FatJet_lsf3"};
  TTreeReaderArray<Int_t> FatJet_jetId = {fReader, "FatJet_jetId"};
  TTreeReaderArray<Int_t> FatJet_subJetIdx1 = {fReader, "FatJet_subJetIdx1"};
  TTreeReaderArray<Int_t> FatJet_subJetIdx2 = {fReader, "FatJet_subJetIdx2"};
  TTreeReaderArray<Int_t> FatJet_electronIdx3SJ = {fReader, "FatJet_electronIdx3SJ"};
  TTreeReaderArray<Int_t> FatJet_muonIdx3SJ = {fReader, "FatJet_muonIdx3SJ"};
  TTreeReaderArray<UChar_t> FatJet_nConstituents = {fReader, "FatJet_nConstituents"};
  */
  //FsrPhoton
  /*
  TTreeReaderValue<UInt_t> nFsrPhoton = {fReader, "nFsrPhoton"};
  TTreeReaderArray<Float_t> FsrPhoton_dROverEt2 = {fReader, "FsrPhoton_dROverEt2"};
  TTreeReaderArray<Float_t> FsrPhoton_eta = {fReader, "FsrPhoton_eta"};
  TTreeReaderArray<Float_t> FsrPhoton_phi = {fReader, "FsrPhoton_phi"};
  TTreeReaderArray<Float_t> FsrPhoton_pt = {fReader, "FsrPhoton_pt"};
  TTreeReaderArray<Float_t> FsrPhoton_relIso03 = {fReader, "FsrPhoton_relIso03"};
  TTreeReaderArray<Int_t> FsrPhoton_muonIdx = {fReader, "FsrPhoton_muonIdx"};
  */
  
  //GenJetAKB (read using fReader_MC)
  /*
  TTreeReaderValue<UInt_t> nGenJetAK8 = {fReader, "nGenJetAK8"};
  TTreeReaderArray<Float_t> GenJetAK8_eta = {fReader, "GenJetAK8_eta"};
  TTreeReaderArray<Float_t> GenJetAK8_mass = {fReader, "GenJetAK8_mass"};
  TTreeReaderArray<Float_t> GenJetAK8_phi = {fReader, "GenJetAK8_phi"};
  TTreeReaderArray<Float_t> GenJetAK8_pt = {fReader, "GenJetAK8_pt"};
  
  //GenJet
  TTreeReaderValue<UInt_t> nGenJet = {fReader, "nGenJet"};
  TTreeReaderArray<Float_t> GenJet_eta = {fReader, "GenJet_eta"};
  TTreeReaderArray<Float_t> GenJet_mass = {fReader, "GenJet_mass"};
  TTreeReaderArray<Float_t> GenJet_phi = {fReader, "GenJet_phi"};
  TTreeReaderArray<Float_t> GenJet_pt = {fReader, "GenJet_pt"};
  */
  //GenParticles (read using fReader_MC)
  /*
  TTreeReaderValue<UInt_t> nGenPart = {fReader, "nGenPart"};
  TTreeReaderArray<Float_t> GenPart_eta = {fReader, "GenPart_eta"};
  TTreeReaderArray<Float_t> GenPart_mass = {fReader, "GenPart_mass"};
  TTreeReaderArray<Float_t> GenPart_phi = {fReader, "GenPart_phi"};
  TTreeReaderArray<Float_t> GenPart_pt = {fReader, "GenPart_pt"};
  TTreeReaderArray<Int_t> GenPart_genPartIdxMother = {fReader, "GenPart_genPartIdxMother"};
  TTreeReaderArray<Int_t> GenPart_pdgId = {fReader, "GenPart_pdgId"};
  TTreeReaderArray<Int_t> GenPart_status = {fReader, "GenPart_status"};
  TTreeReaderArray<Int_t> GenPart_statusFlags = {fReader, "GenPart_statusFlags"};
  */
  //Other MC variables (read using fReader_MC)
  /*
  TTreeReaderValue<UInt_t> nSubGenJetAK8 = {fReader, "nSubGenJetAK8"};
  TTreeReaderArray<Float_t> SubGenJetAK8_eta = {fReader, "SubGenJetAK8_eta"};
  TTreeReaderArray<Float_t> SubGenJetAK8_mass = {fReader, "SubGenJetAK8_mass"};
  TTreeReaderArray<Float_t> SubGenJetAK8_phi = {fReader, "SubGenJetAK8_phi"};
  TTreeReaderArray<Float_t> SubGenJetAK8_pt = {fReader, "SubGenJetAK8_pt"};
  TTreeReaderValue<Float_t> Generator_binvar = {fReader, "Generator_binvar"};
  TTreeReaderValue<Float_t> Generator_scalePDF = {fReader, "Generator_scalePDF"};
  TTreeReaderValue<Float_t> Generator_weight = {fReader, "Generator_weight"};
  TTreeReaderValue<Float_t> Generator_x1 = {fReader, "Generator_x1"};
  TTreeReaderValue<Float_t> Generator_x2 = {fReader, "Generator_x2"};
  TTreeReaderValue<Float_t> Generator_xpdf1 = {fReader, "Generator_xpdf1"};
  TTreeReaderValue<Float_t> Generator_xpdf2 = {fReader, "Generator_xpdf2"};
  TTreeReaderValue<Int_t> Generator_id1 = {fReader, "Generator_id1"};
  TTreeReaderValue<Int_t> Generator_id2 = {fReader, "Generator_id2"};
  TTreeReaderValue<Float_t> GenVtx_x = {fReader, "GenVtx_x"};
  TTreeReaderValue<Float_t> GenVtx_y = {fReader, "GenVtx_y"};
  TTreeReaderValue<Float_t> GenVtx_z = {fReader, "GenVtx_z"};
  TTreeReaderValue<UInt_t> nGenVisTau = {fReader, "nGenVisTau"};
  TTreeReaderArray<Float_t> GenVisTau_eta = {fReader, "GenVisTau_eta"};
  TTreeReaderArray<Float_t> GenVisTau_mass = {fReader, "GenVisTau_mass"};
  TTreeReaderArray<Float_t> GenVisTau_phi = {fReader, "GenVisTau_phi"};
  TTreeReaderArray<Float_t> GenVisTau_pt = {fReader, "GenVisTau_pt"};
  TTreeReaderArray<Int_t> GenVisTau_charge = {fReader, "GenVisTau_charge"};
  TTreeReaderArray<Int_t> GenVisTau_genPartIdxMother = {fReader, "GenVisTau_genPartIdxMother"};
  TTreeReaderArray<Int_t> GenVisTau_status = {fReader, "GenVisTau_status"};
  TTreeReaderValue<Float_t> genWeight = {fReader, "genWeight"};
  TTreeReaderValue<Float_t> LHEWeight_originalXWGTUP = {fReader, "LHEWeight_originalXWGTUP"};
  TTreeReaderValue<UInt_t> nLHEPdfWeight = {fReader, "nLHEPdfWeight"};
  TTreeReaderArray<Float_t> LHEPdfWeight = {fReader, "LHEPdfWeight"};
  TTreeReaderValue<UInt_t> nLHEReweightingWeight = {fReader, "nLHEReweightingWeight"};
  TTreeReaderArray<Float_t> LHEReweightingWeight = {fReader, "LHEReweightingWeight"};
  TTreeReaderValue<UInt_t> nLHEScaleWeight = {fReader, "nLHEScaleWeight"};
  TTreeReaderArray<Float_t> LHEScaleWeight = {fReader, "LHEScaleWeight"};
  TTreeReaderValue<UInt_t> nPSWeight = {fReader, "nPSWeight"};
  TTreeReaderArray<Float_t> PSWeight = {fReader, "PSWeight"};
  */
  //IsoTrack
  /*
  TTreeReaderValue<UInt_t> nIsoTrack = {fReader, "nIsoTrack"};
  TTreeReaderArray<Float_t> IsoTrack_dxy = {fReader, "IsoTrack_dxy"};
  TTreeReaderArray<Float_t> IsoTrack_dz = {fReader, "IsoTrack_dz"};
  TTreeReaderArray<Float_t> IsoTrack_eta = {fReader, "IsoTrack_eta"};
  TTreeReaderArray<Float_t> IsoTrack_pfRelIso03_all = {fReader, "IsoTrack_pfRelIso03_all"};
  TTreeReaderArray<Float_t> IsoTrack_pfRelIso03_chg = {fReader, "IsoTrack_pfRelIso03_chg"};
  TTreeReaderArray<Float_t> IsoTrack_phi = {fReader, "IsoTrack_phi"};
  TTreeReaderArray<Float_t> IsoTrack_pt = {fReader, "IsoTrack_pt"};
  TTreeReaderArray<Float_t> IsoTrack_miniPFRelIso_all = {fReader, "IsoTrack_miniPFRelIso_all"};
  TTreeReaderArray<Float_t> IsoTrack_miniPFRelIso_chg = {fReader, "IsoTrack_miniPFRelIso_chg"};
  TTreeReaderArray<Int_t> IsoTrack_charge = {fReader, "IsoTrack_charge"};
  TTreeReaderArray<Int_t> IsoTrack_fromPV = {fReader, "IsoTrack_fromPV"};
  TTreeReaderArray<Int_t> IsoTrack_pdgId = {fReader, "IsoTrack_pdgId"};
  TTreeReaderArray<Bool_t> IsoTrack_isHighPurityTrack = {fReader, "IsoTrack_isHighPurityTrack"};
  TTreeReaderArray<Bool_t> IsoTrack_isPFcand = {fReader, "IsoTrack_isPFcand"};
  TTreeReaderArray<Bool_t> IsoTrack_isFromLostTrack = {fReader, "IsoTrack_isFromLostTrack"};
  */
  //Jet
  TTreeReaderValue<UInt_t> nJet = {fReader, "nJet"};
  TTreeReaderArray<Float_t> Jet_area = {fReader, "Jet_area"};
  TTreeReaderArray<Float_t> Jet_btagCSVV2 = {fReader, "Jet_btagCSVV2"};
  TTreeReaderArray<Float_t> Jet_btagDeepB = {fReader, "Jet_btagDeepB"};
  TTreeReaderArray<Float_t> Jet_btagDeepCvB = {fReader, "Jet_btagDeepCvB"};
  TTreeReaderArray<Float_t> Jet_btagDeepCvL = {fReader, "Jet_btagDeepCvL"};
  TTreeReaderArray<Float_t> Jet_btagDeepFlavB = {fReader, "Jet_btagDeepFlavB"};
  TTreeReaderArray<Float_t> Jet_btagDeepFlavCvB = {fReader, "Jet_btagDeepFlavCvB"};
  TTreeReaderArray<Float_t> Jet_btagDeepFlavCvL = {fReader, "Jet_btagDeepFlavCvL"};
  TTreeReaderArray<Float_t> Jet_btagDeepFlavQG = {fReader, "Jet_btagDeepFlavQG"};
  TTreeReaderArray<Float_t> Jet_chEmEF = {fReader, "Jet_chEmEF"};
  TTreeReaderArray<Float_t> Jet_chFPV0EF = {fReader, "Jet_chFPV0EF"};
  TTreeReaderArray<Float_t> Jet_chHEF = {fReader, "Jet_chHEF"};
  TTreeReaderArray<Float_t> Jet_eta = {fReader, "Jet_eta"};
  TTreeReaderArray<Float_t> Jet_hfsigmaEtaEta = {fReader, "Jet_hfsigmaEtaEta"};
  TTreeReaderArray<Float_t> Jet_hfsigmaPhiPhi = {fReader, "Jet_hfsigmaPhiPhi"};
  TTreeReaderArray<Float_t> Jet_mass = {fReader, "Jet_mass"};
  TTreeReaderArray<Float_t> Jet_muEF = {fReader, "Jet_muEF"};
  TTreeReaderArray<Float_t> Jet_muonSubtrFactor = {fReader, "Jet_muonSubtrFactor"};
  TTreeReaderArray<Float_t> Jet_neEmEF = {fReader, "Jet_neEmEF"};
  TTreeReaderArray<Float_t> Jet_neHEF = {fReader, "Jet_neHEF"};
  TTreeReaderArray<Float_t> Jet_phi = {fReader, "Jet_phi"};
  TTreeReaderArray<Float_t> Jet_pt = {fReader, "Jet_pt"};
  TTreeReaderArray<Float_t> Jet_puIdDisc = {fReader, "Jet_puIdDisc"};
  TTreeReaderArray<Float_t> Jet_qgl = {fReader, "Jet_qgl"};
  TTreeReaderArray<Float_t> Jet_rawFactor = {fReader, "Jet_rawFactor"};
  TTreeReaderArray<Float_t> Jet_bRegCorr = {fReader, "Jet_bRegCorr"};
  TTreeReaderArray<Float_t> Jet_bRegRes = {fReader, "Jet_bRegRes"};
  TTreeReaderArray<Float_t> Jet_cRegCorr = {fReader, "Jet_cRegCorr"};
  TTreeReaderArray<Float_t> Jet_cRegRes = {fReader, "Jet_cRegRes"};
  TTreeReaderArray<Int_t> Jet_electronIdx1 = {fReader, "Jet_electronIdx1"};
  TTreeReaderArray<Int_t> Jet_electronIdx2 = {fReader, "Jet_electronIdx2"};
  TTreeReaderArray<Int_t> Jet_hfadjacentEtaStripsSize = {fReader, "Jet_hfadjacentEtaStripsSize"};
  TTreeReaderArray<Int_t> Jet_hfcentralEtaStripSize = {fReader, "Jet_hfcentralEtaStripSize"};
  TTreeReaderArray<Int_t> Jet_jetId = {fReader, "Jet_jetId"};
  TTreeReaderArray<Int_t> Jet_muonIdx1 = {fReader, "Jet_muonIdx1"};
  TTreeReaderArray<Int_t> Jet_muonIdx2 = {fReader, "Jet_muonIdx2"};
  TTreeReaderArray<Int_t> Jet_nElectrons = {fReader, "Jet_nElectrons"};
  TTreeReaderArray<Int_t> Jet_nMuons = {fReader, "Jet_nMuons"};
  TTreeReaderArray<Int_t> Jet_puId = {fReader, "Jet_puId"};
  TTreeReaderArray<UChar_t> Jet_nConstituents = {fReader, "Jet_nConstituents"};

  //L1PreFiringWeight
  /*
  TTreeReaderValue<Float_t> L1PreFiringWeight_Dn = {fReader, "L1PreFiringWeight_Dn"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_ECAL_Dn = {fReader, "L1PreFiringWeight_ECAL_Dn"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_ECAL_Nom = {fReader, "L1PreFiringWeight_ECAL_Nom"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_ECAL_Up = {fReader, "L1PreFiringWeight_ECAL_Up"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_Muon_Nom = {fReader, "L1PreFiringWeight_Muon_Nom"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_Muon_StatDn = {fReader, "L1PreFiringWeight_Muon_StatDn"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_Muon_StatUp = {fReader, "L1PreFiringWeight_Muon_StatUp"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_Muon_SystDn = {fReader, "L1PreFiringWeight_Muon_SystDn"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_Muon_SystUp = {fReader, "L1PreFiringWeight_Muon_SystUp"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_Nom = {fReader, "L1PreFiringWeight_Nom"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_Up = {fReader, "L1PreFiringWeight_Up"};
  */
  //LHE (read using fReader_MC)
  /*
  TTreeReaderValue<Float_t> LHE_HT = {fReader, "LHE_HT"};
  TTreeReaderValue<Float_t> LHE_HTIncoming = {fReader, "LHE_HTIncoming"};
  TTreeReaderValue<Float_t> LHE_Vpt = {fReader, "LHE_Vpt"};
  TTreeReaderValue<Float_t> LHE_AlphaS = {fReader, "LHE_AlphaS"};
  TTreeReaderValue<UChar_t> LHE_Njets = {fReader, "LHE_Njets"};
  TTreeReaderValue<UChar_t> LHE_Nb = {fReader, "LHE_Nb"};
  TTreeReaderValue<UChar_t> LHE_Nc = {fReader, "LHE_Nc"};
  TTreeReaderValue<UChar_t> LHE_Nuds = {fReader, "LHE_Nuds"};
  TTreeReaderValue<UChar_t> LHE_Nglu = {fReader, "LHE_Nglu"};
  TTreeReaderValue<UChar_t> LHE_NpNLO = {fReader, "LHE_NpNLO"};
  TTreeReaderValue<UChar_t> LHE_NpLO = {fReader, "LHE_NpLO"};
  TTreeReaderValue<UInt_t> nLHEPart = {fReader, "nLHEPart"};
  TTreeReaderArray<Float_t> LHEPart_pt = {fReader, "LHEPart_pt"};
  TTreeReaderArray<Float_t> LHEPart_eta = {fReader, "LHEPart_eta"};
  TTreeReaderArray<Float_t> LHEPart_phi = {fReader, "LHEPart_phi"};
  TTreeReaderArray<Float_t> LHEPart_mass = {fReader, "LHEPart_mass"};
  TTreeReaderArray<Float_t> LHEPart_incomingpz = {fReader, "LHEPart_incomingpz"};
  TTreeReaderArray<Int_t> LHEPart_pdgId = {fReader, "LHEPart_pdgId"};
  TTreeReaderArray<Int_t> LHEPart_status = {fReader, "LHEPart_status"};
  TTreeReaderArray<Int_t> LHEPart_spin = {fReader, "LHEPart_spin"};
  */
  //LowPtElectron
  /*
  TTreeReaderValue<UInt_t> nLowPtElectron = {fReader, "nLowPtElectron"};
  TTreeReaderArray<Float_t> LowPtElectron_ID = {fReader, "LowPtElectron_ID"};
  TTreeReaderArray<Float_t> LowPtElectron_convVtxRadius = {fReader, "LowPtElectron_convVtxRadius"};
  TTreeReaderArray<Float_t> LowPtElectron_deltaEtaSC = {fReader, "LowPtElectron_deltaEtaSC"};
  TTreeReaderArray<Float_t> LowPtElectron_dxy = {fReader, "LowPtElectron_dxy"};
  TTreeReaderArray<Float_t> LowPtElectron_dxyErr = {fReader, "LowPtElectron_dxyErr"};
  TTreeReaderArray<Float_t> LowPtElectron_dz = {fReader, "LowPtElectron_dz"};
  TTreeReaderArray<Float_t> LowPtElectron_dzErr = {fReader, "LowPtElectron_dzErr"};
  TTreeReaderArray<Float_t> LowPtElectron_eInvMinusPInv = {fReader, "LowPtElectron_eInvMinusPInv"};
  TTreeReaderArray<Float_t> LowPtElectron_embeddedID = {fReader, "LowPtElectron_embeddedID"};
  TTreeReaderArray<Float_t> LowPtElectron_energyErr = {fReader, "LowPtElectron_energyErr"};
  TTreeReaderArray<Float_t> LowPtElectron_eta = {fReader, "LowPtElectron_eta"};
  TTreeReaderArray<Float_t> LowPtElectron_hoe = {fReader, "LowPtElectron_hoe"};
  TTreeReaderArray<Float_t> LowPtElectron_mass = {fReader, "LowPtElectron_mass"};
  TTreeReaderArray<Float_t> LowPtElectron_miniPFRelIso_all = {fReader, "LowPtElectron_miniPFRelIso_all"};
  TTreeReaderArray<Float_t> LowPtElectron_miniPFRelIso_chg = {fReader, "LowPtElectron_miniPFRelIso_chg"};
  TTreeReaderArray<Float_t> LowPtElectron_phi = {fReader, "LowPtElectron_phi"};
  TTreeReaderArray<Float_t> LowPtElectron_pt = {fReader, "LowPtElectron_pt"};
  TTreeReaderArray<Float_t> LowPtElectron_ptbiased = {fReader, "LowPtElectron_ptbiased"};
  TTreeReaderArray<Float_t> LowPtElectron_r9 = {fReader, "LowPtElectron_r9"};
  TTreeReaderArray<Float_t> LowPtElectron_scEtOverPt = {fReader, "LowPtElectron_scEtOverPt"};
  TTreeReaderArray<Float_t> LowPtElectron_sieie = {fReader, "LowPtElectron_sieie"};
  TTreeReaderArray<Float_t> LowPtElectron_unbiased = {fReader, "LowPtElectron_unbiased"};
  TTreeReaderArray<Int_t> LowPtElectron_charge = {fReader, "LowPtElectron_charge"};
  TTreeReaderArray<Int_t> LowPtElectron_convWP = {fReader, "LowPtElectron_convWP"};
  TTreeReaderArray<Int_t> LowPtElectron_pdgId = {fReader, "LowPtElectron_pdgId"};
  TTreeReaderArray<Bool_t> LowPtElectron_convVeto = {fReader, "LowPtElectron_convVeto"};
  TTreeReaderArray<UChar_t> LowPtElectron_lostHits = {fReader, "LowPtElectron_lostHits"};
  */
  //MET
  // TTreeReaderValue<Float_t> GenMET_phi = {fReader, "GenMET_phi"};
  // TTreeReaderValue<Float_t> GenMET_pt = {fReader, "GenMET_pt"};
  TTreeReaderValue<Float_t> MET_MetUnclustEnUpDeltaX = {fReader, "MET_MetUnclustEnUpDeltaX"};
  TTreeReaderValue<Float_t> MET_MetUnclustEnUpDeltaY = {fReader, "MET_MetUnclustEnUpDeltaY"};
  TTreeReaderValue<Float_t> MET_covXX = {fReader, "MET_covXX"};
  TTreeReaderValue<Float_t> MET_covXY = {fReader, "MET_covXY"};
  TTreeReaderValue<Float_t> MET_covYY = {fReader, "MET_covYY"};
  TTreeReaderValue<Float_t> MET_phi = {fReader, "MET_phi"};
  TTreeReaderValue<Float_t> MET_pt = {fReader, "MET_pt"};
  TTreeReaderValue<Float_t> MET_significance = {fReader, "MET_significance"};
  TTreeReaderValue<Float_t> MET_sumEt = {fReader, "MET_sumEt"};
  TTreeReaderValue<Float_t> MET_sumPtUnclustered = {fReader, "MET_sumPtUnclustered"};

  //Muon
  TTreeReaderValue<UInt_t> nMuon = {fReader, "nMuon"};
  TTreeReaderArray<Float_t> Muon_dxy = {fReader, "Muon_dxy"};
  TTreeReaderArray<Float_t> Muon_dxyErr = {fReader, "Muon_dxyErr"};
  TTreeReaderArray<Float_t> Muon_dxybs = {fReader, "Muon_dxybs"};
  TTreeReaderArray<Float_t> Muon_dz = {fReader, "Muon_dz"};
  TTreeReaderArray<Float_t> Muon_dzErr = {fReader, "Muon_dzErr"};
  TTreeReaderArray<Float_t> Muon_eta = {fReader, "Muon_eta"};
  TTreeReaderArray<Float_t> Muon_ip3d = {fReader, "Muon_ip3d"};
  TTreeReaderArray<Float_t> Muon_jetPtRelv2 = {fReader, "Muon_jetPtRelv2"};
  TTreeReaderArray<Float_t> Muon_jetRelIso = {fReader, "Muon_jetRelIso"};
  TTreeReaderArray<Float_t> Muon_mass = {fReader, "Muon_mass"};
  TTreeReaderArray<Float_t> Muon_miniPFRelIso_all = {fReader, "Muon_miniPFRelIso_all"};
  TTreeReaderArray<Float_t> Muon_miniPFRelIso_chg = {fReader, "Muon_miniPFRelIso_chg"};
  TTreeReaderArray<Float_t> Muon_pfRelIso03_all = {fReader, "Muon_pfRelIso03_all"};
  TTreeReaderArray<Float_t> Muon_pfRelIso03_chg = {fReader, "Muon_pfRelIso03_chg"};
  TTreeReaderArray<Float_t> Muon_pfRelIso04_all = {fReader, "Muon_pfRelIso04_all"};
  TTreeReaderArray<Float_t> Muon_phi = {fReader, "Muon_phi"};
  TTreeReaderArray<Float_t> Muon_pt = {fReader, "Muon_pt"};
  TTreeReaderArray<Float_t> Muon_ptErr = {fReader, "Muon_ptErr"};
  TTreeReaderArray<Float_t> Muon_segmentComp = {fReader, "Muon_segmentComp"};
  TTreeReaderArray<Float_t> Muon_sip3d = {fReader, "Muon_sip3d"};
  TTreeReaderArray<Float_t> Muon_softMva = {fReader, "Muon_softMva"};
  TTreeReaderArray<Float_t> Muon_tkRelIso = {fReader, "Muon_tkRelIso"};
  TTreeReaderArray<Float_t> Muon_tunepRelPt = {fReader, "Muon_tunepRelPt"};
  TTreeReaderArray<Float_t> Muon_mvaLowPt = {fReader, "Muon_mvaLowPt"};
  TTreeReaderArray<Float_t> Muon_mvaTTH = {fReader, "Muon_mvaTTH"};
  TTreeReaderArray<Int_t> Muon_charge = {fReader, "Muon_charge"};
  TTreeReaderArray<Int_t> Muon_jetIdx = {fReader, "Muon_jetIdx"};
  TTreeReaderArray<Int_t> Muon_nStations = {fReader, "Muon_nStations"};
  TTreeReaderArray<Int_t> Muon_nTrackerLayers = {fReader, "Muon_nTrackerLayers"};
  TTreeReaderArray<Int_t> Muon_pdgId = {fReader, "Muon_pdgId"};
  TTreeReaderArray<Int_t> Muon_tightCharge = {fReader, "Muon_tightCharge"};
  TTreeReaderArray<Int_t> Muon_fsrPhotonIdx = {fReader, "Muon_fsrPhotonIdx"};
  TTreeReaderArray<UChar_t> Muon_highPtId = {fReader, "Muon_highPtId"};
  TTreeReaderArray<Bool_t> Muon_highPurity = {fReader, "Muon_highPurity"};
  TTreeReaderArray<Bool_t> Muon_inTimeMuon = {fReader, "Muon_inTimeMuon"};
  TTreeReaderArray<Bool_t> Muon_isGlobal = {fReader, "Muon_isGlobal"};
  TTreeReaderArray<Bool_t> Muon_isPFcand = {fReader, "Muon_isPFcand"};
  TTreeReaderArray<Bool_t> Muon_isStandalone = {fReader, "Muon_isStandalone"};
  TTreeReaderArray<Bool_t> Muon_isTracker = {fReader, "Muon_isTracker"};
  TTreeReaderArray<UChar_t> Muon_jetNDauCharged = {fReader, "Muon_jetNDauCharged"};
  TTreeReaderArray<Bool_t> Muon_looseId = {fReader, "Muon_looseId"};
  TTreeReaderArray<Bool_t> Muon_mediumId = {fReader, "Muon_mediumId"};
  TTreeReaderArray<Bool_t> Muon_mediumPromptId = {fReader, "Muon_mediumPromptId"};
  TTreeReaderArray<UChar_t> Muon_miniIsoId = {fReader, "Muon_miniIsoId"};
  TTreeReaderArray<UChar_t> Muon_multiIsoId = {fReader, "Muon_multiIsoId"};
  TTreeReaderArray<UChar_t> Muon_mvaId = {fReader, "Muon_mvaId"};
  TTreeReaderArray<UChar_t> Muon_mvaLowPtId = {fReader, "Muon_mvaLowPtId"};
  TTreeReaderArray<UChar_t> Muon_pfIsoId = {fReader, "Muon_pfIsoId"};
  TTreeReaderArray<UChar_t> Muon_puppiIsoId = {fReader, "Muon_puppiIsoId"};
  TTreeReaderArray<Bool_t> Muon_softId = {fReader, "Muon_softId"};
  TTreeReaderArray<Bool_t> Muon_softMvaId = {fReader, "Muon_softMvaId"};
  TTreeReaderArray<Bool_t> Muon_tightId = {fReader, "Muon_tightId"};
  TTreeReaderArray<UChar_t> Muon_tkIsoId = {fReader, "Muon_tkIsoId"};
  TTreeReaderArray<Bool_t> Muon_triggerIdLoose = {fReader, "Muon_triggerIdLoose"};

  //Photon
  TTreeReaderValue<UInt_t> nPhoton = {fReader, "nPhoton"};
  TTreeReaderArray<Float_t> Photon_dEscaleDown = {fReader, "Photon_dEscaleDown"};
  TTreeReaderArray<Float_t> Photon_dEscaleUp = {fReader, "Photon_dEscaleUp"};
  TTreeReaderArray<Float_t> Photon_dEsigmaDown = {fReader, "Photon_dEsigmaDown"};
  TTreeReaderArray<Float_t> Photon_dEsigmaUp = {fReader, "Photon_dEsigmaUp"};
  TTreeReaderArray<Float_t> Photon_eCorr = {fReader, "Photon_eCorr"};
  TTreeReaderArray<Float_t> Photon_energyErr = {fReader, "Photon_energyErr"};
  TTreeReaderArray<Float_t> Photon_eta = {fReader, "Photon_eta"};
  TTreeReaderArray<Float_t> Photon_hoe = {fReader, "Photon_hoe"};
  TTreeReaderArray<Float_t> Photon_mass = {fReader, "Photon_mass"};
  TTreeReaderArray<Float_t> Photon_mvaID = {fReader, "Photon_mvaID"};
  TTreeReaderArray<Float_t> Photon_mvaID_Fall17V1p1 = {fReader, "Photon_mvaID_Fall17V1p1"};
  TTreeReaderArray<Float_t> Photon_pfRelIso03_all = {fReader, "Photon_pfRelIso03_all"};
  TTreeReaderArray<Float_t> Photon_pfRelIso03_chg = {fReader, "Photon_pfRelIso03_chg"};
  TTreeReaderArray<Float_t> Photon_phi = {fReader, "Photon_phi"};
  TTreeReaderArray<Float_t> Photon_pt = {fReader, "Photon_pt"};
  TTreeReaderArray<Float_t> Photon_r9 = {fReader, "Photon_r9"};
  TTreeReaderArray<Float_t> Photon_sieie = {fReader, "Photon_sieie"};
  TTreeReaderArray<Int_t> Photon_charge = {fReader, "Photon_charge"};
  TTreeReaderArray<Int_t> Photon_cutBased = {fReader, "Photon_cutBased"};
  TTreeReaderArray<Int_t> Photon_cutBased_Fall17V1Bitmap = {fReader, "Photon_cutBased_Fall17V1Bitmap"};
  TTreeReaderArray<Int_t> Photon_electronIdx = {fReader, "Photon_electronIdx"};
  TTreeReaderArray<Int_t> Photon_jetIdx = {fReader, "Photon_jetIdx"};
  TTreeReaderArray<Int_t> Photon_pdgId = {fReader, "Photon_pdgId"};
  TTreeReaderArray<Int_t> Photon_vidNestedWPBitmap = {fReader, "Photon_vidNestedWPBitmap"};
  TTreeReaderArray<Bool_t> Photon_electronVeto = {fReader, "Photon_electronVeto"};
  TTreeReaderArray<Bool_t> Photon_isScEtaEB = {fReader, "Photon_isScEtaEB"};
  TTreeReaderArray<Bool_t> Photon_isScEtaEE = {fReader, "Photon_isScEtaEE"};
  TTreeReaderArray<Bool_t> Photon_mvaID_WP80 = {fReader, "Photon_mvaID_WP80"};
  TTreeReaderArray<Bool_t> Photon_mvaID_WP90 = {fReader, "Photon_mvaID_WP90"};
  TTreeReaderArray<Bool_t> Photon_pixelSeed = {fReader, "Photon_pixelSeed"};
  TTreeReaderArray<UChar_t> Photon_seedGain = {fReader, "Photon_seedGain"};

  //Pileup (read using fReader_MC)
  /*
  TTreeReaderValue<Float_t> Pileup_nTrueInt = {fReader, "Pileup_nTrueInt"};
  TTreeReaderValue<Float_t> Pileup_pudensity = {fReader, "Pileup_pudensity"};
  TTreeReaderValue<Float_t> Pileup_gpudensity = {fReader, "Pileup_gpudensity"};
  TTreeReaderValue<Int_t> Pileup_nPU = {fReader, "Pileup_nPU"};
  TTreeReaderValue<Int_t> Pileup_sumEOOT = {fReader, "Pileup_sumEOOT"};
  TTreeReaderValue<Int_t> Pileup_sumLOOT = {fReader, "Pileup_sumLOOT"};
  */
  
  //PuppiMet
  TTreeReaderValue<Float_t> PuppiMET_phi = {fReader, "PuppiMET_phi"};
  TTreeReaderValue<Float_t> PuppiMET_phiJERDown = {fReader, "PuppiMET_phiJERDown"};
  TTreeReaderValue<Float_t> PuppiMET_phiJERUp = {fReader, "PuppiMET_phiJERUp"};
  TTreeReaderValue<Float_t> PuppiMET_phiJESDown = {fReader, "PuppiMET_phiJESDown"};
  TTreeReaderValue<Float_t> PuppiMET_phiJESUp = {fReader, "PuppiMET_phiJESUp"};
  TTreeReaderValue<Float_t> PuppiMET_phiUnclusteredDown = {fReader, "PuppiMET_phiUnclusteredDown"};
  TTreeReaderValue<Float_t> PuppiMET_phiUnclusteredUp = {fReader, "PuppiMET_phiUnclusteredUp"};
  TTreeReaderValue<Float_t> PuppiMET_pt = {fReader, "PuppiMET_pt"};
  TTreeReaderValue<Float_t> PuppiMET_ptJERDown = {fReader, "PuppiMET_ptJERDown"};
  TTreeReaderValue<Float_t> PuppiMET_ptJERUp = {fReader, "PuppiMET_ptJERUp"};
  TTreeReaderValue<Float_t> PuppiMET_ptJESDown = {fReader, "PuppiMET_ptJESDown"};
  TTreeReaderValue<Float_t> PuppiMET_ptJESUp = {fReader, "PuppiMET_ptJESUp"};
  TTreeReaderValue<Float_t> PuppiMET_ptUnclusteredDown = {fReader, "PuppiMET_ptUnclusteredDown"};
  TTreeReaderValue<Float_t> PuppiMET_ptUnclusteredUp = {fReader, "PuppiMET_ptUnclusteredUp"};
  TTreeReaderValue<Float_t> PuppiMET_sumEt = {fReader, "PuppiMET_sumEt"};

  //RawMet
  /*
  TTreeReaderValue<Float_t> RawMET_phi = {fReader, "RawMET_phi"};
  TTreeReaderValue<Float_t> RawMET_pt = {fReader, "RawMET_pt"};
  TTreeReaderValue<Float_t> RawMET_sumEt = {fReader, "RawMET_sumEt"};

  //RawPuppiMET
  TTreeReaderValue<Float_t> RawPuppiMET_phi = {fReader, "RawPuppiMET_phi"};
  TTreeReaderValue<Float_t> RawPuppiMET_pt = {fReader, "RawPuppiMET_pt"};
  TTreeReaderValue<Float_t> RawPuppiMET_sumEt = {fReader, "RawPuppiMET_sumEt"};
  */
  //fixedGridRhoFastjet
  /*
  TTreeReaderValue<Float_t> fixedGridRhoFastjetAll = {fReader, "fixedGridRhoFastjetAll"};
  TTreeReaderValue<Float_t> fixedGridRhoFastjetCentral = {fReader, "fixedGridRhoFastjetCentral"};
  TTreeReaderValue<Float_t> fixedGridRhoFastjetCentralCalo = {fReader, "fixedGridRhoFastjetCentralCalo"};
  TTreeReaderValue<Float_t> fixedGridRhoFastjetCentralChargedPileUp = {fReader, "fixedGridRhoFastjetCentralChargedPileUp"};
  TTreeReaderValue<Float_t> fixedGridRhoFastjetCentralNeutral = {fReader, "fixedGridRhoFastjetCentralNeutral"};
  */
  //GenDressedLepton (read using fReader_MC)
  /*
  TTreeReaderValue<UInt_t> nGenDressedLepton = {fReader, "nGenDressedLepton"};
  TTreeReaderArray<Float_t> GenDressedLepton_eta = {fReader, "GenDressedLepton_eta"};
  TTreeReaderArray<Float_t> GenDressedLepton_mass = {fReader, "GenDressedLepton_mass"};
  TTreeReaderArray<Float_t> GenDressedLepton_phi = {fReader, "GenDressedLepton_phi"};
  TTreeReaderArray<Float_t> GenDressedLepton_pt = {fReader, "GenDressedLepton_pt"};
  TTreeReaderArray<Int_t> GenDressedLepton_pdgId = {fReader, "GenDressedLepton_pdgId"};
  TTreeReaderArray<Bool_t> GenDressedLepton_hasTauAnc = {fReader, "GenDressedLepton_hasTauAnc"};
  TTreeReaderValue<UInt_t> nGenIsolatedPhoton = {fReader, "nGenIsolatedPhoton"};
  TTreeReaderArray<Float_t> GenIsolatedPhoton_eta = {fReader, "GenIsolatedPhoton_eta"};
  TTreeReaderArray<Float_t> GenIsolatedPhoton_mass = {fReader, "GenIsolatedPhoton_mass"};
  TTreeReaderArray<Float_t> GenIsolatedPhoton_phi = {fReader, "GenIsolatedPhoton_phi"};
  TTreeReaderArray<Float_t> GenIsolatedPhoton_pt = {fReader, "GenIsolatedPhoton_pt"};

  //SoftActivityJet
  TTreeReaderValue<UInt_t> nSoftActivityJet = {fReader, "nSoftActivityJet"};
  TTreeReaderArray<Float_t> SoftActivityJet_eta = {fReader, "SoftActivityJet_eta"};
  TTreeReaderArray<Float_t> SoftActivityJet_phi = {fReader, "SoftActivityJet_phi"};
  TTreeReaderArray<Float_t> SoftActivityJet_pt = {fReader, "SoftActivityJet_pt"};
  TTreeReaderValue<Float_t> SoftActivityJetHT = {fReader, "SoftActivityJetHT"};
  TTreeReaderValue<Float_t> SoftActivityJetHT10 = {fReader, "SoftActivityJetHT10"};
  TTreeReaderValue<Float_t> SoftActivityJetHT2 = {fReader, "SoftActivityJetHT2"};
  TTreeReaderValue<Float_t> SoftActivityJetHT5 = {fReader, "SoftActivityJetHT5"};
  TTreeReaderValue<Int_t> SoftActivityJetNjets10 = {fReader, "SoftActivityJetNjets10"};
  TTreeReaderValue<Int_t> SoftActivityJetNjets2 = {fReader, "SoftActivityJetNjets2"};
  TTreeReaderValue<Int_t> SoftActivityJetNjets5 = {fReader, "SoftActivityJetNjets5"};

  //SubJet
  TTreeReaderValue<UInt_t> nSubJet = {fReader, "nSubJet"};
  TTreeReaderArray<Float_t> SubJet_btagCSVV2 = {fReader, "SubJet_btagCSVV2"};
  TTreeReaderArray<Float_t> SubJet_btagDeepB = {fReader, "SubJet_btagDeepB"};
  TTreeReaderArray<Float_t> SubJet_eta = {fReader, "SubJet_eta"};
  TTreeReaderArray<Float_t> SubJet_mass = {fReader, "SubJet_mass"};
  TTreeReaderArray<Float_t> SubJet_n2b1 = {fReader, "SubJet_n2b1"};
  TTreeReaderArray<Float_t> SubJet_n3b1 = {fReader, "SubJet_n3b1"};
  TTreeReaderArray<Float_t> SubJet_phi = {fReader, "SubJet_phi"};
  TTreeReaderArray<Float_t> SubJet_pt = {fReader, "SubJet_pt"};
  TTreeReaderArray<Float_t> SubJet_rawFactor = {fReader, "SubJet_rawFactor"};
  TTreeReaderArray<Float_t> SubJet_tau1 = {fReader, "SubJet_tau1"};
  TTreeReaderArray<Float_t> SubJet_tau2 = {fReader, "SubJet_tau2"};
  TTreeReaderArray<Float_t> SubJet_tau3 = {fReader, "SubJet_tau3"};
  TTreeReaderArray<Float_t> SubJet_tau4 = {fReader, "SubJet_tau4"};
  */
  //Tau
  TTreeReaderValue<UInt_t> nTau = {fReader, "nTau"};
  TTreeReaderArray<Float_t> Tau_chargedIso = {fReader, "Tau_chargedIso"};
  TTreeReaderArray<Float_t> Tau_dxy = {fReader, "Tau_dxy"};
  TTreeReaderArray<Float_t> Tau_dz = {fReader, "Tau_dz"};
  TTreeReaderArray<Float_t> Tau_eta = {fReader, "Tau_eta"};
  TTreeReaderArray<Float_t> Tau_leadTkDeltaEta = {fReader, "Tau_leadTkDeltaEta"};
  TTreeReaderArray<Float_t> Tau_leadTkDeltaPhi = {fReader, "Tau_leadTkDeltaPhi"};
  TTreeReaderArray<Float_t> Tau_leadTkPtOverTauPt = {fReader, "Tau_leadTkPtOverTauPt"};
  TTreeReaderArray<Float_t> Tau_mass = {fReader, "Tau_mass"};
  TTreeReaderArray<Float_t> Tau_neutralIso = {fReader, "Tau_neutralIso"};
  TTreeReaderArray<Float_t> Tau_phi = {fReader, "Tau_phi"};
  TTreeReaderArray<Float_t> Tau_photonsOutsideSignalCone = {fReader, "Tau_photonsOutsideSignalCone"};
  TTreeReaderArray<Float_t> Tau_pt = {fReader, "Tau_pt"};
  TTreeReaderArray<Float_t> Tau_puCorr = {fReader, "Tau_puCorr"};
  TTreeReaderArray<Float_t> Tau_rawDeepTau2017v2p1VSe = {fReader, "Tau_rawDeepTau2017v2p1VSe"};
  TTreeReaderArray<Float_t> Tau_rawDeepTau2017v2p1VSjet = {fReader, "Tau_rawDeepTau2017v2p1VSjet"};
  TTreeReaderArray<Float_t> Tau_rawDeepTau2017v2p1VSmu = {fReader, "Tau_rawDeepTau2017v2p1VSmu"};
  TTreeReaderArray<Float_t> Tau_rawIso = {fReader, "Tau_rawIso"};
  TTreeReaderArray<Float_t> Tau_rawIsodR03 = {fReader, "Tau_rawIsodR03"};
  TTreeReaderArray<Int_t> Tau_charge = {fReader, "Tau_charge"};
  TTreeReaderArray<Int_t> Tau_decayMode = {fReader, "Tau_decayMode"};
  TTreeReaderArray<Int_t> Tau_jetIdx = {fReader, "Tau_jetIdx"};
  TTreeReaderArray<Bool_t> Tau_idAntiEleDeadECal = {fReader, "Tau_idAntiEleDeadECal"};
  TTreeReaderArray<UChar_t> Tau_idAntiMu = {fReader, "Tau_idAntiMu"};
  TTreeReaderArray<Bool_t> Tau_idDecayModeOldDMs = {fReader, "Tau_idDecayModeOldDMs"};
  TTreeReaderArray<UChar_t> Tau_idDeepTau2017v2p1VSe = {fReader, "Tau_idDeepTau2017v2p1VSe"};
  TTreeReaderArray<UChar_t> Tau_idDeepTau2017v2p1VSjet = {fReader, "Tau_idDeepTau2017v2p1VSjet"};
  TTreeReaderArray<UChar_t> Tau_idDeepTau2017v2p1VSmu = {fReader, "Tau_idDeepTau2017v2p1VSmu"};

  //tkMET
  /*
  TTreeReaderValue<Float_t> TkMET_phi = {fReader, "TkMET_phi"};
  TTreeReaderValue<Float_t> TkMET_pt = {fReader, "TkMET_pt"};
  TTreeReaderValue<Float_t> TkMET_sumEt = {fReader, "TkMET_sumEt"};
  */
  //TrigObj
  TTreeReaderValue<UInt_t> nTrigObj = {fReader, "nTrigObj"};
  TTreeReaderArray<Float_t> TrigObj_pt = {fReader, "TrigObj_pt"};
  TTreeReaderArray<Float_t> TrigObj_eta = {fReader, "TrigObj_eta"};
  TTreeReaderArray<Float_t> TrigObj_phi = {fReader, "TrigObj_phi"};
  TTreeReaderArray<Float_t> TrigObj_l1pt = {fReader, "TrigObj_l1pt"};
  TTreeReaderArray<Float_t> TrigObj_l1pt_2 = {fReader, "TrigObj_l1pt_2"};
  TTreeReaderArray<Float_t> TrigObj_l2pt = {fReader, "TrigObj_l2pt"};
  TTreeReaderArray<Int_t> TrigObj_id = {fReader, "TrigObj_id"};
  TTreeReaderArray<Int_t> TrigObj_l1iso = {fReader, "TrigObj_l1iso"};
  TTreeReaderArray<Int_t> TrigObj_l1charge = {fReader, "TrigObj_l1charge"};
  TTreeReaderArray<Int_t> TrigObj_filterBits = {fReader, "TrigObj_filterBits"};
  //  TTreeReaderValue<Int_t> genTtbarId = {fReader, "genTtbarId"};

  //PV
  /*
  TTreeReaderValue<UInt_t> nOtherPV = {fReader, "nOtherPV"};
  TTreeReaderArray<Float_t> OtherPV_z = {fReader, "OtherPV_z"};
  TTreeReaderValue<Float_t> PV_ndof = {fReader, "PV_ndof"};
  TTreeReaderValue<Float_t> PV_x = {fReader, "PV_x"};
  TTreeReaderValue<Float_t> PV_y = {fReader, "PV_y"};
  TTreeReaderValue<Float_t> PV_z = {fReader, "PV_z"};
  TTreeReaderValue<Float_t> PV_chi2 = {fReader, "PV_chi2"};
  TTreeReaderValue<Float_t> PV_score = {fReader, "PV_score"};
  TTreeReaderValue<Int_t> PV_npvs = {fReader, "PV_npvs"};
  TTreeReaderValue<Int_t> PV_npvsGood = {fReader, "PV_npvsGood"};
  */
  //SV
  /*
  TTreeReaderValue<UInt_t> nSV = {fReader, "nSV"};
  TTreeReaderArray<Float_t> SV_dlen = {fReader, "SV_dlen"};
  TTreeReaderArray<Float_t> SV_dlenSig = {fReader, "SV_dlenSig"};
  TTreeReaderArray<Float_t> SV_dxy = {fReader, "SV_dxy"};
  TTreeReaderArray<Float_t> SV_dxySig = {fReader, "SV_dxySig"};
  TTreeReaderArray<Float_t> SV_pAngle = {fReader, "SV_pAngle"};
  TTreeReaderArray<Int_t> SV_charge = {fReader, "SV_charge"};
  */
  //------------
  //Some gen variables (read using fReader_MC)
  /*
  TTreeReaderArray<Int_t> boostedTau_genPartIdx = {fReader, "boostedTau_genPartIdx"};
  TTreeReaderArray<UChar_t> boostedTau_genPartFlav = {fReader, "boostedTau_genPartFlav"};
  TTreeReaderArray<Int_t> Electron_genPartIdx = {fReader, "Electron_genPartIdx"};
  TTreeReaderArray<UChar_t> Electron_genPartFlav = {fReader, "Electron_genPartFlav"};
  TTreeReaderArray<Int_t> FatJet_genJetAK8Idx = {fReader, "FatJet_genJetAK8Idx"};
  TTreeReaderArray<Int_t> FatJet_hadronFlavour = {fReader, "FatJet_hadronFlavour"};
  TTreeReaderArray<UChar_t> FatJet_nBHadrons = {fReader, "FatJet_nBHadrons"};
  TTreeReaderArray<UChar_t> FatJet_nCHadrons = {fReader, "FatJet_nCHadrons"};
  TTreeReaderArray<Int_t> GenJetAK8_partonFlavour = {fReader, "GenJetAK8_partonFlavour"};
  TTreeReaderArray<UChar_t> GenJetAK8_hadronFlavour = {fReader, "GenJetAK8_hadronFlavour"};
  TTreeReaderArray<Int_t> GenJet_partonFlavour = {fReader, "GenJet_partonFlavour"};
  TTreeReaderArray<UChar_t> GenJet_hadronFlavour = {fReader, "GenJet_hadronFlavour"};
  TTreeReaderValue<Float_t> GenVtx_t0 = {fReader, "GenVtx_t0"};
  TTreeReaderArray<Int_t> Jet_genJetIdx = {fReader, "Jet_genJetIdx"};
  TTreeReaderArray<Int_t> Jet_hadronFlavour = {fReader, "Jet_hadronFlavour"};
  TTreeReaderArray<Int_t> Jet_partonFlavour = {fReader, "Jet_partonFlavour"};
  TTreeReaderArray<Int_t> LowPtElectron_genPartIdx = {fReader, "LowPtElectron_genPartIdx"};
  TTreeReaderArray<UChar_t> LowPtElectron_genPartFlav = {fReader, "LowPtElectron_genPartFlav"};
  TTreeReaderArray<Int_t> Muon_genPartIdx = {fReader, "Muon_genPartIdx"};
  TTreeReaderArray<UChar_t> Muon_genPartFlav = {fReader, "Muon_genPartFlav"};
  TTreeReaderArray<Int_t> Photon_genPartIdx = {fReader, "Photon_genPartIdx"};
  TTreeReaderArray<UChar_t> Photon_genPartFlav = {fReader, "Photon_genPartFlav"};
  TTreeReaderValue<Float_t> MET_fiducialGenPhi = {fReader, "MET_fiducialGenPhi"};
  TTreeReaderValue<Float_t> MET_fiducialGenPt = {fReader, "MET_fiducialGenPt"};
  */
  //CleanMask
  
  TTreeReaderArray<UChar_t> Electron_cleanmask = {fReader, "Electron_cleanmask"};
  TTreeReaderArray<UChar_t> Jet_cleanmask = {fReader, "Jet_cleanmask"};
  TTreeReaderArray<UChar_t> Muon_cleanmask = {fReader, "Muon_cleanmask"};
  TTreeReaderArray<UChar_t> Photon_cleanmask = {fReader, "Photon_cleanmask"};
  TTreeReaderArray<UChar_t> Tau_cleanmask = {fReader, "Tau_cleanmask"};
  //TTreeReaderArray<Int_t> SubJet_hadronFlavour = {fReader, "SubJet_hadronFlavour"};
  // TTreeReaderArray<UChar_t> SubJet_nBHadrons = {fReader, "SubJet_nBHadrons"};
  // TTreeReaderArray<UChar_t> SubJet_nCHadrons = {fReader, "SubJet_nCHadrons"};
  // TTreeReaderArray<Float_t> SV_chi2 = {fReader, "SV_chi2"};
  // TTreeReaderArray<Float_t> SV_eta = {fReader, "SV_eta"};
  // TTreeReaderArray<Float_t> SV_mass = {fReader, "SV_mass"};
  // TTreeReaderArray<Float_t> SV_ndof = {fReader, "SV_ndof"};
  // TTreeReaderArray<Float_t> SV_phi = {fReader, "SV_phi"};
  // TTreeReaderArray<Float_t> SV_pt = {fReader, "SV_pt"};
  // TTreeReaderArray<Float_t> SV_x = {fReader, "SV_x"};
  // TTreeReaderArray<Float_t> SV_y = {fReader, "SV_y"};
  // TTreeReaderArray<Float_t> SV_z = {fReader, "SV_z"};
  //TTreeReaderArray<UChar_t> SV_ntracks = {fReader, "SV_ntracks"};
  
  // TTreeReaderArray<Int_t> Tau_genPartIdx = {fReader, "Tau_genPartIdx"};
  // TTreeReaderArray<UChar_t> Tau_genPartFlav = {fReader, "Tau_genPartFlav"};
  
  /*
    TTreeReaderValue<Bool_t> L1_AlwaysTrue = {fReader, "L1_AlwaysTrue"};
    TTreeReaderValue<Bool_t> L1_BPTX_AND_Ref1_VME = {fReader, "L1_BPTX_AND_Ref1_VME"};
    TTreeReaderValue<Bool_t> L1_BPTX_AND_Ref3_VME = {fReader, "L1_BPTX_AND_Ref3_VME"};
    TTreeReaderValue<Bool_t> L1_BPTX_AND_Ref4_VME = {fReader, "L1_BPTX_AND_Ref4_VME"};
    TTreeReaderValue<Bool_t> L1_BPTX_BeamGas_B1_VME = {fReader, "L1_BPTX_BeamGas_B1_VME"};
    TTreeReaderValue<Bool_t> L1_BPTX_BeamGas_B2_VME = {fReader, "L1_BPTX_BeamGas_B2_VME"};
    TTreeReaderValue<Bool_t> L1_BPTX_BeamGas_Ref1_VME = {fReader, "L1_BPTX_BeamGas_Ref1_VME"};
    TTreeReaderValue<Bool_t> L1_BPTX_BeamGas_Ref2_VME = {fReader, "L1_BPTX_BeamGas_Ref2_VME"};
    TTreeReaderValue<Bool_t> L1_BPTX_NotOR_VME = {fReader, "L1_BPTX_NotOR_VME"};
    TTreeReaderValue<Bool_t> L1_BPTX_OR_Ref3_VME = {fReader, "L1_BPTX_OR_Ref3_VME"};
    TTreeReaderValue<Bool_t> L1_BPTX_OR_Ref4_VME = {fReader, "L1_BPTX_OR_Ref4_VME"};
    TTreeReaderValue<Bool_t> L1_BPTX_RefAND_VME = {fReader, "L1_BPTX_RefAND_VME"};
    TTreeReaderValue<Bool_t> L1_BptxMinus = {fReader, "L1_BptxMinus"};
    TTreeReaderValue<Bool_t> L1_BptxOR = {fReader, "L1_BptxOR"};
    TTreeReaderValue<Bool_t> L1_BptxPlus = {fReader, "L1_BptxPlus"};
    TTreeReaderValue<Bool_t> L1_BptxXOR = {fReader, "L1_BptxXOR"};
    TTreeReaderValue<Bool_t> L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142 = {fReader, "L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142"};
    TTreeReaderValue<Bool_t> L1_DoubleEG8er2p5_HTT260er = {fReader, "L1_DoubleEG8er2p5_HTT260er"};
    TTreeReaderValue<Bool_t> L1_DoubleEG8er2p5_HTT280er = {fReader, "L1_DoubleEG8er2p5_HTT280er"};
    TTreeReaderValue<Bool_t> L1_DoubleEG8er2p5_HTT300er = {fReader, "L1_DoubleEG8er2p5_HTT300er"};
    TTreeReaderValue<Bool_t> L1_DoubleEG8er2p5_HTT320er = {fReader, "L1_DoubleEG8er2p5_HTT320er"};
    TTreeReaderValue<Bool_t> L1_DoubleEG8er2p5_HTT340er = {fReader, "L1_DoubleEG8er2p5_HTT340er"};
    TTreeReaderValue<Bool_t> L1_DoubleEG_15_10_er2p5 = {fReader, "L1_DoubleEG_15_10_er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleEG_20_10_er2p5 = {fReader, "L1_DoubleEG_20_10_er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleEG_22_10_er2p5 = {fReader, "L1_DoubleEG_22_10_er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleEG_25_12_er2p5 = {fReader, "L1_DoubleEG_25_12_er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleEG_25_14_er2p5 = {fReader, "L1_DoubleEG_25_14_er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleEG_27_14_er2p5 = {fReader, "L1_DoubleEG_27_14_er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleEG_LooseIso20_10_er2p5 = {fReader, "L1_DoubleEG_LooseIso20_10_er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleEG_LooseIso22_10_er2p5 = {fReader, "L1_DoubleEG_LooseIso22_10_er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleEG_LooseIso22_12_er2p5 = {fReader, "L1_DoubleEG_LooseIso22_12_er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleEG_LooseIso25_12_er2p5 = {fReader, "L1_DoubleEG_LooseIso25_12_er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleIsoTau32er2p1 = {fReader, "L1_DoubleIsoTau32er2p1"};
    TTreeReaderValue<Bool_t> L1_DoubleIsoTau34er2p1 = {fReader, "L1_DoubleIsoTau34er2p1"};
    TTreeReaderValue<Bool_t> L1_DoubleIsoTau36er2p1 = {fReader, "L1_DoubleIsoTau36er2p1"};
    TTreeReaderValue<Bool_t> L1_DoubleJet100er2p3_dEta_Max1p6 = {fReader, "L1_DoubleJet100er2p3_dEta_Max1p6"};
    TTreeReaderValue<Bool_t> L1_DoubleJet100er2p5 = {fReader, "L1_DoubleJet100er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleJet112er2p3_dEta_Max1p6 = {fReader, "L1_DoubleJet112er2p3_dEta_Max1p6"};
    TTreeReaderValue<Bool_t> L1_DoubleJet120er2p5 = {fReader, "L1_DoubleJet120er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleJet150er2p5 = {fReader, "L1_DoubleJet150er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5 = {fReader, "L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5"};
    TTreeReaderValue<Bool_t> L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5 = {fReader, "L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5"};
    TTreeReaderValue<Bool_t> L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5 = {fReader, "L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5"};
    TTreeReaderValue<Bool_t> L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5 = {fReader, "L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5"};
    TTreeReaderValue<Bool_t> L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5 = {fReader, "L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5"};
    TTreeReaderValue<Bool_t> L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5 = {fReader, "L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5"};
    TTreeReaderValue<Bool_t> L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp = {fReader, "L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp"};
    TTreeReaderValue<Bool_t> L1_DoubleJet40er2p5 = {fReader, "L1_DoubleJet40er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleJet_100_30_DoubleJet30_Mass_Min620 = {fReader, "L1_DoubleJet_100_30_DoubleJet30_Mass_Min620"};
    TTreeReaderValue<Bool_t> L1_DoubleJet_110_35_DoubleJet35_Mass_Min620 = {fReader, "L1_DoubleJet_110_35_DoubleJet35_Mass_Min620"};
    TTreeReaderValue<Bool_t> L1_DoubleJet_115_40_DoubleJet40_Mass_Min620 = {fReader, "L1_DoubleJet_115_40_DoubleJet40_Mass_Min620"};
    TTreeReaderValue<Bool_t> L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28 = {fReader, "L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28"};
    TTreeReaderValue<Bool_t> L1_DoubleJet_120_45_DoubleJet45_Mass_Min620 = {fReader, "L1_DoubleJet_120_45_DoubleJet45_Mass_Min620"};
    TTreeReaderValue<Bool_t> L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28 = {fReader, "L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28"};
    TTreeReaderValue<Bool_t> L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ = {fReader, "L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ"};
    TTreeReaderValue<Bool_t> L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp = {fReader, "L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp"};
    TTreeReaderValue<Bool_t> L1_DoubleJet_80_30_Mass_Min420_Mu8 = {fReader, "L1_DoubleJet_80_30_Mass_Min420_Mu8"};
    TTreeReaderValue<Bool_t> L1_DoubleJet_90_30_DoubleJet30_Mass_Min620 = {fReader, "L1_DoubleJet_90_30_DoubleJet30_Mass_Min620"};
    TTreeReaderValue<Bool_t> L1_DoubleLooseIsoEG22er2p1 = {fReader, "L1_DoubleLooseIsoEG22er2p1"};
    TTreeReaderValue<Bool_t> L1_DoubleLooseIsoEG24er2p1 = {fReader, "L1_DoubleLooseIsoEG24er2p1"};
    TTreeReaderValue<Bool_t> L1_DoubleMu0 = {fReader, "L1_DoubleMu0"};
    TTreeReaderValue<Bool_t> L1_DoubleMu0_Mass_Min1 = {fReader, "L1_DoubleMu0_Mass_Min1"};
    TTreeReaderValue<Bool_t> L1_DoubleMu0_OQ = {fReader, "L1_DoubleMu0_OQ"};
    TTreeReaderValue<Bool_t> L1_DoubleMu0_SQ = {fReader, "L1_DoubleMu0_SQ"};
    TTreeReaderValue<Bool_t> L1_DoubleMu0_SQ_OS = {fReader, "L1_DoubleMu0_SQ_OS"};
    TTreeReaderValue<Bool_t> L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8 = {fReader, "L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8"};
    TTreeReaderValue<Bool_t> L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 = {fReader, "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4"};
    TTreeReaderValue<Bool_t> L1_DoubleMu0er1p5_SQ = {fReader, "L1_DoubleMu0er1p5_SQ"};
    TTreeReaderValue<Bool_t> L1_DoubleMu0er1p5_SQ_OS = {fReader, "L1_DoubleMu0er1p5_SQ_OS"};
    TTreeReaderValue<Bool_t> L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 = {fReader, "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"};
    TTreeReaderValue<Bool_t> L1_DoubleMu0er1p5_SQ_dR_Max1p4 = {fReader, "L1_DoubleMu0er1p5_SQ_dR_Max1p4"};
    TTreeReaderValue<Bool_t> L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4 = {fReader, "L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4"};
    TTreeReaderValue<Bool_t> L1_DoubleMu0er2p0_SQ_dR_Max1p4 = {fReader, "L1_DoubleMu0er2p0_SQ_dR_Max1p4"};
    TTreeReaderValue<Bool_t> L1_DoubleMu10_SQ = {fReader, "L1_DoubleMu10_SQ"};
    TTreeReaderValue<Bool_t> L1_DoubleMu18er2p1 = {fReader, "L1_DoubleMu18er2p1"};
    TTreeReaderValue<Bool_t> L1_DoubleMu3_OS_DoubleEG7p5Upsilon = {fReader, "L1_DoubleMu3_OS_DoubleEG7p5Upsilon"};
    TTreeReaderValue<Bool_t> L1_DoubleMu3_SQ_ETMHF50_HTT60er = {fReader, "L1_DoubleMu3_SQ_ETMHF50_HTT60er"};
    TTreeReaderValue<Bool_t> L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5 = {fReader, "L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5 = {fReader, "L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5 = {fReader, "L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleMu3_SQ_HTT220er = {fReader, "L1_DoubleMu3_SQ_HTT220er"};
    TTreeReaderValue<Bool_t> L1_DoubleMu3_SQ_HTT240er = {fReader, "L1_DoubleMu3_SQ_HTT240er"};
    TTreeReaderValue<Bool_t> L1_DoubleMu3_SQ_HTT260er = {fReader, "L1_DoubleMu3_SQ_HTT260er"};
    TTreeReaderValue<Bool_t> L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8 = {fReader, "L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8"};
    TTreeReaderValue<Bool_t> L1_DoubleMu4_SQ_EG9er2p5 = {fReader, "L1_DoubleMu4_SQ_EG9er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleMu4_SQ_OS = {fReader, "L1_DoubleMu4_SQ_OS"};
    TTreeReaderValue<Bool_t> L1_DoubleMu4_SQ_OS_dR_Max1p2 = {fReader, "L1_DoubleMu4_SQ_OS_dR_Max1p2"};
    TTreeReaderValue<Bool_t> L1_DoubleMu4p5_SQ_OS = {fReader, "L1_DoubleMu4p5_SQ_OS"};
    TTreeReaderValue<Bool_t> L1_DoubleMu4p5_SQ_OS_dR_Max1p2 = {fReader, "L1_DoubleMu4p5_SQ_OS_dR_Max1p2"};
    TTreeReaderValue<Bool_t> L1_DoubleMu4p5er2p0_SQ_OS = {fReader, "L1_DoubleMu4p5er2p0_SQ_OS"};
    TTreeReaderValue<Bool_t> L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18 = {fReader, "L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18"};
    TTreeReaderValue<Bool_t> L1_DoubleMu5Upsilon_OS_DoubleEG3 = {fReader, "L1_DoubleMu5Upsilon_OS_DoubleEG3"};
    TTreeReaderValue<Bool_t> L1_DoubleMu5_SQ_EG9er2p5 = {fReader, "L1_DoubleMu5_SQ_EG9er2p5"};
    TTreeReaderValue<Bool_t> L1_DoubleMu9_SQ = {fReader, "L1_DoubleMu9_SQ"};
    TTreeReaderValue<Bool_t> L1_DoubleMu_12_5 = {fReader, "L1_DoubleMu_12_5"};
    TTreeReaderValue<Bool_t> L1_DoubleMu_15_5_SQ = {fReader, "L1_DoubleMu_15_5_SQ"};
    TTreeReaderValue<Bool_t> L1_DoubleMu_15_7 = {fReader, "L1_DoubleMu_15_7"};
    TTreeReaderValue<Bool_t> L1_DoubleMu_15_7_Mass_Min1 = {fReader, "L1_DoubleMu_15_7_Mass_Min1"};
    TTreeReaderValue<Bool_t> L1_DoubleMu_15_7_SQ = {fReader, "L1_DoubleMu_15_7_SQ"};
    TTreeReaderValue<Bool_t> L1_DoubleTau70er2p1 = {fReader, "L1_DoubleTau70er2p1"};
    TTreeReaderValue<Bool_t> L1_ETM120 = {fReader, "L1_ETM120"};
    TTreeReaderValue<Bool_t> L1_ETM150 = {fReader, "L1_ETM150"};
    TTreeReaderValue<Bool_t> L1_ETMHF100 = {fReader, "L1_ETMHF100"};
    TTreeReaderValue<Bool_t> L1_ETMHF100_HTT60er = {fReader, "L1_ETMHF100_HTT60er"};
    TTreeReaderValue<Bool_t> L1_ETMHF110 = {fReader, "L1_ETMHF110"};
    TTreeReaderValue<Bool_t> L1_ETMHF110_HTT60er = {fReader, "L1_ETMHF110_HTT60er"};
    TTreeReaderValue<Bool_t> L1_ETMHF110_HTT60er_NotSecondBunchInTrain = {fReader, "L1_ETMHF110_HTT60er_NotSecondBunchInTrain"};
    TTreeReaderValue<Bool_t> L1_ETMHF120 = {fReader, "L1_ETMHF120"};
    TTreeReaderValue<Bool_t> L1_ETMHF120_HTT60er = {fReader, "L1_ETMHF120_HTT60er"};
    TTreeReaderValue<Bool_t> L1_ETMHF120_NotSecondBunchInTrain = {fReader, "L1_ETMHF120_NotSecondBunchInTrain"};
    TTreeReaderValue<Bool_t> L1_ETMHF130 = {fReader, "L1_ETMHF130"};
    TTreeReaderValue<Bool_t> L1_ETMHF130_HTT60er = {fReader, "L1_ETMHF130_HTT60er"};
    TTreeReaderValue<Bool_t> L1_ETMHF140 = {fReader, "L1_ETMHF140"};
    TTreeReaderValue<Bool_t> L1_ETMHF150 = {fReader, "L1_ETMHF150"};
    TTreeReaderValue<Bool_t> L1_ETMHF90_HTT60er = {fReader, "L1_ETMHF90_HTT60er"};
    TTreeReaderValue<Bool_t> L1_ETT1200 = {fReader, "L1_ETT1200"};
    TTreeReaderValue<Bool_t> L1_ETT1600 = {fReader, "L1_ETT1600"};
    TTreeReaderValue<Bool_t> L1_ETT2000 = {fReader, "L1_ETT2000"};
    TTreeReaderValue<Bool_t> L1_FirstBunchAfterTrain = {fReader, "L1_FirstBunchAfterTrain"};
    TTreeReaderValue<Bool_t> L1_FirstBunchBeforeTrain = {fReader, "L1_FirstBunchBeforeTrain"};
    TTreeReaderValue<Bool_t> L1_FirstBunchInTrain = {fReader, "L1_FirstBunchInTrain"};
    TTreeReaderValue<Bool_t> L1_FirstCollisionInOrbit = {fReader, "L1_FirstCollisionInOrbit"};
    TTreeReaderValue<Bool_t> L1_FirstCollisionInTrain = {fReader, "L1_FirstCollisionInTrain"};
    TTreeReaderValue<Bool_t> L1_HCAL_LaserMon_Trig = {fReader, "L1_HCAL_LaserMon_Trig"};
    TTreeReaderValue<Bool_t> L1_HCAL_LaserMon_Veto = {fReader, "L1_HCAL_LaserMon_Veto"};
    TTreeReaderValue<Bool_t> L1_HTT120er = {fReader, "L1_HTT120er"};
    TTreeReaderValue<Bool_t> L1_HTT160er = {fReader, "L1_HTT160er"};
    TTreeReaderValue<Bool_t> L1_HTT200er = {fReader, "L1_HTT200er"};
    TTreeReaderValue<Bool_t> L1_HTT255er = {fReader, "L1_HTT255er"};
    TTreeReaderValue<Bool_t> L1_HTT280er = {fReader, "L1_HTT280er"};
    TTreeReaderValue<Bool_t> L1_HTT280er_QuadJet_70_55_40_35_er2p4 = {fReader, "L1_HTT280er_QuadJet_70_55_40_35_er2p4"};
    TTreeReaderValue<Bool_t> L1_HTT320er = {fReader, "L1_HTT320er"};
    TTreeReaderValue<Bool_t> L1_HTT320er_QuadJet_70_55_40_40_er2p4 = {fReader, "L1_HTT320er_QuadJet_70_55_40_40_er2p4"};
    TTreeReaderValue<Bool_t> L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3 = {fReader, "L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3"};
    TTreeReaderValue<Bool_t> L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3 = {fReader, "L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3"};
    TTreeReaderValue<Bool_t> L1_HTT360er = {fReader, "L1_HTT360er"};
    TTreeReaderValue<Bool_t> L1_HTT400er = {fReader, "L1_HTT400er"};
    TTreeReaderValue<Bool_t> L1_HTT450er = {fReader, "L1_HTT450er"};
    TTreeReaderValue<Bool_t> L1_IsoEG32er2p5_Mt40 = {fReader, "L1_IsoEG32er2p5_Mt40"};
    TTreeReaderValue<Bool_t> L1_IsoEG32er2p5_Mt44 = {fReader, "L1_IsoEG32er2p5_Mt44"};
    TTreeReaderValue<Bool_t> L1_IsoEG32er2p5_Mt48 = {fReader, "L1_IsoEG32er2p5_Mt48"};
    TTreeReaderValue<Bool_t> L1_IsoTau40er2p1_ETMHF100 = {fReader, "L1_IsoTau40er2p1_ETMHF100"};
    TTreeReaderValue<Bool_t> L1_IsoTau40er2p1_ETMHF110 = {fReader, "L1_IsoTau40er2p1_ETMHF110"};
    TTreeReaderValue<Bool_t> L1_IsoTau40er2p1_ETMHF120 = {fReader, "L1_IsoTau40er2p1_ETMHF120"};
    TTreeReaderValue<Bool_t> L1_IsoTau40er2p1_ETMHF90 = {fReader, "L1_IsoTau40er2p1_ETMHF90"};
    TTreeReaderValue<Bool_t> L1_IsolatedBunch = {fReader, "L1_IsolatedBunch"};
    TTreeReaderValue<Bool_t> L1_LastBunchInTrain = {fReader, "L1_LastBunchInTrain"};
    TTreeReaderValue<Bool_t> L1_LastCollisionInTrain = {fReader, "L1_LastCollisionInTrain"};
    TTreeReaderValue<Bool_t> L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3 = {fReader, "L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3"};
    TTreeReaderValue<Bool_t> L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3 = {fReader, "L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3"};
    TTreeReaderValue<Bool_t> L1_LooseIsoEG24er2p1_HTT100er = {fReader, "L1_LooseIsoEG24er2p1_HTT100er"};
    TTreeReaderValue<Bool_t> L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3 = {fReader, "L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3"};
    TTreeReaderValue<Bool_t> L1_LooseIsoEG26er2p1_HTT100er = {fReader, "L1_LooseIsoEG26er2p1_HTT100er"};
    TTreeReaderValue<Bool_t> L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3 = {fReader, "L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3"};
    TTreeReaderValue<Bool_t> L1_LooseIsoEG28er2p1_HTT100er = {fReader, "L1_LooseIsoEG28er2p1_HTT100er"};
    TTreeReaderValue<Bool_t> L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3 = {fReader, "L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3"};
    TTreeReaderValue<Bool_t> L1_LooseIsoEG30er2p1_HTT100er = {fReader, "L1_LooseIsoEG30er2p1_HTT100er"};
    TTreeReaderValue<Bool_t> L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3 = {fReader, "L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3"};
    TTreeReaderValue<Bool_t> L1_MinimumBiasHF0_AND_BptxAND = {fReader, "L1_MinimumBiasHF0_AND_BptxAND"};
    TTreeReaderValue<Bool_t> L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6 = {fReader, "L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6"};
    TTreeReaderValue<Bool_t> L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6 = {fReader, "L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6"};
    TTreeReaderValue<Bool_t> L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6 = {fReader, "L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6"};
    TTreeReaderValue<Bool_t> L1_Mu18er2p1_Tau24er2p1 = {fReader, "L1_Mu18er2p1_Tau24er2p1"};
    TTreeReaderValue<Bool_t> L1_Mu18er2p1_Tau26er2p1 = {fReader, "L1_Mu18er2p1_Tau26er2p1"};
    TTreeReaderValue<Bool_t> L1_Mu20_EG10er2p5 = {fReader, "L1_Mu20_EG10er2p5"};
    TTreeReaderValue<Bool_t> L1_Mu22er2p1_IsoTau32er2p1 = {fReader, "L1_Mu22er2p1_IsoTau32er2p1"};
    TTreeReaderValue<Bool_t> L1_Mu22er2p1_IsoTau34er2p1 = {fReader, "L1_Mu22er2p1_IsoTau34er2p1"};
    TTreeReaderValue<Bool_t> L1_Mu22er2p1_IsoTau36er2p1 = {fReader, "L1_Mu22er2p1_IsoTau36er2p1"};
    TTreeReaderValue<Bool_t> L1_Mu22er2p1_IsoTau40er2p1 = {fReader, "L1_Mu22er2p1_IsoTau40er2p1"};
    TTreeReaderValue<Bool_t> L1_Mu22er2p1_Tau70er2p1 = {fReader, "L1_Mu22er2p1_Tau70er2p1"};
    TTreeReaderValue<Bool_t> L1_Mu3_Jet120er2p5_dR_Max0p4 = {fReader, "L1_Mu3_Jet120er2p5_dR_Max0p4"};
    TTreeReaderValue<Bool_t> L1_Mu3_Jet120er2p5_dR_Max0p8 = {fReader, "L1_Mu3_Jet120er2p5_dR_Max0p8"};
    TTreeReaderValue<Bool_t> L1_Mu3_Jet16er2p5_dR_Max0p4 = {fReader, "L1_Mu3_Jet16er2p5_dR_Max0p4"};
    TTreeReaderValue<Bool_t> L1_Mu3_Jet30er2p5 = {fReader, "L1_Mu3_Jet30er2p5"};
    TTreeReaderValue<Bool_t> L1_Mu3_Jet35er2p5_dR_Max0p4 = {fReader, "L1_Mu3_Jet35er2p5_dR_Max0p4"};
    TTreeReaderValue<Bool_t> L1_Mu3_Jet60er2p5_dR_Max0p4 = {fReader, "L1_Mu3_Jet60er2p5_dR_Max0p4"};
    TTreeReaderValue<Bool_t> L1_Mu3_Jet80er2p5_dR_Max0p4 = {fReader, "L1_Mu3_Jet80er2p5_dR_Max0p4"};
    TTreeReaderValue<Bool_t> L1_Mu3er1p5_Jet100er2p5_ETMHF40 = {fReader, "L1_Mu3er1p5_Jet100er2p5_ETMHF40"};
    TTreeReaderValue<Bool_t> L1_Mu3er1p5_Jet100er2p5_ETMHF50 = {fReader, "L1_Mu3er1p5_Jet100er2p5_ETMHF50"};
    TTreeReaderValue<Bool_t> L1_Mu5_EG23er2p5 = {fReader, "L1_Mu5_EG23er2p5"};
    TTreeReaderValue<Bool_t> L1_Mu5_LooseIsoEG20er2p5 = {fReader, "L1_Mu5_LooseIsoEG20er2p5"};
    TTreeReaderValue<Bool_t> L1_Mu6_DoubleEG10er2p5 = {fReader, "L1_Mu6_DoubleEG10er2p5"};
    TTreeReaderValue<Bool_t> L1_Mu6_DoubleEG12er2p5 = {fReader, "L1_Mu6_DoubleEG12er2p5"};
    TTreeReaderValue<Bool_t> L1_Mu6_DoubleEG15er2p5 = {fReader, "L1_Mu6_DoubleEG15er2p5"};
    TTreeReaderValue<Bool_t> L1_Mu6_DoubleEG17er2p5 = {fReader, "L1_Mu6_DoubleEG17er2p5"};
    TTreeReaderValue<Bool_t> L1_Mu6_HTT240er = {fReader, "L1_Mu6_HTT240er"};
    TTreeReaderValue<Bool_t> L1_Mu6_HTT250er = {fReader, "L1_Mu6_HTT250er"};
    TTreeReaderValue<Bool_t> L1_Mu7_EG23er2p5 = {fReader, "L1_Mu7_EG23er2p5"};
    TTreeReaderValue<Bool_t> L1_Mu7_LooseIsoEG20er2p5 = {fReader, "L1_Mu7_LooseIsoEG20er2p5"};
    TTreeReaderValue<Bool_t> L1_Mu7_LooseIsoEG23er2p5 = {fReader, "L1_Mu7_LooseIsoEG23er2p5"};
    TTreeReaderValue<Bool_t> L1_NotBptxOR = {fReader, "L1_NotBptxOR"};
    TTreeReaderValue<Bool_t> L1_QuadJet36er2p5_IsoTau52er2p1 = {fReader, "L1_QuadJet36er2p5_IsoTau52er2p1"};
    TTreeReaderValue<Bool_t> L1_QuadJet60er2p5 = {fReader, "L1_QuadJet60er2p5"};
    TTreeReaderValue<Bool_t> L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0 = {fReader, "L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0"};
    TTreeReaderValue<Bool_t> L1_QuadMu0 = {fReader, "L1_QuadMu0"};
    TTreeReaderValue<Bool_t> L1_QuadMu0_OQ = {fReader, "L1_QuadMu0_OQ"};
    TTreeReaderValue<Bool_t> L1_QuadMu0_SQ = {fReader, "L1_QuadMu0_SQ"};
    TTreeReaderValue<Bool_t> L1_SecondBunchInTrain = {fReader, "L1_SecondBunchInTrain"};
    TTreeReaderValue<Bool_t> L1_SecondLastBunchInTrain = {fReader, "L1_SecondLastBunchInTrain"};
    TTreeReaderValue<Bool_t> L1_SingleEG10er2p5 = {fReader, "L1_SingleEG10er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleEG15er2p5 = {fReader, "L1_SingleEG15er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleEG26er2p5 = {fReader, "L1_SingleEG26er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleEG34er2p5 = {fReader, "L1_SingleEG34er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleEG36er2p5 = {fReader, "L1_SingleEG36er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleEG38er2p5 = {fReader, "L1_SingleEG38er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleEG40er2p5 = {fReader, "L1_SingleEG40er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleEG42er2p5 = {fReader, "L1_SingleEG42er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleEG45er2p5 = {fReader, "L1_SingleEG45er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleEG50 = {fReader, "L1_SingleEG50"};
    TTreeReaderValue<Bool_t> L1_SingleEG60 = {fReader, "L1_SingleEG60"};
    TTreeReaderValue<Bool_t> L1_SingleEG8er2p5 = {fReader, "L1_SingleEG8er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleIsoEG24er1p5 = {fReader, "L1_SingleIsoEG24er1p5"};
    TTreeReaderValue<Bool_t> L1_SingleIsoEG24er2p1 = {fReader, "L1_SingleIsoEG24er2p1"};
    TTreeReaderValue<Bool_t> L1_SingleIsoEG26er1p5 = {fReader, "L1_SingleIsoEG26er1p5"};
    TTreeReaderValue<Bool_t> L1_SingleIsoEG26er2p1 = {fReader, "L1_SingleIsoEG26er2p1"};
    TTreeReaderValue<Bool_t> L1_SingleIsoEG26er2p5 = {fReader, "L1_SingleIsoEG26er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleIsoEG28er1p5 = {fReader, "L1_SingleIsoEG28er1p5"};
    TTreeReaderValue<Bool_t> L1_SingleIsoEG28er2p1 = {fReader, "L1_SingleIsoEG28er2p1"};
    TTreeReaderValue<Bool_t> L1_SingleIsoEG28er2p5 = {fReader, "L1_SingleIsoEG28er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleIsoEG30er2p1 = {fReader, "L1_SingleIsoEG30er2p1"};
    TTreeReaderValue<Bool_t> L1_SingleIsoEG30er2p5 = {fReader, "L1_SingleIsoEG30er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleIsoEG32er2p1 = {fReader, "L1_SingleIsoEG32er2p1"};
    TTreeReaderValue<Bool_t> L1_SingleIsoEG32er2p5 = {fReader, "L1_SingleIsoEG32er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleIsoEG34er2p5 = {fReader, "L1_SingleIsoEG34er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleJet10erHE = {fReader, "L1_SingleJet10erHE"};
    TTreeReaderValue<Bool_t> L1_SingleJet120 = {fReader, "L1_SingleJet120"};
    TTreeReaderValue<Bool_t> L1_SingleJet120_FWD3p0 = {fReader, "L1_SingleJet120_FWD3p0"};
    TTreeReaderValue<Bool_t> L1_SingleJet120er2p5 = {fReader, "L1_SingleJet120er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleJet12erHE = {fReader, "L1_SingleJet12erHE"};
    TTreeReaderValue<Bool_t> L1_SingleJet140er2p5 = {fReader, "L1_SingleJet140er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleJet140er2p5_ETMHF80 = {fReader, "L1_SingleJet140er2p5_ETMHF80"};
    TTreeReaderValue<Bool_t> L1_SingleJet140er2p5_ETMHF90 = {fReader, "L1_SingleJet140er2p5_ETMHF90"};
    TTreeReaderValue<Bool_t> L1_SingleJet160er2p5 = {fReader, "L1_SingleJet160er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleJet180 = {fReader, "L1_SingleJet180"};
    TTreeReaderValue<Bool_t> L1_SingleJet180er2p5 = {fReader, "L1_SingleJet180er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleJet200 = {fReader, "L1_SingleJet200"};
    TTreeReaderValue<Bool_t> L1_SingleJet20er2p5_NotBptxOR = {fReader, "L1_SingleJet20er2p5_NotBptxOR"};
    TTreeReaderValue<Bool_t> L1_SingleJet20er2p5_NotBptxOR_3BX = {fReader, "L1_SingleJet20er2p5_NotBptxOR_3BX"};
    TTreeReaderValue<Bool_t> L1_SingleJet35 = {fReader, "L1_SingleJet35"};
    TTreeReaderValue<Bool_t> L1_SingleJet35_FWD3p0 = {fReader, "L1_SingleJet35_FWD3p0"};
    TTreeReaderValue<Bool_t> L1_SingleJet35er2p5 = {fReader, "L1_SingleJet35er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleJet43er2p5_NotBptxOR_3BX = {fReader, "L1_SingleJet43er2p5_NotBptxOR_3BX"};
    TTreeReaderValue<Bool_t> L1_SingleJet46er2p5_NotBptxOR_3BX = {fReader, "L1_SingleJet46er2p5_NotBptxOR_3BX"};
    TTreeReaderValue<Bool_t> L1_SingleJet60 = {fReader, "L1_SingleJet60"};
    TTreeReaderValue<Bool_t> L1_SingleJet60_FWD3p0 = {fReader, "L1_SingleJet60_FWD3p0"};
    TTreeReaderValue<Bool_t> L1_SingleJet60er2p5 = {fReader, "L1_SingleJet60er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleJet8erHE = {fReader, "L1_SingleJet8erHE"};
    TTreeReaderValue<Bool_t> L1_SingleJet90 = {fReader, "L1_SingleJet90"};
    TTreeReaderValue<Bool_t> L1_SingleJet90_FWD3p0 = {fReader, "L1_SingleJet90_FWD3p0"};
    TTreeReaderValue<Bool_t> L1_SingleJet90er2p5 = {fReader, "L1_SingleJet90er2p5"};
    TTreeReaderValue<Bool_t> L1_SingleLooseIsoEG28er1p5 = {fReader, "L1_SingleLooseIsoEG28er1p5"};
    TTreeReaderValue<Bool_t> L1_SingleLooseIsoEG30er1p5 = {fReader, "L1_SingleLooseIsoEG30er1p5"};
    TTreeReaderValue<Bool_t> L1_SingleMu0_BMTF = {fReader, "L1_SingleMu0_BMTF"};
    TTreeReaderValue<Bool_t> L1_SingleMu0_DQ = {fReader, "L1_SingleMu0_DQ"};
    TTreeReaderValue<Bool_t> L1_SingleMu0_EMTF = {fReader, "L1_SingleMu0_EMTF"};
    TTreeReaderValue<Bool_t> L1_SingleMu0_OMTF = {fReader, "L1_SingleMu0_OMTF"};
    TTreeReaderValue<Bool_t> L1_SingleMu10er1p5 = {fReader, "L1_SingleMu10er1p5"};
    TTreeReaderValue<Bool_t> L1_SingleMu12_DQ_BMTF = {fReader, "L1_SingleMu12_DQ_BMTF"};
    TTreeReaderValue<Bool_t> L1_SingleMu12_DQ_EMTF = {fReader, "L1_SingleMu12_DQ_EMTF"};
    TTreeReaderValue<Bool_t> L1_SingleMu12_DQ_OMTF = {fReader, "L1_SingleMu12_DQ_OMTF"};
    TTreeReaderValue<Bool_t> L1_SingleMu12er1p5 = {fReader, "L1_SingleMu12er1p5"};
    TTreeReaderValue<Bool_t> L1_SingleMu14er1p5 = {fReader, "L1_SingleMu14er1p5"};
    TTreeReaderValue<Bool_t> L1_SingleMu15_DQ = {fReader, "L1_SingleMu15_DQ"};
    TTreeReaderValue<Bool_t> L1_SingleMu16er1p5 = {fReader, "L1_SingleMu16er1p5"};
    TTreeReaderValue<Bool_t> L1_SingleMu18 = {fReader, "L1_SingleMu18"};
    TTreeReaderValue<Bool_t> L1_SingleMu18er1p5 = {fReader, "L1_SingleMu18er1p5"};
    TTreeReaderValue<Bool_t> L1_SingleMu20 = {fReader, "L1_SingleMu20"};
    TTreeReaderValue<Bool_t> L1_SingleMu22 = {fReader, "L1_SingleMu22"};
    TTreeReaderValue<Bool_t> L1_SingleMu22_BMTF = {fReader, "L1_SingleMu22_BMTF"};
    TTreeReaderValue<Bool_t> L1_SingleMu22_EMTF = {fReader, "L1_SingleMu22_EMTF"};
    TTreeReaderValue<Bool_t> L1_SingleMu22_OMTF = {fReader, "L1_SingleMu22_OMTF"};
    TTreeReaderValue<Bool_t> L1_SingleMu25 = {fReader, "L1_SingleMu25"};
    TTreeReaderValue<Bool_t> L1_SingleMu3 = {fReader, "L1_SingleMu3"};
    TTreeReaderValue<Bool_t> L1_SingleMu5 = {fReader, "L1_SingleMu5"};
    TTreeReaderValue<Bool_t> L1_SingleMu6er1p5 = {fReader, "L1_SingleMu6er1p5"};
    TTreeReaderValue<Bool_t> L1_SingleMu7 = {fReader, "L1_SingleMu7"};
    TTreeReaderValue<Bool_t> L1_SingleMu7_DQ = {fReader, "L1_SingleMu7_DQ"};
    TTreeReaderValue<Bool_t> L1_SingleMu7er1p5 = {fReader, "L1_SingleMu7er1p5"};
    TTreeReaderValue<Bool_t> L1_SingleMu8er1p5 = {fReader, "L1_SingleMu8er1p5"};
    TTreeReaderValue<Bool_t> L1_SingleMu9er1p5 = {fReader, "L1_SingleMu9er1p5"};
    TTreeReaderValue<Bool_t> L1_SingleMuCosmics = {fReader, "L1_SingleMuCosmics"};
    TTreeReaderValue<Bool_t> L1_SingleMuCosmics_BMTF = {fReader, "L1_SingleMuCosmics_BMTF"};
    TTreeReaderValue<Bool_t> L1_SingleMuCosmics_EMTF = {fReader, "L1_SingleMuCosmics_EMTF"};
    TTreeReaderValue<Bool_t> L1_SingleMuCosmics_OMTF = {fReader, "L1_SingleMuCosmics_OMTF"};
    TTreeReaderValue<Bool_t> L1_SingleMuOpen = {fReader, "L1_SingleMuOpen"};
    TTreeReaderValue<Bool_t> L1_SingleMuOpen_NotBptxOR = {fReader, "L1_SingleMuOpen_NotBptxOR"};
    TTreeReaderValue<Bool_t> L1_SingleMuOpen_er1p1_NotBptxOR_3BX = {fReader, "L1_SingleMuOpen_er1p1_NotBptxOR_3BX"};
    TTreeReaderValue<Bool_t> L1_SingleMuOpen_er1p4_NotBptxOR_3BX = {fReader, "L1_SingleMuOpen_er1p4_NotBptxOR_3BX"};
    TTreeReaderValue<Bool_t> L1_SingleTau120er2p1 = {fReader, "L1_SingleTau120er2p1"};
    TTreeReaderValue<Bool_t> L1_SingleTau130er2p1 = {fReader, "L1_SingleTau130er2p1"};
    TTreeReaderValue<Bool_t> L1_TOTEM_1 = {fReader, "L1_TOTEM_1"};
    TTreeReaderValue<Bool_t> L1_TOTEM_2 = {fReader, "L1_TOTEM_2"};
    TTreeReaderValue<Bool_t> L1_TOTEM_3 = {fReader, "L1_TOTEM_3"};
    TTreeReaderValue<Bool_t> L1_TOTEM_4 = {fReader, "L1_TOTEM_4"};
    TTreeReaderValue<Bool_t> L1_TripleEG16er2p5 = {fReader, "L1_TripleEG16er2p5"};
    TTreeReaderValue<Bool_t> L1_TripleEG_16_12_8_er2p5 = {fReader, "L1_TripleEG_16_12_8_er2p5"};
    TTreeReaderValue<Bool_t> L1_TripleEG_16_15_8_er2p5 = {fReader, "L1_TripleEG_16_15_8_er2p5"};
    TTreeReaderValue<Bool_t> L1_TripleEG_18_17_8_er2p5 = {fReader, "L1_TripleEG_18_17_8_er2p5"};
    TTreeReaderValue<Bool_t> L1_TripleEG_18_18_12_er2p5 = {fReader, "L1_TripleEG_18_18_12_er2p5"};
    TTreeReaderValue<Bool_t> L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5 = {fReader, "L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5"};
    TTreeReaderValue<Bool_t> L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5 = {fReader, "L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5"};
    TTreeReaderValue<Bool_t> L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5 = {fReader, "L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5"};
    TTreeReaderValue<Bool_t> L1_TripleMu0 = {fReader, "L1_TripleMu0"};
    TTreeReaderValue<Bool_t> L1_TripleMu0_OQ = {fReader, "L1_TripleMu0_OQ"};
    TTreeReaderValue<Bool_t> L1_TripleMu0_SQ = {fReader, "L1_TripleMu0_SQ"};
    TTreeReaderValue<Bool_t> L1_TripleMu3 = {fReader, "L1_TripleMu3"};
    TTreeReaderValue<Bool_t> L1_TripleMu3_SQ = {fReader, "L1_TripleMu3_SQ"};
    TTreeReaderValue<Bool_t> L1_TripleMu_5SQ_3SQ_0OQ = {fReader, "L1_TripleMu_5SQ_3SQ_0OQ"};
    TTreeReaderValue<Bool_t> L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9 = {fReader, "L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9"};
    TTreeReaderValue<Bool_t> L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9 = {fReader, "L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"};
    TTreeReaderValue<Bool_t> L1_TripleMu_5_3_3 = {fReader, "L1_TripleMu_5_3_3"};
    TTreeReaderValue<Bool_t> L1_TripleMu_5_3_3_SQ = {fReader, "L1_TripleMu_5_3_3_SQ"};
    TTreeReaderValue<Bool_t> L1_TripleMu_5_3p5_2p5 = {fReader, "L1_TripleMu_5_3p5_2p5"};
    TTreeReaderValue<Bool_t> L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17 = {fReader, "L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17"};
    TTreeReaderValue<Bool_t> L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17 = {fReader, "L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17"};
    TTreeReaderValue<Bool_t> L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17 = {fReader, "L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17"};
    TTreeReaderValue<Bool_t> L1_TripleMu_5_5_3 = {fReader, "L1_TripleMu_5_5_3"};
    TTreeReaderValue<Bool_t> L1_UnpairedBunchBptxMinus = {fReader, "L1_UnpairedBunchBptxMinus"};
    TTreeReaderValue<Bool_t> L1_UnpairedBunchBptxPlus = {fReader, "L1_UnpairedBunchBptxPlus"};
    TTreeReaderValue<Bool_t> L1_ZeroBias = {fReader, "L1_ZeroBias"};
    TTreeReaderValue<Bool_t> L1_ZeroBias_copy = {fReader, "L1_ZeroBias_copy"};
    TTreeReaderValue<Bool_t> L1_UnprefireableEvent = {fReader, "L1_UnprefireableEvent"};
  */

  //Flags
   TTreeReaderValue<Bool_t> Flag_HBHENoiseFilter = {fReader, "Flag_HBHENoiseFilter"};
   TTreeReaderValue<Bool_t> Flag_HBHENoiseIsoFilter = {fReader, "Flag_HBHENoiseIsoFilter"};
   TTreeReaderValue<Bool_t> Flag_CSCTightHaloFilter = {fReader, "Flag_CSCTightHaloFilter"};
   TTreeReaderValue<Bool_t> Flag_CSCTightHaloTrkMuUnvetoFilter = {fReader, "Flag_CSCTightHaloTrkMuUnvetoFilter"};
   TTreeReaderValue<Bool_t> Flag_CSCTightHalo2015Filter = {fReader, "Flag_CSCTightHalo2015Filter"};
   TTreeReaderValue<Bool_t> Flag_globalTightHalo2016Filter = {fReader, "Flag_globalTightHalo2016Filter"};
   TTreeReaderValue<Bool_t> Flag_globalSuperTightHalo2016Filter = {fReader, "Flag_globalSuperTightHalo2016Filter"};
   TTreeReaderValue<Bool_t> Flag_HcalStripHaloFilter = {fReader, "Flag_HcalStripHaloFilter"};
   TTreeReaderValue<Bool_t> Flag_hcalLaserEventFilter = {fReader, "Flag_hcalLaserEventFilter"};
   TTreeReaderValue<Bool_t> Flag_EcalDeadCellTriggerPrimitiveFilter = {fReader, "Flag_EcalDeadCellTriggerPrimitiveFilter"};
   TTreeReaderValue<Bool_t> Flag_EcalDeadCellBoundaryEnergyFilter = {fReader, "Flag_EcalDeadCellBoundaryEnergyFilter"};
   TTreeReaderValue<Bool_t> Flag_ecalBadCalibFilter = {fReader, "Flag_ecalBadCalibFilter"};
   TTreeReaderValue<Bool_t> Flag_goodVertices = {fReader, "Flag_goodVertices"};
   TTreeReaderValue<Bool_t> Flag_eeBadScFilter = {fReader, "Flag_eeBadScFilter"};
   TTreeReaderValue<Bool_t> Flag_ecalLaserCorrFilter = {fReader, "Flag_ecalLaserCorrFilter"};
   TTreeReaderValue<Bool_t> Flag_trkPOGFilters = {fReader, "Flag_trkPOGFilters"};
   TTreeReaderValue<Bool_t> Flag_chargedHadronTrackResolutionFilter = {fReader, "Flag_chargedHadronTrackResolutionFilter"};
   TTreeReaderValue<Bool_t> Flag_muonBadTrackFilter = {fReader, "Flag_muonBadTrackFilter"};
   TTreeReaderValue<Bool_t> Flag_BadChargedCandidateFilter = {fReader, "Flag_BadChargedCandidateFilter"};
   TTreeReaderValue<Bool_t> Flag_BadPFMuonFilter = {fReader, "Flag_BadPFMuonFilter"};
   TTreeReaderValue<Bool_t> Flag_BadPFMuonDzFilter = {fReader, "Flag_BadPFMuonDzFilter"};
   TTreeReaderValue<Bool_t> Flag_hfNoisyHitsFilter = {fReader, "Flag_hfNoisyHitsFilter"};
   TTreeReaderValue<Bool_t> Flag_BadChargedCandidateSummer16Filter = {fReader, "Flag_BadChargedCandidateSummer16Filter"};
   TTreeReaderValue<Bool_t> Flag_BadPFMuonSummer16Filter = {fReader, "Flag_BadPFMuonSummer16Filter"};
   TTreeReaderValue<Bool_t> Flag_trkPOG_manystripclus53X = {fReader, "Flag_trkPOG_manystripclus53X"};
   TTreeReaderValue<Bool_t> Flag_trkPOG_toomanystripclus53X = {fReader, "Flag_trkPOG_toomanystripclus53X"};
   TTreeReaderValue<Bool_t> Flag_trkPOG_logErrorTooManyClusters = {fReader, "Flag_trkPOG_logErrorTooManyClusters"};
   TTreeReaderValue<Bool_t> Flag_METFilters = {fReader, "Flag_METFilters"};
  // TTreeReaderValue<Bool_t> L1Reco_step = {fReader, "L1Reco_step"};

  //HLT_variables (some of these may not work. Be careful while removing the comments)
  /*
    TTreeReaderValue<Bool_t> HLTriggerFirstPath = {fReader, "HLTriggerFirstPath"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet360_TrimMass30 = {fReader, "HLT_AK8PFJet360_TrimMass30"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet380_TrimMass30 = {fReader, "HLT_AK8PFJet380_TrimMass30"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet400_TrimMass30 = {fReader, "HLT_AK8PFJet400_TrimMass30"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet420_TrimMass30 = {fReader, "HLT_AK8PFJet420_TrimMass30"};
    TTreeReaderValue<Bool_t> HLT_AK8PFHT750_TrimMass50 = {fReader, "HLT_AK8PFHT750_TrimMass50"};
    TTreeReaderValue<Bool_t> HLT_AK8PFHT800_TrimMass50 = {fReader, "HLT_AK8PFHT800_TrimMass50"};
    TTreeReaderValue<Bool_t> HLT_AK8PFHT850_TrimMass50 = {fReader, "HLT_AK8PFHT850_TrimMass50"};
    TTreeReaderValue<Bool_t> HLT_AK8PFHT900_TrimMass50 = {fReader, "HLT_AK8PFHT900_TrimMass50"};
    TTreeReaderValue<Bool_t> HLT_CaloJet500_NoJetID = {fReader, "HLT_CaloJet500_NoJetID"};
    TTreeReaderValue<Bool_t> HLT_CaloJet550_NoJetID = {fReader, "HLT_CaloJet550_NoJetID"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL = {fReader, "HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon = {fReader, "HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon"};
    TTreeReaderValue<Bool_t> HLT_Trimuon5_3p5_2_Upsilon_Muon = {fReader, "HLT_Trimuon5_3p5_2_Upsilon_Muon"};
    TTreeReaderValue<Bool_t> HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon = {fReader, "HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon"};
    TTreeReaderValue<Bool_t> HLT_DoubleEle25_CaloIdL_MW = {fReader, "HLT_DoubleEle25_CaloIdL_MW"};
    TTreeReaderValue<Bool_t> HLT_DoubleEle27_CaloIdL_MW = {fReader, "HLT_DoubleEle27_CaloIdL_MW"};
    TTreeReaderValue<Bool_t> HLT_DoubleEle33_CaloIdL_MW = {fReader, "HLT_DoubleEle33_CaloIdL_MW"};
    TTreeReaderValue<Bool_t> HLT_DoubleEle24_eta2p1_WPTight_Gsf = {fReader, "HLT_DoubleEle24_eta2p1_WPTight_Gsf"};
    TTreeReaderValue<Bool_t> HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350 = {fReader, "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350"};
    TTreeReaderValue<Bool_t> HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350 = {fReader, "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350"};
    TTreeReaderValue<Bool_t> HLT_Ele27_Ele37_CaloIdL_MW = {fReader, "HLT_Ele27_Ele37_CaloIdL_MW"};
    TTreeReaderValue<Bool_t> HLT_Mu27_Ele37_CaloIdL_MW = {fReader, "HLT_Mu27_Ele37_CaloIdL_MW"};
    TTreeReaderValue<Bool_t> HLT_Mu37_Ele27_CaloIdL_MW = {fReader, "HLT_Mu37_Ele27_CaloIdL_MW"};
    TTreeReaderValue<Bool_t> HLT_Mu37_TkMu27 = {fReader, "HLT_Mu37_TkMu27"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu4_3_Bs = {fReader, "HLT_DoubleMu4_3_Bs"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu4_3_Jpsi = {fReader, "HLT_DoubleMu4_3_Jpsi"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu4_JpsiTrk_Displaced = {fReader, "HLT_DoubleMu4_JpsiTrk_Displaced"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu4_LowMassNonResonantTrk_Displaced = {fReader, "HLT_DoubleMu4_LowMassNonResonantTrk_Displaced"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu3_Trk_Tau3mu = {fReader, "HLT_DoubleMu3_Trk_Tau3mu"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu3_TkMu_DsTau3Mu = {fReader, "HLT_DoubleMu3_TkMu_DsTau3Mu"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu4_PsiPrimeTrk_Displaced = {fReader, "HLT_DoubleMu4_PsiPrimeTrk_Displaced"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu4_Mass3p8_DZ_PFHT350 = {fReader, "HLT_DoubleMu4_Mass3p8_DZ_PFHT350"};
    TTreeReaderValue<Bool_t> HLT_Mu3_PFJet40 = {fReader, "HLT_Mu3_PFJet40"};
    TTreeReaderValue<Bool_t> HLT_Mu7p5_L2Mu2_Jpsi = {fReader, "HLT_Mu7p5_L2Mu2_Jpsi"};
    TTreeReaderValue<Bool_t> HLT_Mu7p5_L2Mu2_Upsilon = {fReader, "HLT_Mu7p5_L2Mu2_Upsilon"};
    TTreeReaderValue<Bool_t> HLT_Mu7p5_Track2_Jpsi = {fReader, "HLT_Mu7p5_Track2_Jpsi"};
    TTreeReaderValue<Bool_t> HLT_Mu7p5_Track3p5_Jpsi = {fReader, "HLT_Mu7p5_Track3p5_Jpsi"};
    TTreeReaderValue<Bool_t> HLT_Mu7p5_Track7_Jpsi = {fReader, "HLT_Mu7p5_Track7_Jpsi"};
    TTreeReaderValue<Bool_t> HLT_Mu7p5_Track2_Upsilon = {fReader, "HLT_Mu7p5_Track2_Upsilon"};
    TTreeReaderValue<Bool_t> HLT_Mu7p5_Track3p5_Upsilon = {fReader, "HLT_Mu7p5_Track3p5_Upsilon"};
    TTreeReaderValue<Bool_t> HLT_Mu7p5_Track7_Upsilon = {fReader, "HLT_Mu7p5_Track7_Upsilon"};
    TTreeReaderValue<Bool_t> HLT_Mu3_L1SingleMu5orSingleMu7 = {fReader, "HLT_Mu3_L1SingleMu5orSingleMu7"};
    TTreeReaderValue<Bool_t> HLT_DoublePhoton33_CaloIdL = {fReader, "HLT_DoublePhoton33_CaloIdL"};
    TTreeReaderValue<Bool_t> HLT_DoublePhoton70 = {fReader, "HLT_DoublePhoton70"};
    TTreeReaderValue<Bool_t> HLT_DoublePhoton85 = {fReader, "HLT_DoublePhoton85"};
    TTreeReaderValue<Bool_t> HLT_Ele20_WPTight_Gsf = {fReader, "HLT_Ele20_WPTight_Gsf"};
    TTreeReaderValue<Bool_t> HLT_Ele15_WPLoose_Gsf = {fReader, "HLT_Ele15_WPLoose_Gsf"};
    TTreeReaderValue<Bool_t> HLT_Ele17_WPLoose_Gsf = {fReader, "HLT_Ele17_WPLoose_Gsf"};
    TTreeReaderValue<Bool_t> HLT_Ele20_WPLoose_Gsf = {fReader, "HLT_Ele20_WPLoose_Gsf"};
    TTreeReaderValue<Bool_t> HLT_Ele20_eta2p1_WPLoose_Gsf = {fReader, "HLT_Ele20_eta2p1_WPLoose_Gsf"};
    TTreeReaderValue<Bool_t> HLT_DiEle27_WPTightCaloOnly_L1DoubleEG = {fReader, "HLT_DiEle27_WPTightCaloOnly_L1DoubleEG"};
    TTreeReaderValue<Bool_t> HLT_Ele27_WPTight_Gsf = {fReader, "HLT_Ele27_WPTight_Gsf"};
    TTreeReaderValue<Bool_t> HLT_Ele28_WPTight_Gsf = {fReader, "HLT_Ele28_WPTight_Gsf"};
    TTreeReaderValue<Bool_t> HLT_Ele30_WPTight_Gsf = {fReader, "HLT_Ele30_WPTight_Gsf"};
    TTreeReaderValue<Bool_t> HLT_Ele32_WPTight_Gsf = {fReader, "HLT_Ele32_WPTight_Gsf"};
    TTreeReaderValue<Bool_t> HLT_Ele35_WPTight_Gsf = {fReader, "HLT_Ele35_WPTight_Gsf"};
    TTreeReaderValue<Bool_t> HLT_Ele35_WPTight_Gsf_L1EGMT = {fReader, "HLT_Ele35_WPTight_Gsf_L1EGMT"};
    TTreeReaderValue<Bool_t> HLT_Ele38_WPTight_Gsf = {fReader, "HLT_Ele38_WPTight_Gsf"};
    TTreeReaderValue<Bool_t> HLT_Ele40_WPTight_Gsf = {fReader, "HLT_Ele40_WPTight_Gsf"};
    TTreeReaderValue<Bool_t> HLT_Ele32_WPTight_Gsf_L1DoubleEG = {fReader, "HLT_Ele32_WPTight_Gsf_L1DoubleEG"};
    TTreeReaderValue<Bool_t> HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1 = {fReader, "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1"};
    TTreeReaderValue<Bool_t> HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1 = {fReader, "HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1"};
    TTreeReaderValue<Bool_t> HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1 = {fReader, "HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1"};
    TTreeReaderValue<Bool_t> HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1 = {fReader, "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1"};
    TTreeReaderValue<Bool_t> HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1 = {fReader, "HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1"};
    TTreeReaderValue<Bool_t> HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1 = {fReader, "HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1"};
    TTreeReaderValue<Bool_t> HLT_HT450_Beamspot = {fReader, "HLT_HT450_Beamspot"};
    TTreeReaderValue<Bool_t> HLT_HT300_Beamspot = {fReader, "HLT_HT300_Beamspot"};
    TTreeReaderValue<Bool_t> HLT_ZeroBias_Beamspot = {fReader, "HLT_ZeroBias_Beamspot"};
    TTreeReaderValue<Bool_t> HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1 = {fReader, "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1"};
    TTreeReaderValue<Bool_t> HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1 = {fReader, "HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1"};
    TTreeReaderValue<Bool_t> HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1 = {fReader, "HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1"};
    TTreeReaderValue<Bool_t> HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1 = {fReader, "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1"};
    TTreeReaderValue<Bool_t> HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1 = {fReader, "HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1"};
    TTreeReaderValue<Bool_t> HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1 = {fReader, "HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1"};
    TTreeReaderValue<Bool_t> HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1 = {fReader, "HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1"};
    TTreeReaderValue<Bool_t> HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1 = {fReader, "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1"};
    TTreeReaderValue<Bool_t> HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1 = {fReader, "HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1"};
    TTreeReaderValue<Bool_t> HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1 = {fReader, "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1"};
    TTreeReaderValue<Bool_t> HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1 = {fReader, "HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1"};
    TTreeReaderValue<Bool_t> HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1 = {fReader, "HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1"};
    TTreeReaderValue<Bool_t> HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1 = {fReader, "HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1"};
    TTreeReaderValue<Bool_t> HLT_IsoMu20 = {fReader, "HLT_IsoMu20"};
    TTreeReaderValue<Bool_t> HLT_IsoMu24 = {fReader, "HLT_IsoMu24"};
    TTreeReaderValue<Bool_t> HLT_IsoMu24_eta2p1 = {fReader, "HLT_IsoMu24_eta2p1"};
    TTreeReaderValue<Bool_t> HLT_IsoMu27 = {fReader, "HLT_IsoMu27"};
    TTreeReaderValue<Bool_t> HLT_IsoMu30 = {fReader, "HLT_IsoMu30"};
    TTreeReaderValue<Bool_t> HLT_UncorrectedJetE30_NoBPTX = {fReader, "HLT_UncorrectedJetE30_NoBPTX"};
    TTreeReaderValue<Bool_t> HLT_UncorrectedJetE30_NoBPTX3BX = {fReader, "HLT_UncorrectedJetE30_NoBPTX3BX"};
    TTreeReaderValue<Bool_t> HLT_UncorrectedJetE60_NoBPTX3BX = {fReader, "HLT_UncorrectedJetE60_NoBPTX3BX"};
    TTreeReaderValue<Bool_t> HLT_UncorrectedJetE70_NoBPTX3BX = {fReader, "HLT_UncorrectedJetE70_NoBPTX3BX"};
    TTreeReaderValue<Bool_t> HLT_L1SingleMu18 = {fReader, "HLT_L1SingleMu18"};
    TTreeReaderValue<Bool_t> HLT_L1SingleMu25 = {fReader, "HLT_L1SingleMu25"};
    TTreeReaderValue<Bool_t> HLT_L2Mu10 = {fReader, "HLT_L2Mu10"};
    TTreeReaderValue<Bool_t> HLT_L2Mu10_NoVertex_NoBPTX3BX = {fReader, "HLT_L2Mu10_NoVertex_NoBPTX3BX"};
    TTreeReaderValue<Bool_t> HLT_L2Mu10_NoVertex_NoBPTX = {fReader, "HLT_L2Mu10_NoVertex_NoBPTX"};
    TTreeReaderValue<Bool_t> HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX = {fReader, "HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX"};
    TTreeReaderValue<Bool_t> HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX = {fReader, "HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX"};
    TTreeReaderValue<Bool_t> HLT_L2Mu50 = {fReader, "HLT_L2Mu50"};
    TTreeReaderValue<Bool_t> HLT_L2Mu23NoVtx_2Cha = {fReader, "HLT_L2Mu23NoVtx_2Cha"};
    TTreeReaderValue<Bool_t> HLT_L2Mu23NoVtx_2Cha_CosmicSeed = {fReader, "HLT_L2Mu23NoVtx_2Cha_CosmicSeed"};
    TTreeReaderValue<Bool_t> HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4 = {fReader, "HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4"};
    TTreeReaderValue<Bool_t> HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4 = {fReader, "HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4"};
    TTreeReaderValue<Bool_t> HLT_DoubleL2Mu50 = {fReader, "HLT_DoubleL2Mu50"};
    TTreeReaderValue<Bool_t> HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed = {fReader, "HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed"};
    TTreeReaderValue<Bool_t> HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched = {fReader, "HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched"};
    TTreeReaderValue<Bool_t> HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed = {fReader, "HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed"};
    TTreeReaderValue<Bool_t> HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched = {fReader, "HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched"};
    TTreeReaderValue<Bool_t> HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4 = {fReader, "HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4"};
    TTreeReaderValue<Bool_t> HLT_DoubleL2Mu23NoVtx_2Cha = {fReader, "HLT_DoubleL2Mu23NoVtx_2Cha"};
    TTreeReaderValue<Bool_t> HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched = {fReader, "HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched"};
    TTreeReaderValue<Bool_t> HLT_DoubleL2Mu25NoVtx_2Cha = {fReader, "HLT_DoubleL2Mu25NoVtx_2Cha"};
    TTreeReaderValue<Bool_t> HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched = {fReader, "HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched"};
    TTreeReaderValue<Bool_t> HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4 = {fReader, "HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4"};
    TTreeReaderValue<Bool_t> HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL = {fReader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL"};
    TTreeReaderValue<Bool_t> HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL = {fReader, "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL"};
    TTreeReaderValue<Bool_t> HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = {fReader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"};
    TTreeReaderValue<Bool_t> HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ = {fReader, "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ"};
    TTreeReaderValue<Bool_t> HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 = {fReader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"};
    TTreeReaderValue<Bool_t> HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 = {fReader, "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8"};
    TTreeReaderValue<Bool_t> HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 = {fReader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"};
    TTreeReaderValue<Bool_t> HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 = {fReader, "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8"};
    TTreeReaderValue<Bool_t> HLT_Mu25_TkMu0_Onia = {fReader, "HLT_Mu25_TkMu0_Onia"};
    TTreeReaderValue<Bool_t> HLT_Mu30_TkMu0_Psi = {fReader, "HLT_Mu30_TkMu0_Psi"};
    TTreeReaderValue<Bool_t> HLT_Mu30_TkMu0_Upsilon = {fReader, "HLT_Mu30_TkMu0_Upsilon"};
    TTreeReaderValue<Bool_t> HLT_Mu20_TkMu0_Phi = {fReader, "HLT_Mu20_TkMu0_Phi"};
    TTreeReaderValue<Bool_t> HLT_Mu25_TkMu0_Phi = {fReader, "HLT_Mu25_TkMu0_Phi"};
    TTreeReaderValue<Bool_t> HLT_Mu12 = {fReader, "HLT_Mu12"};
    TTreeReaderValue<Bool_t> HLT_Mu15 = {fReader, "HLT_Mu15"};
    TTreeReaderValue<Bool_t> HLT_Mu20 = {fReader, "HLT_Mu20"};
    TTreeReaderValue<Bool_t> HLT_Mu27 = {fReader, "HLT_Mu27"};
    TTreeReaderValue<Bool_t> HLT_Mu50 = {fReader, "HLT_Mu50"};
    TTreeReaderValue<Bool_t> HLT_Mu55 = {fReader, "HLT_Mu55"};
    TTreeReaderValue<Bool_t> HLT_OldMu100 = {fReader, "HLT_OldMu100"};
    TTreeReaderValue<Bool_t> HLT_TkMu100 = {fReader, "HLT_TkMu100"};
    TTreeReaderValue<Bool_t> HLT_DiPFJetAve40 = {fReader, "HLT_DiPFJetAve40"};
    TTreeReaderValue<Bool_t> HLT_DiPFJetAve60 = {fReader, "HLT_DiPFJetAve60"};
    TTreeReaderValue<Bool_t> HLT_DiPFJetAve80 = {fReader, "HLT_DiPFJetAve80"};
    TTreeReaderValue<Bool_t> HLT_DiPFJetAve140 = {fReader, "HLT_DiPFJetAve140"};
    TTreeReaderValue<Bool_t> HLT_DiPFJetAve200 = {fReader, "HLT_DiPFJetAve200"};
    TTreeReaderValue<Bool_t> HLT_DiPFJetAve260 = {fReader, "HLT_DiPFJetAve260"};
    TTreeReaderValue<Bool_t> HLT_DiPFJetAve320 = {fReader, "HLT_DiPFJetAve320"};
    TTreeReaderValue<Bool_t> HLT_DiPFJetAve400 = {fReader, "HLT_DiPFJetAve400"};
    TTreeReaderValue<Bool_t> HLT_DiPFJetAve500 = {fReader, "HLT_DiPFJetAve500"};
    TTreeReaderValue<Bool_t> HLT_DiPFJetAve60_HFJEC = {fReader, "HLT_DiPFJetAve60_HFJEC"};
    TTreeReaderValue<Bool_t> HLT_DiPFJetAve80_HFJEC = {fReader, "HLT_DiPFJetAve80_HFJEC"};
    TTreeReaderValue<Bool_t> HLT_DiPFJetAve100_HFJEC = {fReader, "HLT_DiPFJetAve100_HFJEC"};
    TTreeReaderValue<Bool_t> HLT_DiPFJetAve160_HFJEC = {fReader, "HLT_DiPFJetAve160_HFJEC"};
    TTreeReaderValue<Bool_t> HLT_DiPFJetAve220_HFJEC = {fReader, "HLT_DiPFJetAve220_HFJEC"};
    TTreeReaderValue<Bool_t> HLT_DiPFJetAve300_HFJEC = {fReader, "HLT_DiPFJetAve300_HFJEC"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet15 = {fReader, "HLT_AK8PFJet15"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet25 = {fReader, "HLT_AK8PFJet25"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet40 = {fReader, "HLT_AK8PFJet40"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet60 = {fReader, "HLT_AK8PFJet60"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet80 = {fReader, "HLT_AK8PFJet80"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet140 = {fReader, "HLT_AK8PFJet140"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet200 = {fReader, "HLT_AK8PFJet200"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet260 = {fReader, "HLT_AK8PFJet260"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet320 = {fReader, "HLT_AK8PFJet320"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet400 = {fReader, "HLT_AK8PFJet400"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet450 = {fReader, "HLT_AK8PFJet450"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet500 = {fReader, "HLT_AK8PFJet500"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet550 = {fReader, "HLT_AK8PFJet550"};
    TTreeReaderValue<Bool_t> HLT_PFJet15 = {fReader, "HLT_PFJet15"};
    TTreeReaderValue<Bool_t> HLT_PFJet25 = {fReader, "HLT_PFJet25"};
    TTreeReaderValue<Bool_t> HLT_PFJet40 = {fReader, "HLT_PFJet40"};
    TTreeReaderValue<Bool_t> HLT_PFJet60 = {fReader, "HLT_PFJet60"};
    TTreeReaderValue<Bool_t> HLT_PFJet80 = {fReader, "HLT_PFJet80"};
    TTreeReaderValue<Bool_t> HLT_PFJet140 = {fReader, "HLT_PFJet140"};
    TTreeReaderValue<Bool_t> HLT_PFJet200 = {fReader, "HLT_PFJet200"};
    TTreeReaderValue<Bool_t> HLT_PFJet260 = {fReader, "HLT_PFJet260"};
    TTreeReaderValue<Bool_t> HLT_PFJet320 = {fReader, "HLT_PFJet320"};
    TTreeReaderValue<Bool_t> HLT_PFJet400 = {fReader, "HLT_PFJet400"};
    TTreeReaderValue<Bool_t> HLT_PFJet450 = {fReader, "HLT_PFJet450"};
    TTreeReaderValue<Bool_t> HLT_PFJet500 = {fReader, "HLT_PFJet500"};
    TTreeReaderValue<Bool_t> HLT_PFJet550 = {fReader, "HLT_PFJet550"};
    TTreeReaderValue<Bool_t> HLT_PFJetFwd15 = {fReader, "HLT_PFJetFwd15"};
    TTreeReaderValue<Bool_t> HLT_PFJetFwd25 = {fReader, "HLT_PFJetFwd25"};
    TTreeReaderValue<Bool_t> HLT_PFJetFwd40 = {fReader, "HLT_PFJetFwd40"};
    TTreeReaderValue<Bool_t> HLT_PFJetFwd60 = {fReader, "HLT_PFJetFwd60"};
    TTreeReaderValue<Bool_t> HLT_PFJetFwd80 = {fReader, "HLT_PFJetFwd80"};
    TTreeReaderValue<Bool_t> HLT_PFJetFwd140 = {fReader, "HLT_PFJetFwd140"};
    TTreeReaderValue<Bool_t> HLT_PFJetFwd200 = {fReader, "HLT_PFJetFwd200"};
    TTreeReaderValue<Bool_t> HLT_PFJetFwd260 = {fReader, "HLT_PFJetFwd260"};
    TTreeReaderValue<Bool_t> HLT_PFJetFwd320 = {fReader, "HLT_PFJetFwd320"};
    TTreeReaderValue<Bool_t> HLT_PFJetFwd400 = {fReader, "HLT_PFJetFwd400"};
    TTreeReaderValue<Bool_t> HLT_PFJetFwd450 = {fReader, "HLT_PFJetFwd450"};
    TTreeReaderValue<Bool_t> HLT_PFJetFwd500 = {fReader, "HLT_PFJetFwd500"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJetFwd15 = {fReader, "HLT_AK8PFJetFwd15"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJetFwd25 = {fReader, "HLT_AK8PFJetFwd25"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJetFwd40 = {fReader, "HLT_AK8PFJetFwd40"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJetFwd60 = {fReader, "HLT_AK8PFJetFwd60"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJetFwd80 = {fReader, "HLT_AK8PFJetFwd80"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJetFwd140 = {fReader, "HLT_AK8PFJetFwd140"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJetFwd200 = {fReader, "HLT_AK8PFJetFwd200"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJetFwd260 = {fReader, "HLT_AK8PFJetFwd260"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJetFwd320 = {fReader, "HLT_AK8PFJetFwd320"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJetFwd400 = {fReader, "HLT_AK8PFJetFwd400"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJetFwd450 = {fReader, "HLT_AK8PFJetFwd450"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJetFwd500 = {fReader, "HLT_AK8PFJetFwd500"};
    TTreeReaderValue<Bool_t> HLT_PFHT180 = {fReader, "HLT_PFHT180"};
    TTreeReaderValue<Bool_t> HLT_PFHT250 = {fReader, "HLT_PFHT250"};
    TTreeReaderValue<Bool_t> HLT_PFHT370 = {fReader, "HLT_PFHT370"};
    TTreeReaderValue<Bool_t> HLT_PFHT430 = {fReader, "HLT_PFHT430"};
    TTreeReaderValue<Bool_t> HLT_PFHT510 = {fReader, "HLT_PFHT510"};
    TTreeReaderValue<Bool_t> HLT_PFHT590 = {fReader, "HLT_PFHT590"};
    TTreeReaderValue<Bool_t> HLT_PFHT680 = {fReader, "HLT_PFHT680"};
    TTreeReaderValue<Bool_t> HLT_PFHT780 = {fReader, "HLT_PFHT780"};
    TTreeReaderValue<Bool_t> HLT_PFHT890 = {fReader, "HLT_PFHT890"};
    TTreeReaderValue<Bool_t> HLT_PFHT1050 = {fReader, "HLT_PFHT1050"};
    TTreeReaderValue<Bool_t> HLT_PFHT500_PFMET100_PFMHT100_IDTight = {fReader, "HLT_PFHT500_PFMET100_PFMHT100_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFHT500_PFMET110_PFMHT110_IDTight = {fReader, "HLT_PFHT500_PFMET110_PFMHT110_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFHT700_PFMET85_PFMHT85_IDTight = {fReader, "HLT_PFHT700_PFMET85_PFMHT85_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFHT700_PFMET95_PFMHT95_IDTight = {fReader, "HLT_PFHT700_PFMET95_PFMHT95_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFHT800_PFMET75_PFMHT75_IDTight = {fReader, "HLT_PFHT800_PFMET75_PFMHT75_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFHT800_PFMET85_PFMHT85_IDTight = {fReader, "HLT_PFHT800_PFMET85_PFMHT85_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFMET110_PFMHT110_IDTight = {fReader, "HLT_PFMET110_PFMHT110_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFMET120_PFMHT120_IDTight = {fReader, "HLT_PFMET120_PFMHT120_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFMET130_PFMHT130_IDTight = {fReader, "HLT_PFMET130_PFMHT130_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFMET140_PFMHT140_IDTight = {fReader, "HLT_PFMET140_PFMHT140_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1 = {fReader, "HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1"};
    TTreeReaderValue<Bool_t> HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1 = {fReader, "HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1"};
    TTreeReaderValue<Bool_t> HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1 = {fReader, "HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1"};
    TTreeReaderValue<Bool_t> HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1 = {fReader, "HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1"};
    TTreeReaderValue<Bool_t> HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1 = {fReader, "HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1"};
    TTreeReaderValue<Bool_t> HLT_PFMET120_PFMHT120_IDTight_PFHT60 = {fReader, "HLT_PFMET120_PFMHT120_IDTight_PFHT60"};
    TTreeReaderValue<Bool_t> HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = {fReader, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60"};
    TTreeReaderValue<Bool_t> HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60 = {fReader, "HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60"};
    TTreeReaderValue<Bool_t> HLT_PFMETTypeOne110_PFMHT110_IDTight = {fReader, "HLT_PFMETTypeOne110_PFMHT110_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFMETTypeOne120_PFMHT120_IDTight = {fReader, "HLT_PFMETTypeOne120_PFMHT120_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFMETTypeOne130_PFMHT130_IDTight = {fReader, "HLT_PFMETTypeOne130_PFMHT130_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFMETTypeOne140_PFMHT140_IDTight = {fReader, "HLT_PFMETTypeOne140_PFMHT140_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFMETNoMu110_PFMHTNoMu110_IDTight = {fReader, "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = {fReader, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFMETNoMu130_PFMHTNoMu130_IDTight = {fReader, "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight"};
    TTreeReaderValue<Bool_t> HLT_PFMETNoMu140_PFMHTNoMu140_IDTight = {fReader, "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight"};
    TTreeReaderValue<Bool_t> HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight = {fReader, "HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight"};
    TTreeReaderValue<Bool_t> HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight = {fReader, "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"};
    TTreeReaderValue<Bool_t> HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight = {fReader, "HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight"};
    TTreeReaderValue<Bool_t> HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight = {fReader, "HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight"};
    TTreeReaderValue<Bool_t> HLT_L1ETMHadSeeds = {fReader, "HLT_L1ETMHadSeeds"};
    TTreeReaderValue<Bool_t> HLT_CaloMHT90 = {fReader, "HLT_CaloMHT90"};
    TTreeReaderValue<Bool_t> HLT_CaloMET80_NotCleaned = {fReader, "HLT_CaloMET80_NotCleaned"};
    TTreeReaderValue<Bool_t> HLT_CaloMET90_NotCleaned = {fReader, "HLT_CaloMET90_NotCleaned"};
    TTreeReaderValue<Bool_t> HLT_CaloMET100_NotCleaned = {fReader, "HLT_CaloMET100_NotCleaned"};
    TTreeReaderValue<Bool_t> HLT_CaloMET110_NotCleaned = {fReader, "HLT_CaloMET110_NotCleaned"};
    TTreeReaderValue<Bool_t> HLT_CaloMET250_NotCleaned = {fReader, "HLT_CaloMET250_NotCleaned"};
    TTreeReaderValue<Bool_t> HLT_CaloMET70_HBHECleaned = {fReader, "HLT_CaloMET70_HBHECleaned"};
    TTreeReaderValue<Bool_t> HLT_CaloMET80_HBHECleaned = {fReader, "HLT_CaloMET80_HBHECleaned"};
    TTreeReaderValue<Bool_t> HLT_CaloMET90_HBHECleaned = {fReader, "HLT_CaloMET90_HBHECleaned"};
    TTreeReaderValue<Bool_t> HLT_CaloMET100_HBHECleaned = {fReader, "HLT_CaloMET100_HBHECleaned"};
    TTreeReaderValue<Bool_t> HLT_CaloMET250_HBHECleaned = {fReader, "HLT_CaloMET250_HBHECleaned"};
    TTreeReaderValue<Bool_t> HLT_CaloMET300_HBHECleaned = {fReader, "HLT_CaloMET300_HBHECleaned"};
    TTreeReaderValue<Bool_t> HLT_CaloMET350_HBHECleaned = {fReader, "HLT_CaloMET350_HBHECleaned"};
    TTreeReaderValue<Bool_t> HLT_PFMET200_NotCleaned = {fReader, "HLT_PFMET200_NotCleaned"};
    TTreeReaderValue<Bool_t> HLT_PFMET200_HBHECleaned = {fReader, "HLT_PFMET200_HBHECleaned"};
    TTreeReaderValue<Bool_t> HLT_PFMET250_HBHECleaned = {fReader, "HLT_PFMET250_HBHECleaned"};
    TTreeReaderValue<Bool_t> HLT_PFMET300_HBHECleaned = {fReader, "HLT_PFMET300_HBHECleaned"};
    TTreeReaderValue<Bool_t> HLT_PFMET200_HBHE_BeamHaloCleaned = {fReader, "HLT_PFMET200_HBHE_BeamHaloCleaned"};
    TTreeReaderValue<Bool_t> HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned = {fReader, "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned"};
    TTreeReaderValue<Bool_t> HLT_MET105_IsoTrk50 = {fReader, "HLT_MET105_IsoTrk50"};
    TTreeReaderValue<Bool_t> HLT_MET120_IsoTrk50 = {fReader, "HLT_MET120_IsoTrk50"};
    TTreeReaderValue<Bool_t> HLT_SingleJet30_Mu12_SinglePFJet40 = {fReader, "HLT_SingleJet30_Mu12_SinglePFJet40"};
    TTreeReaderValue<Bool_t> HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71 = {fReader, "HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71"};
    TTreeReaderValue<Bool_t> HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71 = {fReader, "HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71"};
    TTreeReaderValue<Bool_t> HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71 = {fReader, "HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71"};
    TTreeReaderValue<Bool_t> HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71 = {fReader, "HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71"};
    TTreeReaderValue<Bool_t> HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 = {fReader, "HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71"};
    TTreeReaderValue<Bool_t> HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 = {fReader, "HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71"};
    TTreeReaderValue<Bool_t> HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 = {fReader, "HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71"};
    TTreeReaderValue<Bool_t> HLT_DoublePFJets40_CaloBTagDeepCSV_p71 = {fReader, "HLT_DoublePFJets40_CaloBTagDeepCSV_p71"};
    TTreeReaderValue<Bool_t> HLT_DoublePFJets100_CaloBTagDeepCSV_p71 = {fReader, "HLT_DoublePFJets100_CaloBTagDeepCSV_p71"};
    TTreeReaderValue<Bool_t> HLT_DoublePFJets200_CaloBTagDeepCSV_p71 = {fReader, "HLT_DoublePFJets200_CaloBTagDeepCSV_p71"};
    TTreeReaderValue<Bool_t> HLT_DoublePFJets350_CaloBTagDeepCSV_p71 = {fReader, "HLT_DoublePFJets350_CaloBTagDeepCSV_p71"};
    TTreeReaderValue<Bool_t> HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 = {fReader, "HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71"};
    TTreeReaderValue<Bool_t> HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 = {fReader, "HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71"};
    TTreeReaderValue<Bool_t> HLT_Photon300_NoHE = {fReader, "HLT_Photon300_NoHE"};
    TTreeReaderValue<Bool_t> HLT_Mu8_TrkIsoVVL = {fReader, "HLT_Mu8_TrkIsoVVL"};
    TTreeReaderValue<Bool_t> HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ = {fReader, "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ"};
    TTreeReaderValue<Bool_t> HLT_Mu8_DiEle12_CaloIdL_TrackIdL = {fReader, "HLT_Mu8_DiEle12_CaloIdL_TrackIdL"};
    TTreeReaderValue<Bool_t> HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ = {fReader, "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ"};
    TTreeReaderValue<Bool_t> HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350 = {fReader, "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350"};
    TTreeReaderValue<Bool_t> HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ = {fReader, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"};
    TTreeReaderValue<Bool_t> HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30 = {fReader, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30"};
    TTreeReaderValue<Bool_t> HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30 = {fReader, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30"};
    TTreeReaderValue<Bool_t> HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5 = {fReader, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5"};
    TTreeReaderValue<Bool_t> HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5 = {fReader, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5"};
    TTreeReaderValue<Bool_t> HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL = {fReader, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"};
    TTreeReaderValue<Bool_t> HLT_Mu17_TrkIsoVVL = {fReader, "HLT_Mu17_TrkIsoVVL"};
    TTreeReaderValue<Bool_t> HLT_Mu19_TrkIsoVVL = {fReader, "HLT_Mu19_TrkIsoVVL"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK4DiJet20_Mu5 = {fReader, "HLT_BTagMu_AK4DiJet20_Mu5"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK4DiJet40_Mu5 = {fReader, "HLT_BTagMu_AK4DiJet40_Mu5"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK4DiJet70_Mu5 = {fReader, "HLT_BTagMu_AK4DiJet70_Mu5"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK4DiJet110_Mu5 = {fReader, "HLT_BTagMu_AK4DiJet110_Mu5"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK4DiJet170_Mu5 = {fReader, "HLT_BTagMu_AK4DiJet170_Mu5"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK4Jet300_Mu5 = {fReader, "HLT_BTagMu_AK4Jet300_Mu5"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK8DiJet170_Mu5 = {fReader, "HLT_BTagMu_AK8DiJet170_Mu5"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK8Jet170_DoubleMu5 = {fReader, "HLT_BTagMu_AK8Jet170_DoubleMu5"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK8Jet300_Mu5 = {fReader, "HLT_BTagMu_AK8Jet300_Mu5"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK4DiJet20_Mu5_noalgo = {fReader, "HLT_BTagMu_AK4DiJet20_Mu5_noalgo"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK4DiJet40_Mu5_noalgo = {fReader, "HLT_BTagMu_AK4DiJet40_Mu5_noalgo"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK4DiJet70_Mu5_noalgo = {fReader, "HLT_BTagMu_AK4DiJet70_Mu5_noalgo"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK4DiJet110_Mu5_noalgo = {fReader, "HLT_BTagMu_AK4DiJet110_Mu5_noalgo"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK4DiJet170_Mu5_noalgo = {fReader, "HLT_BTagMu_AK4DiJet170_Mu5_noalgo"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK4Jet300_Mu5_noalgo = {fReader, "HLT_BTagMu_AK4Jet300_Mu5_noalgo"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK8DiJet170_Mu5_noalgo = {fReader, "HLT_BTagMu_AK8DiJet170_Mu5_noalgo"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo = {fReader, "HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo"};
    TTreeReaderValue<Bool_t> HLT_BTagMu_AK8Jet300_Mu5_noalgo = {fReader, "HLT_BTagMu_AK8Jet300_Mu5_noalgo"};
    TTreeReaderValue<Bool_t> HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL = {fReader, "HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL"};
    TTreeReaderValue<Bool_t> HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = {fReader, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"};
    TTreeReaderValue<Bool_t> HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL = {fReader, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"};
    TTreeReaderValue<Bool_t> HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = {fReader, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"};
    TTreeReaderValue<Bool_t> HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL = {fReader, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"};
    TTreeReaderValue<Bool_t> HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL = {fReader, "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"};
    TTreeReaderValue<Bool_t> HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ = {fReader, "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"};
    TTreeReaderValue<Bool_t> HLT_Mu12_DoublePhoton20 = {fReader, "HLT_Mu12_DoublePhoton20"};
    TTreeReaderValue<Bool_t> HLT_TriplePhoton_20_20_20_CaloIdLV2 = {fReader, "HLT_TriplePhoton_20_20_20_CaloIdLV2"};
    TTreeReaderValue<Bool_t> HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL = {fReader, "HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL"};
    TTreeReaderValue<Bool_t> HLT_TriplePhoton_30_30_10_CaloIdLV2 = {fReader, "HLT_TriplePhoton_30_30_10_CaloIdLV2"};
    TTreeReaderValue<Bool_t> HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL = {fReader, "HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL"};
    TTreeReaderValue<Bool_t> HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL = {fReader, "HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL"};
    TTreeReaderValue<Bool_t> HLT_Photon20 = {fReader, "HLT_Photon20"};
    TTreeReaderValue<Bool_t> HLT_Photon33 = {fReader, "HLT_Photon33"};
    TTreeReaderValue<Bool_t> HLT_Photon50 = {fReader, "HLT_Photon50"};
    TTreeReaderValue<Bool_t> HLT_Photon75 = {fReader, "HLT_Photon75"};
    TTreeReaderValue<Bool_t> HLT_Photon90 = {fReader, "HLT_Photon90"};
    TTreeReaderValue<Bool_t> HLT_Photon120 = {fReader, "HLT_Photon120"};
    TTreeReaderValue<Bool_t> HLT_Photon150 = {fReader, "HLT_Photon150"};
    TTreeReaderValue<Bool_t> HLT_Photon175 = {fReader, "HLT_Photon175"};
    TTreeReaderValue<Bool_t> HLT_Photon200 = {fReader, "HLT_Photon200"};
    TTreeReaderValue<Bool_t> HLT_Photon100EB_TightID_TightIso = {fReader, "HLT_Photon100EB_TightID_TightIso"};
    TTreeReaderValue<Bool_t> HLT_Photon110EB_TightID_TightIso = {fReader, "HLT_Photon110EB_TightID_TightIso"};
    TTreeReaderValue<Bool_t> HLT_Photon120EB_TightID_TightIso = {fReader, "HLT_Photon120EB_TightID_TightIso"};
    TTreeReaderValue<Bool_t> HLT_Photon100EBHE10 = {fReader, "HLT_Photon100EBHE10"};
    TTreeReaderValue<Bool_t> HLT_Photon100EEHE10 = {fReader, "HLT_Photon100EEHE10"};
    TTreeReaderValue<Bool_t> HLT_Photon100EE_TightID_TightIso = {fReader, "HLT_Photon100EE_TightID_TightIso"};
    TTreeReaderValue<Bool_t> HLT_Photon50_R9Id90_HE10_IsoM = {fReader, "HLT_Photon50_R9Id90_HE10_IsoM"};
    TTreeReaderValue<Bool_t> HLT_Photon75_R9Id90_HE10_IsoM = {fReader, "HLT_Photon75_R9Id90_HE10_IsoM"};
    TTreeReaderValue<Bool_t> HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3 = {fReader, "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3"};
    TTreeReaderValue<Bool_t> HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3 = {fReader, "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3"};
    TTreeReaderValue<Bool_t> HLT_Photon90_R9Id90_HE10_IsoM = {fReader, "HLT_Photon90_R9Id90_HE10_IsoM"};
    TTreeReaderValue<Bool_t> HLT_Photon120_R9Id90_HE10_IsoM = {fReader, "HLT_Photon120_R9Id90_HE10_IsoM"};
    TTreeReaderValue<Bool_t> HLT_Photon165_R9Id90_HE10_IsoM = {fReader, "HLT_Photon165_R9Id90_HE10_IsoM"};
    TTreeReaderValue<Bool_t> HLT_Photon90_CaloIdL_PFHT700 = {fReader, "HLT_Photon90_CaloIdL_PFHT700"};
    TTreeReaderValue<Bool_t> HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 = {fReader, "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90"};
    TTreeReaderValue<Bool_t> HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95 = {fReader, "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95"};
    TTreeReaderValue<Bool_t> HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55 = {fReader, "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55"};
    TTreeReaderValue<Bool_t> HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55 = {fReader, "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55"};
    TTreeReaderValue<Bool_t> HLT_Photon35_TwoProngs35 = {fReader, "HLT_Photon35_TwoProngs35"};
    TTreeReaderValue<Bool_t> HLT_IsoMu24_TwoProngs35 = {fReader, "HLT_IsoMu24_TwoProngs35"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_Jpsi_L1_NoOS = {fReader, "HLT_Dimuon0_Jpsi_L1_NoOS"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_Jpsi_NoVertexing_NoOS = {fReader, "HLT_Dimuon0_Jpsi_NoVertexing_NoOS"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_Jpsi = {fReader, "HLT_Dimuon0_Jpsi"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_Jpsi_NoVertexing = {fReader, "HLT_Dimuon0_Jpsi_NoVertexing"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_Jpsi_L1_4R_0er1p5R = {fReader, "HLT_Dimuon0_Jpsi_L1_4R_0er1p5R"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R = {fReader, "HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_Jpsi3p5_Muon2 = {fReader, "HLT_Dimuon0_Jpsi3p5_Muon2"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_Upsilon_L1_4p5 = {fReader, "HLT_Dimuon0_Upsilon_L1_4p5"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_Upsilon_L1_5 = {fReader, "HLT_Dimuon0_Upsilon_L1_5"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_Upsilon_L1_4p5NoOS = {fReader, "HLT_Dimuon0_Upsilon_L1_4p5NoOS"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_Upsilon_L1_4p5er2p0 = {fReader, "HLT_Dimuon0_Upsilon_L1_4p5er2p0"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_Upsilon_L1_4p5er2p0M = {fReader, "HLT_Dimuon0_Upsilon_L1_4p5er2p0M"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_Upsilon_NoVertexing = {fReader, "HLT_Dimuon0_Upsilon_NoVertexing"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_Upsilon_L1_5M = {fReader, "HLT_Dimuon0_Upsilon_L1_5M"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_LowMass_L1_0er1p5R = {fReader, "HLT_Dimuon0_LowMass_L1_0er1p5R"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_LowMass_L1_0er1p5 = {fReader, "HLT_Dimuon0_LowMass_L1_0er1p5"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_LowMass = {fReader, "HLT_Dimuon0_LowMass"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_LowMass_L1_4 = {fReader, "HLT_Dimuon0_LowMass_L1_4"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_LowMass_L1_4R = {fReader, "HLT_Dimuon0_LowMass_L1_4R"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_LowMass_L1_TM530 = {fReader, "HLT_Dimuon0_LowMass_L1_TM530"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_Upsilon_Muon_L1_TM0 = {fReader, "HLT_Dimuon0_Upsilon_Muon_L1_TM0"};
    TTreeReaderValue<Bool_t> HLT_Dimuon0_Upsilon_Muon_NoL1Mass = {fReader, "HLT_Dimuon0_Upsilon_Muon_NoL1Mass"};
    TTreeReaderValue<Bool_t> HLT_TripleMu_5_3_3_Mass3p8_DZ = {fReader, "HLT_TripleMu_5_3_3_Mass3p8_DZ"};
    TTreeReaderValue<Bool_t> HLT_TripleMu_10_5_5_DZ = {fReader, "HLT_TripleMu_10_5_5_DZ"};
    TTreeReaderValue<Bool_t> HLT_TripleMu_12_10_5 = {fReader, "HLT_TripleMu_12_10_5"};
    TTreeReaderValue<Bool_t> HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15 = {fReader, "HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15"};
    TTreeReaderValue<Bool_t> HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1 = {fReader, "HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1"};
    TTreeReaderValue<Bool_t> HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15 = {fReader, "HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15"};
    TTreeReaderValue<Bool_t> HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1 = {fReader, "HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu3_DZ_PFMET50_PFMHT60 = {fReader, "HLT_DoubleMu3_DZ_PFMET50_PFMHT60"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu3_DZ_PFMET70_PFMHT70 = {fReader, "HLT_DoubleMu3_DZ_PFMET70_PFMHT70"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu3_DZ_PFMET90_PFMHT90 = {fReader, "HLT_DoubleMu3_DZ_PFMET90_PFMHT90"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass = {fReader, "HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu4_Jpsi_Displaced = {fReader, "HLT_DoubleMu4_Jpsi_Displaced"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu4_Jpsi_NoVertexing = {fReader, "HLT_DoubleMu4_Jpsi_NoVertexing"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu4_JpsiTrkTrk_Displaced = {fReader, "HLT_DoubleMu4_JpsiTrkTrk_Displaced"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu43NoFiltersNoVtx = {fReader, "HLT_DoubleMu43NoFiltersNoVtx"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu48NoFiltersNoVtx = {fReader, "HLT_DoubleMu48NoFiltersNoVtx"};
    TTreeReaderValue<Bool_t> HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL = {fReader, "HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL"};
    TTreeReaderValue<Bool_t> HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL = {fReader, "HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL"};
    TTreeReaderValue<Bool_t> HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL = {fReader, "HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL"};
    TTreeReaderValue<Bool_t> HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL = {fReader, "HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu33NoFiltersNoVtxDisplaced = {fReader, "HLT_DoubleMu33NoFiltersNoVtxDisplaced"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu40NoFiltersNoVtxDisplaced = {fReader, "HLT_DoubleMu40NoFiltersNoVtxDisplaced"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu20_7_Mass0to30_L1_DM4 = {fReader, "HLT_DoubleMu20_7_Mass0to30_L1_DM4"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu20_7_Mass0to30_L1_DM4EG = {fReader, "HLT_DoubleMu20_7_Mass0to30_L1_DM4EG"};
    TTreeReaderValue<Bool_t> HLT_HT425 = {fReader, "HLT_HT425"};
    TTreeReaderValue<Bool_t> HLT_HT430_DisplacedDijet40_DisplacedTrack = {fReader, "HLT_HT430_DisplacedDijet40_DisplacedTrack"};
    TTreeReaderValue<Bool_t> HLT_HT500_DisplacedDijet40_DisplacedTrack = {fReader, "HLT_HT500_DisplacedDijet40_DisplacedTrack"};
    TTreeReaderValue<Bool_t> HLT_HT430_DisplacedDijet60_DisplacedTrack = {fReader, "HLT_HT430_DisplacedDijet60_DisplacedTrack"};
    TTreeReaderValue<Bool_t> HLT_HT400_DisplacedDijet40_DisplacedTrack = {fReader, "HLT_HT400_DisplacedDijet40_DisplacedTrack"};
    TTreeReaderValue<Bool_t> HLT_HT650_DisplacedDijet60_Inclusive = {fReader, "HLT_HT650_DisplacedDijet60_Inclusive"};
    TTreeReaderValue<Bool_t> HLT_HT550_DisplacedDijet60_Inclusive = {fReader, "HLT_HT550_DisplacedDijet60_Inclusive"};
    TTreeReaderValue<Bool_t> HLT_DiJet110_35_Mjj650_PFMET110 = {fReader, "HLT_DiJet110_35_Mjj650_PFMET110"};
    TTreeReaderValue<Bool_t> HLT_DiJet110_35_Mjj650_PFMET120 = {fReader, "HLT_DiJet110_35_Mjj650_PFMET120"};
    TTreeReaderValue<Bool_t> HLT_DiJet110_35_Mjj650_PFMET130 = {fReader, "HLT_DiJet110_35_Mjj650_PFMET130"};
    TTreeReaderValue<Bool_t> HLT_TripleJet110_35_35_Mjj650_PFMET110 = {fReader, "HLT_TripleJet110_35_35_Mjj650_PFMET110"};
    TTreeReaderValue<Bool_t> HLT_TripleJet110_35_35_Mjj650_PFMET120 = {fReader, "HLT_TripleJet110_35_35_Mjj650_PFMET120"};
    TTreeReaderValue<Bool_t> HLT_TripleJet110_35_35_Mjj650_PFMET130 = {fReader, "HLT_TripleJet110_35_35_Mjj650_PFMET130"};
    TTreeReaderValue<Bool_t> HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned = {fReader, "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned"};
    TTreeReaderValue<Bool_t> HLT_Ele28_eta2p1_WPTight_Gsf_HT150 = {fReader, "HLT_Ele28_eta2p1_WPTight_Gsf_HT150"};
    TTreeReaderValue<Bool_t> HLT_Ele28_HighEta_SC20_Mass55 = {fReader, "HLT_Ele28_HighEta_SC20_Mass55"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu20_7_Mass0to30_Photon23 = {fReader, "HLT_DoubleMu20_7_Mass0to30_Photon23"};
    TTreeReaderValue<Bool_t> HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5 = {fReader, "HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5"};
    TTreeReaderValue<Bool_t> HLT_Ele15_IsoVVVL_PFHT450_PFMET50 = {fReader, "HLT_Ele15_IsoVVVL_PFHT450_PFMET50"};
    TTreeReaderValue<Bool_t> HLT_Ele15_IsoVVVL_PFHT450 = {fReader, "HLT_Ele15_IsoVVVL_PFHT450"};
    TTreeReaderValue<Bool_t> HLT_Ele50_IsoVVVL_PFHT450 = {fReader, "HLT_Ele50_IsoVVVL_PFHT450"};
    TTreeReaderValue<Bool_t> HLT_Ele15_IsoVVVL_PFHT600 = {fReader, "HLT_Ele15_IsoVVVL_PFHT600"};
    TTreeReaderValue<Bool_t> HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60 = {fReader, "HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60"};
    TTreeReaderValue<Bool_t> HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60 = {fReader, "HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60"};
    TTreeReaderValue<Bool_t> HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60 = {fReader, "HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60"};
    TTreeReaderValue<Bool_t> HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5 = {fReader, "HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5"};
    TTreeReaderValue<Bool_t> HLT_Mu15_IsoVVVL_PFHT450_PFMET50 = {fReader, "HLT_Mu15_IsoVVVL_PFHT450_PFMET50"};
    TTreeReaderValue<Bool_t> HLT_Mu15_IsoVVVL_PFHT450 = {fReader, "HLT_Mu15_IsoVVVL_PFHT450"};
    TTreeReaderValue<Bool_t> HLT_Mu50_IsoVVVL_PFHT450 = {fReader, "HLT_Mu50_IsoVVVL_PFHT450"};
    TTreeReaderValue<Bool_t> HLT_Mu15_IsoVVVL_PFHT600 = {fReader, "HLT_Mu15_IsoVVVL_PFHT600"};
    TTreeReaderValue<Bool_t> HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight = {fReader, "HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight"};
    TTreeReaderValue<Bool_t> HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight = {fReader, "HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight"};
    TTreeReaderValue<Bool_t> HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight = {fReader, "HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight"};
    TTreeReaderValue<Bool_t> HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight = {fReader, "HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight"};
    TTreeReaderValue<Bool_t> HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight = {fReader, "HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight"};
    TTreeReaderValue<Bool_t> HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight = {fReader, "HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight"};
    TTreeReaderValue<Bool_t> HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight = {fReader, "HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight"};
    TTreeReaderValue<Bool_t> HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight = {fReader, "HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight"};
    TTreeReaderValue<Bool_t> HLT_Dimuon10_PsiPrime_Barrel_Seagulls = {fReader, "HLT_Dimuon10_PsiPrime_Barrel_Seagulls"};
    TTreeReaderValue<Bool_t> HLT_Dimuon20_Jpsi_Barrel_Seagulls = {fReader, "HLT_Dimuon20_Jpsi_Barrel_Seagulls"};
    TTreeReaderValue<Bool_t> HLT_Dimuon12_Upsilon_y1p4 = {fReader, "HLT_Dimuon12_Upsilon_y1p4"};
    TTreeReaderValue<Bool_t> HLT_Dimuon14_Phi_Barrel_Seagulls = {fReader, "HLT_Dimuon14_Phi_Barrel_Seagulls"};
    TTreeReaderValue<Bool_t> HLT_Dimuon18_PsiPrime = {fReader, "HLT_Dimuon18_PsiPrime"};
    TTreeReaderValue<Bool_t> HLT_Dimuon25_Jpsi = {fReader, "HLT_Dimuon25_Jpsi"};
    TTreeReaderValue<Bool_t> HLT_Dimuon18_PsiPrime_noCorrL1 = {fReader, "HLT_Dimuon18_PsiPrime_noCorrL1"};
    TTreeReaderValue<Bool_t> HLT_Dimuon24_Upsilon_noCorrL1 = {fReader, "HLT_Dimuon24_Upsilon_noCorrL1"};
    TTreeReaderValue<Bool_t> HLT_Dimuon24_Phi_noCorrL1 = {fReader, "HLT_Dimuon24_Phi_noCorrL1"};
    TTreeReaderValue<Bool_t> HLT_Dimuon25_Jpsi_noCorrL1 = {fReader, "HLT_Dimuon25_Jpsi_noCorrL1"};
    TTreeReaderValue<Bool_t> HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8 = {fReader, "HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8"};
    TTreeReaderValue<Bool_t> HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ = {fReader, "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ"};
    TTreeReaderValue<Bool_t> HLT_DiMu9_Ele9_CaloIdL_TrackIdL = {fReader, "HLT_DiMu9_Ele9_CaloIdL_TrackIdL"};
    TTreeReaderValue<Bool_t> HLT_DoubleIsoMu20_eta2p1 = {fReader, "HLT_DoubleIsoMu20_eta2p1"};
    TTreeReaderValue<Bool_t> HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx = {fReader, "HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx"};
    TTreeReaderValue<Bool_t> HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx = {fReader, "HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx"};
    TTreeReaderValue<Bool_t> HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx = {fReader, "HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx"};
    TTreeReaderValue<Bool_t> HLT_Mu8 = {fReader, "HLT_Mu8"};
    TTreeReaderValue<Bool_t> HLT_Mu17 = {fReader, "HLT_Mu17"};
    TTreeReaderValue<Bool_t> HLT_Mu19 = {fReader, "HLT_Mu19"};
    TTreeReaderValue<Bool_t> HLT_Mu17_Photon30_IsoCaloId = {fReader, "HLT_Mu17_Photon30_IsoCaloId"};
    TTreeReaderValue<Bool_t> HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30 = {fReader, "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30"};
    TTreeReaderValue<Bool_t> HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30 = {fReader, "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30"};
    TTreeReaderValue<Bool_t> HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30 = {fReader, "HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30"};
    TTreeReaderValue<Bool_t> HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30 = {fReader, "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30"};
    TTreeReaderValue<Bool_t> HLT_Ele8_CaloIdM_TrackIdM_PFJet30 = {fReader, "HLT_Ele8_CaloIdM_TrackIdM_PFJet30"};
    TTreeReaderValue<Bool_t> HLT_Ele17_CaloIdM_TrackIdM_PFJet30 = {fReader, "HLT_Ele17_CaloIdM_TrackIdM_PFJet30"};
    TTreeReaderValue<Bool_t> HLT_Ele23_CaloIdM_TrackIdM_PFJet30 = {fReader, "HLT_Ele23_CaloIdM_TrackIdM_PFJet30"};
    TTreeReaderValue<Bool_t> HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 = {fReader, "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165"};
    TTreeReaderValue<Bool_t> HLT_Ele115_CaloIdVT_GsfTrkIdT = {fReader, "HLT_Ele115_CaloIdVT_GsfTrkIdT"};
    TTreeReaderValue<Bool_t> HLT_Ele135_CaloIdVT_GsfTrkIdT = {fReader, "HLT_Ele135_CaloIdVT_GsfTrkIdT"};
    TTreeReaderValue<Bool_t> HLT_Ele145_CaloIdVT_GsfTrkIdT = {fReader, "HLT_Ele145_CaloIdVT_GsfTrkIdT"};
    TTreeReaderValue<Bool_t> HLT_Ele200_CaloIdVT_GsfTrkIdT = {fReader, "HLT_Ele200_CaloIdVT_GsfTrkIdT"};
    TTreeReaderValue<Bool_t> HLT_Ele250_CaloIdVT_GsfTrkIdT = {fReader, "HLT_Ele250_CaloIdVT_GsfTrkIdT"};
    TTreeReaderValue<Bool_t> HLT_Ele300_CaloIdVT_GsfTrkIdT = {fReader, "HLT_Ele300_CaloIdVT_GsfTrkIdT"};
    TTreeReaderValue<Bool_t> HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5 = {fReader, "HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5"};
    TTreeReaderValue<Bool_t> HLT_PFHT330PT30_QuadPFJet_75_60_45_40 = {fReader, "HLT_PFHT330PT30_QuadPFJet_75_60_45_40"};
    TTreeReaderValue<Bool_t> HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94 = {fReader, "HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94"};
    TTreeReaderValue<Bool_t> HLT_PFHT400_SixPFJet32 = {fReader, "HLT_PFHT400_SixPFJet32"};
    TTreeReaderValue<Bool_t> HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59 = {fReader, "HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59"};
    TTreeReaderValue<Bool_t> HLT_PFHT450_SixPFJet36 = {fReader, "HLT_PFHT450_SixPFJet36"};
    TTreeReaderValue<Bool_t> HLT_PFHT350 = {fReader, "HLT_PFHT350"};
    TTreeReaderValue<Bool_t> HLT_PFHT350MinPFJet15 = {fReader, "HLT_PFHT350MinPFJet15"};
    TTreeReaderValue<Bool_t> HLT_Photon60_R9Id90_CaloIdL_IsoL = {fReader, "HLT_Photon60_R9Id90_CaloIdL_IsoL"};
    TTreeReaderValue<Bool_t> HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL = {fReader, "HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL"};
    TTreeReaderValue<Bool_t> HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15 = {fReader, "HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15"};
    TTreeReaderValue<Bool_t> HLT_ECALHT800 = {fReader, "HLT_ECALHT800"};
    TTreeReaderValue<Bool_t> HLT_DiSC30_18_EIso_AND_HE_Mass70 = {fReader, "HLT_DiSC30_18_EIso_AND_HE_Mass70"};
    TTreeReaderValue<Bool_t> HLT_Physics = {fReader, "HLT_Physics"};
    TTreeReaderValue<Bool_t> HLT_Physics_part0 = {fReader, "HLT_Physics_part0"};
    TTreeReaderValue<Bool_t> HLT_Physics_part1 = {fReader, "HLT_Physics_part1"};
    TTreeReaderValue<Bool_t> HLT_Physics_part2 = {fReader, "HLT_Physics_part2"};
    TTreeReaderValue<Bool_t> HLT_Physics_part3 = {fReader, "HLT_Physics_part3"};
    TTreeReaderValue<Bool_t> HLT_Physics_part4 = {fReader, "HLT_Physics_part4"};
    TTreeReaderValue<Bool_t> HLT_Physics_part5 = {fReader, "HLT_Physics_part5"};
    TTreeReaderValue<Bool_t> HLT_Physics_part6 = {fReader, "HLT_Physics_part6"};
    TTreeReaderValue<Bool_t> HLT_Physics_part7 = {fReader, "HLT_Physics_part7"};
    TTreeReaderValue<Bool_t> HLT_Random = {fReader, "HLT_Random"};
    TTreeReaderValue<Bool_t> HLT_ZeroBias = {fReader, "HLT_ZeroBias"};
    TTreeReaderValue<Bool_t> HLT_ZeroBias_Alignment = {fReader, "HLT_ZeroBias_Alignment"};
    TTreeReaderValue<Bool_t> HLT_ZeroBias_part0 = {fReader, "HLT_ZeroBias_part0"};
    TTreeReaderValue<Bool_t> HLT_ZeroBias_part1 = {fReader, "HLT_ZeroBias_part1"};
    TTreeReaderValue<Bool_t> HLT_ZeroBias_part2 = {fReader, "HLT_ZeroBias_part2"};
    TTreeReaderValue<Bool_t> HLT_ZeroBias_part3 = {fReader, "HLT_ZeroBias_part3"};
    TTreeReaderValue<Bool_t> HLT_ZeroBias_part4 = {fReader, "HLT_ZeroBias_part4"};
    TTreeReaderValue<Bool_t> HLT_ZeroBias_part5 = {fReader, "HLT_ZeroBias_part5"};
    TTreeReaderValue<Bool_t> HLT_ZeroBias_part6 = {fReader, "HLT_ZeroBias_part6"};
    TTreeReaderValue<Bool_t> HLT_ZeroBias_part7 = {fReader, "HLT_ZeroBias_part7"};
    TTreeReaderValue<Bool_t> HLT_AK4CaloJet30 = {fReader, "HLT_AK4CaloJet30"};
    TTreeReaderValue<Bool_t> HLT_AK4CaloJet40 = {fReader, "HLT_AK4CaloJet40"};
    TTreeReaderValue<Bool_t> HLT_AK4CaloJet50 = {fReader, "HLT_AK4CaloJet50"};
    TTreeReaderValue<Bool_t> HLT_AK4CaloJet80 = {fReader, "HLT_AK4CaloJet80"};
    TTreeReaderValue<Bool_t> HLT_AK4CaloJet100 = {fReader, "HLT_AK4CaloJet100"};
    TTreeReaderValue<Bool_t> HLT_AK4CaloJet120 = {fReader, "HLT_AK4CaloJet120"};
    TTreeReaderValue<Bool_t> HLT_AK4PFJet30 = {fReader, "HLT_AK4PFJet30"};
    TTreeReaderValue<Bool_t> HLT_AK4PFJet50 = {fReader, "HLT_AK4PFJet50"};
    TTreeReaderValue<Bool_t> HLT_AK4PFJet80 = {fReader, "HLT_AK4PFJet80"};
    TTreeReaderValue<Bool_t> HLT_AK4PFJet100 = {fReader, "HLT_AK4PFJet100"};
    TTreeReaderValue<Bool_t> HLT_AK4PFJet120 = {fReader, "HLT_AK4PFJet120"};
    TTreeReaderValue<Bool_t> HLT_SinglePhoton10_Eta3p1ForPPRef = {fReader, "HLT_SinglePhoton10_Eta3p1ForPPRef"};
    TTreeReaderValue<Bool_t> HLT_SinglePhoton20_Eta3p1ForPPRef = {fReader, "HLT_SinglePhoton20_Eta3p1ForPPRef"};
    TTreeReaderValue<Bool_t> HLT_SinglePhoton30_Eta3p1ForPPRef = {fReader, "HLT_SinglePhoton30_Eta3p1ForPPRef"};
    TTreeReaderValue<Bool_t> HLT_Photon20_HoverELoose = {fReader, "HLT_Photon20_HoverELoose"};
    TTreeReaderValue<Bool_t> HLT_Photon30_HoverELoose = {fReader, "HLT_Photon30_HoverELoose"};
    TTreeReaderValue<Bool_t> HLT_EcalCalibration = {fReader, "HLT_EcalCalibration"};
    TTreeReaderValue<Bool_t> HLT_HcalCalibration = {fReader, "HLT_HcalCalibration"};
    TTreeReaderValue<Bool_t> HLT_L1UnpairedBunchBptxMinus = {fReader, "HLT_L1UnpairedBunchBptxMinus"};
    TTreeReaderValue<Bool_t> HLT_L1UnpairedBunchBptxPlus = {fReader, "HLT_L1UnpairedBunchBptxPlus"};
    TTreeReaderValue<Bool_t> HLT_L1NotBptxOR = {fReader, "HLT_L1NotBptxOR"};
    TTreeReaderValue<Bool_t> HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142 = {fReader, "HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142"};
    TTreeReaderValue<Bool_t> HLT_CDC_L2cosmic_5_er1p0 = {fReader, "HLT_CDC_L2cosmic_5_er1p0"};
    TTreeReaderValue<Bool_t> HLT_CDC_L2cosmic_5p5_er1p0 = {fReader, "HLT_CDC_L2cosmic_5p5_er1p0"};
    TTreeReaderValue<Bool_t> HLT_HcalNZS = {fReader, "HLT_HcalNZS"};
    TTreeReaderValue<Bool_t> HLT_HcalPhiSym = {fReader, "HLT_HcalPhiSym"};
    TTreeReaderValue<Bool_t> HLT_HcalIsolatedbunch = {fReader, "HLT_HcalIsolatedbunch"};
    TTreeReaderValue<Bool_t> HLT_IsoTrackHB = {fReader, "HLT_IsoTrackHB"};
    TTreeReaderValue<Bool_t> HLT_IsoTrackHE = {fReader, "HLT_IsoTrackHE"};
    TTreeReaderValue<Bool_t> HLT_ZeroBias_FirstCollisionAfterAbortGap = {fReader, "HLT_ZeroBias_FirstCollisionAfterAbortGap"};
    TTreeReaderValue<Bool_t> HLT_ZeroBias_IsolatedBunches = {fReader, "HLT_ZeroBias_IsolatedBunches"};
    TTreeReaderValue<Bool_t> HLT_ZeroBias_FirstCollisionInTrain = {fReader, "HLT_ZeroBias_FirstCollisionInTrain"};
    TTreeReaderValue<Bool_t> HLT_ZeroBias_LastCollisionInTrain = {fReader, "HLT_ZeroBias_LastCollisionInTrain"};
    TTreeReaderValue<Bool_t> HLT_ZeroBias_FirstBXAfterTrain = {fReader, "HLT_ZeroBias_FirstBXAfterTrain"};
    TTreeReaderValue<Bool_t> HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr = {fReader, "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr"};
    TTreeReaderValue<Bool_t> HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90 = {fReader, "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90"};
    TTreeReaderValue<Bool_t> HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100 = {fReader, "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100"};
    TTreeReaderValue<Bool_t> HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110 = {fReader, "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110"};
    TTreeReaderValue<Bool_t> HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120 = {fReader, "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120"};
    TTreeReaderValue<Bool_t> HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130 = {fReader, "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130"};
    TTreeReaderValue<Bool_t> HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140 = {fReader, "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140"};
    TTreeReaderValue<Bool_t> HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr = {fReader, "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr"};
    TTreeReaderValue<Bool_t> HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr = {fReader, "HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr"};
    TTreeReaderValue<Bool_t> HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1 = {fReader, "HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1"};
    TTreeReaderValue<Bool_t> HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1 = {fReader, "HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1"};
    TTreeReaderValue<Bool_t> HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1 = {fReader, "HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1"};
    TTreeReaderValue<Bool_t> HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL = {fReader, "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL"};
    TTreeReaderValue<Bool_t> HLT_Rsq0p35 = {fReader, "HLT_Rsq0p35"};
    TTreeReaderValue<Bool_t> HLT_Rsq0p40 = {fReader, "HLT_Rsq0p40"};
    TTreeReaderValue<Bool_t> HLT_RsqMR300_Rsq0p09_MR200 = {fReader, "HLT_RsqMR300_Rsq0p09_MR200"};
    TTreeReaderValue<Bool_t> HLT_RsqMR320_Rsq0p09_MR200 = {fReader, "HLT_RsqMR320_Rsq0p09_MR200"};
    TTreeReaderValue<Bool_t> HLT_RsqMR300_Rsq0p09_MR200_4jet = {fReader, "HLT_RsqMR300_Rsq0p09_MR200_4jet"};
    TTreeReaderValue<Bool_t> HLT_RsqMR320_Rsq0p09_MR200_4jet = {fReader, "HLT_RsqMR320_Rsq0p09_MR200_4jet"};
    TTreeReaderValue<Bool_t> HLT_IsoMu27_MET90 = {fReader, "HLT_IsoMu27_MET90"};
    TTreeReaderValue<Bool_t> HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg = {fReader, "HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg"};
    TTreeReaderValue<Bool_t> HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg = {fReader, "HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg"};
    TTreeReaderValue<Bool_t> HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg = {fReader, "HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg"};
    TTreeReaderValue<Bool_t> HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg = {fReader, "HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg"};
    TTreeReaderValue<Bool_t> HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg = {fReader, "HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg"};
    TTreeReaderValue<Bool_t> HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg = {fReader, "HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg"};
    TTreeReaderValue<Bool_t> HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg = {fReader, "HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg"};
    TTreeReaderValue<Bool_t> HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg = {fReader, "HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg"};
    TTreeReaderValue<Bool_t> HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1 = {fReader, "HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1"};
    TTreeReaderValue<Bool_t> HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1 = {fReader, "HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1"};
    TTreeReaderValue<Bool_t> HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1 = {fReader, "HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1"};
    TTreeReaderValue<Bool_t> HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50 = {fReader, "HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50"};
    TTreeReaderValue<Bool_t> HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3 = {fReader, "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3"};
    TTreeReaderValue<Bool_t> HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3 = {fReader, "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3"};
    TTreeReaderValue<Bool_t> HLT_PFMET100_PFMHT100_IDTight_PFHT60 = {fReader, "HLT_PFMET100_PFMHT100_IDTight_PFHT60"};
    TTreeReaderValue<Bool_t> HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60 = {fReader, "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60"};
    TTreeReaderValue<Bool_t> HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60 = {fReader, "HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60"};
    TTreeReaderValue<Bool_t> HLT_Mu18_Mu9_SameSign = {fReader, "HLT_Mu18_Mu9_SameSign"};
    TTreeReaderValue<Bool_t> HLT_Mu18_Mu9_SameSign_DZ = {fReader, "HLT_Mu18_Mu9_SameSign_DZ"};
    TTreeReaderValue<Bool_t> HLT_Mu18_Mu9 = {fReader, "HLT_Mu18_Mu9"};
    TTreeReaderValue<Bool_t> HLT_Mu18_Mu9_DZ = {fReader, "HLT_Mu18_Mu9_DZ"};
    TTreeReaderValue<Bool_t> HLT_Mu20_Mu10_SameSign = {fReader, "HLT_Mu20_Mu10_SameSign"};
    TTreeReaderValue<Bool_t> HLT_Mu20_Mu10_SameSign_DZ = {fReader, "HLT_Mu20_Mu10_SameSign_DZ"};
    TTreeReaderValue<Bool_t> HLT_Mu20_Mu10 = {fReader, "HLT_Mu20_Mu10"};
    TTreeReaderValue<Bool_t> HLT_Mu20_Mu10_DZ = {fReader, "HLT_Mu20_Mu10_DZ"};
    TTreeReaderValue<Bool_t> HLT_Mu23_Mu12_SameSign = {fReader, "HLT_Mu23_Mu12_SameSign"};
    TTreeReaderValue<Bool_t> HLT_Mu23_Mu12_SameSign_DZ = {fReader, "HLT_Mu23_Mu12_SameSign_DZ"};
    TTreeReaderValue<Bool_t> HLT_Mu23_Mu12 = {fReader, "HLT_Mu23_Mu12"};
    TTreeReaderValue<Bool_t> HLT_Mu23_Mu12_DZ = {fReader, "HLT_Mu23_Mu12_DZ"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05 = {fReader, "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi = {fReader, "HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi"};
    TTreeReaderValue<Bool_t> HLT_DoubleMu3_DCA_PFMET50_PFMHT60 = {fReader, "HLT_DoubleMu3_DCA_PFMET50_PFMHT60"};
    TTreeReaderValue<Bool_t> HLT_TripleMu_5_3_3_Mass3p8_DCA = {fReader, "HLT_TripleMu_5_3_3_Mass3p8_DCA"};
    TTreeReaderValue<Bool_t> HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1 = {fReader, "HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1"};
    TTreeReaderValue<Bool_t> HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1 = {fReader, "HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1"};
    TTreeReaderValue<Bool_t> HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1 = {fReader, "HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1"};
    TTreeReaderValue<Bool_t> HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2 = {fReader, "HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2"};
    TTreeReaderValue<Bool_t> HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2 = {fReader, "HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2"};
    TTreeReaderValue<Bool_t> HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2 = {fReader, "HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2"};
    TTreeReaderValue<Bool_t> HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2 = {fReader, "HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2"};
    TTreeReaderValue<Bool_t> HLT_QuadPFJet98_83_71_15 = {fReader, "HLT_QuadPFJet98_83_71_15"};
    TTreeReaderValue<Bool_t> HLT_QuadPFJet103_88_75_15 = {fReader, "HLT_QuadPFJet103_88_75_15"};
    TTreeReaderValue<Bool_t> HLT_QuadPFJet105_88_76_15 = {fReader, "HLT_QuadPFJet105_88_76_15"};
    TTreeReaderValue<Bool_t> HLT_QuadPFJet111_90_80_15 = {fReader, "HLT_QuadPFJet111_90_80_15"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17 = {fReader, "HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1 = {fReader, "HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02 = {fReader, "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2 = {fReader, "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2"};
    TTreeReaderValue<Bool_t> HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 = {fReader, "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4"};
    TTreeReaderValue<Bool_t> HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55 = {fReader, "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55"};
    TTreeReaderValue<Bool_t> HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto = {fReader, "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto"};
    TTreeReaderValue<Bool_t> HLT_Mu12_IP6_part0 = {fReader, "HLT_Mu12_IP6_part0"};
    TTreeReaderValue<Bool_t> HLT_Mu12_IP6_part1 = {fReader, "HLT_Mu12_IP6_part1"};
    TTreeReaderValue<Bool_t> HLT_Mu12_IP6_part2 = {fReader, "HLT_Mu12_IP6_part2"};
    TTreeReaderValue<Bool_t> HLT_Mu12_IP6_part3 = {fReader, "HLT_Mu12_IP6_part3"};
    TTreeReaderValue<Bool_t> HLT_Mu12_IP6_part4 = {fReader, "HLT_Mu12_IP6_part4"};
    TTreeReaderValue<Bool_t> HLT_Mu9_IP5_part0 = {fReader, "HLT_Mu9_IP5_part0"};
    TTreeReaderValue<Bool_t> HLT_Mu9_IP5_part1 = {fReader, "HLT_Mu9_IP5_part1"};
    TTreeReaderValue<Bool_t> HLT_Mu9_IP5_part2 = {fReader, "HLT_Mu9_IP5_part2"};
    TTreeReaderValue<Bool_t> HLT_Mu9_IP5_part3 = {fReader, "HLT_Mu9_IP5_part3"};
    TTreeReaderValue<Bool_t> HLT_Mu9_IP5_part4 = {fReader, "HLT_Mu9_IP5_part4"};
    TTreeReaderValue<Bool_t> HLT_Mu7_IP4_part0 = {fReader, "HLT_Mu7_IP4_part0"};
    TTreeReaderValue<Bool_t> HLT_Mu7_IP4_part1 = {fReader, "HLT_Mu7_IP4_part1"};
    TTreeReaderValue<Bool_t> HLT_Mu7_IP4_part2 = {fReader, "HLT_Mu7_IP4_part2"};
    TTreeReaderValue<Bool_t> HLT_Mu7_IP4_part3 = {fReader, "HLT_Mu7_IP4_part3"};
    TTreeReaderValue<Bool_t> HLT_Mu7_IP4_part4 = {fReader, "HLT_Mu7_IP4_part4"};
    TTreeReaderValue<Bool_t> HLT_Mu9_IP4_part0 = {fReader, "HLT_Mu9_IP4_part0"};
    TTreeReaderValue<Bool_t> HLT_Mu9_IP4_part1 = {fReader, "HLT_Mu9_IP4_part1"};
    TTreeReaderValue<Bool_t> HLT_Mu9_IP4_part2 = {fReader, "HLT_Mu9_IP4_part2"};
    TTreeReaderValue<Bool_t> HLT_Mu9_IP4_part3 = {fReader, "HLT_Mu9_IP4_part3"};
    TTreeReaderValue<Bool_t> HLT_Mu9_IP4_part4 = {fReader, "HLT_Mu9_IP4_part4"};
    TTreeReaderValue<Bool_t> HLT_Mu8_IP5_part0 = {fReader, "HLT_Mu8_IP5_part0"};
    TTreeReaderValue<Bool_t> HLT_Mu8_IP5_part1 = {fReader, "HLT_Mu8_IP5_part1"};
    TTreeReaderValue<Bool_t> HLT_Mu8_IP5_part2 = {fReader, "HLT_Mu8_IP5_part2"};
    TTreeReaderValue<Bool_t> HLT_Mu8_IP5_part3 = {fReader, "HLT_Mu8_IP5_part3"};
    TTreeReaderValue<Bool_t> HLT_Mu8_IP5_part4 = {fReader, "HLT_Mu8_IP5_part4"};
    TTreeReaderValue<Bool_t> HLT_Mu8_IP6_part0 = {fReader, "HLT_Mu8_IP6_part0"};
    TTreeReaderValue<Bool_t> HLT_Mu8_IP6_part1 = {fReader, "HLT_Mu8_IP6_part1"};
    TTreeReaderValue<Bool_t> HLT_Mu8_IP6_part2 = {fReader, "HLT_Mu8_IP6_part2"};
    TTreeReaderValue<Bool_t> HLT_Mu8_IP6_part3 = {fReader, "HLT_Mu8_IP6_part3"};
    TTreeReaderValue<Bool_t> HLT_Mu8_IP6_part4 = {fReader, "HLT_Mu8_IP6_part4"};
    TTreeReaderValue<Bool_t> HLT_Mu9_IP6_part0 = {fReader, "HLT_Mu9_IP6_part0"};
    TTreeReaderValue<Bool_t> HLT_Mu9_IP6_part1 = {fReader, "HLT_Mu9_IP6_part1"};
    TTreeReaderValue<Bool_t> HLT_Mu9_IP6_part2 = {fReader, "HLT_Mu9_IP6_part2"};
    TTreeReaderValue<Bool_t> HLT_Mu9_IP6_part3 = {fReader, "HLT_Mu9_IP6_part3"};
    TTreeReaderValue<Bool_t> HLT_Mu9_IP6_part4 = {fReader, "HLT_Mu9_IP6_part4"};
    TTreeReaderValue<Bool_t> HLT_Mu8_IP3_part0 = {fReader, "HLT_Mu8_IP3_part0"};
    TTreeReaderValue<Bool_t> HLT_Mu8_IP3_part1 = {fReader, "HLT_Mu8_IP3_part1"};
    TTreeReaderValue<Bool_t> HLT_Mu8_IP3_part2 = {fReader, "HLT_Mu8_IP3_part2"};
    TTreeReaderValue<Bool_t> HLT_Mu8_IP3_part3 = {fReader, "HLT_Mu8_IP3_part3"};
    TTreeReaderValue<Bool_t> HLT_Mu8_IP3_part4 = {fReader, "HLT_Mu8_IP3_part4"};
    TTreeReaderValue<Bool_t> HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1 = {fReader, "HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1"};
    TTreeReaderValue<Bool_t> HLT_TrkMu6NoFiltersNoVtx = {fReader, "HLT_TrkMu6NoFiltersNoVtx"};
    TTreeReaderValue<Bool_t> HLT_TrkMu16NoFiltersNoVtx = {fReader, "HLT_TrkMu16NoFiltersNoVtx"};
    TTreeReaderValue<Bool_t> HLT_DoubleTrkMu_16_6_NoFiltersNoVtx = {fReader, "HLT_DoubleTrkMu_16_6_NoFiltersNoVtx"};
    TTreeReaderValue<Bool_t> HLTriggerFinalPath = {fReader, "HLTriggerFinalPath"};
    TTreeReaderValue<Bool_t> L1simulation_step = {fReader, "L1simulation_step"};
  */

  //################################################################################################
  
  SingleTopAna(TTree * /*tree*/ =0) { }
  virtual ~SingleTopAna() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();


  //User defined functions are declared here
  
  void SetHstFileName(const char *HstFileName){ _HstFileName = HstFileName;}
  void SetSumFileName(const char *SumFileName){ _SumFileName = SumFileName;}
  void SetSample(int sample){_sample=sample;}
  void SetVerbose(int verbose){ _verbosity = verbose; }
  void SetData(int data){_data=data;}
  void SetYear(int year){_year = year;}
  void SetEra(TString era){_era=era;}
  void SetLep(int lep){_lep=lep;}
  void BookHistograms();
  float delta_phi(float phi1, float phi2);
  float transv_mass(float E_lep, float MET, float dphi);
  int GenMother(int ind, int mom_ind);
  int get_mother(int index);

  double Muon_2018UL_Reco(float pt, float eta);
  double Muon_2018UL_ID(float pt, float eta);
  double Muon_2018UL_Iso(float pt, float eta);
  double TrigEff_2018_IsoMu24_Data(float pt, float eta);
  
public:
  struct Hists {

    //###################################
    //Histograms are declared here.

    TH1F *nmuons_tight;
    TH1F *nele_loose;
    TH1F *cut_based;
    TH1F *cutbased_baseline;
    
    TH1F *n1mu2J;
    TH1F *n1e2J;
    TH1F *n1mu2J_cut;
    TH1F *evt_wt_cut;
    TH1F *evt_wt;
    TH1F *elesize;
    TH1F *lepsize;

    TH1F *nEvtRan;
    TH1F *nEvtTotal;
    
    //FOR LEPTONS
    TH1F *nLep;
    TH1F *leppt;
    TH1F *lepeta;
    TH1F *lepphi;
    TH1F *lep[30];

    //FOR MUONS
    TH1F *nmuons;
    TH1F *mupt;
    TH1F *mueta;
    TH1F *muphi;
    TH1F *mu_iso;
    TH1F *muprop[30];  
    TH1F *mu[30];

    //FOR ELECTRONS
    TH1F *nelectrons;
    TH1F *elept;
    TH1F *eleeta;
    TH1F *elephi;
    TH1F *ele[30];

    //FOR JETS
    TH1F *njets;
    TH1F *jetpt;
    TH1F *jeteta;
    TH1F *jetphi;
    TH1F *jet_score;
    TH1F *jet[30];
    //b jets
    TH1F *bjet_score_deepB;
    TH1F *n_bjets;
    // TH1F *bjet_score_CSV;
    // TH1F *n_bjets_CSV;
    TH1F *b_jetpt;
    TH1F *b_jeteta;
    TH1F *b_jetphi;
    TH1F *b_jet[30];

    // FOR MET
    TH1F *MET[30];
    TH1F *N[30];

    //FOR OBJECTS AFTER SELECTIONS (cut)
    TH1F *mu_cut[30];
    TH1F *jet_cut[30];
    TH1F *b_jet_cut[30];
    TH1F *MET_cut[30];

  };

  struct Lepton {//The struct 'Lepton' can store the following variables:
    TLorentzVector v;
    int id;  int ind;
    float wt;
    int status; 
    int momid;
    int cutbased;
    //1 = tight
    //-1 = loose
  };

  void Sortpt(vector<Lepton> vec);
  
protected:
  Hists h;
  
private:
  //Global variables go here. Make them global only if necessary.
  TFile *_HstFile;
  const char *_HstFileName;
  const char *_SumFileName;
  int _verbosity,_exclude;
  int _data, _lep, _year, _sample;

  //Event counters can be declared here.
  int nEvtTotal, nEvtRan;
  bool GoodEvt, GoodEvt2016, GoodEvt2017, GoodEvt2018; //Flags
  TString _era;
  
  int n1mu2J ,n1mu2J_cut;
  
  Lepton mylep, myj0,myj1,mybjet;
  
  
  
  
  //######################
  // Declare arrays here:
  //######################
  vector<Lepton> goodMu;
  vector<Lepton> goodEle;
  vector<Lepton> leptons;
  vector<Lepton> goodPhoton;
  vector<Lepton> goodJets;
  vector<Lepton> genEle;
  vector<Lepton> genMu;
  vector<Lepton> goodbJets;

  vector<Lepton> goodMu_tight;
  vector<Lepton> goodEle_loose;
  ClassDef(SingleTopAna,0);
};

#endif

#ifdef SingleTopAna_cxx
void SingleTopAna::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  fReader.SetTree(tree);
}

Bool_t SingleTopAna::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}



#endif // #ifdef SingleTopAna_cxx
