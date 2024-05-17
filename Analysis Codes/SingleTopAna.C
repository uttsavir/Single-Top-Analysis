//***************************************************
// Code for Single Top Analysis**********************
//***************************************************

#define SingleTopAna_cxx
#include "SingleTopAna.h"
#include <TH2.h>
#include <TStyle.h>

void SingleTopAna::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  
  TString option = GetOption();
}

void SingleTopAna::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  
  TString option = GetOption();
  
  //Initialization of the counters:
  nEvtRan        = 0;
  nEvtTotal      = 0;
  n1mu2J         = 0;
  n1mu2J_cut  = 0;
  //Other custom counters can be initialized here.
  
  _HstFile = new TFile(_HstFileName,"recreate");
  BookHistograms();
}

void SingleTopAna::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  
  _HstFile->Write();
  _HstFile->Close();

  //The following lines are displayed on the root prompt.
  cout<<"Total events ran = "<<nEvtRan<<endl;
  cout<<"Total good events = "<<nEvtTotal<<endl;
  cout<<"------------------------------------------------------------------"<<endl;
  cout<<"Total events with final state 1Muon and 2jets(>1bjet) = "<<n1mu2J<<endl;
  cout<<"------------------------------------------------------------------"<<endl;
  cout<<"Total events with final state 1Muon and 2jets(>1bjet)-CUT = "<<n1mu2J_cut<<endl;
  
  //The following lines are written on the sum_<process name>.txt file
  ofstream fout(_SumFileName);
  fout<<"Total events ran = "<<nEvtRan<<endl;
  fout<<"Total good events  = "<<nEvtTotal<<endl;
}

void SingleTopAna::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
}

Bool_t SingleTopAna::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  
  // The processing can be stopped by calling Abort().
  
  // Use fStatus to set the return value of TTree::Process().
  
  // The return value is currently not used.

  fReader.SetLocalEntry(entry);
  if(_data == 0)
    fReader_MC.SetLocalEntry(entry);
  if(_data == 1)
    fReader_Data.SetLocalEntry(entry);

  //Verbosity determines the number of processed events after which the root prompt is supposed to display a status update.
  if(_verbosity==0 && nEvtTotal%20000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;      
  else if(_verbosity>0 && nEvtTotal%20000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;
  
  //The following flags throws away some events based on unwanted properties (such as detector problems)
  GoodEvt2018 = (_year==2018 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2017 = (_year==2017 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2016 = (_year==2016 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  
  GoodEvt = GoodEvt2018 && GoodEvt2017 && GoodEvt2016;
  
  nEvtRan++;                             //Total number of events containing everything (including the trash events).
  h.nEvtRan->Fill(1); 
  
  if(GoodEvt){
    nEvtTotal++;                         //Total number of events containing goodEvents
    h.nEvtTotal->Fill(1); 
    //The analysis is done for these good events.                                    

    //###############################
    //Construction of the arrays:
    // 1. goodMu
    // 2. goodEle
    // 3. Lepton
    // 4. goodJets
    //#################################


    //***************************************
    //goodMu array****************************
    //****************************************
    int nmu = 0;                         // This counts the number of muons in each event.
    goodMu.clear(), leptons.clear(), goodMu_tight.clear();
    // Make sure to empty the array from previous event.
    for(unsigned int i=0; i<(*nMuon); i++){
      // This loop runs over all the muon candidates. Some of them will pass our selection criteria.
      // These will be stored in the goodMu array.
      Lepton temp;                       // 'temp' is the i-th candidate.
      temp.v.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],0.105); //the muon mass in GeV is 0.105
      temp.id = -13*Muon_charge[i];    //pdgID for mu- = 13, pdgID for mu+ = -13  
      temp.ind = i; 
      //temp.cutbased;
      
      //These are the flags the 'temp' object i.e. the muon candidate has to pass.
      /*    bool passCuts =temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && Muon_tightId[i];
      passCuts = passCuts && Muon_pfRelIso04_all[i]<0.15;
      //passCuts = passCuts && fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;
      passCuts = passCuts && fabs(Muon_dxy[i])<0.01 && fabs(Muon_dz[i])<0.02;

      bool passCuts_loose = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && Muon_looseId[i] && Muon_pfRelIso04_all[i]<0.15 && fabs(Muon_dxy[i])<0.01 && fabs(Muon_dz[i])<0.02;
      */
      bool passCut =  temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && Muon_pfRelIso04_all[i]<0.2 && fabs(Muon_dxy[i])<0.01 && fabs(Muon_dz[i])<0.02;
      
      bool passCuts_tight =  passCut && Muon_tightId[i] ;
      bool passCuts_loose = passCut && Muon_looseId[i];

      if(passCuts_loose){
	goodMu.push_back(temp);
	leptons.push_back(temp);
      }

      if(passCuts_tight){
	goodMu_tight.push_back(temp); 
      }
    } // This 'for' loop has created a goodMu array.
    
    //Now we sort in decreasing pT order
    Sortpt(goodMu);//loose muons
    Sortpt(goodMu_tight);//tight muons


    
    //**************************************
    //goodEle array*************************
    //**************************************
    int ne = 0;                         
    goodEle.clear();                      
    for(unsigned int i=0; i<(*nElectron); i++){
      Lepton temp;                       
      temp.v.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],0.000511);
      //the electron mass = 0.511 MeV = 0.000511 GeV
      //pdgID for e- = 11, pdgID for e+ = -11  
      temp.id = -11*Electron_charge[i];    
      temp.ind = i;
      
      //These are the flags the 'temp' object i.e. the electron candidate has to pass.
      bool passCuts_tight = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && Electron_cutBased[i]>2;
     
      
      bool isprompt;
      if(fabs(temp.v.Eta())<=1.479){//for barrel
	if(fabs(Electron_dxy[i])<0.05 && fabs(Electron_dz[i])<0.1)
	  isprompt = true;
      }
      if(fabs(temp.v.Eta())>1.479){//endcap
	if(fabs(Electron_dxy[i])<0.1 && fabs(Electron_dz[i])<0.2)
	  isprompt = true;
    }
      
      // cleaning the electrons against muons-------------------
      bool is_muon_close = false;//muon is not close
      //float dR =0 ;
      for(int i=0; i<((int)goodMu.size()); i++){
	float dR = temp.v.DeltaR(goodMu.at(i).v);
	if(dR<0.4){
	  is_muon_close = true;
	}
      }
      bool passCuts_loose = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && Electron_cutBased[i]>0;
      
      if(passCuts_loose && isprompt && !is_muon_close){
	goodEle.push_back(temp);//loose ele
	leptons.push_back(temp);
      }
    }

    //Now we sort  in decreasing pT order
    Sortpt(goodEle);
   
    Sortpt(leptons);


    
    //*******************************************
    //goodJets array*****************************
    //*******************************************
    int njet = 0;                         
    goodJets.clear();
    goodbJets.clear();
    for(unsigned int i=0; i<(*nJet); i++){
      Lepton temp;                            
      temp.v.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i], Jet_mass[i]); 
      temp.ind = i; 
      
      //These are the flags the 'temp' object i.e. the electron candidate has to pass.
      bool passCuts = temp.v.Pt()>30 && fabs(temp.v.Eta())<2.4;
      
      // cleaning the jets against the electrons--------------
      bool is_e_close = false;//ele is not close
      //float dR =0.0 ;
      for(int i=0; i<((int)goodEle.size()); i++){
	float dR = temp.v.DeltaR(goodEle.at(i).v);
	if(dR<0.4){
	  is_e_close = true;
	}
      }
      // cleaning the jets against muons--------------
      bool is_mu_close = false;//muon is not close
      for(int i=0; i<((int)goodMu.size()); i++){
	float dR = temp.v.DeltaR(goodMu.at(i).v);
	if(dR<0.4){
	  is_mu_close = true;
	}
      }
      
      if(passCuts && !is_mu_close && !is_e_close){
	// If 'temp' satisfies all the conditions, it is pushed back into goodJets
	if(_year == 2016 ? Jet_jetId[i]>=1 : Jet_jetId[i]>=2)
	  goodJets.push_back(temp);
	if(Jet_btagDeepFlavB[i]>0.2783)
	  goodbJets.push_back(temp);
      }
    }
    
    Sortpt(goodJets);
    Sortpt(goodbJets);




    
    //-------------------------------------------------------------------------------------------------------------
    //                                                  Analysis Block
    //--------------------------------------------------------------------------------------------------------------
    h.nmuons_tight->Fill((int)goodMu_tight.size());
        

    //For Leptons------------------------------------------------------
    h.nLep->Fill((int)leptons.size());
    for(int i=0; i<(int)leptons.size(); i++){
      h.leppt->Fill(leptons.at(i).v.Pt());
      h.lepeta->Fill(leptons.at(i).v.Eta());
      h.lepphi->Fill(leptons.at(i).v.Phi());
      h.cut_based->Fill(leptons.at(i).cutbased); 
    }
    
    //For muons--------------------------------------------------------  
    h.nmuons->Fill((int)goodMu.size());
    for(int i=0; i<(int)goodMu.size(); i++){
      h.mupt->Fill(goodMu.at(i).v.Pt());
      h.mueta->Fill(goodMu.at(i).v.Eta());
      h.muphi->Fill(goodMu.at(i).v.Phi());
    }
    
    //For electrons---------------------------------------------------------  
    h.nelectrons->Fill((int)goodEle.size());
    
    for(int i=0; i<(int)goodEle.size(); i++){
      h.elept->Fill(goodEle.at(i).v.Pt());
      h.eleeta->Fill(goodEle.at(i).v.Eta());
      h.elephi->Fill(goodEle.at(i).v.Phi());
    }

    if((int)goodEle.size()>1){
      h.ele[0]->Fill(goodEle.at(0).v.Pt());
      h.ele[1]->Fill(goodEle.at(0).v.Eta());
      h.ele[2]->Fill(goodEle.at(0).v.Phi());
      h.ele[3]->Fill(goodEle.at(1).v.Pt());
      h.ele[4]->Fill(goodEle.at(1).v.Eta());
      h.ele[5]->Fill(goodEle.at(1).v.Phi());
    }
   
     
    //For Jets--------------------------------------------------------------------
    
    h.njets->Fill((int)goodJets.size());
    
    for(int i=0; i<(int)goodJets.size(); i++){
      h.jetpt->Fill(goodJets.at(i).v.Pt());
      h.jeteta->Fill(goodJets.at(i).v.Eta());
      h.jetphi->Fill(goodJets.at(i).v.Phi());
     
      float jetscore_deepB = Jet_btagDeepB[goodJets.at(i).ind];
      h.jet_score -> Fill(jetscore_deepB);
    }
    //b jets--------------------------------
    h.n_bjets -> Fill((int)goodbJets.size());
    for(int i=0; i<(int)goodbJets.size(); i++){
      h.b_jetpt->Fill(goodbJets.at(i).v.Pt());
      h.b_jeteta->Fill(goodbJets.at(i).v.Eta());
      h.b_jetphi->Fill(goodbJets.at(i).v.Phi());
    }

    
    //----------------------------------------------------------------------------------------------------------
    //             EVENT SELECTION (1mu2J(bJet>0) - inclusive on jets only
    //---------------------------------------------------------------------------------------------------------
    
    bool   trigger1 = false;
    bool is_1mu2j      = false;
    bool trigger2 = false;
        
    if((int)goodMu.size()==1&& (int)goodJets.size()>1){// exactly muon and atleast 2 jets
      // cout<< "werking1"<<endl;
      if((int)goodMu_tight.size() ==1){//the 1 mu has to be a tight muon
      	//cout<<"working2" <<endl;  
	if((int)goodEle.size()==0){//no electrons
	  //cout<<"working3" <<endl;  
	  mylep=goodMu.at(0);
	  myj0=goodJets.at(0);
	  myj1=goodJets.at(1);
	  if(goodMu.at(0).v.Pt() > 26){
	    // cout<<"working1" <<endl;  
	    trigger1 = true;
	    if(fabs(Muon_pfRelIso04_all[mylep.ind])<0.15){
	       //&& fabs(Muon_sip3d[mylep.ind])<4){
	      //cout<<"working2" <<endl;  
	      trigger2 = true;
	      if(trigger1 && trigger2){
		//cout<<"working3" <<endl;  
		is_1mu2j =true;
	      }
	    }
	  }
	}
      }
    }

//--------------------------------------------------------------------------------------------
    if(is_1mu2j){
      //  cout<<"working7" <<endl;  
      //corrections-----------     
      double evt_wt=1.0;
  
      if(_data==0 && (int)goodMu.size()>0){	//MC
	
	//scale factor---------------
	//mu0--------------------
	//double SFReco_mu0= Muon_2018UL_Reco(goodMu.at(0).v.Pt(),goodMu.at(0).v.Eta());
	//double SFID_mu0  = Muon_2018UL_ID(goodMu.at(0).v.Pt(),goodMu.at(0).v.Eta());
	//double SFIso_mu0 = Muon_2018UL_Iso(goodMu.at(0).v.Pt(),goodMu.at(0).v.Eta());
	//----------
	//double SF_mu0 = SFReco_mu0 * SFID_mu0 * SFIso_mu0 ;
	/*
	//mu1------------------
	double SFReco_mu1= Muon_2018UL_Reco(goodMu.at(1).v.Pt(),goodMu.at(1).v.Eta());
	double SFID_mu1  = Muon_2018UL_ID(goodMu.at(1).v.Pt(),goodMu.at(1).v.Eta());
	double SFIso_mu1 = Muon_2018UL_Iso(goodMu.at(1).v.Pt(),goodMu.at(1).v.Eta());
	//--------
	double SF_mu1 = SFReco_mu1 * SFID_mu1 * SFIso_mu1 ;
	
	//mu2------------------
	double SFReco_mu2= Muon_2018UL_Reco(goodMu.at(2).v.Pt(),goodMu.at(2).v.Eta());
	double SFID_mu2  = Muon_2018UL_ID(goodMu.at(2).v.Pt(),goodMu.at(2).v.Eta());
	double SFIso_mu2 = Muon_2018UL_Iso(goodMu.at(2).v.Pt(),goodMu.at(2).v.Eta());
	//--------
	double SF_mu2 = SFReco_mu2 * SFID_mu2 * SFIso_mu2 ;
	
	//TOTAL SCALE FACTOR
	double ScaleFactor = SF_mu0 * SF_mu1 * SF_mu2 ; 
	*/
	/*double ScaleFactor = SF_mu0 ;
	  
	//trigger effi----------------------------
	double trigger_effi_mu0 = TrigEff_2018_IsoMu24_Data(goodMu.at(0).v.Pt(),goodMu.at(0).v.Eta());
	//	double trigger_effi_mu1 = TrigEff_2018_IsoMu24_Data(goodMu.at(1).v.Pt(),goodMu.at(1).v.Eta());
	//	double trigger_effi_mu2 = TrigEff_2018_IsoMu24_Data(goodMu.at(2).v.Pt(),goodMu.at(2).v.Eta());
	//TOTAL SCALE FACTOR
	//	double TriggerEfficiency =1- ((1- trigger_effi_mu0) * (1- trigger_effi_mu1) * (1- trigger_effi_mu2));
	double TriggerEfficiency = trigger_effi_mu0;
	
	//eventweight = SF times TriggerEff
	evt_wt = ScaleFactor * TriggerEfficiency;*/
	evt_wt = 1.0;
      }//data0
      
      if(_data==1){	//Data
	evt_wt = 1.0;
      }//data1
      
      h.evt_wt->Fill(evt_wt);
      n1mu2J++;
      h.n1mu2J->Fill(1);//filling no of 1l2j events in a histogram
      
      
      h.elesize->Fill((int)goodEle.size(),evt_wt);
      h.lepsize->Fill((int)leptons.size(),evt_wt);
  
      for(int i=0; i<(int)leptons.size(); i++){
	h.cutbased_baseline->Fill(leptons.at(i).cutbased); 
      }
      
      //muons
      h.mu[0]->Fill((int)goodMu.size(),evt_wt);
      h.mu[1]->Fill(mylep.v.Pt(),evt_wt);
      h.mu[2]->Fill(mylep.v.Eta(),evt_wt);
      h.mu[3]->Fill(mylep.v.Phi(),evt_wt);
      
      //if((int)leptons.size()>0){
      //here if we use if leptons.size()>1, if will give us an incorrect answer because it only selects event which have atleast 2 leps- and neglects events which have 1 lep. but in case of jets using goodJet.size()>1, is alright because all our events have atleast 2 jets(becuse of our event selection)
      
      //Ploting invariant mass
      h.mu[4]-> Fill((mylep.v + myj0.v).M(),evt_wt);
      h.mu[5]-> Fill((mylep.v + myj1.v).M(),evt_wt);

      //isolation
      h.mu[9]->Fill( Muon_pfRelIso04_all[mylep.ind],evt_wt );
      //IMPACT PARAMETER---------------
      h.mu[10]->Fill(Muon_dxy[mylep.ind],evt_wt);	
      // dz var
      h.mu[11]->Fill(Muon_dz[mylep.ind],evt_wt);
      //3d iso
      h.mu[12]->Fill(Muon_sip3d[mylep.ind],evt_wt);
      //JETS
      h.jet[0]->Fill((int)goodJets.size(),evt_wt);
      h.jet[1]->Fill(myj0.v.Pt(),evt_wt);
      h.jet[2]->Fill(myj0.v.Eta(),evt_wt);
      h.jet[3]->Fill(myj0.v.Phi(),evt_wt);
	  
      h.jet[4]->Fill(myj1.v.Pt(),evt_wt);
      h.jet[5]->Fill(myj1.v.Eta(),evt_wt);
      h.jet[6]->Fill(myj1.v.Phi(),evt_wt);
    
      
      
      //  h.jet[8] ->Fill(myj0.v.DeltaPhi(myj1.v),evt_wt); // angle between the leading and subleading jet
      float dphi_J0J1 = delta_phi((myj0.v.Phi()) , (myj1.v.Phi()));
      h.jet[8] ->Fill(dphi_J0J1, evt_wt); // angle between the leading and subleading jet
      h.jet[9] ->Fill(myj0.v.DeltaR(myj1.v),evt_wt);// delta R (angular distance) between the leading and subleading jets
      h.jet[10]->Fill((myj0.v + myj1.v).M(),evt_wt);//Ploting invariant mass of leading & subleading jet in each event	
      h.jet[11]->Fill((myj0.v + myj1.v).Pt(),evt_wt); //dijet pt
      

      //2jet system
      float pT_J0J1  = (myj0.v + myj1.v).Pt(); //pt of the 2 jet system
      float phi_J0J1 = (myj0.v + myj1.v).Phi(); //phi of the 2 jet system
      
      float dphi_J0J1_mu0= delta_phi(phi_J0J1, (mylep.v.Phi()));
      h.jet[17] ->Fill(dphi_J0J1_mu0,evt_wt);
      
      //dr and dphi
      float dR_mu0J0= mylep.v.DeltaR(myj0.v);
      h.jet[12] ->Fill(dR_mu0J0,evt_wt);
      float dR_mu0J1= mylep.v.DeltaR(myj1.v);
      h.jet[13] ->Fill(dR_mu0J1,evt_wt);
      //      float dPhi_mu0J0= mylep.v.DeltaPhi(myj0.v);
      float dPhi_mu0J0= delta_phi((mylep.v.Phi()) , (myj0.v.Phi()));
      h.jet[14] ->Fill(dPhi_mu0J0,evt_wt);
      //      float dPhi_mu0J1= mylep.v.DeltaPhi(myj1.v);
      float dPhi_mu0J1= delta_phi((mylep.v.Phi()) , (myj1.v.Phi()));
      h.jet[15] ->Fill(dPhi_mu0J1,evt_wt);
      
            
      //jet score
      float jetscore_sum= Jet_btagDeepFlavB[myj0.ind] + Jet_btagDeepFlavB[myj1.ind];
      float jetscore_product= Jet_btagDeepFlavB[myj0.ind] * Jet_btagDeepFlavB[myj1.ind];
      
      h.jet[18] -> Fill(jetscore_sum,evt_wt);     
      h.jet[19] -> Fill(jetscore_product,evt_wt);
     

      //bJets
      h.b_jet[0]->Fill((int)goodbJets.size(),evt_wt);
      if((int)goodbJets.size()>0){
	h.b_jet[1]->Fill(goodbJets.at(0).v.Pt(),evt_wt);
	h.b_jet[2]->Fill(goodbJets.at(0).v.Eta(),evt_wt);
	h.b_jet[3]->Fill(goodbJets.at(0).v.Phi(),evt_wt);
	

	h.b_jet[4]-> Fill((mylep.v + goodbJets.at(0).v).M(),evt_wt);
	//-----------------------------------------------------------------------------------------------------------
	float dR_mu0bJ0 = mylep.v.DeltaR(goodbJets.at(0).v);
	h.b_jet[5] ->Fill(dR_mu0bJ0,evt_wt);
	
	float dphi_mu0bjet = delta_phi((mylep.v.Phi()) , (goodbJets.at(0).v.Phi()));
	h.b_jet[6] ->Fill(dphi_mu0bjet , evt_wt);
	
	//
	float dphi_bJ0met = delta_phi( (goodbJets.at(0).v.Phi()) , (*PuppiMET_phi));
	float mT_bjet =transv_mass((goodbJets.at(0).v.Pt()),(*PuppiMET_pt) , dphi_bJ0met);
	h.MET[5] -> Fill(mT_bjet,evt_wt);
	h.MET[9] ->Fill(dphi_bJ0met,evt_wt);
      }
           
      //For MET ---------------------------------------
      h.MET[0]->Fill (*PuppiMET_pt,evt_wt);
      h.MET[1]->Fill (*PuppiMET_phi,evt_wt);

      
      float dphi_mu0met = delta_phi((mylep.v.Phi()),(*PuppiMET_phi));
      float mT_mu0 =transv_mass((mylep.v.Pt()),(*PuppiMET_pt) , dphi_mu0met);
      h.MET[2] -> Fill(mT_mu0,evt_wt);
      h.MET[6] -> Fill(dphi_mu0met,evt_wt);
      
      
      //jets
      float dphi_J0met = delta_phi((myj0.v.Phi()) , (*PuppiMET_phi));
      float mT_J0 =transv_mass((myj0.v.Pt()),(*PuppiMET_pt) , dphi_J0met);
      h.MET[3] -> Fill(mT_J0,evt_wt);
      h.MET[7] ->Fill(dphi_J0met,evt_wt);
      
      float dphi_J1met = delta_phi((myj1.v.Phi()) , (*PuppiMET_phi));
      float mT_J1 =transv_mass((myj1.v.Pt()),(*PuppiMET_pt) , dphi_J1met);
      h.MET[4] ->Fill(mT_J1,evt_wt);
      h.MET[8] ->Fill(dphi_J1met,evt_wt); 
      
      //2 jet system 
      float dphi_J0J1met= delta_phi(phi_J0J1, (*PuppiMET_phi));
      h.MET[10] ->Fill(dphi_J0J1met,evt_wt);
      
      float mT_J0J1met =transv_mass((pT_J0J1),(*PuppiMET_pt) , dphi_J0met);
      h.MET[11] -> Fill(mT_J0J1met,evt_wt);
    
      
      //Sum of pT--------------
      //ST= sum of HT, MET, LT
      
      float Ht = 0.0;
      float St = 0.0;
      float Lt = 0.0;
      
      //if((int)goodJets.size()>0 && (int)goodMu.size()>0){
      for(int j = 0;j< (int)goodJets.size();j++){
	Ht = Ht + goodJets.at(j).v.Pt();//scalar sum of pT of all the jets
      }
      for(int i=0; i<(int)goodMu.size(); i++){
	Lt =  leptons.at(i).v.Pt();//scalar sum of pT of all leptons
	St = Ht +Lt;
      }
      
      h.jet[7]->Fill(Ht,evt_wt);
      h.mu[6]->Fill(Lt,evt_wt);
      
      St = St + *PuppiMET_pt;
      //scalar sum of pT
      h.mu[7]->Fill(St,evt_wt);
	
      //Ht of 2 jets
      float Ht_j0j1 =0.0;
      for(int j=0; j<2; j++){
	Ht_j0j1 = Ht_j0j1 + goodJets.at(j).v.Pt(); 
      }
      h.jet[16]->Fill(Ht_j0j1,evt_wt);
      
      //HT50  HT70 
      float Ht50 =0.0;
      float Ht70 =0.0;
      
      for(int j = 0;j< (int)goodJets.size();j++){
	if(goodJets.at(j).v.Pt() > 50)
	  Ht50 = Ht50 + goodJets.at(j).v.Pt();//scalar sum of pT of all the jets having pT>50
	if(goodJets.at(j).v.Pt() > 70)
	  Ht70 = Ht70 + goodJets.at(j).v.Pt();//scalar sum of pT of all the jets having pT>70
      }
      h.jet[20]->Fill(Ht50,evt_wt);
      h.jet[21]->Fill(Ht70,evt_wt);
      
      //momentum in the z direction (Pz)
      float vector_add_Pz = (mylep.v + myj0.v + myj1.v).Pz();
      h.mu[8]->Fill(vector_add_Pz,evt_wt);
      
      
      //---------------------------------------------------------------------------------------------------------------------
      //                             CUTS
      //----------------------------------------------------------------------------------------------------------------------
      
      if(Ht>150 &&
	 jetscore_sum>1 &&
	 Ht70>20){

	n1mu2J_cut++;
	//	 cout<<"test"<<n1mu2J_cut<<endl;
	h.n1mu2J_cut->Fill(2);//filling no of 1l2j events after cut in a histogram
	
	h.evt_wt_cut->Fill(evt_wt);
	
	//muons
	h.mu_cut[0]->Fill((int)goodMu.size(),evt_wt);
	h.mu_cut[1]->Fill(mylep.v.Pt(),evt_wt);
	h.mu_cut[2]->Fill(mylep.v.Eta(),evt_wt);
	h.mu_cut[3]->Fill(mylep.v.Phi(),evt_wt);
           
	//Ploting invariant mass
	h.mu_cut[4]-> Fill((mylep.v + myj0.v).M(),evt_wt);
	h.mu_cut[5]-> Fill((mylep.v + myj1.v).M(),evt_wt);

	//isolation
	h.mu_cut[9]->Fill( Muon_pfRelIso04_all[mylep.ind],evt_wt );
	//IMPACT PARAMETER---------------
	h.mu_cut[10]->Fill(Muon_dxy[mylep.ind],evt_wt);	
	// dz var
	h.mu_cut[11]->Fill(Muon_dz[mylep.ind],evt_wt);
	//3d iso
	h.mu_cut[12]->Fill(Muon_sip3d[mylep.ind],evt_wt);


	//JETS
	h.jet_cut[0]->Fill((int)goodJets.size(),evt_wt);
	h.jet_cut[1]->Fill(myj0.v.Pt(),evt_wt);
	h.jet_cut[2]->Fill(myj0.v.Eta(),evt_wt);
	h.jet_cut[3]->Fill(myj0.v.Phi(),evt_wt);
	  
	h.jet_cut[4]->Fill(myj1.v.Pt(),evt_wt);
	h.jet_cut[5]->Fill(myj1.v.Eta(),evt_wt);
	h.jet_cut[6]->Fill(myj1.v.Phi(),evt_wt);
	
	h.jet_cut[8] ->Fill(dphi_J0J1, evt_wt); // angle between the leading and subleading jet
	h.jet_cut[9] ->Fill(myj0.v.DeltaR(myj1.v),evt_wt);// delta R (angular distance) between the leading and subleading jets
	h.jet_cut[10]->Fill((myj0.v + myj1.v).M(),evt_wt);//Ploting invariant mass of leading & subleading jet in each event	
	h.jet_cut[11]->Fill((myj0.v + myj1.v).Pt(),evt_wt); //dijet pt
      

	//2jet system
	h.jet_cut[17] ->Fill(dphi_J0J1_mu0,evt_wt);
      
	//dr and dphi
	h.jet_cut[12] ->Fill(dR_mu0J0,evt_wt);
	h.jet_cut[13] ->Fill(dR_mu0J1,evt_wt);
	
	h.jet_cut[14] ->Fill(dPhi_mu0J0,evt_wt);
	
	h.jet_cut[15] ->Fill(dPhi_mu0J1,evt_wt);
                  
	//jet score
	h.jet_cut[18] -> Fill(jetscore_sum,evt_wt);     
	h.jet_cut[19] -> Fill(jetscore_product,evt_wt);
     
	
	//bJets
	h.b_jet_cut[0]->Fill((int)goodbJets.size(),evt_wt);
	if((int)goodbJets.size()>0){
	  h.b_jet_cut[1]->Fill(goodbJets.at(0).v.Pt(),evt_wt);
	  h.b_jet_cut[2]->Fill(goodbJets.at(0).v.Eta(),evt_wt);
	  h.b_jet_cut[3]->Fill(goodbJets.at(0).v.Phi(),evt_wt);
	  h.b_jet_cut[4]-> Fill((mylep.v + goodbJets.at(0).v).M(),evt_wt);
	
	  float dR_mu0bJ0 = mylep.v.DeltaR(goodbJets.at(0).v);
	  h.b_jet_cut[5] ->Fill(dR_mu0bJ0,evt_wt);
	
	  float dphi_mu0bjet = delta_phi((mylep.v.Phi()) , (goodbJets.at(0).v.Phi()));
	  h.b_jet_cut[6] ->Fill(dphi_mu0bjet , evt_wt);
	
	  //met
	  float dphi_bJ0met = delta_phi( (goodbJets.at(0).v.Phi()) , (*PuppiMET_phi));
	  float mT_bjet =transv_mass((goodbJets.at(0).v.Pt()),(*PuppiMET_pt) , dphi_bJ0met);
	  h.MET_cut[5] -> Fill(mT_bjet,evt_wt);
	  h.MET_cut[9] ->Fill(dphi_bJ0met,evt_wt);
	}
	
	//For MET ---------------------------------------
	h.MET_cut[0]->Fill (*PuppiMET_pt,evt_wt);
	h.MET_cut[1]->Fill (*PuppiMET_phi,evt_wt);

	//muon
      	h.MET_cut[2] -> Fill(mT_mu0,evt_wt);
	h.MET_cut[6] -> Fill(dphi_mu0met,evt_wt);
            
	//jets
	h.MET_cut[3] -> Fill(mT_J0,evt_wt);
	h.MET_cut[7] ->Fill(dphi_J0met,evt_wt);
      
	h.MET_cut[4] ->Fill(mT_J1,evt_wt);
	h.MET_cut[8] ->Fill(dphi_J1met,evt_wt); 
      
	//2 jet system 
	h.MET_cut[10] ->Fill(dphi_J0J1met,evt_wt);
      	h.MET_cut[11] -> Fill(mT_J0J1met,evt_wt);

      
	//Sum of pT--------------
	//ST= sum of HT, MET, LT
      	h.jet_cut[7]->Fill(Ht,evt_wt);
	h.mu_cut[6]->Fill(Lt,evt_wt);
	
	//scalar sum of pT
	h.mu_cut[7]->Fill(St,evt_wt);
	
	//Ht of 2 jets
	h.jet_cut[16]->Fill(Ht_j0j1,evt_wt);
	
	//HT50  HT70 
	h.jet_cut[20]->Fill(Ht50,evt_wt);
	h.jet_cut[21]->Fill(Ht70,evt_wt);
	
	//momentum in the z direction (Pz)
	//float vector_add_Pz = (mylep.v + myj0.v + myj1.v).Pz();
	h.mu_cut[8]->Fill(vector_add_Pz,evt_wt);
	
      }//cuts
    }//is_1mu2j	
    
    //############################END#######################################
    
  }//GoodEvt
  return kTRUE;  
  
}//Process



//---------------------------------------------------------------------------------------------------------------------
//                                                          Functions
//---------------------------------------------------------------------------------------------------------------------
//######################################
//USER DEFINED FUNCTIONS
//1. Sortpt
//2. get_mother
//3. GenMother
//4. delta_phi
//5. transv_mass
//6. scale factor muon reco
//7. scale factor muon id
//8. scale factor muon iso
//9. trigger_efficiency
//######################################

void SingleTopAna::Sortpt(vector<Lepton> vec)
{
  
  for(int i=0; i<(int)vec.size()-1; i++){
    for(int j=i+1; j<(int)vec.size(); j++){
      if( vec[i].v.Pt() < vec[j].v.Pt() ) swap(vec.at(i), vec.at(j));
    }
  }
}

/*
int SingleTopAna::get_mother(int i)
{
  int pid = GenPart_pdgId[i];
  int motherid = GenPart_pdgId[GenPart_genPartIdxMother[i]];
  
  return motherid;
}

int SingleTopAna::GenMother(int ind, int mom_ind)
{
  int p_id = GenPart_pdgId[ind];
  int m_id = GenPart_pdgId[mom_ind];
  while(p_id==m_id){
    ind = mom_ind;
    mom_ind = GenPart_genPartIdxMother[ind];
    p_id = GenPart_pdgId[ind];
    m_id = GenPart_pdgId[mom_ind];
  }
  return m_id;
}
*/
float SingleTopAna::delta_phi(float phi1, float phi2)
{
  //The correct deltaPhi falls in the interval [0 , pi]
  phi1 = TVector2::Phi_0_2pi(phi1);
  phi2 = TVector2::Phi_0_2pi(phi2);
  float dphi = fabs(phi1 - phi2);
  if(dphi>TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
  return dphi;
}

float SingleTopAna::transv_mass(float E_lep, float MET, float dphi)
{
  //The inputs are the Energy of the lepton, MET and dPhi between the lepton and MET
  float mT = sqrt(2* E_lep * MET *(1-cos(dphi)));
  return mT;
}



//Corrections-------------------------------------------------
//Muon Reco--------------------------------


double SingleTopAna::Muon_2018UL_Reco(float pt, float eta){
  double scale_factor = 1.0;//default value
  if(fabs(eta)<=0.9){
    if( 15<pt && pt<=20)         scale_factor =  1.04000400946163;
    else if( 20<pt && pt<=25)    scale_factor =  0.989014677272758;
    else if( 25<pt && pt<=30)    scale_factor =  0.9989055590645546;
    else if( 30<pt && pt<=40)    scale_factor =  0.999760853514322;
    else if( 40<pt && pt<=50)    scale_factor =  0.999565467845588;
    else if( 50<pt && pt<=60)    scale_factor =  0.9998204415944014;
    else if( 60<pt && pt<=120)   scale_factor =  0.9975015485587808;
    else if( 120<pt && pt<=3000) scale_factor =  0.9975015485587808;   
  }
  else if(0.9<fabs(eta) && fabs(eta)<=1.2){
    if( 15<pt && pt<=20)         scale_factor =  1.0334239991476872;
    else if( 20<pt && pt<=25)    scale_factor =  0.9832043508246824;
    else if( 25<pt && pt<=30)    scale_factor =  0.9962074764783356;
    else if( 30<pt && pt<=40)    scale_factor =  0.999066550372675;
    else if( 40<pt && pt<=50)    scale_factor =  0.999297799609978;
    else if( 50<pt && pt<=60)    scale_factor =  0.999898284881032;
    else if( 60<pt && pt<=120)   scale_factor =  0.9978985301114488;
    else if( 120<pt && pt<=3000) scale_factor =  0.9978985301114488;
  }
  else if(1.2<fabs(eta) && fabs(eta)<=2.1){
    if( 15<pt && pt<=20)         scale_factor =  0.9052340457842656;
    else if( 20<pt && pt<=25)    scale_factor =  0.9932902480824264;
    else if( 25<pt && pt<=30)    scale_factor =  0.9941744704916416;
    else if( 30<pt && pt<=40)    scale_factor =  0.9981802152409986 ;
    else if( 40<pt && pt<=50)    scale_factor =  0.9992437146078612;
    else if( 50<pt && pt<=60)    scale_factor =  0.9988904857161484;
    else if( 60<pt && pt<=120)   scale_factor =  0.9961937512394232;
    else if( 120<pt && pt<=3000) scale_factor =  0.9961937512394232;
  }
  else if(2.1<fabs(eta) && fabs(eta)<=2.4){
    if( 15<pt && pt<=20)         scale_factor =  0.9606433668106512;
    else if( 20<pt && pt<=25)    scale_factor =  0.9911896930756448;
    else if( 25<pt && pt<=30)    scale_factor =  0.9910856248132396;
    else if( 30<pt && pt<=40)    scale_factor =  0.9907553768818946;
    else if( 40<pt && pt<=50)    scale_factor =  0.9933622417972257;
    else if( 50<pt && pt<=60)    scale_factor =  0.9972680454815808;
    else if( 60<pt && pt<=120)   scale_factor =  0.9665749388673832;
    else if( 120<pt && pt<=3000) scale_factor =  0.9665749388673832;
  }
  return scale_factor;
}


//Muon id---------------------------------

double SingleTopAna::Muon_2018UL_ID(float pt, float eta){
  double scale_factor = 1.0;//default value
  if(fabs(eta)<=0.9){
    if( 15<pt && pt<=20)         scale_factor =  0.997503630932556;
    else if( 20<pt && pt<=25)    scale_factor =  0.9970098193760408;
    else if( 25<pt && pt<=30)    scale_factor =  0.9967038731014324;
    else if( 30<pt && pt<=40)    scale_factor =  0.9968067958281034;
    else if( 40<pt && pt<=50)    scale_factor =  0.9967508548887346;
    else if( 50<pt && pt<=60)    scale_factor =  0.9966229487883712;
    else if( 60<pt && pt<=120)   scale_factor =  0.9960404563074258;
    else if( 120<pt && pt<=3000) scale_factor =  0.9960404563074258;   
  }
  
  else if(0.9<fabs(eta) && fabs(eta)<=1.2){
    if( 15<pt && pt<=20)         scale_factor =  0.9946687083751804;
    else if( 20<pt && pt<=25)    scale_factor =  0.9945901104821656;
    else if( 25<pt && pt<=30)    scale_factor =  0.9949123046983688;
    else if( 30<pt && pt<=40)    scale_factor =  0.9956851566638726;
    else if( 40<pt && pt<=50)    scale_factor =  0.995312757370098;
    else if( 50<pt && pt<=60)    scale_factor =  0.9953453170396492;
    else if( 60<pt && pt<=120)   scale_factor =  0.9952502092400578;
    else if( 120<pt && pt<=3000) scale_factor =  0.9952502092400578;
  }
  
  else if(1.2<fabs(eta) && fabs(eta)<=2.1){
    if( 15<pt && pt<=20)         scale_factor =  0.9935059740562387;
    else if( 20<pt && pt<=25)    scale_factor =  0.9938964039703309;
    else if( 25<pt && pt<=30)    scale_factor =  0.9940771317309104;
    else if( 30<pt && pt<=40)    scale_factor =  0.9941876412222066;
    else if( 40<pt && pt<=50)    scale_factor =  0.9947739560711352;
    else if( 50<pt && pt<=60)    scale_factor =  0.9944640815519576;
    else if( 60<pt && pt<=120)   scale_factor =  0.9944321015185944;
    else if( 120<pt && pt<=3000) scale_factor =  0.9944321015185944;
  }
  
  else if(2.1<fabs(eta) && fabs(eta)<=2.4){
    if( 15<pt && pt<=20)         scale_factor =  0.9741187787139896;
    else if( 20<pt && pt<=25)    scale_factor =  0.9727131632621572;
    else if( 25<pt && pt<=30)    scale_factor =  0.9739795198691632;
    else if( 30<pt && pt<=40)    scale_factor =  0.9750724343730348;
    else if( 40<pt && pt<=50)    scale_factor =  0.9747976871311668;
    else if( 50<pt && pt<=60)    scale_factor =  0.9743746068723314;
    else if( 60<pt && pt<=120)   scale_factor =  0.970208378910226;
    else if( 120<pt && pt<=3000) scale_factor =  0.970208378910226;
  }
  
  return scale_factor;
}

//Muon iso-----------------------
double SingleTopAna::Muon_2018UL_Iso(float pt, float eta){
  double scale_factor = 1.0;//default value
  if(fabs(eta)<=0.9){
    if( 15<pt && pt<=20)         scale_factor =  0.9855709248030112;
    else if( 20<pt && pt<=25)    scale_factor =  0.9915374149428736;
    else if( 25<pt && pt<=30)    scale_factor =  0.9892802138927208;
    else if( 30<pt && pt<=40)    scale_factor =  0.99429319088059;
    else if( 40<pt && pt<=50)    scale_factor =  0.9960381793609222;
    else if( 50<pt && pt<60)     scale_factor =  0.9966122162214248;
    else if( 60<pt && pt<=120)   scale_factor =  0.9982512434902386;
    else if( 120<pt && pt<=3000) scale_factor =  0.9982512434902386;
  }
  
  else if(0.9<fabs(eta) && fabs(eta)<=1.2){
    if( 15<pt && pt<=20)         scale_factor =  0.9704305474456316;
    else if( 20<pt && pt<=25)    scale_factor =  0.979495631900224;
    else if( 25<pt && pt<=30)    scale_factor =  0.986109906870666;
    else if( 30<pt && pt<=40)    scale_factor =  0.9882799134496284;
    else if( 40<pt && pt<=50)    scale_factor =  0.9920317410868984;
    else if( 50<pt && pt<=60)    scale_factor =  0.9944502663616424;
    else if( 60<pt && pt<=120)   scale_factor =  0.9964834155289044;
    else if( 120<pt && pt<=3000) scale_factor =  0.9964834155289044;
  }
  
  else if(1.2<fabs(eta) && fabs(eta)<=2.1){
    if( 15<pt && pt<=20)         scale_factor =  1.0033355765204894;
    else if( 20<pt && pt<=25)    scale_factor =  1.0015352730781717;
    else if( 25<pt && pt<=30)    scale_factor =  0.999623464660824;
    else if( 30<pt && pt<=40)    scale_factor =  0.9994357007740622;
    else if( 40<pt && pt<=50)    scale_factor =  0.99882097063955;
    else if( 50<pt && pt<=60)    scale_factor =  0.999185450924234;
    else if( 60<pt && pt<=120)   scale_factor =  0.9995171516432572;
    else if( 120<pt && pt<=3000) scale_factor =  0.9995171516432572;
  }
  
  else if(2.1<fabs(eta) && fabs(eta)<=2.4){
    if( 15<pt && pt<=20)         scale_factor =  1.0252542402601108;
    else if( 20<pt && pt<=25)    scale_factor =  1.0198418523288777;
    else if( 25<pt && pt<30)     scale_factor =  1.012454313839175;
    else if( 30<pt && pt<=40)    scale_factor =  1.008064837308479;
    else if( 40<pt && pt<=50)    scale_factor =  1.0037288955534218;
    else if( 50<pt && pt<=60)    scale_factor =  1.0028456983691447;
    else if( 60<pt && pt<=120)   scale_factor =  1.0029888816740478;
    else if( 120<pt && pt<=3000) scale_factor =  1.0029888816740478;
  }
  
  return scale_factor;
}


//trigger efficiency---------------------------
double SingleTopAna::TrigEff_2018_IsoMu24_Data(float pt, float eta){
  double eff = 0.0; //default value
  if( pt<10 || eta>2.4 ) return 0.0;
  else if(fabs(eta)<=1.479) eff = 0.5*0.950463*(1.0+TMath::Erf((pt-23.9593)/(2.0*0.375996)));
  else if(fabs(eta)>1.479)  eff = 0.5*0.953162*(1.0+TMath::Erf((pt-23.9459)/(2.0*0.457351)));
  return eff;
}


//#############################
//BOOKING OF HISTOGRAMS
//############################

void SingleTopAna::BookHistograms()
{
  //These histograms are stored in the hst_<process name>.root file in the same order.
  //Example : new TH1F ("hst_name", "hst title", NBins, startVal, EndVal);

  //#################################
  //         basic plots
  //#################################
  h.nmuons_tight   = new TH1F ("n_Mu_tight",           "Number of mu - tight",                               10,0,10);
  
  h.nEvtRan        =  new TH1F("n_EvtRan" ,                "No of events ran" ,                                     10,0,10);
  h.nEvtTotal      =  new TH1F("n_EvtTotal" ,              "No of good events " ,                                 10,0,10);
  h.nLep              =  new TH1F ("n_leptons",               "Number of Leptons",                                  10,0,10);

  h.leppt             = new TH1F("lep_pt",                       "Lepton p_{T}",                                           1000,0,1000);
  h.lepeta           = new TH1F("lep_Eta",                     "Lepton Eta",                                                64,-3.2,3.2);
  h.lepphi           = new TH1F("lep_Phi",                     "Lepton Phi",                                                64,-3.2,3.2);
  h.cut_based    = new TH1F ("lep_cutbased",           "Lepton Cutbased",                                     6,-3,3);

     
  h.nmuons       =   new TH1F("n_muons",                  "Number of Muons",                                    10, 0, 10);
  h.mupt           =   new TH1F("mu_pT",                       "Muon p_{T}",                                            1000,0,1000);
  h.mueta         =   new TH1F("mu_Eta",                     "Muon Eta",                                                 64,-3.2,3.2);
  h.muphi         =   new TH1F("mu_Phi",                      "Muon Phi",                                                64,-3.2,3.2);

  h.nelectrons = new TH1F("n_ele",                            "Number of Electrons",                             10, 0, 10);
  h.elept          = new TH1F("ele_pT",                          "Electron p_{T}",                                      1000,0,1000);
  h.eleeta        = new TH1F("ele_Eta",                         "Electron Eta",                                          64,-3.2,3.2);
  h.elephi        = new TH1F("ele_Phi",                         "Electron Phi",                                          64,-3.2,3.2);


  h.ele[0]       = new TH1F("ele0_pT",                         "Leading Electron p_{T}",                      1000,0,1000);
  h.ele[1]       = new TH1F("ele0_Eta",                        "Leading Electron Eta",                          64,-3.2,3.2);
  h.ele[2]       = new TH1F("ele0_Phi",                        "Leading Electron Phi",                          64,-3.2,3.2);
  h.ele[3]       = new TH1F("ele1_pT",                         "Sub Leading Electron p_{T}",             1000,0,1000);
  h.ele[4]       = new TH1F("ele1_Eta",                        "Sub Leading Electron Eta",                 64,-3.2,3.2);
  h.ele[5]       = new TH1F("ele1_Phi",                        "Sub Leading Electron Phi",                  64,-3.2,3.2);

  h.njets        = new TH1F("n_jets",                             "Number of Jets",                                    10,0,10);

  h.jetpt        = new TH1F("jet_pT",                              "Jet p_{T}",                                            1000,0,1000);
  h.jeteta      = new TH1F("jet_Eta",                             "Jet Eta",                                                64,-3.2,3.2);
  h.jetphi      = new TH1F("jet_Phi",                             "Jet Phi",                                                64,-3.2,3.2);
  h.jet_score = new TH1F("jet_score",                         "Jet Score",                                            1000,-2,2);


  h.n_bjets   = new TH1F("n_bjets",                            "Number of b Jets",                                10,0,10);
  h.b_jetpt   = new TH1F("b_jet_pT",                           "b Jet p_{T}",                                        1000,0,1000);
  h.b_jeteta = new TH1F("b_jet_Eta",                          "b Jet Eta",                                             64,-3.2,3.2);
  h.b_jetphi = new TH1F("b_jet_Phi",                          "b Jet Phi",                                              64,-3.2,3.2);

  //#################################
  //         baseline event selections
  //#################################
  h.n1mu2J  = new TH1F("n_1mu2j_baseline" ,           "No of events with 1Mu 2J" ,                10,0,10);
  h.evt_wt    = new TH1F("evt_wt_baseline",                "Event Weight",                                   2000, 0.85,1);

  
  // FOR LEPTONS---------------------------------------
 
  h.elesize =   new TH1F("n_ele_baseline",      "Number of ele- baseline",                            10, 0, 10);//this is after evt selection
  h.lepsize =  new TH1F("n_lep_baseline",      "Number of lep- baseline",                            10, 0, 10);//this is after event selection

  h.cutbased_baseline  = new TH1F ("lep_cutbased_baseline",   "Lepton Cutbased- baseline",                            6,-3,3);

  //FOR MUONS----------------------------------------------------------
  
  h.mu[0] = new TH1F("nMuon_baseline",        "No of Muons- baseline",                                10,0,10);
  h.mu[1] = new TH1F("mu0_pT_baseline",            "Leading Muon p_{T}- baseline",                                   1000,0,1000);
  h.mu[2] = new TH1F("mu0_eta_baseline",           "Leading Muon Eta- baseline",                                        64,-3.2,3.2);
  h.mu[3] = new TH1F("mu0_phi_baseline",           "Leading Muon Phi- baseline",                                        64,-3.2,3.2);

  h.mu[4] = new TH1F("invar_mass_mu0j0_baseline",  "Invariant Mass of Leading Muon and Leading Jet- baseline",            1000,0,1000);
  h.mu[5] = new TH1F("invar_mass_mu0j1_baseline",  "Invariant Mass of Leading Muon and Subleading Jet- baseline",       1000,0,1000);


  
  h.mu[6] = new TH1F("Lt_baseline" ,               "Scalar sum of pT of all Muon- baseline" ,                          1600,0,1600);
  h.mu[7] = new TH1F("ST_baseline" ,              "Scalar sum of pT- baseline" ,                                              1600,0,1600);


  h.mu[8] = new TH1F("objects_pz_baseline",            "P_{z} of all the Objects in the EveNt Selection- baseline",     1600,0,1600);


  h.mu[9] = new TH1F("Iso_mu0_baseline" ,               "Muon Isolation- baseline " ,            600,0,0.2);
  h.mu[10] = new TH1F("mu_dxy_baseline" ,               "Muon dxy- baseline" ,                           1000,-0.06,0.06);
  h.mu[11] = new TH1F("mu_dz_baseline" ,               "Muon dz- baseline" ,                          1000,-0.09,0.09);

  h.mu[12] = new TH1F("mu_sip3d_baseline" ,               "Muon 3D impact parameter significance- baseline" ,                          500,0,10);


  
  //FOR JETS------------------------------------------
 
    h.jet[0]  = new TH1F("n_jets_baseline",         "Number of  Jets- baseline ",                                           10,0,10);
  h.jet[1] = new TH1F("jet0_pT_baseline",             "Leading Jet p_{T}- baseline",                                               1000,0,1000);
  h.jet[2] = new TH1F("jet0_Eta_baseline",            "Leading Jet Eta- baseline",                                                 64,-3.2,3.2);
  h.jet[3] = new TH1F("jet0_Phi_baseline",            "Leading Jet Phi- baseline",                                                 64,-3.2,3.2);

  h.jet[4] = new TH1F("jet1_pT_baseline",             "Sub Leading Jet p_{T}- baseline",                                           1000,0,1000);
  h.jet[5] = new TH1F("jet1_Eta_baseline",            "Sub Leading Jet Eta- baseline",                                             64,-3.2,3.2);
  h.jet[6] = new TH1F("jet1_Phi_baseline",            "Sub Leading Jet Phi- baseline",                                             64,-3.2,3.2);
  
  
  h.jet[7] = new TH1F("Ht_baseline",                  "Scalar sum of p_{T} of Jets (H_{t})- baseline",                            1600,0,1600);

  
  h.jet[16] = new TH1F("Ht_j0j1_baseline"  ,          "Scalar sum of p_{T} of Leading and Subleading Jets- baseline",             1600,0,1600);
  h.jet[20] = new TH1F("Ht_50_baseline"  ,            "Scalar sum of p_{T} of Jets (p_{T}>50)- baseline",                         1600,0,1600);
  h.jet[21] = new TH1F("Ht_70_baseline"  ,            "Scalar sum of p_{T} of Jets (p_{T}>70)- baseline",                         1600,0,1600);

  

  h.jet[8] = new TH1F("dphi_jets_baseline",           "Angle between Leading and Subleading Jets- baseline"  ,                    64,-3.2,3.2);
  h.jet[9] = new TH1F("dR_jets_baseline",             "Ang. distance between Leading and Subleading Jets- baseline"  ,            200,0,10);
  h.jet[10] = new TH1F("invar_jet_mass_baseline",     "Invariant Mass of Jets- baseline",                                         1000,0,1000);
  h.jet[11] = new TH1F("dijet_pt_baseline",           "Dijet p_{T}- baseline",                                                    1000,0,1000);

  
  h.jet[12] = new TH1F("dR_mu0j0_baseline",           "Ang. Distance between Leading Muon and Leading Jet- baseline",              200,0,10);
  h.jet[13] = new TH1F("dR_mu0j1_baseline",           "Ang. Distance between Leading Muon and Subleading Jet- baseline",           200,0,10);

  h.jet[14] = new TH1F("dphi_mu0j0_baseline",         "Angle between Leading Lepton and Leading jet- baseline"  ,                  64,-3.2,3.2);
  h.jet[15] = new TH1F("dphi_mu0j1_baseline",         "Angle between Leading Lepton and Subleading jet- baseline"  ,               64,-3.2,3.2);
  h.jet[17] = new TH1F("dphi_J0J1_mu0_baseline",      "Angle between the 2 Jet System and Leading Lepton- baseline",               64,-3.2,3.2);
 
  h.jet[18] =new TH1F("jet_score_sum_baseline",       "Jet Score Sum- baseline",                                                   1000,0,2);
  h.jet[19] =new TH1F("jet_score_product_baseline",   "Jet Score Product- baseline",                                               1000,0,1);

  
  //bjets----------------------------------------
 
  h.b_jet[0]  = new TH1F("n_bjets_baseline",            "Number of b Jets- baseline",                                                10,0,10);
  h.b_jet[1] = new TH1F("bjet0_pT_baseline",                "Leading bJet p_{T}- baseline",                                                    1000,0,1000);
  h.b_jet[2] = new TH1F("bjet0_Eta_baseline",               "Leading bJet Eta- baseline",                                                      64,-3.2,3.2);
  h.b_jet[3] = new TH1F("bjet0_Phi_baseline",               "Leading bJet Phi- baseline",                                                      64,-3.2,3.2);

  h.b_jet[4] = new TH1F("invar_mu0bjet_mass_baseline",       "Invariant Mass of Leading Muon &bJet- baseline",                                 1000,0,1000);
  h.b_jet[5] = new TH1F("dR_mu0bj0_baseline",                "Ang. distance between leading Muon and subleading jet- baseline",                200,0,10);
  h.b_jet[6] = new TH1F("dphi_mu0bj0_baseline",              "Angle between Leading Muon and subleading jet- baseline",                        64,-3.2,3.2);

 
  //FOR MET------------------------------
  h.MET[0]= new TH1F("puppi_met_pT_baseline",                "MET p_{T}- baseline",                                                     1000,0,1000);
  h.MET[1]= new TH1F("puppi_met_phi_baseline",               "MET phi- baseline",                                                       64,-3.2,3.2);
  
  h.MET[2]= new TH1F("mT_mu0_baseline",                      "Transverse Mass of MET and Leading Muon- baseline",                       1000,0,1000);
  h.MET[3]= new TH1F("mT_Jet0_baseline",                     "Transverse Mass of MET and leading Jet- baseline",                        1000,0,1000);
  h.MET[4]= new TH1F("mT_Jet1_baseline",                     "Transverse Mass of MET and sub leading Jet- baseline",                    1000,0,1000);
  h.MET[5]= new TH1F("mT_bJet_baseline",                     "Transverse Mass of MET and b Jet- baseline",                              1000,0,1000);


  h.MET[6] = new TH1F("dphi_mu0MET_baseline",                "Angle between leading Muon and MET- baseline",                             64,-3.2,3.2);
  h.MET[7] = new TH1F("dphi_J0MET_baseline",                 "Angle between leading Jet and MET- baseline",                              64,-3.2,3.2);
  h.MET[8] = new TH1F("dphi_J1MET_baseline",                 "Angle between Subleading Jet and MET- baseline",                           64,-3.2,3.2);
  h.MET[9] = new TH1F("dphi_bJ0MET_baseline",                "Angle between leading bJet and MET- baseline",                             64,-3.2,3.2);

  
  h.MET[10] = new TH1F("dphi_J0J1MET_baseline",              "Angle between the 2 Jet System and MET- baseline",                         64,-3.2,3.2);
  h.MET[11] = new TH1F("mT_J0J1_baseline",                   "Transverse mass of MET and 2 Jet System- baseline",                         1000,0,1000);

  

  //----------------------------------------------------------------------------------------------------------------------------
  //-------------------------------------------------CUTS--------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------------------------------------
  
  h.n1mu2J_cut  = new TH1F("n_1mu2j_cut1" ,           "No of events with 1Mu 2J- cut1" ,                         10,0,10);
  h.evt_wt_cut = new TH1F("evt_wt_cut1",                  "Event Weight- cut1",  2000, 0.85,1);

  
  // FOR LEPTONS---------------------------------------
  /*
    h.elesize =   new TH1F("n_ele_baseline",      "Number of ele- baseline",                            10, 0, 10);//this is after evt selection
    h.lepsize =  new TH1F("n_lep_baseline",      "Number of lep- baseline",                            10, 0, 10);//this is after event selection
    
    h.cutbased_baseline  = new TH1F ("lep_cutbased_baseline",   "Lepton Cutbased- baseline",                            6,-3,3);
  */

  
  //FOR MUONS----------------------------------------------------------
  
  h.mu_cut[0] = new TH1F("nMuon_cut1" ,    "No of Muons- cut1",                                10,0,10);
  h.mu_cut[1] = new TH1F("mu0_pT_cut1",         "Leading Muon p_{T}- cut1",                                   1000,0,1000);
  h.mu_cut[2] = new TH1F("mu0_eta_cut1",       "Leading Muon Eta- cut1",                                        64,-3.2,3.2);
  h.mu_cut[3] = new TH1F("mu0_phi_cut1",        "Leading Muon Phi- cut1",                                        64,-3.2,3.2);

  h.mu_cut[4] = new TH1F("invar_mass_mu0j0_cut1", "Invariant Mass of Leading Muon and Leading Jet- cut1",            1000,0,1000);
  h.mu_cut[5] = new TH1F("invar_mass_mu0j1_cut1",  "Invariant Mass of Leading Muon and Subleading Jet- cut1",       1000,0,1000);


  
  h.mu_cut[6] = new TH1F("Lt_cut1" ,               "Scalar sum of pT of all Muon- cut1" ,                          1600,0,1600);
  h.mu_cut[7] = new TH1F("ST_cut1" ,              "Scalar sum of pT- cut1" ,                                              1600,0,1600);


  h.mu_cut[8] = new TH1F("objects_pz_cut1",            "P_{z} of all the Objects in the EveNt Selection- cut1",     1600,0,1600);


  h.mu_cut[9] = new TH1F("Iso_mu0_cut1" ,               "Muon Isolation- cut1 " ,            600,0,0.2);
  h.mu_cut[10] = new TH1F("mu_dxy_cut1" ,               "Muon dxy- cut1" ,                           1000,-0.06,0.06);
  h.mu_cut[11] = new TH1F("mu_dz_cut1" ,               "Muon dz- cut1" ,                          1000,-0.09,0.09);

  h.mu_cut[12] = new TH1F("mu_sip3d_cut1" ,               "Muon 3D impact parameter significance- cut1" ,                          500,0,10);

  
  //FOR JETS------------------------------------------
 
  
  h.jet_cut[0]  = new TH1F("n_jets_cut1",         "Number of  Jets- cut1 ",                                           10,0,10);
  h.jet_cut[1] = new TH1F("jet0_pT_cut1",             "Leading Jet p_{T}- cut1",                                               1000,0,1000);
  h.jet_cut[2] = new TH1F("jet0_Eta_cut1",            "Leading Jet Eta- cut1",                                                 64,-3.2,3.2);
  h.jet_cut[3] = new TH1F("jet0_Phi_cut1",            "Leading Jet Phi- cut1",                                                 64,-3.2,3.2);

  h.jet_cut[4] = new TH1F("jet1_pT_cut1",             "Sub Leading Jet p_{T}- cut1",                                           1000,0,1000);
  h.jet_cut[5] = new TH1F("jet1_Eta_cut1",            "Sub Leading Jet Eta- cut1",                                             64,-3.2,3.2);
  h.jet_cut[6] = new TH1F("jet1_Phi_cut1",            "Sub Leading Jet Phi- cut1",                                             64,-3.2,3.2);
  
  
  h.jet_cut[7] = new TH1F("Ht_cut1",                  "Scalar sum of p_{T} of Jets (H_{t})- cut1",                            1600,0,1600);

  
  h.jet_cut[16] = new TH1F("Ht_j0j1_cut1"  ,          "Scalar sum of p_{T} of Leading and Subleading Jets- cut1",             1600,0,1600);
  h.jet_cut[20] = new TH1F("Ht_50_cut1"  ,            "Scalar sum of p_{T} of Jets (p_{T}>50) -cut1",                         1600,0,1600);
  h.jet_cut[21] = new TH1F("Ht_70_cut1"  ,            "Scalar sum of p_{T} of Jets (p_{T}>70)- cut1",                         1600,0,1600);

  

  h.jet_cut[8] = new TH1F("dphi_jets_cut1",           "Angle between Leading and Subleading Jets- cut1"  ,                    64,-3.2,3.2);
  h.jet_cut[9] = new TH1F("dR_jets_cut1",             "Ang. distance between Leading and Subleading Jets - cut1"  ,            200,0,10);
  h.jet_cut[10] = new TH1F("invar_jet_mass_cut1",     "Invariant Mass of Jets- cut1",                                         1000,0,1000);
  h.jet_cut[11] = new TH1F("dijet_pt_cut1",           "Dijet p_{T}- cut1",                                                    1000,0,1000);

  
  h.jet_cut[12] = new TH1F("dR_mu0j0_cut1",           "Ang. Distance between Leading Muon and Leading Jet- cut1",              200,0,10);
  h.jet_cut[13] = new TH1F("dR_mu0j1_cut1",           "Ang. Distance between Leading Muon and Subleading Jet- cut1",           200,0,10);

  h.jet_cut[14] = new TH1F("dphi_mu0j0_cut1",         "Angle between Leading Lepton and Leading jet- cut1"  ,                  64,-3.2,3.2);
  h.jet_cut[15] = new TH1F("dphi_mu0j1_cut1",         "Angle between Leading Lepton and Subleading jet- cut1"  ,               64,-3.2,3.2);
  h.jet_cut[17] = new TH1F("dphi_J0J1_mu0_cut1",      "Angle between the 2 Jet System and Leading Lepton- cut1",               64,-3.2,3.2);
 
  h.jet_cut[18] =new TH1F("jet_score_sum_cut1",       "Jet Score Sum- cut1",                                                   1000,0,2);
  h.jet_cut[19] =new TH1F("jet_score_product_cut1",   "Jet Score Product- cut1",                                               1000,0,1);

  
  //bjets----------------------------------------
 
  h.b_jet_cut[0]  = new TH1F("n_bjets_cut1",            "Number of b Jets- cut1",                                                10,0,10);
  h.b_jet_cut[1] = new TH1F("bjet0_pT_cut1",                "Leading bJet p_{T}- cut1",                                                    1000,0,1000);
  h.b_jet_cut[2] = new TH1F("bjet0_Eta_cut1",               "Leading bJet Eta- cut1",                                                      64,-3.2,3.2);
  h.b_jet_cut[3] = new TH1F("bjet0_Phi_cut1",               "Leading bJet Phi- cut1",                                                      64,-3.2,3.2);

  h.b_jet_cut[4] = new TH1F("invar_mu0bjet_mass_cut1",       "Invariant Mass of Leading Muon &bJet- cut1",                                 1000,0,1000);
  h.b_jet_cut[5] = new TH1F("dR_mu0bj0_cut1",                "Ang. distance between leading Muon and subleading jet- cut1",                200,0,10);
  h.b_jet_cut[6] = new TH1F("dphi_mu0bj0_cut1",              "Angle between Leading Muon and subleading jet- cut1",                        64,-3.2,3.2);

 
  //FOR MET------------------------------
  h.MET_cut[0]= new TH1F("puppi_met_pT_cut1",                "MET p_{T}- cut1",                                                     1000,0,1000);
  h.MET_cut[1]= new TH1F("puppi_met_phi_cut1",               "MET phi- cut1",                                                       64,-3.2,3.2);
  
  h.MET_cut[2]= new TH1F("mT_mu0_cut1",                      "Transverse Mass of MET and Leading Muon- cut1",                       1000,0,1000);
  h.MET_cut[3]= new TH1F("mT_Jet0_cut1",                     "Transverse Mass of MET and leading Jet- cut1",                        1000,0,1000);
  h.MET_cut[4]= new TH1F("mT_Jet1_cut1",                     "Transverse Mass of MET and sub leading Jet- cut1",                    1000,0,1000);
  h.MET_cut[5]= new TH1F("mT_bJet_cut1",                     "Transverse Mass of MET and b Jet- cut1",                              1000,0,1000);


  h.MET_cut[6] = new TH1F("dphi_mu0MET_cut1",                "Angle between leading Muon and MET- cut1",                             64,-3.2,3.2);
  h.MET_cut[7] = new TH1F("dphi_J0MET_cut1",                 "Angle between leading Jet and MET- cut1",                              64,-3.2,3.2);
  h.MET_cut[8] = new TH1F("dphi_J1MET_cut1",                 "Angle between Subleading Jet and MET- cut1",                           64,-3.2,3.2);
  h.MET_cut[9] = new TH1F("dphi_bJ0MET_cut1",                "Angle between leading bJet and MET- cut1",                             64,-3.2,3.2);

  
  h.MET_cut[10] = new TH1F("dphi_J0J1MET_cut1",              "Angle between the 2 Jet System and MET- cut1",                         64,-3.2,3.2);
  h.MET_cut[11] = new TH1F("mT_J0J1_cut1",                   "Transverse mass of MET and 2 Jet System- cut1",                         1000,0,1000);


}//book hist
















