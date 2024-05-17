
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
/*
This is a driver script.
It decides which code to run over which sample, the names
of output files and so on.
*/


void ana(int sample=0){
  const char *hstfilename, *sumfilename;
  //Declare a chain for input files.
  TChain *chain = new TChain("Events"); //"Events"
  //Declare an instance of our code class
  SingleTopAna m_selec;
  
  //SingleTop s-channel
  //-------------------------DATA------------------------------------------
  if(sample==0){
    chain->Add("Data/*");
    hstfilename = "hst_singletop_partdata.root";
    sumfilename = "sum_singletop_partdata.txt";
    m_selec.SetData(1); //MC=0, data=1
    m_selec.SetYear(2018);
  }

  
  //-------------------------MC---------------------------------
  if(sample==1){
    //Add one file to chain. This is the input file.
    // chain->Add("samples/schannel/VLL_SingleTop_s-channel_LeptonDecays_10.root");
    chain -> Add("samples/schannel/*");
    
    //Set Names of outputfiles
    hstfilename = "hst_singletop_schannel.root";
    sumfilename = "sum_singletop_schannel.txt";
    //Set some options
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2018);
  }
  
  //SingleTop t-channel 
  
  if(sample==2){
    chain->Add("samples/tchannel_top/*");
    hstfilename = "hst_singletop_tchannel_top.root";
    sumfilename = "sum_singletop_tchannel_top.txt";
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2018);
  }
  
  //Single antiTop t-channel
  
  if(sample==3){
    chain->Add("samples/tchannel_antitop/*");
    hstfilename = "hst_singletop_tchannel_antitop.root";
    sumfilename = "sum_singletop_tchannel_antitop.txt";
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2018);
  }
  
  //Single Top tW-channel
  
  if(sample==4){
    chain->Add("samples/tWchannel_top/*");
    hstfilename = "hst_singletop_tWchannel_top.root";
    sumfilename = "sum_singletop_tWchannel_top.txt";
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2018);
  }
  
  //Single antiTop tW-channel
  
  if(sample==5){
    chain->Add("samples/tWchannel_antitop/*");
    hstfilename = "hst_singletop_tWchannel_antitop.root";
    sumfilename = "sum_singletop_tWchannel_antitop.txt";
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2018);
  }

  //BACKGROUNDS---------------------------

  //ttbar
   if(sample==6){
    chain->Add("backgrounds/ttbar/*");
    hstfilename = "hst_singletop_background_ttbar.root";
    sumfilename = "sum_singletop_background_ttbar.txt";
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2018);
  }

   //wjets
   if(sample==7){
    chain->Add("backgrounds/wjets/*");
    hstfilename = "hst_singletop_background_wjets.root";
    sumfilename = "sum_singletop_background_wjets.txt";
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2018);
  }

    //ww
   if(sample==8){
    chain->Add("backgrounds/ww/*");
    hstfilename = "hst_singletop_background_ww.root";
    sumfilename = "sum_singletop_background_ww.txt";
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2018);
  }
  
  
  std::cout<<"Output files are "<<hstfilename<<" and "<<sumfilename<<std::endl;
  // Set some more options.. set the output file names.
  m_selec.SetHstFileName(hstfilename);
  m_selec.SetSumFileName(sumfilename);
  m_selec.SetVerbose(10);//set verbosity level for output.
  // Call the process function which runs the code.
  chain->Process(&m_selec);
  
}
