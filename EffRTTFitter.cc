#include "EffRTTFitter.hh"
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TLorentzVector.h>

#include <cassert>
#include <sstream>
#include <iomanip>

#include "CPlot.hh"
#include "KStyle.hh"
//#include "ZSignals.hh"
#include "CEffUser1D.hh"


// RooFit headers
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
//#include "RooDataHist.h"
//#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
//#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"


#include "RooAbsPdf.h"
//#include "RooHistPdf.h"




// bin size constants
// #define BIN_SIZE_PASS 20
// #define BIN_SIZE_FAIL 20 //was 20



//--------------------------------------------------------------------------------------------------
EffRTTFitter::EffRTTFitter():
fIsInitialized(false),
fSigPass      (0),
fBkgPass      (0),
fSigFail      (0),
fBkgFail      (0),
fOutputDir    ("."),
topmvacut     (-999.),
maketemplates     (true),
fFitLo (0),
fFitHi (1000),
BIN_SIZE_PASS (20),
BIN_SIZE_FAIL (20),
SIGFR12MASK(0x111),
FR2SF(0.0),
FR2ERR(0.0),
MODEL(1),
VERBOSE(true),
Nsigexp(0.), 
Nbkg1exp(0.),
Nbkg2exp(0.),
nsig_obs(0.),
nbkg1_obs(0.),
nbkg2_obs(0.),
eff_sig_exp(0.),
eff_fr1_exp(0.),
eff_fr2_exp(0.),
eff_sig_obs(0.), 
effl_sig_obs(0.), 
effh_sig_obs(0.),
fr1_sig_obs(0.), 
fr1l_sig_obs(0.), 
fr1h_sig_obs(0.),
fr2_sig_obs(0.), 
fr2l_sig_obs(0.), 
fr2h_sig_obs(0.),
nsigpass_exp(0.), 
nbkg1pass_exp(0.), 
nbkg2pass_exp(0.), 
nsigfail_exp(0.), 
nbkg1fail_exp(0.), 
nbkg2fail_exp(0.), 
npass_data(0.),
nfail_data(0.),
nsigpass_obs(0.), 
nsigpass_err_obs(0.),
nbkg1pass_obs(0.), 
nbkg1pass_err_obs(0.),
nbkg2pass_obs(0.), 
nbkg2pass_err_obs(0.),
nsigfail_obs(0.), 
nsigfail_err_obs(0.),
nbkg1fail_obs(0.), 
nbkg1fail_err_obs(0.),
nbkg2fail_obs(0.), 
nbkg2fail_err_obs(0.),
chisqpass(0.),
chisqfail(0.)

{

  fout = 0;
  
  histfile_signal = 0;
  histfile_bkg1 = 0;  
  histfile_bkg2 = 0;

      
  fPassTree = 0; 
  fFailTree = 0;

  Nbkg1Pass = 0; 
  Nbkg1Fail = 0; 
  Nbkg2Pass = 0; 
  Nbkg2Fail = 0;
  NsigPass  = 0 ;
  NsigFail  = 0 ;
  h_pass_data = 0;
  h_fail_data = 0;
  h_fail_diff_af = 0;
  h_pass_diff_af = 0;

  dataCombined=0;

  sigModPass      = 0;
  bkg1ModPass     = 0;
  bkg2ModPass     = 0;
  sigModFail      = 0;
  bkg1ModFail     = 0;
  bkg2ModFail     = 0;
  sigModPass_rdh  = 0; 
  bkg1ModPass_rdh = 0;
  bkg2ModPass_rdh = 0;  
  sigModFail_rdh  = 0;
  bkg1ModFail_rdh = 0;
  bkg2ModFail_rdh = 0;   
  sigModPass2     = 0;
  bkg1ModPass2    = 0;
  bkg2ModPass2    = 0;
  sigModFail2     = 0; 
  bkg1ModFail2    = 0;
  bkg2ModFail2    = 0;
  modelPass       = 0;
  modelFail       = 0;
  dataPass        = 0;
  dataFail        = 0;
  fitResult       = 0;

  h_pass_data     = 0; 
  h_fail_data     = 0; 
  h_pass_signal   = 0; 
  h_pass_bkg1     = 0;  
  h_pass_bkg2     = 0;  
  h_fail_signal   = 0; 
  h_fail_bkg1     = 0;  
  h_fail_bkg2     = 0;  
  // h_pass_signal_n = 0;  
  // h_pass_bkg1_n   = 0;  
  // h_pass_bkg2_n   = 0; 
  // h_fail_signal_n = 0;  
  // h_fail_bkg1_n   = 0;  
  // h_fail_bkg2_n   = 0; 
  h_pass_signal_af= 0;  
  h_pass_bkg1_af  = 0;  
  h_pass_bkg2_af  = 0;  //after fit (normalized to fitted fractios) 
  h_fail_signal_af= 0; 
  h_fail_bkg1_af  = 0;  
  h_fail_bkg2_af  = 0;  


}

//--------------------------------------------------------------------------------------------------
EffRTTFitter::~EffRTTFitter()
{
  if(fout){
    fout->Close();
    delete fout;
    fout = 0;
  }
  if(fPassTree){
    delete fPassTree;     
    fPassTree=0;  
  }
  if(fFailTree){
    delete fFailTree;     
    fFailTree=0;  
  }

  if(Nbkg1Pass){
    delete Nbkg1Pass; 
    Nbkg1Pass = 0;
  }
  if(Nbkg2Pass){
    delete Nbkg2Pass; 
    Nbkg2Pass = 0;
  }
  if(NsigPass){
    delete NsigPass;
    NsigPass = 0;
  }
  if(Nbkg1Fail){
    delete Nbkg1Fail; 
    Nbkg1Fail = 0 ;
  }
  if(Nbkg2Fail){
    delete Nbkg2Fail; 
    Nbkg2Fail = 0;
  }
  if(NsigFail){
    delete NsigFail;
    NsigFail = 0;
  }
  if(modelPass){
    delete modelPass;
    modelPass = 0;
  }
  if(modelFail){
    delete modelFail;  
    modelFail = 0;
  }
  if(dataCombined){
    delete dataCombined;
    dataCombined = 0;
  }
  if(dataPass){
    delete dataPass;
    dataPass = 0;
  }
  if(dataFail){
    delete dataFail;
    dataFail = 0;
  }
  if(sigModPass_rdh){
    delete sigModPass_rdh;
    sigModPass_rdh = 0;
  }
  if(bkg1ModPass_rdh){
    delete bkg1ModPass_rdh;
    bkg1ModPass_rdh = 0;
  }
  if(bkg2ModPass_rdh){
    delete bkg2ModPass_rdh;
    bkg2ModPass_rdh = 0;
  }
  if(sigModFail_rdh){
    delete sigModFail_rdh;
    sigModFail_rdh = 0;
  }
  if(bkg1ModFail_rdh){
    delete bkg1ModFail_rdh;
    bkg1ModFail_rdh = 0;
  }
  if(bkg2ModFail_rdh){
    delete bkg2ModFail_rdh;
    bkg2ModFail_rdh = 0;
  }
  if(sigModPass){
    delete sigModPass;
    sigModPass = 0;
  }
  if(sigModPass2){
    delete sigModPass2;
    sigModPass2 = 0;
  }
  if(bkg1ModPass){
    delete bkg1ModPass;  
    bkg1ModPass = 0;
  }
  if(bkg2ModPass){
    delete bkg2ModPass;
    bkg2ModPass = 0;
  }
  if(sigModFail){
    delete sigModFail;
    sigModFail = 0;
  }
  if(sigModFail2){
    delete sigModFail2;
    sigModFail2 = 0;
  }
  if(bkg1ModFail){
    delete bkg1ModFail;
    bkg1ModFail = 0;
  }
  if(bkg2ModFail){
    delete bkg2ModFail;
    bkg2ModFail = 0;
  }

  if(histfile_signal){
    delete histfile_signal;  
    histfile_signal = 0;
  }
  if(histfile_bkg1){
    delete histfile_bkg1;  
    histfile_bkg1 = 0;
  }
  if(histfile_bkg2){
    delete histfile_bkg2;
    histfile_bkg2 = 0;
  }  
}

//--------------------------------------------------------------------------------------------------
void EffRTTFitter::initialize(const std::string infname, const std::string outdir, const std::string temfname1, const std::string temfname2, const std::string temfname3, const double TOPMVACUT, const bool MAKETEMPLATES, const float FITLO, const float FITHI, const float BINSIZEPASS, const float BINSIZEFAIL, const unsigned int SIGNALFAKERATE12MASK, const float FR2SCALEFACTOR, const float FR2ERROR)
{
  std::cout << "   [EffRTTFitter] Initializing... " << std::endl;


  fFitLo = FITLO;
  fFitHi = FITHI;
  topmvacut = TOPMVACUT;
  maketemplates = MAKETEMPLATES;
  BIN_SIZE_PASS = BINSIZEPASS;
  BIN_SIZE_FAIL = BINSIZEFAIL;
  SIGFR12MASK = SIGNALFAKERATE12MASK;
  FR2SF = FR2SCALEFACTOR;
  FR2ERR = FR2ERROR;


  // set up output directory
  fOutputDir = outdir;
  gSystem->mkdir(fOutputDir.c_str(),true);

  if(VERBOSE)
    std::cout<<"SIGFR12MASK "<<SIGFR12MASK<<std::endl;


  //------------------------------------------------------------------------------------------------
  // MAKE TEMPLATES (MC)
  //===============================================================================================
  if(maketemplates){
    if((SIGFR12MASK & 0x001) == 0x001)
      makeBinnedTemplates(temfname1, "signal");
    if((SIGFR12MASK & 0x010) == 0x010)
      makeBinnedTemplates(temfname2, "bkg1");
    if((SIGFR12MASK & 0x100) == 0x100)
      makeBinnedTemplates(temfname3, "bkg2");
  }

  
  //------------------------------------------------------------------------------------------------
  // RETRIEVE INFO FROM INPTUS TO MAKE PASS/FAIL TREES (DATA)
  //===============================================================================================

  //unsigned int runNum, lumiSec, evtNum;   // event ID
  float        evtWeight; // // LUMI*scale1fb*PU*btgag*lep ID*Trig turnon*others*kfact*topPt*NORMTODATA
  float topmass;
  float res_topmva;

  float varforfit; //will be computed on the fly


  
  TFile *infile = new TFile(infname.c_str());    assert(infile);
  TTree *intree = (TTree*)infile->Get("Events"); assert(intree);

  //from Reformat..TnP.C
  // intree->SetBranchAddress("runNum",    &runNum);
  // intree->SetBranchAddress("lumiSec",   &lumiSec);
  // intree->SetBranchAddress("evtNum",    &evtNum);
  intree->SetBranchAddress("evtWeight", &evtWeight);
  intree->SetBranchAddress("topmass", &topmass);
  intree->SetBranchAddress("res_topmva",   &res_topmva);


  unsigned int pass;
  
  char tname[50];
  float wgt;

  {
    sprintf(tname,"pass");
    fPassTree = new TTree(tname,"");
    fPassTree->Branch("m",&varforfit,"m/F");
    fPassTree->Branch("w",&wgt, "w/F");
    fPassTree->SetDirectory(0);
    sprintf(tname,"fail");
    fFailTree = new TTree(tname,"");
    fFailTree->Branch("m",&varforfit,"m/F");
    fFailTree->Branch("w",&wgt, "w/F");
    fFailTree->SetDirectory(0);
  }


  //
  // loop over probes
  //
  for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    
    varforfit = topmass;
    

    //FUTHER SELECTION
    if(varforfit < fFitLo) continue;
    if(varforfit > fFitHi) continue;

    
    if(res_topmva > topmvacut) pass = 1;
    else pass = 0;



    wgt = 1; //scale1fb;
    

    //
    // Fill trees
    //
    if(pass) {
      fPassTree->Fill();
    
    } else {
      fFailTree->Fill();
      
    }
  }
  delete infile;
  infile=0, intree=0;
  
  
  fIsInitialized = true;
}

//--------------------------------------------------------------------------------------------------
void EffRTTFitter::computeEff()
{
  assert(fIsInitialized);
  
  std::cout << "   [EffRTTFitter] Computing efficiencies..." << std::endl;


  //------------------------------------------------------------------------------------------------
  // Efficiency calculation
  //================================================================================================
  CreateOutFile(false); 
  
  //  double eff, errl, errh;  
  // performFit(eff, errl, errh,
  // 	     fPassTree, fFailTree); 
  performFit(); 

  //Store templates in out.root for later fun
  WriteHistos(false);

  //Print fit results straight from Roofit
  printRooFitresults(false); 

  //Print efficiencies, SF, etc..
  printresults(true,true);//print to file
  printresults(false,false); //print to screen

  //DO NOT PUT PRINTOUTS HERE-> it will crash. NOT SURE WHY
  //std::cout<<"---"<<std::endl;


}


//--------------------------------------------------------------------------------------------------
void EffRTTFitter::makeBinnedTemplates(const std::string temfname, string name)
{
  std::cout << "   [EffRTTFitter] Creating binned templates... "; std::cout.flush();
  char hname[50];
  
  
  TH1D* h_pass;
  TH1D* h_fail;

  sprintf(hname,"h_pass_%s",name.c_str());
  h_pass = new TH1D(hname,"",int(fFitHi-fFitLo)/BIN_SIZE_PASS,fFitLo,fFitHi);
  h_pass->SetDirectory(0);
  sprintf(hname,"h_fail_%s",name.c_str());
  h_fail = new TH1D(hname,"",int(fFitHi-fFitLo)/BIN_SIZE_FAIL,fFitLo,fFitHi);
  h_fail->SetDirectory(0);

  
    
  //------------------------------------------------------------------------------------------------
  // RETRIEVE INFO FROM INPTUS TO MAKE PASS/FAIL TREES (MC)
  //===============================================================================================

  //unsigned int runNum, lumiSec, evtNum;   // event ID
  float        evtWeight; // // LUMI*scale1fb*PU*btgag*lep ID*Trig turnon*others*kfact*topPt*NORMTODATA
  float topmass;
  float res_topmva;

  float varforfit; //will be computed on the fly


  
  TFile *infile = new TFile(temfname.c_str());    assert(infile);
  TTree *intree = (TTree*)infile->Get("Events"); assert(intree);

  //from Reformat..TnP.C
  // intree->SetBranchAddress("runNum",    &runNum);
  // intree->SetBranchAddress("lumiSec",   &lumiSec);
  // intree->SetBranchAddress("evtNum",    &evtNum);
  intree->SetBranchAddress("evtWeight", &evtWeight);
  intree->SetBranchAddress("topmass", &topmass);
  intree->SetBranchAddress("res_topmva",   &res_topmva);

  

  unsigned int pass;
  
  for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);


    varforfit = topmass;

    double weight = evtWeight; //evtWeight contains all weights (NORMTODATA included) 
    

    if(res_topmva > topmvacut) pass = 1;
    else pass = 0;
    


    if(pass) { 
      h_pass->Fill(varforfit,weight);
    } else {
      h_fail->Fill(varforfit,weight);
    }    
  }
  infile->Close();
 
  string str_tmp = fOutputDir+"/binnedTemplates_"+name+".root";
  TFile outfile(str_tmp.c_str(), "RECREATE");

  h_pass->Write();
  h_fail->Write();

  if(VERBOSE){
    cout<<"Yields"<<endl;
    cout<<"Pass "<<h_pass->Integral()<<" Fail "<<h_fail->Integral()<<" Pass+Fail "<<h_pass->Integral()+h_fail->Integral()<<endl;
  }

  delete h_pass;
  delete h_fail;

  outfile.Write();
  outfile.Close(); 

  cout << "Done!" << endl;
}




//--------------------------------------------------------------------------------------------------
//void EffRTTFitter::performFit(double &resEff, double &resErrl, double &resErrh,
//                           TTree *passTree, TTree *failTree)
void EffRTTFitter::performFit()

{
  // string str_fname = fOutputDir+"/out.root";
  // TFile *fout = new TFile(str_fname.c_str(),"RECREATE");



  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  std::cout << " ...performFit..." << std::endl; 
  
  string str_title;
  str_title = "Top Mass [GeV/c^2]";
  RooRealVar m("m",str_title.c_str(),fFitLo,fFitHi);

  
  m.setBins(10000);
  

  //retrieve templates from input files
  // TFile *histfile_signal = 0;
  // TFile *histfile_bkg1 = 0;  TFile *histfile_bkg2 = 0;

  histfile_signal = TFile::Open(((TString)fOutputDir)+"/binnedTemplates_signal.root");
  if((SIGFR12MASK & 0x001) == 0x001) assert(histfile_signal);
  histfile_bkg1 = new TFile(((TString)fOutputDir)+"/binnedTemplates_bkg1.root");  
  if((SIGFR12MASK & 0x010) == 0x010) assert(histfile_bkg1);
  histfile_bkg2 = new TFile(((TString)fOutputDir)+"/binnedTemplates_bkg2.root");  
  if((SIGFR12MASK & 0x100) == 0x100) assert(histfile_bkg2);

  //pass
  if((SIGFR12MASK & 0x001) == 0x001){
    char hname_signal[50];    
    sprintf(hname_signal,"h_pass_%s","signal");
    h_pass_signal = (TH1D*)histfile_signal->Get(hname_signal);
    cout<<"looking for histo "<<hname_signal<<" in file "<<histfile_signal->GetName()<<std::endl;
    assert(h_pass_signal);
  }
  if((SIGFR12MASK & 0x010) == 0x010){
    char hname_bkg1[50];
    sprintf(hname_bkg1,"h_pass_%s","bkg1");
    h_pass_bkg1 = (TH1D*)histfile_bkg1->Get(hname_bkg1);
    assert(h_pass_bkg1);
  }
  if((SIGFR12MASK & 0x100) == 0x100){
    char hname_bkg2[50];
    sprintf(hname_bkg2,"h_pass_%s","bkg2");
    h_pass_bkg2 = (TH1D*)histfile_bkg2->Get(hname_bkg2);
    assert(h_pass_bkg2);
  }

  //fail
  if((SIGFR12MASK & 0x001) == 0x001){
    char hname_signal[50];
    sprintf(hname_signal,"h_fail_%s","signal");
    h_fail_signal = (TH1D*)histfile_signal->Get(hname_signal);
    assert(h_fail_signal);
  }
  if((SIGFR12MASK & 0x010) == 0x010){
    char hname_bkg1[50];
    sprintf(hname_bkg1,"h_fail_%s","bkg1");
    h_fail_bkg1 = (TH1D*)histfile_bkg1->Get(hname_bkg1);
    assert(h_fail_bkg1);
  }
  if((SIGFR12MASK & 0x100) == 0x100){
    char hname_bkg2[50];
    sprintf(hname_bkg2,"h_fail_%s","bkg2");
    h_fail_bkg2 = (TH1D*)histfile_bkg2->Get(hname_bkg2);
    assert(h_fail_bkg2);
  }


  
  // Define Pass/Fail categories in Data
  RooCategory sample("sample","");
  sample.defineType("Pass",1);
  sample.defineType("Fail",2);
  
  TH1D histPass("histPass","",int(fFitHi-fFitLo)/BIN_SIZE_PASS,fFitLo,fFitHi); 
  TH1D histFail("histFail","",int(fFitHi-fFitLo)/BIN_SIZE_FAIL,fFitLo,fFitHi);
  //  RooAbsData *dataCombined=0;
  
  std::cout<<"passtree entries "<<fPassTree->Draw("m>>histPass","w")<<std::endl;
  std::cout<<"failtree entries "<<fFailTree->Draw("m>>histFail","w")<<std::endl;
  dataPass = new RooDataHist("dataPass","dataPass",RooArgSet(m),&histPass);
  dataFail = new RooDataHist("dataFail","dataFail",RooArgSet(m),&histFail);
  //m.setBins(100);  

  dataCombined = new RooDataHist("dataCombined","dataCombined",RooArgList(m),
				 RooFit::Index(sample),
				 RooFit::Import("Pass",*((RooDataHist*)dataPass)),
				 RooFit::Import("Fail",*((RooDataHist*)dataFail)));
  
  


  // Make Pass PDF from MC histograms
  if((SIGFR12MASK & 0x001) == 0x001){
    sigModPass_rdh= new RooDataHist("sigModPass_rdh","sigModPass_rdh",RooArgSet(m),h_pass_signal);
    sigModPass = new RooHistPdf("sigModPass","sigModPass",RooArgSet(m),*sigModPass_rdh);
    sigModPass2 = new CMCTemplateConvGaussian(m,h_pass_signal,0);

    
  }
  if((SIGFR12MASK & 0x010) == 0x010){
    bkg1ModPass_rdh = new RooDataHist("bkg1ModPass_rdh","bkg1ModPass_rdh",RooArgSet(m),h_pass_bkg1);
    bkg1ModPass = new RooHistPdf("bkg1ModPass","bkg1ModPass",RooArgSet(m),*bkg1ModPass_rdh);
    bkg1ModPass2 = new CMCTemplateConvGaussian(m,h_pass_bkg1,1);
  }
  if((SIGFR12MASK & 0x100) == 0x100){
    bkg2ModPass_rdh = new RooDataHist("bkg2ModPass_rdh","bkg2ModPass_rdh",RooArgSet(m),h_pass_bkg2);
    bkg2ModPass = new RooHistPdf("bkg2ModPass","bkg2ModPass",RooArgSet(m),*bkg2ModPass_rdh);
    bkg2ModPass2 = new CMCTemplateConvGaussian(m,h_pass_bkg2,2);
  }
  

  if(VERBOSE){
    if(sigModPass)
      cout<<"sigModPass "<<sigModPass<<endl;
    if(sigModPass2)
      cout<<"sigModPass2 "<<sigModPass2<<endl;
    if(bkg1ModPass)
      cout<<"bkg1ModPass "<<bkg1ModPass<<endl;
    if(bkg1ModPass2)
      cout<<"bkg1ModPass2 "<<bkg1ModPass2<<endl;
    if(bkg2ModPass)
      cout<<"bkg2ModPass "<<bkg2ModPass<<endl;
    if(bkg2ModPass2)
      cout<<"bkg2ModPass2 "<<bkg2ModPass2<<endl;

    cout<<"sigModPass_rdh "<<sigModPass_rdh<<endl;
    cout<<"bkg1ModPass_rdh "<<bkg1ModPass_rdh<<endl;
    cout<<"bkg2ModPass_rdh "<<bkg2ModPass_rdh<<endl;

  }

  // Make Fail PDF from MC histograms      
  if((SIGFR12MASK & 0x001) == 0x001){
    sigModFail_rdh = new RooDataHist("sigModFail_rdh","sigModFail_rdh",RooArgSet(m),h_fail_signal);
    sigModFail = new RooHistPdf("sigModFail","sigModFail",RooArgSet(m),*sigModFail_rdh);
    sigModFail2 = new CMCTemplateConvGaussian(m,h_fail_signal,0);
  }
  if((SIGFR12MASK & 0x010) == 0x010){
    bkg1ModFail_rdh = new RooDataHist("bkg1ModFail_rdh","bkg1ModFail_rdh",RooArgSet(m),h_fail_bkg1);
    bkg1ModFail = new RooHistPdf("bkg1ModFail","bkg1ModFail",RooArgSet(m),*bkg1ModFail_rdh);
    bkg1ModFail2 = new CMCTemplateConvGaussian(m,h_fail_bkg1,1);
  }
  if((SIGFR12MASK & 0x100) == 0x100){
    bkg2ModFail_rdh = new RooDataHist("bkg2ModFail_rdh","bkg2ModFail_rdh",RooArgSet(m),h_fail_bkg2);
    bkg2ModFail = new RooHistPdf("bkg2ModFail","bkg2ModFail",RooArgSet(m),*bkg2ModFail_rdh);
    bkg2ModFail2 = new CMCTemplateConvGaussian(m,h_fail_bkg2,2);
  }

  if(VERBOSE){
    if(sigModFail)
      cout<<"sigModFail "<<sigModFail<<endl;
    if(sigModFail2)
      cout<<"sigModFail2 "<<sigModFail2<<endl;
    if(bkg1ModFail)
      cout<<"bkg1ModFail "<<bkg1ModFail<<endl;
    if(bkg1ModFail2)
      cout<<"bkg1ModFail2 "<<bkg1ModFail2<<endl;
    if(bkg2ModFail)
      cout<<"bkg2ModFail "<<bkg2ModFail<<endl;
    if(bkg2ModFail2)
      cout<<"bkg2ModFail2 "<<bkg2ModFail2<<endl;
    

    cout<<"sigModFail_rdh "<<sigModFail_rdh<<endl;
    cout<<"bkg1ModFail_rdh "<<bkg1ModFail_rdh<<endl;
    cout<<"bkg2ModFail_rdh "<<bkg2ModFail_rdh<<endl;

  }


  //expected yields
  if(h_pass_signal)  nsigpass_exp = h_pass_signal->Integral();
  if(h_pass_bkg1)    nbkg1pass_exp = h_pass_bkg1->Integral();
  if(h_pass_bkg2)    nbkg2pass_exp = h_pass_bkg2->Integral();

  if(h_fail_signal)  nsigfail_exp = h_fail_signal->Integral();
  if(h_fail_bkg1)    nbkg1fail_exp = h_fail_bkg1->Integral();
  if(h_fail_bkg2)    nbkg2fail_exp = h_fail_bkg2->Integral();

  Nsigexp = nsigpass_exp+nsigfail_exp;
  Nbkg1exp = nbkg1pass_exp+nbkg1fail_exp;
  Nbkg2exp = nbkg2pass_exp+nbkg2fail_exp;


  //observed pass/fail yields
  npass_data = histPass.Integral();
  nfail_data = histFail.Integral();

  //Expected Efficiencies
  if((SIGFR12MASK & 0x001) == 0x001 && (nsigpass_exp+nsigfail_exp)>0)
    eff_sig_exp = nsigpass_exp/(nsigpass_exp+nsigfail_exp);
  if((SIGFR12MASK & 0x010) == 0x010 && (nbkg1pass_exp+nbkg1fail_exp) > 0)
    eff_fr1_exp = nbkg1pass_exp/(nbkg1pass_exp+nbkg1fail_exp);
  if((SIGFR12MASK & 0x100) == 0x100 && (nbkg2pass_exp+nbkg2fail_exp) > 0)
    eff_fr2_exp = nbkg2pass_exp/(nbkg2pass_exp+nbkg2fail_exp);


  ////////////////////
  //Get ready to Fit
  ///////////////////


  // Define nuisance parameters
  RooRealVar Nsig("Nsig","Signal Yield",Nsigexp*0.9,0,Nsigexp*1.1);
  RooRealVar Nbkg1("Nbkg1","Bkg1 Yield",Nbkg1exp*0.9,0,Nbkg1exp*1.1);
  RooRealVar Nbkg2("Nbkg2","Bkg2 Yield",Nbkg2exp*0.9,0,Nbkg2exp*1.1);

  RooRealVar eff("eff","Efficiency",eff_sig_exp,eff_sig_exp/2,1);
  RooRealVar FR1("FR1","FR Bkg1",eff_fr1_exp,eff_fr1_exp/2,1);
  double fr2min = eff_fr2_exp*FR2SF - FR2ERR;
  double fr2max = eff_fr2_exp*FR2SF + FR2ERR;
  std::cout<<"&&&&&&&& Fr2exp*FR2SF, fr2min, fr2max, "<<eff_fr2_exp*FR2SF<<" "<<fr2min<<" "<<fr2max<<std::endl;
  RooRealVar FR2("FR2","FR Bkg2",eff_fr2_exp*FR2SF,fr2min,fr2max);

  
  // Define actual fit parameters (function of the above nuisances)
  if((SIGFR12MASK & 0x001) == 0x001)
    NsigPass = new RooFormulaVar("NsigPass","eff*Nsig",RooArgList(eff,Nsig));
  if((SIGFR12MASK & 0x010) == 0x010)
    Nbkg1Pass = new RooFormulaVar("Nbkg1Pass","FR1*Nbkg1",RooArgList(FR1,Nbkg1)); 
  //  Nbkg2Pass = new RooFormulaVar("Nbkg2Pass","FR2overFR1*FR1*Nbkg2",RooArgList(FR2overFR1,FR1,Nbkg2)); 
  if((SIGFR12MASK & 0x100) == 0x100)
    Nbkg2Pass = new RooFormulaVar("Nbkg2Pass","FR2*Nbkg2",RooArgList(FR2,Nbkg2)); 

  if((SIGFR12MASK & 0x001) == 0x001)
    NsigFail = new RooFormulaVar("NsigFail","(1.0-eff)*Nsig",RooArgList(eff,Nsig));
  if((SIGFR12MASK & 0x010) == 0x010)
    Nbkg1Fail = new RooFormulaVar("Nbkg1Fail","(1.0-FR1)*Nbkg1",RooArgList(FR1,Nbkg1));
  //    Nbkg2Fail = new RooFormulaVar("Nbkg2Fail","(1.0-FR2overFR1*FR1)*Nbkg2",RooArgList(FR2overFR1,FR1,Nbkg2));
  if((SIGFR12MASK & 0x100) == 0x100)
    Nbkg2Fail = new RooFormulaVar("Nbkg2Fail","(1.0-FR2)*Nbkg2",RooArgList(FR2,Nbkg2));

  // Define models
  if(NsigPass  && Nbkg1Pass  && Nbkg2Pass){
    if(MODEL==1 && sigModPass && bkg1ModPass && bkg2ModPass)
      modelPass = new RooAddPdf("modelPass","Model for PASS sample",RooArgList(*sigModPass,*bkg1ModPass,*bkg2ModPass),RooArgList(*NsigPass,*Nbkg1Pass,*Nbkg2Pass));
    else if(MODEL==2 && sigModPass2 && bkg1ModPass2 && bkg2ModPass2){
      modelPass = new RooAddPdf("modelPass","Model for PASS sample",RooArgList(*(sigModPass2->model),*(bkg1ModPass2->model),*(bkg2ModPass2->model)),RooArgList(*NsigPass,*Nbkg1Pass,*Nbkg2Pass));
    }
  }
  else if (bkg2ModPass && Nbkg2Pass)
    modelPass = new RooAddPdf("modelPass","Model for PASS sample",RooArgList(*bkg2ModPass),RooArgList(*Nbkg2Pass));

  if(NsigFail && Nbkg1Fail&& Nbkg2Fail){
    if(MODEL==1 && sigModFail && bkg1ModFail && bkg2ModFail )
      modelFail = new RooAddPdf("modelFail","Model for FAIL sample",RooArgList(*sigModFail,*bkg1ModFail,*bkg2ModFail),RooArgList(*NsigFail,*Nbkg1Fail,*Nbkg2Fail)); 
    else if(MODEL==2 && sigModFail2 && bkg1ModFail2 && bkg2ModFail2 ){
      modelFail = new RooAddPdf("modelFail","Model for FAIL sample",RooArgList(*(sigModFail2->model),*(bkg1ModFail2->model),*(bkg2ModFail2->model)),RooArgList(*NsigFail,*Nbkg1Fail,*Nbkg2Fail)); 

    }
  }
  else if (bkg2ModFail && Nbkg2Fail)
    modelFail = new RooAddPdf("modelFail","Model for PASS sample",RooArgList(*bkg2ModFail),RooArgList(*Nbkg2Fail));


  if(VERBOSE){  
    cout<<"*****before setting totalPdf"<<endl;
    cout<<"modelPass "<<modelPass<<endl;
    cout<<"modelFail "<<modelFail<<endl;
  }

  //Add models to the overall PDF  
  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(*modelPass,"Pass"); 
  totalPdf.addPdf(*modelFail,"Fail");

  //keep a copy of templates to be normalized to fit results
  CloneTemplates(true);

  SetTemplatesColors(h_pass_bkg1, h_pass_bkg2, h_pass_signal, false);
  SetTemplatesColors(h_fail_bkg1, h_fail_bkg2, h_fail_signal, false);
  SetTemplatesColors(h_pass_bkg1_af, h_pass_bkg2_af, h_pass_signal_af, false);
  SetTemplatesColors(h_fail_bkg1_af, h_fail_bkg2_af, h_fail_signal_af, false);



  // Plot the imported histogram(s)                                                                               
  PlotDataVsMC(h_pass_bkg1, h_pass_bkg2, h_pass_signal, histPass, false, "TemplatesandData_Passsample_beforefit_mt", true);   
  PlotDataVsMC(h_fail_bkg1, h_fail_bkg2, h_fail_signal, histFail, false,  "TemplatesandData_Failsample_beforefit_mt", false);


  TLine * a = new TLine();
  a->SetLineColor(2);
  a->SetLineWidth(2);
  a->SetLineStyle(1);

  TLegend *legend = new TLegend(0.65,0.54,0.99,0.99);
  legend->SetFillColor(0);
  legend->AddEntry(&histPass,"Data" , "lp");
  if(h_pass_signal)
    legend->AddEntry(h_pass_signal, "Signal" , "lf");
  if(h_pass_bkg2)
    legend->AddEntry(h_pass_bkg2, "Non-tt1l BG" , "lf");
  if(h_pass_bkg1)
    legend->AddEntry(h_pass_bkg1, "tt1l Combinatorial BG" , "lf");



  RooPlotDataVsMC(m, legend, "TemplatesandData_Passsample_beforefit.png","TemplatesandData_Failsample_beforefit.png", true);   






  //MT Mar 7 : commented -> need to understand it
  // std::cout<<"plotting smeared V non smeared pass templates"<<std::endl;
  // RooPlot* mframePass_smearnotsmear_bf = m.frame(Bins(int(fFitHi-fFitLo)/BIN_SIZE_PASS));
  // mframePass_smearnotsmear_bf->SetTitle("");
  // dataPass->plotOn(mframePass_smearnotsmear_bf,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));//needed for the normalization
  // if(MODEL==1)
  //   modelPass->plotOn(mframePass_smearnotsmear_bf,Components("sigModPass"),LineStyle(1),LineColor(kGreen));
  // else if(MODEL==2){
  //   modelPass->plotOn(mframePass_smearnotsmear_bf,Components("signalPass"),LineStyle(1),LineColor(kGreen));
  // }
  // TCanvas* cpass_smearnotsmear_bf = new TCanvas("cpass_smearnotsmear_bf","cpass_smearnotsmear_bf",800,1200);
  // cpass_smearnotsmear_bf->SetWindowPosition(cpass_smearnotsmear_bf->GetWindowTopX()+cpass_smearnotsmear_bf->GetBorderSize()+800,0);
  // cpass_smearnotsmear_bf->SetTickx(1);
  // cpass_smearnotsmear_bf->SetTicky(1);
  // //    mframePass_smearnotsmear_bf->Draw();
  // h_pass_signal->SetFillColor(0);
  // mframePass_smearnotsmear_bf->Draw();
  // h_pass_signal->Draw("same");



  // //std::cout<<(sigModPass2->model)->evaluate()<<std::cout;



  // str_tmp = fOutputDir+"/TemplatesandData_Passsample_smearnotsmear_beforefit.png"; 
  // cpass_smearnotsmear_bf->SaveAs(str_tmp.c_str());

  // std::cout<<"plotting smeared V non smeared fail templates"<<std::endl;
  // RooPlot* mframeFail_smearnotsmear_bf = m.frame(Bins(int(fFitHi-fFitLo)/BIN_SIZE_FAIL));
  // mframeFail_smearnotsmear_bf->SetTitle("");
  // dataFail->plotOn(mframeFail_smearnotsmear_bf,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));//needed for the normalization
  // if(MODEL==1)
  //   modelFail->plotOn(mframeFail_smearnotsmear_bf,Components("sigModFail"),LineStyle(1),LineColor(kGreen));
  // else if(MODEL==2)
  //   modelFail->plotOn(mframeFail_smearnotsmear_bf,Components("signalFail"),LineStyle(1),LineColor(kBlue));
  // TCanvas* cfail_smearnotsmear_bf = new TCanvas("cfail_smearnotsmear_bf","cfail_smearnotsmear_bf",800,1200);
  // cfail_smearnotsmear_bf->SetWindowPosition(cfail_smearnotsmear_bf->GetWindowTopX()+cfail_smearnotsmear_bf->GetBorderSize()+800,0);
  // cfail_smearnotsmear_bf->SetTickx(1);
  // cfail_smearnotsmear_bf->SetTicky(1);

    
  // h_fail_signal->SetFillColor(0);
  // mframeFail_smearnotsmear_bf->Draw();
  // h_fail_signal->Draw("same");
    
  // //      legend->Draw();
    
    
  // // if(MODEL==2){
  // // 	RooRealVar *meanfail = (RooRealVar*)((sigModPass2->model)->mean);
  // // 	std::cout<<"mean "<<meanfail->Eval();
  // // }
    
  // str_tmp = fOutputDir+"/TemplatesandData_Failsample_smearnotsmear_beforefit.png"; 
  // cfail_smearnotsmear_bf->SaveAs(str_tmp.c_str());
  
    

  //  }

  if(VERBOSE)
    cout<<"---Fit---"<<endl;


  ////////////////  
  //Fit
  ///////////////
  fitResult = totalPdf.fitTo(*dataCombined,
                             RooFit::Extended(),
                             //RooFit::Strategy(2),
                             //RooFit::Minos(RooArgSet(eff)),
                             RooFit::NumCPU(4),
                             RooFit::Save());
  
  if(VERBOSE)
    cout<<"---Done Fitting: retrieving results---"<<endl;



  // Refit w/o MINOS if MINOS errors are strange...
  if((fabs(eff.getErrorLo())<5e-5) || (eff.getErrorHi()<5e-5))
    fitResult = totalPdf.fitTo(*dataCombined, RooFit::Extended(), RooFit::Strategy(1), RooFit::Save());


  //Observed Efficiencies
  if((SIGFR12MASK & 0x001) == 0x001){
    eff_sig_obs  = eff.getVal(); 
    effl_sig_obs = fabs(eff.getErrorLo());
    effh_sig_obs = eff.getErrorHi();
  }
  fr1_sig_obs = FR1.getVal();
  fr1l_sig_obs = FR1.getErrorLo();
  fr1h_sig_obs = FR1.getErrorHi();
  fr2_sig_obs = FR2.getVal();
  fr2l_sig_obs = FR2.getErrorLo();
  fr2h_sig_obs = FR2.getErrorHi();
  nsig_obs = Nsig.getVal();
  nbkg1_obs = Nbkg1.getVal();
  nbkg2_obs = Nbkg2.getVal();

  if(NsigPass){  
    nsigpass_obs = NsigPass->getVal();
    nsigpass_err_obs = NsigPass->getPropagatedError(*fitResult);
  }
  if(Nbkg1Pass){ 
    nbkg1pass_obs = Nbkg1Pass->getVal();
    nbkg1pass_err_obs = Nbkg1Pass->getPropagatedError(*fitResult);
  }
  if(Nbkg2Pass){ 
    nbkg2pass_obs = Nbkg2Pass->getVal();
    nbkg2pass_err_obs = Nbkg2Pass->getPropagatedError(*fitResult);
  }

  if(NsigFail){
    nsigfail_obs = NsigFail->getVal();
    nsigfail_err_obs = NsigFail->getPropagatedError(*fitResult);
  }
  if(Nbkg1Fail){
    nbkg1fail_obs = Nbkg1Fail->getVal();
    nbkg1fail_err_obs = Nbkg1Fail->getPropagatedError(*fitResult);    
  }
  if(Nbkg2Fail){
    nbkg2fail_obs = Nbkg2Fail->getVal();
    nbkg2fail_err_obs = Nbkg2Fail->getPropagatedError(*fitResult);
  }




  if(VERBOSE)
    cout<<"---Scale Templates according to Fit Results---"<<endl;

  //Scale template to fitted values
  if(h_pass_signal_af)
    h_pass_signal_af->Scale(NsigPass->getVal()/h_pass_signal_af->Integral());
  if(h_pass_bkg1_af)
    h_pass_bkg1_af->Scale(Nbkg1Pass->getVal()/h_pass_bkg1_af->Integral());
  if(h_pass_bkg2_af)
    h_pass_bkg2_af->Scale(Nbkg2Pass->getVal()/h_pass_bkg2_af->Integral());
  if(h_fail_signal_af)
    h_fail_signal_af->Scale(NsigFail->getVal()/h_fail_signal_af->Integral());
  if(h_fail_bkg1_af)
    h_fail_bkg1_af->Scale(Nbkg1Fail->getVal()/h_fail_bkg1_af->Integral());
  if(h_fail_bkg2_af)
    h_fail_bkg2_af->Scale(Nbkg2Fail->getVal()/h_fail_bkg2_af->Integral());

  //Calc postfit chisq
  chisqpass = CalcChisq(h_pass_bkg1_af,h_pass_bkg2_af,h_pass_signal_af,histPass,false);
  chisqfail = CalcChisq(h_fail_bkg1_af,h_fail_bkg2_af,h_fail_signal_af,histFail,false);


  // Plot the scaled histogram(s) Vs data
  PlotDataVsMC(h_pass_bkg1_af, h_pass_bkg2_af, h_pass_signal_af, histPass, false, "TemplatesandData_Passsample_afterfit_mt", true);   
  PlotDataVsMC(h_fail_bkg1_af, h_fail_bkg2_af, h_fail_signal_af, histFail, false, "TemplatesandData_Failsample_afterfit_mt", false);




  //plot postfit distributions Vs data via Roofit methods (TO BE INSERTED IN A FUCNTION)
  if(MODEL==1){
  RooPlot *mframePass = m.frame(Bins(int(fFitHi-fFitLo)/BIN_SIZE_PASS));
  mframePass->SetTitle("");
  dataPass->plotOn(mframePass,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));    

  modelPass->plotOn(mframePass);
  if(MODEL==1 && sigModPass)
    modelPass->plotOn(mframePass,Components("sigModPass"),LineStyle(1),LineColor(kGreen));
  if(MODEL==2 && sigModPass2)
    modelPass->plotOn(mframePass,Components("signalPass"),LineStyle(1),LineColor(kGreen));
  if(bkg1ModPass)
    modelPass->plotOn(mframePass,Components("bkg1ModPass"),LineStyle(1),LineColor(kRed));
  if(bkg2ModPass)
    modelPass->plotOn(mframePass,Components("bkg2ModPass"),LineStyle(1),LineColor(kMagenta));

  
  RooPlot *mframeFail; 

  mframeFail = m.frame(Bins(int(fFitHi-fFitLo)/BIN_SIZE_FAIL));
  mframeFail->SetTitle("");
  dataFail->plotOn(mframeFail,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
  modelFail->plotOn(mframeFail);
  if(MODEL==1 && sigModFail)
    modelFail->plotOn(mframeFail,Components("sigModFail"),LineStyle(1),LineColor(kGreen));
  if(MODEL==2 && sigModFail2)
    modelFail->plotOn(mframeFail,Components("signalFail"),LineStyle(1),LineColor(kGreen));
  if(bkg1ModFail)
    modelFail->plotOn(mframeFail,Components("bkg1ModFail"),LineStyle(1),LineColor(kRed));
  if(bkg2ModFail)
    modelFail->plotOn(mframeFail,Components("bkg2ModFail"),LineStyle(1),LineColor(kMagenta));
  

  TCanvas *cpass = MakeCanvas("cpass","cpass",800,1200);
  cpass->SetWindowPosition(cpass->GetWindowTopX()+cpass->GetBorderSize()+800,0);
  cpass->SetTickx(1);
  cpass->SetTicky(1);
  mframePass->Draw();

  legend->Draw();

  string str_tmp = fOutputDir+"/TemplatesandData_Passsample_afterfit.png"; 
  cpass->SaveAs(str_tmp.c_str());
  








  TCanvas *cfail; 
  
  cfail = MakeCanvas("cfail","cfail",800,1200); 
  cfail->SetWindowPosition(cfail->GetWindowTopX()+cfail->GetBorderSize()+800,0);
  cfail->SetTickx(1);
  cfail->SetTicky(1);  
  mframeFail->Draw();

  legend->Draw();
  str_tmp = fOutputDir+"/TemplatesandData_Failsample_afterfit.png"; 
  cfail->SaveAs(str_tmp.c_str());
  }
  




  // if(fout)
  //   fout->Close();

  h_pass_data = (TH1D*)histPass.Clone("h_pass_data");
  h_fail_data = (TH1D*)histFail.Clone("h_fail_data");


  delete legend;



}

void EffRTTFitter::printresults(bool PrintToFile, bool verbose)
{



  if(verbose)
    cout<<"----printresults---"<<endl;

  //  cout<<"fitResult "<<fitResult<<endl;


  ofstream outfile;
  ofstream txtfile;

  if(PrintToFile){
    char txtfname[100];
    sprintf(txtfname,"%s/fitres.txt",fOutputDir.c_str());
    txtfile.open(txtfname);
    assert(txtfile.is_open());
  }

  ostream & outFile = (PrintToFile ? txtfile : cout);

  if(verbose)
    cout<<"----calc chisq pass---"<<endl;



  if(verbose)
    cout<<"----Printing SF and Obs Eff---"<<endl;


  outFile<<"++++++++++++++++++++"<<std::endl;
  outFile<<"++++++++++SF++++++++++"<<std::endl;
  outFile<<"++++++++++++++++++++"<<std::endl;
  if((SIGFR12MASK & 0x001) == 0x001)
    outFile<<"SF (Eff) "<<eff_sig_obs/eff_sig_exp<<" - "<<effl_sig_obs/eff_sig_exp<<" + "<<effh_sig_obs/eff_sig_exp<<std::endl;
  if((SIGFR12MASK & 0x010) == 0x010)
    outFile<<"SF (FR1) "<<fr1_sig_obs/eff_fr1_exp<<" - "<<fabs(fr1l_sig_obs)/eff_fr1_exp<<" + "<<fr1h_sig_obs/eff_fr1_exp<<std::endl;
  if((SIGFR12MASK & 0x100) == 0x100)
    outFile<<"SF (FR2) "<<fr2_sig_obs/eff_fr2_exp<<" - "<<fabs(fr2l_sig_obs)/eff_fr2_exp<<" + "<<fr2h_sig_obs/eff_fr2_exp<<std::endl;
  

  outFile<<"++++++++++Eff (Fit)++++++++++"<<std::endl;
  if((SIGFR12MASK & 0x001) == 0x001)
    outFile<<"Eff "<<eff_sig_obs<<" - "<<effl_sig_obs<<" + "<<effh_sig_obs<<std::endl;
  if((SIGFR12MASK & 0x010) == 0x010)
    outFile<<"FR1 "<<fr1_sig_obs<<" - "<<fabs(fr1l_sig_obs)<<" + "<<fr1h_sig_obs<<std::endl;

    
  if((SIGFR12MASK & 0x100) == 0x100)
    outFile<<"FR2 "<<fr2_sig_obs<<" - "<<fabs(fr2l_sig_obs)<<" + "<<fr2h_sig_obs<<std::endl; 


  if(verbose)
    cout<<"----Printing Exp Eff---"<<endl;


  outFile<<"++++++++++Eff (Expect)++++++++++"<<std::endl; //compute errors
  if((nsigpass_exp+nsigfail_exp) > 0)
    outFile<<"Eff "<<eff_sig_exp<<std::endl;
  if((nbkg1pass_exp+nbkg1fail_exp) > 0)
    outFile<<"FR1 "<<eff_fr1_exp<<std::endl;
  if((nbkg2pass_exp+nbkg2fail_exp) > 0)
    outFile<<"FR2 "<<eff_fr2_exp<<std::endl;

  //    outFile<<"#chi^{2}/dof (Pass) "<<mframePass->chiSquare()<<std::endl; //still need be to interpreted
  //outFile<<"#chi^{2}/dof (Fail) "<<mframeFail->chiSquare()<<std::endl;
    
  outFile<<"#chi^{2}/dof (PassMarco) "<<chisqpass<<std::endl;
  outFile<<"#chi^{2}/dof (FailMarco) "<<chisqfail<<std::endl;


  if(verbose)
    cout<<"----Printing other Obs Values---"<<endl;

    
    
  outFile<<"++++Other Fitted Values++++"<<endl;
  if((SIGFR12MASK & 0x001) == 0x001)
    outFile<<"Nsig "<<nsig_obs<<" (expect: "<<Nsigexp<<")"<<endl;
  if((SIGFR12MASK & 0x010) == 0x010)
    outFile<<"Nbkg1 "<<nbkg1_obs<<" (expect: "<<Nbkg1exp<<")"<<endl;
  if((SIGFR12MASK & 0x100) == 0x100)
    outFile<<"Nbkg2 "<<nbkg2_obs<<" (expect: "<<Nbkg2exp<<")"<<endl;
  outFile<<"*********Fitted Pass**********"<<std::endl;
  outFile<<"Yields (data) "<<(int)fPassTree->GetEntries()<<std::endl;
  if(NsigPass)
    outFile<<"Nsig "<<nsigpass_obs<<" +/- "<<nsigpass_err_obs<<std::endl;
  if(Nbkg1Pass)
    outFile<<"Nbkg1 "<<nbkg1pass_obs<<" +/- "<<nbkg1pass_err_obs<<std::endl;
  if(Nbkg2Pass)
    outFile<<"Nbkg2 "<<nbkg2pass_obs<<" +/- "<<nbkg2pass_err_obs<<std::endl;
  //propagate error fraction  
  if(NsigPass && Nbkg1Pass && Nbkg2Pass){
    outFile<<"FracNsig, FracNbkg1, FracNbkg2 "<<nsigpass_obs/npass_data<<" "<<nbkg1pass_obs/npass_data<<" "<<nbkg1pass_obs/npass_data<<std::endl;
      
    outFile<<"*********EXPECTED FRACMC Pass**********"<<std::endl;
    outFile<<"Signal, Bkg1, Bkg2 "<<nsigpass_exp/npass_data<<" "<<nbkg1pass_exp/npass_data<<" "<<nbkg2pass_exp/npass_data<<std::endl;
  }
  outFile<<"*********Fitted Fail**********"<<std::endl;
  outFile<<"Yields (data) "<<(int)fFailTree->GetEntries()<<std::endl;
  if(NsigFail)
    outFile<<"Nsig "<<nsigfail_obs<<" +/- "<<nsigfail_err_obs<<std::endl;
  if(Nbkg1Fail)
    outFile<<"Nbkg1 "<<nbkg1fail_obs<<" +/- "<<nbkg1fail_err_obs<<std::endl;
  if(Nbkg2Fail)
    outFile<<"Nbkg2 "<<nbkg2fail_obs<<" +/- "<<nbkg2fail_err_obs<<std::endl;
  //propagate error fraction  
  if(NsigFail && Nbkg1Fail && Nbkg2Fail){
    outFile<<"FracNsig, FracNbkg1, FracNbkg2 "<<nsigfail_obs/nfail_data<<" "<<nbkg1fail_obs/nfail_data<<" "<<nbkg1fail_obs/nfail_data<<std::endl;
      
    outFile<<"*********EXPECTED FRACMC Fail**********"<<std::endl;
    outFile<<"Signal, Bkg1, Bkg2 "<<nsigfail_exp/nfail_data<<" "<<nbkg1fail_exp/nfail_data<<" "<<nbkg2fail_exp/nfail_data<<std::endl;
  }
  
    
  outFile<<"&&&&&&&&&&&PREFIT YIELDS&&&&&&&&&&&&&&"<<std::endl;
    
  outFile<<"**********Fail**********"<<std::endl;
  outFile<<"h_signal, h_bkg1, h_bkg2 "<<nsigfail_exp<<" "<<nbkg1fail_exp<<" "<<nbkg2fail_exp<<std::endl;
    
  outFile<<"**********Pass********"<<std::endl;
  outFile<<"h_signal, h_bkg1, h_bkg2 "<<nsigpass_exp<<" "<<nbkg1pass_exp<<" "<<nbkg2pass_exp<<std::endl;
  outFile<<"*********TOT**********"<<std::endl;
  outFile<<"h_signal, h_bkg1, h_bkg2 "<<nsigfail_exp+nsigpass_exp<<" "<<nbkg1fail_exp+nbkg1pass_exp<<" "<<nbkg2fail_exp+nbkg2pass_exp<<std::endl;
  outFile<<"*********TOTMC**********"<<std::endl;
  outFile<<"Pass, Fail TOT "<<nsigpass_exp+nbkg1pass_exp+nbkg2pass_exp<<" "<<nsigfail_exp+nbkg1fail_exp+nbkg2fail_exp<<" "<<nsigpass_exp+nbkg1pass_exp+nbkg2pass_exp+nsigfail_exp+nbkg1fail_exp+nbkg2fail_exp<<std::endl;
  outFile<<"*******Data*******"<<std::endl;
  outFile<<"Pass, Fail, Tot "<<npass_data<<" "<<nfail_data<<" "<<npass_data+nfail_data<<std::endl;
  outFile<<"**************"<<std::endl;
    
  outFile<<"&&&&&&&&&&&POSTFIT YIELDS&&&&&&&&&&&&&&"<<std::endl;
  outFile<<"Yields (data) "<<(int)fPassTree->GetEntries()+(int)fFailTree->GetEntries()<<std::endl;
  if(NsigPass && NsigFail)
    outFile<<"Nsig "<<nsigpass_obs+nsigfail_obs<<" +/- "<<TMath::Sqrt(nsigpass_err_obs*nsigpass_err_obs+nsigfail_err_obs*nsigfail_err_obs)<<std::endl;
  if(Nbkg1Pass && Nbkg1Fail)
    outFile<<"Nbkg1 "<<nbkg1pass_obs+nbkg1fail_obs<<" +/- "<<TMath::Sqrt(nbkg1pass_err_obs*nbkg1pass_err_obs+nbkg1fail_err_obs*nbkg1fail_err_obs)<<std::endl;
  if(Nbkg2Pass && Nbkg2Fail)
    outFile<<"Nbkg2 "<<nbkg2pass_obs+nbkg2fail_obs<<" +/- "<<TMath::Sqrt(nbkg2pass_err_obs*nbkg2pass_err_obs+nbkg2fail_err_obs*nbkg2fail_err_obs)<<std::endl;
  if(NsigPass && NsigFail && Nbkg1Pass && Nbkg1Fail && Nbkg2Pass && Nbkg2Fail)
    outFile<<"TOTMC "<<nsigpass_obs+nsigfail_obs+nbkg1pass_obs+nbkg1fail_obs+nbkg2pass_obs+nbkg2fail_obs<<endl;
   
  outFile << endl;

  if(PrintToFile) txtfile.close();

  if(verbose)
    std::cout<<"Done printresults"<<endl;

}

void EffRTTFitter::printRooFitresults(bool verbose)
{
  if(verbose)
    cout<<"----printRooFitresults---"<<endl;

  //PRINT FITTED NUISANCE< CHISQ, ETC...
  ofstream txtfile2;
  char txtfname2[100];
  sprintf(txtfname2,"%s/fitres2.txt",fOutputDir.c_str());
  txtfile2.open(txtfname2);
  assert(txtfile2.is_open());

  txtfile2 <<"DONE WITH PRINTING MY RESULTS"<< endl;
  txtfile2 << endl;  txtfile2 << endl;
  txtfile2 <<"NOW DUMPING fit results from fitResult->printStream"<< endl;
  fitResult->printStream(txtfile2,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile2 <<"PRINT VALUE"<< endl;
  fitResult->printValue(txtfile2);
  txtfile2.close();

  if(verbose)
    cout<<"----Done printRooFitresults2---"<<endl;

}

void EffRTTFitter::CloneTemplates(bool verbose)
{
  if(verbose)
    cout<<"---CloneTemplates---"<<endl;

  if(h_pass_signal){
    h_pass_signal_af = (TH1D*)h_pass_signal->Clone("h_pass_signal_af");
    h_pass_signal_af->Sumw2();
  }
  if(h_pass_bkg1){
    h_pass_bkg1_af = (TH1D*)h_pass_bkg1->Clone("h_pass_bkg1_af");
    h_pass_bkg1_af->Sumw2();
  }
  if(h_pass_bkg2){
    h_pass_bkg2_af = (TH1D*)h_pass_bkg2->Clone("h_pass_bkg2_af");
    h_pass_bkg2_af->Sumw2();
  }
  if(h_fail_signal){
    h_fail_signal_af = (TH1D*)h_fail_signal->Clone("h_fail_signal_af");
    h_fail_signal_af->Sumw2();
  }
  if(h_fail_bkg1){
    h_fail_bkg1_af = (TH1D*)h_fail_bkg1->Clone("h_fail_bkg1_af");
    h_fail_bkg1_af->Sumw2();
  }
  if(h_fail_bkg2){
    h_fail_bkg2_af = (TH1D*)h_fail_bkg2->Clone("h_fail_bkg2_af");
    h_fail_bkg2_af->Sumw2();
  }

  if(verbose)
    cout<<"---Done CloneTemplates---"<<endl;

}

void EffRTTFitter::SetTemplatesColors(TH1D *hMC_1, TH1D *hMC_2, TH1D *hMC_3, bool verbose){

  if(verbose)
    cout<<"---SetTemplatesColors---"<<endl;

  if(hMC_1){
    hMC_1->SetLineColor(kRed+2);    hMC_1->SetFillColor(kRed);     hMC_1->SetFillStyle(kRed);
  }
  if(hMC_2){
    hMC_2->SetLineColor(kMagenta+2);     hMC_2->SetFillColor(kMagenta);     hMC_2->SetFillStyle(kMagenta);
  }
  if(hMC_3){
    hMC_3->SetLineColor(kGreen+2);     hMC_3->SetFillColor(kGreen);     hMC_3->SetFillStyle(kGreen);
  }


  if(verbose)
    cout<<"---Done SetTemplatesColors---"<<endl;

}

// void EffRTTFitter::PlotDataVsMC(TH1D *hMC_1, TH1D *hMC_2, TH1D *hMC_3, TH1D hData, std::string pictname, bool verbose){
//   //hMC_1 = h_pass_bkg1
//   //hMC_2 = h_pass_bkg2
//   //hMC_3= h_pass_signal

//   if(verbose)
//     cout<<"----PlotDataVsMC---"<<endl;

//   THStack *h_tot = new THStack("h_tot",""); 
//   if(hMC_1)
//     h_tot->Add(hMC_1);
//   if(hMC_2)
//     h_tot->Add(hMC_2);
//   if(hMC_3)
//     h_tot->Add(hMC_3);

//   TH1D *h_tot_ = (TH1D *)hMC_2->Clone("h_tot_");
//   h_tot_->Sumw2();
//   if(hMC_1)
//     h_tot_->Add(hMC_1);
//   if(hMC_3)
//     h_tot_->Add(hMC_3);

//   TH1D *h_ratio = (TH1D *)hData.Clone("h_tot_");
//   h_ratio->Sumw2();
//   h_ratio->Divide(h_tot_);  


//   TLine * a = new TLine();
//   a->SetLineColor(2);
//   a->SetLineWidth(2);
//   a->SetLineStyle(1);

//   TLegend *legend = new TLegend(0.65,0.54,0.99,0.99);
//   legend->SetFillColor(0);
//   legend->AddEntry(&hData,"Data" , "lp");
//   if(hMC_3)
//     legend->AddEntry(hMC_3, "Signal" , "lf");
//   if(hMC_2)
//     legend->AddEntry(hMC_2, "Non-tt1l BG" , "lf");
//   if(hMC_1)
//     legend->AddEntry(hMC_1, "tt1l Combinatorial BG" , "lf");


//   std::cout<<"plotting data and model prefit via TH1"<<std::endl;
    
//   TCanvas* c = new TCanvas("c","c",800,1200);
//   //    c->SetWindowPosition(c->GetWindowTopX()+c->GetBorderSize()+800,0);
//   //c->SetTickx(1);
//   //c->SetTicky(1);
//   TPad *pad1 = NULL; TPad *pad2 = NULL;
//   pad1 = new TPad("pad1","",0.01,0.41,0.99,1,0);
//   pad2 = new TPad("pad2","Data-Exp",0.01,0.,0.99,0.4);
//   pad2->SetLeftMargin(0.059);
//   pad2->SetRightMargin(0.02);
//   pad1->SetRightMargin(0.02);
//   pad1->SetLeftMargin(0.059);
//   pad1->Draw(); pad2->Draw();

//   pad1->cd();
    
//   cout<<"h_tot "<<h_tot<<endl;

//   if(h_tot->GetMaximum()<hData.GetMaximum())
//     h_tot->SetMaximum(hData.GetMaximum()*1.1);
//   h_tot->Draw("hist");
  

//   hData.Draw("Esame");

//     //    h_tot->GetXaxis()->SetTitle("Top Mass [GeV/c]"); //crash...figure out why

//     //    c->Update();


//   legend->Draw();

//   pad2->cd();

//   h_ratio->GetYaxis()->SetRangeUser(0.5,1.5);
//   h_ratio->Draw("E");
//   a->DrawLine(fFitLo,1.,fFitHi,1.);

//   string str_tmp = fOutputDir+"/"+pictname;
//   c->SaveAs(str_tmp.c_str());

//   if(pad1)  delete pad1;
//   if(pad2)  delete pad2;
//   if(c)     delete c;
//   if(legend)delete legend;
//   if(h_tot) delete h_tot;
//   if(a)     delete a;

//   if(verbose)
//     cout<<"Done"<<endl;
// }

void EffRTTFitter::PlotDataVsMC(TH1D *hMC_1, TH1D *hMC_2, TH1D *hMC_3, TH1D hData, bool doLogY, std::string picname, bool verbose){
  //hMC_1 = h_pass_bkg1
  //hMC_2 = h_pass_bkg2
  //hMC_3= h_pass_signal

  if(verbose)
    cout<<"----PlotDataVsMC---"<<endl;

  //make Data/MC ratio

  TH1D *h_tot_ = (TH1D *)hMC_2->Clone("h_tot_");
  h_tot_->Sumw2();
  if(hMC_1)
    h_tot_->Add(hMC_1);
  if(hMC_3)
    h_tot_->Add(hMC_3);

  TH1D *hRatio = (TH1D *)hData.Clone("h_tot_");
  hRatio->Sumw2();
  hRatio->Divide(h_tot_);  

  //plot
  TString xaxislabel = "Top Mass [GeV/C]";
  const float xmin=0;  const float xmax=1000; 

  TCanvas *c = MakeCanvas("c","c",800,600);
  c->Divide(1,2,0,0);
  c->cd(1)->SetPad(0,0.3,1.0,1.0);
  c->cd(1)->SetTopMargin(0.1);
  c->cd(1)->SetBottomMargin(0.01);
  c->cd(1)->SetLeftMargin(0.15);
  c->cd(1)->SetRightMargin(0.07);
  c->cd(1)->SetTickx(1);
  c->cd(1)->SetTicky(1);
  c->cd(2)->SetPad(0,0,1.0,0.3);
  c->cd(2)->SetTopMargin(0.05);
  c->cd(2)->SetBottomMargin(0.45);
  c->cd(2)->SetLeftMargin(0.15);
  c->cd(2)->SetRightMargin(0.07);
  c->cd(2)->SetTickx(1);
  c->cd(2)->SetTicky(1);


  CPlot plot(picname,"",xaxislabel,"Events");
  plot.AddHist1D(&hData,"Data","E");
  plot.AddToStack(hMC_2,"Other Background",kRed-3,kRed-2);
  plot.AddToStack(hMC_1,"t#bar{t}(1l) Combinatorial",kOrange-3,kOrange-2);
  plot.AddToStack(hMC_3,"t#bar{t}(1l) Hadronic Matched",kAzure-3,kAzure-2);
  plot.AddTextBox("35.87 fb^{-1} (13 TeV)",0.66,0.99,0.95,0.925,0,kBlack);
  plot.AddTextBox("CMS",0.20,0.88,0.30,0.82,0,kBlack,62);
  plot.AddTextBox("Preliminary",0.20,0.82,0.37,0.77,0,kBlack,52);
  plot.SetYRange(0.001,1.4*hData.GetBinContent(hData.GetMaximumBin()));
  if(doLogY)
    plot.SetLogy();

  hRatio->SetMarkerSize(0.8);  
  CPlot plotRatio(picname,"",xaxislabel,"Data / MC");
  plotRatio.AddHist1D(hRatio,"EX0",kBlack);
  plotRatio.SetYRange(0.6,1.5);
  plotRatio.AddLine(xmin,1.0,xmax,1.0,kBlack,3);



  plot.TransLegend(-0.03,-0.05);
  plot.Draw(c,false,"png",1);
  plotRatio.Draw(c,true,"png",2); 
  plot.Draw(c,false,"pdf",1);
  plotRatio.Draw(c,true,"pdf",2); 

  //  gSystem->Rename(picname+".png",fOutputDir+"/"+picname+".png");
  string orig   = picname+".png";
  string destin = fOutputDir+"/"+orig;
  gSystem->Rename(orig.c_str(), destin.c_str());
  orig   = picname+".pdf";
  destin = fOutputDir+"/"+orig;
  gSystem->Rename(orig.c_str(), destin.c_str());


  if(verbose)
    cout<<"Done"<<endl;
}


void EffRTTFitter::RooPlotDataVsMC(RooRealVar m, TLegend *legend, string pictnamepass, string pictnamefail, bool verbose){
  //TO BE IMPROVED: it would be nice to plot either pass or fail to cut the # lines by two
  if(verbose)
    cout<<"---RooPlotDataVsMC---"<<endl;



  //pass
  RooPlot* mframePass_bf = m.frame(Bins(int(fFitHi-fFitLo)/BIN_SIZE_PASS));
  mframePass_bf->SetTitle("");
  dataPass->plotOn(mframePass_bf,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));    
  modelPass->plotOn(mframePass_bf);
  if(MODEL==1 && sigModPass)
    modelPass->plotOn(mframePass_bf,Components("sigModPass"),LineStyle(1),LineColor(kGreen));
  if(MODEL==2 && sigModPass2)
    modelPass->plotOn(mframePass_bf,Components("signalPass"),LineStyle(1),LineColor(kGreen));
  if(bkg1ModPass)
    modelPass->plotOn(mframePass_bf,Components("bkg1ModPass"),LineStyle(1),LineColor(kRed));
  if(bkg2ModPass)
    modelPass->plotOn(mframePass_bf,Components("bkg2ModPass"),LineStyle(1),LineColor(kMagenta));
  TCanvas* cpass_bf = new TCanvas("cpass_bf","cpass_bf",800,1200);
  cpass_bf->SetWindowPosition(cpass_bf->GetWindowTopX()+cpass_bf->GetBorderSize()+800,0);
  cpass_bf->SetTickx(1);
  cpass_bf->SetTicky(1);
  mframePass_bf->Draw();
  legend->Draw(); 

  string str_tmp = fOutputDir+"/"+pictnamepass;
  cpass_bf->SaveAs(str_tmp.c_str());

  //fail
    RooPlot* mframeFail_bf; 
    mframeFail_bf = m.frame(Bins(int(fFitHi-fFitLo)/BIN_SIZE_PASS));
    mframeFail_bf->SetTitle("");
    dataFail->plotOn(mframeFail_bf,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));    
    modelFail->plotOn(mframeFail_bf);
    if(MODEL==1 && sigModFail)
      modelFail->plotOn(mframeFail_bf,Components("sigModFail"),LineStyle(1),LineColor(kGreen));
    if(MODEL==2 && sigModFail2)
      modelFail->plotOn(mframeFail_bf,Components("signalFail"),LineStyle(1),LineColor(kGreen));
    
    modelFail->plotOn(mframeFail_bf,Components("bkg1ModFail"),LineStyle(1),LineColor(kRed));
    modelFail->plotOn(mframeFail_bf,Components("bkg2ModFail"),LineStyle(1),LineColor(kMagenta));
    TCanvas* cfail_bf = new TCanvas("cfail_bf","cfail_bf",800,1200);
    cfail_bf->SetWindowPosition(cfail_bf->GetWindowTopX()+cfail_bf->GetBorderSize()+800,0);
    cfail_bf->SetTickx(1);
    cfail_bf->SetTicky(1);
    mframeFail_bf->Draw();
    
    
    legend->Draw();
    
    str_tmp = fOutputDir+"/"+pictnamefail; 
    cfail_bf->SaveAs(str_tmp.c_str());




  if(verbose)
    cout<<"---Done RooPlotDataVsMC---"<<endl;
}

double EffRTTFitter::CalcChisq(TH1D *hMC_1, TH1D *hMC_2, TH1D *hMC_3, TH1D hData, bool verbose){
  if(verbose)
    cout<<"---CalcChisq---"<<endl;

  double chisq_tmp;

  TH1D *h_tot = (TH1D *)hMC_2->Clone("h_tot");
  if(hMC_1)
    h_tot->Add(hMC_1);
  if(hMC_3)
    h_tot->Add(hMC_3);

  TH1D *h_diff = (TH1D *)hData.Clone("h_diff");
  h_diff->Add(h_tot,-1);

  for(int ib=0; ib<h_diff->GetNbinsX(); ib++){
    double bc = h_diff->GetBinContent(ib);
    double berr = h_diff->GetBinError(ib);
    if(berr != 0)
      chisq_tmp += bc*bc/(berr*berr);
  } 


  if(verbose)
    cout<<"---Done CalcChisq---"<<endl;


  return chisq_tmp;
}

void EffRTTFitter::CreateOutFile(bool verbose){
  if(verbose)
    cout<<"---CreateOutFile---"<<endl;



  string str_fname = fOutputDir+"/out.root";

  if(verbose)
    cout<<"str_fname "<<str_fname<<endl;

  fout = new TFile(str_fname.c_str(),"RECREATE");



  if(verbose)
    cout<<"---CreateOutFile---"<<endl;
}

void EffRTTFitter::WriteHistos(bool verbose){
  if(verbose)
    cout<<"---WriteHistos---"<<endl;
  
  fout->cd();
  
  
  h_pass_data->Write();
  h_fail_data->Write();
  if(h_pass_bkg1)
    h_pass_bkg1->Write();
  if(h_pass_bkg2)
    h_pass_bkg2->Write();
  if(h_pass_signal)
    h_pass_signal->Write();
  
  if(h_fail_bkg1)
    h_fail_bkg1->Write();
  if(h_fail_bkg2)
    h_fail_bkg2->Write();
  if(h_fail_signal)
    h_fail_signal->Write();



  if(h_pass_bkg1_af)
    h_pass_bkg1_af->Write();
  if(h_pass_bkg2_af)
    h_pass_bkg2_af->Write();
  if(h_pass_signal_af)
    h_pass_signal_af->Write();

  if(h_fail_bkg1_af)
    h_fail_bkg1_af->Write();
  if(h_fail_bkg2_af)
    h_fail_bkg2_af->Write();
  if(h_fail_signal_af)
    h_fail_signal_af->Write();



  if(verbose)
    cout<<"---WriteHistos---"<<endl;
}

