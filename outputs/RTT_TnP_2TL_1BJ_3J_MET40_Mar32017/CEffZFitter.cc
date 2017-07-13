#include "CEffZFitter.hh"
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
#include "ZSignals.hh"
#include "CEffUser1D.hh"


// RooFit headers
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"


#include "RooAbsPdf.h"
#include "RooHistPdf.h"




// bin size constants
// #define BIN_SIZE_PASS 20
// #define BIN_SIZE_FAIL 20 //was 20
//MODEL=1 nominal, MODEL=2 IF you wanna smear with gaussian (for systematic)


//--------------------------------------------------------------------------------------------------
CEffZFitter::CEffZFitter():
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
VERBOSE(true)
{}

//--------------------------------------------------------------------------------------------------
CEffZFitter::~CEffZFitter()
{
  delete fPassTree;     fPassTree=0;  
  
}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::initialize(const std::string infname, const std::string outdir, const std::string temfname1, const std::string temfname2, const std::string temfname3, const double TOPMVACUT, const bool MAKETEMPLATES, const float FITLO, const float FITHI, const float BINSIZEPASS, const float BINSIZEFAIL, const unsigned int SIGNALFAKERATE12MASK, const float FR2SCALEFACTOR, const float FR2ERROR)
{
  std::cout << "   [CEffZFitter] Initializing... " << std::endl;


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
void CEffZFitter::computeEff()
{
  assert(fIsInitialized);
  
  std::cout << "   [CEffZFitter] Computing efficiencies..." << std::endl;


  //------------------------------------------------------------------------------------------------
  // Efficiency calculation
  //================================================================================================

  
  double eff, errl, errh;  
  performFit(eff, errl, errh,
	     fPassTree, fFailTree); 


  std::cout<<"eff, errl, errh"<<errh<<" "<<errl<<" "<<errh<<std::endl;  
  
}


//--------------------------------------------------------------------------------------------------
void CEffZFitter::makeBinnedTemplates(const std::string temfname, string name)
{
  std::cout << "   [CEffZFitter] Creating binned templates... "; std::cout.flush();
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
void CEffZFitter::performFit(double &resEff, double &resErrl, double &resErrh,
                             TTree *passTree, TTree *failTree)

{
  string str_fname = fOutputDir+"/out.root";
  TFile *fout = new TFile(str_fname.c_str(),"RECREATE");

  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  std::cout << " ...performFit..." << std::endl; 
  
  string str_title;
  str_title = "Top Mass [GeV/c^2]";
  RooRealVar m("m",str_title.c_str(),fFitLo,fFitHi);

  
  m.setBins(10000);
  
  
  TFile *histfile_signal = 0;
  TFile *histfile_bkg1 = 0;  TFile *histfile_bkg2 = 0;
  {
    string str_tmp = fOutputDir+"/binnedTemplates_signal.root";
    histfile_signal = new TFile(str_tmp.c_str());
    assert(histfile_signal);

   
    str_tmp = fOutputDir+"/binnedTemplates_bkg1.root";
    histfile_bkg1 = new TFile(str_tmp.c_str());  
    str_tmp = fOutputDir+"/binnedTemplates_bkg2.root";
    histfile_bkg2 = new TFile(str_tmp.c_str());
    assert(histfile_bkg1);  assert(histfile_bkg2);
  }
  
  // Define categories
  RooCategory sample("sample","");
  sample.defineType("Pass",1);
  sample.defineType("Fail",2);
  
  RooAbsData *dataPass=0;
  RooAbsData *dataFail=0;
  TH1D histPass("histPass","",int(fFitHi-fFitLo)/BIN_SIZE_PASS,fFitLo,fFitHi); 
  TH1D histFail("histFail","",int(fFitHi-fFitLo)/BIN_SIZE_FAIL,fFitLo,fFitHi);
  RooAbsData *dataCombined=0;
  
  std::cout<<"passtree entries "<<passTree->Draw("m>>histPass","w")<<std::endl;
  std::cout<<"failtree entries "<<failTree->Draw("m>>histFail","w")<<std::endl;
  dataPass = new RooDataHist("dataPass","dataPass",RooArgSet(m),&histPass);
  dataFail = new RooDataHist("dataFail","dataFail",RooArgSet(m),&histFail);
  //m.setBins(100);  

  dataCombined = new RooDataHist("dataCombined","dataCombined",RooArgList(m),
				 RooFit::Index(sample),
				 RooFit::Import("Pass",*((RooDataHist*)dataPass)),
				 RooFit::Import("Fail",*((RooDataHist*)dataFail)));
  
  

  std::cout<<"qui"<<std::endl;
  RooHistPdf *sigModPass = 0;
  CSignalModel     *sigModPass2 = 0; //smeared with a gaussian
  RooHistPdf *bkg1ModPass = 0;
  CSignalModel     *bkg1ModPass2 = 0; //smeared with a gaussian
  RooHistPdf *bkg2ModPass = 0;
  CSignalModel     *bkg2ModPass2 = 0; //smeared with a gaussian
  RooHistPdf *sigModFail = 0;
  CSignalModel     *sigModFail2 = 0; //smeared with a gaussian
  RooHistPdf *bkg1ModFail = 0;
  CSignalModel     *bkg1ModFail2 = 0; //smeared with a gaussian
  RooHistPdf *bkg2ModFail = 0;
  CSignalModel     *bkg2ModFail2 = 0; //smeared with a gaussian
  RooDataHist *sigModPass_rdh = 0; 
  RooDataHist *bkg1ModPass_rdh = 0;
  RooDataHist *bkg2ModPass_rdh = 0;  
  RooDataHist *sigModFail_rdh = 0;
  RooDataHist *bkg1ModFail_rdh = 0;
  RooDataHist *bkg2ModFail_rdh = 0;   

  TH1D *h_pass_signal = 0; 
  TH1D *h_pass_bkg1  = 0;  
  TH1D *h_pass_bkg2  = 0;  
  TH1D *h_fail_signal  = 0; 
  TH1D *h_fail_bkg1  = 0;  
  TH1D *h_fail_bkg2  = 0;  
  TH1D *h_pass_signal_n  = 0;  
  TH1D *h_pass_bkg1_n  = 0;  
  TH1D *h_pass_bkg2_n  = 0; 
  TH1D *h_fail_signal_n  = 0;  
  TH1D *h_fail_bkg1_n  = 0;  
  TH1D *h_fail_bkg2_n  = 0; 
  TH1D *h_pass_signal_af  = 0;  
  TH1D *h_pass_bkg1_af  = 0;  
  TH1D *h_pass_bkg2_af  = 0;  //after fit (normalized to fitted fractios) 
  TH1D *h_fail_signal_af = 0; 
  TH1D *h_fail_bkg1_af  = 0;  
  TH1D *h_fail_bkg2_af  = 0;  
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

  
    // Make PDF from MC histograms
  
  
  std::cout<<"qui3"<<std::endl;
  

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
                                                                                                              
  // Make PDF from MC histograms      
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


  float nfail_sig = 0; float nfail_bkg1 = 0; float nfail_bkg2 = 0; float npass_data = 0;
  float npass_sig = 0; float npass_bkg1 = 0; float npass_bkg2 = 0; float nfail_data = 0;
  
  if(h_fail_signal)
    nfail_sig = h_fail_signal->Integral();
  if(h_fail_bkg1)
    nfail_bkg1 = h_fail_bkg1->Integral();
  if(h_fail_bkg2)
    nfail_bkg2 = h_fail_bkg2->Integral();

  if(h_pass_signal)
    npass_sig = h_pass_signal->Integral();
  if(h_pass_bkg1)
    npass_bkg1 = h_pass_bkg1->Integral();
  if(h_pass_bkg2)
    npass_bkg2 = h_pass_bkg2->Integral();

  npass_data = histPass.Integral();
  nfail_data = histFail.Integral();

  if(h_pass_signal){
    h_pass_signal_n = (TH1D*)h_pass_signal->Clone("h_pass_signal_n");
    h_pass_signal_n->Sumw2();
    h_pass_signal_n->Scale(npass_data/h_pass_signal_n->Integral());
  }
  if(h_pass_bkg1){
    h_pass_bkg1_n = (TH1D*)h_pass_bkg1->Clone("h_pass_bkg1_n");
    h_pass_bkg1_n->Sumw2();
    h_pass_bkg1_n->Scale(npass_data/h_pass_bkg1_n->Integral());
  }
  if(h_pass_bkg2){
    h_pass_bkg2_n = (TH1D*)h_pass_bkg2->Clone("h_pass_bkg2_n");
    h_pass_bkg2_n->Sumw2();
    h_pass_bkg2_n->Scale(npass_data/h_pass_bkg2_n->Integral());
  }
  if(h_fail_signal){
    h_fail_signal_n = (TH1D*)h_fail_signal->Clone("h_fail_signal_n");
    h_fail_signal_n->Sumw2();
    h_fail_signal_n->Scale(nfail_data/h_fail_signal_n->Integral());
  }
  if(h_fail_bkg1){
    h_fail_bkg1_n = (TH1D*)h_fail_bkg1->Clone("h_fail_bkg1_n");
    h_fail_bkg1_n->Sumw2();
    h_fail_bkg1_n->Scale(nfail_data/h_fail_bkg1_n->Integral());
  }
  if(h_fail_bkg2){
    h_fail_bkg2_n = (TH1D*)h_fail_bkg2->Clone("h_fail_bkg2_n");
    h_fail_bkg2_n->Sumw2();
    h_fail_bkg2_n->Scale(nfail_data/h_fail_bkg2_n->Integral());
  }
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
  
  
  std::cout<<"PREFIT YIELDS"<<std::endl;
  std::cout<<"**********Fail**********"<<std::endl;
  std::cout<<"h_signal, h_bkg1, h_bkg2 "<<nfail_sig<<" "<<nfail_bkg1<<" "<<nfail_bkg2<<std::endl;

  std::cout<<"**********Pass********"<<std::endl;
  std::cout<<"h_signal, h_bkg1, h_bkg2 ratiobkg1/bkg2 "<<npass_sig<<" "<<npass_bkg1<<" "<<npass_bkg2<<" "<<npass_bkg2/npass_bkg1<<std::endl;
  
  std::cout<<"*********TOT**********"<<std::endl;
  std::cout<<"h_signal, h_bkg1, h_bkg2 "<<nfail_sig+npass_sig<<" "<<nfail_bkg1+npass_bkg1<<" "<<nfail_bkg2+npass_bkg2<<std::endl;
  std::cout<<"*********TOTMC**********"<<std::endl;
  std::cout<<"Pass, Fail TOT "<<npass_sig+npass_bkg1+npass_bkg2<<" "<<nfail_sig+nfail_bkg1+nfail_bkg2<<" "<<npass_sig+npass_bkg1+npass_bkg2+nfail_sig+nfail_bkg1+nfail_bkg2<<std::endl;
  std::cout<<"*******Data*******"<<std::endl;
  std::cout<<"Pass, Fail, Tot "<<npass_data<<" "<<nfail_data<<" "<<npass_data+nfail_data<<std::endl;
  std::cout<<"**************"<<std::endl;


  double Nsigexp = npass_sig+nfail_sig;
  double Nbkg1exp = npass_bkg1+nfail_bkg1;
  double Nbkg2exp = npass_bkg2+nfail_bkg2;

  double Effexp = npass_sig/Nsigexp;  
  double Fr1exp = npass_bkg1/Nbkg1exp;
  double Fr2exp = npass_bkg2/Nbkg2exp;

  
  // Define free parameters


  RooRealVar Nsig("Nsig","Signal Yield",Nsigexp*0.9,0,Nsigexp*1.1);
  RooRealVar Nbkg1("Nbkg1","Bkg1 Yield",Nbkg1exp*0.9,0,Nbkg1exp*1.1);
  RooRealVar Nbkg2("Nbkg2","Bkg2 Yield",Nbkg2exp*0.9,0,Nbkg2exp*1.1);


  RooRealVar eff("eff","Efficiency",Effexp,Effexp/2,1);
  RooRealVar FR1("FR1","FR Bkg1",Fr1exp,Fr1exp/2,1);
  double fr2min = Fr2exp*FR2SF - FR2ERR;
  double fr2max = Fr2exp*FR2SF + FR2ERR;
  std::cout<<"&&&&&&&& Fr2exp*FR2SF, fr2min, fr2max, "<<Fr2exp*FR2SF<<" "<<fr2min<<" "<<fr2max<<std::endl;
  RooRealVar FR2("FR2","FR Bkg2",Fr2exp*FR2SF,fr2min,fr2max);


  cout<<"qui2"<<endl;

  

  RooFormulaVar *Nbkg1Pass = 0; 
  RooFormulaVar *Nbkg1Fail = 0; 
  RooFormulaVar *Nbkg2Pass = 0; 
  RooFormulaVar *Nbkg2Fail = 0;
  RooFormulaVar  *NsigPass = 0;
  RooFormulaVar  *NsigFail = 0;

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



  cout<<"qui3"<<endl;
  


 
  RooAddPdf *modelPass=0, *modelFail=0;

  if(NsigPass  && Nbkg1Pass  && Nbkg2Pass){
    if(MODEL==1 && sigModPass && bkg1ModPass && bkg2ModPass)
      modelPass = new RooAddPdf("modelPass","Model for PASS sample",RooArgList(*sigModPass,*bkg1ModPass,*bkg2ModPass),RooArgList(*NsigPass,*Nbkg1Pass,*Nbkg2Pass));
    else if(MODEL==2 && sigModPass2 && bkg1ModPass2 && bkg2ModPass2){
      //      modelPass = new RooAddPdf("modelPass","Model for PASS sample",RooArgList(*(sigModPass2->model),*(bkg1ModPass2->model),*(bkg2ModPass2->model)),RooArgList(*NsigPass,*Nbkg1Pass,*Nbkg2Pass)); //MT Jun 1

      modelPass = new RooAddPdf("modelPass","Model for PASS sample",RooArgList(*(sigModPass2->model),*(bkg1ModPass2->model),*bkg2ModPass),RooArgList(*NsigPass,*Nbkg1Pass,*Nbkg2Pass));
      //modelPass = new RooAddPdf("modelPass","Model for PASS sample",RooArgList(*(sigModPass2->model),*bkg1ModPass,*(bkg2ModPass2->model)),RooArgList(*NsigPass,*Nbkg1Pass,*Nbkg2Pass));
    }
  }
  else if (bkg2ModPass && Nbkg2Pass)
    modelPass = new RooAddPdf("modelPass","Model for PASS sample",RooArgList(*bkg2ModPass),RooArgList(*Nbkg2Pass));

  if(NsigFail && Nbkg1Fail&& Nbkg2Fail){
    if(MODEL==1 && sigModFail && bkg1ModFail && bkg2ModFail )
      modelFail = new RooAddPdf("modelFail","Model for FAIL sample",RooArgList(*sigModFail,*bkg1ModFail,*bkg2ModFail),RooArgList(*NsigFail,*Nbkg1Fail,*Nbkg2Fail)); 
    else if(MODEL==2 && sigModFail2 && bkg1ModFail2 && bkg2ModFail2 ){
      //	modelFail = new RooAddPdf("modelFail","Model for FAIL sample",RooArgList(*(sigModFail2->model),*(bkg1ModFail2->model),*(bkg2ModFail2->model)),RooArgList(*NsigFail,*Nbkg1Fail,*Nbkg2Fail)); //MT Jun 1
      
      modelFail = new RooAddPdf("modelFail","Model for FAIL sample",RooArgList(*(sigModFail2->model),*(bkg1ModFail2->model),*bkg2ModFail),RooArgList(*NsigFail,*Nbkg1Fail,*Nbkg2Fail)); 
      //modelFail = new RooAddPdf("modelFail","Model for FAIL sample",RooArgList(*(sigModFail2->model),*bkg1ModFail,*(bkg2ModFail2->model)),RooArgList(*NsigFail,*Nbkg1Fail,*Nbkg2Fail)); 
    }
  }
  else if (bkg2ModFail && Nbkg2Fail)
    modelFail = new RooAddPdf("modelFail","Model for PASS sample",RooArgList(*bkg2ModFail),RooArgList(*Nbkg2Fail));


  
  cout<<"*****before setting totalPdf"<<endl;
  cout<<"modelPass "<<modelPass<<endl;
  cout<<"modelFail "<<modelFail<<endl;
  
  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(*modelPass,"Pass"); 
  totalPdf.addPdf(*modelFail,"Fail");

  cout<<"qui41"<<endl;

  if(h_pass_signal){
    h_pass_signal->SetLineColor(kGreen+2);     h_pass_signal->SetFillColor(kGreen);     h_pass_signal->SetFillStyle(kGreen);
  }
  if(h_pass_bkg1){
    h_pass_bkg1->SetLineColor(kRed+2);    h_pass_bkg1->SetFillColor(kRed);     h_pass_bkg1->SetFillStyle(kRed);
  }
  if(h_pass_bkg2){
    h_pass_bkg2->SetLineColor(kMagenta+2);     h_pass_bkg2->SetFillColor(kMagenta);     h_pass_bkg2->SetFillStyle(kMagenta);
  }
    cout<<"qui42"<<endl;
  if(h_fail_signal){
    h_fail_signal->SetLineColor(kGreen+2);     h_fail_signal->SetFillColor(kGreen); h_fail_signal->SetFillStyle(kGreen);
  }
  if(h_fail_bkg1){
    h_fail_bkg1->SetLineColor(kRed+2);    h_fail_bkg1->SetFillColor(kRed); h_fail_bkg1->SetFillStyle(kRed);
  }
  if(h_fail_bkg2){
    h_fail_bkg2->SetLineColor(kMagenta+2);     h_fail_bkg2->SetFillColor(kMagenta); h_fail_bkg2->SetFillStyle(kMagenta);
  }

  cout<<"qui43"<<endl;

  if(h_pass_signal_n)
    h_pass_signal_n->SetLineColor(kGreen+2);
  if(h_pass_bkg1_n)
    h_pass_bkg1_n->SetLineColor(kRed+2);    
  if(h_pass_bkg2_n)
    h_pass_bkg2_n->SetLineColor(kMagenta+2);

  cout<<"qui44"<<endl;
  if(h_pass_signal_n)
    h_pass_signal_n->SetLineColor(kGreen+2);
  if(h_pass_bkg1_n)
    h_pass_bkg1_n->SetLineColor(kRed+2);    
  if(h_pass_bkg2_n)
    h_pass_bkg2_n->SetLineColor(kMagenta+2);
  if(h_fail_signal_n)
    h_fail_signal_n->SetLineColor(kGreen+2);
  if(h_fail_bkg1_n)
    h_fail_bkg1_n->SetLineColor(kRed+2);    
  if(h_fail_bkg2_n)
    h_fail_bkg2_n->SetLineColor(kMagenta+2);

  cout<<"qui45"<<endl;
  if(h_pass_signal_af){
    h_pass_signal_af->SetLineColor(kGreen+2);     h_pass_signal_af->SetFillColor(kGreen);//     h_pass_signal_af->SetFillStyle(kGreen);
  }
  if(h_pass_bkg1_af){
    h_pass_bkg1_af->SetLineColor(kRed+2);    h_pass_bkg1_af->SetFillColor(kRed);    // h_pass_bkg1_af->SetFillStyle(kRed);
  }
  if(h_pass_bkg2_af){
    h_pass_bkg2_af->SetLineColor(kMagenta+2);     h_pass_bkg2_af->SetFillColor(kMagenta);  //   h_pass_bkg2_af->SetFillStyle(kMagenta);
  }
  if(h_fail_signal_af){
    h_fail_signal_af->SetLineColor(kGreen+2);     h_fail_signal_af->SetFillColor(kGreen); //h_fail_signal_af->SetFillStyle(kGreen);
  }
  if(h_fail_bkg1_af){
    h_fail_bkg1_af->SetLineColor(kRed+2);    h_fail_bkg1_af->SetFillColor(kRed); //h_fail_bkg1_af->SetFillStyle(kRed);
  }
  if(h_fail_bkg2_af){
    h_fail_bkg2_af->SetLineColor(kMagenta+2);     h_fail_bkg2_af->SetFillColor(kMagenta);// h_fail_bkg2_af->SetFillStyle(kMagenta);
  }


  cout<<"qui5"<<endl;


  // Plot the imported histogram(s)                                                                                  
  THStack *h_pass_tot = new THStack("h_pass_tot",""); 
  if(h_pass_bkg1)
    h_pass_tot->Add(h_pass_bkg1);
  if(h_pass_bkg2)
  h_pass_tot->Add(h_pass_bkg2);
  if(h_pass_signal)
  h_pass_tot->Add(h_pass_signal);


  cout<<"qui6"<<endl;

  TH1D *h_pass_tot_ = (TH1D *)h_pass_bkg2->Clone("h_pass_tot_");
  h_pass_tot_->Sumw2();
  if(h_pass_bkg1)
    h_pass_tot_->Add(h_pass_bkg1);
  if(h_pass_signal)
    h_pass_tot_->Add(h_pass_signal);

  TH1D *h_pass_ratio = (TH1D *)histPass.Clone("h_pass_tot_");
  h_pass_ratio->Sumw2();
  h_pass_ratio->Divide(h_pass_tot_);  

  //  TH1D *h_fail_tot; 
  THStack *h_fail_tot = new THStack("h_fail_tot",""); 
  TH1D *h_fail_tot_ = NULL;
  TH1D *h_fail_ratio = NULL;
  if(h_fail_bkg1)
    h_fail_tot->Add(h_fail_bkg1);
  if(h_fail_bkg2)
    h_fail_tot->Add(h_fail_bkg2);
  if(h_fail_signal)
    h_fail_tot->Add(h_fail_signal);
  
  h_fail_tot_ = (TH1D *)h_fail_bkg2->Clone("h_fail_tot_");
  h_fail_tot_->Sumw2();
  if(h_fail_bkg1)
    h_fail_tot_->Add(h_fail_bkg2);
  if(h_fail_signal)
    h_fail_tot_->Add(h_fail_signal);


  h_fail_ratio = (TH1D *)histFail.Clone("h_fail_tot_");
  h_fail_ratio->Sumw2();
  h_fail_ratio->Divide(h_fail_tot_);  
  

  cout<<"****after fail histos**"<<endl;

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

  cout<<"qui5"<<endl;
  
  {
    std::cout<<"plotting data and model prefit via TH1"<<std::endl;
    
    TCanvas* cpass_bf_mt = new TCanvas("cpass_bf_mt","cpass_bf_mt",800,1200);
    //    cpass_bf_mt->SetWindowPosition(cpass_bf_mt->GetWindowTopX()+cpass_bf_mt->GetBorderSize()+800,0);
    //cpass_bf_mt->SetTickx(1);
    //cpass_bf_mt->SetTicky(1);
    TPad *pad1pass_bf_mt = NULL; TPad *pad2pass_bf_mt = NULL;
    pad1pass_bf_mt = new TPad("pad1pass_bf_mt","",0.01,0.41,0.99,1,0);
    pad2pass_bf_mt = new TPad("pad2pass_bf_mt","Data-Exp",0.01,0.,0.99,0.4);
    pad2pass_bf_mt->SetLeftMargin(0.059);
    pad2pass_bf_mt->SetRightMargin(0.02);
    pad1pass_bf_mt->SetRightMargin(0.02);
    pad1pass_bf_mt->SetLeftMargin(0.059);
    pad1pass_bf_mt->Draw(); pad2pass_bf_mt->Draw();

    pad1pass_bf_mt->cd();
    
    cout<<"h_pass_tot "<<h_pass_tot<<endl;

    if(h_pass_tot->GetMaximum()<histPass.GetMaximum())
      h_pass_tot->SetMaximum(histPass.GetMaximum()*1.1);
    h_pass_tot->Draw("hist");
  
    cout<<"qui51"<<endl;
    histPass.Draw("Esame");

    //    h_pass_tot->GetXaxis()->SetTitle("Top Mass [GeV/c]"); //crash...figure out why

    //    cpass_bf_mt->Update();

    cout<<"qui52"<<endl;
    legend->Draw();

    pad2pass_bf_mt->cd();

    h_pass_ratio->Draw("E");
    a->DrawLine(fFitLo,1.,fFitHi,1.);

    fout->cd();
    histPass.Write();
    if(h_pass_bkg1)
      h_pass_bkg1->Write();
    if(h_pass_bkg2)
      h_pass_bkg2->Write();
    if(h_pass_signal)
      h_pass_signal->Write();

    cout<<"qui53"<<endl;

    string str_tmp = fOutputDir+"/TemplatesandData_Passsample_beforefit_mt.png"; 
    cpass_bf_mt->SaveAs(str_tmp.c_str());
    // str_tmp = fOutputDir+"/TemplatesandData_Passsample_beforefit_mt.C"; 
    // cpass_bf_mt->SaveAs(str_tmp.c_str());

 

    if(h_pass_signal_n && h_pass_bkg1_n && h_pass_bkg2_n){
      TCanvas* cpass_bf_norm_mt = new TCanvas("cpass_bf_norm_mt","cpass_bf_norm_mt",800,1200);
      cpass_bf_norm_mt->SetWindowPosition(cpass_bf_norm_mt->GetWindowTopX()+cpass_bf_norm_mt->GetBorderSize()+800,0);
      cpass_bf_norm_mt->SetTickx(1);
      cpass_bf_norm_mt->SetTicky(1);
      h_pass_signal_n->Draw("hist");
      h_pass_bkg1_n->Draw("samehist");
      h_pass_bkg2_n->Draw("samehist");
      histPass.Draw("Esame");
      
      legend->Draw();
    cout<<"qui54"<<endl;      
      str_tmp = fOutputDir+"/TemplatesandData_Passsample_beforefit_norm_mt.png"; 
      cpass_bf_norm_mt->SaveAs(str_tmp.c_str());
      // str_tmp = fOutputDir+"/TemplatesandData_Passsample_beforefit_norm_mt.C"; 
      // cpass_bf_norm_mt->SaveAs(str_tmp.c_str());
    }

    cout<<"qui6"<<endl;

    TCanvas* cfail_bf_mt; 
    cfail_bf_mt = new TCanvas("cfail_bf_mt","cfail_bf_mt",800,1200);
    TPad *pad1fail_bf_mt = NULL; TPad *pad2fail_bf_mt = NULL;
    pad1fail_bf_mt = new TPad("pad1fail_bf_mt","",0.01,0.41,0.99,1,0);
    pad2fail_bf_mt = new TPad("pad2fail_bf_mt","Data-Exp",0.01,0.,0.99,0.4);
    pad2fail_bf_mt->SetLeftMargin(0.059);
    pad2fail_bf_mt->SetRightMargin(0.02);
    pad1fail_bf_mt->SetRightMargin(0.02);
    pad1fail_bf_mt->SetLeftMargin(0.059);
    pad1fail_bf_mt->Draw(); pad2fail_bf_mt->Draw();

    pad1fail_bf_mt->cd();
    if(h_fail_tot->GetMaximum()<histFail.GetMaximum())
      h_fail_tot->SetMaximum(histFail.GetMaximum()*1.1);
    h_fail_tot->Draw("hist");      
    histFail.Draw("Esame");
    
    legend->Draw();
    
    pad2fail_bf_mt->cd();
    
    h_fail_ratio->Draw("E");
    a->DrawLine(fFitLo,1.,fFitHi,1.);
    
    
    fout->cd();
    histFail.Write();
    if(h_fail_bkg1)
      h_fail_bkg1->Write();
    if(h_fail_bkg2)
      h_fail_bkg2->Write();
    if(h_fail_signal)
      h_fail_signal->Write();
    
      
    str_tmp = fOutputDir+"/TemplatesandData_Failsample_beforefit_mt.png"; 
    cfail_bf_mt->SaveAs(str_tmp.c_str());
    // str_tmp = fOutputDir+"/TemplatesandData_Failsample_beforefit_mt.C"; 
    // cfail_bf_mt->SaveAs(str_tmp.c_str());
    
    
    
    
    if(h_fail_signal_n && h_fail_bkg1_n && h_fail_bkg2_n){
      TCanvas* cfail_bf_norm_mt = new TCanvas("cfail_bf_norm_mt","cfail_bf_norm_mt",800,1200);
      cfail_bf_norm_mt->SetWindowPosition(cfail_bf_norm_mt->GetWindowTopX()+cfail_bf_norm_mt->GetBorderSize()+800,0);
      cfail_bf_norm_mt->SetTickx(1);
      cfail_bf_norm_mt->SetTicky(1);
      
      h_fail_signal_n->Draw("hist");
      h_fail_bkg1_n->Draw("samehist");
      h_fail_bkg2_n->Draw("samehist");
      histFail.Draw("Esame");
      
      legend->Draw();
      
      str_tmp = fOutputDir+"/TemplatesandData_Failsample_beforefit_norm_mt.png"; 
      cfail_bf_norm_mt->SaveAs(str_tmp.c_str());
      // str_tmp = fOutputDir+"/TemplatesandData_Failsample_beforefit_norm_mt.C"; 
      // cfail_bf_norm_mt->SaveAs(str_tmp.c_str());
    }
      
    
  }
  
  
  
  
  {
    std::cout<<"plotting data and model prefit"<<std::endl;
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

    string str_tmp = fOutputDir+"/TemplatesandData_Passsample_beforefit.png"; 
    cpass_bf->SaveAs(str_tmp.c_str());
    // str_tmp = fOutputDir+"/TemplatesandData_Passsample_beforefit.C"; 
    // cpass_bf->SaveAs(str_tmp.c_str());



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
    
    str_tmp = fOutputDir+"/TemplatesandData_Failsample_beforefit.png"; 
    cfail_bf->SaveAs(str_tmp.c_str());
    // str_tmp = fOutputDir+"/TemplatesandData_Failsample_beforefit.C"; 
    // cfail_bf->SaveAs(str_tmp.c_str());
    



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
  
    

  }

  

  
  RooFitResult *fitResult=0;
  std::cout<<"before fit"<<std::endl;
  fitResult = totalPdf.fitTo(*dataCombined,
                             RooFit::Extended(),
                             //RooFit::Strategy(2),
                             //RooFit::Minos(RooArgSet(eff)),
                             RooFit::NumCPU(4),
                             RooFit::Save());
  
  std::cout<<"after fit"<<std::endl;
  
  // Refit w/o MINOS if MINOS errors are strange...
  if((fabs(eff.getErrorLo())<5e-5) || (eff.getErrorHi()<5e-5))
    fitResult = totalPdf.fitTo(*dataCombined, RooFit::Extended(), RooFit::Strategy(1), RooFit::Save());

  if((SIGFR12MASK & 0x001) == 0x001){
    resEff  = eff.getVal(); 
    resErrl = fabs(eff.getErrorLo());
    resErrh = eff.getErrorHi();
  }
  else{
    resEff  = 0.;
    resErrl  = 0.;
    resErrh  = 0.;
  }

  //  float Fr2err = FR1.getVal()*FR2overFR1.getVal()*TMath::Sqrt((FR2overFR1.getPropagatedError(*fitResult)*FR2overFR1.getPropagatedError(*fitResult))/(FR2overFR1.getVal()*FR2overFR1.getVal())+(FR1.getPropagatedError(*fitResult)*FR1.getPropagatedError(*fitResult))/(FR1.getVal()*FR1.getVal()));




  if(h_pass_signal_af)
    h_pass_signal_af->Scale(NsigPass->getVal()/h_pass_signal_af->Integral());
  if(h_pass_bkg1_af)
    h_pass_bkg1_af->Scale(Nbkg1Pass->getVal()/h_pass_bkg1_af->Integral());
  if(h_pass_bkg2_af)
    h_pass_bkg2_af->Scale(Nbkg2Pass->getVal()/h_pass_bkg2_af->Integral());
  THStack *h_pass_tot_af = new THStack("h_pass_tot_af",""); 
  if(h_pass_bkg1_af)
    h_pass_tot_af->Add(h_pass_bkg1_af);
  if(h_pass_bkg2_af)
    h_pass_tot_af->Add(h_pass_bkg2_af);
  if(h_pass_signal_af)
    h_pass_tot_af->Add(h_pass_signal_af);
  TH1D *h_pass_tot_af_ = (TH1D *)h_pass_bkg2_af->Clone("h_pass_tot_af_");
  h_pass_tot_af_->Sumw2();
  if(h_pass_bkg1_af)
    h_pass_tot_af_->Add(h_pass_bkg1_af);
  if(h_pass_signal_af)
    h_pass_tot_af_->Add(h_pass_signal_af);

  TH1D *h_pass_ratio_af = (TH1D *)histPass.Clone("h_pass_tot_af_");
  TH1D *h_pass_diff_af = (TH1D *)histPass.Clone("h_pass_tot_af_");
  h_pass_ratio_af->Sumw2();
  h_pass_diff_af->Sumw2();
  h_pass_ratio_af->Divide(h_pass_tot_af_);  
  h_pass_diff_af->Add(h_pass_tot_af_,-1);  

  // for(int ib=0; ib<=h_pass_ratio_af->GetNbinsX(); ib++){
  //   cout<<"bincenter "<<h_pass_ratio_af->GetBinCenter(ib)<<endl;
  //   // double errreldata = histPass.GetBinError(ib)/histPass.GetBinContent(ib);
  //   // double errrelmc = h_pass_tot_af_->GetBinError(ib)/h_pass_tot_af_->GetBinContent(ib);
  //   // double errrel = TMath::Sqrt(errreldata*errreldata+errrelmc*errrelmv);
  //   //    h_pass_ratio_af->SetBinError(errrel*h_pass_ratio_af->GetBinContent(ib));
  //   cout<<"h_pass_ratio_af bc, berr "<<h_pass_ratio_af->GetBinContent(ib)<<" "<<h_pass_ratio_af->GetBinError(ib)<<endl;
  //  cout<<"histPass bc, berr "<<histPass.GetBinContent(ib)<<" "<<histPass.GetBinError(ib)<<endl;
  //  cout<<"h_pass_tot_af_ bc, berr "<<h_pass_tot_af_->GetBinContent(ib)<<" "<<h_pass_tot_af_->GetBinError(ib)<<endl;   
  // }


  if(h_fail_signal_af)
    h_fail_signal_af->Scale(NsigFail->getVal()/h_fail_signal_af->Integral());
  if(h_fail_bkg1_af)
    h_fail_bkg1_af->Scale(Nbkg1Fail->getVal()/h_fail_bkg1_af->Integral());
  if(h_fail_bkg2_af)
    h_fail_bkg2_af->Scale(Nbkg2Fail->getVal()/h_fail_bkg2_af->Integral());
  THStack *h_fail_tot_af = new THStack("h_fail_tot_af",""); 
  if(h_fail_bkg1_af)
    h_fail_tot_af->Add(h_fail_bkg1_af);
  if(h_fail_bkg2_af)
    h_fail_tot_af->Add(h_fail_bkg2_af);
  if(h_fail_signal_af)
    h_fail_tot_af->Add(h_fail_signal_af);
  TH1D *h_fail_tot_af_ = (TH1D *)h_fail_bkg2_af->Clone("h_fail_tot_af_");
  h_fail_tot_af_->Sumw2();
  if(h_fail_bkg1_af)
    h_fail_tot_af_->Add(h_fail_bkg1_af);
  if(h_fail_signal_af)
    h_fail_tot_af_->Add(h_fail_signal_af);

  TH1D *h_fail_ratio_af = (TH1D *)histFail.Clone("h_fail_tot_af_");
  TH1D *h_fail_diff_af = (TH1D *)histFail.Clone("h_fail_tot_af_");
  h_fail_ratio_af->Sumw2();
  h_fail_diff_af->Sumw2();
  h_fail_ratio_af->Divide(h_fail_tot_af_);  
  h_fail_diff_af->Add(h_fail_tot_af_,-1);  




  //plot

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
  // str_tmp = fOutputDir+"/TemplatesandData_Passsample_afterfit.C"; 
  // cpass->SaveAs(str_tmp.c_str());

  cout<<"qui10"<<endl;

  TCanvas *cpass_af_mt = MakeCanvas("cpass_af_mt","cpass_af_mt",800,1200);
  TPad *pad1pass_af_mt = NULL; TPad *pad2pass_af_mt = NULL;
  pad1pass_af_mt = new TPad("pad1pass_af_mt","",0.05,0.41,0.99,1,0);
  pad2pass_af_mt = new TPad("pad2pass_af_mt","Data-Exp",0.05,0.,0.99,0.4);
  pad2pass_af_mt->SetLeftMargin(0.059);
  pad2pass_af_mt->SetRightMargin(0.02);
  pad1pass_af_mt->SetRightMargin(0.02);
  pad1pass_af_mt->SetLeftMargin(0.059);
  pad1pass_af_mt->Draw(); pad2pass_af_mt->Draw();
  
  pad1pass_af_mt->cd();
  if(h_pass_tot_af->GetMaximum()<histPass.GetMaximum())
    h_pass_tot_af->SetMaximum(histPass.GetMaximum()*1.1);
  h_pass_tot_af->Draw("hist");
  
  cout<<"qui11"<<endl;

  histPass.Draw("Esame");

  cout<<"qui12"<<endl;
  legend->Draw();
  
  pad2pass_af_mt->cd();
  
  h_pass_ratio_af->SetMaximum(2);
  h_pass_ratio_af->SetMinimum(2);
  h_pass_ratio_af->Draw("E");
  //h_pass_diff_af->Draw("E");
  a->DrawLine(fFitLo,1.,fFitHi,1.);
  
  cout<<"qui13"<<endl;

  fout->cd();
  //histPass->Write();
  if(h_pass_bkg1_af)
    h_pass_bkg1_af->Write();
  if(h_pass_bkg2_af)
    h_pass_bkg2_af->Write();
  if(h_pass_signal_af)
    h_pass_signal_af->Write();

  cout<<"qui14"<<endl;

  str_tmp = fOutputDir+"/TemplatesandData_Passsample_afterfit_mt.png"; 
  cpass_af_mt->SaveAs(str_tmp.c_str());






  TCanvas *cfail; 
  
  cfail = MakeCanvas("cfail","cfail",800,1200); 
  cfail->SetWindowPosition(cfail->GetWindowTopX()+cfail->GetBorderSize()+800,0);
  cfail->SetTickx(1);
  cfail->SetTicky(1);  
  mframeFail->Draw();

  legend->Draw();
  str_tmp = fOutputDir+"/TemplatesandData_Failsample_afterfit.png"; 
  cfail->SaveAs(str_tmp.c_str());
  // str_tmp = fOutputDir+"/TemplatesandData_Failsample_afterfit.C"; 
  // cfail->SaveAs(str_tmp.c_str());
  
  
  TCanvas *cfail_af_mt = MakeCanvas("cfail_af_mt","cfail_af_mt",800,1200);
  TPad *pad1fail_af_mt = NULL; TPad *pad2fail_af_mt = NULL;
  pad1fail_af_mt = new TPad("pad1fail_af_mt","",0.05,0.41,0.99,1,0);
  pad2fail_af_mt = new TPad("pad2fail_af_mt","Data-Exp",0.05,0.,0.99,0.4);
  pad2fail_af_mt->SetLeftMargin(0.059);
  pad2fail_af_mt->SetRightMargin(0.02);
  pad1fail_af_mt->SetRightMargin(0.02);
  pad1fail_af_mt->SetLeftMargin(0.059);
  pad1fail_af_mt->Draw(); pad2fail_af_mt->Draw();
  
  pad1fail_af_mt->cd();
  if(h_fail_tot_af->GetMaximum()<histFail.GetMaximum())
    h_fail_tot_af->SetMaximum(histFail.GetMaximum()*1.1);
  h_fail_tot_af->Draw("hist");
  
  histFail.Draw("Esame");

  
  legend->Draw();
  
  pad2fail_af_mt->cd();
  
  h_fail_ratio_af->Draw("E");
    //h_fail_diff_af->Draw("E");
  a->DrawLine(fFitLo,1.,fFitHi,1.);
  
  fout->cd();
  //histFail->Write();
  if(h_fail_bkg1_af)
    h_fail_bkg1_af->Write();
  if(h_fail_bkg2_af)
    h_fail_bkg2_af->Write();
  if(h_fail_signal_af)
    h_fail_signal_af->Write();
  

  str_tmp = fOutputDir+"/TemplatesandData_Failsample_afterfit_mt.png"; 
  cfail_af_mt->SaveAs(str_tmp.c_str());
  




  

  
  std::cout<<"++++++++++Fit++++++++++"<<std::endl;
  if((SIGFR12MASK & 0x001) == 0x001)
    std::cout<<"Eff "<<resEff<<" - "<<resErrl<<" + "<<resErrh<<std::endl;
  if((SIGFR12MASK & 0x010) == 0x010)
    std::cout<<"FR1 "<<FR1.getVal()<<" - "<<fabs(FR1.getErrorLo())<<" + "<<FR1.getErrorHi()<<std::endl;
  //    std::cout<<"FR2 "<<FR2overFR1.getVal()*FR1.getVal()<<" + "<<Fr2err<<" - "<<Fr2err<<std::endl;
  if((SIGFR12MASK & 0x100) == 0x100){
    std::cout<<"FR2 "<<FR2.getVal()<<" - "<<fabs(FR2.getErrorLo())<<" + "<<FR2.getErrorHi()<<std::endl;
    std::cout<<"Nbkg2 "<<Nbkg2.getVal()<<" - "<<fabs(Nbkg2.getErrorLo())<<" + "<<Nbkg2.getErrorHi()<<std::endl;
  }
  
  
  std::cout<<"++++++++++Expect++++++++++"<<std::endl; //compute errors
  if((npass_sig+nfail_sig) > 0)
    std::cout<<"Eff "<<npass_sig/(npass_sig+nfail_sig)<<std::endl;
  if((npass_bkg1+nfail_bkg1) > 0)
      std::cout<<"FR1 "<<npass_bkg1/(npass_bkg1+nfail_bkg1)<<std::endl;
  if((npass_bkg2+nfail_bkg2) > 0)
    std::cout<<"FR2 "<<npass_bkg2/(npass_bkg2+nfail_bkg2)<<std::endl;
  
  std::cout<<"++++++++++Fit (extras)++++++++++"<<std::endl;
  std::cout<<"#chi^{2}/dof (Pass)"<<mframePass->chiSquare()<<std::endl;
  std::cout<<"#chi^{2}/dof (Fail)"<<mframeFail->chiSquare()<<std::endl;

  std::cout<<"Extra"<<std::endl;
  //    std::cout<<"FR2overFR1.getVal() "<<FR2overFR1.getVal()<<std::endl;
  std::cout<<"*********Postfit Total Yields*********"<<std::endl;
  std::cout<<"Yields (data) "<<(int)passTree->GetEntries()+(int)failTree->GetEntries()<<std::endl;
  if(NsigPass && NsigFail)
    std::cout<<"Nsig "<<NsigPass->getVal()+NsigFail->getVal()<<" +/- "<<TMath::Sqrt(NsigPass->getPropagatedError(*fitResult)*NsigPass->getPropagatedError(*fitResult)+NsigFail->getPropagatedError(*fitResult)*NsigFail->getPropagatedError(*fitResult))<<std::endl;
  if(Nbkg1Pass && Nbkg1Fail)
    std::cout<<"Nbkg1 "<<Nbkg1Pass->getVal()+Nbkg1Fail->getVal()<<" +/- "<<TMath::Sqrt(Nbkg1Pass->getPropagatedError(*fitResult)*Nbkg1Pass->getPropagatedError(*fitResult)+Nbkg1Fail->getPropagatedError(*fitResult)*Nbkg1Fail->getPropagatedError(*fitResult))<<std::endl;
  if(Nbkg2Pass && Nbkg2Fail)
    std::cout<<"Nbkg2 "<<Nbkg2Pass->getVal()+Nbkg2Fail->getVal()<<" +/- "<<TMath::Sqrt(Nbkg2Pass->getPropagatedError(*fitResult)*Nbkg2Pass->getPropagatedError(*fitResult)+Nbkg2Pass->getPropagatedError(*fitResult)*Nbkg2Pass->getPropagatedError(*fitResult))<<std::endl;    


  
  // float errFitratio = TMath::Sqrt(Nbkg2Pass->getPropagatedError(*fitResult)*Nbkg2Pass->getPropagatedError(*fitResult)/(Nbkg2Pass->getVal()*Nbkg2Pass->getVal())+Nbkg1Pass->getPropagatedError(*fitResult)*Nbkg1Pass->getPropagatedError(*fitResult)/(Nbkg1Pass->getVal()*Nbkg1Pass->getVal()))*(Nbkg2Pass->getVal()/Nbkg1Pass->getVal());

  std::cout<<"*********Fitted Pass**********"<<std::endl;
  std::cout<<"Yields (data) "<<(int)passTree->GetEntries()<<std::endl;
  if(NsigPass)
    std::cout<<"Nsig "<<NsigPass->getVal()<<" +/- "<<NsigPass->getPropagatedError(*fitResult)<<std::endl;
  if(Nbkg1Pass)
    std::cout<<"Nbkg1 "<<Nbkg1Pass->getVal()<<" +/- "<<Nbkg1Pass->getPropagatedError(*fitResult)<<std::endl;
  if(Nbkg2Pass)
    std::cout<<"Nbkg2 "<<Nbkg2Pass->getVal()<<" +/- "<<Nbkg2Pass->getPropagatedError(*fitResult)<<std::endl;
  //  std::cout<<"Nbkg2/Nbkg1 "<<Nbkg2Pass->getVal()/Nbkg1Pass->getVal()<<" +/- "<<errFitratio<<std::endl;
  //propagate error fraction 
  if(NsigPass && Nbkg1Pass && Nbkg2Pass){
    std::cout<<"FracNsig, FracNbkg1, FracNbkg2 "<<NsigPass->getVal()/npass_data<<" "<<Nbkg1Pass->getVal()/npass_data<<" "<<Nbkg2Pass->getVal()/npass_data<<std::endl;
    std::cout<<"*********EXPECTED FRACMC Pass**********"<<std::endl;
    std::cout<<"Signal, Bkg1, Bkg2 "<<npass_sig/npass_data<<" "<<npass_bkg1/npass_data<<" "<<npass_bkg2/npass_data<<std::endl;
  }

  std::cout<<"*********Fitted Fail**********"<<std::endl;
  std::cout<<"Yields (data) "<<(int)failTree->GetEntries()<<std::endl;
  if(NsigFail)
    std::cout<<"Nsig "<<NsigFail->getVal()<<" +/- "<<NsigFail->getPropagatedError(*fitResult)<<std::endl;
  if(Nbkg1Fail)
    std::cout<<"Nbkg1 "<<Nbkg1Fail->getVal()<<" +/- "<<Nbkg1Fail->getPropagatedError(*fitResult)<<std::endl;
  if(Nbkg2Fail)
    std::cout<<"Nbkg2 "<<Nbkg2Fail->getVal()<<" +/- "<<Nbkg2Fail->getPropagatedError(*fitResult)<<std::endl;
  //propagate error fraction  
    if(NsigFail && Nbkg1Fail && Nbkg2Fail){
      std::cout<<"FracNsig, FracNbkg1, FracNbkg2 "<<NsigFail->getVal()/nfail_data<<" "<<Nbkg1Fail->getVal()/nfail_data<<" "<<Nbkg2Fail->getVal()/nfail_data<<std::endl;
      std::cout<<"*********EXPECTED FRACMC Fail**********"<<std::endl;
      std::cout<<"Signal, Bkg1, Bkg2 "<<nfail_sig/nfail_data<<" "<<nfail_bkg1/nfail_data<<" "<<nfail_bkg2/nfail_data<<std::endl;
    }
  
 

  //
  // Write fit results
  //
  ofstream txtfile;
  char txtfname[100];    
  sprintf(txtfname,"%s/fitres.txt",fOutputDir.c_str());
  txtfile.open(txtfname);
  assert(txtfile.is_open());


  double chisqpasstmp=0;
  for(int ib=0; ib<h_pass_diff_af->GetNbinsX(); ib++){
    double bc = h_pass_diff_af->GetBinContent(ib);
    double berr = h_pass_diff_af->GetBinError(ib);
    if(berr != 0)
      chisqpasstmp += bc*bc/(berr*berr);
  } 


  
  txtfile<<"++++++++++Fit++++++++++"<<std::endl;
  if((SIGFR12MASK & 0x001) == 0x001)
    txtfile<<"Eff "<<resEff<<" - "<<resErrl<<" + "<<resErrh<<std::endl;
  if((SIGFR12MASK & 0x010) == 0x010)
    txtfile<<"FR1 "<<FR1.getVal()<<" - "<<fabs(FR1.getErrorLo())<<" + "<<FR1.getErrorHi()<<std::endl;
  //    txtfile<<"FR2 "<<FR2overFR1.getVal()*FR1.getVal()<<" + "<<Fr2err<<" - "<<Fr2err<<std::endl;
  if((SIGFR12MASK & 0x100) == 0x100)
    txtfile<<"FR2 "<<FR2.getVal()<<" - "<<fabs(FR2.getErrorLo())<<" + "<<FR2.getErrorHi()<<std::endl;
  

  txtfile<<"++++++++++Expect++++++++++"<<std::endl; //compute errors
  if((npass_sig+nfail_sig) > 0)
    txtfile<<"Eff "<<npass_sig/(npass_sig+nfail_sig)<<std::endl;
  if((npass_bkg1+nfail_bkg1) > 0)
    txtfile<<"FR1 "<<npass_bkg1/(npass_bkg1+nfail_bkg1)<<std::endl;
  if((npass_bkg2+nfail_bkg2) > 0)
    txtfile<<"FR2 "<<npass_bkg2/(npass_bkg2+nfail_bkg2)<<std::endl;
  //    txtfile<<"#chi^{2}/dof (Pass) "<<mframePass->chiSquare()<<std::endl; //still need be to interpreted
  //txtfile<<"#chi^{2}/dof (Fail) "<<mframeFail->chiSquare()<<std::endl;
  
  double chisqfailtmp=0;
  for(int ib=0; ib<h_fail_diff_af->GetNbinsX(); ib++){
    double bc = h_fail_diff_af->GetBinContent(ib);
    double berr = h_fail_diff_af->GetBinError(ib);
    if(berr != 0)
      chisqfailtmp += bc*bc/(berr*berr);
  } 
  txtfile<<"#chi^{2}/dof (PassMarco) "<<chisqpasstmp<<std::endl;
  txtfile<<"#chi^{2}/dof (FailMarco) "<<chisqfailtmp<<std::endl;
  

  txtfile<<"++++Other Fitted Values++++"<<endl;
  if((SIGFR12MASK & 0x001) == 0x001)
    txtfile<<"Nsig "<<Nsig.getVal()<<" (expect: "<<Nsigexp<<")"<<endl;
  if((SIGFR12MASK & 0x010) == 0x010)
    txtfile<<"Nbkg1 "<<Nbkg1.getVal()<<" (expect: "<<Nbkg1exp<<")"<<endl;
  if((SIGFR12MASK & 0x100) == 0x100)
    txtfile<<"Nbkg2 "<<Nbkg2.getVal()<<" (expect: "<<Nbkg2exp<<")"<<endl;
    
  //  txtfile<<"FR2overFR1.getVal() "<<FR2overFR1.getVal()<<std::endl;
  txtfile<<"*********Fitted Pass**********"<<std::endl;
  txtfile<<"Yields (data) "<<(int)passTree->GetEntries()<<std::endl;
  if(NsigPass)
    txtfile<<"Nsig "<<NsigPass->getVal()<<" +/- "<<NsigPass->getPropagatedError(*fitResult)<<std::endl;
  if(Nbkg1Pass)
    txtfile<<"Nbkg1 "<<Nbkg1Pass->getVal()<<" +/- "<<Nbkg1Pass->getPropagatedError(*fitResult)<<std::endl;
  if(Nbkg2Pass)
    txtfile<<"Nbkg2 "<<Nbkg2Pass->getVal()<<" +/- "<<Nbkg2Pass->getPropagatedError(*fitResult)<<std::endl;
  //  txtfile<<"Nbkg2/Nbkg1 "<<Nbkg2Pass->getVal()/Nbkg1Pass->getVal()<<" +/- "<<errFitratio<<std::endl;
  //propagate error fraction  
  if(NsigPass && Nbkg1Pass && Nbkg2Pass){
    txtfile<<"FracNsig, FracNbkg1, FracNbkg2 "<<NsigPass->getVal()/npass_data<<" "<<Nbkg1Pass->getVal()/npass_data<<" "<<Nbkg2Pass->getVal()/npass_data<<std::endl;
    txtfile<<"*********EXPECTED FRACMC Pass**********"<<std::endl;
    txtfile<<"Signal, Bkg1, Bkg2 "<<npass_sig/npass_data<<" "<<npass_bkg1/npass_data<<" "<<npass_bkg2/npass_data<<std::endl;
  }

  txtfile<<"*********Fitted Fail**********"<<std::endl;
  txtfile<<"Yields (data) "<<(int)failTree->GetEntries()<<std::endl;
  if(NsigFail)
    txtfile<<"Nsig "<<NsigFail->getVal()<<" +/- "<<NsigFail->getPropagatedError(*fitResult)<<std::endl;
  if(Nbkg1Fail)
    txtfile<<"Nbkg1 "<<Nbkg1Fail->getVal()<<" +/- "<<Nbkg1Fail->getPropagatedError(*fitResult)<<std::endl;
  if(Nbkg2Fail)
    txtfile<<"Nbkg2 "<<Nbkg2Fail->getVal()<<" +/- "<<Nbkg2Fail->getPropagatedError(*fitResult)<<std::endl;
  //propagate error fraction  
  if(NsigFail && Nbkg1Fail && Nbkg2Fail){
    txtfile<<"FracNsig, FracNbkg1, FracNbkg2 "<<NsigFail->getVal()/nfail_data<<" "<<Nbkg1Fail->getVal()/nfail_data<<" "<<Nbkg2Fail->getVal()/nfail_data<<std::endl;
    
    txtfile<<"*********EXPECTED FRACMC Fail**********"<<std::endl;
    txtfile<<"Signal, Bkg1, Bkg2 "<<nfail_sig/nfail_data<<" "<<nfail_bkg1/nfail_data<<" "<<nfail_bkg2/nfail_data<<std::endl;
  }
  

  txtfile<<"&&&&&&&&&&&PREFIT YIELDS&&&&&&&&&&&&&&"<<std::endl;
  
  txtfile<<"**********Fail**********"<<std::endl;
  txtfile<<"h_signal, h_bkg1, h_bkg2 "<<nfail_sig<<" "<<nfail_bkg1<<" "<<nfail_bkg2<<std::endl;
  
  txtfile<<"**********Pass********"<<std::endl;
  txtfile<<"h_signal, h_bkg1, h_bkg2 "<<npass_sig<<" "<<npass_bkg1<<" "<<npass_bkg2<<std::endl;
  txtfile<<"*********TOT**********"<<std::endl;
  txtfile<<"h_signal, h_bkg1, h_bkg2 "<<nfail_sig+npass_sig<<" "<<nfail_bkg1+npass_bkg1<<" "<<nfail_bkg2+npass_bkg2<<std::endl;
  txtfile<<"*********TOTMC**********"<<std::endl;
  txtfile<<"Pass, Fail TOT "<<npass_sig+npass_bkg1+npass_bkg2<<" "<<nfail_sig+nfail_bkg1+nfail_bkg2<<" "<<npass_sig+npass_bkg1+npass_bkg2+nfail_sig+nfail_bkg1+nfail_bkg2<<std::endl;
  txtfile<<"*******Data*******"<<std::endl;
  txtfile<<"Pass, Fail, Tot "<<npass_data<<" "<<nfail_data<<" "<<npass_data+nfail_data<<std::endl;
  txtfile<<"**************"<<std::endl;
    
  txtfile<<"&&&&&&&&&&&POSTFIT YIELDS&&&&&&&&&&&&&&"<<std::endl;
  txtfile<<"Yields (data) "<<(int)passTree->GetEntries()+(int)failTree->GetEntries()<<std::endl;
  if(NsigPass && NsigFail)
    txtfile<<"Nsig "<<NsigPass->getVal()+NsigFail->getVal()<<" +/- "<<TMath::Sqrt(NsigPass->getPropagatedError(*fitResult)*NsigPass->getPropagatedError(*fitResult)+NsigFail->getPropagatedError(*fitResult)*NsigFail->getPropagatedError(*fitResult))<<std::endl;
  if(Nbkg1Pass && Nbkg1Fail)
    txtfile<<"Nbkg1 "<<Nbkg1Pass->getVal()+Nbkg1Fail->getVal()<<" +/- "<<TMath::Sqrt(Nbkg1Pass->getPropagatedError(*fitResult)*Nbkg1Pass->getPropagatedError(*fitResult)+Nbkg1Fail->getPropagatedError(*fitResult)*Nbkg1Fail->getPropagatedError(*fitResult))<<std::endl;
  if(Nbkg2Pass && Nbkg2Fail)
    txtfile<<"Nbkg2 "<<Nbkg2Pass->getVal()+Nbkg2Fail->getVal()<<" +/- "<<TMath::Sqrt(Nbkg2Pass->getPropagatedError(*fitResult)*Nbkg2Pass->getPropagatedError(*fitResult)+Nbkg2Pass->getPropagatedError(*fitResult)*Nbkg2Pass->getPropagatedError(*fitResult))<<std::endl;    
  if(NsigPass && NsigFail && Nbkg1Pass && Nbkg1Fail && Nbkg2Pass && Nbkg2Fail)
    txtfile<<"TOTMC "<<NsigPass->getVal()+NsigFail->getVal()+Nbkg1Pass->getVal()+Nbkg1Fail->getVal()+Nbkg2Pass->getVal()+Nbkg2Fail->getVal()<<endl;

  

  
  txtfile <<"DONE WITH PRINTING MY RESULTS"<< endl;
  txtfile << endl;  txtfile << endl;
  txtfile <<"NOW DUMPING fit results from fitResult->printStream"<< endl;
  fitResult->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile <<"PRINT VALUE"<< endl;
  fitResult->printValue(txtfile);
  txtfile << endl;
  txtfile.close();

  if(fout)
    fout->Close();


  //
  // Clean up
  //

  
  delete modelPass;
  delete modelFail;  
  delete dataCombined;
  delete dataPass;
  delete dataFail;
  if(sigModPass_rdh)
    delete sigModPass_rdh;
  if(bkg1ModPass_rdh)
    delete bkg1ModPass_rdh;
  if(bkg2ModPass_rdh)
    delete bkg2ModPass_rdh;
  if(sigModFail_rdh)
    delete sigModFail_rdh;
  if(bkg1ModFail_rdh)
    delete bkg1ModFail_rdh;
  if(bkg2ModFail_rdh)
    delete bkg2ModFail_rdh;
  if(sigModPass)
    delete sigModPass;
  if(sigModPass2)
    delete sigModPass2;
  if(bkg1ModPass)
    delete bkg1ModPass;  
  if(bkg2ModPass)
    delete bkg2ModPass;
  if(sigModFail)
    delete sigModFail;
  if(sigModFail2)
    delete sigModFail2;
  if(bkg1ModFail)
    delete bkg1ModFail;
  if(bkg2ModFail)
    delete bkg2ModFail;
  if(Nbkg1Pass)
  delete Nbkg1Pass; 
  if(Nbkg2Pass)
    delete Nbkg2Pass; 
  if(NsigPass)
    delete NsigPass;
  if(Nbkg1Fail)
    delete Nbkg1Fail; 
  if(Nbkg2Fail)
    delete Nbkg2Fail; 
  if(NsigFail)
    delete NsigFail;
  delete histfile_signal;  delete histfile_bkg1;  delete histfile_bkg2;
  delete legend;
}

float CEffZFitter::calcmass(TLorentzVector *obj1, TLorentzVector *obj2, TLorentzVector *obj3){
  float invmass = -999.;
  TLorentzVector jetSumLV;

  jetSumLV = (*obj1)+(*obj2)+(*obj3);

  invmass = jetSumLV.Mag();

  return invmass;

}

 
float CEffZFitter::calcwmass(TLorentzVector *obj1, TLorentzVector *obj2){
  //calc the W mass as the mass of the two jets (out of three) w/ lowest csv score                                                                                                                                                            
  float invmass = -999.;
  TLorentzVector jet1Wmass, jet2Wmass, jetSumWmass;

  jet1Wmass = (*obj1);
  jet2Wmass = (*obj2);


  jetSumWmass = jet1Wmass + jet2Wmass;

  invmass = jetSumWmass.Mag();

  return invmass;

}
