#ifndef EFFRTTFITTER_HH
#define EFFRTTFITTER_HH

//================================================================================================
//
// Signal Extraction
//-------------------
//  0: probe counting
//  1: Breit-Wigner convolved with Crystal Ball function
//  2: MC template convolved with Gaussian
//  3: Phil's Crystal Ball based "Voigtian" shape
//  4: Unbinned MC data convolved with Gaussian
//
// Background Model
//------------------
//  0: no background
//  1: exponential model
//  2: erfc*exp model
//  3: double exponential model
//  4: linear*exp model
//  5: quadratic*exp model
//
//________________________________________________________________________________________________

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <TLorentzVector.h>
#include "TLegend.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"

#include "ZSignals.hh"

#include "CPlot.hh"                   // helper class for plots



class TTree;
class TCanvas;
class TGraphAsymmErrors;
class TH1D;
class TH2D;
class RooFitResult;

class EffRTTFitter
{
public:
  EffRTTFitter();
  ~EffRTTFitter();
  
  
 
  void initialize(const std::string infname, const std::string outdir, const std::string temfname1, const std::string temfname2, const std::string temfname3, const double TOPMVACUT, const bool MAKETEMPLATES, const float FITLO, const float FITHI, const float BINSIZEPASS, const float BINSIZEFAIL, const unsigned int SIGNALFAKERATE12MASK, const float FR2SCALEFACTOR, const float FR2ERROR);


  void computeEff();

  void SetMODEL(unsigned int model) { MODEL = model;};  
  void SetVERBOSE(bool verbose) { VERBOSE = verbose;};
    
protected:
  
  void makeBinnedTemplates(const std::string temfname, std::string name); //MT
  
  
  
  
  // void performFit(double &resEff, double &resErrl, double &resErrh,
  // 		  TTree *passTree, TTree *failTree);
  void performFit();

  
  //  float calcmass(TLorentzVector *obj1, TLorentzVector *obj2, TLorentzVector *obj3);  
  //float calcwmass(TLorentzVector *obj1, TLorentzVector *obj2);

  void printresults(bool PrintToFile, bool verbose);
  void CloneTemplates(bool verbose);
  void SetTemplatesColors(TH1D *hMC_1, TH1D *hMC_2, TH1D *hMC_3, bool verbose);
  void PlotDataVsMC(TH1D *hMC_1, TH1D *hMC_2, TH1D *hMC_3, TH1D hData, bool doLogY, std::string picname, bool verbose);
  void RooPlotDataVsMC(RooRealVar m, TLegend *legend, std::string pictnamepass, std::string pictnamefail, bool verbose);
  double CalcChisq(TH1D *hMC_1, TH1D *hMC_2, TH1D *hMC_3, TH1D hData, bool verbose);
  void printRooFitresults(bool verbose);
  void CreateOutFile(bool verbose);
  void WriteHistos(bool verbose);


  ///// data members /////
  
  bool fIsInitialized;
  bool maketemplates;  
  
  //output root file where templates are stored
  TFile *fout;

  int fSigPass, fBkgPass, fSigFail, fBkgFail;
  
  double fLo, fHi;        // signal extraction fit variable window  
  double fFitLo, fFitHi;  // fit variable window
  
  double topmvacut;

  float BIN_SIZE_PASS;
  float BIN_SIZE_FAIL;
  unsigned int SIGFR12MASK;
  float FR2SF; //SF w/ expectations
  float FR2ERR; //RANGE AROUND FR2SF*EXP allowed in the fit

  unsigned int MODEL; //=1 fit with templates from MC, =2 smear MC with Gaussian (CMCTemplateConvGaussian)

  bool VERBOSE;
  
  // output directory for results
  std::string fOutputDir;
      
  TTree *fPassTree, *fFailTree;

  RooFormulaVar *Nbkg1Pass; 
  RooFormulaVar *Nbkg1Fail; 
  RooFormulaVar *Nbkg2Pass; 
  RooFormulaVar *Nbkg2Fail;
  RooFormulaVar  *NsigPass;
  RooFormulaVar  *NsigFail;

  RooAbsData *dataCombined; //Pass+fail data

  //files containing template (created by makebinnedtemplates)
  TFile *histfile_signal;
  TFile *histfile_bkg1;  
  TFile *histfile_bkg2;

  //data histos pass/fail
  TH1D *h_pass_data;
  TH1D *h_fail_data;
  //Templates
  TH1D *h_pass_signal         ;
  TH1D *h_pass_bkg1  	      ;
  TH1D *h_pass_bkg2  	      ;
  TH1D *h_fail_signal 	      ;
  TH1D *h_fail_bkg1  	      ;
  TH1D *h_fail_bkg2  	      ;
  //Templates normalized to data
  // TH1D *h_pass_signal_n       ;
  // TH1D *h_pass_bkg1_n         ;
  // TH1D *h_pass_bkg2_n  	      ;
  // TH1D *h_fail_signal_n       ;
  // TH1D *h_fail_bkg1_n         ;
  // TH1D *h_fail_bkg2_n  	      ;
  //Templates to fitted values
  TH1D *h_pass_signal_af      ;
  TH1D *h_pass_bkg1_af        ;
  TH1D *h_pass_bkg2_af        ;
  TH1D *h_fail_signal_af      ;
  TH1D *h_fail_bkg1_af        ;
  TH1D *h_fail_bkg2_af        ;


  //Pass/Fail PDF
  RooHistPdf *sigModPass ;
  RooHistPdf *bkg1ModPass;
  RooHistPdf *bkg2ModPass;
  RooHistPdf *sigModFail ;
  RooHistPdf *bkg1ModFail;
  RooHistPdf *bkg2ModFail;

  //Pass/Fail smeared with a gaussian PDF 
  CSignalModel     *sigModFail2 ;
  CSignalModel     *bkg1ModFail2;
  CSignalModel     *bkg2ModFail2;
  CSignalModel     *sigModPass2 ;
  CSignalModel     *bkg1ModPass2;
  CSignalModel     *bkg2ModPass2;


  //Pass/Fail DataHist
  RooDataHist *sigModPass_rdh ;
  RooDataHist *bkg1ModPass_rdh;
  RooDataHist *bkg2ModPass_rdh;  
  RooDataHist *sigModFail_rdh ;
  RooDataHist *bkg1ModFail_rdh;
  RooDataHist *bkg2ModFail_rdh;   


  RooAddPdf *modelPass;
  RooAddPdf *modelFail;

  RooAbsData *dataPass;
  RooAbsData *dataFail;


  RooFitResult *fitResult;



  //Efficiencies
  double Nsigexp; 
  double Nbkg1exp;
  double Nbkg2exp;
  double nsig_obs;
  double nbkg1_obs;
  double nbkg2_obs;


  float eff_sig_exp;
  float eff_fr1_exp;
  float eff_fr2_exp;
  double eff_sig_obs, effl_sig_obs, effh_sig_obs;
  double fr1_sig_obs, fr1l_sig_obs, fr1h_sig_obs;
  double fr2_sig_obs, fr2l_sig_obs, fr2h_sig_obs;


  //Expected yields
  double nsigpass_exp;
  double nbkg1pass_exp; 
  double nbkg2pass_exp;

  double nsigfail_exp; 
  double nbkg1fail_exp; 
  double nbkg2fail_exp; 

  //Observed yields
  double npass_data;
  double nfail_data;



  double nsigpass_obs, nsigpass_err_obs;
  double nbkg1pass_obs, nbkg1pass_err_obs;
  double nbkg2pass_obs, nbkg2pass_err_obs;
  double nsigfail_obs, nsigfail_err_obs;
  double nbkg1fail_obs, nbkg1fail_err_obs;
  double nbkg2fail_obs, nbkg2fail_err_obs;

  double chisqpass, chisqfail;                            //chisq for pass/fail (postfit)

  //
  TH1D *h_fail_diff_af;
  TH1D *h_pass_diff_af;
};

#endif
