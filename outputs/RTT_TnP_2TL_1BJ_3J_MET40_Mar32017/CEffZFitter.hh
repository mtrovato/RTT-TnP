#ifndef CEFFZFITTER_HH
#define CEFFZFITTER_HH

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

class TTree;
class TCanvas;
class TGraphAsymmErrors;
class TH1D;
class TH2D;
class RooFitResult;

class CEffZFitter
{
public:
  CEffZFitter();
  ~CEffZFitter();
  
  
 
  void initialize(const std::string infname, const std::string outdir, const std::string temfname1, const std::string temfname2, const std::string temfname3, const double TOPMVACUT, const bool MAKETEMPLATES, const float FITLO, const float FITHI, const float BINSIZEPASS, const float BINSIZEFAIL, const unsigned int SIGNALFAKERATE12MASK, const float FR2SCALEFACTOR, const float FR2ERROR);


  void computeEff();

  
    
protected:

  void makeBinnedTemplates(const std::string temfname, std::string name); //MT
  
  
  
  
  void performFit(double &resEff, double &resErrl, double &resErrh,
		  TTree *passTree, TTree *failTree);


  
  float calcmass(TLorentzVector *obj1, TLorentzVector *obj2, TLorentzVector *obj3);  
  float calcwmass(TLorentzVector *obj1, TLorentzVector *obj2);
  ///// data members /////
  
  bool fIsInitialized;
  bool maketemplates;  
  

  int fSigPass, fBkgPass, fSigFail, fBkgFail;
  
  double fLo, fHi;        // signal extraction fit variable window  
  double fFitLo, fFitHi;  // fit variable window
  
  double topmvacut;

  float BIN_SIZE_PASS;
  float BIN_SIZE_FAIL;
  unsigned int SIGFR12MASK;
  float FR2SF; //SF w/ expectations
  float FR2ERR; //RANGE AROUND FR2SF*EXP allowed in the fit

  int MODEL; //=1 fit with templates from MC, =2 smear MC

  bool VERBOSE;
  
  // output directory for results
  std::string fOutputDir;
      
  TTree *fPassTree, *fFailTree;
  
  
};

#endif
