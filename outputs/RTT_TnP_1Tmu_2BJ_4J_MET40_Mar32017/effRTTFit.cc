#include "EffRTTFitter.hh"
//#include "CEffPlotter.hh"
#include "KStyle.hh"

#include <iostream>
#include <string>
#include <cstdlib>

int main(int argc, char **argv)
{
  std::cout<<"START: effRTTFit"<<std::endl;
  //--------------------------------------------------------------------------------------------------------------
  // Settings and constants
  //==============================================================================================================

  // handle input arguments
  //MT

  const std::string infname  = argv[1];         // input ROOT file of probes
  const std::string outdir   = argv[2];         // output directory
  const std::string temfname1 = argv[3];        // ROOT file for generating MC-based signal templates
  const std::string temfname2 = argv[4];        // ROOT file for generating MC-based bg1 (Wrong comb tt1l) templates
  const std::string temfname3 = argv[5];        // ROOT file for generating MC-based (nontt1l bg) templates
  const float         topmvacut  = atof(argv[6]); //topmvacut to split inclusive sample in pass and fail  
  const bool         maketemplates  = atoi(argv[7]); //decide whether to make templates or not (they may already exist)
  const float fitLo = atof(argv[8]);
  const float fitHi = atof(argv[9]);
  const float BIN_SIZE_PASS = atof(argv[10]); 
  const float BIN_SIZE_FAIL = atof(argv[11]); 
  const unsigned int SIGNALFAKERATE12MASK = atoll(argv[12]); //0x001 sig, 0x010 bkg1, 0x100 bk2
  const float FR2SCALEFACTOR = atof(argv[13]);
  const float FR2ERROR = atof(argv[14]);
  const unsigned int MODEL = atoll(argv[15]);

 
  // other settings
  //  const double       fitLo = 0;
  //const double       fitHi = 1000;
  
  std::cout<<"****effRTTFit****"<<std::endl;
  std::cout << std::endl;
  std::cout << " <> Processing probes file: " << infname << std::endl;
  std::cout << " outdir: " << outdir << std::endl;
  std::cout << " outdir: " << outdir << std::endl;
  std::cout << " temfname1, temfname2, temfname3: " << temfname1 <<" "<<temfname2<<" "<<temfname3<< std::endl;
  std::cout << "topmvacut "<<topmvacut<<std::endl;
  std::cout << "maketemplates "<<maketemplates<<std::endl;
  std::cout << "fitLo, fitHi "<<fitLo<<" "<<fitHi<<std::endl;
  std::cout << "BIN_SIZE_PASS BIN_SIZE_FAIL "<<BIN_SIZE_PASS<<" "<<BIN_SIZE_FAIL<<" "<<std::endl;
  std::cout << "SIGNALFAKERATE12MASK "<<SIGNALFAKERATE12MASK<<std::endl; 
  std::cout << "FR2SCALEFACTOR "<<FR2SCALEFACTOR<<std::endl;
  std::cout << "FR2ERROR "<<FR2ERROR<<std::endl;
  std::cout << "MODEL "<<MODEL<<std::endl;

  std::cout << std::endl;

  KStyle();  
  
  EffRTTFitter fitter;

  std::cout<<"initialize"<<std::endl;
  fitter.initialize(infname, outdir, temfname1,temfname2,temfname3,
		   topmvacut, maketemplates,
		   fitLo, fitHi, BIN_SIZE_PASS, BIN_SIZE_FAIL, SIGNALFAKERATE12MASK, FR2SCALEFACTOR, FR2ERROR);

  fitter.SetMODEL(2);

  fitter.computeEff();
  
  std::cout << " <> Output saved in " << outdir << std::endl;






  return 0;
}
