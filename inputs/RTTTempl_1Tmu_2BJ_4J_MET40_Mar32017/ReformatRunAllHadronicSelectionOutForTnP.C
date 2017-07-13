#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TBenchmark.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include "CPlot.hh"                   // helper class for plots

#include "LeptonEffUtils.hh"
#include "METEffUtils.hh"
#include "kfactors.hh"
#include "Utils.C"

#endif





// Main macro function
// //--------------------------------------------------------------------------------------------------
void ReformatRunAllHadronicSelectionOutForTnP(const char* input, const char* outfilename, const char* LUMI_c, const char* NORMTODATA_c, string processforTnP="", bool isData=false, bool isSemiLepSelection=true) 
{
  ////singlelepdata =>    //(0= 0 lep; 1 = 1 muo, 2=1 ele, 3=2 muo; 4=2 ele; 5=1ele, 1muo will be decide below according to the filename

  const double LUMI = atof(LUMI_c);
  double NORMTODATA = atof(NORMTODATA_c); //if you wanna print the MC normalized to data 
  unsigned int singlelepdata = 1; //by default it's the 1 mu

  bool useflatkfact=true;
  

  int semilep = -1; //for tt events request at gen level 1 lep (semilep=1) or !=1 lep (semilep=0); for other events no request (semilep=-1)



  cout<<"input "<<input<<endl;
  cout<<"outfilename "<<outfilename<<endl;
  cout<<"LUMI "<<LUMI<<endl;
  cout<<"NORMTODATA "<<NORMTODATA<<endl;
  cout<<"processforTnP "<<processforTnP<<endl;
  cout<<"isData "<<isData<<endl;
  cout<<"isSemiLepSelection "<<isSemiLepSelection<<endl;


  bool verbose = false;


  gBenchmark->Start("ReformatRunAllHadronicSelectionOutForTnP");

  vector<string> infilenames;
  
  // 
  // parse input file
  //  
  ifstream ifs;
  ifs.open(input); 
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) { infilenames.push_back(line); }
  ifs.close();
  
  TTree::SetMaxTreeSize(kMaxLong64);

  //
  //Save Yields in a Yields.txt
  //
  std::ofstream yieldfile;
  string yieldstring = "./Yields.txt";
  yieldfile.open (yieldstring.c_str(), std::ofstream::out | std::ofstream::app);


  //--------------------------------------------------------------------------------------------------------------
  // Settings
  //==============================================================================================================
  
  
  const bool useMuonMETturnon = true; //false is we want to use ele parameterization



  //--------------------------------------------------------------------------------------------------------------
  // variables to read in bacon bits
  //============================================================================================================== 
  //TO BE CHANGED ACCORDING TO runALLHadronicPreselection.cc
  unsigned int runNum, lumiSec, evtNum;                          // event ID
  unsigned int metfilter;                                        // MET filter bits
  unsigned int npv;                                              // number of PV
  float npu;                                                     // mean expected PU
  unsigned int njets;                                            // jet multiplicity
  unsigned int nbjetsL, nbjetsM, nbjetsT;                        // b-tag multiplicity
  unsigned int nele, nmuo, ntau, nlep;                           // loose/tight lepton multiplicity depending on LOOSE FLAG = loose/tight
  unsigned int nele_oth, nmuo_oth, nlep_oth;                     // tight/loose lepton multiplicity depending on LOOSE FLAG = loose/tight
  unsigned int nlpho;                                            // loose pho (only if domonojet)
  float scale1fb;                                                // cross section scale factor per 1/fb
  float lhew;                                                    // LHE weight
  float ttWgt;                                                   // event weight for top pT reweighting

  float btag0Wgt, btag0WgtUp, btag0WgtDown;                      // event weight from b-tag scale factors
  float btag1Wgt, btag1WgtUp, btag1WgtDown;
  float btag2Wgt, btag2WgtUp, btag2WgtDown;

  float pfmetraw, pfmetphiraw, dphijetmetraw;    // raw PF MET
  float pfmet, pfmetphi, dphijetmet;                // PF MET
  float dphijetmet4;                                             //up to 4 jets used
  float calomet, calometphi;                                     // Calo MET
  float recoil;

  //additional lep var
  TClonesArray electron("TLorentzVector"); 
  TClonesArray muon("TLorentzVector");  
  TClonesArray tau("TLorentzVector");

  int lep1charge=-999; int lep2charge=-999;

  //trigger
  bool passtrigger;


  // generator level info
  float genVPt;
  float qcdScaleWgtUp, qcdScaleWgtDown;
  int             genId11, genId12, genId13;
  int             genId21, genId22, genId23;
  TLorentzVector *genpar11=0, *genpar12=0, *genpar13=0;
  TLorentzVector *genpar21=0, *genpar22=0, *genpar23=0;
  bool ishf=0;
  unsigned int nlepgenfortt=0; //not produced yet at .cc 
  unsigned int numb=0; //not produced yet at .cc

  // resolved or inclusive variables
  float           res_wmva, res_topmva;
  float           res_chisq, res_prob, res_cost, res_fitmass, res_fitmassW;
  float res_bdt_dphij1b_out, res_bdt_dphij2b_out, res_bdt_drj1b_out, res_bdt_drj2b_out;

  float           res_wmva2, res_topmva2;
  float           res_chisq2, res_prob2, res_cost2, res_fitmass2, res_fitmassW2;
  float res_bdt_dphij4b_out, res_bdt_dphij5b_out, res_bdt_drj4b_out, res_bdt_drj5b_out;


  float           res_jet1csv,     res_jet2csv,     res_jet3csv;
  float           res_jet4csv,     res_jet5csv,     res_jet6csv;
  float           res_jet1qgid,    res_jet2qgid,    res_jet3qgid;
  float           res_jet4qgid,    res_jet5qgid,    res_jet6qgid;
  TVector2       *res_jet1pull=0, *res_jet2pull=0, *res_jet3pull=0;
  TLorentzVector *res_jet1=0,     *res_jet2=0,     *res_jet3=0;
  TLorentzVector *res_jet4=0,     *res_jet5=0,     *res_jet6=0;

  unsigned int cutflow;                                            //keep track of which cut have been passed
 //END: TO BE CHANGED ACCORDING TO runAllHadronicPreselection.cc

  //VARIABLE TO BE COMPUTED IN THIS MACRO
  float evtWeight;
  float isHadrTopHighestMVA, isHadrTopHighestMVA2, isHadrTopHighestMVA_2, isHadrTopHighestMVA2_2; //0,1 is semileptonic ttbar event, -1 not semileptonic ttbar event; 0: highest (or second highest) MVA score is not top hadronic, 1 is top hadronic. N.B: in event with 2 hadr tops I should sum isHadrTopHighestMVA+isHadrTopHighestMVA_2 to find out if jet triplet with the highest matched to first hadronic top; and isHadrTopHighestMV2A+isHadrTopHighestMVA2_2 to find out if jet triplet with the secon score matched to one hadornic top

 float dphitop1top2; //deltaphi between the first 3-jet system and the second 3-jet system: they are the two hadronic top if tagger is used
 float topmass, topmass2, wmass, wmass2; //take all the jets in the event above threshold and build the mass from the combination of jets (3 and 2 respectively) which give the mass closer to the nominal top and W mass respectively
 float toppt, toppt2;


  float pfmt;
  //END: VARIABLE TO BE COMPUTED IN THIS MACRO


  //additional variables for support
  TClonesArray *electron_arr  = 0; TBranch *electron_br = 0; //don't remove the 0 or it wil crash
  TClonesArray *muon_arr = 0; TBranch *muon_br = 0;
  TClonesArray *tau_arr = 0; TBranch *tau_br = 0;

  float btagWgt, btagWgtp, btagWgtm;             // event weight from b-tag scale factors

  vector<TLorentzVector> ele;
  vector<TLorentzVector> muo;
  vector<TLorentzVector> ta;




  TFile *outFile = new TFile(outfilename,"RECREATE");
  TH1D hTotalEvents("TotalEvents","TotalEvents",1,-10,10);
  TTree *outTree = new TTree("Events","Events");

  outTree->Branch("runNum",         &runNum,         "runNum/i");
  outTree->Branch("lumiSec",        &lumiSec,        "lumiSec/i");
  outTree->Branch("evtNum",         &evtNum,         "evtNum/i");
  // resolved variables
  outTree->Branch("res_wmva",     &res_wmva,     "res_wmva/F");
  outTree->Branch("res_topmva",   &res_topmva,   "res_topmva/F");
  outTree->Branch("topmass",   &topmass,   "topmass/F");
  //END:   //TO BE CHANGED ACCORDING TO runALLHadronicPreselection.cc
  outTree->Branch("evtWeight", &evtWeight, "evtWeight/F");
  //--------------------------------------------------------------------------------------------------------------
  // Merge input files
  //============================================================================================================== 
  
  TFile *infile=0;
  TTree *intree=0;

  double neventsv=0.;
  int nentries=0;

  //MT Mar 2: I dont to store total events
  //
  // First sum up TotalEvents histogram from all files in order to properly compute 
  //
  // for(unsigned int ifile=0; ifile<infilenames.size(); ifile++) {
  //   cout<<"infilenames["<<ifile<<"]"<<infilenames[ifile]<<endl;
  //   infile = new TFile(infilenames[ifile].c_str());
  //   hTotalEvents.Add((TH1*)infile->Get("TotalEvents"));
  //   infile->Close();
  //   delete infile;
  //   infile=0;
  //   intree=0;
  // }

  //
  // Loop through input files
  //

  //CHANGE ACCORDiNG TO RUNALLHADRONICPRESELECTIOn
  for(UInt_t ifile=0; ifile<infilenames.size(); ifile++) {
    cout << "Adding " << infilenames[ifile] << endl;
    infile = new TFile(infilenames[ifile].c_str()); assert(infile);
    intree = (TTree*)infile->Get("Events"); assert(intree);


    //    cout<<"processforTnP "<<processforTnP<<" infilenames[ifile] "<<infilenames[ifile].find("TT_powheg")<<endl;

    if((processforTnP=="Signal" || processforTnP=="Bkg") && infilenames[ifile].find("TT_") != string::npos)
      semilep = 1;
    else if(infilenames[ifile].find("TT_") != string::npos)
      semilep = 0;
    else
      semilep = -1;
    
    cout<<"semilep "<<semilep<<endl;

    if(infilenames[ifile].find("1Tmu") != string::npos)
      singlelepdata = 1;
    if(infilenames[ifile].find("2Tmu") != string::npos)
      singlelepdata = 3;
    else if(infilenames[ifile].find("2Te") != string::npos)
      singlelepdata = 4;
    if(infilenames[ifile].find("1Tmu1Tele") != string::npos)
      singlelepdata = 5;
    
    cout<<"singlelepdata "<<singlelepdata<<endl;
    
    intree->SetBranchAddress("runNum",         &runNum                                             );   
     intree->SetBranchAddress("lumiSec",        &lumiSec         			                 );
     intree->SetBranchAddress("evtNum",         &evtNum          			                 );
     intree->SetBranchAddress("metfilter",      &metfilter       			                 );
     intree->SetBranchAddress("npv",            &npv             			                 );
     intree->SetBranchAddress("npu",            &npu             			                 );
     intree->SetBranchAddress("njets",          &njets           			                 );
     intree->SetBranchAddress("nbjetsL",        &nbjetsL         			                 );
     intree->SetBranchAddress("nbjetsM",        &nbjetsM         			                 );
     intree->SetBranchAddress("nbjetsT",        &nbjetsT         			                 );
     intree->SetBranchAddress("nele",           &nele            			                 );
     intree->SetBranchAddress("nmuo",           &nmuo            			                 );
     intree->SetBranchAddress("ntau",           &ntau            			                 );
     intree->SetBranchAddress("nlep",           &nlep            			                 );
     intree->SetBranchAddress("nlep_oth",       &nlep_oth        			                 );
     intree->SetBranchAddress("nele_oth",       &nele_oth        			                 );
     intree->SetBranchAddress("nmuo_oth",       &nmuo_oth        			                 );
     intree->SetBranchAddress("nlpho",          &nlpho        			                 );
     intree->SetBranchAddress("scale1fb",       &scale1fb        			                 );
     intree->SetBranchAddress("lhew",           &lhew);
     intree->SetBranchAddress("btag0Wgt",       &btag0Wgt        			                 );
     intree->SetBranchAddress("btag0WgtUp",     &btag0WgtUp      			                 );
     intree->SetBranchAddress("btag0WgtDown",   &btag0WgtDown    			                 );
     intree->SetBranchAddress("btag1Wgt",       &btag1Wgt        			                 );
     intree->SetBranchAddress("btag1WgtUp",     &btag1WgtUp      			                 );
     intree->SetBranchAddress("btag1WgtDown",   &btag1WgtDown    			                 );
     intree->SetBranchAddress("btag2Wgt",       &btag2Wgt        			                 );
     intree->SetBranchAddress("btag2WgtUp",     &btag2WgtUp      			                 );
     intree->SetBranchAddress("btag2WgtDown",   &btag2WgtDown    			                 );
     intree->SetBranchAddress("ttWgt",          &ttWgt           			                 );
     intree->SetBranchAddress("pfmetraw",       &pfmetraw        			                 );
     intree->SetBranchAddress("pfmetphiraw",    &pfmetphiraw     			                 );
     intree->SetBranchAddress("dphijetmetraw",  &dphijetmetraw   			                 );
     intree->SetBranchAddress("pfmet",          &pfmet           			                 );
     intree->SetBranchAddress("pfmetphi",       &pfmetphi        			                 );
     intree->SetBranchAddress("dphijetmet",     &dphijetmet      			                 );
     intree->SetBranchAddress("dphijetmet4",    &dphijetmet4      			                 );
     intree->SetBranchAddress("calomet",        &calomet         			                 );
     intree->SetBranchAddress("calometphi",     &calometphi                                         );
     intree->SetBranchAddress("recoil",         &recoil                                         );

     // lepton variables
     //  outTree->Branch("lepId", &lepId, "lepId/I");
     intree->SetBranchAddress("electron",    &electron_arr,    &electron_br               );
     intree->SetBranchAddress("muon",    &muon_arr,    &muon_br		                 );
     intree->SetBranchAddress("tau",    &tau_arr,    &tau_br  		                 );
									                 
     intree->SetBranchAddress("lep1charge", &lep1charge				                 );
     intree->SetBranchAddress("lep2charge", &lep2charge                                             ); 
									                 
									                 
     // generator level info						                 
     intree->SetBranchAddress("genVPt",          &genVPt                                            );
     intree->SetBranchAddress("qcdScaleWgtUp",   &qcdScaleWgtUp			                 );
     intree->SetBranchAddress("qcdScaleWgtDown", &qcdScaleWgtDown			                 );
     intree->SetBranchAddress("genId11", &genId11					                 );
     intree->SetBranchAddress("genId12", &genId12					                 );
     intree->SetBranchAddress("genId13", &genId13					                 );
     intree->SetBranchAddress("genId21", &genId21					                 );
     intree->SetBranchAddress("genId22", &genId22					                 );
     intree->SetBranchAddress("genId23", &genId23					                 );
     intree->SetBranchAddress("genpar11", &genpar11		                 );
     intree->SetBranchAddress("genpar12", &genpar12		                 );
     intree->SetBranchAddress("genpar13", &genpar13		                 );
     intree->SetBranchAddress("genpar21", &genpar21		                 );
     intree->SetBranchAddress("genpar22", &genpar22		                 );
     intree->SetBranchAddress("genpar23", &genpar23                               );  
     intree->SetBranchAddress("ishf",   &ishf                                                       );
     intree->SetBranchAddress("nlepgenfortt",   &nlepgenfortt                                                       );
     intree->SetBranchAddress("numb",   &numb                                                       );
    
     // resolved variables
     intree->SetBranchAddress("res_wmva",     &res_wmva                                             );
     intree->SetBranchAddress("res_topmva",   &res_topmva						 );
    
     intree->SetBranchAddress("res_wmva2",     &res_wmva2						 );
     intree->SetBranchAddress("res_topmva2",   &res_topmva2					 );

     intree->SetBranchAddress("res_jet1csv",  &res_jet1csv						 );
     intree->SetBranchAddress("res_jet2csv",  &res_jet2csv						 );
     intree->SetBranchAddress("res_jet3csv",  &res_jet3csv						 );
     intree->SetBranchAddress("res_jet4csv",  &res_jet4csv						 );
     intree->SetBranchAddress("res_jet5csv",  &res_jet5csv						 );
     intree->SetBranchAddress("res_jet6csv",  &res_jet6csv						 );
     intree->SetBranchAddress("res_jet1qgid", &res_jet1qgid					 );
     intree->SetBranchAddress("res_jet2qgid", &res_jet2qgid					 );
     intree->SetBranchAddress("res_jet3qgid", &res_jet3qgid					 );
     intree->SetBranchAddress("res_jet4qgid", &res_jet4qgid					 );
     intree->SetBranchAddress("res_jet5qgid", &res_jet5qgid					 );
     intree->SetBranchAddress("res_jet6qgid", &res_jet6qgid					 );
     intree->SetBranchAddress("res_jet1pull",  &res_jet1pull				 );
     intree->SetBranchAddress("res_jet2pull",  &res_jet2pull				 );
     intree->SetBranchAddress("res_jet3pull",  &res_jet3pull				 );
     intree->SetBranchAddress("res_jet1",  &res_jet1				 );
     intree->SetBranchAddress("res_jet2",  &res_jet2				 );
     intree->SetBranchAddress("res_jet3",  &res_jet3				 );
     intree->SetBranchAddress("res_jet4",  &res_jet4				 );
     intree->SetBranchAddress("res_jet5",  &res_jet5				 );
     intree->SetBranchAddress("res_jet6",  &res_jet6				 );
     
     intree->SetBranchAddress("res_prob",     &res_prob						 );
     intree->SetBranchAddress("res_chisq",    &res_chisq                                          );
     intree->SetBranchAddress("res_cost",     &res_cost						 );
     intree->SetBranchAddress("res_fitmass",  &res_fitmass					 );
     intree->SetBranchAddress("res_fitmassW", &res_fitmassW					 );
     intree->SetBranchAddress("res_bdt_dphij1b", &res_bdt_dphij1b_out				 );
     intree->SetBranchAddress("res_bdt_dphij2b", &res_bdt_dphij2b_out				 );
     intree->SetBranchAddress("res_bdt_drj1b", &res_bdt_drj1b_out				 );
     intree->SetBranchAddress("res_bdt_drj2b", &res_bdt_drj2b_out				 );
     
     intree->SetBranchAddress("res_prob2",     &res_prob2					 );
     intree->SetBranchAddress("res_chisq2",    &res_chisq2					 );
     intree->SetBranchAddress("res_cost2",     &res_cost2					 );
     intree->SetBranchAddress("res_fitmass2",  &res_fitmass2					 );
     intree->SetBranchAddress("res_fitmassW2", &res_fitmassW2					 );
     intree->SetBranchAddress("res_bdt_dphij4b", &res_bdt_dphij4b_out				 );
     intree->SetBranchAddress("res_bdt_dphij5b", &res_bdt_dphij5b_out				 );
     intree->SetBranchAddress("res_bdt_drj4b", &res_bdt_drj4b_out				 );
     intree->SetBranchAddress("res_bdt_drj5b", &res_bdt_drj5b_out				 );
    
     //trigger										 
     intree->SetBranchAddress("passtrigger",       &passtrigger					 );
											 
     //cutflow										 
     intree->SetBranchAddress("cutflow",    &cutflow                                                );
     //END: CHANGE ACCORDING TO runALLHadronicPreselection.cc



    const string cmssw_base = getenv("CMSSW_BASE");
    const string puWeightFilename = cmssw_base + ("/src/DMSAna/Utils/data/PUWeights_Summer16.root");
    TFile puWeightFile(puWeightFilename.c_str());

    TH1D *hPUWeights = (TH1D*)puWeightFile.Get("PUWeights");
    hPUWeights->SetDirectory(0);
    puWeightFile.Close();




	
    cout<<"intree->GetEntries() "<<intree->GetEntries()<<endl;
    for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++) {

    //    for(unsigned int ientry=0; ientry<10; ientry++) {
      if(verbose)
	if(ientry%1000 == 0)
	  cout<<"******ientry "<<ientry<<endl;

      //      if(ientry>10000) continue;

      intree->GetEntry(ientry);



      //decide the btag SF to be used
      if(isSemiLepSelection){
	btagWgt = btag2Wgt;
	btagWgtp = btag2WgtUp;
	btagWgtm = btag2WgtDown;
      }
      else{
	btagWgt = btag1Wgt;
	btagWgtp = btag1WgtUp;
	btagWgtm = btag1WgtDown;
      }


     if(electron_arr){
       electron_arr->Clear(); 
       electron_br->GetEntry(ientry);
     }
     else{
       cout<<"*******WARNING****** electron branch not set"<<endl;
     }
     //     cout<<"AAA muon_arr->GetEntriesFast() "<<muon_arr->GetEntriesFast()<<endl;
     if(muon_arr){
       muon_arr->Clear(); 
       //cout<<"BBB muon_arr->GetEntriesFast() "<<muon_arr->GetEntriesFast()<<endl;
       muon_br->GetEntry(ientry);
       //cout<<"CCC muon_arr->GetEntriesFast() "<<muon_arr->GetEntriesFast()<<endl;
     }
     else{
       cout<<"*******WARNING****** muon branch not set"<<endl;
     }


     ele.erase(ele.begin(),ele.begin()+ele.size());
     muo.erase(muo.begin(),muo.begin()+muo.size());
     ele.clear();
     muo.clear();

     if(electron_arr){ 
       //       for(unsigned int ie=0; ie<electron_arr->GetEntriesFast(); ie++) { //FOUND BUG WITH GetEntriesFast(). in runAllHadronic.cc I have a check nele != vEle.sizee() exit(1)
       for(unsigned int ie=0; ie<nele; ie++) {
	 TLorentzVector *electron_tmp = (TLorentzVector*)electron_arr->At(ie);
	 ele.push_back(*electron_tmp);
	 if(verbose)
	   cout<<"electron Pt "<<electron_tmp->Pt()<<endl;
       }      
     }//if(electron_arr)
     std::sort (ele.begin(), ele.begin()+ele.size(), lorentzvecptorder);

     if(nele != ele.size()){
       cout<<"nele "<<nele<<" ele.size() "<<ele.size()<<" are different..gonna cause problems"<<endl;
       cout<<"runNum, evtNum "<<runNum<<" "<<evtNum<<endl;
       for (std::vector<TLorentzVector>::iterator it = ele.begin() ; it < ele.end(); ++it){
	 cout<<"Pt, eta "<<it->Pt()<<" "<<it->Eta()<<endl;
       }
       exit(1);
     }

     if(ele.size()>=2){
       if(ele[0].Pt()<ele[1].Pt()){ 
	 cerr<<"AAAAAAA Wrong ordering "<<ele.size()<<" ele 0, 1 "<<ele[0].Pt()<<" "<<ele[1].Pt()<<endl;
	 exit(1);
       }
     }

     if(muon_arr){
       //       for(unsigned int im=0; im<muon_arr->GetEntriesFast(); im++) { //FOUND BUG WITH GetEntriesFast(). nmuo != vMuo.size() exit(1): for event with nmuo=1 returned two muons, where the second muon ws the second muon of the previous stored event

       for(unsigned int im=0; im<nmuo; im++) {
	 TLorentzVector *muon_tmp = (TLorentzVector*)muon_arr->At(im);
	 //	 cout<<"muon Pt "<<muon_tmp->Pt()<<endl;
	 muo.push_back(*muon_tmp);
	 if(verbose)
	   cout<<"muon Pt "<<muon_tmp->Pt()<<endl;
       }
     }//if(muon_arr)
     std::sort (muo.begin(), muo.begin()+muo.size(), lorentzvecptorder);

     if(nmuo != muo.size()){
       cout<<"nmuo "<<nmuo<<" muo.size() "<<muo.size()<<" are different..gonna cause problems"<<endl;
       cout<<"runNum, evtNum "<<runNum<<" "<<evtNum<<endl;
       for (std::vector<TLorentzVector>::iterator it = muo.begin() ; it < muo.end(); ++it){
	 cout<<"Pt, eta "<<it->Pt()<<" "<<it->Eta()<<endl;
       }
       exit(1);
     }

     if(muo.size()>=2){
       if(muo[0].Pt()<muo[1].Pt()){ 
	 cerr<<"AAAAAAA Wrong ordering "<<muo.size()<<" muo 0, 1 "<<muo[0].Pt()<<" "<<muo[1].Pt()<<endl;
	 exit(1);
       }
     }


      //END TO BE CHANGED ACCORDING TO runALLHadronicPreselection.cc (for TClonesArray only)
      

     //TMVA top tagger: assuming that res_jet1, 2, 3 are the triplet associated to highest MVA score
     //assuming the genpar1, genpar2, genpar3 are from the hadronic top
     //count the number of match -> if three found hadronic top
      float isHadrTopHighestMVA = -999;
      //      cout<<"infilenames[ifile] "<<infilenames[ifile]<<endl;
      //cout<<"infilenames[ifile].find(TT "<<infilenames[ifile].find("TT_powheg")<<endl;


      if(infilenames[ifile].find("TT_") != string::npos && (njets>2 && genpar11->Pt() >0.001 && genpar12->Pt()  > 0.001 && genpar13->Pt()>0.001))
	 isHadrTopHighestMVA = isHadronicTop(res_jet1, res_jet2, res_jet3, genpar11, genpar12, genpar13, genId11, genId12, genId13, verbose);
      else
	isHadrTopHighestMVA = -1.;


     //     cout<<"isHadrTopHighestMVA "<<isHadrTopHighestMVA<<endl;
     //


     //MT 2 Mar: commented

     // TLorentzVector vLep, vLep2; //highest, second highest energy lep
     // vLep.SetPtEtaPhiM(0.,0,0,0);
     // vLep2.SetPtEtaPhiM(0.,0,0,0);
     // //     TLorentzVector vEle[nele], vMuo[nmuo];
     // TLorentzVector vEle[10], vMuo[10];

     // for(unsigned int ie=0; ie<nele; ie++){
     //   //cout<<"ele "<<ie<<" Pt, Eta "<<ele[ie].Pt()<<", "<<ele[ie].Eta()<<endl;
     //   TLorentzVector* vEle_tmp = (TLorentzVector*)electron_arr->At(ie);
     //   vEle[ie] = *vEle_tmp;
     //   if(vEle[ie].Pt()>vLep.Pt() || vLep.Pt()==0){
     // 	 vLep2 = vLep;
     // 	   vLep = vEle[ie];	   
     //   }
     //   else if (vEle[ie].Pt()>vLep2.Pt() || vLep2.Pt()==0)	   
     // 	 vLep2 = vEle[ie];
     // }
     // for(unsigned int im=0; im<nmuo; im++){
     //   TLorentzVector* vMuo_tmp = (TLorentzVector*)muon_arr->At(im);
     //   vMuo[im] = *vMuo_tmp;

     //   if(vMuo[im].Pt()>vLep.Pt() || vLep.Pt()==0){
     // 	 vLep2 = vLep;
     // 	 vLep = vMuo[im];
     //   }
     //   else if (vMuo[im].Pt()>vLep2.Pt() || vLep2.Pt()==0)
     // 	 vLep2 = vMuo[im]; 
     // }





     // if(nele==1 || nmuo==1){//should be nlep but be careful about ntau
     //   pfmt = computeMt(vLep.Pt(), vLep.Phi(), pfmet, pfmetphi);
     // }
     // else{
     //   pfmt = -999.;
     // }
     

     //

     if(res_jet1 != 0 && res_jet2 !=0 && res_jet3 !=0 && njets>=3){
       toppt = calctoppt(res_jet1,res_jet2,res_jet3);
       topmass = calcmass(res_jet1,res_jet2,res_jet3); //Calctopmass(njets, res_jet1,res_jet2,res_jet3,res_jet4,res_jet5,res_jet6);
     }
     else{
       toppt = -999.;
       topmass = -999.;
     }


     //MT 2 Mar commented

     // if((nele>=1 || nmuo>=1) && njets>=3){//njets>=3 since I need this variable for the top tagger which requires >=3jets
     //   dR_blep = vLep.DeltaR(*res_jet3);
     // }
     // else
     //   dR_blep = -999.;




     if(verbose)
       cout<<"Start selection"<<endl;
      

     if(!passtrigger && isData) continue;

      //selection

     //     cout<<"metfilter "<<metfilter<<" passtrigger "<<passtrigger<<" njets "<<njets<<" nbjets "<<nbjets<<" pfmet "<<pfmet<<" nmuo "<<nmuo<<" pfmt "<<pfmt<<" isHadrTopHighestMVA "<<isHadrTopHighestMVA<<endl;


     if(metfilter & 2)    continue;  // beam halo filter
     if(metfilter & 1)    continue;  // HBHE noise filter
     if(metfilter & 1024) continue;  // HBHEiso noise filter
     if(metfilter & 8)    continue;  // ECAL TP filter
     if(metfilter & 32)   continue;  // ee badSC noise filter
     if(metfilter & 8192) continue;  // badMuon
     if(metfilter & 4096) continue;  // badCharged hadron


     //1lep region
     if(isSemiLepSelection){
       if(njets<4) continue;
       if(pfmet<40) continue;
       if(nbjetsM < 2) continue;
       if(nmuo != 1) continue; 
       if(nlep_oth != 1) continue; 
       if(res_jet1csv<0) continue;
       if(res_jet2csv<0) continue;
     }
     else{
       if(njets<3) continue;
       if(pfmet<40) continue;
       if(nbjetsM < 1) continue;
       if(infilenames[ifile].find("2Tmu") != string::npos){
	 if(nmuo != 2) continue; 
       }
       else if(infilenames[ifile].find("2Tele") != string::npos){
	 if(nele != 2) continue; 
       }
       else if(infilenames[ifile].find("1Tmu1Tele") != string::npos){
	 if(nele != 1 || nmuo != 1) continue; 
       }
       if(lep1charge*lep2charge == 1) continue;
       if(nlep_oth != 2) continue; 
       if(res_jet1csv<0) continue;
       if(res_jet2csv<0) continue;

     }
     //     cout<<"+++++++++infilenames["<<ifile<<"] "<<infilenames[ifile]<<endl;

     //2lep region
     // if(njets  < 3)  continue;
     // if(infilenames[ifile].find("2tightmuon") != string::npos){
     //   //cout<<"qui"<<endl;
     //   if(nmuo != 2) continue; 
     // }
     // if(infilenames[ifile].find("2tightele") != string::npos){ 
     //   if(nele != 2) continue; 
     //   if (fabs((vEle[0].Eta())>2.1 && fabs(vEle[1].Eta())>2.1)) continue;
     // }
     // else if(infilenames[ifile].find("tighteletightmuon") != string::npos){
     //   if(nmuo != 1) continue;
     //   if(nele != 1) continue;
     //   if (fabs(vEle[0].Eta())>2.1) continue;
     // }
     // if(nlep != 2) continue;
     // if(nbjets < 1) continue;
     // if(lep1charge*lep2charge == 1) continue;
     // if(res_jet1csv<0) continue;
     // if(res_jet2csv<0) continue;
     // if(topmass<0 || topmass>1000) continue;

     
     
     //     if(toppt>100) continue;
     //if(toppt<=100 || toppt>200) continue;
     //     if(toppt<=200) continue;

     //     cout<<"isHadrTopHighestMVA "<<isHadrTopHighestMVA<<" nlepgenfortt "<<nlepgenfortt<<endl;


     if(verbose)
       cout<<"Event passed detector level cuts"<<endl;

     if(isSemiLepSelection && !isData){
       if(semilep==1){
	 if(nlepgenfortt==0){//all hadronic is assumed to be signal
	   if(processforTnP!="Signal") continue;
	 }//if(nlepgenfortt==0)
	 else if(nlepgenfortt==1){ //semilep
	   if(processforTnP=="Signal" && isHadrTopHighestMVA != 1) continue;
	   else if(processforTnP=="Bkg" && isHadrTopHighestMVA != 0) continue;
	 }//if(nlepgenfortt==1)
	 else
	   continue;
       }//if(semilep==1)
       else if(semilep==0){ //tdil
	 //	 if(nlepgenfortt==1)
	 if(nlepgenfortt!=2 && infilenames[ifile].find("TT_powheg") != string::npos)
	   continue;
       }
     }//if(isSemiLepSelection)


     //     if(verbose)
     //       cout<<"Event passed gen level cuts"<<endl;


      //
      
	//     cout<<"Event Passed selection evtNum "<<evtNum<<endl;

       
     double weight = 1.;
     double weight_[10] = {1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};
     //       if(verbose==2)
     //	 cout<<"isData "<<isData<<endl;

     if(!isData) {     
       weight_[0] = LUMI*scale1fb;//*8077.00/10573.75;
       	   
       weight_[1] *= hPUWeights->GetBinContent(hPUWeights->FindBin(npu));

       if(btagWgt!=btagWgt) weight_[2] *= 1;
       else                 weight_[2] *= btagWgt;
	 

       //     cout<<"muo size, ele size "<<muo.size()<<" "<<ele.size()<<endl;
       if(singlelepdata==0){
	 weight_[4] *= METturnon13TeV(pfmet,true,useMuonMETturnon)/METturnon13TeV(pfmet,false,useMuonMETturnon);
       }
       else if((singlelepdata==1 || singlelepdata==3) && nmuo>0){
	 weight_[3] *= getMuonTrkSF(muo[0].Eta())*getTightMuonSF(muo[0].Pt(),muo[0].Eta());
	 if(singlelepdata==1)
	   weight_[4] *= getMuonTrigSF(muo[0].Pt(),muo[0].Eta());
	 else if(singlelepdata==3){//implies two muons
	   weight_[3] *= getMuonTrkSF(muo[1].Eta())*getTightMuonSF(muo[1].Pt(), muo[1].Eta());
	   
	   double ineff1Data=0, ineff2Data=0, ineff1MC=0, ineff2MC=0;
	   
	   ineff1Data = 1. - getMuonTrigDataEff(muo[0].Pt(), muo[0].Eta(), LUMI);
	   ineff1MC   = 1. - getMuonTrigMCEff(muo[0].Pt(), muo[0].Eta());
	   
	   ineff2Data = 1. - getMuonTrigDataEff(muo[1].Pt(), muo[1].Eta(), LUMI);
	   ineff2MC   = 1. - getMuonTrigMCEff(muo[1].Pt(), muo[1].Eta());
	     
	   weight_[4] *= (1. - ineff1Data*ineff2Data)/(1. - ineff1MC*ineff2MC);

	 }

       }//else if((singlelepdata==1 || singlelepdata==3) && nmuo>0)
       else if((singlelepdata==2 || singlelepdata==4) && nele>0){
	 weight_[3] *= getTightEleSF(ele[0].Pt(),ele[0].Eta());
	 if(singlelepdata==2)
	   weight_[4] *= getEleTrigSF(ele[0].Pt(), ele[0].Eta());
	 else if (singlelepdata==4){
	   weight_[3] *= getTightEleSF(ele[1].Pt(), ele[1].Eta());
	   double ineff1Data=0, ineff2Data=0, ineff1MC=0, ineff2MC=0;
	     
	   ineff1Data = 1. - getEleTrigDataEff(ele[0].Pt(), ele[0].Eta());
	   ineff1MC   = 1. - getEleTrigMCEff(ele[0].Pt(), ele[0].Eta());

	   ineff2Data = 1. - getEleTrigDataEff(ele[1].Pt(), ele[1].Eta());
	   ineff2MC   = 1. - getEleTrigMCEff(ele[1].Pt(), ele[1].Eta());
	   
	   weight_[4] *= (1. - ineff1Data*ineff2Data)/(1. - ineff1MC*ineff2MC);
	 }

       }//       else if((singlelepdata==2 || singlelepdata==4) && nele>0){
       else if (singlelepdata==5 && nele>0 && nmuo>0){
	 weight_[3] *= getTightEleSF(ele[0].Pt(), ele[0].Eta());
	 weight_[3] *= getTightMuonSF(muo[0].Pt(), muo[0].Eta());

	 double ineff1Data=0, ineff2Data=0, ineff1MC=0, ineff2MC=0;
	     
	     
	 ineff1Data = 1. - getMuonTrigDataEff(muo[0].Pt(), muo[0].Eta(), LUMI);
	 ineff1MC   = 1. - getMuonTrigMCEff(muo[0].Pt(), muo[0].Eta());
	 
	 ineff2Data = 1. - getEleTrigDataEff(ele[0].Pt(), muo[0].Eta());
	 ineff2MC   = 1. - getEleTrigMCEff(ele[0].Pt(), ele[0].Eta());
	 
	 weight_[4] *= (1. - ineff1Data*ineff2Data)/(1. - ineff1MC*ineff2MC);
       }//if (singlelepdata==5 && nele>0 && nmuo>0)
       
 
       if((infilenames[ifile].find("WJetsToLNu") != string::npos)){
	   
	 if(genVPt>=0){	   
	   weight_[6] *= getQCDW(genVPt);
	   weight_[6] *= getEWKW(genVPt);
	 }
       }
       if((infilenames[ifile].find("DYJetsToLL") != string::npos) || (infilenames[ifile].find("ZJetsToNuNu") != string::npos)){

	 if(genVPt>=0){	   
	   weight_[6] *= getQCDZ(genVPt);
	   weight_[6] *= getEWKZ(genVPt);
	 }
	   
       }//if(infilenames[ifile].find("Zjets") != string::npos)

       if(infilenames[ifile].find("TT_powheg") != string::npos)
	 weight_[7] *= ttWgt;
	 
	 


     }//if(!isData)

       
     weight *= weight_[0]*weight_[1]*weight_[2]*weight_[3]*weight_[4]*weight_[5]*weight_[6]*weight_[7]*weight_[8]*weight_[9];
     if(!isData) weight *= NORMTODATA;     
     evtWeight = weight;
     


     //if(verbose){			    
     //       cout<<"rN, eN, weight; met, genVPt "<<runNum<<" "<<evtNum<<" "<<weight<<"; "<<pfmet<<" "<<genVPt<<endl;
     //cout<<"Single weights: LUMI*scale1fb, PU, btgag, lep ID, Trig turnon, others, kfact, topPt, RTT Eff SF, FR2  "<<weight_[0]<<", "<<weight_[1]<<", "<<weight_[2]<<", "<<weight_[3]<<", "<<weight_[4]<<", "<<weight_[5]<<", "<<weight_[6]<<", "<<weight_[7]<<", "<<weight_[8]<<", "<<weight_[9]<<endl;


       //}


     outTree->Fill();
      


     nentries++;
     neventsv +=weight;
    }

    infile->Close();
    delete infile;
    infile=0;
    intree=0;
  }

  outFile->Write();
  outFile->Close();

  cout << endl;
  cout << " <> " << outfilename << " created!" << endl;
  cout << endl;

  cout<<" nentries "<<nentries<<endl;
  cout << " neventsv "<< setprecision(2) << fixed << setw(10) << neventsv<<endl;
  yieldfile << " neventsv "<< setprecision(2) << fixed << setw(10) << neventsv<<endl;

  yieldfile.close();
  gBenchmark->Show("ReformatRunAllHadronicSelectionOutForTnP");
}

