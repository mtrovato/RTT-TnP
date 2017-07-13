#!/bin/bash

#run ./effZFit -999 SingleMu_2012-22Jan2013_allhadrbits_TnP.root muon_selection_data Summer12_TTJets_SemiLeptMGDecays_allhadrbits_Signal_TnP.root Summer12_TTJets_SemiLeptMGDecays_allhadrbits_Bkg1_TnP.root Summer12_TTJets_SemiLeptMGDecays_allhadrbits_Bkg2_TnP.root
#and make sure that Data and MC have more or less the same yields. If not may have a problem estimating parameters



#Signal.root and Bkg1.root Bkg2.root is where I get the templates from (make it a list)

#Signal.root = tt1l w/ isHadronictop=1; Bkg1.root=w/ isHadronictop=0; Bkg2.root all the non tt1l bg. Write a function which starts from runAllHadronicSelection output and reformat the output as descibed above (i.e: simpler version of mergeAllHadronic.root)

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`



#evtWeight in the trees carries all weights (NORMTODATA included)

#general settings
fitLo=0
fitHi=1000
maketemplates=1

#1muo Settings
#BIN_SIZE_PASS=20 
#BIN_SIZE_FAIL=20
#SIGNALFAKERATE12MASK=273 #1=0x001 w/ signal, 16=0x010 w/ bkg1, 256=0x100 w/ bkg2 all 3 = 273
#suffix=Jan132016_1muon1looselepnjets4nbjets2met40_topptmore200

#2l Settings
BIN_SIZE_PASS=1000 
BIN_SIZE_FAIL=1000
SIGNALFAKERATE12MASK=256
suffix=RTTTempl_2TL_1BJ_3J_MET40_Mar32017









dirin=inputs/${suffix}/


dir=RTT_TnP_2TL_1BJ_3J_MET40_Mar32017



TOPMVACUT="0.0" 
FR2SF="1.00"
FR2ERR="1.00"

echo "******"; 
echo "MVACUT $TOPMVACUT"; 
echo "FR2SF $FR2SF"; 
echo "FR2ERR $FR2ERR"; 
echo "******";


./effZFit ${dirin}/Data_TnP${suffix}.root ${dir} ${dirin}/Signal_TnP${suffix}.root ${dirin}/Bkg1_TnP${suffix}.root ${dirin}/Bkg2_TnP${suffix}.root ${TOPMVACUT} ${maketemplates} ${fitLo} ${fitHi} ${BIN_SIZE_PASS} ${BIN_SIZE_FAIL} ${LUMI} ${SIGNALFAKERATE12MASK} ${FR2SF} ${FR2ERR}  #>& ${dir}.log



cp *.cc ${dir}
cp *.hh ${dir}
cp run_example.sh ${dir}
mv ${dir}.log ${dir}
mv ${dir} outputs
echo "done"

