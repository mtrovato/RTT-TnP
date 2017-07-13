#!/bin/bash

#run ./effRTTFit -999 SingleMu_2012-22Jan2013_allhadrbits_TnP.root muon_selection_data Summer12_TTJets_SemiLeptMGDecays_allhadrbits_Signal_TnP.root Summer12_TTJets_SemiLeptMGDecays_allhadrbits_Bkg1_TnP.root Summer12_TTJets_SemiLeptMGDecays_allhadrbits_Bkg2_TnP.root
#and make sure that Data and MC have more or less the same yields. If not may have a problem estimating parameters



#Signal.root and Bkg1.root Bkg2.root is where I get the templates from (make it a list)

#Signal.root = tt1l w/ isHadronictop=1; Bkg1.root=w/ isHadronictop=0; Bkg2.root all the non tt1l bg. Write a function which starts from runAllHadronicSelection output and reformat the output as descibed above (i.e: simpler version of mergeAllHadronic.root)

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`



#evtWeight in the trees carries all weights (NORMTODATA included)
################
#general settings
################
fitLo=0
fitHi=1000
maketemplates=1
TOPMVACUT="0.0" 


#1muo Settings
BIN_SIZE_PASS=20 
BIN_SIZE_FAIL=20
SIGNALFAKERATE12MASK=273 #1=0x001 w/ signal, 16=0x010 w/ bkg1, 256=0x100 w/ bkg2 all 3 = 273
FR2SF="0.990"
FR2ERR="0.004"
MODEL=2
suffix=RTTTempl_1Tmu_2BJ_4J_MET40_MODEL2_Mar32017
dirout=RTT_TnP_1Tmu_2BJ_4J_MET40_Mar32017

#2l Settings
# BIN_SIZE_PASS=1000 
# BIN_SIZE_FAIL=1000
# SIGNALFAKERATE12MASK=256
# FR2SF="1.00"
# FR2ERR="1.00"
# suffix=RTTTempl_2TL_1BJ_3J_MET40_Mar32017
# dirout=RTT_TnP_2TL_1BJ_3J_MET40_Mar32017
###################

rm -rf outputs/${dirout}_old
mv outputs/${dirout} outputs/${dirout}_old

dirin=inputs/${suffix}/


./effRTTFit ${dirin}/Data_TnP${suffix}.root ${dirout} ${dirin}/Signal_TnP${suffix}.root ${dirin}/Bkg1_TnP${suffix}.root ${dirin}/Bkg2_TnP${suffix}.root ${TOPMVACUT} ${maketemplates} ${fitLo} ${fitHi} ${BIN_SIZE_PASS} ${BIN_SIZE_FAIL} ${LUMI} ${SIGNALFAKERATE12MASK} ${FR2SF} ${FR2ERR} ${MODEL}  #>& ${dirout}.log



cp *.cc ${dirout}
cp *.hh ${dirout}
cp run.sh ${dirout}
mv ${dirout}.log ${dirout}
mv ${dirout} outputs
echo "done"

