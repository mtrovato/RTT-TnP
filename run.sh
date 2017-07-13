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
dateout="Jul122017"
fitLo=0
fitHi=1000
maketemplates=1
TOPMVACUT=$2 #"0.0" 
VERBOSE=1
isSyst=$1 #0
doSingleMu=$3 #1 #1=1mu; !=1=dil
################

##############
if [ $doSingleMu -eq 1 ]
then
#1muo Settings
##############
    BIN_SIZE_PASS=20 
    BIN_SIZE_FAIL=20
    SIGNALFAKERATE12MASK=273 #1=0x001 w/ signal, 16=0x010 w/ bkg1, 256=0x100 w/ bkg2 all 3 = 273
############
    if [ $isSyst -eq 0 ] #nominal
    then
	FR2SF=$4 # "0.990" #from dil TnP
	FR2ERR=$5 # "0.004"
	MODEL=1
	suffix=RTTTempl_1Tmu_2BJ_4J_MET40_Mar32017
	dirout=RTT_TnP_1Tmu_2BJ_4J_MET40_${dateout}
    elif [ $isSyst -eq 1 ]
    then
	FR2SF=$4 # "0.917"
	FR2ERR=$5 # "0.004"
	MODEL=1
	suffix=RTTTempl_1Tmu_2BJ_4J_MET40_Herwig_Mar32017
	dirout=RTT_TnP_1Tmu_2BJ_4J_MET40_Herwig_${dateout}
    elif [ $isSyst -eq 2 ]
    then
	FR2SF=$4 # "0.970"
	FR2ERR=$5 # "0.004"
	MODEL=1
	suffix=RTTTempl_1Tmu_2BJ_4J_MET40_isrup_Mar32017
	dirout=RTT_TnP_1Tmu_2BJ_4J_MET40_isrup_${dateout}
    elif [ $isSyst -eq 3 ]
    then
	FR2SF=$4 # "1.006"
	FR2ERR=$5 # "0.004"
	MODEL=1
	suffix=RTTTempl_1Tmu_2BJ_4J_MET40_isrdown_Mar32017
	dirout=RTT_TnP_1Tmu_2BJ_4J_MET40_isrdown_${dateout}
    elif [ $isSyst -eq 4 ]
    then
	FR2SF=$4 # "1.015"
	FR2ERR=$5 # "0.004"
	MODEL=1
	suffix=RTTTempl_1Tmu_2BJ_4J_MET40_fsrup_Mar32017
	dirout=RTT_TnP_1Tmu_2BJ_4J_MET40_fsrup_${dateout}
    elif [ $isSyst -eq 5 ]
    then
	FR2SF=$4 # "0.983"
	FR2ERR=$5 # "0.004"
	MODEL=1
	suffix=RTTTempl_1Tmu_2BJ_4J_MET40_fsrdown_Mar32017
	dirout=RTT_TnP_1Tmu_2BJ_4J_MET40_fsrdown_${dateout}
    elif [ $isSyst -eq 6 ]
    then
	FR2SF=$4 # "0.990"
	FR2ERR=$5 # "0.004"
	MODEL=1
	suffix=RTTTempl_1Tmu_2BJ_4J_MET40_hdampUP_Mar32017
	dirout=RTT_TnP_1Tmu_2BJ_4J_MET40_hdampUP_${dateout}
    elif [ $isSyst -eq 7 ]
    then
	FR2SF=$4 # "0.997"
	FR2ERR=$5 # "0.004"
	MODEL=1
	suffix=RTTTempl_1Tmu_2BJ_4J_MET40_hdampDOWN_Mar32017
	dirout=RTT_TnP_1Tmu_2BJ_4J_MET40_hdampDOWN_${dateout}
    elif [ $isSyst -eq 8 ]
    then
	#same inputs as nominal
	FR2SF=$4 # "0.990"
	FR2ERR=$5 # "0.004"
	MODEL=1
	suffix=RTTTempl_1Tmu_2BJ_4J_MET40_dR015_Mar32017
	dirout=RTT_TnP_1Tmu_2BJ_4J_MET40_dR015_${dateout}
    elif [ $isSyst -eq 9 ]
    then
	#same inputs as nominal
	FR2SF=$4 # "0.990"
	FR2ERR=$5 # "0.004"
	MODEL=1
	suffix=RTTTempl_1Tmu_2BJ_4J_MET40_dR045_Mar32017
	dirout=RTT_TnP_1Tmu_2BJ_4J_MET40_dR045_${dateout}
    elif [ $isSyst -eq 10 ]
    then
	FR2SF=$4 # "0.991"
	FR2ERR=$5 # "0.004"
	MODEL=1
	suffix=RTTTempl_1Tmu_2BJ_4J_MET40_qcdscaleup_Mar32017
	dirout=RTT_TnP_1Tmu_2BJ_4J_MET40_qcdscaleup_${dateout}
    elif [ $isSyst -eq 11 ]
    then
	FR2SF=$4 # "0.986"
	FR2ERR=$5 # "0.004"
	MODEL=1
	suffix=RTTTempl_1Tmu_2BJ_4J_MET40_qcdscaledown_Mar32017
	dirout=RTT_TnP_1Tmu_2BJ_4J_MET40_qcdscaledown_${dateout}
    elif [ $isSyst -eq 12 ]
    then
	#same inputs as nominal
	FR2SF=$4 # "0.990"
	FR2ERR=$5 # "0.004"
	MODEL=2
	suffix=RTTTempl_1Tmu_2BJ_4J_MET40_Mar32017
	dirout=RTT_TnP_1Tmu_2BJ_4J_MET40_GaussSmearTempl_${dateout}
    fi #$isSyst -eq 0
else
#2l Settings
    BIN_SIZE_PASS=1000 
    BIN_SIZE_FAIL=1000
    SIGNALFAKERATE12MASK=256
    FR2SF="1.00"
    FR2ERR="1.00"
    MODEL=1

    if [ $isSyst -eq 0 ] #nominal
    then
	suffix=RTTTempl_2TL_1BJ_3J_MET40_Mar32017
	dirout=RTT_TnP_2TL_1BJ_3J_MET40_${dateout}
    elif [ $isSyst -eq 1 ]
	then
	suffix=RTTTempl_2TL_1BJ_3J_MET40_Herwig_Mar32017
	dirout=RTT_TnP_2TL_1BJ_3J_MET40_Herwig_${dateout}
    elif [ $isSyst -eq 2 ]
	then
	suffix=RTTTempl_2TL_1BJ_3J_MET40_isrup_Mar32017
	dirout=RTT_TnP_2TL_1BJ_3J_MET40_isrup_${dateout}
    elif [ $isSyst -eq 3 ]
	then
	suffix=RTTTempl_2TL_1BJ_3J_MET40_isrdown_Mar32017
	dirout=RTT_TnP_2TL_1BJ_3J_MET40_isrdown_${dateout}
    elif [ $isSyst -eq 4 ]
	then
	suffix=RTTTempl_2TL_1BJ_3J_MET40_fsrup_Mar32017
	dirout=RTT_TnP_2TL_1BJ_3J_MET40_fsrup_${dateout}
    elif [ $isSyst -eq 5 ]
	then
	suffix=RTTTempl_2TL_1BJ_3J_MET40_fsrdown_Mar32017
	dirout=RTT_TnP_2TL_1BJ_3J_MET40_fsrdown_${dateout}
    elif [ $isSyst -eq 6 ]
	then
	suffix=RTTTempl_2TL_1BJ_3J_MET40_hdampUP_Mar32017
	dirout=RTT_TnP_2TL_1BJ_3J_MET40_hdampUP_${dateout}
    elif [ $isSyst -eq 7 ]
	then
	suffix=RTTTempl_2TL_1BJ_3J_MET40_hdampDOWN_Mar32017
	dirout=RTT_TnP_2TL_1BJ_3J_MET40_hdampDOWN_${dateout}
    elif [ $isSyst -eq 8 ] || [ $isSyst -eq 9 ] 
        then
	echo "USE FR2, ERRFR2 NOMINAL VALUES"
	exit 
    elif [ $isSyst -eq 10 ]
	then
	suffix=RTTTempl_2TL_1BJ_3J_MET40_qcdscaleup_Mar32017
	dirout=RTT_TnP_2TL_1BJ_3J_MET40_qcdscaleup_${dateout}
    elif [ $isSyst -eq 11 ]
	then
	suffix=RTTTempl_2TL_1BJ_3J_MET40_qcdscaledown_Mar32017
	dirout=RTT_TnP_2TL_1BJ_3J_MET40_qcdscaledown_${dateout}
    elif [ $isSyst -eq 12 ] 
    then
	echo "USE FR2, ERRFR2 NOMINAL VALUES"
	exit
    fi 
fi #if [ $doSingleMu -eq 1 ]


###################

rm -rf outputs/${dirout}_old
mv outputs/${dirout} outputs/${dirout}_old

dirin=inputs/${suffix}/


####print settings
echo "*****************************"
echo "dateout=${dateout}"
echo "fitLo=${fitLo} fitHi=${fitHi}"
echo "maketemplates=${maketemplates}"
echo "TOPMVACUT=${TOPMVACUT}"
echo "VERBOSE=${VERBOSE}"
echo "isSyst=${isSyst}"
echo "doSingleMu=${doSingleMu}" 
echo "FR2SF=${FR2SF}"
echo "FR2ERR=${FR2ERR}"
echo "MODEL=${MODEL}"
echo "suffix=${suffix}"
echo "dirout=${dirout}"
echo "BIN_SIZE_PASS=${BIN_SIZE_PASS} BIN_SIZE_FAIL=${BIN_SIZE_FAIL}"
echo "SIGNALFAKERATE12MASK=${SIGNALFAKERATE12MASK}"
echo "*****************************"
###

./effRTTFit ${dirin}/Data_TnP${suffix}.root ${dirout} ${dirin}/Signal_TnP${suffix}.root ${dirin}/Bkg1_TnP${suffix}.root ${dirin}/Bkg2_TnP${suffix}.root ${TOPMVACUT} ${maketemplates} ${fitLo} ${fitHi} ${BIN_SIZE_PASS} ${BIN_SIZE_FAIL} ${LUMI} ${SIGNALFAKERATE12MASK} ${FR2SF} ${FR2ERR} ${MODEL} ${VERBOSE}  >& ${dirout}.log



cp *.cc ${dirout}
cp *.hh ${dirout}
cp run.sh ${dirout}
mv ${dirout}.log ${dirout}
mv ${dirout} outputs
echo "done"

