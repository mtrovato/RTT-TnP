#! /bin/bash

#INPUT: output of update.sh

######settings###
dir=../../baconbits/baconbits_1Tmu_2BJ_4J_RTT_Feb282017
suffix=RTTTempl_1Tmu_2BJ_4J_MET40_Mar32017
LUMI=35.867
NORMTODATA=0.8805565355; #253655.00/288062.14
#################

####processes####
d0_dib=("Summer16_WWToLNuQQ_powheg_allhadrbits" "Summer16_WWTo2L2Nu_powheg_allhadrbits" "Summer16_WWTo4Q_powheg_allhadrbits" "Summer16_WZTo1L1Nu2Q_amcatnlo_allhadrbits" "Summer16_WZTo1L3Nu_amcatnlo_allhadrbits" "Summer16_WZTo2L2Q_amcatnlo_allhadrbits" "Summer16_WZTo3LNu_powheg_allhadrbits" "Summer16_ZZTo2L2Nu_powheg_allhadrbits" "Summer16_ZZTo2L2Q_amcatnlo_allhadrbits" "Summer16_ZZTo2Q2Nu_amcatnlo_allhadrbits" "Summer16_ZZto4L_amcatnlo_allhadrbits" "Summer16_ZZTo4Q_amcatnlo_allhadrbits")

d0_st=("Summer16_ST_s_channel_4f_leptonDecays_allhadrbits" "Summer16_ST_t_channel_top_4f_inclusiveDecays_TuneCUETP8M2T4_allhadrbits" "Summer16_ST_t_channel_antitop_4f_inclusiveDecays_TuneCUETP8M2T4_allhadrbits" "Summer16_ST_tW_antitop_5f_NoFullyHadronicDecays_allhadrbits" "Summer16_ST_tW_top_5f_NoFullyHadronicDecays_allhadrbits")

d0_zj=("Summer16_ZJetsToNuNu_HT-100To200_allhadrbits" "Summer16_ZJetsToNuNu_HT-200To400_allhadrbits" "Summer16_ZJetsToNuNu_HT-400To600_allhadrbits" "Summer16_ZJetsToNuNu_HT-600To800_allhadrbits" "Summer16_ZJetsToNuNu_HT-800To1200_allhadrbits" "Summer16_ZJetsToNuNu_HT-1200To2500_allhadrbits" "Summer16_ZJetsToNuNu_HT-2500ToInf_allhadrbits")

d0_dyll=("Summer16_DYJetsToLL_M-50_HT-70to100_allhadrbits" "Summer16_DYJetsToLL_M-50_HT-100to200_allhadrbits" "Summer16_DYJetsToLL_M-50_HT-200to400_allhadrbits" "Summer16_DYJetsToLL_M-50_HT-400to600_allhadrbits" "Summer16_DYJetsToLL_M-50_HT-600to800_allhadrbits" "Summer16_DYJetsToLL_M-50_HT-800to1200_allhadrbits" "Summer16_DYJetsToLL_M-50_HT-1200to2500_allhadrbits" "Summer16_DYJetsToLL_M-50_HT-2500toInf_allhadrbits")

d0_wjln=("Summer16_WJetsToLNu_HT-70To100_allhadrbits" "Summer16_WJetsToLNu_HT-100To200_allhadrbits" "Summer16_WJetsToLNu_HT-200To400_allhadrbits" "Summer16_WJetsToLNu_HT-400To600_allhadrbits" "Summer16_WJetsToLNu_HT-600To800_allhadrbits" "Summer16_WJetsToLNu_HT-800To1200_allhadrbits" "Summer16_WJetsToLNu_HT-1200To2500_allhadrbits" "Summer16_WJetsToLNu_HT-2500ToInf_allhadrbits")

ttonly=Summer16_TT_powheg_allhadrbits
d0_tt=(${ttonly} "Summer16_TTGJets_allhadrbits" "Summer16_TTZToLLNuNu_M-10_allhadrbits" "Summer16_TTZToQQ_allhadrbits" "Summer16_TTWJetsToLNu_allhadrbits" "Summer16_TTWJetsToQQ_allhadrbits")


d0_qcd=("Summer16_QCD_HT100to200_allhadrbits" "Summer16_QCD_HT200to300_allhadrbits" "Summer16_QCD_HT300to500_allhadrbits" "Summer16_QCD_HT500to700_allhadrbits" "Summer16_QCD_HT700to1000_allhadrbits" "Summer16_QCD_HT1000to1500_allhadrbits" "Summer16_QCD_HT1500to2000_allhadrbits" "Summer16_QCD_HT2000toInf_allhadrbits")  

processData=SingleMuon_Run2016-ReReco_allhadrbits
processSignal=${ttonly}
processBkg1=${ttonly}
processBkg2=("${d0_dib[@]}" "${d0_st[@]}" "${d0_zj[@]}" "${d0_dyll[@]}" "${d0_wjln[@]}" "${d0_tt[@]}" "${d0_qcd[@]}")



# echo "Data: $processData"
# echo "Signal: $processSignal"
# echo "Bkg1: $processBkg1"
# echo "Bkg2: "
# for (( i=0; i<${#processBkg2[@]}; i++ ))
# do
#     process_tmp=${processBkg2[i]}
#     echo $process_tmp
# done



################


rm listData.txt listSignal.txt listBkg1.txt listBkg2.txt

#Data
ls ../${dir}/${processData}.root > listData.txt;    
root -l -q -b ReformatRunAllHadronicSelectionOutForTnP.C+\(\"listData.txt\",\"Data_TnP${suffix}.root\",\"$LUMI\",\"$NORMTODATA\",\"\",true\); 


#Signal
ls ../${dir}/${processSignal}.root > listSignal.txt;    
root -l -q -b ReformatRunAllHadronicSelectionOutForTnP.C+\(\"listSignal.txt\",\"Signal_TnP${suffix}.root\",\"$LUMI\",\"$NORMTODATA\",\"Signal\",false\); # I save all input events w/ isAllHadronicTop=1

#Bkg1
ls ../${dir}/${processBkg1}.root > listBkg1.txt;    
root -l -q -b ReformatRunAllHadronicSelectionOutForTnP.C+\(\"listBkg1.txt\",\"Bkg1_TnP${suffix}.root\",\"$LUMI\",\"$NORMTODATA\",\"Bkg\",false\); #all input events w/ isAllHadronicTop=0

#Bkg2
for (( i=0; i<${#processBkg2[@]}; i++ ))
do
    process_tmp=${processBkg2[i]}
    ls ../${dir}/${process_tmp}.root >> listBkg2.txt; 
done

root -l -q -b ReformatRunAllHadronicSelectionOutForTnP.C+\(\"listBkg2.txt\",\"Bkg2_TnP${suffix}.root\",\"$LUMI\",\"$NORMTODATA\",\"\",false\);


rm listData.txt listSignal.txt listBkg1.txt listBkg2.txt

dir=../../../ReformatforTnPoutput/${suffix}
mkdir $dir
mv Data_TnP${suffix}.root ${dir}
mv Signal_TnP${suffix}.root ${dir}
mv Bkg1_TnP${suffix}.root ${dir}
mv Bkg2_TnP${suffix}.root ${dir}
cp ReformatRunAllHadronicSelectionOutForTnP.sh ${dir}
cp ReformatRunAllHadronicSelectionOutForTnP.C ${dir}





