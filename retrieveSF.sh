#!/bin/bash
####Retrieve SF, stat, syst uncertainties from TnP directories (outputs of runrun.sh). Compute the total uncertainty, by considering stat uncertainty when the given syst varied SF is not significantly different than nominal

abs()
{
#abs
#https://stackoverflow.com/questions/8654051/how-to-compare-two-floating-point-numbers-in-bash
    if (( $(echo $1'<'0 | bc -l) ))
    then
	value="$(echo "${Sys_Herwig}*-1" | bc)"
    else
	value=$1
    fi

    echo $value

}

CompareCentralAndUnc()
{
    cent=$1
    unc=$2
    if  (( $(echo $cent'<'$unc | bc -l) ))
then
   value=$unc
else
   value=$cent
fi

    echo $value

}


################
TOPMVACUT="0.4"
##########

dirin="outputs/RTT04/"
ARRAYSYST=("" "Herwig_" "isrup_" "isrdown_" "fsrup_" "fsrdown_" "hdampUP_" "hdampDOWN_" "dR015_" "dR045_" "qcdscaleup_" "qcdscaledown_" "GaussSmearTempl_")

echo "----------------------------"
echo "dirin=$dirin"

#STEP1
i=0
while [ $i -lt ${#ARRAYSYST[@]} ]; do
    subdir=$dirin$"RTT_TnP_1Tmu_2BJ_4J_MET40_"${ARRAYSYST[$i]}"Jul122017/"
#    echo "file=$file"

    #mean
    grep "SF (Eff)" $subdir/"fitres.txt" | awk '{print $3}' > ${subdir}SF_Eff.txt
    grep "SF (FR1)" $subdir/"fitres.txt" | awk '{print $3}' > ${subdir}SF_FR1.txt
    grep "SF (FR2)" $subdir/"fitres.txt" | awk '{print $3}' > ${subdir}SF_FR2.txt
    #uncertainty (assuming symmetric uncertainties)
    grep "SF (Eff)" $subdir/"fitres.txt" | awk '{print $5}' > ${subdir}SFUnc_Eff.txt
    grep "SF (FR1)" $subdir/"fitres.txt" | awk '{print $5}' > ${subdir}SFUnc_FR1.txt
    grep "SF (FR2)" $subdir/"fitres.txt" | awk '{print $5}' > ${subdir}SFUnc_FR2.txt
    
    let i=i+1
done

#STEP2

option=$1 #Eff, FR1, FR2

i=0
declare -a SF_Eff
declare -a SFUnc_Eff
while [ $i -lt ${#ARRAYSYST[@]} ]; do
    subdir=$dirin$"RTT_TnP_1Tmu_2BJ_4J_MET40_"${ARRAYSYST[$i]}"Jul122017/"
    SF[$i]=`cat ${subdir}/SF_${option}.txt`
    SFUnc[$i]=`cat ${subdir}/SFUnc_${option}.txt`

 #   echo "SF ${ARRAYSYST[$i]} = ${SF[$i]}"
 #   echo "SFUnc ${ARRAYSYST[$i]} = ${SFUnc[$i]}"

    let i=i+1
done

#compute deltas


#nominal
Nom=${SF[0]}
dNom=${SFUnc[0]}

#Herwig
Herwig=${SF[1]}
dHerwig=${SFUnc[1]}
Sys_Herwig="$(echo "${Herwig} - ${Nom}" | bc)"
SysUnc_Herwig="$(echo "sqrt(${dHerwig}*${dHerwig}+${dNom}*${dNom})" | bc)"

#ISR
isrup=${SF[2]}
disrup=${SFUnc[2]}
isrdown=${SF[3]}
disrdown=${SFUnc[3]}
Sys_isr="$(echo "0.5*${isrup} - 0.5*${isrdown}" | bc)"
SysUnc_isr="$(echo "0.5*sqrt(${disrup}*${disrup}+${disrdown}*${disrdown})" | bc)"

#FSR
fsrup=${SF[4]}
dfsrup=${SFUnc[4]}
fsrdown=${SF[5]}
dfsrdown=${SFUnc[5]}
Sys_fsr="$(echo "0.5*${fsrup} - 0.5*${fsrdown}" | bc)"
SysUnc_fsr="$(echo "0.5*sqrt(${dfsrup}*${dfsrup}+${dfsrdown}*${dfsrdown})" | bc)"

#HDAMP
hdampup=${SF[6]}
dhdampup=${SFUnc[6]}
hdampdown=${SF[7]}
dhdampdown=${SFUnc[7]}
Sys_hdamp="$(echo "0.5*${hdampup} - 0.5*${hdampdown}" | bc)"
SysUnc_hdamp="$(echo "0.5*sqrt(${dhdampup}*${dhdampup}+${dhdampdown}*${dhdampdown})" | bc)"

#dR
dRup=${SF[8]}
ddRup=${SFUnc[8]}
dRdown=${SF[9]}
ddRdown=${SFUnc[9]}
Sys_dR="$(echo "0.5*${dRup} - 0.5*${dRdown}" | bc)"
SysUnc_dR="$(echo "0.5*sqrt(${ddRup}*${ddRup}+${ddRdown}*${ddRdown})" | bc)"

#qcdscale
qcdscaleup=${SF[8]}
dqcdscaleup=${SFUnc[8]}
qcdscaledown=${SF[9]}
dqcdscaledown=${SFUnc[9]}
Sys_qcdscale="$(echo "0.5*${qcdscaleup} - 0.5*${qcdscaledown}" | bc)"
SysUnc_qcdscale="$(echo "0.5*sqrt(${dqcdscaleup}*${dqcdscaleup}+${dqcdscaledown}*${dqcdscaledown})" | bc)"

#jecjer
jecjer=${SF[10]}
djecjer=${SFUnc[10]}
Sys_jecjer="$(echo "${jecjer} - ${Nom}" | bc)"
SysUnc_jecjer="$(echo "sqrt(${djecjer}*${djecjer}+${dNom}*${dNom})" | bc)"


DSF[0]=${DSF_TMP[0]}
#absolute value
Sys_Herwig=$(abs $Sys_Herwig)
SYS_Herwig=$(CompareCentralAndUnc $Sys_Herwig $SysUnc_Herwig)
Sys_isr=$(abs $Sys_isr)
SYS_isr=$(CompareCentralAndUnc $Sys_isr $SysUnc_isr)
Sys_fsr=$(abs $Sys_fsr)
SYS_fsr=$(CompareCentralAndUnc $Sys_fsr $SysUnc_fsr)
Sys_hdamp=$(abs $Sys_hdamp)
SYS_hdamp=$(CompareCentralAndUnc $Sys_hdamp $SysUnc_hdamp)
Sys_dR=$(abs $Sys_dR)
SYS_dR=$(CompareCentralAndUnc $Sys_dR $SysUnc_dR)
Sys_qcdscale=$(abs $Sys_qcdscale)
SYS_qcdscale=$(CompareCentralAndUnc $Sys_qcdscale $SysUnc_qcdscale)
Sys_jecjer=$(abs $Sys_jecjer)
SYS_jecjer=$(CompareCentralAndUnc $Sys_jecjer $SysUnc_jecjer)


#no herwig
TOT_UNC="$(echo "sqrt(${dNom}*${dNom}+${SYS_jecjer}*${SYS_jecjer}+${SYS_isr}*${SYS_isr}+${SYS_fsr}*${SYS_fsr}+${SYS_hdamp}*${SYS_hdamp}+${SYS_dR}*${SYS_dR}+${SYS_qcdscale}*${SYS_qcdscale})" | bc)"





echo "$Nom"
echo "$dNom $SYS_jecjer $SYS_Herwig $SYS_isr $SYS_fsr $SYS_hdamp $SYS_dR $SYS_qcdscale    $TOT_UNC"








