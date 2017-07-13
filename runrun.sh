#!/bin/bash
################
TOPMVACUT="0.4"
##########

#singleMu
ARRAYISSSYS=(0 1 2 3 4 5 6 7 8 9 10 11 12)
doSingleMu=1
FR2SF=(0.996 0.898 0.972 1.016 1.044 0.971 0.989 1.004 0.996 0.996 0.998 0.991 0.996)
FR2ERR=(0.006 0.006 0.006 0.007 0.007 0.006 0.006 0.007 0.006 0.006 0.006 0.006 0.006)

#DIL
#ARRAYISSSYS=(0 1 2 3 4 5 6 7 10 11)
#doSingleMu=0


echo "----------------------------"
echo "doSingleMu=$doSingleMu"

i=0
while [ $i -lt ${#ARRAYISSSYS[@]} ]; do
    echo "isSyst=${ARRAYISSSYS[$i]}"
    echo "FR2SF=${FR2SF[$i]}"
    echo "FR2ERR=${FR2ERR[$i]}"
    echo "----------------------------"
    source run.sh ${ARRAYISSSYS[$i]} $TOPMVACUT $doSingleMu ${FR2SF[$i]} ${FR2ERR[$i]}

    let i=i+1
done






