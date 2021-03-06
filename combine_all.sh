#!/bin/bash

for A in dum # 3H 3He 2H 1H 
do
    for kin in slow mid fast  # mid slow fast #slow  mid
    do 
	if [ -e skim/${A}_${kin} ]
	then

	    if [ `ls skim/${A}_${kin} | wc -l` -eq "0" ]
	    then
		echo "$A $kin has no runs to combine."
	    else 
		echo "Working on $A $kin ..."
		./combiner skim/${A}_${kin}/*.root
		mv combine_out.root combiner_out/${A}_${kin}.root
	    fi
	fi
    done
done
