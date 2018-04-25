#!/bin/bash

for kin in slow
do
    for A in 3He 
    do
	echo "Working on isotope $A in $kin kinematics ..."

	# Test if a runlist exists
	runlist=infiles/${A}_${kin}_kin.dat

	if [ -e $runlist ]
	then
	    for run in {100000..103000}
	    do
		
		# Test if that run number exists in this runlist
		if grep -Fq "$run" $runlist
		then
		    echo "Working on run $run ..."

		    # Create the skim output file
		    skim_out=/chafs1/work1/tritium/coinc/inclusive_skim/${A}_${kin}/skim_${run}.root
		    # Test if such a skimmed file already exists
		    if [ -e $skim_out ]
		    then
			echo "File has already been skimmed."
		    else
			# Test if at least one replayed file exists
			if [ -e /chafs1/work1/tritium/Rootfiles/tritium_${run}.root ]
			then
			    ./incl_skimmer $run "$A" "$kin" /chafs1/work1/tritium/Rootfiles/tritium_${run}*.root
			else
			    echo "Replayed file is missing. Skipping."
			fi
		    fi
		fi
	    done
	else
	    echo "No runlist file exists. Moving on..."
	fi
    done
done
