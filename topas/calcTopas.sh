#!/bin/bash
# heterogeneity

export TOPAS_G4_DATA_DIR=~/G4Data

directoriesString=$(find $PWD/$1 -maxdepth 3 -type d | sort -n) #get names of subdirectories
directoriesUnsorted=($directoriesString) #save as array
counter=0
for f in "${directoriesUnsorted[@]}"
do
	res=($f/*.txt) #check if folder contains txt file to run	
	if [[ "$res" != *"*"* ]]; then
		directories[counter]=$f
		let counter=counter+1
	fi
done

unset directoriesString
unset directoriesUnsorted
echo "Submitted calculation of ${#directories[@]} folders."
counter=1

for directory in "${directories[@]}"
do	
	cd $directory
	seconds1=$(date '+%s')
	fileTest=(*run5.txt)
	if [ -n "$fileTest" ]; then
		subFiles=(matRad_plan_field*.txt)
		counterSub=1
	fi

	for s in $directory/matRad_plan_*.txt
	do
		# Start second timer if multiple runs are present (long simulation with substeps)
		if [ -n "$fileTest" ]; then
			secondsSub1=$(date '+%s')
		fi
		
		# Start TOPAS calculation of specific txt file
		echo "Starting calculation for file $(basename "$s") in folder $directory"
		~/topas/bin/topas ${s%.*}.txt > ${s%.*}.out > ${s%.*}.err

		# check if produced bin file contains values or is empty
		file=${s##*/}
		if [ -s $directory/score_${file%.*}_physicalDose.bin ]; then
        		echo "Calculation successful.\n"
		else
			seconds2=$(date '+%s')
			let difference[counter]=seconds2-seconds1
        		printf "#######################################\nERROR:\n${HOSTNAME:11:7}: No results in file ${s%.*} of folder '${1:0:-1}'! Finished in ${difference[counter]} seconds.\n#######################################" | telegram-send --stdin
			continue 2
		fi

		# output for multiple subfolders/runs to give updates on the length
		if [ -n "$fileTest" ]; then
			secondsSub2=$(date '+%s')
			let differenceSub[counterSub]=secondsSub2-secondsSub1
			sumSub=$(IFS=+; echo "$((${differenceSub[*]}))")
			let estimatedTimeSub=(${#subFiles[@]}-counterSub)
			let estimatedTimeSub=estimatedTimeSub*sumSub/counterSub

			if [ "$estimatedTimeSub" -ne "0" ]; then
				finishDateSub=$(date '+%a, %d.%m., %H:%Mh' --date="@$((secondsSub2+ estimatedTimeSub))")
				printf "${HOSTNAME:11:7}: Simulation of subFile $file ($counterSub/${#subFiles[@]}) in '$(basename $directory)' finished in ${differenceSub[counterSub]} seconds.\nEstimated date finished: $finishDateSub." | telegram-send --silent --stdin
			else
				printf "${HOSTNAME:11:7}: Simulation of subFile $file ($counterSub/${#subFiles[@]}) in '$(basename $directory)' finished in ${differenceSub[counterSub]} seconds." | telegram-send --silent --stdin	
			fi
			let counterSub=counterSub+1
		fi

	done

	seconds2=$(date '+%s')
	let difference[counter]=seconds2-seconds1
	sum=$(IFS=+; echo "$((${difference[*]}))")
	let estimatedTime=(${#directories[@]}-counter)
	let estimatedTime=estimatedTime*sum/counter

	if [ "$estimatedTime" -ne "0" ]; then
		finishDate=$(date '+%a, %d.%m., %H:%Mh' --date="@$((seconds2 + estimatedTime))")
		printf "${HOSTNAME:11:7}: Simulation of folder $counter/${#directories[@]} in '${1:0:-1}' finished in $(date -d@${difference[counter]} -u +%Hh%Mm).\nEstimated date finished: $finishDate." | telegram-send --silent --stdin
	else
		printf "${HOSTNAME:11:7}: Simulation of folder $counter/${#directories[@]} in '${1:0:-1}' finished in $(date -d@${difference[counter]} -u +%Hh%Mm)." | telegram-send --silent --stdin	
	fi
	let counter=counter+1
done

cd $PWD
printf "***************************************\n${HOSTNAME:11:7}: Simulation of folder '${1:0:-1}' completed. Total time: $(date -d@$sum -u +%Hh%Mm).\n***************************************" | telegram-send --stdin
echo "All done!"