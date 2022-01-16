#!/usr/bin/bash
# script compares all sequence files in a directory against each other using mash distance
# dependent on having mash distance estimation tool

#When running script this is argument 1 after script name

file_directory="sequence_files" 
file_extension=".fa"
kmer_size="21"
sketch_size="5000"
available_cores=4 # available cores for parallel processing. Modify this depending on how many your computer has.
output_file="output_file_example.tsv"


cd "${file_directory}" # changes directory to where seuqnece files are stored, this is system argument 1
mkdir temp # creates temporary directory
ls *"${file_extension}" > temp/genome_names.txt # lists all sequence files into a text file


readarray -t genome_file_array < temp/genome_names.txt # generates an array out of the sequence files to run one by one
start_time=$(date +%s) # finds current time seconds for start
cores_used=0 # 
    
for a in "${genome_file_array[@]}" # loops over sequence file array
do


    if [[ "${cores_used}" -eq "${available_cores}" ]] # checks how many cores are currently being used, if all cores being used then statement activtaed
    	
        then
        	cores_used=0 # resets cores used

        	wait # waits till all cores have been used before continuing
 
    	fi
    
    
    cores_used=$((cores_used+1))
    
	mash sketch -k "${kmer_size}" -s "${sketch_size}" ${a} > /dev/null 2>&1 # runs mash sketch of sequence file, can be modified if differing sketch size/distance wanted


done

wait
end_time=$(date +%s) # finds current time in seconds
elapsed=$(( end_time - start_time )) # calculates time took in seconds
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')" # runs command to give overall time took


start_time=$(date +%s) # same as before, finds starting time
echo -e "genome_1\tgenome_2\tmash distance\tp_val\tmatching kmers" > "${output_file}" # system argument 2 is results file, recommend it to be a tsv eg. "mash_calc_result.tsv"

for a in "${genome_file_array[@]}" # loops over sequence file array
do
    sed -i '1d' temp/genome_names.txt # 
    readarray -t second_genome_file_array < temp/genome_names.txt



	cores_used=0 # initial cores used value for running mash distance
	
    for b in "${second_genome_file_array[@]}" # loops over the genome file again
	do
    	
        
        if [[ ${a} == ${b} ]] # means that it doesn't compare the sequence file to itself by checking the two sequence files aren't the same
    	
        then
        	
            continue # so continues onto next item in loop if this is the case
    	
        fi
   	 
    	if [[ "${cores_used}" -eq "${available_cores}" ]] # checks how many cores are currently being used, if all cores being used then statement activtaed
    	
        then
        	cores_used=0 # resets cores used
        	wait # waits till all cores have been used before continuing
 
    	fi
   	 

    	cores_used=$((cores_used+1))


    	mash dist "${a}.msh" "${b}.msh" >> "${output_file}" & # runs mash distance against each other using msh files already produced



	done
	wait # makes sure all mash distancing for this loop has been done
done

mkdir mash_outputs # creates a mash outputs directory, this is in the initial directory of the sequence files
mv *.msh mash_outputs # moves the .msh sketch files produced to there, maybe not necessary so can always rm *.msh instead

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"

rm -r temp