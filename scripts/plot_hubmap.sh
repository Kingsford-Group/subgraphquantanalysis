#!/bin/bash

script_folder=$0
quant_folder=$1
salmon_bootstrap_file_path=$2

hubmap_folder=${quant_folder}/HumanBodyMap

while read -r line; do
	ID=${line}
	# convert salmon bootstrap output
	python ${salmon_bootstrap_file_path} ${hubmap_folder}/${ID}/ ${hubmap_folder}/${ID}
done < ${script_folder}/../data/Metadata_hubmap.txt

# plotting
python3 ${script_folder}/evaluate_nonidentifiable_degree.py ${hubmap_folder} ${quant_folder}
