#!/bin/bash

script_folder=$0
folder=$1

while read -r line; do
	IFS=',' read -ra x <<< "${line}"

	echo "fastq-dump --split-files --gzip ${x[0]}"
	fastq-dump --split-files --gzip -O ${folder} ${x[0]}
done < ${script_folder%/*}/Metadata_immune29.txt
