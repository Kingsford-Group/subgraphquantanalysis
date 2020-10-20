#!/bin/bash

script_folder=$0
quant_folder=$1

# compute I value
python ${script_folder}/ComputeIValue_mcf10.py ${quant_folder}

# plot example
Rscript scripts/draw_figure_mcf10.R ${quant_folder}
