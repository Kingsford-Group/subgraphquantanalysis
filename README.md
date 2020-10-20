# Overview
This repo contains the analysis code for the uncertainty of expression estimates due to non-identifiability. The algorithm for the bounding the range of uncertain expression estimates is implemented and combined in the graph quantification repo ([https://github.com/Kingsford-Group/subgraphquant](https://github.com/Kingsford-Group/subgraphquant)). Specifically, bounding the range of uncertain expression under graph quantification (assuming reference transcripts are incomplete) is in [flow_graph.py](https://github.com/Kingsford-Group/subgraphquant/blob/master/src/flow_graph.py), and under reference-transcript quantification (assuming reference transcripts are complete) is in [lp_fixed_transcripts.py](https://github.com/Kingsford-Group/subgraphquant/blob/master/src/lp_fixed_transcripts.py). This analysis is based on the output of these codes.

This analysis contains two parts: evaluating the degree of expression estimation uncertainty across multiple tissues using Human Body Map data, evaluating the detected differentially expressed (DE) transcripts in a MCF10 dataset and in a T-cell dataset.

The analysis depends on the following python and R packages:
python packages:
+ [matplotlib](https://matplotlib.org/)
+ [matplotlib-venn](https://pypi.org/project/matplotlib-venn/)
+ [pandas](https://pandas.pydata.org/)
+ [seaborn](https://seaborn.pydata.org/)
R packages:
+ [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html)
+ [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)


# Preparing for all datasets and running quantification
The SRA accession of three datasets can be found in `data/Metadata_*.txt`. The following command will download required reference files and RNA-seq datasets, and quantify the expression using [Salmon](https://salmon.readthedocs.io/en/latest/) and [Graph Salmon](https://github.com/Kingsford-Group/subgraphquant). Salmon generates the expression estimates under the assumption that the reference transcripts are complete; Graph Salmon estimates the abundances of splice graph edges assuming that the reference transcripts are incomplete and that any full paths in splice graphs can express.
```
./script/metarun.sh <output folder>
```

# Evaluating the degree of expression estimation uncertainty in Human Body Map dataset
Using Human Body Map dataset, we evaluate the number of transcripts for which the uncertainty in unidentifiable expression estimates is so large that the ranking of expression between the transcript and its sibling isoforms cannot be determined. The following script will generate a figure that reproduces the figure in manuscript. The figure is named `percentage_flipped_transcripts.pdf` and located in the Human Body Map quantification directory.
```
python3 scripts/evaluate_nonidentifiable_degree.py <path to Human Body Map quantification output> <path to data downloading output>
```

# Evaluating the reliability of differentially expressed (DE) transcripts by uncertainty of expression estimates
We detect differentially expressed (DE) transcripts on a MCF10 and a T-cell dataset using Salmon quantification, tximport, and DESeq2. We then evaluate the reliability of the detected DE transcripts by comparing the ranges of uncertain expression estimates between DE groups. The following command will run DE detection and reproduce the figure for MCF10 dataset. The generated figure will be located in the MCF10 quantification folder under name `IValue_mcf10.pdf`.
```
python scripts/ComputeIValue.py <path to data downloading folder>
Rscript scripts/draw_figure_mcf10.R <path to data downloading folder>
```

The following command will run DE detection and reproduce the figure for T-cell dataset. The generated figure will be located in the T-cell quantification folder under name `IValue_immune29.pdf`.
```
python scripts/ComputeIValue.py <path to data downloading folder>
Rscript scripts/draw_figure_immune29.R <path to data downloading folder>
```
