# StarterKit
A set of simple R scripts for basic, stand-alone analysis of HuProt (gpr) result files.


## Overview
The following StartKit provides standalone scripts for the generation of basic QC plots and Analytic Results tables from raw gpr files provided by CDI Lab. They constitute a subset of the techniques used in-house and are provided for the benefit of the users who would like to independently (re)-produce the core results in CDI-provided reports.


## Required .tsv files
The inputs to a given StartKit script (of which there are two, namely `QC.R` and `pipeline.R`) are provided in the `input` folder. The only input files required are the .gpr files


## Running the pipeline
Running the StartKit scripts can be done by either in RStudio or by running the scripts from the command line, *keeping in my that the scripts must be places within the same parent folder as the input and output folders!*

```shell
Rscript QC.R
Rscript pipeline.R
```

### Output Files
Once the scripts have been run the output folder will contain the resulting QC and pipeline outputs. For example the following plot shows that the duplicate spots on the slide rarely diverge extemely in values (i.e. beyond a fold-change of 2):

![alt text](https://github.com/cdi-lab/SimplerStarterKit/blob/main/output/red_outs.png "divergent replicates are rare")
