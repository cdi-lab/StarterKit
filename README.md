# StarterKit
A set of simple R scripts for basic, stand-alone analysis of HuProt (gpr) result files.


## Overview
The following StartKit provides standalone scripts for the generation of basic QC plots and Analytic Results tables from raw gpr files provided by CDI Lab. They constitute a subset of the techniques used in-house and are provided for the benefit of the users who would like to independently (re)-produce the core results in CDI-provided reports.


## Required .tsv files
The inputs to a given StartKit script (of which there are two, namely `QC.R` and `pipeline.R`) are provided in the `input` folder. The key input files are, naturally, the .gpr files, but there are two additional files that will be expected. The first is `channels.txt` which specifies the channels used in the experiment. As simple example of this file is shown below:

| Label       | Secondary |
| ------------- | ------------- |
| red   | IgG |
| green   | IgA |

An additional file may be provided which specifies the A-vs-B comparisons which will be made by the `pipeline.R` script. An examples is provided in the repository and looks like this:

|Sample|Treatment_vs_Control|Subgroup_vs_Control|
|------|--------------------|-------------------|
|01_CDI7890911.gpr|case||
|02_CDI7890912.gpr|case||
|03_CDI7890972.gpr|case||
|04_CDI7890997.gpr|case|case|
|05_CDI7890974.gpr|case||  
|06_CDI7890975.gpr|case||  
|07_CDI7890982.gpr|case||
|08_CDI7890983.gpr|case|case|
|09_CDI7890984.gpr|case|case|
|10_CDI7890985.gpr|case||
|11_CDI7890986.gpr|case|case|
|12_CDI7890987.gpr|case|| 
|13_CDI7890988.gpr|case||
|14_CDI7890989.gpr|case|case|
|15_CDI7890990.gpr|case||
|16_CDI7890991.gpr|control|control|
|17_CDI7890993.gpr|control|control|
|18_CDI7890994.gpr|control|control|
|19_CDI7890995.gpr|control|control|
|20_CDI7890996.gpr|control|control|

## Running the pipeline

```shell
```

### Output Files
