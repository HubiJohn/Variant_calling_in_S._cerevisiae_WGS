# Variant_calling_in_S._cerevisiae_WGS
The aim fo this project was to build a semi-automated pipeline for variant calling of Candian S. cerevisiae

This project was a part of a graded assignment for the "Analysis of NGS data" course during my MEng studies. 
The main aim was to construct a minimalistic a semi-automated pipeline which result in some quirks described here.
A a great part of the assignment was to assess the biological cosequences of detected variants (SNPs) using open source databases
and present everything to the course group. The presentation translated to English is provided in the depository.

The bioproject used for the assignment is PRJNA838724
(associated article: https://academic.oup.com/g3journal/article/13/8/jkad130/7194280#413067566)

Functioalities and quirks include:
- The main script always creates a directory Project with all its subdirectories.
- When a bioproject number is provided, all FASTQ files associated with it are downloaded.
- When a run number is provided, only FASTQ files associated with it are downloaded
- Adapter sequences for trimming are defined within the bash code which is not optimal and against best programming practices.

What could be improved:
- conda environments could be utilised
- whole script could be broken up into smaller scripts that are called by the a main script
- the main script could take more arguments, e.g.:
    -  the name for the project directory,
    -  adapter seq files,
    -  a file with arguments for all the steps of the pipeline, 
- functionality that would allow to start the analysis from a certain step and stop the analysis at any failed step
