# Tempus Coding Challenge
Having trouble running any part of this code? Please contact me at bbowles1@alumni.nd.edu. In a hurry? Run the script in "test" mode (`python vcf_annot.py -i <input_file> -t`) for faster code execution.

**The challenge:** Prototype a variant interpretation tool to label the following.
1. Depth of sequence coverage at the site of variation.
2. Number of reads supporting the variant.
3. Percentage of reads supporting the variant versus those supporting reference reads.
4. Using the VEP hgvs API, get the:
  + gene of the variant
  + type of variation (substitution, insertion, CNV, etc.)
  + variant effect (missense, silent, intergenic, etc.)
5. The minor allele frequency of the variant if available.

This repository contains the following files:
+ *vcf_annot.py*: The main annotation script for VCF 4.0 files.
+ *annot.tsv*: An annotated list of variants initially given in the `test_vcf_data.txt`
+ *Tempus_Demo.ipynb*: a Jupyter notebook walking through some of the design decisions for the variant interpreter.
+ *Tempus_Demo.html*: saved .html output from Tempus_Demo.ipynb.
+ *env.yml*: a full list of the dev environment packages for reference.

**Dependencies**
This code is designed to run with minimal dependencies. You will require an installation of the `requests` and `pandas` python libraries (preferably pandas 1.4). A full list of dependencies and my python version can be found in the env.yml file.

**Running the Annotation Pipeline**

The general syntax for annotating an input VCF file is to run `python vcf_annot.py -i <input_file> -o <output_file>`. Users can also specify the following options:
-t: Run in "test mode," where a random subset of 250 input variants is selected to demo the pipeline functionality. 
-v: Verbose. If verbose, print summary statistics for read depth, variant supporting reads, and variant supporting read percent.

**Notes**
+ This script requires VCF 4.0 format as input with 'TR' and 'TC' fields in the INFO column.
+ This script slightly expands dataframe size versus the initial input, as variants A) map to multiple unique genes, and B) some input variants contain multiple ALT calls at a single locus, which are expanded to separate rows.