# SomatoSim

## Description
SomatoSim is a tool that simulates single nucleotide variants at any variant allele fraction and depth of coverage.
SomatoSim takes an analysis ready BAM file as input and SomatoSim outputs a BAM file containing simulated variants.
Use SomatoSim to evaluate variant calling tools, algorithms, and pipelines.

## Dependencies
1. Python version 3.6.8 (https://www.python.org/downloads/)
2. Numpy version 1.16.2 (https://numpy.org/)
3. Pandas version 0.25.1 (https://pandas.pydata.org/)
4. Matplotlib version 3.1.1 (https://matplotlib.org/)
5. Pysam version 0.15.0 (https://pysam.readthedocs.io/en/latest/)
6. SAMtools version 1.9 (http://www.htslib.org/)


## Installation
```text
$ git clone https://github.com/BieseckerLab/SomatoSim.git
$ cd SomatoSim
$ python3 -m pip install . 
```

## Usage
```text
# The following are the command line options for SomatoSim
$ somatosim \
  -i <input bam file path> \
  -b <input bed file path> \
  -o <output directory path> \
  --vaf-low <number> \
  --vaf-high <number> \
  --output-prefix <output files prefix> \
  --number-snv <number> \
  --random-seed <number> \
  --down-sample <number> \
  --target-coverage <number> \
  --coverage-tolerance <number> \
  --coverage-MQ <number> \
  --coverage-BQ <number> \
  --minimum-separation <number> \
  <--verbose> \
  --read-min-MQ <number> \
  --position-min-BQ <number>
```

### Attribution
This is software was developed by Marwan A. Hawari, Celine S. Hong, and Leslie G. Biesecker at the National Human Genome Research Institute (NHGRI), National Institutes of Health (NIH). Please include proper attribution of the NHGRI as the developer of this program and include a link to the following [https://github.com/BieseckerLab/SomatoSim] in all publications and other public disclosures that reference the program and/or include data or research results that were generated utilizing the program.

### Public Domain Notice
This software is a United States Government Work. Anyone may use the software on a worldwide and royalty-free basis for any purpose and anyone may reproduce and prepare derivative works without restriction. Although all reasonable efforts have been taken to ensure the accuracy and reliability of the software, the National Human Genome Research Institute (NHGRI), National Institutes of Health (NIH) and the U.S. Government do not and cannot warrant the performance or any results that may be obtained by using this software. NHGRI, NIH and the U.S. Government disclaim all warranties as to performance, merchantability or fitness for any particular purpose. No indemnification is intended or provided by the US government.

