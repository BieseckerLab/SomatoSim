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

## Tutorial

1. Download and install dependencies:

Make sure you have the latest version of SAMtools and pip installed.

Refer to the SAMtools documentation (http://www.htslib.org/download/) for SAMtools installation instructions.

When you install SomatoSim, pip will also automatically download and install any python packages you need to use SomatoSim. You can upgrade your current version of pip using the following command:

```text
$ python3 -m pip install --upgrade pip
```

2. Download and install SomatoSim:


```text
$ git clone https://github.com/BieseckerLab/SomatoSim.git
$ cd SomatoSim
$ python3 -m pip install .
```

3. Run SomatoSim:

Here we will simulate 100 SNVs with VAFs between 0.01 and 0.05. The input BAM and BED files are supplied in this SomatoSim repository. 

```text
$ somatosim \
  -i test_data/test_BAM.bam \
  -b test_data/test_BED.bed \
  -o ./output_dir \
  --vaf-low 0.01 \
  --vaf-high 0.05 \
  --number-snv 100 \
  --random-seed 0
```

4. Explore SomatoSim's outputs

The output_dir directory should contain 4 files: 
* A new BAM file: `test_BAM.somatosim.bam`
* A new BAM index file: `test_BAM.somatosim.bam.bai`
* A results file: `simulation_output.txt`
* A log file: `simulation_log.txt`

Below is an example of the first 10 lines of the results file. The first 3 columns are the SNV coordinates, the `input_VAF` column is the target VAF, the `input_coverage` column is the original coverage at the position, the `output_VAF` column is the VAF calculated after making the mutations, the `output_coverage` column is the coverage after making the mutations (and should be the same as the `input_VAF` column), the `ref_allele` column is the original allele at the position, and the `alt_allele` is the new allele introduced by SomatoSim.

```text
chromosome position position input_VAF input_coverage output_VAF output_coverage ref_allele alt_allele
1           109898015   109898016   0.05        329         0.04863     329         C           A          
1           211571650   211571651   0.05        285         0.04912     285         C           T          
1           14042085    14042086    0.04        380         0.03947     380         C           T          
1           156127891   156127892   0.05        269         0.04833     269         G           T          
1           198222215   198222216   0.05        399         0.05013     399         G           C          
1           197009494   197009495   0.02        335         0.0209      335         A           G          
1           153148814   153148815   0.01        289         0.01038     289         T           A          
1           92729264    92729265    0.03        425         0.03059     425         T           G          
2           89999382    89999383    0.03        242         0.02893     242         T           G 
```

Below is an example of the log file. It reports which versions of the dependencies are currently in use, the user's input parameters, and metrics related to the three stages (variant selection, variant simulation, and variant evaluation) of SomatoSim. 

```text
  ____                        _       ____  _        
 / ___|  ___  _ __ ___   __ _| |_ ___/ ___|(_)_ __ ___  
 \___ \ / _ \| '_ ` _ \ / _` | __/ _ \___ \| | '_ ` _ \ 
  ___) | (_) | | | | | | (_| | || (_) |__) | | | | | | | 
 |____/ \___/|_| |_| |_|\__,_|\__\___/____/|_|_| |_| |_| 
 
Mon Oct 26 16:43:06 2020
 
SomatoSim version 1.0.0
Python version 3.6.10
Numpy version 1.19.2
Pandas version 1.1.3
Matplotlib version 3.3.2
Pysam version 0.16.0.1
Samtools version 1.11
 
Log file for: test_BAM
 
Simulation run parameters:
--------------------------
BAM file path: test_data/test_BAM.bam
BED file path: test_data/test_BED.bed
BED file input contains no pre-defined VAF or alternate allele
 
VAF range: 0.01 - 0.05
Output directory: ./output_dir
Using default output prefix
Number of positions to mutate: 100
Random seed: 0
Samtools mpileup coverage calculation using mapping quality (MQ) = 20 and base quality (BQ) = 20
Minimum base pair separation = 1
Verbose mode: Off
 
Variant selection stage set-up:
-------------------------------
Using input BED format
 
Average coverage of original BAM file: 289.784
 
 
Variant selection stage metrics: 
--------------------------------
Total number of genomic positions in the BED file: 159709
Number of genomic ranges in the BED file: 288
Total number of genomic positions iterated through: 137
Final number of input positions selected to mutate: 100
Average coverage of the final selected input positions: 284.2
 
VAF distribution (input):
[[ 0.    9.  ]
 [ 0.01 20.  ]
 [ 0.02 17.  ]
 [ 0.03 26.  ]
 [ 0.04 20.  ]
 [ 0.05  8.  ]
 [ 0.06  0.  ]]
 
Alternate allele distribution (input):
[['T' '26']
 ['G' '27']
 ['A' '24']
 ['C' '23']]
 
Variant selection stage runtime: 9.046 seconds
 
Variant simulation stage metrics:
---------------------------------
Minimum mapping quality for reads to mutate: 20
Minimum base quality for variant allele: 20
 
Number of mutations expected from the variant selection stage: 880
Number of reads selected for variant simulation: 880
Total number of mutations: 880
Total reads in the input BAM file: 341240
 
Variant simulation stage runtime: 44.823 seconds
 
Variant evaluation stage metrics:
-------------------------------
Final number of unique spike-in positions inside the target VAF value range: 100
Final mean coverage at spike-in positions inside the target VAF value range: 284.2
 
VAF distribution (output):
[[ 0.    9.  ]
 [ 0.01 20.  ]
 [ 0.02 17.  ]
 [ 0.03 26.  ]
 [ 0.04 20.  ]
 [ 0.05  8.  ]
 [ 0.06  0.  ]]
 
Alternate allele distribution (output):
[['T' '26']
 ['G' '27']
 ['A' '24']
 ['C' '23']]
 
####################
Simulation complete!
####################
 
Variant evaluation stage runtime: 1.814 seconds
 
Total runtime: 55.684 seconds
```

## Attribution
This is software was developed by Marwan A. Hawari, Celine S. Hong, and Leslie G. Biesecker at the National Human Genome Research Institute (NHGRI), National Institutes of Health (NIH). Please include proper attribution of the NHGRI as the developer of this program and include a link to the following [https://github.com/BieseckerLab/SomatoSim] in all publications and other public disclosures that reference the program and/or include data or research results that were generated utilizing the program.

## Public Domain Notice
This software is a United States Government Work. Anyone may use the software on a worldwide and royalty-free basis for any purpose and anyone may reproduce and prepare derivative works without restriction. Although all reasonable efforts have been taken to ensure the accuracy and reliability of the software, the National Human Genome Research Institute (NHGRI), National Institutes of Health (NIH) and the U.S. Government do not and cannot warrant the performance or any results that may be obtained by using this software. NHGRI, NIH and the U.S. Government disclaim all warranties as to performance, merchantability or fitness for any particular purpose. No indemnification is intended or provided by the US government.

