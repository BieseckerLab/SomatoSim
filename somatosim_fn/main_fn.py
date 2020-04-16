#!/usr/bin/env python3

global np
import numpy as np
global pd
import pandas as pd
global matplotlib
import matplotlib
global pysam
import pysam
global os
import os
global sys
import sys
global time
import time
global datetime
import datetime

global print_both,convert_BED,vaf_alt_distribution,count_alleles,select_reads
from somatosim_fn.minor_fn import print_both,convert_BED,vaf_alt_distribution,count_alleles,select_reads,check_BED_columns

def variant_selection(bam_directory_in,bed_in,output_directory,vaf_low_in,vaf_high_in,output_prefix,number_positions_to_mutate,seed_in,desired_downsampled_coverage,target_coverage,coverage_tolerance,pileup_MQ,pileup_BQ,min_separation,verbose_mode_in):

    np.set_printoptions(threshold=sys.maxsize)

    global start_time
    start_time=time.time()

    global SomatoSim_version
    SomatoSim_version = "1.0.0"

    #define bam file path (should also contain the bai file)
    global rawbamname
    rawbamname = bam_directory_in

    if not os.path.exists(rawbamname):
        raise Exception('BAM file path does not exist')


    #extract the prefix of the bam file
    global bamprefix
    bamprefix = rawbamname.split('/')[-1][:-4]

    #define bed file path
    global bedfilepath
    bedfilepath = bed_in

    if not os.path.exists(bedfilepath):
        raise Exception('BED file path does not exist')

    #define simulation output directory 
    global fullrunpath_user
    fullrunpath_user = output_directory


    #Declare the bounds for VAF
    global vaf_low
    global vaf_high

    if not (0 < vaf_low_in < 1) or not (0 < vaf_high_in < 1):
        raise ValueError('Use decimals for VAF, not percentages')
    else:
        vaf_low = vaf_low_in*100
        vaf_high = vaf_high_in*100 + 1

    #Declare the random seed
    global randseed
    randseed = seed_in


    #Define the MQ and BQ to be used in samtools mpileup coverage calculation
    global coverage_MQ
    coverage_MQ = pileup_MQ # q
    global coverage_BQ
    coverage_BQ = pileup_BQ # Q

    #Define the prefix for the output file
    global out_prefix
    out_prefix = output_prefix


    #if choosing to downsample the BAM file, define the desired down-sampled coverage
    global desired_ds_cov
    desired_ds_cov = desired_downsampled_coverage


    #define the target coverage
    target_cov = target_coverage


    #Declare tolerance on the target coverage:
    global cov_tolerance
    cov_tolerance = coverage_tolerance


    #Define the minimum number of bases in between simulated variants
    global minimum_separation
    minimum_separation = min_separation
    
    global verbose_mode
    verbose_mode = verbose_mode_in

    global haplo_status
    haplo_status = False

    ################################################################################################################
    ########################################## END OF DECLARATIONS #################################################
    ################################################################################################################


    #Define the full run path
    global fullrunpath
    if fullrunpath_user[-1] == '/':
        fullrunpath = fullrunpath_user[:-1]
    else:
        fullrunpath = fullrunpath_user
    #Create a temporary directory
    global tmpdir
    if not os.path.exists(fullrunpath):
        os.system('mkdir '+ fullrunpath)
    if out_prefix:
        tmpdir = fullrunpath + '/' + out_prefix + '_tmp'
    else:
        tmpdir = fullrunpath + '/tmp'
    if not os.path.exists(tmpdir):
        os.system('mkdir '+ tmpdir)

    #Create a log file
    if out_prefix:
        if os.path.exists(fullrunpath + '/' + out_prefix + '_simulation_log.txt'):
            os.system('rm ' + fullrunpath + '/' + out_prefix + '_simulation_log.txt')
        logfile = open(fullrunpath + '/' + out_prefix + '_simulation_log.txt','w')
        out_prefix_status = True
    else:
        if os.path.exists(fullrunpath + '/simulation_log.txt'):
            os.system('rm ' + fullrunpath + '/simulation_log.txt')
        logfile = open(fullrunpath + '/simulation_log.txt','w')
        out_prefix_status = False




    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$     


    global BED
    global oriBED
    # read in target bed file
    readbed = pd.read_csv(bedfilepath, delim_whitespace=True, header=None, dtype=str)
    BED = readbed.values
    oriBED = readbed.values


    #Create many global variables
    global originalbamname
    global newbamname
    global newbamname_sorted
    global avgcov_value
    global bed_select
    global forpileup


    #Define the number of positions to simulate
    global number_to_mutate
    if number_positions_to_mutate == 0:
        number_to_mutate = len(BED)
    else:
        number_to_mutate = number_positions_to_mutate



    #if the original BED file has a third column
    global user_vaf
    global user_alt
    global user_vaf_status
    global user_alt_status




    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 

    print_both("  ____                        _       ____  _        ",logfile)
    print_both(" / ___|  ___  _ __ ___   __ _| |_ ___/ ___|(_)_ __ ___  ",logfile)
    print_both(" \___ \ / _ \| '_ ` _ \ / _` | __/ _ \___ \| | '_ ` _ \ ",logfile)
    print_both("  ___) | (_) | | | | | | (_| | || (_) |__) | | | | | | | ",logfile)
    print_both(" |____/ \___/|_| |_| |_|\__,_|\__\___/____/|_|_| |_| |_| ",logfile)

    print_both(" ",logfile)

    print_both(datetime.datetime.now().strftime("%c"),logfile)

    print_both(" ",logfile)

    print_both("SomatoSim version " + SomatoSim_version,logfile)
    print_both("Python version " + sys.version[:6],logfile)

    print_both("Numpy version " + np.__version__,logfile)
    print_both("Pandas version " + pd.__version__,logfile)
    print_both("Matplotlib version " + matplotlib.__version__,logfile)

    print_both("Pysam version " + pysam.__version__,logfile)
    print_both("Samtools version " + os.popen('samtools --version').read().split()[1],logfile)
    print_both(" ",logfile)

    print('Log file for: ' + str(bamprefix),file=logfile)
    print(" ",file=logfile)


    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    #Print the user's arguments to stdout and to the log file

    print("Started simulation run:")
    print("-----------------------")
    print("User inputs:")

    print("Simulation run parameters:", file=logfile)
    print("--------------------------", file=logfile)

    print_both("BAM file path: " + rawbamname,logfile)
    print_both("BED file path: " + bedfilepath,logfile)

    #*****************************
    user_vaf,user_alt,user_vaf_status,user_alt_status = check_BED_columns(oriBED,logfile)
    #*****************************


    if user_vaf_status == True:
        print_both("VAF range: " + str(float(np.around(min(user_vaf.astype(float)),decimals=5))) + " - " + str(float(np.around(max(user_vaf.astype(float)),decimals=5))),logfile)
    else:
        print_both("VAF range: " + str(np.around((vaf_low/100),decimals=5)) + " - " + str(np.around((vaf_high)/100-0.01,decimals=5)),logfile)
    print_both("Output directory: " + fullrunpath,logfile)

    if out_prefix:
        print_both("Output prefix: " + out_prefix,logfile)
    else:
        print_both("Using default output prefix",logfile)
    print_both("Number of positions to mutate: " + str(number_to_mutate),logfile)
    print_both("Random seed: " + str(randseed),logfile)
    if desired_ds_cov != 0 and type(desired_ds_cov) == int:
        print_both("Desired downsampled coverage: " + str(desired_ds_cov),logfile)
    if target_cov != 0 and type(target_cov) == int:
        print_both("Target coverage: " + str(target_cov),logfile)
        print_both("Coverage tolerance: " + str(cov_tolerance),logfile)

    print_both("Samtools mpileup coverage calculation using mapping quality (MQ) = " + str(coverage_MQ) + " and base quality (BQ) = " + str(coverage_BQ),logfile)
    print_both("Minimum base pair separation = " + str(minimum_separation),logfile)
    if verbose_mode == True:
        print_both("Verbose mode: On",logfile)
    else:
        print_both("Verbose mode: Off",logfile)

    print_both(" ",logfile)
    print("Started variant selection stage:")
    print("--------------------------------")
    print("Variant selection stage set-up:", file=logfile)
    print("-------------------------------", file=logfile)




    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$     



    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    if desired_ds_cov == 0:
        #if NOT downsampling, then calculate the average coverage
        originalbamname = rawbamname
        try:
            #try using the format of BED that was originally given
            avgcov_value = os.popen("samtools mpileup -x -q " + str(coverage_MQ) + " -Q " + str(coverage_BQ) + " -l "+ bedfilepath + " " + rawbamname + " 2> /dev/null | awk -F '\t' '{sum+=$4} END { print sum/NR}' 2>/dev/null").read()[:-1]
            attempt = str(float(avgcov_value))
            print_both("Using input BED format",logfile)
            print_both(" ",logfile)

        except:
            #but if it doesn't work with the given BAM file, then use the converted format BED file
            new_BED = convert_BED(BED,tmpdir)
            BED = new_BED
            avgcov_value = os.popen("samtools mpileup -x -q " + str(coverage_MQ) + " -Q " + str(coverage_BQ) + " -l "+ tmpdir + "/new_BED.bed " + rawbamname + " 2> /dev/null | awk -F '\t' '{sum+=$4} END { print sum/NR}' 2>/dev/null").read()[:-1]
            attempt = str(float(avgcov_value))
            print_both("Using converted BED format",logfile)
            print_both(" ",logfile)

        print_both('Average coverage of original BAM file: ' + str(float(avgcov_value)), logfile)
        print_both(" ",logfile)

        if out_prefix:
            newbamname = fullrunpath + "/" + out_prefix + ".somatosim.notsorted.bam"
            newbamname_sorted = fullrunpath + "/" + out_prefix + ".somatosim.bam"
        else:
            newbamname = fullrunpath + "/" + bamprefix + ".somatosim.notsorted.bam"
            newbamname_sorted = fullrunpath + "/" + bamprefix + ".somatosim.bam"


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    elif desired_ds_cov != 0 and type(desired_ds_cov) == int:

        try:
            #try using the format of BED that was originally given
            avgcov_value = os.popen("samtools mpileup -x -q " + str(coverage_MQ) + " -Q " + str(coverage_BQ) + " -l "+ bedfilepath + " " + rawbamname + " 2> /dev/null | awk -F '\t' '{sum+=$4} END { print sum/NR}' 2>/dev/null").read()[:-1]
            if out_prefix:
                originalbamname = fullrunpath + "/" + out_prefix + '.ds.bam'
            else:
                originalbamname = fullrunpath + "/" + bamprefix + '.ds.bam'
            #define desired output coverage for BAM and then downsample the BAM
            ds_percentage_multiplier = int(round(desired_ds_cov/float(avgcov_value) * 100))
            #bash command to downsample BAM
            if len(str(ds_percentage_multiplier)) < 2:
                avgcov_value_ds = os.popen("samtools view -bs " + str(randseed) + ".0" + str(ds_percentage_multiplier) + " " + rawbamname + " > " + originalbamname + " ; samtools index " + originalbamname + " ; samtools mpileup -x -q " + str(coverage_MQ) + " -Q " + str(coverage_BQ) + " -l " + bedfilepath + " " + originalbamname + " 2> /dev/null | awk -F '\t' '{sum+=$4} END { print sum/NR}' 2>/dev/null").read()[:-1]    
            else:
                avgcov_value_ds = os.popen("samtools view -bs " + str(randseed) + "." + str(ds_percentage_multiplier) + " " + rawbamname + " > " + originalbamname + " ; samtools index " + originalbamname + " ; samtools mpileup -x -q " + str(coverage_MQ) + " -Q " + str(coverage_BQ) + " -l " + bedfilepath + " " + originalbamname + " 2> /dev/null | awk -F '\t' '{sum+=$4} END { print sum/NR}' 2>/dev/null").read()[:-1]     
            # Calculate average coverage for the bamfile of interest after downsampling
            attempt = str(float(avgcov_value))
            attempt_ds = str(float(avgcov_value_ds))
            print_both("Using input BED format",logfile)
            print_both(" ",logfile)

        except:
            new_BED = convert_BED(BED,tmpdir)
            BED=new_BED
            #but if it doesn't work with the given BAM file, then use the converted format BED file
            avgcov_value = os.popen("samtools mpileup -x -q " + str(coverage_MQ) + " -Q " + str(coverage_BQ) + " -l "+ tmpdir + "/new_BED.bed " + rawbamname + " 2> /dev/null | awk -F '\t' '{sum+=$4} END { print sum/NR}'").read()[:-1]
            if out_prefix:
                originalbamname = fullrunpath + "/" + out_prefix + '.ds.bam'
            else:
                originalbamname = fullrunpath + "/" + bamprefix + '.ds.bam'
            #define desired output coverage for BAM and then downsample the BAM
            ds_percentage_multiplier = int(round(desired_ds_cov/float(avgcov_value) * 100)) 
            #bash command to downsample BAM
            if len(str(ds_percentage_multiplier)) < 2:
                avgcov_value_ds = os.popen("samtools view -bs " + str(randseed) + ".0" + str(ds_percentage_multiplier) + " " + rawbamname + " > " + originalbamname + " ; samtools index " + originalbamname + " ; samtools mpileup -x -q " + str(coverage_MQ) + " -Q " + str(coverage_BQ) + " -l "+ tmpdir + "/new_BED.bed " + originalbamname + " 2> /dev/null | awk -F '\t' '{sum+=$4} END { print sum/NR}' 2>/dev/null").read()[:-1]    
            else:
                avgcov_value_ds = os.popen("samtools view -bs " + str(randseed) + "." + str(ds_percentage_multiplier) + " " + rawbamname + " > " + originalbamname + " ; samtools index " + originalbamname + " ; samtools mpileup -x -q " + str(coverage_MQ) + " -Q " + str(coverage_BQ) + " -l "+ tmpdir + "/new_BED.bed " + originalbamname + " 2> /dev/null | awk -F '\t' '{sum+=$4} END { print sum/NR}' 2>/dev/null").read()[:-1]    
            # Calculate average coverage for the bamfile of interest after downsampling
            attempt = str(float(avgcov_value))
            attempt_ds = str(float(avgcov_value_ds))
            print_both("Using converted BED format",logfile)
            print_both(" ",logfile)

        print_both('Average coverage of original BAM file: '+ str(float(avgcov_value)),logfile)
        print_both('Average coverage of down-sampled BAM file: '+ str(float(avgcov_value_ds)),logfile)
        print_both(" ",logfile)

        if out_prefix:
            newbamname = fullrunpath + "/" + out_prefix + ".ds.somatosim.notsorted.bam"
            newbamname_sorted = fullrunpath + "/" + out_prefix + ".ds.somatosim.bam"
        else:
            newbamname = fullrunpath + "/" + bamprefix + ".ds.somatosim.notsorted.bam"
            newbamname_sorted = fullrunpath + "/" + bamprefix + ".ds.somatosim.bam"


    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    ################################################################################################################
    ########################################## START OF VARIANT SELECTION ##########################################
    ################################################################################################################




    no_coverage_bucket = []         
    allele_status_bucket = []
    in_target_cov_bucket = []
    too_close_bucket = []
    ref_alt_status_bucket = []
    var_nonzero_bucket = []


    #Define and initialize variables
    BED_positions = np.array(BED[:,1:3])
    bases_per_range = (BED_positions[:,1].astype(int) - BED_positions[:,0].astype(int))
    total_number_bases = np.sum(bases_per_range)
    potential = []
    potential_range = []
    #k represents the line in the BED file = the genomic range in the BED file
    k=0
    #q represents the bases iterated through
    q=0
    #d represents the number of BED ranges that have been exhausted
    d=0
    rando_bank = []
    np.random.seed(randseed)
    possible_allele = np.array(['T','G','A','C'])
    exhausted = []
    run_count = np.zeros((len(BED),1))

    #Create a data structure that will enable selection of random positions from each BED genomic range
    random_position_tracker = np.zeros((len(BED),1))
    entire_BED_random_positions = []
    for r in range(len(BED)):
        np.random.seed(randseed)
        total_positions_BED = np.linspace(float(BED[r,1]),float(BED[r,2]),int(float(BED[r,2])-float(BED[r,1]))+1)[1:]
        total_positions_BED_random = np.random.choice(total_positions_BED,int(float(BED[r,2])-float(BED[r,1])),replace=False)
        entire_BED_random_positions.append(total_positions_BED_random)

    print("Started selecting positions for simulation")

    if number_to_mutate < total_number_bases:
        length_checker = number_to_mutate
    else:
        length_checker = total_number_bases


    #Start while loop to search for base positions that are suitable for variant simulation
    while True:

        #if the BED file range has not been exhausted 
        if k not in exhausted:
            pass
        #if the BED file range is exhausted, then move on to the next BED file range, and if it is out of range, then reset and go back to the first line
        elif k in exhausted and d != len(BED):
            k+=1
            if k>=len(BED):
                k=0
            continue
        #if you've exhausted the current BED file range and the number of exhausted BED file ranges = total ranges in the BED file, then break out of loop -- variant selection is done
        elif k in exhausted and d == len(BED):
            break

        #select the position from the array of randomly distributed positions
        selected_range = np.array(entire_BED_random_positions)[k]
        rando = int(selected_range[int(random_position_tracker[k])])
        random_position_tracker[k] +=1

        #check if this selection will exhaust the BED file range, if it is, then add the index of the bed file range to an array called "exhausted"
        num_line_positions = int(BED[k,2])-int(BED[k,1])
        if int(random_position_tracker[k]) == num_line_positions:
            exhausted.append(k)
            d+=1


        #Assign a vaf, check whether or not the user has input a vaf, otherwise just randomly assign a vaf
        if user_vaf_status == True:
            vaf = float(user_vaf[k][0])
        else:
            #if both inputs are integer percentages, then only assign integer percentage VAFS, otherwise can assign non-integer percentage VAF
            if int(vaf_low) == vaf_low and int(vaf_high) == vaf_high:
                vaf = np.random.randint(vaf_low,vaf_high)/100 #**********#
            else:
                vaf = np.random.uniform(vaf_low, vaf_high-1)/100


        if "chr" in str(BED[k,0]):
            line = str(str(BED[k,0])) + ":" + str(int(float(rando))) + "-" + str(int(float(rando)))
        else:
            try:
                line = str(str(int(float(BED[k,0]))) + ":" + str(int(float(rando))) + "-" + str(int(float(rando))))
            except:
                line = str(str(BED[k,0])) + ":" + str(int(float(rando))) + "-" + str(int(float(rando)))


        #Calculate depth of coverage
        value_singlecovval = os.popen('samtools mpileup -x -q ' + str(coverage_MQ) + ' -Q ' + str(coverage_BQ) + ' -r ' + line + " " + originalbamname + ' 2> /dev/null | cut -f 4').read()[:-1]

        #check if that coverage value exists
        if len(value_singlecovval) == 0:
            value_singlecovval = 0
            allelecount = [0,0,0,0,0]

        else:
            allelecount = count_alleles(coverage_MQ,coverage_BQ,line,originalbamname)


        #check if monoallelic
        allele_status = np.size(np.array(allelecount[:4]).nonzero()) == 1

        #check if position is in coverage range
        if value_singlecovval != 0:
            #if the user selected a target coverage value that is an integer
            if target_cov != 0 and type(target_cov) == int:
                #if the coverage of this position is within that target ccoverage and tolerance, then add a "True" to the end of the array
                if int(value_singlecovval) < float(target_cov)*(1 + cov_tolerance) and int(value_singlecovval) > float(target_cov)*(1 - cov_tolerance):
                    in_target_cov = True
                #otherwise, add a False
                else:
                    in_target_cov = False
            #if the user did not select a target coverage, then just add a True anyway
            else:
                in_target_cov = True
        else:
            in_target_cov = False

        #assign ref and alternate allele. If user has decided the alternate, then use that, else figure out the reference and randomly assign an alt that is not equal to reference
        if user_alt_status == True:
            alt_allele = user_alt[k][0]
            allelecount_poly_filt = np.array(allelecount)[:4].astype(float)
            if allele_status == True:
                allelenumeric = allelecount_poly_filt.nonzero()[0]
                if int(allelenumeric) == 0:
                    ref_allele = 'T'
                if int(allelenumeric) == 1:
                    ref_allele = 'G'
                if int(allelenumeric) == 2:
                    ref_allele = 'A'
                if int(allelenumeric) == 3:
                    ref_allele = 'C'
            else:
                allelenumeric = 0
                ref_allele = 'T'

        else:
            allelecount_poly_filt = np.array(allelecount)[:4].astype(float)
            if allele_status == True:
                allelenumeric = allelecount_poly_filt.nonzero()[0]
                if int(allelenumeric) == 0:
                    alt_allele = np.random.choice(possible_allele[1:])
                    ref_allele = 'T'
                if int(allelenumeric) == 1:
                    alt_allele = np.random.choice(possible_allele[[0,2,3]])
                    ref_allele = 'G'
                if int(allelenumeric) == 2:
                    alt_allele = np.random.choice(possible_allele[[0,1,3]])
                    ref_allele = 'A'
                if int(allelenumeric) == 3:
                    alt_allele = np.random.choice(possible_allele[[0,1,2]])
                    ref_allele = 'C'

            else:
                allelenumeric = 0
                alt_allele = 'T'
                ref_allele = 'T'


        #check that ref is not the same as the alternate
        if ref_allele != alt_allele:
            ref_alt_status = True
        else:
            ref_alt_status = False


        #assign the number of reads that must be mutated to achieve the randomly assigned VAF

        var_depth = int(round(float(np.array(int(value_singlecovval)) * float(np.array(vaf)))))

        #check that variant depth not <= 0 and that the resulting output VAF = var_depth / total depth will be within the the coverage range  ONLY IF DOING RANDOM MUTATION
        if user_vaf_status == False:  
            if (var_depth != 0) and not (var_depth < 0) and ((vaf_low/100 - 0.005) < (var_depth/int(value_singlecovval)) < (vaf_high/100 - 0.005)):
                var_nonzero = True
            else:
                var_nonzero = False

        elif user_vaf_status == True:
            if (var_depth != 0) and not (var_depth < 0) and (min(user_vaf.astype(float))-0.005 < (var_depth/int(value_singlecovval)) < max(user_vaf.astype(float))+0.005):
                var_nonzero = True
            else:
                var_nonzero = False

        bed_in_cov_with_allele_vardepth = [BED[k,0], #chromosome [0]
                                           rando, #position [1]
                                           rando, #position [2]
                                           vaf, #vaf [3]
                                           int(value_singlecovval), #DOC [4]
                                           allelecount[0], #T count [5]
                                           allelecount[1], #G count [6]
                                           allelecount[2], #A count [7]
                                           allelecount[3], #C count [8]
                                           allelecount[4], #Sum of allele counts (includes deletions) [9]
                                           allele_status, #True if mono-allelic [10]
                                           in_target_cov, #True if within the coverage range. If no target selected, then True. [11]
                                           ref_allele, #reference allele [12]
                                           alt_allele, #alternate allele [13]
                                           ref_alt_status, #True if ref is not the same as alt [14]
                                           var_depth, #Number of reads to mutate in order to achieve the desired VAF [15]
                                           var_nonzero] #True if var_depth is >=0 and (var_depth/DOC) is within the input VAF range (minVAF-0.005 to maxVAF + 0.005) [16]

        positions_range = np.arange(rando-minimum_separation,rando+minimum_separation+1)
        positions_range_with_chrom = []
        for min_sep_position in range(len(positions_range)):
            positions_range_with_chrom.append(str(BED[k,0])+":"+str(positions_range[min_sep_position]))

        if allele_status == True and in_target_cov == True and ref_alt_status == True and var_nonzero == True and (str(BED[k,0])+":"+str(rando)) not in np.array(potential_range).reshape(1,-1).tolist()[0]:
            potential.append(bed_in_cov_with_allele_vardepth)
            potential_range.append(positions_range_with_chrom)
            if len(potential) == int(length_checker/10):
                print("    Selected 10% of desired positions...")
            if len(potential) == int(length_checker/4):
                print("    Selected 25% of desired positions...")
            if len(potential) == int(length_checker/2):
                print("    Selected 50% of desired positions...")
            if len(potential) == int(3*length_checker/4):
                print("    Selected 75% of desired positions...")
            if len(potential) == int(9*length_checker/10):
                print("    Selected 90% of desired positions...")
            if len(potential) == int(length_checker):
                print("    Selected 100% of desired positions...")

        else:
            if value_singlecovval == 0:
                no_coverage_bucket.append(str(BED[k,0])+":"+str(rando-1)+"-"+str(rando))
            elif value_singlecovval != 0:
                if allele_status != True:
                    allele_status_bucket.append(str(BED[k,0])+":"+str(rando-1)+"-"+str(rando))
                elif in_target_cov != True:
                    in_target_cov_bucket.append(str(BED[k,0])+":"+str(rando-1)+"-"+str(rando))
                elif (str(BED[k,0])+":"+str(rando)) in np.array(potential_range).reshape(1,-1).tolist()[0]:
                    too_close_bucket.append(str(BED[k,0])+":"+str(rando-1)+"-"+str(rando))
                elif ref_alt_status != True:
                    ref_alt_status_bucket.append(str(BED[k,0])+":"+str(rando-1)+"-"+str(rando))
                elif var_nonzero != True:
                    var_nonzero_bucket.append(str(BED[k,0])+":"+str(rando-1)+"-"+str(rando))


        k += 1
        q += 1

        #reset position in BED file if reach the end
        if k >= len(BED):
            k=0

        #break if you get the number you were looking for
        if len(potential) == number_to_mutate:
            break

        #break if you have iterated through every base in the BED file
        if q == total_number_bases:
            break

    if verbose_mode == True:
        if len(no_coverage_bucket) != 0:
            print_both(" ",logfile)
            print_both("WARNING: The following " + str(len(no_coverage_bucket)) + " positions were not selected because there are no reads in the BAM file at this position: ",logfile)
            print_both(np.array(no_coverage_bucket),logfile)
        if len(allele_status_bucket) != 0:
            print_both(" ",logfile)
            print_both("WARNING: The following " + str(len(allele_status_bucket)) + " positions were not selected because there was an existing variant: ",logfile)
            print_both(np.array(allele_status_bucket),logfile)
        if len(in_target_cov_bucket) != 0:
            print_both(" ",logfile)
            print_both("WARNING: The following " + str(len(in_target_cov_bucket)) + " positions were not selected because it was outside the target DOC: ",logfile)
            print_both(np.array(in_target_cov_bucket),logfile)
        if len(too_close_bucket) != 0:
            print_both(" ",logfile)
            print_both("WARNING: The following " + str(len(too_close_bucket)) + " positions were not selected because it was too close to a previously selected position: ",logfile)
            print_both(np.array(too_close_bucket),logfile)
        if len(ref_alt_status_bucket) != 0:
            print_both(" ",logfile)
            print_both("WARNING: The following " + str(len(ref_alt_status_bucket)) + " positions were not selected because the user-specified alternate allele is the same as the reference: ",logfile)
            print_both(np.array(ref_alt_status_bucket),logfile)
        if len(var_nonzero_bucket) != 0:
            print_both(" ",logfile)
            print_both("WARNING: The following " + str(len(var_nonzero_bucket)) + " positions were not selected because it would require 0 reads to be mutated for the selected VAF and DOC: ",logfile)
            print_both(np.array(var_nonzero_bucket),logfile) 


    #rename variant depth and alternate alleles 
    VAR_DEPTH_FIX = np.reshape(np.array(potential)[:,15],(-1,1))

    ref_allele_reshape = np.reshape(np.array(potential)[:,12],(-1,1))
    alt_allele_reshape = np.reshape(np.array(potential)[:,13],(-1,1))


    #Round the VAF (important for non-integer VAF percentages) 
    rounded_input_VAF = np.reshape(np.around(np.array(potential)[:,3].astype(float),decimals=5),(-1,1))
    reshaped_DOC = np.reshape(np.array(potential)[:,4],(-1,1))
    simulation_input_array_rounded_0based = np.concatenate((np.array(potential)[:,0].reshape(-1,1),
                                                            (np.array(potential)[:,1].astype(int)-1).reshape(-1,1),
                                                            np.array(potential)[:,2].reshape(-1,1),
                                                            rounded_input_VAF,
                                                            reshaped_DOC,
                                                            VAR_DEPTH_FIX,
                                                            np.array(potential)[:,5:11],
                                                            ref_allele_reshape,alt_allele_reshape),1)


    #Put everything together into the simulation input array -- this is the final output of the variant selection stage
    global simulation_input_array
    simulation_input_array = np.concatenate((np.array(potential)[:,:5],VAR_DEPTH_FIX,np.array(potential)[:,5:11],ref_allele_reshape,alt_allele_reshape),1)



    print(" ",file=logfile)
    print('Variant selection stage metrics: ', file=logfile)
    print("--------------------------------", file=logfile)
    print("Total number of genomic positions in the BED file:",total_number_bases,file=logfile)
    print('Number of genomic ranges in the BED file:', len(BED), file=logfile)
    print("Total number of genomic positions iterated through:", q,file=logfile)
    print('Final number of input positions selected to mutate:', len(simulation_input_array), file=logfile)
    print('Average coverage of the final selected input positions:',str(np.mean(simulation_input_array[:,4].astype(float))), file=logfile)




    ################################################################################################################
    ########################################## END OF VARIANT SELECTION ############################################
    ################################################################################################################



    vaf_alt_distribution(simulation_input_array[:,5].astype(float)/simulation_input_array[:,4].astype(float),
                         simulation_input_array[:,13],
                         logfile,
                         "(input)")



    print(" ")
    print("Completed variant selection stage!")
    print(" ")



    print_both("Variant selection stage runtime: " + str(round((time.time() - start_time),3)) + " seconds", logfile)
    print_both(" ",logfile)


    logfile.close()


def variant_simulation(min_MQ,min_BQ):

    start_time2=time.time()
    if out_prefix:
        logfile = open(fullrunpath + '/' + out_prefix + '_simulation_log.txt','a')
    else:
        logfile = open(fullrunpath + '/simulation_log.txt','a')

    ######## Variant simulation STAGE ##########
    print("Started variant simulation stage:")
    print("---------------------------------")
    #Min MQ when selecting reads
    mapping_quality_score = min_MQ
    print("Minimum mapping quality for reads to mutate:",str(mapping_quality_score))

    #Define the minimum base quality score at the position to check for. Reads with a base quality at the specified position less than this value won't be used.
    base_quality_score = min_BQ
    print("Minimum base quality for variant allele:",str(base_quality_score))

    ################################################################################################################
    ########################################## START OF READ SELECTION #############################################
    ################################################################################################################

    print(" ")
    print("Started searching for suitable reads")
    #Load in the bam file using pysam
    samfile = pysam.AlignmentFile(originalbamname, "rb" )
    #initialize the sense variable which will be used to alternate between selecting positive and negative sense strands
    sense=False
    #initialize a bank of reads that are used to avoid using duplicate reads
    all_reads = []
    #initialize a variable to keep track of which position we are currently working on
    j=0 #j is the position index
    #initialize a list to contain information about reads as well as the reads themselves
    var_series_intermediate2 = []
    #Start by iterating through the positions determined by simulation_input_array
    np.random.seed(randseed)
    for row in simulation_input_array:
        #initialize a variant_depth_counter to keep track of how many reads have been selected to mutate for a given position
        variant_depth_counter=0
        #set variable chrom to the 0th index of the simulation_input_array, which is the chromosome value
        if "chr" in str(row[0]):
            chrom = row[0]
        else:
            try:
                chrom = str(int(float(row[0])))
            except:
                chrom = row[0]
        #set variable start to the 1st index of the simulation_input_array, which is the position start value
        start = float(row[1])
        #set variable var_depth to the 5th index of the simulation_input_array, which is the number of reads to mutate at the given position
        var_depth = float(row[5])
        #Now iterate through all the reads at the given position
        read_bank_original = []
        for read in samfile.fetch(chrom, start - 1, start):
            read_bank_original.append([read,start])
        read_bank = np.array(read_bank_original)[np.random.choice(len(read_bank_original), len(read_bank_original), replace=False)]


        read_bank_pos = []
        read_bank_neg = []
        for st,read in enumerate(read_bank[:,0]):
            try:
                read_refs = read.get_reference_positions(full_length=True)
                start = read_bank[st,1]-1
                read_index = [h for h, s in enumerate(read_refs) if start == s]
                base_quality_check = read.query_qualities[read_index[0]]
            except:
                pass
            #Check if the read is properly paired, not a duplicate, not unmapped, doesn't fail QC, and is not a secondary read
            if read.is_proper_pair and not read.is_duplicate and not read.is_unmapped and not read.is_qcfail and not read.is_secondary and (read.mapping_quality >= mapping_quality_score) and (base_quality_check >= base_quality_score):
                #Check if position is in the read reference positions
                read_refs = read.get_reference_positions(full_length=True)
                #so read_index is the actual index within the read that contains the position of interest (it generally won't be the first position in that read_refs list because generally the position we are mutating occurs somewhere in the middle of a read)
                read_index = [h for h, s in enumerate(read_refs) if start == s]
                if len(read_index)!=0 and len(read_refs) == len(read.query_sequence):
                    if read.is_reverse == False:
                        read_bank_pos.append(read)
                    if read.is_reverse == True:
                        read_bank_neg.append(read)

        select_reads(read_bank_pos,
                         read_bank_neg,
                         var_depth,
                         all_reads,
                         haplo_status,
                         variant_depth_counter,
                         simulation_input_array,
                         j,
                         var_series_intermediate2)


        #register that you worked through this read position
        j+=1

    #convert list to array
    global var_series
    var_series = np.array(var_series_intermediate2)
    #turn the read data and position number into a dictionary such that read data = key and position number = value
    var_dict = dict(var_series[:,[1,2]])


    print("Started writing new BAM file")

    #Create a new sam file to write reads to, use the old bam file as a template
    new_samfile = pysam.AlignmentFile(newbamname, 'wb', template=samfile) # create file to write


    ################################################################################################################
    ########################################## END OF READ SELECTION ###############################################
    ################################################################################################################

    #initialize some counting variables
    bam_read_counter=0
    hit_counter=0
    mutated_counter=0
    unique_read_counter=0

    #Start iterating through ALL reads in the input (down-sampled) BAM file

    for read in samfile.fetch():
        if read in var_dict:
            bank_index = np.where(var_series[:,1]==read)[0]
            for bank_index_position in range(len(bank_index)):
                #save the read and quality in variables for later
                new_read = read
                new_sequence = read.query_sequence
                tempq = read.query_qualities
                #then save the position number (j) (index) of that position
                var_index = var_series[bank_index[bank_index_position],2]
                #then isolate just the one row from the simulation_input_array for that exact position number (j) (index)
                bam_row = simulation_input_array[var_index]
                #set variable start to the 1st index of the simulation_input_array, which is the position start value
                start = float(bam_row[1]) - 1
                #set variable alt to the 12th index of the simulation_input_array, which is the alternate allele value
                alt = bam_row[13]
                #read.get_reference_positions() is an array of the read containing base-pair position number (like the start variable) for every position in that read (so it will be a linear sequence from beginning to end)
                read_refs = read.get_reference_positions(full_length=True)
                #so read_index is the actual index within the read that contains the position of interest (it generally won't be the first position in that read_refs list because generally the position we are mutating occurs somewhere in the middle of a read)
                read_index = [h for h, s in enumerate(read_refs) if start == s]

                #Check if the read_index exists 
                if len(read_index)!=0 and len(read_refs) == len(read.query_sequence) and len(alt) == 1:
                    #Replace query_sequence of the read at the alt_allele position (read_index)
                    new_read.query_sequence = new_sequence[:read_index[0]] + str(alt) + new_sequence[(read_index[0]+1):]

                    #Then re-write read qualities
                    new_read.query_qualities = tempq

                    #MH write the NEW read to the new sam file
                    if bank_index_position == len(bank_index)-1:
                        new_samfile.write(new_read)

                    #keep track of number of positions (j index) replaced 
                    mutated_counter += 1
                    if mutated_counter == int(len(var_series)/10):
                        print("    10% of selected reads were mutated...")
                    if mutated_counter == int(len(var_series)/4):
                        print("    25% of selected reads were mutated...")
                    if mutated_counter == int(len(var_series)/2):
                        print("    50% of selected reads were mutated...")
                    if mutated_counter == int(3*len(var_series)/4):
                        print("    75% of selected reads were mutated...")
                    if mutated_counter == int(9*len(var_series)/10):
                        print("    90% of selected reads were mutated...")
                    if mutated_counter == int(len(var_series)):
                        print("    100% of selected reads were mutated...")


                #Register how many times a read in the input ds BAM file is in our previously generated list of reads to mutate
                hit_counter += 1

            unique_read_counter +=1

        else:
            #if don't hit a read that is to be edited, then write the OLD read to the new sam file
            new_samfile.write(read)



        bam_read_counter +=1    



    #close, sort, and index the new bam file
    samfile.close()
    new_samfile.close()
    pysam.sort('-o',newbamname_sorted,newbamname)
    pysam.index(newbamname_sorted)





    print("Variant simulation stage metrics:", file=logfile)
    print("---------------------------------", file=logfile)
    print("Minimum mapping quality for reads to mutate:",str(mapping_quality_score),file=logfile)
    print("Minimum base quality for variant allele:",str(base_quality_score),file=logfile)
    print(" ",file=logfile)
    print('Number of mutations expected from the variant selection stage:', int(np.sum(simulation_input_array[:,5].astype(float))), file=logfile)
    print("Number of reads selected for variant simulation: " + str(len(var_series)), file=logfile)
    print("Total number of mutations: " + str(mutated_counter), file=logfile)    #mutated_counter
    print("Total reads in the input BAM file: " + str(bam_read_counter), file=logfile) #bam_read_counter

    print(" ")
    print("Completed variant simulation stage!")


    print_both(" ",logfile)
    print_both("Variant simulation stage runtime: " + str(round((time.time() - start_time2),3)) + " seconds", logfile)
    print_both(" ",logfile)

    logfile.close()



def variant_evaluation():


    start_time3=time.time()
    if out_prefix:
        logfile = open(fullrunpath + '/' + out_prefix + '_simulation_log.txt','a')
    else:
        logfile = open(fullrunpath + '/simulation_log.txt','a')

    print('Started variant evaluation stage:')
    print("-------------------------------")

    #~~~~~~~~~~~~~~~~~~~~
    global out_filtpileup
    # prepare reformatted input bed file for mpileup for just the filtered positions 
    out_filtpileup = []
    for p in range(len(simulation_input_array)):
        if "chr" in str(simulation_input_array[0,0]):
            fline = str(simulation_input_array[p,0]) + ":" + str(simulation_input_array[p,1]) + "-" + str(simulation_input_array[p,1])
            out_filtpileup.append(fline)
        else:
            try:
                fline = str(int(float(simulation_input_array[p,0]))) + ":" + str(simulation_input_array[p,1]) + "-" + str(simulation_input_array[p,1])
                out_filtpileup.append(fline)
            except:
                fline = str(simulation_input_array[p,0]) + ":" + str(simulation_input_array[p,1]) + "-" + str(simulation_input_array[p,1])
                out_filtpileup.append(fline)

    #########################

    #~~~~~~~~~~~~~~~~~~~~

    # prepare reformatted input bed file for mpileup for just the filtered positions 
    print("Started counting alleles in output BAM")
    out_allelecount_list = []
    for f in range(len(out_filtpileup)):
    #     out_allelecount_inter = count_alleles(coverage_MQ,coverage_BQ,out_filtpileup[f],newbamname_sorted,tmpdir)
        out_allelecount_inter = count_alleles(coverage_MQ,coverage_BQ,out_filtpileup[f],newbamname_sorted)
        out_allelecount_list.append(out_allelecount_inter)
        if f == int(len(out_filtpileup)/4):
            print("    25% complete...")
        if f==int(len(out_filtpileup)/2):
            print("    50% complete...")
        if f==int(3*len(out_filtpileup)/4):
            print("    75% complete...")
        if f==int(len(out_filtpileup)-1):
            print("    100% complete...")


    out_allelecount = np.array(out_allelecount_list)   

    global raw_simulation_array
    raw_simulation_array = np.concatenate((simulation_input_array, out_allelecount[:,:5]),1)




    #########################  

    #~~~~~~~~~~~~~~~~~~~~

    # filter out false positives (0% VAF) and polymorphic sites
    out_allelecount_filt=np.zeros((len(out_allelecount),1))
    out_allelecount_nomutation = np.zeros((len(out_allelecount),1))
    out_allelecount_polymutation = np.zeros((len(out_allelecount),1))

    for b in range(len(out_allelecount)):
        out_allelecount_filt[b]=np.size(out_allelecount[b,:4].nonzero()) == 2
        out_allelecount_nomutation[b]=np.size(out_allelecount[b,:4].nonzero()) < 2
        out_allelecount_polymutation[b]=np.size(out_allelecount[b,:4].nonzero()) > 2

    out_bool_ac_f = np.array(out_allelecount_filt,dtype=bool)
    global out_finarray
    out_finarray = np.concatenate((raw_simulation_array, out_bool_ac_f),1)



    #~~~~~~~~~~~~~~~~~~~~
    #calculate variant allele fraction
    out_allelecount_poly_filt = out_finarray[:,14:18].astype(float)
    out_allele_numeric = []
    for jj in range(len(out_allelecount_poly_filt)):
        if len(out_allelecount_poly_filt[jj].nonzero()[0]) != 1:
            out_allele_numeric.append(out_allelecount_poly_filt[jj].nonzero()[0])
        else:
            monoallele=out_allelecount_poly_filt[jj].nonzero()[0][0]
            nonexist= np.delete(np.linspace(0,3,4).astype(int),out_allelecount_poly_filt[jj].nonzero()[0][0])[0]
            out_allele_numeric.append(np.array([nonexist,monoallele]))

    out_allele_numeric_fix=np.reshape(out_allele_numeric,(np.size(out_allele_numeric,0),2))

    out_allele_alt_ref = []
    for mm in range(np.size(out_allele_numeric_fix,0)):
        out_allele_alt_ref.append(out_allelecount_poly_filt[mm,out_allele_numeric_fix[mm]])
    out_allele_alt_ref_fix=np.reshape(out_allele_alt_ref,(np.size(out_allele_alt_ref,0),2))

    out_allelefrac = np.zeros((np.size(out_allele_alt_ref_fix,0),1))
    for uu in range(np.size(out_allele_alt_ref_fix,0)):
        if out_allele_alt_ref_fix[uu,0] < out_allele_alt_ref_fix[uu,1]:
            out_allelefrac[uu] = round(out_allele_alt_ref_fix[uu,0]/int(out_finarray[uu,4]),5)
        elif out_allele_alt_ref_fix[uu,1] < out_allele_alt_ref_fix[uu,0]:
            out_allelefrac[uu] = round(out_allele_alt_ref_fix[uu,1]/int(out_finarray[uu,4]),5)


    #add calculated VAF to output array
    out_finarray_allelefrac = np.concatenate((out_finarray[:,0:5], out_allelefrac, out_finarray[:,18].reshape((-1,1)),out_finarray[:,12:14]),1)


    ######################### 


    inrange_list = []
    #check for failed mutations
    incomplete_mutation_list = []
    for ik in range(len(out_finarray_allelefrac)):
        if user_vaf_status == False:
            if not ((vaf_low/100 - 0.005) < float(out_finarray_allelefrac[ik,5]) < (vaf_high/100 - 0.005)):
                incomplete_mutation_list.append(out_finarray[ik])
            else:
                inrange_list.append(out_finarray_allelefrac[ik])
        elif user_vaf_status == True:
            if not (min(user_vaf.astype(float))-0.005 < float(out_finarray_allelefrac[ik,5]) < max(user_vaf.astype(float))+0.005):
                incomplete_mutation_list.append(out_finarray[ik])
            else:
                inrange_list.append(out_finarray_allelefrac[ik])

    incomplete_mutation_array = np.array(incomplete_mutation_list)
    inrange_array = np.array(inrange_list)

    if len(incomplete_mutation_array) != 0:
        rounded_incomplete_input_VAF = np.reshape(np.around(incomplete_mutation_array[:,3].astype(float),decimals=5),(-1,1))
        incomplete_mutation_array_rounded_0based = np.concatenate((incomplete_mutation_array[:,0].reshape(-1,1),
                                                            (incomplete_mutation_array[:,1].astype(int)-1).reshape(-1,1),
                                                            incomplete_mutation_array[:,2].reshape(-1,1),
                                                            rounded_incomplete_input_VAF,
                                                            incomplete_mutation_array[:,4:-1]),1)


        if out_prefix:
            np.savetxt(fullrunpath + '/' + out_prefix + '_failed_mutation.txt', incomplete_mutation_array_rounded_0based, fmt="%-11s",delimiter=" ",header="chromosome position position input_VAF coverage variant_coverage count_T_in count_G_in count_A_in count_C_in count_total_in monoallelic_in ref_allele alt_allele count_T_out count_G_out count_A_out count_C_out count_total_out",comments='')              
        else:
            np.savetxt(fullrunpath + '/failed_mutation.txt', incomplete_mutation_array_rounded_0based, fmt="%-11s",delimiter=" ",header="chromosome position position input_VAF coverage variant_coverage count_T_in count_G_in count_A_in count_C_in count_total_in monoallelic_in ref_allele alt_allele count_T_out count_G_out count_A_out count_C_out count_total_out",comments='')


    ######################### 

    #includes failed mutations
    simulation_output_array_total = np.reshape(out_finarray_allelefrac, (np.size(out_finarray_allelefrac,0),np.size(out_finarray_allelefrac,1)))


    #only includes in-range mutations
    simulation_output_array = np.reshape(inrange_array, (np.size(inrange_array,0),np.size(inrange_array,1)))




    rounded_output_array_input_VAF = np.reshape(np.around(simulation_output_array[:,3].astype(float),decimals=5),(-1,1))
    simulation_output_array_rounded_0based = np.concatenate((simulation_output_array[:,0].reshape(-1,1),
                                                            (simulation_output_array[:,1].astype(int)-1).reshape(-1,1),
                                                            simulation_output_array[:,2].reshape(-1,1),
                                                             rounded_output_array_input_VAF,
                                                             simulation_output_array[:,4:]),1)


    if out_prefix:
        np.savetxt(fullrunpath + '/' + out_prefix + '_simulation_output.txt', simulation_output_array_rounded_0based, fmt="%-11s", delimiter=" ",header="chromosome position position input_VAF input_coverage output_VAF output_coverage ref_allele alt_allele", comments='')
    else:    
        np.savetxt(fullrunpath + '/simulation_output.txt', simulation_output_array_rounded_0based, fmt="%-11s", delimiter=" ",header="chromosome position position input_VAF input_coverage output_VAF output_coverage ref_allele alt_allele", comments='')


    print("Variant evaluation stage metrics:", file=logfile)
    print("-------------------------------", file=logfile)
    if len(incomplete_mutation_array) !=0:
        print('WARNING: Number of failed mutations (positions with a final VAF value outside of the target VAF value range):', len(incomplete_mutation_array), file=logfile)
    print('Final number of unique spike-in positions inside the target VAF value range:', len(inrange_array), file=logfile)
    print("Final mean coverage at spike-in positions inside the target VAF value range:", np.mean(inrange_array[:,6].astype(float)), file=logfile)
    print(" ")


    #plot includes both failed and in-range mutations
    vaf_alt_distribution(simulation_output_array_total[:,5],simulation_output_array_total[:,8],logfile,"(output)")


    # remove the temporary directory
    os.system('rm ' + newbamname)
    os.system('rm -r ' + tmpdir)
    if desired_ds_cov != 0 and type(desired_ds_cov) == int:
        os.system('rm ' + originalbamname)
        os.system('rm ' + originalbamname + '.bai')

    print("Completed variant evaluation stage!")
    print(" ")


    print_both("####################",logfile)
    print_both("Simulation complete!",logfile)
    print_both("####################",logfile)
    print_both(" ",logfile)



    print_both("Variant evaluation stage runtime: " + str(round((time.time() - start_time3),3)) + " seconds", logfile)
    print_both(" ",logfile)
    print_both("Total runtime: " + str(round((time.time() - start_time),3)) + " seconds", logfile)
    print_both(" ",logfile)



    logfile.close()


