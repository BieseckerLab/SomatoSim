#!/usr/bin/env python3

global np
import numpy as np
global matplotlib
import matplotlib
global plt
import matplotlib.pyplot as plt
global figure
from matplotlib.pyplot import figure
global os
import os


def print_both(input_string,savefile):
    print(input_string,file=savefile)
    print(input_string)

def convert_BED(original_BED,temporary_directory):
    #create a BED file with the opposite chromosome notation (ex: "chr5" instead of "5" or vice versa)
    if "chr" in str(original_BED[0,0]):
        new_chrom = []
        new_position_1 = []
        new_position_2 = []
        for q in range(len(original_BED)):
            new_chrom.append(original_BED[q,0][3:])
            new_position_1.append(int(float(original_BED[q,1])))
            new_position_2.append(int(float(original_BED[q,2])))
        new_BED = np.concatenate((np.reshape(np.array(new_chrom),(-1,1)),np.reshape(np.array(new_position_1),(-1,1)),np.reshape(np.array(new_position_2),(-1,1))),1)
        np.savetxt(temporary_directory + '/new_BED.bed', new_BED, fmt='%s',delimiter='\t')
    else:
        new_chrom = []
        new_position_1 = []
        new_position_2 = []
        for q in range(len(original_BED)):
            try:
                new_chrom.append("chr" + str(int(float(original_BED[q,0]))))
            except:
                new_chrom.append("chr" + str((original_BED[q,0])))
            new_position_1.append(int(float(original_BED[q,1])))
            new_position_2.append(int(float(original_BED[q,2])))
        new_BED = np.concatenate((np.reshape(np.array(new_chrom),(-1,1)),np.reshape(np.array(new_position_1),(-1,1)),np.reshape(np.array(new_position_2),(-1,1))),1)
        np.savetxt(temporary_directory + '/new_BED.bed', new_BED, fmt='%s',delimiter='\t')
        
    return new_BED


def vaf_alt_distribution(vaf_values,alt_values,save_file,description_string=''):
    ######## PLOT VAF DISTRIBUTION ############################
    np.set_printoptions(suppress=True)
    figure(num=2, figsize=(8, 6), dpi=300, facecolor='w', edgecolor='k')

    histo_vaf = vaf_values.astype(float) ##############
    histo_start = 0

    histo_end = round((np.max(histo_vaf) + 0.02),2)
    histo_num_bins = int((histo_end * 100) + 1)
    VAF_bin_equal=np.linspace(histo_start,histo_end,histo_num_bins)

    histo = plt.hist(histo_vaf,bins=VAF_bin_equal)
    bincount = np.reshape(histo[0], (-1,1))
    binval = np.reshape(histo[1][:-1],(-1,1))
    NumPosPerVAF = np.concatenate((binval,bincount),1)
    plt.title('VAF distribution')
    plt.xlabel('VAF')
    plt.ylabel('Count')
    plt.close()
    
    print(" ", file=save_file)
    print("Vaiant allele fraction distribution " + description_string + ":", file=save_file)  ##############
    print(NumPosPerVAF, file=save_file)
    print(" ", file=save_file)
    
    ######## PLOT ALTERNATE ALLELE DISTRIBUTION ##############
    count_T = 0
    count_G = 0
    count_A = 0
    count_C = 0


    for c in range(len(alt_values)): ##############
        if alt_values[c] == 'T': ##############
            count_T += 1
        if alt_values[c] == 'G': ##############
            count_G +=1
        if alt_values[c] == 'A': ##############
            count_A +=1
        if alt_values[c] == 'C': ##############
            count_C +=1

    count_all_alleles=[count_T,count_G,count_A,count_C]

    allele_labels = ['T', 'G', 'A', 'C']
    allele_labels_numeric = np.arange(len(allele_labels))

    # Create plot
    figure(num=1, figsize=(8, 6), dpi=300, facecolor='w', edgecolor='k')
    plt.bar(allele_labels_numeric, count_all_alleles)

    # Add title and axis names
    plt.title('Alternate allele distribution')
    plt.xlabel('Alternate allele')
    plt.ylabel('Count')

    # Create names
    plt.xticks(allele_labels_numeric, allele_labels)

    plt.close()


    allele_count_array = np.concatenate((np.reshape(np.array(allele_labels),(-1,1)),np.reshape(np.array(count_all_alleles),(-1,1))),1) ##############

    print("Alternate allele distribution " + description_string + ":", file=save_file) ##############
    print(allele_count_array, file=save_file) ##############
    print(" ", file=save_file)
    ###############################################




def count_alleles(mapping_qual,base_qual,position_check,bam_check):
    
    pileup_sequence_read = os.popen('samtools mpileup -x -q ' + str(mapping_qual) + ' -Q ' + str(base_qual) + ' -r '+ position_check + ' ' + bam_check + ' 2> /dev/null | cut -f 5 | tr "[a-z]" "[A-Z]"').read()[:-1]

    pileup_sequence_value = np.array(list(pileup_sequence_read))

    caret = np.where(pileup_sequence_value=='^')[0]
    if len(caret) != 0:
        caret_plus_one = np.append(caret,caret+1)
        cleaned_seq = np.delete(pileup_sequence_value,caret_plus_one)
        cleaned_string = ''.join(cleaned_seq)
    else:
        cleaned_string = ''.join(pileup_sequence_value)

    cleaned_array = np.array(list(cleaned_string))

    ins_index = np.where(cleaned_array=='+')[0]

    if len(ins_index) != 0:

        all_ins_delete_index = []
        for iii in range(len(ins_index)):
            ins_value = int(cleaned_array[ins_index[iii]+1])
            ins_delete_index = np.linspace(ins_index[iii]+2,ins_index[iii]+ins_value+1,ins_value).astype(int)
            all_ins_delete_index.append(ins_delete_index)

        all_ins_delete_index_reshape = np.concatenate(all_ins_delete_index, axis=0)
        remove_ins = np.delete(cleaned_array,all_ins_delete_index_reshape)
        remove_ins_string = ''.join(remove_ins)

    else:
        remove_ins = cleaned_array
        remove_ins_string = cleaned_string


    del_index = np.where(remove_ins=='-')[0]

    if len(del_index) != 0:
        all_del_delete_index = []
        for iii in range(len(del_index)):
            del_value = int(remove_ins[del_index[iii]+1])
            del_delete_index = np.linspace(del_index[iii]+2,del_index[iii]+del_value+1,del_value).astype(int)
            all_del_delete_index.append(del_delete_index)


        all_del_delete_index_reshape = np.concatenate(all_del_delete_index, axis=0)

        remove_del = np.delete(remove_ins,all_del_delete_index_reshape)
        remove_del_string = ''.join(remove_del)


    else:
        remove_del = remove_ins
        remove_del_string = remove_ins_string


    # [T G A C total del ins]
    allelecount_fn = [0,0,0,0,0,0,0,0]
    for ch in range(len(remove_del_string)):
        if remove_del_string[ch] == 'T':
            allelecount_fn[0] +=1
        if remove_del_string[ch] == 'G':
            allelecount_fn[1] +=1
        if remove_del_string[ch] == 'A':
            allelecount_fn[2] +=1
        if remove_del_string[ch] == 'C':
            allelecount_fn[3] +=1
        if remove_del_string[ch] == '*':
            allelecount_fn[5] +=1
        if remove_del_string[ch] == '+':
            allelecount_fn[6] +=1
        if remove_del_string[ch] == '-':
            allelecount_fn[7] +=1


    allelecount_fn[4] = np.sum(np.append(allelecount_fn[:4],allelecount_fn[5])) 
    
    return allelecount_fn


def select_reads(read_bank_pos_fn,
                 read_bank_neg_fn,
                 var_depth_fn,
                 all_reads_fn,
                 haplo_status_fn,
                 variant_depth_counter_fn,
                 simulation_input_array_fn,
                 j_fn,
                 var_series_intermediate2_fn):



    fraction_strand_pos_fn = len(read_bank_pos_fn)/(len(read_bank_pos_fn) + len(read_bank_neg_fn))
    balanced_pos = round(fraction_strand_pos_fn * var_depth_fn)

    fraction_strand_neg_fn = len(read_bank_neg_fn)/(len(read_bank_pos_fn) + len(read_bank_neg_fn))
    balanced_neg = round(fraction_strand_neg_fn * var_depth_fn)

    num_pos = 0
    num_neg= 0

    succ_pos = 0
    succ_neg = 0
    current_reads = []
    while  num_pos < len(read_bank_pos_fn) or num_neg< len(read_bank_neg_fn):
        while num_pos < len(read_bank_pos_fn):
            if haplo_status_fn == True:
                if read_bank_pos_fn[num_pos] not in current_reads and succ_pos < balanced_pos:
                    all_reads_fn.append(read_bank_pos_fn[num_pos])
                    current_reads.append(read_bank_pos_fn[num_pos])
                    variant_depth_counter_fn += 1
                    var_series_intermediate1 = [read_bank_pos_fn[num_pos].query_name,read_bank_pos_fn[num_pos],j_fn,read_bank_pos_fn[num_pos].is_reverse]
                    var_series_intermediate2_fn.append(var_series_intermediate1) # this contains the query name, the actual query, and the position index (j_fn)
                    if len(var_series_intermediate2_fn) == int(np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    10% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(2.5*np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    25% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(5*np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    50% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(7.5*np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    75% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(9*np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    90% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(np.sum(simulation_input_array_fn[:,5].astype(float))):
                        print("    100% of expected reads were selected...")
                    num_pos += 1
                    succ_pos +=1
                    break
                else:
                    num_pos += 1
                    continue
            else:
                if read_bank_pos_fn[num_pos] not in current_reads and read_bank_pos_fn[num_pos] not in all_reads_fn and succ_pos < balanced_pos:
                    all_reads_fn.append(read_bank_pos_fn[num_pos])
                    current_reads.append(read_bank_pos_fn[num_pos])
                    variant_depth_counter_fn += 1
                    var_series_intermediate1 = [read_bank_pos_fn[num_pos].query_name,read_bank_pos_fn[num_pos],j_fn,read_bank_pos_fn[num_pos].is_reverse]
                    var_series_intermediate2_fn.append(var_series_intermediate1) # this contains the query name, the actual query, and the position index (j_fn)
                    if len(var_series_intermediate2_fn) == int(np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    10% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(2.5*np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    25% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(5*np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    50% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(7.5*np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    75% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(9*np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    90% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(np.sum(simulation_input_array_fn[:,5].astype(float))):
                        print("    100% of expected reads were selected...")
                    num_pos += 1
                    succ_pos +=1
                    break
                else:
                    num_pos += 1
                    continue
        if variant_depth_counter_fn == var_depth_fn:
            break
        while num_neg< len(read_bank_neg_fn):
            if haplo_status_fn == True:
                if read_bank_neg_fn[num_neg] not in current_reads and succ_neg < balanced_neg:
                    all_reads_fn.append(read_bank_neg_fn[num_neg])
                    current_reads.append(read_bank_neg_fn[num_neg])
                    variant_depth_counter_fn += 1
                    var_series_intermediate1 = [read_bank_neg_fn[num_neg].query_name,read_bank_neg_fn[num_neg],j_fn,read_bank_neg_fn[num_neg].is_reverse]
                    var_series_intermediate2_fn.append(var_series_intermediate1) # this contains the query name, the actual query, and the position index (j_fn)
                    if len(var_series_intermediate2_fn) == int(np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    10% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(2.5*np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    25% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(5*np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    50% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(7.5*np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    75% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(9*np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    90% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(np.sum(simulation_input_array_fn[:,5].astype(float))):
                        print("    100% of expected reads were selected...")
                    num_neg+= 1
                    succ_neg+=1
                    break
                else:
                    num_neg+= 1
                    continue
            else:
                if read_bank_neg_fn[num_neg] not in current_reads and read_bank_neg_fn[num_neg] not in all_reads_fn and succ_neg < balanced_neg:
                    all_reads_fn.append(read_bank_neg_fn[num_neg])
                    current_reads.append(read_bank_neg_fn[num_neg])
                    variant_depth_counter_fn += 1
                    var_series_intermediate1 = [read_bank_neg_fn[num_neg].query_name,read_bank_neg_fn[num_neg],j_fn,read_bank_neg_fn[num_neg].is_reverse]
                    var_series_intermediate2_fn.append(var_series_intermediate1) # this contains the query name, the actual query, and the position index (j_fn)
                    if len(var_series_intermediate2_fn) == int(np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    10% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(2.5*np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    25% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(5*np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    50% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(7.5*np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    75% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(9*np.sum(simulation_input_array_fn[:,5].astype(float)/10)):
                        print("    90% of expected reads were selected...")
                    if len(var_series_intermediate2_fn) == int(np.sum(simulation_input_array_fn[:,5].astype(float))):
                        print("    100% of expected reads were selected...")
                    num_neg+= 1
                    succ_neg+=1
                    break
                else:
                    num_neg+= 1
                    continue                
        if variant_depth_counter_fn == var_depth_fn:
            break

def check_BED_columns(oriBED_fn,logfile_fn):
    #if the original BED file has a third column
    try:
        oriBED_fn[:,3]
        #see if it is a number between 0 and 1 in the third column
        try:
            float(oriBED_fn[0,3])
            #if yes, store those values in user_vaf
            if 0 < float(oriBED_fn[0,3]) < 1:
                user_vaf_fn = np.reshape(oriBED_fn[:,3],(-1,1))
                user_vaf_status_fn = True
                try:
                    #see if there's a 4th column in the BED file
                    oriBED_fn[:,4]
                    #see if that 4th column has alternate alleles, and if it does, store them is user_alt
                    if (oriBED_fn[0,4] == 'T') or (oriBED_fn[0,4] == 'G') or (oriBED_fn[0,4] == 'A') or (oriBED_fn[0,4] == 'C'):
                        user_alt_fn = np.reshape(oriBED_fn[:,4],(-1,1))
                        user_alt_status_fn = True
                        print_both('BED file input given with VAF and alternate allele',logfile_fn)
                        print_both(" ",logfile_fn)
                    #otherwise if the 4th column isn't alternate alleles, then just store the VAF values
                    else:
                        print_both('BED file input given with VAF',logfile_fn)
                        print_both(" ",logfile_fn)
                        user_alt_status_fn = False
                #otherwise if there is not 4th column, then just store the VAF values
                except:
                    print_both('BED file input given with VAF',logfile_fn)
                    print_both(" ",logfile_fn)
                    user_alt_status_fn = False
            else:
                print_both('BED file input contains no pre-defined VAF or alternate allele',logfile_fn)
                print_both(" ",logfile_fn)
                user_vaf_status_fn = False
                user_alt_status_fn = False
        #if the third column isn't a number, then see if it is an alternate allele
        except:
            #if yes, store those values in user_alt
            if (oriBED_fn[0,3] == 'T') or (oriBED_fn[0,3] == 'G') or (oriBED_fn[0,3] == 'A') or (oriBED_fn[0,3] == 'C'):
                user_alt_fn = np.reshape(oriBED_fn[:,3],(-1,1))
                user_alt_status_fn = True
                try:
                    #see if there's a 4th column in the BED file
                    float(oriBED_fn[0,4])
                    #if the 4th column is a number between 0 and 1, save that column as the VAF
                    if 0 < float(oriBED_fn[0,4]) < 1:
                        user_vaf_fn = np.reshape(oriBED_fn[:,4],(-1,1))
                        user_vaf_status_fn = True
                        print_both('BED file input given with alternate allele and VAF',logfile_fn)
                        print_both(" ",logfile_fn)
                    #otherwise if the 4th column is not a number between 0 and 1, then just store the alt allele value
                    else:
                        print_both('BED file input given with alternate allele',logfile_fn)
                        print_both(" ",logfile_fn)
                        user_vaf_status_fn = False
                #otherwise if the 4th column isn't there (or isn't a float), then just store the alternate allele value
                except:
                    print_both('BED file input given with alternate allele',logfile_fn)
                    print_both(" ",logfile_fn)   
                    user_vaf_status_fn = False
            #if the third column isn't a number and it isn't an alternate allele, then go to random mode
            else:
                print_both('BED file input contains no pre-defined VAF or alternate allele',logfile_fn)
                print_both(" ",logfile_fn)
                user_vaf_status_fn = False
                user_alt_status_fn = False            
    except:
        print_both('BED file input contains no pre-defined VAF or alternate allele',logfile_fn)
        print_both(" ",logfile_fn)
        user_vaf_status_fn = False
        user_alt_status_fn = False
        
    if user_vaf_status_fn == False:
        user_vaf_fn = np.array([])
    if user_alt_status_fn == False:
        user_alt_fn = np.array([])
    
    return user_vaf_fn,user_alt_fn,user_vaf_status_fn,user_alt_status_fn


