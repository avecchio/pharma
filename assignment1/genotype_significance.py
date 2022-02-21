from scipy.stats import chisquare, chi2_contingency

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

import argparse
import sys

###################################
#File Processing
###################################

# helper function
def parse_tsv_entry(tsv_line):
    return tsv_line.strip().split("\t")

# read all phenotype entries from text file and store in a dictionary
# dictionary format is { subject_id : phenotype_response}

# input: list of variants. Each variant is in the form of a dictionary
#        where the column matches the column in the VCF file
# output: (int) count of variants (length of array)

def read_phenotype_responses_from_file(phenotype_filepath):
    phenotype_responses = {}
    try:
        with open(phenotype_filepath) as phenotype_file:
            for line in phenotype_file.readlines()[1:]:
                phenotype_response_entry = parse_tsv_entry(line)
                subject_id = phenotype_response_entry[0]
                phenotype_response = phenotype_response_entry[1]
                if subject_id not in phenotype_responses:
                    phenotype_responses[subject_id] = phenotype_response
                else:
                    print(f'Error, duplicate entry {subject_id}')
                    return
        return phenotype_responses
    except:
        print(f'Error reading or locating file {phenotype_filepath}')
        sys.exit()

# input: list of variants. Each variant is in the form of a dictionary
#        where the column matches the column in the VCF file
# output: (int) count of variants (length of array)

def extract_vcf_entry_info(unparsed_info):
    info_dict = {}
    info_data = unparsed_info.split(";")
    for info_entry in info_data:
        try:
            info_key, info_value = info_entry.split("=")
            info_dict[info_key] = info_value
        except:
            pass
    return info_dict
# input: list of variants. Each variant is in the form of a dictionary
#        where the column matches the column in the VCF file
# output: (int) count of variants (length of array)

def extract_vcf_entries(vcf_filepath):
    vcf_header = None
    vcf_entries = []
    try:
        with open(vcf_filepath) as vcf_file:
            for line in vcf_file.readlines():
                # track the header
                if line.startswith('#CHROM'):
                    vcf_header = parse_tsv_entry(line[1:])
                # track all lines that are not commented
                elif line[0] != '#':
                    vcf_entry = parse_tsv_entry(line)
                    vcf_entries.append(vcf_entry)
        vcf_data = []
        for vcf_entry in vcf_entries:
            vcf_entry_dict = {}
            for index in range(0, len(vcf_header)):
                if vcf_header[index] == 'INFO':
                    entry_info_dict = extract_vcf_entry_info(vcf_entry[index])
                    vcf_entry_dict[vcf_header[index]] = entry_info_dict
                else:
                    vcf_entry_dict[vcf_header[index]] = vcf_entry[index]
            vcf_data.append(vcf_entry_dict)
        return vcf_data
    except:
        print(f'Error reading or locating file {vcf_filepath}')
        sys.exit()


###################################
#Statistical methods for question 1
###################################

# input: list of variants. Each variant is in the form of a dictionary
#        where the column matches the column in the VCF file
# output: (int) count of variants (length of array)
def count_vcf_entries(variants):
    return len(variants)

# input: list of variants. Each variant is in the form of a dictionary
#        where the column matches the column in the VCF file
# output: (int) count of variants have more then one alternate allele
def count_multiple_alternate_alleles(variants):
    ac_counter = 0
    for variant in variants:
        # alternate allele counts exists within INFO column
        variant_info = variant['INFO']
        # two numbers might exist for the alternate allele count
        # in the format of #,#
        # check to see if the first number (the alternate allele count) is greater then one 
        if int(variant_info['AC'].split(",")[0]) > 1:
            ac_counter += 1
        # only a single number exists for the AC value
        # check if alternate allele count is greater then one
        elif int(variant_info['AC']) > 1:
            ac_counter += 1
    return ac_counter

# input: list of variants. Each variant is in the form of a dictionary
#        where the column matches the column in the VCF file
# output: (int) count of variants that passed all filters
def count_passing_variants(variants):
    passing_variant_counter = 0
    for variant in variants:
        # if the variant passed all filter (recorded in the vcf file)
        # increase the counter by one
        if variant['FILTER'] == 'PASS':
            passing_variant_counter += 1
    return passing_variant_counter

# input: list of variants. Each variant is in the form of a dictionary
#        where the column matches the column in the VCF file
# output: (int) count of variants that are single nucleotide polymorphisms
def count_snp_variants(variants):
    snp_entries_counter = 0
    for variant in variants:
        # if the variant is a SNP, increase the counter by one
        if variant['INFO']['VT'] == 'SNP':
            snp_entries_counter += 1
    return snp_entries_counter

# input: list of variants. Each variant is in the form of a dictionary
#        where the column matches the column in the VCF file
# output: (int) count of variants that are indels
def count_indel_variants(variants):
    snp_indels_counter = 0
    for variant in variants:
        # if the variant is an indel, increase the counter by one
        if variant['INFO']['VT'] == 'INDEL':
            snp_indels_counter += 1
    return snp_indels_counter

def count_hwe_significant_variants(snp_phenotype_p_values, statistical_significance_threshold):
    hwe_sig = 0
    for p_value_stats in snp_phenotype_p_values:
        if p_value_stats['hwe_p_value'] <= statistical_significance_threshold:
            hwe_sig += 1

    return hwe_sig

#############################################
# Core processing or statistical calculations
#############################################


# input: genotype in the form of .|.
# output: coded genotype [aa, Aa, AA]
# convert genotype in vcf into readable form
# the following mapping pattern is used
# 0|0 => AA
# 1|1, 2|2, 3|3, etc. => aa
# 1|0, 0|2, 1|2, 3|1, etc. => Aa
def convert_coded_genotype_to_readable_form(genotype):
    # split the genotype into the following array format [ . , .]
    genotype_frac = genotype.split("|")
    # if all the values in the genotype are zeros (dominant)
    # return the code for homozygous dominant: AA
    if genotype_frac.count('0') == len(genotype_frac):
        return 'AA'
    elif genotype_frac.count(genotype_frac[0]) == len(genotype_frac):
        return 'aa'
    # all other 
    else:
        return 'Aa'

# input: list of variants. Each variant is in the form of a dictionary
#        where the column matches the column in the VCF file
# output: (dictionary) the contingency table of allele counts grouped by case or control
# constructs a contingency table for the counts of each type of allele
# that is determined by convert_coded_genotype_to_readable_form
# the allele counts are grouped based on 'case' or 'control' which
# is a flag ('1') set in the phenotype responses for each subject
def count_phenotype_responses_for_variant(variant, phenotype_responses):
        genotype_freq_table = {
            'cases' : {'aa': 0, 'Aa': 0, 'AA': 0},
            'controls': {'aa': 0, 'Aa': 0, 'AA': 0}
        }
        for subject in phenotype_responses:
            allele = convert_coded_genotype_to_readable_form(variant[subject])
            if phenotype_responses[subject] == '1':
                genotype_freq_table['cases'][allele] += 1
            else:
                genotype_freq_table['controls'][allele] += 1

        variant_with_genotype_freq_table = variant.copy()
        variant_with_genotype_freq_table['genotype_freq_table'] = genotype_freq_table
        return variant_with_genotype_freq_table


# input: (dictionary) the contingency table of allele counts grouped by case or control
# output: chi square (float) and p_value (float) statistic for significance of variant
#         if any column (aa, Aa, AA) sums to 0 in the contingency table, the resulting values
#         will be 1 (chi square) and 0 (p value) respectively
def variant_chi_square_analysis(genome_freq_table):

    df = pd.DataFrame(genome_freq_table).transpose()

    column_sums = list(df.sum(axis=0).values)
    # if any of the columns sums to zero
    # return the p value as 1
    for val in column_sums:
        if val == 0:
            return 0, 1

    # using chi2_contingency to determine if there exists a
    # statistical difference between any of the groups (aa, Aa, AA)
    chi2, p, dof, ex = chi2_contingency(df, correction=False)
    return chi2, p

# input: frequeny of alleles: aa (int), Aa (int) AA (int)
#        where the column matches the column in the VCF file
# output: the hwe chi square (int) and p_value (float) statistics
#         calculate hwe chi square and p value statistic for each variant
#         If any column (aa, Aa, AA) sums to 0 in the contingency table, the resulting values
#         will be 1 (chi square) and 0 (p value) respectively
def calculate_hwe(aa, Aa, AA):

    # calculate the expected values
    total = aa + Aa + AA
    p = (Aa + (AA*2)) / (total * 2)
    q = 1 - p

    expected_aa = (q * q) * total
    expected_Aa = 2 * (p * q) * total
    expected_AA = (p * p) * total

    observed_values = [aa, Aa, AA]
    expected_values = [
        expected_aa,
        expected_Aa,
        expected_AA
    ]

    # if any of the columns (aa, Aa, or AA) sum to zero,
    # chi square is set to 0 and the p value is set to 1 as the actual
    # chi square equation cannot properly calcuate the real value
    for column_index in range(0, len(observed_values)):
        if observed_values[column_index] + expected_values[column_index] == 0:
            return 0, 1

    # using scipy chisquare to determine if there is a statistical difference
    # between all observed vs expected values
    res = chisquare(observed_values, f_exp=expected_values)
    chi_square, p_value = res

    return chi_square, p_value


###################################
#Output generators
###################################

def record_significant_p_values(variants, id_key, p_value_key, hwe_p_value_key, significance, output_filepath):
    record_entries = []
    record_entries.append(f'VAR_ID P_VALUE HWE_P_VALUE')
    for variant in variants:
        if variant[p_value_key] < significance:
            p_val = variant[p_value_key]
            hwe_p_value = variant[hwe_p_value_key]
            id = variant[id_key]
            record_entries.append(f'{id} {p_val} {hwe_p_value}')

    with open(output_filepath, 'w') as f:
        f.write('\n'.join(record_entries))

def generate_manhattan_plot(
        data, chromosome_key, p_value_key, export_path = None, significance = 5e-8,
        colors = ['#E24E42', '#008F95'], refSNP = False
    ):

    data['-log10(p_value)'] = -np.log10(data[p_value_key])
    data[chromosome_key] = data[chromosome_key].astype('category')
    data['ind'] = range(len(data))
    data_grouped = data.groupby((chromosome_key))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(data_grouped):
        group.plot(kind='scatter', x='ind', y='-log10(p_value)', color=colors[num % len(colors)], ax=ax, s= 10000/len(data))
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2))

    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, len(data)])
    ax.set_ylim([0, data['-log10(p_value)'].max() + 1])
    ax.set_xlabel('Chromosome')
    plt.axhline(y=significance, color='black', linestyle='-', linewidth = 1)
    plt.xticks(fontsize=8, rotation=60)
    plt.yticks(fontsize=8)

    if refSNP:
        for index, row in data.iterrows():
            if row['-log10(p_value)'] >= significance:
                ax.annotate(str(row[refSNP]), xy = (index, row['-log10(p_value)']))

    if export_path:
        plt.savefig(export_path)

    plt.show()

# Input:
# variants:
# phenotype_responses:
# This is considered to be the main method
def evaluate_variant_phenotype_significance(variants, phenotype_responses):

    # calculate statistics for Question 1
    count_multiple_alternate_alleles(variants)
    count_passing_variants(variants)
    count_snp_variants(variants)
    count_indel_variants(variants)

    
    snp_phenotype_p_values = []
    
    # calculate the 
    for variant in variants:
        snp_phenotype_frequencies = count_phenotype_responses_for_variant(variant, phenotype_responses)
        chi_square, p_value = variant_chi_square_analysis(snp_phenotype_frequencies['genotype_freq_table'])

        gft = snp_phenotype_frequencies['genotype_freq_table']['controls']

        hwe_chi_square, hwe_p_value = calculate_hwe(gft['aa'], gft['Aa'], gft['AA'])

        snp_phenotype_p_value = snp_phenotype_frequencies.copy()
        del snp_phenotype_p_value['genotype_freq_table']
        snp_phenotype_p_value['p_value'] = p_value
        snp_phenotype_p_value['hwe_p_value'] = hwe_p_value
        snp_phenotype_p_values.append(snp_phenotype_p_value)

    chromosome_key = 'CHROM'
    p_value_key = 'p_value'
    hwe_p_value_key = "hwe_p_value"

    statistical_significance_threshold = 5e-8

    hwe_significant_variants_count = count_hwe_significant_variants(snp_phenotype_p_value, statistical_significance_threshold)
    print(f'hwe significant variant count: {hwe_significant_variants_count}')

    record_significant_p_values(
        snp_phenotype_p_values,
        chromosome_key,
        p_value_key,
        hwe_p_value_key,
        statistical_significance_threshold,
        "significant_variants.txt"
    )
    
    generate_manhattan_plot(
        pd.DataFrame(snp_phenotype_p_values),
        chromosome_key=chromosome_key,
        p_value_key=p_value_key,
        significance= -math.log10(statistical_significance_threshold),
        export_path="genotype_significance.png"
    )


###################################
# MAIN
###################################

if __name__ == '__main__':
    # setup the arg parser
    parser = argparse.ArgumentParser(description='GWAS Analysis')

    # add two required arguments. The file paths for the phenotype responses and VCF, respectively
    parser.add_argument('--phen_responses', dest='phen', required=True,
                    help='a two column (no header) text file containing phenotype responses per subject')
    parser.add_argument('--vcf', dest='vcf', required=True,
                    help='The VCF file that contains genotype data for variants')
    args = parser.parse_args()

    # read the phenotype response file and extract all entries
    phenotypes = read_phenotype_responses_from_file(args.phen)
    # read the vcf file and extract all entries
    vcf_data = extract_vcf_entries(args.vcf)
    # calculate all of the statistics and determine statistically significant variants
    evaluate_variant_phenotype_significance(vcf_data, phenotypes)
