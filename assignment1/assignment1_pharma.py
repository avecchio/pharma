
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import argparse

from scipy.stats import chi2_contingency

def parse_tsv_entry(tsv_line):
    return tsv_line.strip().split("\t")

# read all phenotype entries from text file and store in a dictionary
# dictionary format is { subject_id : phenotype_response}
def read_phenotype_responses_from_file(phenotype_filepath):
    phenotype_responses = {}
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

def extract_vcf_entry_info(unparsed_info):
    info_dict = {}
    info_data = unparsed_info.split(";")
    for info_entry in info_data:
        try:
            info_key, info_value = info_entry.split("=")
            info_dict[info_key] = info_value
        except:
            pass
            #print(info_data)
            #print(info_entry)
    return info_dict

def extract_vcf_entries(vcf_filepath):
    vcf_header = None
    vcf_entries = []
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

def count_vcf_entries(variants):
    return len(variants)

def count_multiple_alternate_alleles(variants):
    ac_counter = 0
    for variant in variants:
        variant_info = variant['INFO']
        if int(variant_info['AC'].split(",")[0]) > 1:
            ac_counter += 1
        elif int(variant_info['AC']) > 1:
            ac_counter += 1
    return ac_counter

def count_passing_variants(variants):
    passing_variant_counter = 0
    for variant in variants:
        if variant['FILTER'] == 'PASS':
            passing_variant_counter += 1
    return passing_variant_counter

def count_snp_variants(variants):
    snp_entries_counter = 0
    for variant in variants:
        if variant['INFO']['VT'] == 'SNP':
            snp_entries_counter += 1
    return snp_entries_counter

def count_indel_variants(variants):
    snp_indels_counter = 0
    for variant in variants:
        if variant['INFO']['VT'] == 'INDEL':
            snp_indels_counter += 1
    return snp_indels_counter

def convert_coded_genotype_to_readable_form(genotype):
    if '|' not in genotype:
        print('error')
        print(genotype)
        return
    elif genotype == '0|0': #in ['0|0', '0/0']:
        return 'AA'
    elif '0' in genotype: #genotype_data[0] == '0' or genotype_data[1] == '0':
        return 'Aa'
    else:
        return 'aa'
    #genotype_data = genotype.split("|")
    
#import scipy.stats as stats
from scipy.stats import chisquare
import numpy as np
import math
def variant_chi_square_analysis(freq_data, observed_key, expected_key):
    df = pd.DataFrame(freq_data).transpose()
    observed_values = list(freq_data[observed_key].values())
    expected_values = list(freq_data[expected_key].values())

    if 0 in expected_values or 0.0 in expected_values:
        print(df)
        return 0, 1

    column_sums = list(df.sum(axis=0).values)
    # if any of the columns sums to zero
    # return the p value as 1
    for val in column_sums:
        if val == 0:
            return 0, 1

    #print(observed_values, expected_values)
    res = chisquare(observed_values, f_exp=expected_values)
    chi_square, p_value = res
    if math.isinf(chi_square):
        p_value = 0.1
    return chi_square, p_value #res[0], res[1]

def count_phenotype_responses_for_variant(variant, phenotype_responses):
        genotype_freq_table = {
            'cases' : {'aa': 0, 'Aa': 0, 'AA': 0},
            'controls': {'aa': 0, 'Aa': 0, 'AA': 0}
        }
        for subject in phenotype_responses:
            allele = convert_coded_genotype_to_readable_form(variant[subject])
            if phenotype_responses[subject] == '0':
                genotype_freq_table['cases'][allele] += 1
            else:
                genotype_freq_table['controls'][allele] += 1

        variant_with_genotype_freq_table = variant.copy()
        variant_with_genotype_freq_table['genotype_freq_table'] = genotype_freq_table
        return variant_with_genotype_freq_table

def record_significant_p_values(variants, id_key, p_value_key, significance, output_filepath):
    record_entries = []
    for variant in variants:
        if variant[p_value_key] < significance:
            pval = variant[p_value_key]
            id = variant[id_key]
            record_entries.append(f'{id} {pval}')

    with open(output_filepath, 'w') as f:
        f.write('\n'.join(record_entries))

def calculate_hwe(aa, Aa, AA):
    total = aa + Aa + AA
    p = (Aa + (AA*2)) / (total * 2)
    q = 1 - p

    expected_aa = (q * q) * total
    expected_Aa = 2 * (p * q) * total
    expected_AA = (p * p) * total

    genotype_freq_table = {
        'observed' : {'aa': aa, 'Aa': Aa, 'AA': AA},
        'expected': {
            'aa': expected_aa,
            'Aa': expected_Aa,
            'AA': expected_AA
        }
    }

    chi_square, p_value = variant_chi_square_analysis(genotype_freq_table, 'observed', 'expected')
    print(chi_square, p_value)

def generate_manhattan_plot(
        data, chromosome_key, p_value_key, export_path = None, significance = 6,
        colors = ['#E24E42', '#008F95'], refSNP = False
    ):

    data.to_csv('print.csv')

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

def evaluate_variant_phenotype_significance(variants, phenotype_responses):
    snp_phenotype_p_values = []
    for variant in variants:
        snp_phenotype_frequencies = count_phenotype_responses_for_variant(variant, phenotype_responses)
        chi_square, p_value = variant_chi_square_analysis(snp_phenotype_frequencies['genotype_freq_table'], 'cases', 'controls')
        print(chi_square, p_value)
        snp_phenotype_p_value = snp_phenotype_frequencies.copy()
        del snp_phenotype_p_value['genotype_freq_table']
        snp_phenotype_p_value['p_value'] = p_value
        snp_phenotype_p_values.append(snp_phenotype_p_value)

    chromosome_key = 'CHROM'
    p_value_key = 'p_value'
    statistical_significance_threshold = 5e-8
    record_significant_p_values(
        snp_phenotype_p_values,
        chromosome_key,
        p_value_key,
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


if __name__ == '__main__':

    #calculate_hwe(200, 400, 400)
    #parser = argparse.ArgumentParser(description='Genotype Analysis')
    #parser.add_argument('--phen_responses', dest='phen',
    #                help='a two column (no header) text file containing phenotype responses')
    #parser.add_argument('--vcf', dest='vcf',
    #                help='a vcf file')

    #args = parser.parse_args()
    #print(args)
    phenotypes = read_phenotype_responses_from_file("assignment1_phen.txt")
    vcf_data = extract_vcf_entries("assignment1_geno.vcf")
    evaluate_variant_phenotype_significance(vcf_data, phenotypes)
