
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy.stats import chi2_contingency

def clean_tsv_entry(tsv_line):
    return tsv_line.strip().split("\t")

# read all phenotype entries from text file and store in a dictionary
# dictionary format is { subject_id : phenotype_response}
def read_phenotype_responses_from_file(phenotype_filepath):
    phenotype_responses = {}
    with open(phenotype_filepath) as phenotype_file:
        for line in phenotype_file.readlines()[1:]:
            phenotype_response_entry = clean_tsv_entry(line)
            subject_id = phenotype_response_entry[0]
            phenotype_response = phenotype_response_entry[1]
            if subject_id not in phenotype_responses:
                phenotype_responses[subject_id] = phenotype_response
            else:
                print(f'Error, duplicate entry {subject_id}')
                return
    return phenotype_responses

def extract_vcf_entries(vcf_filepath):
    vcf_header = None
    vcf_entries = []
    with open(vcf_filepath) as vcf_file:
        for line in vcf_file.readlines():
            # track the header
            if line.startswith('#CHROM'):
                vcf_header = clean_tsv_entry(line[1:])
            # track all lines that are not commented
            elif line[0] != '#':
                vcf_entry = clean_tsv_entry(line)
                vcf_entries.append(vcf_entry)
    vcf_data = []
    for vcf_entry in vcf_entries:
        vcf_entry_dict = {}
        for index in range(0, len(vcf_header)):
            vcf_entry_dict[vcf_header[index]] = vcf_entry[index]
        vcf_data.append(vcf_entry_dict)
    return vcf_data

    #    print(len(vcf_entry), len(vcf_header))

                #entry_info = vcf_entry[7].split(";")
                #ac = entry_info[0].split("=")[1]
                #if ac.split(","):
                #    if int(ac.split(",")[0]) > 1:
                #        print(f'hi {ac.split(",")}')
                #        ac_counter += 1
                #elif int(ac) > 1:
                #    ac_counter += 1
                ##print(vcf_entry)
                #counter += 1
        #print(f'Vcf entries: {counter}')
        #print(f'Alternate variant counts: {ac_counter}')


#def is_reactive(subject):
#    return phenotypes[subject] == '1'

def convert_genotype(genotype):
    if genotype == '0|0':
        return 'AA'
    elif genotype == '0|1' or genotype == '1|0':
        return 'Aa'
    else:
        return 'aa'

def chi_square_analysis(data_table):
    df = pd.DataFrame(data_table).transpose()
    column_sums = list(df.sum(axis=0).values)
    # if any of the columns sums to zero
    # return the p value as 1
    for val in column_sums:
        if val == 0:
            return 1
    # calculate the p value of the contingency table
    c, p, dof, expected = chi2_contingency(df)
    return p

def count_phenotype_responses_for_variant(variant, phenotype_responses):
        genotype_freq_table = {
            'cases' : {'aa': 0, 'Aa': 0, 'AA': 0},
            'controls': {'aa': 0, 'Aa': 0, 'AA': 0}
        }
        #print(variant)
        for subject in phenotype_responses:
            allele = convert_genotype(variant[subject])
            if phenotype_responses[subject] == '1':
                genotype_freq_table['cases'][allele] += 1
            else:
                genotype_freq_table['controls'][allele] += 1

        variant_with_genotype_freq_table = variant.copy()
        variant_with_genotype_freq_table['genotype_freq_table'] = genotype_freq_table
        return variant_with_genotype_freq_table
        #return {
        #    'chromosome': variant['CHROM'],
        #    'position': variant['POS'],
        #    'id': variant['ID'],
        #    'data_table': data_table
        #}


def count_phenotype_responses_foreach_variant(variants, phenotype_responses):
    frequencies = []
    for variant in variants:
        snp_phenotype_frequencies = count_phenotype_responses_for_variant(variant, phenotype_responses)
        p_value = chi_square_analysis(snp_phenotype_frequencies['genotype_freq_table'])
        snp_phenotype_frequencies['p_value'] = p_value
        #frequencies.append(variant_data)
    return frequencies

def chi_square(frequencies):
    p_value_data = []
    for freq in frequencies:
        pval = chi_square_analysis(freq['data_table'])
        #print(pval)
        chr = freq['chromosome']
        pval_data = {
            'chromosome': f'chr{chr}',
            'position': freq['position'],
            'id': freq['id'],
            'p_value': pval
        }
        p_value_data.append(pval_data)
    return p_value_data

def FormatData(data, chromosome = 'chr', p_value = 'p_wald'):
    data['-log10(p_value)'] = -np.log10(data[p_value])
    data[chromosome] = data[chromosome].astype('category')
    data['ind'] = range(len(data))
    data_grouped = data.groupby((chromosome))
    return data, data_grouped

def GenerateManhattan(pyhattan_object, export_path = None, significance = 6, colors = ['#E24E42', '#008F95'], refSNP = False):
    data = pyhattan_object[0]
    data_grouped = pyhattan_object[1]

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
    plt.axhline(y=significance, color='gray', linestyle='-', linewidth = 0.5)
    plt.xticks(fontsize=8, rotation=60)
    plt.yticks(fontsize=8)

    if refSNP:
        for index, row in data.iterrows():
            if row['-log10(p_value)'] >= significance:
                ax.annotate(str(row[refSNP]), xy = (index, row['-log10(p_value)']))

    if export_path:
        plt.savefig(export_path)

    plt.show()


phenotypes = read_phenotype_responses_from_file("assignment1_phen.txt")
vcf_data = extract_vcf_entries("assignment1_geno.vcf")
frequencies = count_frequencies(vcf_data, phenotypes)
p_values = chi_square(frequencies)

chr_max_mins = {}
for p_val in p_values:
    if p_val['chromosome'] not in chr_max_mins:
        chr_max_mins[p_val['chromosome']] = []
    chr_max_mins[p_val['chromosome']].append(int(p_val['position']))

chr_diffs = {}
for chr in chr_max_mins:
    positions = chr_max_mins[chr]
    max_pos = max(positions)
    min_pos = min(positions)
    chr_diffs[chr] = min_pos - 1

lines = []
for p_val in p_values:
    if p_val['p_value'] < 5e-8:
        pval = p_val['p_value']
        id = p_val['id']
        lines.append(f'{id} {pval}')

with open('significant_variants.txt', 'w') as f:
    f.write('\n'.join(lines))

formatted_manhattan_data = []
for v in p_values:
    formatted_manhattan_data.append({
        "#CHROM": v['chromosome'],
        "POS": int(v['position']), # - chr_diffs[v['chromosome']],
        "P": v['p_value'],
        "ID": v['id']
    })

data = pd.DataFrame(formatted_manhattan_data)
formatted_data = FormatData(data, chromosome='#CHROM', p_value='P')
GenerateManhattan(formatted_data, export_path="manhattan.png")
