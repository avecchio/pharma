
from assignment1_pharma import *

def test_chi_square_nonzero_columns():
    data_table = {
        'cases' : {'aa': 1, 'Aa': 2, 'AA': 3},
        'controls': {'aa': 0, 'Aa': 4, 'AA': 0}
    }
    chi_square, p_value = variant_chi_square_analysis(data_table, 'cases', 'controls')
    print(p_value)

def test_chi_square_zero_sum_column():
    data_table = {
        'cases' : {'aa': 0, 'Aa': 2, 'AA': 3},
        'controls': {'aa': 0, 'Aa': 4, 'AA': 0}
    }
    chi_square, p_value = variant_chi_square_analysis(data_table, 'cases', 'controls')
    print(p_value)

def test_chi_square_full_values():
    data_table = {
        'cases' : {'aa': 1, 'Aa': 2, 'AA': 3},
        'controls': {'aa': 1, 'Aa': 4, 'AA': 1}
    }
    chi_square, p_value = variant_chi_square_analysis(data_table, 'cases', 'controls')
    print(p_value)

#test_chi_square_nonzero_columns()
#test_chi_square_zero_sum_column()
#test_chi_square_full_values()
print(convert_coded_genotype_to_readable_form('0|0'))
print(convert_coded_genotype_to_readable_form('0|1'))
print(convert_coded_genotype_to_readable_form('1|1'))
print(convert_coded_genotype_to_readable_form('2|0'))
print(convert_coded_genotype_to_readable_form('0|2'))
print(convert_coded_genotype_to_readable_form('1|2'))
print(convert_coded_genotype_to_readable_form('2|1'))