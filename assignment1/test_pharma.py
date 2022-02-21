
from genotype_significance import *

def test_chi_square_nonzero_columns():
    data_table = {
        'cases' : {'aa': 1, 'Aa': 2, 'AA': 3},
        'controls': {'aa': 0, 'Aa': 4, 'AA': 0}
    }
    chi_square, p_value = variant_chi_square_analysis(data_table)
    assert p_value == 0.10836802322189587

def test_chi_square_zero_sum_column():
    data_table = {
        'cases' : {'aa': 0, 'Aa': 2, 'AA': 3},
        'controls': {'aa': 0, 'Aa': 4, 'AA': 0}
    }
    chi_square, p_value = variant_chi_square_analysis(data_table)
    assert p_value == 1

def test_chi_square_full_values():
    data_table = {
        'cases' : {'aa': 1, 'Aa': 2, 'AA': 3},
        'controls': {'aa': 1, 'Aa': 4, 'AA': 1}
    }
    chi_square, p_value = variant_chi_square_analysis(data_table)
    assert p_value == 0.4345982085070783

def test_genotype_conversion():
    assert convert_coded_genotype_to_readable_form('0|0') == 'AA' 
    assert convert_coded_genotype_to_readable_form('0|1') == 'Aa'
    assert convert_coded_genotype_to_readable_form('1|1') == 'aa'
    assert convert_coded_genotype_to_readable_form('2|0') == 'Aa'
    assert convert_coded_genotype_to_readable_form('0|2') == 'Aa'
    assert convert_coded_genotype_to_readable_form('1|2') == 'Aa'
    assert convert_coded_genotype_to_readable_form('2|1') == 'Aa'
    assert convert_coded_genotype_to_readable_form('2|2') == 'aa'
