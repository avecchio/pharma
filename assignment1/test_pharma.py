
from assignment1_pharma import *

def test_chi_square_nonzero_columns():
    data_table = {
        'cases' : {'aa': 0, 'Aa': 2, 'AA': 3},
        'controls': {'aa': 0, 'Aa': 4, 'AA': 0}
    }
    p_value = chi_square_analysis(data_table)
    print(p_value == 1)

def test_chi_square_zero_sum_column():
    data_table = {
        'cases' : {'aa': 1, 'Aa': 2, 'AA': 3},
        'controls': {'aa': 0, 'Aa': 4, 'AA': 0}
    }
    p_value = chi_square_analysis(data_table)
    print(p_value)

test_chi_square_nonzero_columns()
test_chi_square_zero_sum_column()