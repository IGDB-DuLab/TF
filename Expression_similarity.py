# -*- coding: utf-8 -*-
"""
Created on Mon May  3 20:58:07 2021

@author: ZZG
"""



import scipy.stats as stats

def hypergeom_test(totoal_num_in_population, total_num_with_condition_in_pupulation, num_of_subset, num_with_condition_in_subset):
    
    less_pvalue = str(stats.hypergeom.cdf(int(num_with_condition_in_subset) ,int(totoal_num_in_population),int(total_num_with_condition_in_pupulation),int(num_of_subset)))
    great_pvalue = str(stats.hypergeom.sf(int(num_with_condition_in_subset) - 1,int(totoal_num_in_population),int(total_num_with_condition_in_pupulation),int(num_of_subset)))
    
    return less_pvalue, great_pvalue

#traced_cell_list : traced cell number
#TF1_exp_cell_list : expression cell list of TF1
#TF2_exp_cell_list : expression cell list of TF2

co_exp_cell_num = len(set(TF1_exp_cell_list) & set(TF2_exp_cell_list))


similairy = ( float(co_exp_cell_num) * traced_cell_list ) / ( len(TF1_exp_cell_list) * len(TF2_exp_cell_list) )
less_pvalue, great_pvalue = hypergeom_test(traced_cell_list,  len(TF1_exp_cell_list), len(TF2_exp_cell_list), co_exp_cell_num)

# output
# similairy : similairy of two TF expression pattern
# great_pvalue : pvalue