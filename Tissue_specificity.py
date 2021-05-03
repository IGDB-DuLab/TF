# -*- coding: utf-8 -*-
"""
Created on Mon May  3 18:58:06 2021

@author: ZZG
"""



import scipy.stats as stats

def hypergeom_test(totoal_num_in_population, total_num_with_condition_in_pupulation, num_of_subset, num_with_condition_in_subset):
    
    less_pvalue = str(stats.hypergeom.cdf(int(num_with_condition_in_subset) ,int(totoal_num_in_population),int(total_num_with_condition_in_pupulation),int(num_of_subset)))
    great_pvalue = str(stats.hypergeom.sf(int(num_with_condition_in_subset) - 1,int(totoal_num_in_population),int(total_num_with_condition_in_pupulation),int(num_of_subset)))
    
    return less_pvalue, great_pvalue


#all_traced_cell_number :  all traced cell number
#expression_cell_number :  detected expression cell number
#expected_tissue_cell_number : tissue composition of exprected cell based on stage of detected expression cell 
#observed_expression_tissue_cell_number: tissue composition of detected expression cell 

enrichement = observed_expression_tissue_cell_number / expected_tissue_cell_number

less_pvalue, great_pvalue = hypergeom_test(all_traced_cell_number, expected_tissue_cell_number, expression_cell_number, observed_expression_tissue_cell_number)            


# output
# enrichement : enrichment score of specific tissue
# great_pvalue : pvalue

