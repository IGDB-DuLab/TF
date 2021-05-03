# -*- coding: utf-8 -*-
"""
Created on Mon May  3 19:08:36 2021

@author: ZZG
"""



import scipy.stats as stats



def hypergeom_test(totoal_num_in_population, total_num_with_condition_in_pupulation, num_of_subset, num_with_condition_in_subset):
    
    less_pvalue = str(stats.hypergeom.cdf(int(num_with_condition_in_subset) ,int(totoal_num_in_population),int(total_num_with_condition_in_pupulation),int(num_of_subset)))
    great_pvalue = str(stats.hypergeom.sf(int(num_with_condition_in_subset) - 1,int(totoal_num_in_population),int(total_num_with_condition_in_pupulation),int(num_of_subset)))
    
    return less_pvalue, great_pvalue


#generation_cell_number :  all traced cell number at a certain generation
#generation_expression_cell_number :  detected expression cell number at a certain generation
#sublineage_cell_number : cell number which cells belong to sublineage at a certain generation
#observed_expression_lineage_cell_number: expression cell number which cells belong to sublineage at a certain generation

enrichement = (observed_expression_lineage_cell_number/sublineage_cell_number) / (expected_lineage_cell_number/generation_cell_number)
less_pvalue, great_pvalue = hypergeom_test(generation_cell_number, generation_expression_cell_number, sublineage_cell_number, observed_expression_lineage_cell_number)            

# output
# enrichement : enrichment score of specific lineage
# great_pvalue : pvalue

