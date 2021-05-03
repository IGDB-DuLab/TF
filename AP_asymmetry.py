# -*- coding: utf-8 -*-
"""
Created on Mon May  3 19:27:27 2021

@author: ZZG
"""


from scipy.stats import binom_test


#a_count : number of anterior lineage satisfied by thresholds of difference(0.25, 0.5 or 0.75) in expression frequency between anterior and posterior (AP) lineages
#p_count : number of posterior lineage satisfied by thresholds of difference(0.25, 0.5 or 0.75) in expression frequency between anterior and posterior (AP) lineages

if a_count > p_count:
    ratio = a_count / p_count
else:
    ratio = p_count / a_count
    
pvalue = binom_test(a_count, a_count+p_count, 0.5)


# output
# ratio : the asymmetry ratio
# pvalue : pvalue
