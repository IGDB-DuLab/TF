# -*- coding: utf-8 -*-
"""
Created on Mon May  3 21:08:27 2021

@author: ZZG
"""

from scipy.spatial.distance import jensenshannon

#expression_matrix : cellular expression matrix row(TF) * col(cell)
#cell1 : cell1 lineage name
#cell2 : cell2 lineage name

jsd = jensenshannon(expression_matrix[cell1], expression_matrix[cell2])

# output
# jsd : jensen-shannon divergence between two cells
