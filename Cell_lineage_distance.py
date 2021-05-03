# -*- coding: utf-8 -*-
"""
Created on Mon May  3 18:52:55 2021

@author: ZZG
"""

import os
import pandas as pd

def except_cell(cell_name):
    
    transform_name = None
    if cell_name[0] == 'L':
        len_cell_name = len(cell_name)
        for i in range(1,len_cell_name+1):
            if cell_name[:-i] in list(relpace_name_to_cell_name_sheet.index.values):
                transform_name = (relpace_name_to_cell_name_sheet.loc[[cell_name[:-i]]].values[0])[0]
                break
            
        append_name = cell_name[len_cell_name-i:]
        
        for j in append_name:
            if j == '0':
                transform_name+='a'
            else:
                transform_name+='p'
                
    else:
        len_cell_name = len(cell_name)
        for i in range(1,len_cell_name+1):
            if cell_name[:-i] in list(cell_name_to_relpace_name_sheet.index.values):
                transform_name = (cell_name_to_relpace_name_sheet.loc[[cell_name[:-i]]].values[0])[0]
                break
            
        append_name = cell_name[len_cell_name-i:]
        
        for j in append_name: 
            transform_name += transform[j]
            
    return transform_name

def cell_name_transfer(cell_list):
    
    transfer_cell_list = []
    if len(cell_list) == 0:
        pass
        
    elif type(cell_list) != list:
        if cell_list[0] == "L":
            if cell_list in list(relpace_name_to_cell_name_sheet.index.values):
                transfer_cell_list = (relpace_name_to_cell_name_sheet.loc[[cell_list]].values[0])[0]
            else:
                transfer_cell_list = except_cell(cell_list)
        else:
            if cell_list in list(cell_name_to_relpace_name_sheet.index.values):
                transfer_cell_list = (cell_name_to_relpace_name_sheet.loc[[cell_list]].values[0])[0]
            else:
                transfer_cell_list = except_cell(cell_list)
                
    else:
        if cell_list[0][0] == "L":
            for i in cell_list:
                if i in list(relpace_name_to_cell_name_sheet.index.values):
                    transfer_cell_list.append((relpace_name_to_cell_name_sheet.loc[[i]].values[0])[0])
                else:
                    transfer_cell_list.append(except_cell(i))
            
        else:
            for i in cell_list:
                if i in list(cell_name_to_relpace_name_sheet.index.values):
                    transfer_cell_list.append((cell_name_to_relpace_name_sheet.loc[[i]].values[0])[0])
                else:
                    transfer_cell_list.append(except_cell(i))

    return transfer_cell_list

def cell_lineage_distance(cell_1, cell_2):
    
    if cell_1[0] == 'L':
        pass
    else:
        cell_1 = cell_name_transfer(cell_1)
        cell_2 = cell_name_transfer(cell_2)
    
    cell_1_list = list(cell_1)
    cell_2_list = list(cell_2)
    
    len_cell_1_list = len(cell_1_list)
    len_cell_2_list = len(cell_2_list)
    
    min_len = min(len_cell_1_list, len_cell_2_list)
    
    lowest_common_ancestor=''
    for i in range(min_len):
        if cell_1_list[i] == cell_2_list[i]:
            lowest_common_ancestor += cell_1_list[i]
        else:
            break
            
    distance_1 = len_cell_1_list - len(lowest_common_ancestor)
    distance_2 = len_cell_2_list - len(lowest_common_ancestor)
        
    distance = distance_1+distance_2
    
    return distance



file_path = r'./data'

# prepare data
relpace_name_sheet = pd.read_csv(os.path.join(file_path, 'binary_sheet.txt'), sep="\t")

cell_name_to_relpace_name_sheet = relpace_name_sheet.set_index("cell_name")
relpace_name_to_cell_name_sheet = relpace_name_sheet.set_index("replace_name")
funder_cell = {'P0':'L', 'AB':'L0', 'P1':'L1', 'EMS':'L10', 'P2':'L11', 'MS':'L100', 'E':'L101', 'C':'L110', 'D':'L1110', 'P4':'L1111', 'Z2':'L11110', 'Z3':'L11111'}
transform = {'a':'0', 'p':'1', 'l':'0', 'r':'1', 'd':'0', 'v':'1'}


# Calculate cell lineage distancce

# cell 1 : CellLineageName or cell list
# cell 2 : CellLineageName or cell list

cld = cell_lineage_distance(cell_1, cell_2)

# output
# cld : cell lineage distane between two cells
