# -*- coding: utf-8 -*-
"""
Created on Mon May  3 19:37:20 2021

@author: ZZG
"""
import pandas as pd
import os

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


def get_terminal_cell(cell_list):
    
    if cell_list[0][0] == "L":
        cell_list = cell_name_transfer(cell_list)
    else:
        pass
    
    unique_cell_set = set(cell_list)
    unique_cell_list = list(unique_cell_set)
    transfer_cell_list = cell_name_transfer(unique_cell_list)
    transfer_cell_set = set(transfer_cell_list)
    
    not_terminal_list = []
    for i in transfer_cell_list:
        len_i = len(i)
        for j in transfer_cell_list:
            if len(j) <= len_i:
                pass
            else:
                if j[:len_i] == i:
                    not_terminal_list.append(i)
                    break
    
    
    terminal_cell_set = transfer_cell_set - set(not_terminal_list)

    terminal_cell_list = list(terminal_cell_set)
    terminal_cell_list = cell_name_transfer(terminal_cell_list)

    return terminal_cell_list


file_path = r'./data'

# prepare data
relpace_name_sheet = pd.read_csv(os.path.join(file_path, 'binary_sheet.txt'), sep="\t")

cell_name_to_relpace_name_sheet = relpace_name_sheet.set_index("cell_name")
relpace_name_to_cell_name_sheet = relpace_name_sheet.set_index("replace_name")
funder_cell = {'P0':'L', 'AB':'L0', 'P1':'L1', 'EMS':'L10', 'P2':'L11', 'MS':'L100', 'E':'L101', 'C':'L110', 'D':'L1110', 'P4':'L1111', 'Z2':'L11110', 'Z3':'L11111'}
transform = {'a':'0', 'p':'1', 'l':'0', 'r':'1', 'd':'0', 'v':'1'}
    


# input 
#traced_cell_list : traced cell list
#exp_cell_list : expression cell list


traced_terminal_cell_list = get_terminal_cell(traced_cell_list)
exp_cell_list_replace = cell_name_transfer(exp_cell_list)
traced_terminal_cell_list_replace = cell_name_transfer(traced_terminal_cell_list)
    
transient_exp = []
for i in exp_cell_list_replace:
    not_continue_exp = False
    is_not_exist = True 

    if i not in traced_terminal_cell_list_replace: 
        if i+'0' not in exp_cell_list_replace or i+'1' not in exp_cell_list_replace:
            transient_exp.append(i)


# output 
# transient_exp : mother cell of transient expression site


