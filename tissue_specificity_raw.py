# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 23:01:53 2019



@author: ZZG
"""

import os
import copy
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
from collections import Counter
from scipy.stats import fisher_exact
from statsmodels.sandbox.stats.multicomp import multipletests


def hypergeom_test(totoal_num_in_population, total_num_with_condition_in_pupulation, num_of_subset, num_with_condition_in_subset):
    
    less_pvalue = str(stats.hypergeom.cdf(int(num_with_condition_in_subset) ,int(totoal_num_in_population),int(total_num_with_condition_in_pupulation),int(num_of_subset)))
    great_pvalue = str(stats.hypergeom.sf(int(num_with_condition_in_subset) - 1,int(totoal_num_in_population),int(total_num_with_condition_in_pupulation),int(num_of_subset)))
    
    return less_pvalue, great_pvalue

def pvalue_adjust(pvalue_matrix):

    mask_all = np.isfinite(pvalue_matrix.values)
    pval_corrected_all = np.full(pvalue_matrix.shape, np.nan)
    pval_corrected_all[mask_all] = multipletests(pvalue_matrix.values[mask_all], method='fdr_bh')[1]
    pval_corrected_all_df = pd.DataFrame(pval_corrected_all, index=pvalue_matrix.index, columns = pvalue_matrix.columns)
    return pval_corrected_all_df


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

# calculating cell lineage distance
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


def all_cell_lineage_order(cell_list):
    
    if cell_list[0][0] == "L":
        cell_list = cell_name_transfer(cell_list)
    else:
        pass
    
    cell_list_replace = cell_name_transfer(cell_list)
    cell_list_replace.sort()
    cell_list_replace = sorted(cell_list_replace, reverse=False, key=len)
    cell_list = cell_name_transfer(cell_list_replace)
    
    AB = []
    E = []
    MS = []
    D = []
    P = []
    C = []
    
    for i in cell_list:
        if i[0] == "A":
            AB.append(i)
        elif i[0] == "C":
            C.append(i)
        elif i[0] == "D":
            D.append(i)
        elif i[0] == "E":
            E.append(i)
        elif i[0] == "M":
            MS.append(i)  
        else:
            P.append(i)
               
    lineage_cell_list = AB+MS+E+C+D+P  
    
    return lineage_cell_list


def layer_lineage_order(cell_list):
    
    if cell_list[0][0] == "L":
        cell_list = cell_name_transfer(cell_list)
    else:
        pass
    
    cell_list_replace = cell_name_transfer(cell_list)
    cell_list_replace.sort()
    cell_list_replace = sorted(cell_list_replace, reverse=False, key=len)
    cell_list = cell_name_transfer(cell_list_replace)
    
    AB = []
    E = []
    MS = []
    D = []
    P = []
    C = []
    
    for i in cell_list:
        if i[0] == "A":
            AB.append(i)
        elif i[0] == "C":
            C.append(i)
        elif i[0] == "D":
            D.append(i)
        elif i[0] == "E":
            E.append(i)
        elif i[0] == "M":
            MS.append(i)  
        else:
            P.append(i)
    AB.sort()
    E.sort()
    MS.sort()
    D.sort()
    P.sort()
    C.sort()    
    lineage_cell_list = AB+MS+E+C+D+P  
    
    return lineage_cell_list
    

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
    terminal_cell_list = layer_lineage_order(terminal_cell_list)
    return terminal_cell_list


#file_path = r'./data'
file_path = r'I:\TF\Transcription-Factors\DU\NM\github\data'

# prepare data
relpace_name_sheet = pd.read_csv(os.path.join(file_path, 'binary_sheet.txt'), sep="\t")

cell_name_to_relpace_name_sheet = relpace_name_sheet.set_index("cell_name")
relpace_name_to_cell_name_sheet = relpace_name_sheet.set_index("replace_name")
funder_cell = {'P0':'L', 'AB':'L0', 'P1':'L1', 'EMS':'L10', 'P2':'L11', 'MS':'L100', 'E':'L101', 'C':'L110', 'D':'L1110', 'P4':'L1111', 'Z2':'L11110', 'Z3':'L11111'}
transform = {'a':'0', 'p':'1', 'l':'0', 'r':'1', 'd':'0', 'v':'1'}

Cell = pd.read_csv(os.path.join(file_path, r'all_embryonic_cell.txt'), sep='\t', index_col=0)
tissue = Cell['Tissue']        # fate file
death_cell_list = list(tissue[tissue== 'Dea'].index)

tissue_list = ['Neuron', 'Pharynx', 'Skin', 'Muscle', 'Int', 'Germ',  'Dea', 'other']
tissue_list_withot_death = ['Neuron', 'Pharynx', 'Skin', 'Muscle', 'Int', 'Germ', 'other']

cell_list = list(tissue.index.values)
terminal_cell = list(Cell[Cell['Terminal'] == 'terminal'].index.values)
cell_list_replace = cell_name_transfer(cell_list)
terminal_cell_replace = cell_name_transfer(terminal_cell)



#even lineage
new_fate = []
new_terminal_cell_list = []
for each_cell in tissue.index.values:
    
    cell_replace = cell_name_transfer(each_cell)
    
    if cell_replace in terminal_cell_replace:
    
        if cell_replace[-1] == '1':
            sister_replace = cell_replace[:-1]+'0'
        else:
            sister_replace = cell_replace[:-1]+'1'
            
        childern_list = []
        for i in cell_list_replace:
            if sister_replace != i and sister_replace in i:
                childern_list.append(i)
        
        if len(childern_list) >= 2:
            new_fate.extend([each_cell, each_cell])
            new_terminal_cell_list.extend([each_cell, each_cell])
        else:
            new_fate.append(each_cell)
            new_terminal_cell_list.extend([each_cell])
    else:
        new_fate.append(each_cell)

new_fate_pd = tissue.loc[new_fate]


# output even lineage
new_fate_pd_output = pd.DataFrame(new_fate_pd)
new_fate_pd_output.columns = ['cell fate']



new_fate_cell_list = list(new_fate_pd.index.values)
new_fate_cell_list_replace = cell_name_transfer(new_fate_cell_list)
new_fate_pd.index = new_fate_cell_list_replace


new_terminal_cell_replace_name_list = cell_name_transfer(new_terminal_cell_list)


# each cell has terminal cell 
terminal_cell_contain = {}
for each_cell in new_fate_cell_list_replace:
    
    terminal_cell_list_ = []
    
    for each_cell_ in new_terminal_cell_replace_name_list:
        
        if each_cell in each_cell_:
            terminal_cell_list_.append(each_cell_)
    
    terminal_cell_contain[each_cell] = terminal_cell_list_
    





tissue.index = cell_list_replace
# each cell has terminal cell fate
terminal_cell_contain_fate = {}
for each_cell, cell_list_ in terminal_cell_contain.items():
    
    terminal_cell_contain_fate[each_cell] = list(tissue.loc[cell_list_].values)
    
    
    
    
    
# each cell has terminal cell fate num
terminal_cell_contain_fate_count = {}    
for each_cell, cell_fate_ in terminal_cell_contain_fate.items():
    
    terminal_cell_contain_fate_count[each_cell] = dict(Counter(cell_fate_))
    


fate_contain = []
for each_cell, each_cell_fate_contain in terminal_cell_contain_fate_count.items():
    
    terminal_cell_list_num = 0.0
    
    fate_each_cell_ = []
    
    for each_fate in tissue_list:
        if each_fate not in each_cell_fate_contain.keys():
            fate_each_cell_.append(0)
        else:
            fate_each_cell_.append(each_cell_fate_contain[each_fate])
            
    fate_contain.append([cell_name_transfer(each_cell)] + fate_each_cell_)

fate_contain_pd = pd.DataFrame(fate_contain)
fate_contain_pd.columns = ['cell_name'] + tissue_list
fate_contain_pd = fate_contain_pd.set_index('cell_name')


# new
death_pd = fate_contain_pd
death_pd_percent = fate_contain_pd.div(fate_contain_pd.sum(axis=1), axis=0)


##no_death_pd = fate_contain_pd[['neuron', 'glia', 'gut', 'hypodermis', 'pharynx',  'muscle', 'other', 'germ-line']]
#no_death_pd = fate_contain_pd[['neuron', 'glia', 'gut', 'hypodermis', 'muscle', 'other', 'germ-line']]
no_death_pd = fate_contain_pd[tissue_list_withot_death]
no_death_pd_percent = no_death_pd.div(no_death_pd.sum(axis=1), axis=0)


fate_probability = pd.concat([no_death_pd_percent, death_pd_percent[['Dea']]], axis=1)
fate_probability = fate_probability.fillna(0)





stage_cell = pd.read_csv(os.path.join(file_path, r'stage_cell.txt'), sep='\t', index_col=1)
stage_list = list(set(stage_cell.index.values))
stage_list.sort()
stage_cell_list = {}
for each in stage_list:
    stage_cell_list[each] = list(stage_cell.loc[each]['Cell_name'].values)


#fate_probability = pd.read_csv(os.path.join(output_path, r'fate_percent_new_ignore_death.txt'), sep='\t', index_col=0)



#将每个细胞距离终末分化细胞的谱系距离作为其所处的发育阶段，作为选取对照细胞的标准。对照细胞的选取，选取和自己相同发育阶段的细胞
layer_dict = {}
cell_list_for_layer_calculate = copy.deepcopy(cell_list)
cell_list_for_layer_calculate = all_cell_lineage_order(list(set(list(cell_list_for_layer_calculate))))

for cell in cell_list_for_layer_calculate:
    
    cell_replace = cell_name_transfer(cell)
    
    terminal_cell_ = []
    for i in terminal_cell_replace:
        if cell_replace in i:
            terminal_cell_.append(i)
            
    cld_ = []
    for each_cell in terminal_cell_:
        cld_.append(cell_lineage_distance(cell, cell_name_transfer(each_cell)))
    
    layer_dict.setdefault(min(cld_), [])
    layer_dict[min(cld_)].append(cell)
    




#exp_matrix = pd.read_csv(r'Z:\zgzhao\work project\MEMBERS\ZZG\2019-03-19-filter_TF_for_analysis\#exp_matrix_binary_all_linear_specific.txt', sep='\t', index_col=0)
exp_matrix = pd.read_csv(r'I:\TF\Transcription-Factors\DU\NM\exp_matrix_266.txt', sep='\t', index_col=0)
exp_matrix[exp_matrix>0]=1
cell_list = list(exp_matrix.index.values)
terminal_cell_list_given = get_terminal_cell(list(exp_matrix.index.values))
cell_list_replace = cell_name_transfer(cell_list)
terminal_cell_list_given_replace = cell_name_transfer(terminal_cell_list_given)

terminal_cell_list_del_death = list(set(terminal_cell_list_given) - set(death_cell_list))

stage_cell_list_my = {}
for key, values in stage_cell_list.items():
    stage_cell_list_my[key] = list(set(values) & set(cell_list))

stage_cell_list_r = {}
for key, values in stage_cell_list_my.items():
    stage_cell_list_r[key] = cell_name_transfer(values)


# input control
layer_dict = stage_cell_list_my

exp_matrix = exp_matrix.fillna(0)



output_pvalue_hyp = []
output_flag = []

count = 0
for each_TF in exp_matrix.columns.values:
    
    
    gene_ = exp_matrix[each_TF].loc[cell_list].dropna()
    gene_exp_cell_ = list(gene_[gene_==1].index.values)   # exp cell list
    exp_cell_num_ = len(gene_exp_cell_)                   # exp cell num
    
    if exp_cell_num_ != 0:
        
        count_death = 0
        for z in gene_exp_cell_:
            if z in death_cell_list:
                count_death += 1
        
        if count_death == exp_cell_num_:
       
            output_pvalue_ = [each_TF, exp_cell_num_] + [np.nan] * len(tissue_list)
            output_flag_ = [each_TF, exp_cell_num_] + [np.nan] * len(tissue_list)     
        else:
            
            #ALL
            gene_ = exp_matrix[each_TF].dropna()
            gene_exp_cell_ = list(gene_[gene_==1].index.values)   # exp cell list
            exp_cell_num_ = len(gene_exp_cell_)                   # exp cell num 

            output_pvalue_hyp_ = [each_TF, exp_cell_num_]
            output_flag_ = [each_TF, exp_cell_num_]
            
            for each_fate in tissue_list:
                
                fate_pd_ = []
                for each_exp_cell in gene_exp_cell_:
                    
                    if each_exp_cell in death_cell_list and each_fate != 'Dea':
                        pass
                    else:
                        probability_ = fate_probability.loc[each_exp_cell][each_fate]
                        
                        control_cell_list = []
                        for key, values in stage_cell_list_my.items():
                            if each_exp_cell in values:
                                control_cell_list.extend(values)
                        control_cell_list = list(set(control_cell_list))
                        
                        exp_pd_ = fate_probability.loc[control_cell_list][each_fate]
                        ctr_mean_ = exp_pd_.mean()
           
                        fate_pd_.append([each_exp_cell, probability_, ctr_mean_])

                
                fate_pd_pd = pd.DataFrame(fate_pd_)
                fate_pd_pd.columns = ['cell_name', 'obs', 'ctr']
                fate_pd_pd = fate_pd_pd.set_index('cell_name')
                
                obs_sum_ = fate_pd_pd.sum()['obs']
                exp_sum_ = fate_pd_pd.sum()['ctr']
                
                obs_sum_ = np.round(obs_sum_,0)
                exp_sum_ = np.round(exp_sum_,0)
           
            
                
                # O/E
                if exp_sum_ == 0 and obs_sum_ == 0:
                    f_stats = np.nan
                elif exp_sum_ == 0 and obs_sum_ != 0:
                    f_stats = np.inf
                else:
                    f_stats = obs_sum_ / exp_sum_
               
                # pvalue
                fisher_exp_obs = np.round(fate_pd_pd.sum()['obs'], 0)
                fisher_exp_exp = np.round(fate_pd_pd.sum()['ctr'], 0)
                
                ctr_sum_ = fate_probability.loc[cell_list][each_fate].sum()
                

                ctr_sum_ = fate_pd_pd.sum()['ctr']*len(cell_list)/len(gene_exp_cell_)
                ctr_sum_ = np.round(ctr_sum_,0)
                less_pvalue, great_pvalue = hypergeom_test(len(cell_list), ctr_sum_, \
                                                      len(gene_exp_cell_), fisher_exp_obs)            

                output_pvalue_hyp_.append(great_pvalue)
                
                #oddsratio_ = np.nanmean(fate_pd_pd['obs']) / np.nanmean(fate_pd_pd['ctr'])
                oddsratio_ = f_stats
                output_flag_.append(oddsratio_)
    else:

    
        output_pvalue_hyp_ = [each_TF, exp_cell_num_] + [np.nan] * len(tissue_list)
        output_flag_ = [each_TF, exp_cell_num_] + [np.nan] * len(tissue_list)
    

    output_pvalue_hyp.append(output_pvalue_hyp_)
    output_flag.append(output_flag_)
    
    count+=1
    print(count)
    
    
#####################################          filter embs        ###########################################    
# remain_embs = list(set(exp_matrix.columns.values) - set(emb_filter))
remain_embs = list(set(exp_matrix.columns.values))
remain_embs.sort()
#####################################          filter embs        ###########################################    






# ALL
all_odd_pd = pd.DataFrame(output_flag)
all_odd_pd.columns = ['gene_name', 'exp_cell_num'] + tissue_list
all_odd_pd = all_odd_pd.set_index(['gene_name', 'exp_cell_num'])
#output_flag_pd = output_flag_pd.replace(np.inf, np.nan)


all_pvalue_hyp_pd = pd.DataFrame(output_pvalue_hyp)
all_pvalue_hyp_pd.columns = ['gene_name', 'exp_cell_num'] + tissue_list
all_pvalue_hyp_pd = all_pvalue_hyp_pd.set_index(['gene_name',  'exp_cell_num'])


all_odd_pd = all_odd_pd.loc(axis=0)[remain_embs, :]
all_pvalue_hyp_pd = all_pvalue_hyp_pd.loc(axis=0)[remain_embs, :]




all_pvalue_hyp_pd = all_pvalue_hyp_pd.astype(float)
all_qvalue_hyp_pd = pvalue_adjust(all_pvalue_hyp_pd)


qvalue_cutoff = 0.01
odd_cutoff =2
sig_pd = copy.deepcopy(all_qvalue_hyp_pd)
sig_fate = (all_odd_pd>=odd_cutoff) & (all_qvalue_hyp_pd<qvalue_cutoff)
sig_fate = sig_fate.astype(int)

    

fate_list_output = ['Neuron', 'Pharynx', 'Skin', 'Muscle', 'Int', 'Germ', 'Dea', 'other']

all_qvalue_hyp_pd = all_qvalue_hyp_pd.astype(float)
all_odd_pd = all_odd_pd.astype(float)

qvalue_cutoff = 0.05
odd_cutoff =2
sig_pd = copy.deepcopy(all_qvalue_hyp_pd)
sig_fate = (all_odd_pd>=odd_cutoff) & (all_qvalue_hyp_pd<qvalue_cutoff)
sig_fate = sig_fate.astype(int)




mask = sig_fate[fate_list_output]


#step1 = list(mask[(mask.iloc[:, :6].sum(axis=1) != 0)].index)

emb_filter = []

add_exclusive = []
des_exclusive = []
count = 0

exp_percent_fate = []
for each_TF in exp_matrix.columns.values:
    
    if each_TF in emb_filter:
        pass
    else:
    
        gene_ = exp_matrix[each_TF]
        gene_exp_cell_ = list(gene_[gene_==1].index.values)   # exp cell list
        
        exp_cell_num_ = len(gene_exp_cell_)
        
        exp_cell_time = list(set(terminal_cell_list_del_death) & set(gene_exp_cell_))
        
        exp_fate_num = fate_probability.loc[gene_exp_cell_].sum(axis=0)
        exp_cell_death_and_other = sum([exp_fate_num.loc['Dea'], exp_fate_num.loc['other']])
        
        exp_cell_del_death_and_other = exp_fate_num.sum() - exp_cell_death_and_other
        
        des_flag = True
        
        exp_percent_ = []
        for i in tissue_list:
            
            ratio_ = exp_fate_num.loc[i] / exp_cell_del_death_and_other
            
            if exp_fate_num.loc[i] < 15 and i not in ['other', 'Dea', 'Germ']:
                des_flag = False
            
            
            if ratio_ > 0.8 and exp_cell_del_death_and_other >= 2:
                
                if each_TF in mask[i][mask[i]==1].index.values:
                    pass
                else:
                    add_exclusive.append([each_TF, i])
            exp_percent_.append(ratio_)
        if des_flag:
            des_exclusive.append(each_TF)
            
        if len(exp_cell_time) == 0 and len(exp_fate_num[exp_fate_num!=0]) > 1:
            des_exclusive.append(each_TF)
            
        exp_percent_fate.append(exp_percent_)
        
add_exclusive_pd = pd.DataFrame(add_exclusive)
add_exclusive_pd.columns = ['TF', 'fate']
add_exclusive_pd = add_exclusive_pd.drop_duplicates()

exp_percent_fate_pd = pd.DataFrame(exp_percent_fate)
exp_percent_fate_pd.columns = tissue_list
exp_percent_fate_pd.index = exp_matrix.columns.values


for i in add_exclusive_pd.values:
    
    tf_ = i[0]
    fate_ = i[1]
    
    mask.loc[tf_,fate_] = 1


for i in des_exclusive:
    
    mask.loc[i] = 0




mask.reset_index(inplace=True)
mask = mask.set_index('gene_name')

mask = mask.drop_dupliate()






