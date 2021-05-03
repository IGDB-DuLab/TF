# Single-cell Transcription Factor Analysis

This repository contains the source code for a brief description of the analyses of our study titled "A 4D Single-Cell Protein Atlas of Transcription Factors Delineates  Spatiotemporal Patterning During Embryogenesis" by Xuehua Ma#, Zhiguang Zhao#, Long Xiao#, Weina Xu, Yahui Kou, Yanping Zhang, Gang Wu, Yangyang Wang, and Zhuo Du*  

## Description of files
- [Lineage_specificity.py](https://github.com/IGDB-DuLab/Zhao-TF/blob/main/Lineage_specificity.py)  Brief description for identifying lineage specific TF.
- [AP_asymmetry.py](https://github.com/IGDB-DuLab/Zhao-TF/blob/main/AP_asymmetry.py) Brief description for identifying AP asymmetry TF.
- [Tissue_specificity.py](https://github.com/IGDB-DuLab/Zhao-TF/blob/main/Tissue_specificity.py) Brief description for identifying tissue specific TF.
- [Transient_expression.py](https://github.com/IGDB-DuLab/Zhao-TF/blob/main/Transient_expression.py) Brief description for identifying transient expression TF.
- [Expression_similarity.py](https://github.com/IGDB-DuLab/Zhao-TF/blob/main/Expression_similarity.py) Performing analyses for calculating similarity between two TFs.
- [Umap_and_louvain.py](https://github.com/IGDB-DuLab/Zhao-TF/blob/main/Umap_and_louvain.py) Performing analyses related to ontologically-related cells.
- [Calculate_cell_lineage_distance.py](https://github.com/IGDB-DuLab/Zhao-TF/blob/main/Calculate_cell_lineage_distance.py) Calculating cell lineage distance between cells.
- [Calculate_jensenshannon.py](https://github.com/IGDB-DuLab/Zhao-TF/blob/main/Calculate_jensenshannon.py) Calculating cellular state divergences .
- [data](https://github.com/IGDB-DuLab/Zhao-TF/blob/main/data) folder contains the raw data needed for various analysis.

## Prerequisites
All the dependencies are listed in [requirements.txt](https://github.com/genetics-dulab/scCAL/blob/main/requirements.txt). 

Most analyses were performed under python 3.7.1

Additional Python packages required for analysis include:

- pandas (>= 1.0.3)
- scipy (>= 1.2.1)
- umap (>= 0.4.6)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](https://github.com/genetics-dulab/scCAL/blob/main/LICENSE) file for details.

