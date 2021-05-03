# Single-cell Transcription Factor Analysis

This repository contains the source code for a brief description of the analyses of our study titled "A 4D Single-Cell Protein Atlas of Transcription Factors Delineates  Spatiotemporal Patterning During Embryogenesis" by Xuehua Ma#, Zhiguang Zhao#, Long Xiao#, Weina Xu, Yahui Kou, Yanping Zhang, Gang Wu, Yangyang Wang, and Zhuo Du*  

## Description of files
- [Lineage_specificity.py](https://github.com/IGDB-DuLab/Zhao-TF/lineage_specificity.py)  Compensation for depth-dependent attenuation of fluorescence intensity
- [AP_asymmetry.py](https://github.com/genetics-dulab/scCAL/blob/main/data_stats.py) Calculation of data distribution, CV, infomation content and variability.
- [Tissue_specificity.py](https://github.com/genetics-dulab/scCAL/blob/main/CAL_lineage.py) Performing analyses related to lineage fate patterning and anterior-posterior fate asymmetry of CALs .
- [Transient_expression.py](https://github.com/genetics-dulab/scCAL/blob/main/CAL_tissue.py) Performing analyses related to tissue convergence and tissue heterogeneity .
- [Expression_similarity.py](https://github.com/genetics-dulab/scCAL/blob/main/CAL_symmetry.py) Performing analyses related to symmetry predetermination.
- [Umap_and_louvain.py](https://github.com/genetics-dulab/scCAL/blob/main/CDCA.py) Performing analyses related to chromatin co-dynamic region.
- [Calculate_cell_lineage_distance.py](https://github.com/genetics-dulab/scCAL/blob/main/CDCA.py) Performing analyses related to chromatin co-dynamic region.
- [Calculate_jensenshannonpy](https://github.com/genetics-dulab/scCAL/blob/main/CDCA.py) Performing analyses related to chromatin co-dynamic region.
- [data](https://github.com/genetics-dulab/scCAL/blob/main/data) folder contains the raw data needed for various analysis.

## Prerequisites
All the dependencies are listed in [requirements.txt](https://github.com/genetics-dulab/scCAL/blob/main/requirements.txt). 

Most analyses were performed under python 3.7.1

Additional Python packages required for analysis include:

- pandas (>= 1.0.3)
- scipy (>= 1.2.1)
- umap (>= 0.4.6)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](https://github.com/genetics-dulab/scCAL/blob/main/LICENSE) file for details.

