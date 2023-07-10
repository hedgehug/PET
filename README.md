# Pathway Ensemble Tool: PET
<img width="781" alt="image" src="https://user-images.githubusercontent.com/16437494/207137637-32dec909-145c-4a3a-9421-57f62189dfb2.png">

## Requirement
* Python >= 3.6
* Numpy
* Scipy
* GSEA [Download](http://www.gsea-msigdb.org/gsea/downloads.jsp)
* DESeq2 [install](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
* statsmodels [install](https://www.statsmodels.org/stable/install.html)
* seaborn, optional, for plotting new method evaluation result only [install](https://seaborn.pydata.org/installing.html)

## Installation
If all dependencies already installed, please follow the [Tutorial](https://github.com/hedgehug/PET/blob/main/run_PET_tutorial.ipynb) to run PET.

To install PET from scratch, please refer to [Installation](https://github.com/hedgehug/PET/blob/main/Installation.md).

## Required files
* Pathway file, please provide pathway file in .gmt format, refer to: [GSEA gene set format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Gene_Set_Database_Formats). Please check [example pathway file](https://github.com/hedgehug/PET/blob/main/example/c2.cp.kegg.v2023.1.Hs.symbols.gmt). 
* Expression matrix file: tab-delimited text file. **Raw read count** is strongly recommended. First column is gene name, rest columns are sample expression. See exmaple in [example_data.txt](https://github.com/hedgehug/PET/tree/main/example/example_data.txt) folder. 


## Functions

### Run PET
PET ensembles pathway analysis results from three underlying methods, and by default, runs with the best practice pathway analysis settings.
We provide a tutorial PET analysis pipeline in the jupyter notebook [Tutorial](https://github.com/hedgehug/PET/blob/main/run_PET_tutorial.ipynb).

### Evaluate pathway analysis method using Benchmark

The data for constructing the benchmark is in [data](https://github.com/hedgehug/PET/tree/main/data) folder. We curated the raw read count from [ENCODE](https://www.encodeproject.org/) for each target and the corresponding controls. 


| Assay     | Source                                     | Number of targets | Curated target pathways                                                                                                                                                                                                                                                                                                                                                                    |
|-----------|--------------------------------------------|-------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| shRNA knockdown followed by RNA-seq| HepG2, K562                                |197| [HepG2_up](https://github.com/hedgehug/PET/blob/main/data/ENCODE_HepG2_RNA_up.gmt), [HepG2_down](https://github.com/hedgehug/PET/blob/main/data/ENCODE_HepG2_RNA_down.gmt), [K562_up](https://github.com/hedgehug/PET/blob/main/data/ENCODE_K562_RNA_up.gmt), [K562_down](https://github.com/hedgehug/PET/blob/main/data/ENCODE_K562_RNA_down.gmt)                                         |
| eCLIP-seq | HepG2, K562                                |73| [HepG2_pvalue](https://github.com/hedgehug/PET/blob/main/data/ENCODE_HepG2_eCLIP_pval.gmt), [HepG2_sigval](https://github.com/hedgehug/PET/blob/main/data/ENCODE_HepG2_eCLIP_signal_value.gmt), [K562_pvalue](https://github.com/hedgehug/PET/blob/main/data/ENCODE_K562_eCLIP_pval.gmt), [K562_sigval](https://github.com/hedgehug/PET/blob/main/data/ENCODE_K562_eCLIP_signal_value.gmt) |                                                                                                                                                   |
| ChIP-seq  | HepG2, K562, CH12.LX, MEL, GM12878, A549, HEK293 |374| [ChIP_qvalue](https://github.com/hedgehug/PET/blob/main/data/ENCODE_ChIP_seq_qvalue.gmt), [ChIP_score](https://github.com/hedgehug/PET/blob/main/data/ENCODE_ChIP_seq_peak_score.gmt)                                                                                                                                                                                                                                                                                   |

#### Run methods evaluated by Benchmark
We provided template scripts for running other methods

#### Evaluate new method using Benchmark

To evaluate new methods, please run the pathway analysis methods with expression profiles from one cell line, e.g. K562 RNA-seq, and the pathway file from another cell line, e.g. HepG2 RNA-seq. 

We provided a tutorial for how to evaluate a new pathway based on our benchmark here [evaluating new method tutorial](https://github.com/hedgehug/PET/blob/main/evaluate_new_method.ipynb)

