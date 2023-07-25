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

Notes: 
* For RNA-seq, the pathway files contain top 200 up/down-regulated differentially expressed genes after shRNA knockdown with specific target.
* For eCLIP-seq and ChIP-seq, the pathway files contain annotated genes near the top peaks, sorted based on peak signal value/p-value/q-value.

#### Run methods evaluated by Benchmark
We provided template scripts for running other methods in [template scripts](https://github.com/hedgehug/PET/tree/main/template_script).

#### Evaluate new method using Benchmark

To evaluate new methods, please run the pathway analysis methods with expression profiles from one cell line, e.g. K562 RNA-seq, and the pathway file from another cell line, e.g. HepG2 RNA-seq. 

We provided a tutorial for how to evaluate a new pathway based on our benchmark here [evaluating new method tutorial](https://github.com/hedgehug/PET/blob/main/evaluate_new_method.ipynb)

## Copyright
Copyright (c) 2023 Purdue University All rights reserved.
Developed by:  Kazemian lab (Luopin Wang), College of Science, Purdue University, https://kazemianlab.com/
Permission to use, copy, modify, and distribute this software and associated documentation (the “Software”) for educational, research and non-profit purposes, in each case by a government agency, educational institution, or other non-profit organization, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
 
Permission to make commercial use of this Software or permission for any use by a for-profit entity may be requested by contacting:
Purdue Research Foundation
Office of Technology Commercialization
The Convergence Center
101 Foundry Dr., Suite 2500
West Lafayette, IN 47906
ATTN: Business Development Manager
otcip@prf.org
 
DISCLAIMER. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.
 
LIMITATION OF LIABILITY. TO THE MAXIMUM EXTENT PERMITTED BY LAW, IN NO EVENT WILL PURDUE OR PURDUE RESEARCH FOUNDATION BE LIABLE TO ANY USER OF THIS CODE FOR ANY INCIDENTAL, CONSEQUENTIAL, EXEMPLARY OR PUNITIVE DAMAGES OF ANY KIND, LOST GOODWILL, LOST PROFITS, LOST BUSINESS AND/OR ANY INDIRECT ECONOMIC DAMAGES WHATSOEVER, REGARDLESS OF WHETHER SUCH DAMAGES ARISE FROM CLAIMS BASED UPON CONTRACT, NEGLIGENCE, TORT (INCLUDING STRICT LIABILITY OR OTHER LEGAL THEORY), A BREACH OF ANY WARRANTY OR TERM OF THIS AGREEMENT, AND REGARDLESS OF WHETHER PURDUE OR PURDUE RESEARCH FOUNDATION WAS ADVISED OR HAD REASON TO KNOW OF THE POSSIBILITY OF INCURRING SUCH DAMAGES IN ADVANCE.