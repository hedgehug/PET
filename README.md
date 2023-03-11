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
## File format
* Pathway file format, please refer to: [GSEA gene set format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Gene_Set_Database_Formats). Now supporting .gmt format. 
* Expression matrix file: tab-delimited text file. **Raw read count** is strongly recommended. First column is gene name, rest columns are sample expression. See exmaple in [example_data.txt](https://github.com/hedgehug/PET/tree/main/example) folder. 

## Installation
To install PET, simply run
```
git clone https://github.com/hedgehug/PET.git
```

## Functions

### Run PET
PET ensembles pathway analysis results from three underlying methods, and by default, runs with the best practice pathway analysis settings.
We provide a tutorial PET analysis pipeline in the jupyter notebook [Tutorial](https://github.com/hedgehug/PET/blob/main/run_PET_tutorial.ipynb).

### Evaluate new methods

The data for constructing the benchmark is in ./data folder. We curated the raw read count from [ENCODE](https://www.encodeproject.org/) for each target and the corresponding controls. 

To evaluate new methods, please run the pathway analysis methods with expression profiles from one cell line, e.g. K562, and the pathway file from another cell line, e.g. HepG2. The paired comparison are listed below:

| Assay     | Source                                     | Number of targets |
|-----------|--------------------------------------------|-------------------|
| shRNA knockdown followed by RNA-seq| HepG2, K562                                |197|
| eCLIP-seq | HepG2, K562                                |73|
| ChIP-seq  | HepG2, K562, CH12.LX, MEL, GM12878, A549, HEK293 |374|

We provided a tutorial for how to evaluate a new pathway based on our benchmark here [evaluating new method tutorial](https://github.com/hedgehug/PET/blob/main/evaluate_new_method.ipynb)
