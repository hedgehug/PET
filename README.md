# Pathway Ensemble Tool: PET
<img width="781" alt="image" src="https://user-images.githubusercontent.com/16437494/207137637-32dec909-145c-4a3a-9421-57f62189dfb2.png">

## Requirement
* Python >= 3.6
* Numpy
* Scipy
* GSEA [Download](http://www.gsea-msigdb.org/gsea/downloads.jsp)
* DESeq2 [install](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

## File format
* Pathway file format, please refer to: [GSEA gene set format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Gene_Set_Database_Formats). Now supporting .gmx, .gmt, .grp. 
* Expression matrix file: tab-delimited text file. **Raw read count** is strongly recommended. First column is gene name, rest columns are sample expression. See exmaple in [data](https://github.com/hedgehug/PET/tree/main/data) folder. 

## Functions

### Run PET
PET is an ensemble tool for three pathway analysis methods: [GSEA](http://www.gsea-msigdb.org/gsea/index.jsp), Fisher Exact Test and [Enrichr](https://maayanlab.cloud/Enrichr/).

To run PET:
```
python PET.py -g GSEA_dir -o OUT_dir -c configuration.config -p Pathway_file.gmt
```

### Evaluate new methods

The data for constructing the benchmark is in ./data folder. We curated the raw read count from [ENCODE](https://www.encodeproject.org/) for each target and the corresponding controls. 

To evaluate new methods, please run the pathway analysis methods with expression profiles from one cell line, e.g. K562, and the pathway file from another cell line, e.g. HepG2. The paired comparison are listed below:

| Assay  |  Source |
|---|---|
|RNA-seq   |  HepG2 |
|  RNA-seq |K562   |
| ChIP-seq  |  HepG2 |
| ChIP-seq  |  K562 |

### Pathway p-value/FDR look-up

During the study, we noticed certain pathways are more oftenly assigned a higher significance, for example, a lower p-value or false positive rate, in pathway analysis methods. Here, we provided a quick look-up script for pathways in different tissues and cancer typpes.

To look up pathway significance:
```
python pathway_significane_look_up.py -p pathway_name -c BLCA -o output_file.txt
```



