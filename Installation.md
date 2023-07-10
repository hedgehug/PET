# Installation 

We strongly recommend Conda to install all packages required for running PET. 

* Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html), with Python version >=3.6.

* Download PET repository through Git (in terminal) or directly downloaded as [Zipped file](https://github.com/hedgehug/PET/archive/refs/heads/main.zip).
```
git clone https://github.com/hedgehug/PET.git
```
* Enter the downloaded folder, create a new conda environment for PET, this step will install all Python and R dependencies.
```
cd PET/
# for MacOS
conda env create -f environment_mac.yml
# for Linux, tested in Ubuntu
conda env create -f environment_linux.yml
```
* Activate the newly created PET conda environment
```
conda activate pet
```

We will also need GSEA to run PET. Please download [GSEA for the
command line (all platforms)](http://www.gsea-msigdb.org/gsea/downloads.jsp). JDK 11 is required for running GSEA command line tool, [download here](https://www.oracle.com/java/technologies/downloads/).

After installation, please follow [Tutorial](https://github.com/hedgehug/PET/blob/main/run_PET_tutorial.ipynb) to run PET. 

To open jupyter notebook, run following command in PET/:
```
jupyter notebook
```