# Installation 

We strongly recommend Conda to install all packages required for running PET. 

* Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html), with Python version >=3.7.

* Download PET repository.
```
git clone https://github.com/hedgehug/PET.git
```
* Enter the downloaded folder, create a new conda environment for PET (in terminal):
```
conda env create -f environment.yml
```
* Install jupyter notebook and add environment to jupyter notebook :
```
conda install -c anaconda jupyter
python -m ipykernel install --user --name=pet
```
We will also need GSEA to run PET. Please download [GSEA for the
command line (all platforms)](http://www.gsea-msigdb.org/gsea/downloads.jsp).

After installation, please follow [Tutorial](https://github.com/hedgehug/PET/blob/main/run_PET_tutorial.ipynb) to run PET.