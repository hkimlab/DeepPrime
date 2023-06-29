<div align="center">
  
  <img src="https://github.com/hkimlab/DeepPrime/blob/main/docs/images/DeepPrime_logo.png?raw=true" width="600"/>

**Developed by Hyongbum Henry Kim's lab** </br>  
[![Python 3.6](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/release/python-360/)
<a href="https://doi.org/10.1016/j.cell.2023.03.034"><img src="https://img.shields.io/badge/Published-Yu et al., Cell, 2023-blue"></a>

<div align="left">

# Table of Contents
- [About](#about)
- [DeepPrime webtool](#deepprime-webtool)
- [Python package for using DeepPrime: GenET](#python-package-for-using-deepprime-genet)
- [Installation from source code](#installation-from-source-code)
- [Usage](#usage)


## About:

DeepPrime is a deep-learning-based prime editing efficiency prediction tool developed in Laboratory of Genome Editing, Yonsei University. 

It greatly expands upon the previous start-of-the-art pegRNA activity prediction model, DeepPE, which was limited to a specific set of edit type and length combinations.

DeepPrime is developed to predict efficiencies of a nearly all feasible combinations of pegRNA designs. We integrated CNN & RNN to extract inter-sequence features between target DNA and corresponding pegRNA. DeepPrime was trained using 259K pegRNAs with PBS lengths ranging from 1 to 17, RT lengths ranging from 1 to 50, Edit positions ranging from 1 to 30, and editing lengths ranging from 1 to 3.

## DeepPrime webtool
The [webtool app](http://deepcrispr.info/DeepPrime/) can accommodate most applications using default parameters and using the most appropriate primed editing (PE) model for your experimental conditions. It can evaluate all possible prime editing guide RNAs (pegRNAs) for a given target according to the predicted the prime editing efficiency, DeepPrime score.

## Python package for using DeepPrime: GenET
[GenET (Genome Editing Toolkit)](https://github.com/Goosang-Yu/genet) is a library of various python functions for the purpose of analyzing and evaluating data from genome editing experiments.

### Installation

```python
# Create virtual env for genet.
# python 3.8 was tested. 
conda create -n genet python=3.8
conda activate genet

# install genet package in your env.
pip install genet -f https://download.pytorch.org/whl/cu113/torch_stable.html git+https://github.com/Goosang-Yu/genet-models.git

# install ViennaRNA package for prediction module
conda install viennarna
```

### How to use DeepPrime using GenET
```python
from genet import predict as prd

# Place WT sequence and Edited sequence information, respectively.
# And select the edit type you want to make and put it in.
#Input seq: 60bp 5' context + 1bp center + 60bp 3' context (total 121bp)

seq_wt   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'
seq_ed   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'
alt_type = 'sub1'

df_pe = prd.pe_score(seq_wt, seq_ed, alt_type)
df_pe.head()
```
output:
|    | ID     | WT74_On                                                                    | Edited74_On                                                                |   PBSlen |   RTlen |   RT-PBSlen |   Edit_pos |   Edit_len |   RHA_len |   type_sub |   type_ins |   type_del |     Tm1 |     Tm2 |   Tm2new |      Tm3 |     Tm4 |      TmD |   nGCcnt1 |   nGCcnt2 |   nGCcnt3 |   fGCcont1 |   fGCcont2 |   fGCcont3 |   MFE3 |   MFE4 |   DeepSpCas9_score |   PE2max_score |
|---:|:-------|:---------------------------------------------------------------------------|:---------------------------------------------------------------------------|---------:|--------:|------------:|-----------:|-----------:|----------:|-----------:|-----------:|-----------:|--------:|--------:|---------:|---------:|--------:|---------:|----------:|----------:|----------:|-----------:|-----------:|-----------:|-------:|-------:|-------------------:|------------------:|
|  0 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxxxxCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx |        7 |      35 |          42 |         34 |          1 |         1 |          1 |          0 |          0 | 16.191  | 62.1654 |  62.1654 | -277.939 | 58.2253 | -340.105 |         5 |        16 |        21 |    71.4286 |    45.7143 |    50      |  -10.4 |   -0.6 |            45.9675 |         0.0202249 |
|  1 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxxxCCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx |        8 |      35 |          43 |         34 |          1 |         1 |          1 |          0 |          0 | 30.1995 | 62.1654 |  62.1654 | -277.939 | 58.2253 | -340.105 |         6 |        16 |        22 |    75      |    45.7143 |    51.1628 |  -10.4 |   -0.6 |            45.9675 |         0.0541608 |
|  2 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxxACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx |        9 |      35 |          44 |         34 |          1 |         1 |          1 |          0 |          0 | 33.7839 | 62.1654 |  62.1654 | -277.939 | 58.2253 | -340.105 |         6 |        16 |        22 |    66.6667 |    45.7143 |    50      |  -10.4 |   -0.6 |            45.9675 |         0.051455  |
|  3 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxCACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx |       10 |      35 |          45 |         34 |          1 |         1 |          1 |          0 |          0 | 38.5141 | 62.1654 |  62.1654 | -277.939 | 58.2253 | -340.105 |         7 |        16 |        23 |    70      |    45.7143 |    51.1111 |  -10.4 |   -0.6 |            45.9675 |         0.0826205 |
|  4 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx |       11 |      35 |          46 |         34 |          1 |         1 |          1 |          0 |          0 | 40.8741 | 62.1654 |  62.1654 | -277.939 | 58.2253 | -340.105 |         7 |        16 |        23 |    63.6364 |    45.7143 |    50      |  -10.4 |   -0.6 |            45.9675 |         0.0910506 |

  

# Installation from source code:
The webtool app can accommodate most applications by choosing the most appropriate model parameters for your experimental conditions. 

For processing large number of pegRNAs, researchers can download zipped source code, install the necessary python packages, and run DeepPrime on their local systems. We recommend using a Linux-based OS.

### 1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### 2. Create and activate virtual environment
```python
conda create -n dprime python=3.8
conda activate dprime
```


### 3. Install Required Python Packages
```
pip install tensorflow==2.8.0     #Use pip linked to the above python installation
pip install torch==1.10.0+cu113 torchvision==0.11.1+cu113 torchaudio===0.10.0+cu113 -f https://download.pytorch.org/whl/cu113/torch_stable.html
pip install biopython==1.78 
pip install pandas regex silence-tensorflow 
```	
### 4. [Install ViennaRNA](https://www.tbi.univie.ac.at/RNA/documentation.html#install)
```
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_5_x/ViennaRNA-2.5.1.tar.gz
tar -zxvf ViennaRNA-2.5.1.tar.gz
cd ViennaRNA-2.5.1
./configure --with-python3	
make
make install

- OR -

conda install -c bioconda viennarna # using Miniconda

- OR -

pip install ViennaRNA
```

### 5. Download Source Code
```
wget https://github.com/hkimlab/DeepPrime/archive/main.zip
unzip main.zip
```

# Usage:

## Input format (.csv file)
ID, Unedited sequences (121 bp), Unedited sequences (121bp), alt_type (sub1, sub2, sub3, ins1, ... , del3)
| ID                     | RefSeq                                                                                                                    | Edited Seq                                                                                                                | EditType |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- | -------- |
| BRCA1e17_pos34_tat_CAT | AATCCTTTGAGTGTTTTTCATTCTGCAGATGCTGAGTTTGTGTGTGAACGGACACTGAAATATTTTCTAGGAATTGCGGGAGGAAAATGGGTAGTTAGCTATTTCTGTAAGTATAATACTA | AATCCTTTGAGTGTTTTTCATTCTGCAGATGCTGAGTTTGTGTGTGAACGGACACTGAAACATTTTCTAGGAATTGCGGGAGGAAAATGGGTAGTTAGCTATTTCTGTAAGTATAATACTA | sub1     |
| BRCA1e17_pos34_tat_CCA | AATCCTTTGAGTGTTTTTCATTCTGCAGATGCTGAGTTTGTGTGTGAACGGACACTGAAATATTTTCTAGGAATTGCGGGAGGAAAATGGGTAGTTAGCTATTTCTGTAAGTATAATACTA | AATCCTTTGAGTGTTTTTCATTCTGCAGATGCTGAGTTTGTGTGTGAACGGACACTGAAACCATTTCTAGGAATTGCGGGAGGAAAATGGGTAGTTAGCTATTTCTGTAAGTATAATACTA | sub3     |
| BRCA1e17_pos34_tat_CCC | AATCCTTTGAGTGTTTTTCATTCTGCAGATGCTGAGTTTGTGTGTGAACGGACACTGAAATATTTTCTAGGAATTGCGGGAGGAAAATGGGTAGTTAGCTATTTCTGTAAGTATAATACTA | AATCCTTTGAGTGTTTTTCATTCTGCAGATGCTGAGTTTGTGTGTGAACGGACACTGAAACCCTTTCTAGGAATTGCGGGAGGAAAATGGGTAGTTAGCTATTTCTGTAAGTATAATACTA | sub3     |

## Run Command:
```md
python DeepPrime.py [-h] [-f INPUT_FILE] [-n NAME] [-p {PE2,PE2max,PE2max-e,PE4max,PE4max-e,NRCH_PE2,NRCH_PE2max,NRCH_PE4max}] [--cell_type {HEK293T,A549,DLD1,HCT116,HeLa,MDA-MB-231,NIH3T3}] [--pbs_min PBS_MIN] [--pbs_max PBS_MAX] [--jobs JOBS] [--progress]
```
Basic command: python DeepPrime.py -f [filename]

```bash
# example_input
python DeepPrime.py -f ./example_input/dp_core_test.csv
```
```bash
# example_input & choose PE4max system
python DeepPrime.py -f ./example_input/dp_core_test.csv -p PE4max
```
```bash
# example_input & choose PE4max system, cell type, and number of cores
python DeepPrime.py -f ./example_input/dp_core_test.csv -p PE2max --cell_type DLD1 --jobs 4
```
## Optional arguments
-h or --help: show a help message  
-f or --input_file: input path (.csv file)  
-n or --name: name tag of run (results directory name)  
-p or --pe_type: PE system. Choose one of the available PE system (PE2,PE2max,PE2max-e,PE4max,PE4max-e,NRCH_PE2,NRCH_PE2max,NRCH_PE4max). Some cell types support limited PE systems.  
--cell_type: Cell type. Choose one of the available cell line.  
--pbs_min: Minimum length of PBS. (1=<)  
--pbs_max: Maximum length of PBS (=<17)  
--jobs: Number of cores for computing  
--progress: Show processing message


	
## Current available PE models:
### On-target
| Cell type  | PE system   | Model                                                             |
| ---------- | ----------- | ----------------------------------------------------------------- |
| HEK293T    | PE2         | DeepPrime_base                                                    |
| HEK293T    | NRCH_PE2    | DeepPrime-FT: HEK293T, NRCH-PE2 with Optimized scaffold           |
| HEK293T    | NRCH_PE2max | DeepPrime-FT: HEK293T, NRCH-PE2max with Optimized scaffold        |
| HEK293T    | PE2         | DeepPrime-FT: HEK293T, PE2 with Conventional scaffold             |
| HEK293T    | PE2max-e    | DeepPrime-FT: HEK293T, PE2max with Optimized scaffold and epegRNA |
| HEK293T    | PE2max      | DeepPrime-FT: HEK293T, PE2max with Optimized scaffold             |
| HEK293T    | PE4max-e    | DeepPrime-FT: HEK293T, PE4max with Optimized scaffold and epegRNA |
| HEK293T    | PE4max      | DeepPrime-FT: HEK293T, PE4max with Optimized scaffold             |
| A549       | PE2max-e    | DeepPrime-FT: A549, PE2max with Optimized scaffold and epegRNA    |
| A549       | PE2max      | DeepPrime-FT: A549, PE2max with Optimized scaffold                |
| A549       | PE4max-e    | DeepPrime-FT: A549, PE4max with Optimized scaffold and epegRNA    |
| A549       | PE4max      | DeepPrime-FT: A549, PE4max with Optimized scaffold                |
| DLD1       | NRCH_PE4max | DeepPrime-FT: DLD1, NRCH-PE4max with Optimized scaffold           |
| DLD1       | PE2max      | DeepPrime-FT: DLD1, PE2max with Optimized scaffold                |
| DLD1       | PE4max      | DeepPrime-FT: DLD1, PE4max with Optimized scaffold                |
| HCT116     | PE2         | DeepPrime-FT: HCT116, PE2 with Optimized scaffold                 |
| HeLa       | PE2max      | DeepPrime-FT: HeLa, PE2max with Optimized scaffold                |
| MDA-MB-231 | PE2         | DeepPrime-FT: MDA-MB-231, PE2 with Optimized scaffold             |
| NIH3T3     | NRCH_PE4max | DeepPrime-FT: NIH3T3, NRCH-PE4max with Optimized scaffold         |

### Off-target (currently writing the manual)
| Cell type  | PE system   | Model                                                             |
| ---------- | ----------- | ----------------------------------------------------------------- |
| HEK293T    | PE2-off     | DeepPrime-Off: PE2 with conventional scaffold in HEK293T cells    |

		
For off-target analysis:
Currently, only the model trained on PE2 with conventional scaffold in HEK293T cells is capable of running an additional analysis to predict off-target levels for specific pegRNAs.

On the webtool:
First select the Off-target compatible, PE2_Conv, and run your inputs. On the results page, use the check box indicating that you are currently running the off-target compatible analysis. Selecting individual rows will auto-fill the pegRNA IDs and the off-target sequences can be added to the text area in 74bp long formats.

On the source code:
Create input file (ex: dp_off_test.csv), and run
```
python DeepPrime.py off_run <filename>

ex)
python DeepPrime.py -f ./example_input/dp_off_test.csv -p PE2-off
```
