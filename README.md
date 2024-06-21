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
conda create -n genet python=3.10
conda activate genet

# Install genet
pip install genet
```

### How to use DeepPrime using GenET
```python
from genet.predict import DeepPrime

seq = 'CCGAGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCACGCTCCATTATC(C/T)AGCCCCAAAGCGCAACAAGCCCACTGTCTATGGTGTGTCCCCCAACTACGACAAGTGGGA'

pegrna = DeepPrime(seq)

# check designed pegRNAs
pegrna.features.head()
```
output:

|   | ID         | Spacer               | RT-PBS                                            | PBS_len | RTT_len | RT-PBS_len | Edit_pos | Edit_len | RHA_len | Target                                            | ... | deltaTm_Tm4-Tm2 | GC_count_PBS | GC_count_RTT | GC_count_RT-PBS | GC_contents_PBS | GC_contents_RTT | GC_contents_RT-PBS | MFE_RT-PBS-polyT | MFE_Spacer | DeepSpCas9_score |
| - | ---------- | -------------------- | ------------------------------------------------- | ------- | ------- | ---------- | -------- | -------- | ------- | ------------------------------------------------- | --- | --------------- | ------------ | ------------ | --------------- | --------------- | --------------- | ------------------ | ---------------- | ---------- | ---------------- |
| 0 | SampleName | GGTTCATCATCATTCAACGG | TAGATAATGGAGCGTGGTGATGAGCCCGTCGGCCACCGTTGAATG     | 7       | 38      | 45         | 37       | 1        | 1       | AGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCAC... | ... | \-510.285       | 2            | 23           | 25              | 28.57143        | 60.52632        | 55.55556           | \-12.7           | 0          | 76.43662         |
| 1 | SampleName | GGTTCATCATCATTCAACGG | TAGATAATGGAGCGTGGTGATGAGCCCGTCGGCCACCGTTGAATGA    | 8       | 38      | 46         | 37       | 1        | 1       | AGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCAC... | ... | \-510.285       | 2            | 23           | 25              | 25              | 60.52632        | 54.34783           | \-11.4           | 0          | 76.43662         |
| 2 | SampleName | GGTTCATCATCATTCAACGG | TAGATAATGGAGCGTGGTGATGAGCCCGTCGGCCACCGTTGAATGAT   | 9       | 38      | 47         | 37       | 1        | 1       | AGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCAC... | ... | \-510.285       | 2            | 23           | 25              | 22.22222        | 60.52632        | 53.19149           | \-11.4           | 0          | 76.43662         |
| 3 | SampleName | GGTTCATCATCATTCAACGG | TAGATAATGGAGCGTGGTGATGAGCCCGTCGGCCACCGTTGAATGATG  | 10      | 38      | 48         | 37       | 1        | 1       | AGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCAC... | ... | \-510.285       | 3            | 23           | 26              | 30              | 60.52632        | 54.16667           | \-11.2           | 0          | 76.43662         |
| 4 | SampleName | GGTTCATCATCATTCAACGG | TAGATAATGGAGCGTGGTGATGAGCCCGTCGGCCACCGTTGAATGATGA | 11      | 38      | 49         | 37       | 1        | 1       | AGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCAC... | ... | \-510.285       | 3            | 23           | 26              | 27.27273        | 60.52632        | 53.06122           | \-11.2           | 0          | 76.43662         |


  

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
