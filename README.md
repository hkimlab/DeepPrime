# Table of Contents
- About
- Installation
- Usage


## About:

DeepPrime is a deep-learning-based prime editing efficiency prediction tool developed in Laboratory of Genome Editing, Yonsei University. 

It greatly expands upon the previous start-of-the-art pegRNA activity prediction model, DeepPE, which was limited to a specific set of edit type and length combinations.

DeepPrime is developed to predict efficiencies of a nearly all feasible combinations of pegRNA designs. We integrated CNN & RNN to extract inter-sequence features between target DNA and corresponding pegRNA. DeepPrime was trained using 259K pegRNAs with PBS lengths ranging from 1 to 17, RT lengths ranging from 1 to 50, Edit positions ranging from 1 to 30, and editing lengths ranging from 1 to 3.


## Installation:

The webtool app can accommodate most applications by choosing the most appropriate model parameters for your experimental conditions. 

For processing large number of pegRNAs, researchers can download zipped source code, install the necessary python packages, and run DeepPrime on their local systems. We recommend using a Linux-based OS.


### Linux (CentOS/Ubuntu) commands:

	- Install Python OR Miniconda (https://docs.conda.io/en/latest/miniconda.html)
		
		- Python -
		
		wget https://www.python.org/ftp/python/3.8.12/Python-3.8.12.tgz
		tar xzf Python-3.8.12.tgz
		./configure --enable-optimizations
		make altinstall

		- Miniconda -
		
		wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh
		bash Miniconda3-py38_4.12.0-Linux-x86_64.sh

	- Install Required Python Packages -

		pip install tensorflow==2.8.0     #Use pip linked to the above python installation
		pip install torch==1.10.0+cu113 torchvision==0.11.1+cu113 torchaudio===0.10.0+cu113 -f https://download.pytorch.org/whl/cu113/torch_stable.html
		pip install biopython==1.78 
		pip install pandas regex 
		
	- Install ViennaRNA https://www.tbi.univie.ac.at/RNA/documentation.html#install -

		wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_5_x/ViennaRNA-2.5.1.tar.gz
		tar -zxvf ViennaRNA-2.5.1.tar.gz
		cd ViennaRNA-2.5.1
		./configure --with-python3	
		make
		make install

		- OR -
		
		conda install -c bioconda viennarna  *using Miniconda
		
		- OR -
		 pip install ViennaRNA

	- Download Source Code -
		wget http://deepcrispr.info/DeepPrime/src/DeepPrime_src_wExamples.zip
		unzip DeepPrime_src_wExamples.zip


## Usage:
	The main script uses the above directory map:
	The data directory functions as the main I/O path, with filename corresponding to designated analysis/experiment tags.

	Run Command:

	Input WT/Edited Sequences:
        python main_src.py main_run <filename>

        ex)
        python main_src.py main_run Analysis_Example

        Input:
            seq.txt 	  --> Unedited sequences (121 bp), Unedited sequences (121-124 bp)
            options.txt   --> <PE model>, <PBS min>, <PBS max>, <RTT max>, <Edit type, length>

            ex)
            seq.txt 	  --> ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT,ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT
            options.txt   --> PE2,1,17,40,sub1 #single line

        Output:
            pegRNA#.csv 			#all pegRNA designs for each input in seq.txt
            pegRNA#_top_designs.csv #top 10 pegRNAs per input with spacer and extension oligo sequences


    Input ClinVar ID:
        python main_src.py clinvar_check  <filename>;

        ex)
        python main_src.py clinvar_check Analysis_Example

        Input:
            clinvar_run.txt   --> <ClinVar ID>, <model or therapy>  #ClinVar ID can be VCV ID or variant ID

        Output:
            seq.txt 	  --> Unedited sequences (121 bp), Unedited sequences (121-124 bp)
            options.txt   --> <PE model>, <PBS min>, <PBS max>, <RTT max>, <Edit type, length>

    After clinvar_check: input for main_run is automatically generated.

    The main script uses the following directory map:
    The data directory functions as the main I/O path, with filename corresponding to designated analysis or experiment tags.
	
   Current available PE models:
        ----------On-target----------
        PE2		                Baseline: PE2 with conventional scaffold in HEK293T cells

        ---Fine-Tuned Models---
        PE2_Conv 	              		PE2 with conventional scaffold in HEK293T cells
        NRCH_PE2_HEK293T	     	    NRCH_PE2 with conventional scaffold in HEK293T cells

        PE2_Opti_HCT116                 PE2 with optimized scaffold in HCT116 cells
        PE2_Opti_MDA                    PE2 with optimized scaffold in MDA cells

        PE2max_Opti_HEK239T             PE2max with optimized scaffold in HEK293T cells
        PE2max_Opti_HeLa                PE2max with optimized scaffold in HeLa cells
        PE2max_Opti_A549                PE2max with optimized scaffold in A549 cells
        PE2max_Opti_DLD1                PE2max with optimized scaffold in DLD1 cells

        NRCH-PE2max_Opti_HEK293T        NRCH-PE2max with optimized scaffold in HEK293T cells

        PE4max_Opti_HEK293T             PE4max with optimized scaffold in HEK293T cells
        PE4max_Opti_A549                PE4max with optimized scaffold in A549 cells
        PE4max_Opti_DLD1                PE4max with optimized scaffold in DLD1 cells

        NRCH-PE4max_Opti_DLD1           NRCH-PE4max with optimized scaffold in DLD1 cells
        NRCH-PE4max_Opti_NIH            NRCH-PE4max with optimized scaffold in NIH cells

        PE2max_epegRNA_Opti_HEK293T     PE2max combined epegRNA with optimized scaffold in HEK293T cells
        PE2max_epegRNA_Opti_A549        PE2max combined epegRNA with optimized scaffold in A549 cells

        PE4max_epegRNA_Opti_HEK293T     PE4max combined epegRNA with optimized scaffold in HEK293T cells
        PE4max_epegRNA_Opti_A549        PE4max combined epegRNA with optimized scaffold in A549 cells


        ----------Off-target compatible----------
        PE2_Conv		        		PE2 with conventional scaffold in HEK293T cells


	working dir
	|---data
		|---Analysis_Example
			|---seqs.txt
			|---options.txt
			|---clinvarrun.txt
		|--output

	|----models
		|---DeepPrime_base
		|---DeepPrime_off
		|---DeepPrime_var
		|---DeepSpCas9

	|----src
		|---biofeat.py
		|---dprime.py
		|---dspcas9py
		|---utils.py
		
		
For off-target analysis:
Currently, only the model trained on PE2 with conventional scaffold in HEK293T cells is capable of running an additional analysis to predict off-target levels for specific pegRNAs.

On the webtool:
First select the Off-target compatible, PE2_Conv, and run your inputs. On the results page, use the check box indicating that you are currently running the off-target compatible analysis. Selecting individual rows will auto-fill the pegRNA IDs and the off-target sequences can be added to the text area in 74bp long formats.

On the source code:
Create two input files, offseq.txt and pegRNA.txt and run

    
    python main_src.py off_run <filename>

    ex)
    python main_src.py off_run Analysis_Example

