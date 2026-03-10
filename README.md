# VIP2B
VIrome Profiler with type IIB restriction sites for WMS data

## Installation
 
 ### System requirements
 
 #### Dependencies
All scripts in VIP2B are programmed by Perl and Python, and execution of VIP2B is recommended in a conda environment. This program could work properly in the Unix systems, or Mac OSX, as all required packages can be appropreiately download and installed.  
 #### Memory usage
**> 2G RAM** is required to run this pipeline.  
 ### Download the pipeline
 * Clone the latest version from GitHub (recommended):  
 
   `git clone https://github.com/sunzhengCDNM/VIP2B/`  
   `cd VIP2B`
   
    This makes it easy to update the software in the future using `git pull` as bugs are fixed and features are added.
 * Alternatively, directly download the whole GitHub repo without installing GitHub:
 
   `wget https://github.com/sunzhengCDNM/VIP2B/archive/master.zip`  
   `unzip master.zip`
   `cd VIP2B-master`
   
 ### Install VIP2B in a conda environment 
 * Conda installation  
   [Miniconda](https://docs.conda.io/en/latest/miniconda.html) provides the conda environment and package manager, and is the recommended way to install VIP2B. 
 * Create a conda environment for VIP2B pipeline:  
   After installing Miniconda and opening a new terminal, make sure you’re running the latest version of conda:
   
   `conda update conda`
   
   Once you have conda installed, create a conda environment with the yml file `config/conda-20250709.yml`.
   
   `conda env create -n VIP2B.1.1 --file config/conda-20250709.yml`
   
 * Activate the VIP2B conda environment by running the following command:
 
   `conda activate VIP2B.1.1` or `source activate VIP2B.1.1`
   
   Make sure the conda environment of VIP2B has been activated by running the above command before you run VIP2B everytime.  

   Now, everything is ready for VIP2B :), Let's get started.
 
## Using VIP2B
 
### Quick start
VIP2B is a highly automatic pipeline, and only a few parameters are required for the pipeline.
* We prepared a test pair-end sequencing data:  
 
   `cd example`
   `mkdir -p data/`  
   `wget -t 3 -O data/test_seq.R1.fq.gz https://figshare.com/ndownloader/files/52717946`  
   `wget -t 3 -O data/test_seq.R2.fq.gz https://figshare.com/ndownloader/files/52717949`
 
* After downloading the sequencing data, we can finally run VIP2B:  
 
   `python3 ../bin/VIP2B.py -i data.list`

    In `data.list` you can learn how to prepare your input data, both single-end and paired-end data can be used as input.  
    
```
sample1 <tab> shotgun1_left.fastq(.gz) <tab> shotgun1_right.fastq(.gz)
sample2 <tab> shotgun2.fastq(.gz)
sample3 ...
```

### Parameters
The main program is `bin/VIP2B.py` in this repo. You can check out the usage by printing the help information via `python3 bin/VIP2B.py -h`.

```
usage: VIP2B.py [-h] -i INPUT [-o OUTPUT]
                [-l {Class,Order,Family,Genus,Species}] [-e ENZYME]
                [-d DATABASE] [-p PROCESSES] [-t THRESHOLD] [-c CUTOFF]
                [--intersection]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT              The filepath of the sample list. Each line includes an input sample ID and the file path of corresponding DNA sequence data where each field should be separated by <tab>. A line in this file that begins with # will be ignored. like 
                          sample <tab> shotgun.1.fq(.gz) (<tab> shotgun.2.fq.gz)
  -o OUTPUT             Output directory, default ./VIP2B_result
  -l {Class,Order,Family,Genus,Species}
                        taxo level, choose from Class/Order/Family/Genus/Species, default Species
  -e ENZYME             Enzyme, One or more of the following enzymes: AlfI/AloI/BaeI/BcgI/BplI/BsaXI/BslFI/Bsp24I/CjeI/CjePI/CspCI/FalI/HaeIV/Hin4I/PpiI/PsrI. Multiple enzymes should be separated by commas. If you want to select all enzymes, you can simply provide "all". It should be noted that the selection of enzymes needs to correspond with the database. The default combination is AlfI,BcgI,BslFI,CjeI,CjePI,FalI,HaeIV,Hin4I.
  -d DATABASE           Database for pipeline, default /data/USER/liujiang/tmp/sunzheng/20240905_viral/VIP2B.1.1/database/8Enzyme
  -p PROCESSES          Number of processes, note that more threads may require more memory, default 1
  -t THRESHOLD          Threshold for species identification, G5 means using gscore > 5 and M0.5 means using ML probability > 0.5 as a filtering parameter, G2/G5/M0.1/M0.5 are a few commonly used options, default M0.5
  -c CUTOFF             cut off for database, default 30000
  --intersection        intersection or union of tags between genomes, default union

author: Zheng Sun, Jiang Liu
mail: sunzheng0618@gmail.com
date: 2026/3/6 22:02:47
version:  1.1

```
