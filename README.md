# CBTOPE2: Prediction of interactivity of antigens in protein sequence

This repository contains the standalone code for CBTOPE2 prediction tool. CBTOPE2 is an antibody interacting residues predictor. To support the scientific community, we have developed a standalone software and web server for the tool. For a user-friendly access to CBTOPE2, please use https://webs.iiitd.edu.in/raghava/cbtope2/.

## Installation (dependencies)

To install dependency - SPOT-1D-Single, use following command in the same order:
1. ```bash
   git clone https://github.com/jas-preet/SPOT-1D-Single.git
2. ```bash
   cd SPOT-1D-Single

To download the model check points:  
3. ```bash
   wget https://apisz.sparks-lab.org:8443/downloads/Resource/Protein/2_Protein_local_structure_prediction/jits.tar.xz  
4. ```bash
   tar -xvf jits.tar.xz

To install dependencies and create conda environment:  
5. ```bash
   conda create -n cbtope2_env python=3.7  
6. ```bash
   conda activate cbtope2_env

7. 
  if GPU computer: ```bash
     conda install pytorch==1.6.0 torchvision==0.7.0 cudatoolkit=10.1 -c pytorch
  for CPU only 7: ```bash
     conda install pytorch==1.6.0 torchvision==0.7.0 cpuonly -c pytorch

8. ```bash
   conda install pandas=1.1.1

(refer to https://github.com/jas-preet/SPOT-1D-Single/tree/master)

## Installation (package)

9. ```bash
   pip install cbtope-2

## Execute

To install the CBTOPE2 package:
```bash
   python3 -m cbtope_2.standalone -i [filename.fasta] -t [probability threshold = 0.5] -m [1,2]

   option '-m'==1: for PSSM based model
       '-m'==2" for PSSM + RSA based model
    

