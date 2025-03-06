# cbtope_2

System Requirments/Dependencies:

refer to https://github.com/jas-preet/SPOT-1D-Single/tree/master

Hardware Requirments: SPOT-1D-single predictor requires only a standard computer with approximately 16 GB RAM to support the in-memory operations.

Python3.7
Anaconda
CUDA 10.1 (Optional if using GPU)
cuDNN (>= 7.4.1) (Optional if using GPU)

Installation:

To install SPOT-1D-Single and it's dependencies following commands can be used in terminal:

git clone https://github.com/jas-preet/SPOT-1D-Single.git
cd SPOT-1D-Single
To download the model check points from the dropbox use the following commands in the terminal:

wget https://apisz.sparks-lab.org:8443/downloads/Resource/Protein/2_Protein_local_structure_prediction/jits.tar.xz
tar -xvf jits.tar.xz
To install the dependencies and create a conda environment use the following commands

conda create -n spot_single_env python=3.7
conda activate spot_single_env
if GPU computer: 7. conda install pytorch==1.6.0 torchvision==0.7.0 cudatoolkit=10.1 -c pytorch

for CPU only 7. conda install pytorch==1.6.0 torchvision==0.7.0 cpuonly -c pytorch
conda install pandas=1.1.1

Execute:

To run CBTOPE2 standalone use the following command:

python3 standalone.py -i input.fasta -t [threshold probability] -m [1,2,3]
