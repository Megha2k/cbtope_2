# CBTOPE2: Prediction of interactivity of antigens in protein sequence

This repository contains the standalone code for CBTOPE2 prediction tool. CBTOPE2 is an antibody interacting residues predictor. To support the scientific community, we have developed a standalone software and web server for the tool. For a user-friendly access to CBTOPE2, please use https://webs.iiitd.edu.in/raghava/cbtope2/.

## Installation (dependencies)

## 1. Standalone

Download standalone zip file from https://webs.iiitd.edu.in/raghava/cbtope2/standalone.html

1. Unzip the zip file you just downloaded
   
2. ```bash
   pip install Bio joblib sklearn openpyxl
   ```
   
3. **Execute**:
   ```bash
   python standalone.py [-i INPUT] [filename.fasta] [-t probability threshold = 0.5] [-m {1,2}]
   ```

   optional arguments:
      -i INPUT, --input INPUT
                           Input: protein or peptide sequence(s) in FASTA format
                           or single sequence per line in single letter code
     *Input File:* It allow users to provide input in the FASTA format.

     -t THRESHOLD, --threshold THRESHOLD
                           Threshold: Value between 0 to 1, most preferably 0.5
     *Threshold:* User should provide threshold between 0 and 1, by default its 0.5.

     -m {1,2}, --model {1,2}
                        Model Type: Model Type: 1: PSSM based, 2: RSA + PSSM ensemble model (Best Model).
     *Model *: User is allowed to choose between two different models, such as 1 for PSSM based and 2 for RSA + PSSM ensemble 2.

### OR

## 2. PIP package

To install the dependency - SPOT-1D-Single, use the following commands **in the same order:**

1. Clone the repository:
    ```bash
    git clone https://github.com/jas-preet/SPOT-1D-Single.git
    ```

2. Move into the directory:
    ```bash
    cd SPOT-1D-Single
    ```

3. Download the model checkpoints:
    ```bash
    wget https://apisz.sparks-lab.org:8443/downloads/Resource/Protein/2_Protein_local_structure_prediction/jits.tar.xz
    ```

4. Extract the files:
    ```bash
    tar -xvf jits.tar.xz
    ```

5. Create a conda environment:
    ```bash
    conda create -n cbtope2_env python=3.7
    ```

6. Activate the environment:
    ```bash
    conda activate cbtope2_env
    ```

7. Install PyTorch:

    - **If using GPU:**
      ```bash
      conda install pytorch==1.6.0 torchvision==0.7.0 cudatoolkit=10.1 -c pytorch
      ```

    - **For CPU only:**
      ```bash
      conda install pytorch==1.6.0 torchvision==0.7.0 cpuonly -c pytorch
      ```

8. Install pandas:
    ```bash
    conda install pandas=1.1.1
    ```

*(Refer to https://github.com/jas-preet/SPOT-1D-Single/tree/master for details.)*

9. Install PSSM dependency:
    You can install system-specific ncbi psi-blast files from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

10. Install CBTOPE2 package:
    ```bash
    pip install cbtope-2
    ```
