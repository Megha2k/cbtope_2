# CBTOPE2: Prediction of interactivity of antigens in protein sequence

This repository contains the standalone code for CBTOPE2 prediction tool. CBTOPE2 is an antibody interacting residues predictor. To support the scientific community, we have developed a standalone software and web server for the tool. For a user-friendly access to CBTOPE2, please use https://webs.iiitd.edu.in/raghava/cbtope2/.

## Installation (dependencies)

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

---

## Installation (package)

9. Install CBTOPE2:
    ```bash
    pip install cbtope-2
    ```

---

## Execute

To run the CBTOPE2 package:
```bash
python3 -m cbtope_2.standalone -i [filename.fasta] -t [probability threshold = 0.5] -m [1,2]
