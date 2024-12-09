# Artifact Appendix

Paper title: **Enhancing Metric Privacy With a Shuffler**

Artifacts HotCRP Id: **9** (2025.2)

Requested Badge: **Available**, **Functional**, **Reproduced**

## Description
This repository includes the code for the experiments of our paper: Enhancing Metric Privacy With a Shuffler

1. exp1 and exp2 folders

   Experiments of Section 6 

2. calculate_boost.py 

   Calculation of the privacy boost of Geo-Shuffle (Section 4.2, Figure 3)

3. corrupted_shuffler.py

   Plots the privacy of the mechanisms when the shuffler is corrupted (Figure 6)

4. plot_conj.py

   Plots the Conjecture B.3 (Figures 9-12)

5. plot_gdl.py (minor script)

   Visualization of the SGDL distribution (Figure 2)

### Security/Privacy Issues and Ethical Concerns (All badges)
No ethical concerns or security/privacy issues. 

## Basic Requirements (Only for Functional and Reproduced badges)
python3 is required 

### Hardware Requirements
No special hardware is required. 

### Software Requirements
python3 is required with the following commonly used libraries:
1. pandas
2. numpy
3. scipy
4. matplotlib
5. mpmath

##### Datasets:
For the first experiment, create_random_data.py creates a (uniformly) random synthetic dataset. The file "synthetic" is an example created by this script.
You don’t need to run the script again; just use the "synthetic" dataset as it is.

For the second experiment, the dataset is stored in austin.csv. The dataset is taken from:
https://data.world/tronovan/austin-address-with-gps-coords. The dataset is already in the exp2 folder; you do not need to download anything.

### Estimated Time and Storage Consumption
No extra storage space is taken upon execution of the experiments. The experiment scripts take only 35MB of space.

##### Estimated Time:
###### Less than a minute
1. plot_gdl.py
2. calculate_boost.py 

###### Approx. 1-2 hours
1. exp1.py
2. exp2.py
3. plot_conj.py

[For plot_conj.py to compute faster (less than 5 minutes), reduce the parameter D_MAX to 20 to produce results only when d=1,2,...,20 (LINE 43)].

###### Approx. 10 hours

1. corrupted_shuffler.py

[For corrupted_shuffler.py to compute faster (approx. 30 minutes), reduce the parameter N_MAX to 100 to produce results only when there are 1,2,...,100 users (LINE 46)].

## Environment 
python3 and the above commonly used libraries should be installed.

### Accessibility (All badges)
https://github.com/andathan/metric_privacy_shuffle_model

TO-DO: Add Tag

### Set up the environment (Only for Functional and Reproduced badges)
On a Linux machine:

```bash
git clone https://github.com/andathan/metric_privacy_shuffle_model && cd metric_privacy_shuffle_model
sudo apt update && sudo apt install -y python3 python3-pip && pip3 install pandas numpy scipy matplotlib mpmath
