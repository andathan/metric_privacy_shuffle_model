# Artifact Appendix

Paper title: **Enhancing Metric Privacy With a Shuffler**

Artifacts HotCRP Id: **9** (2025.2)

Requested Badge: **Available**, **Functional**, **Reproduced**

## Description
This repository includes the code for the experiments of our paper: Enhancing Metric Privacy With a Shuffler


1. experiments 

	Experiments of Section 6 

2. calculate_boost.py 

    Calculation of the privacy boost of Geo-Shuffle (Section 4.2. | Figure 3)

3. corrupted_shuffler.py

	Plots the privacy of the mechanisms when the shuffler is corrupted (Figure 6)

4. plot_conj.py

	Plots the Conjecture B.3 (Figures 9-12)

5. 	plot_gdl.py

	Visualization of the SGDL distribution  (Figure 2)



### Security/Privacy Issues and Ethical Concerns (All badges)
No ethical concerns or security/privacy issues. 

## Basic Requirements (Only for Functional and Reproduced badges)
python3 is required 

### Hardware Requirements
No special hardware is required. 

### Software Requirements
python3 is required with the following libraries:
1.	pandas
2.	numpy
3.	scipy
4.	matplotlib
5.	mpmath

##### Datasets:
For the first experiment, create_random_data.py creates a (uniformly) random synthetic dataset. The file "synthetic" is an example, created by this script.
You dont need to run the script again, just use the "synthetic" dataset as it is.

For the second experiment, the dataset is stored in austin.csv. The dataset is taken from:
https://data.world/tronovan/austin-address-with-gps-coords



### Estimated Time and Storage Consumption
No extra storage space is taken upon execution of the experiments. The experiment scripts take only 30MB of space.


##### Estimated Time:
###### Less than a minute
1. plot_gdl.py
2. calculate_boost.py 

###### Aprox. 1-2 hours
1. exp1.py
2. exp2.py
3.  plot_conj.py
[To compute faster (less than 5 minutes) reduce the parameter D_MAX to 20 to produce results only when d=1,2,...,20 (LINE 43)]


###### Aprox. 10 hours

2. corrupted_shuffler.py
[To compute faster (approx 30 minutes) reduce the parameter N_MAX to 100 to produce results only when there are 1,2...,100 users  (LINE 46)]

2. Experiments of Section 6 



## Environment 
python3 and the above commonly used libraries should be installed.


### Accessibility (All badges)
https://github.com/andathan/metric_privacy_shuffle_model

TO-DO: Add Tag

### Set up the environment (Only for Functional and Reproduced badges)
On a Linux machine:

TO-DO: ADD GIT CLONE
```bash 
sudo apt update && sudo apt install -y python3 python3-pip && pip3 install pandas numpy scipy matplotlib mpmath
```


### Testing the Environment (Only for Functional and Reproduced badges)
First try running the scripts which do not take much time:
```bash
python3 plot_gdl.py
python3 calculate_boost.py 
```

If it runs smoothly, you are good to go!

## Artifact Evaluation (Only for Functional and Reproduced badges)

### Main Results and Claims

#### Main Result 1: SGDL-Shuffle offers the best utility (MAE), followed by Geo-Shuffle and then by RR-Shuffle. 


#### Main Result 2: Geo-Shuffle amplifies the privacy of Geo-Local (Figure 3)
Consider the Geometric Mechanism in the Local Model (Geo-Local) with epsilon = epsilon_geo. Adding a shuffler to this mechanism (Section 4) results in better privacy (formally calculated in Theorem 4.2.), shown in Figure 3.

#### Main Result 3: The privacy of SGDL-Shuffle diminishes when the shuffler is compromised; Geo-Shuffle still maintains a reasonable privacy.

Consider the Geometric Mechanism in the Local Model (Geo-Local) with epsilon = epsilon_geo. Adding a shuffler to this mechanism (Section 4) results in better privacy (formally calculated in Theorem 4.2.), shown in Figure 3.


#### Main Result 4: Validity of Conjecture B.3: d=1 maximizes the ratio (4)
The complex formula of SGDL makes an analytical solution difficult. However we experimentally show that d=1 maximizes ratio (4) with every combination of parameters (Figures 9-12).

This is intuitive as setting the number of users to 1 results to the worst privacy amplification of Geo-Shuffle's (shuffling only 1 user is pointless). Still, in this worst case, Geo-Shuffle's privacy can be described by Geo-Local (the standard Geometric Mechanism) where it is known that this ratio is independent of d. 


### Experiments 
List each experiment the reviewer has to execute. Describe:
 - How to execute it in detailed steps.
 - What the expected result is.
 - How long it takes and how much space it consumes on disk. (approximately)
 - Which claim and results does it support, and how.

#### Experiment 1: Synthetic data experiment (Section 6.2., Figure 4)

```bash
python3 exp1/exp1.py confg 0.1 0.001
```

Takes approximately 1-2 hours to run.

The experiment tests the utility (MAE) of the proposed mechanisms on a synthetic dataset.

The parameters of the experiment are loaded from the confg file. 

Run the above command and compare the result to Figure 4. 

We use this experiment for Main Result 1, i.e. that SGDL-Shuffle has lower utility loss than every other mechanism, followed by Geo-Shuffle and then RR-shuffle (when n > aprox. 60). 

#### Experiment 2: Real world location experiment (Section 6.3., Figure 5)

```bash
python3 exp2/exp2.py 0.15 4000 1000 600
```

Takes approximately 1-2 hours to run.

The experiment tests the utility (MAE) of the proposed mechanisms on real-world location data.

The parameters of the experiment are defined in parameters.py.

Run the above command and compare the result to Figure 5. 

We use this experiment for Main Result 1, i.e. that SGDL-Shuffle has lower utility loss than every other mechanism, followed by Geo-Shuffle and then RR-shuffle (when n > aprox. 60). 



#### Experiment 3: Geo-Shuffle amplifies the privacy of Geo-Local (Figure 3) 
Figure 3 shows the final epsilon of Geo-Shuffle when the starting epsilon_geo (of Geo-Local) is set. 
For the blue line run:
```bash
python3 calculate_boost.py 0.1 0.001 
```
and for the rest lines:
```bash
python3 calculate_boost.py 0.2 0.001 
python3 calculate_boost.py 0.3 0.001 
python3 calculate_boost.py 0.4 0.001 
python3 calculate_boost.py 0.8 0.001 
```

Compare the results with Figure 3.

We use this experiment to conclude in Main Result 2 i.e. that the resulting epsilon is notably smaller than the initial one (epsilon_geo), in all cases.

#### Experiment 4: Validity of Conjecture B.3: d=1 maximizes the ratio (4)
Figures 10-12 experimentally show the validity of conjecture for any combination of parameters. In each figure, the script outputs the result of the ratio (4), which we name as the function K(r,d) as d increases. Multiple values are shown for r (an integer which corresponds to alpha of Theorem 4.2.).

Figure 10:
```bash
python3 plot_conj.py 0.2 0.001 100
```

Figure 11:
```bash
python3 plot_conj.py 0.8 0.001 100 
```

Figure 12:
```bash
python3 plot_conj.py 0.8 0.001 10000
```

Observe that in all cases d=1 maximizes the ratio. 


## Limitations (Only for Functional and Reproduced badges)
All scripts to produce the  main-body results are included.
Some scripts for the Figures of the Appendix are excluded as they are minor scripts. These are Figure 7, Figure 8, Figure 13, Figure 14 and Figure 22. 

The extra experiments in the Appendix can be recreated by properly adjusting the confg for exp1 and "parameters.py" for exp2. But they do not provide any new interesting results.

## Notes on Reusability (Only for Functional and Reproduced badges)
To increase the accuracy of the experiments, increase the parameters MAX_RUNS (specifies the number of times to run the experiment) and sum_elements (how many elements from the SGDL distribution are sampled).

You can use your own dataset on exp1 by changing the line "databases" in confg.

You can use your own dataset on exp2 by changing the line "db_name" in "parameters.py".

Note that n is the number of users, each having an integer value from 0 to k. 


The code is licensed under GNU General Public License v3.0.


