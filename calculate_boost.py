"""
Visualizes the boost of Geo-Shuffle compared to the standard Geometric Mechanism

python3 calculate_boost.py 

e.g. for Figure 3:
python3 calculate_boost.py 
"""


import sys
import scipy
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import mpmath as mp
from datetime import *
import sys
import pandas as pd
import itertools
from matplotlib import rc
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl



mpl.rcParams.update(mpl.rcParamsDefault)
sum_elements = 10000 #the number of elements to sum from the SGDL distribution (10000 should be enough to produce an accurate visualization, but increase for better accuracy)
mp.dps = 500  #depth of mpmath library since we are handling extremely small numbers
d = 1
N_MAX = 500 #maximum number of users

e1 = 0.1
e2 = 0.2
e3 = 0.3
e4 = 0.4
e5 = 0.8

delta = 0.001


def geometric (epsilon,value,max_value,c):
	"""
	standard implementation of the geometric mechanism
	"""
	p = 1 - np.exp(-epsilon)
	Z = np.random.geometric(p) - np.random.geometric(p)
	Z_prime = value + Z
	
	if (c!=-1):
		if (Z_prime > max_value+c):
			return max_value+c
		elif (Z_prime < -c):
			return 0
		else:
			return Z_prime + c
	else:
		return Z_prime



def gdl_pmf(m,b, p): 
	"""
	calculates the PMF of SGDL on some point m
	"""
	global sum_elements
	assert(p<=1 and p>=0)
	par1 = (mp.mpf(mp.mpf(1-p)**mp.mpf(b))) * (mp.mpf(mp.mpf(1-p)**mp.mpf(b)))
	par2 = 0
	summand = 0 
	prev_summand = 0 
	for k in range (abs(m), abs(m)+sum_elements):
		probabilities = mp.mpf(p)**(mp.mpf(k)) * mp.mpf(p)**mp.mpf(k-abs(m))
		binomials = mp.binomial(b+k-1,k)*mp.binomial(mp.mpf(b)+mp.mpf(k)-abs(m)-1,mp.mpf(k)-abs(m))
		prev_summand = par2
		par2 +=  mp.mpf(binomials)*mp.mpf(probabilities)
		summand = par2
	return_value = mp.mpf(par1)*mp.mpf(par2)
	return return_value



def find_e_of_sgdl(delta,geo_epsilon,b):
	"""
	Finds the resulting epsilon of Geo-Shuffle (Theorem 4.1)
	"""
	t = geo_epsilon/(math.sqrt(b))
	p =np.exp(-geo_epsilon)
	try:
		moment_generating_function = ((mp.mpf(1-p)*mp.mpf(1-p))  / ( mp.mpf(1-p*mp.mpf(math.e**mp.mpf(t))  ) *mp.mpf(1-p*mp.mpf(math.e**mp.mpf(-t))  ) )) ** mp.mpf(b)
	except:
		return geo_epsilon

	frac = (delta/4) * mp.mpf(1/mp.mpf(moment_generating_function))
	if (frac < 0):
		return geo_epsilon

	a = int(-mp.log(frac)/t)
	gdl1= gdl_pmf(-a,b,p)
	gdl2 = gdl_pmf(-a-d,b,p)#
	ratio2 = mp.mpf(gdl1)/mp.mpf(gdl2)
	if (mp.mpf(ratio2)!= 0 and 1/mp.mpf(ratio2) > mp.mpf(ratio2)):
		ratio2 = 1/mp.mpf(ratio2)
	max_ratio = ratio2
	epsilon = math.log(max_ratio)/d
	return epsilon






geo_shuffle_e1 = []
geo_shuffle_e2 = []
geo_shuffle_e3 = []
geo_shuffle_e4 = []
geo_shuffle_e5 = []
n_list = []
n_iter_list = []


#decide which users to run on:
for i in range (1, 10, 1):
	n_iter_list.append(i)
for i in range (10, 50, 5):
	n_iter_list.append(i)
for i in range (50, 100, 10):
	n_iter_list.append(i)
for i in range (100, N_MAX, 50):
	n_iter_list.append(i)


#calculate the resulting privacy (aka boost)
for n in n_iter_list:
	n_list.append(n)
	e1_boost = find_e_of_sgdl(delta,e1,n)
	e2_boost = find_e_of_sgdl(delta,e2,n)
	e3_boost = find_e_of_sgdl(delta,e3,n)
	e4_boost = find_e_of_sgdl(delta,e4,n)
	e5_boost = find_e_of_sgdl(delta,e5,n)
	geo_shuffle_e1.append(e1_boost)
	geo_shuffle_e2.append(e2_boost)
	geo_shuffle_e3.append(e3_boost)
	geo_shuffle_e4.append(e4_boost)
	geo_shuffle_e5.append(e5_boost)



#create a nice plot

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 25})
plt.rcParams["figure.figsize"] = (20,12)


plt.figure().gca().xaxis.set_major_locator(MaxNLocator(integer=True))
epsilon_geo_latex = r'$\varepsilon_{geo}$'
epsilon_latex = r'$\varepsilon$'

plt.gca().set_ylabel(r'$\varepsilon$', fontsize = 42)
plt.ylabel(epsilon_latex)
plt.xlabel("Number of users", fontsize = 40)
plt.plot(n_list, geo_shuffle_e1, label = epsilon_geo_latex + " = " + str(e1))
plt.plot(n_list, geo_shuffle_e2, label = epsilon_geo_latex + " = " + str(e2))
plt.plot(n_list, geo_shuffle_e3, label = epsilon_geo_latex + " = " + str(e3))
plt.plot(n_list, geo_shuffle_e4, label = epsilon_geo_latex + " = " + str(e4))
plt.plot(n_list, geo_shuffle_e5, label = epsilon_geo_latex + " = " + str(e5))


plt.legend()
name_of_file = "calculate_boost.pdf"
plt.savefig(name_of_file, format="pdf", bbox_inches="tight")

