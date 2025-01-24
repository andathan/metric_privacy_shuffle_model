"""
Plots the SGDL distribution 

python3 plot_gdl.py [beta] [p]

e.g. for blue dotted line of Figure 2:
python3 plot_gdl.py 10 0.5 

and the rest lines with:
python3 plot_gdl.py 10 0.8
python3 plot_gdl.py 50 0.5
python3 plot_gdl.py 50 0.8
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
sum_elements = 10000 #the number of elements to sum from the SGDL distribution (10000 should be enough to produce an accurate visualization  but increase for better accuracy)
mp.dps = 500  #depth of mpmath library since we are handling extremely small numbers


R_MAX = 100 #range of results, will show values from -R_MAX up to +R_MAX

assert(len(sys.argv)==3)
b = int(sys.argv[1]) #will correspond to number of users
p = float(sys.argv[2]) #relates to epsilon


def gdl_pmf(m,b, p):
	"""
	calculates the PMF of SGDL on some point m
	"""

	m = math.floor(m)
	assert(p<=1 and p>=0)
	#split it in parts: par1, par2
	#no need for par1 (aka outside the sum), will be omitted on division
	#par1 = (mp.mpf(mp.mpf(1-p)**mp.mpf(b))) * (mp.mpf(mp.mpf(1-p)**mp.mpf(b)))
	par1=1
	par2 = 0
	summand = 0 
	prev_summand = 0 
	for k in range (abs(m), abs(m)+sum_elements):
		probabilities = mp.mpf(p)**(mp.mpf(k))      * mp.mpf(p)**mp.mpf(k-abs(m))
		binomials = mp.binomial(b+k-1,k)*mp.binomial(mp.mpf(b)+mp.mpf(k)-abs(m)-1,mp.mpf(k)-abs(m))
		prev_summand = par2
		par2 +=  mp.mpf(binomials)*mp.mpf(probabilities)
		summand = par2
	return_value = mp.mpf((1-p))**mp.mpf((2*b)) *  mp.mpf(par1)*mp.mpf(par2)
	return return_value



results_b_p = []
r_list = []
for r in range(-R_MAX,R_MAX):
	r_list.append(r)
	results_b_p.append(gdl_pmf(r,b,p))


#create a nice plot:
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 25})
plt.rcParams["figure.figsize"] = (20,12)
plt.figure().gca().xaxis.set_major_locator(MaxNLocator(integer=True))

beta_latex = r'$\beta$'
sgdl_latex = r'SGDL($\beta$,p)'
p_latex = 'P'
y_axis_label =  p_latex + "[y =  " + sgdl_latex + "]"

plt.gca().set_ylabel(y_axis_label, fontsize = 40)
plt.plot(r_list,results_b_p, '.', markersize = 7, label= "" + beta_latex + " = " + str(b) + ", p = " + str(p))

plt.xlabel("y", fontsize = 40)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.grid()
plt.show()
