"""
Plots Conjecture B.3
For every combination of parameters the upper bound of the privacy should be given when d=1

python3 plot_conj.py [geo_epsilon] [delta] [n] 

where:
geo_epsilon: the starting epsilon of Geo-Local, before shuffling
delta: delta of metric privacy
n: number of users

e.g. for Figure 10:
python3 plot_conj.py 0.2 0.001 100 

Figure 11:
python3 plot_conj.py 0.8 0.001 100 

Figure 12:
python3 plot_conj.py 0.8 0.001 10000

"""

import scipy
import numpy as np
import matplotlib.pylab as plt
import random
import math
import mpmath as mp
from datetime import *
import sys
import pandas as pd
import itertools
from mpl_toolkits import mplot3d
from matplotlib import rc
from matplotlib.ticker import MaxNLocator


assert (len(sys.argv) == 4)

geo_epsilon = float(sys.argv[1])
delta = float(sys.argv[2])
n = int(sys.argv[3])
D_MAX = 500 #the maximum d to go up to and calculate the ratio, in the figures of thepaper its 500

sum_elements = 10000 #the number of elements to sum from the SGDL distribution 10000 should be enough

mp.dps = 500 #depth of mpmath library since we are handling extremely small numbers


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

def privacy_of_geo_shuffle(r,d):
	"""
	calculates the privacy of Geo-Shuffle: 
	the ratio between two PMF of SGDL
	one on some point r-d/2 and the other one on r+d/2(Theorem 4.1)
	"""
	global geo_epsilon 
	global n 
	b = n 

	p =np.exp(-geo_epsilon)
	gdl1= gdl_pmf(math.floor(r-(d/2)),b,p)
	gdl2 = gdl_pmf(math.floor(r+(d/2)),b,p)#
	ratio1 = mp.mpf(gdl1)/mp.mpf(gdl2)
	if (mp.mpf(ratio1)!= 0 and 1/mp.mpf(ratio1) > mp.mpf(ratio1)):
		ratio1 = 1/mp.mpf(ratio1)
	max_ratio = ratio1
	return max_ratio

#the points in which to calculate the ratio (correspond to alpha of Theorem 4.2)
r1 = -1
r2 = -15
r3 = -30
r4 = -45 
r5 = -60 
r6 = -75

r1_results = []
r2_results = []
r3_results = []
r4_results = []
r5_results = []
r6_results = []
d_list = []

#calculate the privacy of Geo-Shuffle for d from 1 to D_MAX
for d in range (1, D_MAX):
	print("Running for d =",d,"/",D_MAX)
	d_list.append(d)
	r1_results.append(math.log(privacy_of_geo_shuffle(r1,d))/d )
	r2_results.append(math.log(privacy_of_geo_shuffle(r2,d))/d )
	r3_results.append(math.log(privacy_of_geo_shuffle(r3,d))/d )
	r4_results.append(math.log(privacy_of_geo_shuffle(r4,d))/d )
	r5_results.append(math.log(privacy_of_geo_shuffle(r5,d))/d )
	r6_results.append(math.log(privacy_of_geo_shuffle(r6,d))/d )
	


#make a nice plot
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 25})
plt.rcParams["figure.figsize"] = (20,12)

plt.plot(d_list,r1_results, '.', markersize = 2, label= "r =" + str(r1))
plt.plot(d_list,r2_results, '.',markersize = 2,  label= "r =" + str(r2))
plt.plot(d_list,r3_results, '.',markersize = 2,  label= "r =" + str(r3))
plt.plot(d_list,r4_results, '.',markersize = 2, label= "r =" + str(r4))
plt.plot(d_list,r5_results,'.',markersize = 2,  label= "r =" + str(r5))
plt.plot(d_list,r6_results,'.', markersize = 2, label= "r =" + str(r6))

plt.xlabel("d", fontsize = 40)
plt.ylabel("K(r,d)", size= 40 )
plt.legend(loc = "center", bbox_to_anchor=[0.5, 1.05], 
          ncol=3, fancybox=True, shadow=True, fontsize = 30, markerscale = 11)
plt.grid()
name_of_file = "plot_conj.pdf"
plt.savefig(name_of_file, format="pdf", bbox_inches="tight")
