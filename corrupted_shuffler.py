"""
Plots Theorem 7.1, Proposition 5 and Theorem 7.2 (Figure 6)
Calculates the privacy of the local model (LDP)
when the shuffler is compromised
while the three mechanisms were running with epsilon1 or epsilon2 privacy levels

NOTE: some values for RR-shuffle are missing
because the mechanism cannot achieve the desired privacy for this combination of parameters (explained in Appendix A)

python3 corrupted_shuffler.py [epsilon1] [epsilon2] [delta] 

e.g. for Figure 6:
python3 corrupted_shuffler.py 0.05 0.2 0.001
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
import statistics
sum_elements = 10000 #the number of elements to sum from the SGDL distribution (10000 should be enough for this, but increase for better accuracy)
mp.dps = 500  #depth of mpmath library since we are handling extremely small numbers

assert(len(sys.argv)==4)

epsilon = float(sys.argv[1])
epsilon_2 = float(sys.argv[2])
delta = float(sys.argv[3])

d=1 #distance between datasets, should be set to 1
e_threshold = 0.001 #how close the binary search is ok to get to epsilon to consider that it has found it
k = 1000 #maximum value each user can have
N_MAX = 1000 #maximum number of users
recursion_counter = 0 #program crashes if the recursion depth is above 100


def use_binary_search (target_e,delta,n,e_threshold,MAX_E,MIN_E):
	"""
	Geo-Shuffle:
	To avoid dealing with the complex Theorem 4.2.
	we use binary search to find the correct epsilon that satisfies Theorem 4.2.
	"""
	global recursion_counter
	if (recursion_counter>=100):
		print("Error: Binary search failed...")
		exit()

	recursion_counter+=1
	middle = (MAX_E+MIN_E)/2
	geo_epsilon = middle
	geo_shuffle_epsilon = find_e_of_sgdl(delta,geo_epsilon,n)	
	if (abs(geo_shuffle_epsilon - target_e ) <= e_threshold):
		return geo_epsilon


	if (target_e>geo_shuffle_epsilon):
		return use_binary_search (target_e,delta,n,e_threshold,MAX_E,middle)
	else:
		return use_binary_search (target_e,delta,n,e_threshold,middle,MIN_E)

def find_e_of_sgdl(delta,geo_epsilon,b):
	"""
	SGDL-Shuffle:
	Finds the resulting epsilon of SGDL-Shuffle (Theorem 5.1)
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

def gdl_pmf(m,b, p):
	"""
	Geo-Shuffle:
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

def compute_privacy_sgdl_shuffle(geo_epsilon,b,a):
	"""
	SGDL-Shuffle:
	Finds the privacy of SGDL-Shuffle when the shuffler is corrupted (Theorem 7.2)
	"""
	p =np.exp(-geo_epsilon)

	gdl1= gdl_pmf(a,b,p)
	gdl2 = gdl_pmf(a-d,b,p)#
	ratio1 = mp.mpf(gdl1)/mp.mpf(gdl2)
	if (mp.mpf(ratio1)!= 0 and 1/mp.mpf(ratio1) > mp.mpf(ratio1)):
		ratio1 = 1/mp.mpf(ratio1)
	gdl1= gdl_pmf(-a,b,p)
	gdl2 = gdl_pmf(-a-d,b,p)#
	ratio2 = mp.mpf(gdl1)/mp.mpf(gdl2)
	if (mp.mpf(ratio2)!= 0 and 1/mp.mpf(ratio2) > mp.mpf(ratio2)):
		ratio2 = 1/mp.mpf(ratio2)
	if (mp.mpf(ratio1)>mp.mpf(ratio2)):
		max_ratio = ratio1 
	else:
		max_ratio = ratio2
	max_ratio = ratio2
	epsilon = math.log((max_ratio))/d
	return epsilon


def check_if_possible_rr (max_n, k, epsilon , delta):
	"""
	RR-Shuffle:
	Checks if the desired privacy is possible
	(satisfies restrictions of Theorem 3.1.)
	"""
	if ((max_n*k)<= 14 * math.log10(4/delta)):
		return 0 

	try:
		maximum_possible_epsilon = math.sqrt(  (32*math.log10(4/delta))  /  ((max_n*k)-math.sqrt(2*(max_n*k)*math.log10(2/delta)))) 
		if ( maximum_possible_epsilon  >= epsilon):

			return 0 
	except: 
		return 1
	
	return 1


def find_rr_shuffle_local_epsilon (n, k, target_e,e_threshold,delta):
	"""
	RR-Shuffle:
	Finds the privacy of RR-Shuffle when the shuffler is corrupted (Theorem 7.1)
	"""
	if (check_if_possible_rr (n, k, target_e , delta))==0:
		return None

	l = calculate_lambda (n, target_e,k,e_threshold,delta)
	p = l/(n*k)
	if(p*k <=14 * math.log10(4/delta)):
		return None

	try:
		outside_sqrt =  (32*math.log10(4/delta))  /  ((p*k)-math.sqrt(2*(p*k)*math.log10(2/delta)))
		result =  math.sqrt(  (32*math.log10(4/delta))  /  ((p*k)-math.sqrt(2*(p*k)*math.log10(2/delta)))  )
	except:
		return None
	return result

def calculate_lambda (n, target_e,k,e_threshold,delta):
	"""
	RR-Shuffle:
	Calculates the parameter lambda of RR-shuffle
	"""
	for l in range(1,n*k):
		try:#check if sqrt does not take a negative input
			e0 =  math.sqrt(  (32*math.log10(4/delta))  /  ((l)-math.sqrt(2*(l)*math.log10(2/delta)))  )
		except:
			continue	
		e = e0 
		if (abs(e-target_e)<= e_threshold): 
			return l
	#nothing found...
	#search again with greater theshold
	return calculate_lambda (n, target_e,k,0.001+e_threshold,delta)



results_sgdl = []
results_geo = []
results_sgdl_2 = []
results_geo_2 = []
results_sgdl_3 = []
results_geo_3 = []
results_rr = []
results_rr_2=[]
a_list = []

#Find the privacy (Local DP) when the shuffler gets corrupted
#For RR-Shuffle, Geo-Shuffle and SGDL -Shuffle
#for both epsilons

for n in range(1,N_MAX,5):
	recursion_counter=0
	print("Running for ",n,"/",N_MAX)
	results_sgdl.append(compute_privacy_sgdl_shuffle(epsilon,1/n,0))
	results_geo.append(use_binary_search (epsilon,delta,n,e_threshold,3,0))
	results_rr.append(find_rr_shuffle_local_epsilon (n,k,  epsilon,e_threshold,delta))


	recursion_counter=0
	results_sgdl_2.append(compute_privacy_sgdl_shuffle(epsilon_2,1/n,0))
	results_geo_2.append(use_binary_search (epsilon_2,delta,n,e_threshold,3,0))
	results_rr_2.append(find_rr_shuffle_local_epsilon (n,k,  epsilon_2,e_threshold,delta))

	a_list.append(n)



#make a nice plot
epsilon_latex = r'$\varepsilon$'
p_latex = r'$\mathbb{P}$'
not_equal = r'$\neq$'
beta = r'$\beta$'
epsilon_lambda =  r'$\varepsilon_L$'
epsilon_sigma =  r'$\varepsilon_S$'



rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amssymb}')
plt.rcParams.update({'font.size': 25})
plt.rcParams["figure.figsize"] = (20,12)
plt.figure().gca().xaxis.set_major_locator(MaxNLocator(integer=True))


none_list = []
for i in range(len(a_list)):
	none_list.append(None)



plt.plot(a_list,results_sgdl, color="purple", label="SGDL-Shuffle")
plt.plot(a_list,results_geo, color ="cyan", label="Geo-Shuffle")
plt.plot(a_list,results_rr, color ="red", label="RR-Shuffle")

plt.plot(a_list,results_sgdl_2, '--', color="purple")
plt.plot(a_list,results_geo_2,  '--',color ="cyan")
plt.plot(a_list,results_rr_2,  '--',color ="red")


plt.plot(a_list,none_list,'--', color="black", label='_nolegend_' )
plt.plot(a_list,none_list, color="black",  label='_nolegend_')

ax = plt.gca()
ax.legend_ = None

lines = plt.gca().get_lines()
include = [0,1,2]
legend1 = plt.legend([lines[i] for i in include],[lines[i].get_label() for i in include], loc=1)
legend2 = plt.legend([lines[i] for i in [7,6]],[str(epsilon_sigma) + " = " + str(epsilon),str(epsilon_sigma) + " = " + str(epsilon_2)], loc="upper left")
plt.gca().add_artist(legend1)
plt.gca().add_artist(legend2)

plt.ylabel(epsilon_lambda, fontsize = 35)
plt.xlabel("Number of Users ", fontsize = 35)
name_of_file = "corrupted_shuffler.pdf"
plt.savefig(name_of_file, format="pdf", bbox_inches="tight")

plt.show()