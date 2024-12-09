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


from parameters import sum_elements
d=1 


def geometric (epsilon,value,max_value,c):
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


def use_binary_search (target_e,delta,n,e_threshold,MAX_E,MIN_E):
	global sum_elements
	middle = (MAX_E+MIN_E)/2
	if (MAX_E == MIN_E):
		print("MAX_E = MIN_E. Binary search failed...")
		exit()
	geo_epsilon = middle


	geo_shuffle_epsilon = find_e_of_sgdl(delta,geo_epsilon,n)	
	if (geo_shuffle_epsilon>geo_epsilon):
		#print("Error: Needs more summands...")
		sum_elements = 2*sum_elements
		return use_binary_search (target_e,delta,n,e_threshold,MAX_E,MIN_E)

	
	if (abs(geo_shuffle_epsilon - target_e ) <= e_threshold):
		
		return geo_epsilon


	if (target_e>geo_shuffle_epsilon):
		return use_binary_search (target_e,delta,n,e_threshold,MAX_E,middle)
	else:
		return use_binary_search (target_e,delta,n,e_threshold,middle,MIN_E)

	print("Binary search failed...")
	exit()


def find_geo_shuffle_epsilon(target_e,delta,n,e_threshold):
	if (n<=10):
		return target_e
	else:#a lot of users
		MAX_E = 5
		MIN_E = target_e
		return use_binary_search(target_e,delta,n,e_threshold,MAX_E,MIN_E)

def find_c(delta,epsilon,n):
	c = int(- math.log(	delta*(1+math.exp(-epsilon))/(4*n))/epsilon)
	#print ("c is ",c)
	return c

def debias_geo_shuffle (sum_of_all_bits,n,k,c):
	debiased = sum_of_all_bits - (n * c)
	if (debiased<0):
		return 0
	elif (debiased>n*k):
		return n*k
	else:
		return debiased


def gdl_pmf(m,b, p): #symmetric version P[Y=m], GDL with parameters b and 
	global sum_elements
	#print("GDL PMF CALLED WITH PARAM, m =",m,"b =",b,"and p=",p)
	assert(p<=1 and p>=0)
	#split it in parts: par1, par2
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
	t = 2* geo_epsilon/(math.sqrt(b))
	p =np.exp(-geo_epsilon)
	moment_generating_function = ((mp.mpf(1-p)*mp.mpf(1-p))  / ( mp.mpf(1-p*mp.mpf(math.e**mp.mpf(t))  ) *mp.mpf(1-p*mp.mpf(math.e**mp.mpf(-t))  ) )) ** mp.mpf(b)
	frac = (delta/4) * mp.mpf(1/mp.mpf(moment_generating_function))
	a = int(-mp.log(frac)/t)
	gdl1= gdl_pmf(-a,b,p)
	gdl2 = gdl_pmf(-a-d,b,p)#
	ratio2 = mp.mpf(gdl1)/mp.mpf(gdl2)
	if (mp.mpf(ratio2)!= 0 and 1/mp.mpf(ratio2) > mp.mpf(ratio2)):
		ratio2 = 1/mp.mpf(ratio2)
	max_ratio = ratio2
	epsilon = math.log(max_ratio)/d
	return epsilon

