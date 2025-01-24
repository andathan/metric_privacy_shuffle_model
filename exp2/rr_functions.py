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


def calculate_lambda (n, target_e,k,e_threshold,delta):
	global rr_lambda
	#print(e_threshold)
	#print("CALLED WITH e threshold, ", e_threshold)
	for l in range(1,n*k,1):
		try:#check if sqrt does not take a negative input
			e0 =  math.sqrt(  (32*math.log(4/delta))  /  (l-math.sqrt(2*(l)*math.log(2/delta)))  )
		except:
			continue
		e = e0 
		if (abs(e-target_e)<= e_threshold): 
			#print("Will get",e,"with this l=",l, "(target e was",target_e,")")
			rr_lambda = l 
			return l
	#nothing found...
	#search again with greater theshold
	return calculate_lambda (n, target_e,k,0.0001+ e_threshold ,delta)


def rr_shuffle_quick(number_of_ones, number_of_zeros, l , debias,r ):
	#pick r-l  bits that will be truthfully 
	#the_list = number_of_ones * [1]
	#the_list_2 = number_of_zeros * [0]
	#the_list.append(the_list_2)
	#random.shuffle(the_list)
	#the_list = the_list[:r-l]
	#reuslt = sum(the_list)
	result  = np.random.binomial(n=r-l, p=number_of_ones/r, size=None)
	#and l
	result += l/2
	if (debias == 1):
		return debias_function(r,n,l,result)
	else:
		print("Carefull, not using debias!")
		return result

def rr_shuffle(secret, l, r, n , debias):
	#l bits [ARE EXPECTED TO BE] random, the others will be true
	result = 0
	for i in range(r):
		coin = ber(l/r)
		if (coin==0):
			result+=secret[i]
		else:
			result+=ber(0.5)	
	if (debias == 1):
		return debias_function(r,n,l,result)
	else:
		print("Carefull, not using debias!")
	return result

def ber(p):
	return np.random.binomial(1,p)


def debias_function(r,n,l,sum_all_bits):
	#real debias
	#z = (1/r)*((n)/(n-l))*sum_all_bits - ((l*r)/2)
	#unary debias
	try:
		z = (r/(r - (l))) * (sum_all_bits - ((l)/2))
	except:
		return sum_all_bits
	return z
