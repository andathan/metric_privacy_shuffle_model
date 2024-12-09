"""
Measures the utility (MAE) of the proposed mechanism
in a synthetic dataset of uniformly random integers

Usage: exp1.py [CONF FILE] [epsilon] [delta]

e.g. for Figure 4:
python3 exp1.py confg 0.1 0.001

"""

import math
import numpy as np
import random
from matplotlib import pyplot as plt
from datetime import *
import sys
import mpmath as mp
from matplotlib import rc
from matplotlib.ticker import MaxNLocator
import sys

sum_elements = 2000 #the number of elements to sum from the SGDL distribution (2000 should be enough for this, but increase for better accuracy)
mp.dps = 500  #depth of mpmath library since we are handling extremely small numbers
d =1 

if (len(sys.argv)!=4):
	print("Wrong usage...")
	print("exp1.py [conf file] [epsilon] [delta]")
	exit()



matched_geo_shuffle_epsilon = -1
def parse_confg(file_name):
	MAX_RUNS = n =  k =  debias = truncate = query = -1 
	conf_file = open(file_name,"r")
	lines = conf_file.readlines()
	for line in  lines:
		splited = line.split("=")
		command = splited[0].strip()
		value = splited[1].strip()
		if (command == "MAX_RUNS"):
			MAX_RUNS = int(value)
		elif (command == "n"):
			n = int(value)
		elif (command == "k"):
			k = int(value)
		elif (command == "query"):
			query = value
			if (query!="avg" and query!="sum"):
				print("Error, query must be [avg] or [sum]!")
				exit()
		elif (command == "truncate_geo"): #truncate the output of the geometric mechanism or not, default is 0
			truncate = int(value)
		elif (command == "RR-mechanism"): #set to d to run RR in metric privacy (as in the paper)
			type_of_rr = value
			if (type_of_rr != "2d" and type_of_rr!="d"):
				print("Type of rr mechanism must be either d privacy or 2d privacy")
				exit()
		elif (command == "debias"):
			debias = int(value)
		elif (command == "databases"):
			databases_list = value.split(",")
	if (MAX_RUNS == -1 or truncate == -1 or n == -1 or k == -1 or debias == -1 or len(databases_list) == 0):
		print("confg file is insufficient!")
		exit()
	return MAX_RUNS, n, k, debias, databases_list , type_of_rr,truncate, query



def debias_function(r,n,l,sum_all_bits):
	"""
	Debias function of RR-Shuffle (operation of the analyst to improve accuracy)
	"""
	try:
		z = (r/(r - (l))) * (sum_all_bits - ((l)/2))
	except:
		return sum_all_bits
	return z

def rr_shuffle(secret, l, r, n , debias):
	#l bits [ARE EXPECTED TO BE] random, the others will be true
	result = 0
	for i in range(r):
		coin = ber((l)/r)
		if (coin==0):
			result+=secret[i]
		else:
			result+=ber(0.5)	
	if (debias == 1):
		return debias_function(r,n,l,result)
	else:
		return result

def avg(list): #returns avg of list
	return sum(list)/len(list)

def ber(p):
	return np.random.binomial(1,p)


def geometric (epsilon,value,max_value,c):
	"""
	Geometric Mechanism
	"""
	p = 1 - np.exp(-epsilon)
	Z = np.random.geometric(p) - np.random.geometric(p)
	Z_prime = value + Z
	
	return Z_prime





def use_binary_search (target_e,delta,n,e_threshold,MAX_E,MIN_E):
	"""
	do a binary search to find the correct epsilon of Geo-Shuffle
	"""
	middle = (MAX_E+MIN_E)/2
	geo_epsilon = middle

	geo_shuffle_epsilon = find_e_of_sgdl(delta,geo_epsilon,n)	
	if (abs(geo_shuffle_epsilon - target_e ) <= e_threshold):
		return geo_epsilon


	if (target_e>geo_shuffle_epsilon):
		return use_binary_search (target_e,delta,n,e_threshold,MAX_E,middle)
	else:
		return use_binary_search (target_e,delta,n,e_threshold,middle,MIN_E)

	print("Binary search failed...")
	exit()



def find_geo_shuffle_epsilon(target_e,delta,n,e_threshold):
	"""
	finds the resulting epsilon of geo-shuffle
	using epsilon_geo as the starting point (i.e. Geo-Local without shuffler)
	Uses binary search to avoid computing the complex 4.1. 
	"""
	if (n==1):
		return target_e

	MAX_E = 5
	MIN_E = target_e
	return use_binary_search(target_e,delta,n,e_threshold,MAX_E,MIN_E)
	
def find_c(delta,epsilon,n):
	"""
	Calculates c of Algorithm 3
	"""
	c = int(- math.log(	delta*(1+math.exp(-epsilon))/(4*n))/epsilon)
	return c


def calculate_c(n,delta):
	w = calculate_w(n, delta)
	a = 1
	b = w - 2
	c = 1
	discriminant = b**2 - 4*a*c
	if discriminant < 0:
        # No real roots, the inequality is always non-negative
 		return None
	elif discriminant == 0:
        # One real root
		root = -b / (2*a)
		return [root]
	else:
        # Two real roots
		root1 = (-b + math.sqrt(discriminant)) / (2*a)
		root2 = (-b - math.sqrt(discriminant)) / (2*a)
	c = root1 * (1-np.exp(-epsilon))/n
	return math.ceil(c)



def find_c_sgdl(delta,epsilon,n):
	w = 2*n * math.log((1 - math.pow(1 - delta, 1/n)) / 2)
	b = w-2
	discriminant = b**2 - 4
	root1 = (-b + math.sqrt(discriminant)) /2
	#print ("c is ",c)
	c = root1 * (1-np.exp(-epsilon))/n
	return math.ceil(c)


def debias_geo_shuffle (sum_of_all_bits,n,k,c):
	"""
	Debias function of RR-Shuffle (operation of the analyst to improve accuracy)
	"""
	debiased = sum_of_all_bits - (n * c)
	print(debiased)
	if (debiased<0):
		return 0
	elif (debiased>n*k):
		return n*k
	else:
		return debiased


def gdl_pmf(m,b, p): 
	"""
	Geo-Shuffle:
	calculates the PMF of SGDL on some point m
	"""
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
	"""
	SGDL-Shuffle:
	Finds the resulting epsilon of SGDL-Shuffle (Theorem 5.1)
	"""

	t = 2* geo_epsilon/(math.sqrt(b))
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


def find_avg_utility (mechanism_list, target_secret):
	sum_result = 0 
	for item in mechanism_list:
		sum_result += abs(item - target_secret)
	return sum_result/len(mechanism_list)


def init_file(n,k,MAX_RUNS,delta,epsilon,debias,truncate):
	#open file for writing

	file_name =  "geo_vs_rr_results_var_n" +  datetime.now().strftime("%d_%m_%Y_%H_%M")
	f = open(file_name, "w")
	#print basic information to file
	f.write("*****VAR N************************\n")
	f.write("MAX n =  " +  str(n) + "\n")
	f.write("k =  " +  str(k) + "\n")
	f.write("delta =  " +  str(delta) + "\n")
	f.write("epsilon =  " +  str(epsilon) + "\n")
	f.write("MAX_RUNS =  " +  str(MAX_RUNS) + "\ƒa =n")
	f.write("Debias RR:  " +  str(debias) + "\n")
	f.write("Truncate GEO:  " +  str(truncate) + "\n")
	f.write("**********************************\n")
	return f

def write_autosave(m, result_rr, result_geo, target_secret):
	#print to console
	print("Current run is", m )
	print("Avg rr is", avg(result_rr))
	print("Avg geo is " , avg(result_geo))
	#write to file
	f.write("-------------------------------\n")
	f.write("Current run is " +  str(m) + "\n")
	f.write("Avg Utility RR is " +  str((abs(avg(result_rr) - (target_secret)))) + "\n")
	f.write("Avg Utility GEO is  "+ str((abs(avg(result_geo) - (target_secret))))+  "\n") 
	f.write("------------------------------- \n")



def calculate_lambda (n, target_e,k,e_threshold,delta):
	for l in range(1,n*k):
		try:#check if sqrt does not take a negative input
			e0 =  math.sqrt(  (32*math.log(4/delta))  /  ((l)-math.sqrt(2*(l)*math.log(2/delta)))  )
		except:
			continue	
		e = e0 
		if (abs(e-target_e)<= e_threshold): 
			return l
	#nothing found...
	#search again with greater theshold
	return calculate_lambda (n, target_e,k,0.001+e_threshold,delta)



def calculate_utility(mechanism_name, database, MAX_RUNS, n, k, epsilon,delta,debias,truncate,query):
	"""
	Finds the utility of a mechanism
	"""
	global  matched_geo_shuffle_epsilon
	if (query=="avg"):
		divide = n 
	elif (query=="sum"):
		divide =1 
	else: 
		print("Wrong query=",query)
		exit()
	results = []
	target_secret = (sum(database))
	#RR
	if (mechanism_name == "RR"):
		if ( math.sqrt(  (32*math.log(4/delta))  /  ((n*k)-math.sqrt(2*(n*k)*math.log(2/delta)))  )   >= epsilon):
			l = n * k
		else:
			l = int(calculate_lambda (n, epsilon, k, 0.001,delta))
		r = n*k #we will run the binary RR-mechanism with n*k number of bits
		#now make a list with size n*k with {target_secret} ones and the rest zeros
		secret = []
		for s in range(target_secret):
			secret.append(1)
		for s in range(r-target_secret):
			secret.append(0)
		#random.shuffle(secret) #shuffle it to make it random //maybe this is not necessary 
		#verify that the conversion is correct:
		assert (sum(secret)==target_secret)
		#calculate noise:
		for m in range(MAX_RUNS):
			results.append(rr_shuffle(secret,l,r,n,debias)/divide)
	#GEO Local
	elif (mechanism_name == "GEO"):
		for m in range(MAX_RUNS):
			result_geo_value = 0 
			for u in range(n):#run for each user:
				result_geo_value +=geometric(epsilon,database[u],k,-1)
			results.append(result_geo_value/divide)
	######### GEO SHUFFLE ##########
	elif (mechanism_name == "GEO-SHUFFLE"):
		#find the corresponding epsilon
		if (matched_geo_shuffle_epsilon==-1):
			e_threshold =  0.01
			epsilon = find_geo_shuffle_epsilon(epsilon,delta,n,e_threshold)
			matched_geo_shuffle_epsilon= epsilon
		else:
			epsilon = matched_geo_shuffle_epsilon
		c = find_c (delta,epsilon,n) # using geo_epsilon for c 
		for m in range(MAX_RUNS):
			result_geo_shuffle = 0 

			for u in range(n):#run for each user:
				reported_value = geometric(epsilon,database[u],k,c)
				result_geo_shuffle +=reported_value
				#if you add c to the above line, add this command as well
			#debias_result = debias_geo_shuffle(result_geo_shuffle,n,k,c)
			debias_result = result_geo_shuffle
			results.append(debias_result/divide)
	#########  SGDL SHUFFLE  ##########
	elif (mechanism_name ==  "SGDL-SHUFFLE"): 
		if (divide == n):
			epsilon = epsilon*n #avg query
		for m in range(MAX_RUNS):
			result_geo_central = geometric(epsilon,target_secret,k,-1)
			results.append(result_geo_central/divide)
	else:
		print("Wrong mechanism name! Exiting...")
		exit()

#find utility
	utility = find_avg_utility (results, target_secret/divide)
	return utility 




def list_from_file(data):
	data_list = data.split(",")
	#pop \n 
	data_list.pop()
	#convert to int 
	data_list = [int(x) for x in data_list]
	return data_list 

def print_parameters(MAX_RUNS, n, k, delta, epsilon, debias, databases_list,truncate):
	print("###################")
	print("Running with parameters:")
	print("MAX_RUNS =",MAX_RUNS)
	print("MAX N n =",n)
	print("k =",k)
	print("epsilon =",epsilon)
	print("delta =",delta)
	print("debias =",debias)
	print("Truncate geo =",truncate)
	print("Databases = ", databases_list)

	print("###################")


def write_list_to_file (caption, lst,f):
	f.write (caption)
	f.write(str(lst[0]))
	lst.pop(0)
	for item in lst:
		f.write("," + str(item))
	f.write("\n")


def check_if_possible_rr (max_n, k, epsilon , delta):
	if ((max_n*k)<= 14 * math.log(4/delta)):
		return 0 

	try:
		if ( math.sqrt(  (32*math.log(4/delta))  /  ((max_n*k)-math.sqrt(2*(max_n*k)*math.log(2/delta))))   >= epsilon):
			return 0 
	except: 
		print("ERROR")
		return 0
	
	return 1

	
def main():
	global matched_geo_shuffle_epsilon
	conf_filename = sys.argv[1]
	epsilon = float(sys.argv[2])
	delta = float(sys.argv[3])
	#read config
	MAX_RUNS, n, k, debias, databases_list, type_of_rr, truncate, query = parse_confg(conf_filename)
	print_parameters(MAX_RUNS, n, k, delta, epsilon, debias, databases_list,truncate)

	#Initialize File 
	#f = init_file(n,k,MAX_RUNS,delta,epsilon,debias,truncate)
	for db_filename in databases_list:
		n_list = []
		rr_utility_list = []
		geo_utility_list = []
		geo_shuffle_utility_list = []
		#geo_central_utility_list = []
		sgdl_shuffle_utility_list = []
		#run for each databases indicated at the confg file:
		print("Running for ", db_filename)
		#f.write("********************\n")
		#f.write("Running for db " +db_filename +"\n")
		#f.write("********************\n")

		#calculate the steps to run
		j_list = []
		for j in range(2,50,1):
			j_list.append(j)

		for j in range(50,150,15):
			j_list.append(j)

		for j in range(150,n,30):
			j_list.append(j)


		for j in j_list:
			matched_geo_shuffle_epsilon = -1
			database = []
			n_list.append(j)
			print("Running for n = ", j)
			#f.write("--------------------\n")
			#f.write("Running for n = " + str(j)  +"\n")
			#f.write("--------------------\n")
			#Initialize Database
			db_file = open(db_filename, "r")		
			database = list_from_file(db_file.read())
			#keep only the first j values
			database = database[0:j]
			assert(len(database)==j) #check if database is equal to the number of users
			#calculate utility

			
			if (check_if_possible_rr (j, k, epsilon , delta))==0:
				rr_utility = None
			else:
				rr_utility = calculate_utility("RR", database, MAX_RUNS, j, k, epsilon,delta,debias,truncate,query)
			
			sgdl_shuffle_utility = calculate_utility("SGDL-SHUFFLE", database, MAX_RUNS, j, k, epsilon,delta,debias,truncate,query)

			geo_utility =  calculate_utility("GEO", database, MAX_RUNS, j, k, epsilon,delta,debias,truncate,query)
			geo_shuffle_utility = calculate_utility("GEO-SHUFFLE", database, MAX_RUNS, j, k, epsilon,delta,debias,truncate,query)
			#geo_central_utility = calculate_utility("GEO-CENTRAL", database, MAX_RUNS, j, k, epsilon,delta,debias,truncate,query)
		
			
			#f.write("RR Utility =  " +  str(rr_utility) + "\n")
			#f.write("GEO Utility = " +  str(geo_utility) + "\n")
			#f.write("GEO Shuffle = " +  str(geo_shuffle_utility) + "\n")
			#f.write("GEO Central = " +  str(geo_central_utility) + "\n")
			rr_utility_list.append(rr_utility)
			geo_utility_list.append(geo_utility)
			geo_shuffle_utility_list.append(geo_shuffle_utility)
			#geo_central_utility_list.append(geo_central_utility)
			sgdl_shuffle_utility_list.append(sgdl_shuffle_utility)

		#f.write("-------------------\n")

		#write_list_to_file("RR Utility List = ", rr_utility_list,f)
		#write_list_to_file("GEO Utility List = ", geo_utility_list,f)
		#write_list_to_file("GEO Shuffle List = ", geo_shuffle_utility_list,f)
		#write_list_to_file("GEO Central List = ", geo_central_utility_list,f)
		#write_list_to_file("SGDL shuffle List = ", sgdl_shuffle_utility_list,f)

		#write_list_to_file("N List = ", n_list,f)
		#f.write("-------------------\n")
	#f.close()
	
	#make a nice plot

	rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
	rc('text', usetex=True)
	plt.rcParams.update({'font.size': 25})
	plt.rcParams["figure.figsize"] = (20,12)
	plt.plot(n_list,geo_utility_list, color ="orange",label ="Geo local")
	plt.plot(n_list,rr_utility_list, color="red",label ="RR-Shuffle")

	plt.plot(n_list,geo_shuffle_utility_list, color = "cyan",label ="Geo-Shuffle")
	#plt.plot(n_list,geo_central_utility_list, color = "purple",label ="SGDL-Shuffle")
	plt.plot(n_list,sgdl_shuffle_utility_list, color = "m",label ="SGDL-Shuffle")
	

	plt.xlabel("Number of Users")
	plt.ylabel('Utility Loss')
	plt.legend(loc ="upper right" )
	plt.ylim([-0.1, 10])

	name_of_file = "exp1_max_n_"+str(n)+".pdf"
	plt.savefig(name_of_file, format="pdf", bbox_inches="tight")

	plt.show()
main()