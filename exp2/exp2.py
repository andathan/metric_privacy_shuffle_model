"""
Measures the utility (MAE) of the proposed mechanism
in real world location dataset

Usage: exp2.py [desired_epsilon] [n] [size_of_grid] [privacy radius (in squares)]


e.g. for Figure 5:
python3 exp2.py 0.15 4000 1000 600

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

from rr_functions import *
from geo_functions import *
from grid_functions import *
from parameters import *
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)


'''
#usage:

'''
if (len(sys.argv)!=5):
	print("Wrong usage...")
	print("exp2.py [desired_epsilon] [n] [size_of_grid] [privacy radius (in meters)]")
	exit()


print("-----------------")
print("Database = ",db_name)
print("MAX RUNS = ",MAX_RUNS)
print("n is",number_of_rows)
#print("SIZE OF GRID IS",size_of_grid, "enlarged by ", INCREASE_GRID)
print("SIZE OF GRID IS",size_of_grid)
print("[Privacy RADIUS:]}",PRIVACY_RADIUS,"meters")


def write_list_to_file (caption, lst,f):
	f.write (caption)
	f.write(str(lst[0]))
	lst.pop(0)
	for item in lst:
		f.write("," + str(item))
	f.write("\n")



def apply_mechanism(mechanism_name, database,n,k,epsilon,delta):
	global quick 
	global matched_geo_shuffle_epsilon
	global rr_lambda
	global rr_array
	#print("n is ",n)
	#print("k is ",k)
	target_secret = sum(database) 
	divide = n #avg query
	#print("TARGET SECRET IS ",target_secret)
	#RR
	r = n*k #we will run the binary RR-mechanism with n*k number of bits

	######## RR ########
	if (mechanism_name == "RR"):
		if (rr_lambda==-1):
			l = int(calculate_lambda (n, epsilon, k, e_threshold,delta))
			rr_lambda = l
			secret = rr_array
		else:
			l = rr_lambda
		#now make a list with size n*k with {target_secret} ones and the rest zeros

		secret_ones = [1] * target_secret
		secret_zeros = [0] * (r-target_secret)
		
		secret = secret_ones + secret_zeros

		#random.shuffle(secret) #shuffle it to make it random //maybe this is not necessary 
		#verify that the conversion is correct:
		#assert (sum(secret)==target_secret)
		#calculate noise:
		return (rr_shuffle(secret,l,r,n,1)/divide)

	######### GEO LOCAL ##########
	elif (mechanism_name == "geo"):
		result_geo = 0 
		for u in range(n):#run for each user:
			result_geo +=geometric(epsilon,database[u],k,-1)
		return (result_geo/divide)

	######### GEO SHUFFLE ##########
	elif (mechanism_name == "geo_shuffle"):
		#now change epsilon 
		if (matched_geo_shuffle_epsilon==-1):
			epsilon = find_geo_shuffle_epsilon(epsilon,delta,n,e_threshold)
			matched_geo_shuffle_epsilon= epsilon
		else:
			epsilon = matched_geo_shuffle_epsilon
		c = find_c (delta,epsilon,n) # using geo_epsilon for c 
		result_geo_shuffle = 0 
		for u in range(n):#run for each user:
			result_geo_shuffle +=geometric(epsilon,database[u],k,c)
		debias_result = debias_geo_shuffle(result_geo_shuffle,n,k,c)
		return (debias_result/divide)
		
	######### GEO CENTRAL ##########
	elif (mechanism_name ==  "sgdl_shuffle"): 
		if (divide == n):
			epsilon = epsilon*n #avg query
		result_sgdl_shuffle = geometric(epsilon,target_secret,k,-1)
		return (result_sgdl_shuffle/divide)
	else:
		print("Wrong mechanism name! Exiting...")
		exit()


def find_utility_of_mechanism(mechanism_name, user_x_pos_db, user_y_pos_db,n,k,epsilon,MAX_RUNS,grid_centroid_x_pos,grid_centroid_y_pos):
	#print("Computing ", mechanism_name, "...")
	results = []
	for i in range(MAX_RUNS):
		x_centroid = trim_invalid(apply_mechanism(mechanism_name, user_x_pos_db,n,k,epsilon,delta),k)   
		y_centroid=  trim_invalid(apply_mechanism(mechanism_name, user_y_pos_db,n,k,epsilon,delta),k) 
		distance_to_actual =  find_euclenian_distance(grid_centroid_x_pos, grid_centroid_y_pos, x_centroid, y_centroid)
		results.append(distance_to_actual)


	return (results)



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

if __name__ == '__main__':
	#read dataset
	sparse = 0 #set to 1 to run the experiment with a sparce dataset

	if (sparse==1):
		print("!!!! Sparse db: Removing users in downtown!!!!")
		db = pd.read_csv(db_name,usecols = ["LATITUDE","LONGITUDE"],sep=SEPARATOR)
		db = db.drop(db[(db['LATITUDE'] > 30.16737) & (db['LATITUDE'] < 30.45908) & (db['LONGITUDE'] >  -97.877248) & (db['LONGITUDE'] <  -97.627248)].index)
		db_name = "austin_sparse.csv"
		db.to_csv(db_name, sep=SEPARATOR)

		#show how it looks?
		#fig = px.scatter_mapbox(db, lat="LATITUDE", lon="LONGITUDE", zoom=10, mapbox_style='open-street-map')
		#fig.show()
	
	#Select random rows?
	if (RANDOM_USERS==1): 
		#print("Picking random users...")
		random.seed(20)
		db_read_for_size = pd.read_csv(db_name)
		skip = sorted(random.sample(range(1,len(db_read_for_size)+1),len(db_read_for_size)-number_of_rows)) #the 0-indexed header will not be included in the skip list
	else:
		skip = None

	
	main_db_lat = pd.read_csv(db_name,usecols = ["LATITUDE"],sep=SEPARATOR, nrows = number_of_rows,skiprows =skip  ).values.tolist()
	main_db_lat = list(itertools.chain(*main_db_lat))


	main_db_long = pd.read_csv(db_name,usecols = ["LONGITUDE"],sep=SEPARATOR,nrows = number_of_rows,skiprows =skip ).values.tolist()
	main_db_long = list(itertools.chain(*main_db_long))

	x_grid, y_grid,square_size = make_grid(main_db_lat,main_db_long) 

	#### plot on a map? ####
	#df=pd.read_csv(db_name,usecols = ["LATITUDE","LONGITUDE"],nrows = number_of_rows,skiprows =skip)
	#fig = px.scatter_mapbox(df, lat="LATITUDE", lon="LONGITUDE", zoom=10, mapbox_style='open-street-map')

	#fig.show()

	###################


	NUMBER_OF_SQUARES =  round(PRIVACY_RADIUS/square_size)

	epsilon = desired_epsilon/NUMBER_OF_SQUARES  #epsilon that will be used for our mechanisms
	e_threshold =   epsilon/10   #used in RR and geo shuffle when trying to find epsilon. how close to the target epsilon is close enough?


	total_geo_local = []
	total_geo_shuffle = []
	total_sgdl_shuffle = []
	total_rr = []
	n_list=[]
	df_rr_shuffle = pd.DataFrame()

	rr_sep_produced_results = [275.8349212104782, 164.37161253875394, 126.81369575103822, 108.47266869906707, 96.53622348190946, 90.32776589394855, 79.21551575843169, 79.8366103939013, 91.51487636381461, 68.78805700320788, 69.82572140333295, 78.9392871154992, 88.577842428972, 70.28988457583904, 85.30478630757678, 51.9427122490086, 66.75123998626952, 95.38359615826145, 82.12204217220288, 66.67914569542322]

	epsilon = desired_epsilon/(NUMBER_OF_SQUARES * math.sqrt(2) )#epsilon that will be used for our mechanisms

	e_threshold =  epsilon/10   #used in RR and geo shuffle when trying to find epsilon. how close to the target epsilon is close enough?



	n_iter_list = []

	n_iter_for_rr = []



	for n_iter in range(500,number_of_rows+1,500):
		n_iter_list.append(n_iter)

	geo_local_error = []
	rr_error = []
	sgdl_shuffle_error = []
	geo_shuffle_error = []
	
	for n_iter in n_iter_list:

		print("\n**********************")
		print("        n=",n_iter)
		print("**********************\n")
		n_list.append(n_iter)
		db_lat = main_db_lat[:n_iter]
		db_long = main_db_long[:n_iter] 
		rr_lambda = -1

		#Find actual centroid
		actual_centroid_x = sum(db_lat)/len(db_lat)
		actual_centroid_y = sum(db_long)/len(db_long)
		#print("Actual centroid is in",actual_centroid_x,",",actual_centroid_y)

		#map users to grid
		user_x_pos_db = []
		user_y_pos_db = []
		for i in range(len(db_lat)):
			user_x_pos,user_y_pos = user_to_grid(db_lat[i],db_long[i],x_grid,y_grid)
			user_x_pos_db.append(user_x_pos)
			user_y_pos_db.append(user_y_pos)



		# n, k notation as used in theory
		# n number of users
		# each having a value in [0,k]
		n = len(user_x_pos_db)
		assert(n==n_iter)

		k = size_of_grid

		###############
		#GRID#
		###############
		grid_centroid_x_pos = (sum(user_x_pos_db)/len(user_x_pos_db))
		grid_centroid_y_pos = (sum(user_y_pos_db)/len(user_y_pos_db))

		###############
		#Geo - Local#
		###############
		geo_local_error.append(find_utility_of_mechanism("geo", user_x_pos_db, user_y_pos_db,n,k,epsilon,MAX_RUNS,grid_centroid_x_pos,grid_centroid_y_pos))

		###############
		#RR-shuffle#
		###############
		if (check_if_possible_rr (n_iter, k, epsilon , delta))==0:
				#print ("**** RR cannot achieve epsilon = ", epsilon, "when n = ", n , "*****")
				rr_error.append([-1000] * MAX_RUNS)
		else:
			rr_error.append(find_utility_of_mechanism("RR", user_x_pos_db, user_y_pos_db,n,k,epsilon,MAX_RUNS,grid_centroid_x_pos,grid_centroid_y_pos))

		###############
		#Geo - Shuffle#
		###############
		geo_shuffle_error.append(find_utility_of_mechanism("geo_shuffle", user_x_pos_db, user_y_pos_db,n,k,epsilon,MAX_RUNS,grid_centroid_x_pos,grid_centroid_y_pos))
		###############
		#SGDL - Shuffle#
		###############
		sgdl_shuffle_error.append(find_utility_of_mechanism("sgdl_shuffle", user_x_pos_db, user_y_pos_db,n,k,epsilon,MAX_RUNS,grid_centroid_x_pos,grid_centroid_y_pos))

		#print("-------- Error (meters) ----------")
		#print("GEO Local:", geo_local_error)
		#print("RR:", rr_error)
		#print("GEO Shuffle:", geo_shuffle_error)
		#print("GEO Central:", sgdl_shuffle_error)

		total_geo_local.append(geo_local_error)
		total_rr.append(rr_error)
		total_geo_shuffle.append(geo_shuffle_error)
		total_sgdl_shuffle.append(sgdl_shuffle_error)
		matched_geo_shuffle_epsilon = -1

		column_name = str(n_iter)
		df_rr_shuffle[column_name] = rr_error[0]



	#box plot
	rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
	rc('text', usetex=True)
	plt.rcParams.update({'font.size': 30})
	plt.rcParams["figure.figsize"] = (20,12)



	plt.figtext(0.69, 0.80,  "RR-Shuffle",
	            backgroundcolor="red", color='black', weight='roman',
	            size='medium')


	plt.figtext(0.69, 0.84,  "Geo Local",
	            backgroundcolor="orange", color='black', weight='roman',
	            size='medium')
	
	plt.figtext(0.79, 0.84,  "SGDL-Shuffle",
	            backgroundcolor="purple", color='white', weight='roman',
	            size='medium')

	

	plt.figtext(0.79, 0.80,  "Geo-Shuffle",
	            backgroundcolor="cyan", color='black', weight='roman',
	            size='medium')
	
	c="orange"
	plt.boxplot(geo_local_error,positions = n_iter_list,widths = 100, patch_artist=True,
	            boxprops=dict(facecolor=c, color=c),
	            capprops=dict(color=c),
	            whiskerprops=dict(color=c),
	            flierprops=dict(color=c, markeredgecolor=c),
	            medianprops=dict(color=c))

	c="red"

	plt.boxplot(rr_error,positions = n_iter_list,widths = 100, patch_artist=True,
	            boxprops=dict(facecolor=c, color=c),
	            capprops=dict(color=c),
	            whiskerprops=dict(color=c),
	            flierprops=dict(color=c, markeredgecolor=c),
	            medianprops=dict(color=c))


	c="cyan"

	plt.boxplot(geo_shuffle_error,positions = n_iter_list,widths = 100, patch_artist=True,
	            boxprops=dict(facecolor=c, color=c),
	            capprops=dict(color=c),
	            whiskerprops=dict(color=c),
	            flierprops=dict(color=c, markeredgecolor=c),
	            medianprops=dict(color=c))



	c="purple"

	plt.boxplot(sgdl_shuffle_error,positions = n_iter_list,widths = 100, patch_artist=True,
	            boxprops=dict(facecolor=c, color=c),
	            capprops=dict(color=c),
	            whiskerprops=dict(color=c),
	            flierprops=dict(color=c, markeredgecolor=c),
	            medianprops=dict(color=c))




	#plt.yscale("log")
	#plt.xscale("log")
	#plt.plot(n_list,total_geo_shuffle , color="cyan")
	#plt.plot(n_list,total_sgdl_shuffle, color ="blue")
	#plt.plot(n_list,total_rr, color ="red")
	plt.ylim(-50, 1500)


	plt.xlabel("Number of Users ")
	plt.ylabel('Utility Loss (meters)' )
	file_name =  "exp2_box_plot_e="+ str(epsilon)+"_delta="+str(delta)+"_" +  datetime.now().strftime("%d_%m_%Y_%H_%M") + ".txt"
	#f = open(file_name, "w")
	#f.write("-------------------\n")
	#write_list_to_file("N List = ", n_iter_list,f)
	#write_list_to_file("RR Utility List = ", rr_error,f)
	#write_list_to_file("GEO Local Utility List = ", geo_local_error,f)
	#write_list_to_file("GEO Shuffle List = ", geo_shuffle_error,f)
	#write_list_to_file("GEO Shuffle List = ", sgdl_shuffle_error,f)
	#f.write("-------------------\n")
	#f.close()
	plt.savefig("exp2_box_plot.pdf", format="pdf", bbox_inches="tight")
	plt.show()
