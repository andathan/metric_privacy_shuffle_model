import sys
import mpmath as mp

"""
parametrizes the experiment
"""
############## PARAMETERS ####################
db_name = "exp2/austin.csv" #gowalla.txt or austin.csv
#to_do: For gowalla.txt change at make_grid() gps coordinates of corners

desired_epsilon = float(sys.argv[1]) #0.001
MAX_RUNS = 1000
number_of_rows = int(sys.argv[2]) #1000 # how many rows to read. this will be the number of users n

size_of_grid = int(sys.argv[3]) #10000 # megetheros kathe square: 120 meters # size_of_grid x size_of_grid. this will be k s.t. users have value in [0,k]	
PRIVACY_RADIUS = int(sys.argv[4]) #5 #privacy radius measured in meters 



INCREASE_GRID = 0.00056 # increase the grid size (0.1 = 10% increase)
RANDOM_USERS = 1 #0 = pick the first n of each dataset. 1: pick n random users from each dataset




matched_geo_shuffle_epsilon = -1 #if set, the geo_shuffle will skip trying to find the correct epsilon and use this one


delta = 0.001
quick = 0 #the RR will run quickly but with an estimation



###### INIT ###########


invalid = 0
if number_of_rows =="ALL":
	number_of_rows = 479529
if(db_name == "exp2/austin.csv"):
	SEPARATOR = ','
elif (db_name == "gowalla.txt"):
	SEPARATOR = '\t'

d =1 
sum_elements = 30000
mp.dps = 500

size_of_square = -1 
rr_lambda = -1 
rr_array = []
max_possible_x = -1 
max_possible_y = -1
########################################


