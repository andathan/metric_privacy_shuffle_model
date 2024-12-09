import sys
import numpy as np
import random
####
#Usage: create_database.py [n] [k]
###
assert(len(sys.argv)==3)
n = int(sys.argv[1]) 
k = int(sys.argv[2])

#0 default: 50 | 1 random 0,2k3, | 2: normal scale 10 | 3: random 0,k | 4: normal scale 5 |  5: default: 25
def initialiaze_database (random_values,k,n):
	input_array = []
	if (random_values == 0 ): 
		default_secret = 50
		print("All users have the same value of  = ",default_secret)
		input_array = [default_secret] * (n+1)
	elif (random_values == 1):
		print("Running with random values from 0 to 2(k/3)")
		for i in range(n+1):
			random_value = random.randint(0,int(2*k/3)) #[0 to 2/3 k]
			input_array.append(random_value)
	elif (random_values == 2):
		print("Running with normal distribution as input, with loc = k/2 and scale = 10 ")
		normal_distr = np.random.normal(loc = k/2, scale = 10.0, size = n+1)
		input_array = normal_distr.tolist()
		#convert to int
		input_array = [int(x) for x in input_array]
	elif (random_values == 3):
		print("Running with random values from 0 to (k)")
		for i in range(n+1):
			random_value = random.randint(0,int(k)) #[0 to k]
			input_array.append(random_value)
	elif (random_values == 4):
		print("Running with normal distribution as input, with loc = k/2 and scale = 5 ")
		normal_distr = np.random.normal(loc = k/2, scale = 5.0, size = n+1)
		input_array = normal_distr.tolist()
		#convert to int
		input_array = [int(x) for x in input_array]
	elif (random_values == 5 ): 
		default_secret = 25
		print("All users have the same value of  = ",default_secret)
		input_array = [default_secret] * (n+1)
	elif (random_values == 6 ): 
		default_secret = 33
		print("All users have the same value of  = ",default_secret)
		input_array = [default_secret] * (n+1)
	elif (random_values == 7 ): 
		default_secret = 75
		print("All users have the same value of  = ",default_secret)
		input_array = [default_secret] * (n+1)
	elif (random_values == 8 ): 
		default_secret = 0
		print("All users have the same value of  = ",default_secret)
		input_array = [default_secret] * (n+1)
	elif (random_values == 9 ): 
		default_secret = 100
		print("All users have the same value of  = ",default_secret)
		input_array = [default_secret] * (n+1)

	return input_array

def init_file(db_type):
	#open file for writing
	if (db_type == 0):
		database_name = "ALL_50"
	elif (db_type == 1):
		database_name = "RANDOM_0_2K3"
	elif (db_type == 2):
		database_name = "NORMAL_SCALE_10"
	elif (db_type == 3):
		database_name = "RANDOM_0_K"
	elif (db_type == 4):
		database_name = "NORMAL_SCALE_5"
	elif (db_type == 5):
		database_name = "ALL_25"
	elif (db_type == 6):
		database_name = "ALL_33"
	elif (db_type == 7):
		database_name = "ALL_75"
	elif (db_type == 8):
		database_name = "ALL_0"
	elif (db_type == 9):
		database_name = "ALL_100"
	file_name = database_name + "_N_" + str(n) +"_K_" + str(k)
	f = open(file_name, "w")
	return f

def write_list_to_file(write_list, file):
	file.write(str(round(write_list[0],2)))
	skip_first = 0
	for item in write_list:
		#just to omit the first one to avoid having double commas
		if(skip_first == 0):
			skip_first = 1
			continue
		file.write(",")
		file.write(str(round(item,2)))
	file.write("\n")


db_type = int(input("Select type of initial database: \n 0: All users have 50s \n 1: Random values (0, 2k3)\n 2: Normal Distribution with parameters loc = k/2 scale = 10.0 \n 3: Random Values (0,k) \n 4: Normal Distribution with parameters loc = k/2 scale = 5.0 \n 5: All users have a default value of 25\n 6: All users have default value of 33\n 7: All users have default value of 75 \n 8: all users have default value of 0 \n9: All users have default value of 100 \n" )) # 0 default: 50 , 1 random 0,2k3, 2: normal 1 , 3: random 0,k 4: normal 2 5: default: 25 6:default 33
f = init_file (db_type)
#Create Database
input_array = initialiaze_database(db_type,k,n)
write_list_to_file(input_array,f)
f.close()
