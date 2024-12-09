import numpy as np
import math  
from parameters import *



def find_euclenian_distance(target_x, target_y, x, y):
	global size_of_square
	assert(size_of_square>0)
	x_distance = abs (target_x-x) * size_of_square
	y_distance = abs(target_y - y) *size_of_square

	eu_distance = math.sqrt(x_distance**2 + y_distance**2)
	return eu_distance



def find_avg_grid_coordinates(label,center_x, center_y,x_grid,y_grid):
	if (center_x != len(x_grid)):
		grid_centroid_x = (x_grid[center_x] + x_grid[center_x+1])/2
	else:
		grid_centroid_x = (x_grid[center_x] + max_possible_x)/2

	if(center_y!=len(y_grid)):
		grid_centroid_y = (y_grid[center_y] + y_grid[center_y+1])/2
	else:
		grid_centroid_y = (y_grid[center_y] + max_possible_y)/2
	print("Centroid when using ",label, "is at GPS:", grid_centroid_x,",",grid_centroid_y)



def degrees_to_meters(lat1, lon1, lat2, lon2):
    R = 6378.137 # Radius of earth in KM
    dLat = lat2 * math.pi / 180 - lat1 * math.pi / 180;
    dLon = lon2 * math.pi / 180 - lon1 * math.pi / 180;
    a = math.sin(dLat/2) * math.sin(dLat/2) + math.cos(lat1 * math.pi / 180) * math.cos(lat2 * math.pi / 180) * math.sin(dLon/2) * math.sin(dLon/2);
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a));
    d = R * c #
    return d * 1000 # meters


def make_grid(db_lat, db_long):
	global max_possible_x,max_possible_y
	global size_of_square
	global INCREASE_GRID 

	min_x =  30.091  - INCREASE_GRID * 30.091
	max_x =  30.95534 + INCREASE_GRID *   30.95534
	min_y = -98.27048 - INCREASE_GRID *  98.27048
	max_y = -97.38191 + INCREASE_GRID *  97.38191


	max_possible_x = max_x 
	max_possible_y = max_y
	square_size = degrees_to_meters(max_x,min_y,min_x,min_y)/size_of_grid
	print("Each square will have a size of ",square_size, "meters")
	size_of_square = square_size
	x_step = (max_x - min_x)/size_of_grid
	y_step = (max_y - min_y)/size_of_grid
	#make a table with x coordinates
	x_grid = []
	pos = min_x
	for i in range (size_of_grid):
		x_grid.append(pos)
		pos += x_step
	#make a table with y coordinates	
	y_grid = []
	pos = min_y
	for i in range (size_of_grid):
		y_grid.append(pos)
		pos += y_step
	return x_grid, y_grid,square_size


def user_to_grid (user_lat, user_long,x_grid,y_grid):
	user_x_pos = 0 
	found = 0
	for i in range(1,len(x_grid)):
		if (user_lat<x_grid[i]):
			user_x_pos = i -1 
			found = 1
			break
	if (found==0):
		user_x_pos = len(x_grid)
	user_y_pos = 0 
	found = 0
	for i in range(1,len(y_grid)):
		if (user_long< y_grid[i]):
			user_y_pos = i -1  
			found = 1
			break 
	if (found==0):
		user_y_pos = len(y_grid)

	return user_x_pos,user_y_pos

def trim_invalid (location,k):
	global invalid 
	if (location>=k):
		invalid +=1
		return k 
	elif (location<0):
		invalid +=1
		return 0 
	else:
		return location


def avg(alist):
	return sum(alist)/len(alist)
