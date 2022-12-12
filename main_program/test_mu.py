import leo_satellites as sat
import network_implementation as net_utils
import networkx as nx
import pandas as pd
import numpy as np
import sys
import copy
import scipy.io as sio
from statsmodels.distributions.empirical_distribution import ECDF
import random
import math
import matplotlib.pyplot as plt
from astropy import constants as const
import os
from IPython.display import clear_output
from datetime import datetime
from multiprocessing import Pool
import scipy
import geometry_functions as geom
from matplotlib import font_manager
import matplotlib
import time
import ephem
import math
import pickle

NUMBER_OF_NODES = sat.SATELLITES_PER_ORBIT * sat.NUMBER_OF_ORBITS
NUMBER_OF_ORBITS = sat.NUMBER_OF_ORBITS
SATELLITES_PER_ORBIT = sat.SATELLITES_PER_ORBIT
NUMBER_OF_LEVELS_list = [1,2,3,4,5]
MAX_HOPS = sys.float_info.max
ALTITUDE = 550
EARTH_RADIUS = const.R_earth.value / 1000
MAX_DISTANCE_BETWEEN_SATS = 2*math.sqrt((ALTITUDE**2+2*EARTH_RADIUS*ALTITUDE))
MAX_NUM_LASER_LINK = 5

#graph = sat.create_full_mesh()
EPOCH = '2022/9/21 00:00:00'
OBSERVATION_DATE = '2022/9/21 00:00:00'

SIMULATION_RANGE = 1200

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import pandas as pd
import seaborn as sns
sns.set_style('darkgrid')

import matplotlib.pyplot as plt

SIMULATION_RANGE = 900

# citys = ['London', 'Boston', 'Shanghai', 'Hong Kong', 'Los Angeles',
#              'Paris', 'New York', 'Tokyo', 'Chicago', 'Singapore',
#              'San Francisco', 'Sydney', 'Toronto', 'Mexico City',
#              'Taipei', 'Washington', 'Beijing', 'Rome', 'Berlin',
#              'Dublin', 'Sao Paulo', 'Moscow', 'Seoul', 'Osaka']

citys = ['London', 'Boston', 'Shanghai', 'Hong Kong', 'Los Angeles']

print(10 * math.log10(600*600), 10 * math.log10(800*800))

import numpy as np
noise_reduction_mean = -2
noise_reduction_width = 1
num_user_each_location = 10

# generate users 
# sat.find_valid_ground_station(24, SIMULATION_RANGE, num_gs = num_user_each_location)

# load users 
with open('./ground_station_location_epoch_2022_09_21_duration_24_sr_1200.json', 'rb') as f:
    valid_gs_all = pickle.load(f)


for t in range(24):

    BD = str(ephem.date(ephem.date(OBSERVATION_DATE) + t/24))
    BD = BD.replace(' ', '-')
    BD = BD.replace(':', '-')
    BD = BD.replace('/', '-')
    print(BD)
    over_all_dict = {}
    lon_all_dict = {}
    lat_all_dict = {}

    for city_id, city in enumerate(citys):

        city_gs_coords = valid_gs_all[city]

        for i in range(SIMULATION_RANGE):
            
            DATA = str(ephem.date(ephem.date(OBSERVATION_DATE) + t/24 + i/24/60/60))
            ori_graph, alt_graph, lon_graph = sat.create_spaceX_graph_with_ground_station_distance_gs_loc(DATA, city, city_gs_coords)
            result = []

            for groundstation, first_level in enumerate(ori_graph):
                if groundstation not in over_all_dict:
                    over_all_dict[groundstation] = {}
                    lon_all_dict[groundstation] = {}
                    lat_all_dict[groundstation] = {}
                my_dict = over_all_dict[groundstation]
                lon_dict = lon_all_dict[groundstation]
                lat_dict = lat_all_dict[groundstation]
                for satellite, second_level in enumerate(first_level):    
                    distance = ori_graph[groundstation][satellite]

                    if distance[1] == False:
                        continue

                    d = distance[0] 
                    alt = alt_graph[groundstation][satellite]
                    lon = lon_graph[groundstation][satellite]
                    result.append((groundstation, satellite, distance[0]))
                    if satellite in my_dict:
                        my_dict[satellite][i] = d
                        lat_dict[satellite][i] = alt
                        lon_dict[satellite][i] = lon
                    else:
                        my_dict[satellite] = {i: d}
                        lat_dict[satellite] = {i: alt}
                        lon_dict[satellite] = {i: lon}
                        

        for gs in range(num_user_each_location):
            my_dict = over_all_dict[gs]
            lon_dict = lon_all_dict[gs]
            lat_dict = lat_all_dict[gs]
            results = {}
            results['time'] = list(range(SIMULATION_RANGE))
            vis = {}
            vis['time'] = list(range(SIMULATION_RANGE))
            dis = {}
            dis['time'] = list(range(SIMULATION_RANGE))
            angle = {}
            angle['time'] = list(range(SIMULATION_RANGE))
            az = {}
            az['time'] = list(range(SIMULATION_RANGE))
            for key in my_dict:
                distance = [  my_dict[key][i] if i in my_dict[key] else 0  for i in range(SIMULATION_RANGE)]
                value = [( 100 * 550 / my_dict[key][i] / 8 +  np.random.normal(noise_reduction_mean, noise_reduction_width))*8 if i in my_dict[key] else 0  for i in range(SIMULATION_RANGE)]
                ang = [  lat_dict[key][i] if i in lat_dict[key] else 0  for i in range(SIMULATION_RANGE)]
                lon = [  lon_dict[key][i] if i in lon_dict[key] else 0  for i in range(SIMULATION_RANGE)]
                # if i in my_dict[key]:
                #     print( my_dict[key][i])
                #ang = [  90 - math.asin(550.0 / my_dict[key][i]) /math.pi * 180 if i in my_dict[key] else 0  for i in range(SIMULATION_RANGE)]
                # ang = [  my_dict[key][i] if i in my_dict[key] else 0  for i in range(SIMULATION_RANGE)]
                if sum(value) == 0:
                    continue
                result = []
                for i in range(SIMULATION_RANGE):
                    if i not in my_dict[key]:
                        result.append(0)
                    else:
                        count = 0
                        while  i+count in my_dict[key]:
                            count+=1
                        result.append(float(count))
                results[key] = value
                vis[key] = result
                dis[key] = distance
                angle[key] = ang
                az[key] = lon

            # data = pd.DataFrame.from_dict(results)
            # data.to_csv('output/rss_' + citys[city] + '_' + BD + '.csv')  
            # data = pd.DataFrame.from_dict(vis)
            # data.to_csv('visibility_' + citys[city] + '_' + BD + '.csv')  
            data = pd.DataFrame.from_dict(dis)
            data.to_csv('output/distance_' + city + '_user_' + str(gs) + '_' + BD + '.csv')  
            data = pd.DataFrame.from_dict(angle)
            data.to_csv('output/alt_' + city + '_user_' + str(gs) + '_' + BD + '.csv')  
            data = pd.DataFrame.from_dict(az)
            data.to_csv('output/az_' + city + '_user_' + str(gs) + '_' + BD + '.csv')  