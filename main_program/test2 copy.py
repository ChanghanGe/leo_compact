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
OBSERVATION_DATE = '2022/9/21 00:00:00'

SIMULATION_RANGE = 300

citys = ['London', 'Boston', 'Shanghai', 'Hong Kong', 'Los Angeles']

results = []
my_dict = {}
for i in range(SIMULATION_RANGE):
    
    DATA = str(ephem.date(ephem.date(OBSERVATION_DATE) + i/24/60/60))
    ori_graph = sat.create_spaceX_graph_with_ground_station_distance(DATA)
    
    result = []
    for groundstation, first_level in enumerate(ori_graph):
        max_rss = -100.0
        for satellite, second_level in enumerate(first_level):    
            distance = ori_graph[groundstation][satellite]
            if distance[1] == False:
                continue
            d = -43 - 40 * math.log( distance[0] / 1000, 10)
            # result.append((groundstation, satellite, distance[0]))
            # if satellite in my_dict:
            #     my_dict[satellite][i] = distance[0]
            # else:
            #     my_dict[satellite] = {i: distance[0]}
            max_rss = max(max_rss, d)
        if groundstation not in my_dict:
            my_dict[groundstation] = [max_rss]
            my_dict[groundstation] = [max_rss]
        else:
            my_dict[groundstation].append(max_rss)
            my_dict[groundstation] = [max_rss]

print(len(my_dict))

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import pandas as pd
import seaborn as sns
sns.set_style('darkgrid')

import matplotlib.pyplot as plt

results = {}
# for key in my_dict:
#     result = [ -43 - 40 * math.log( my_dict[key][i] / 1000, 10) if i in my_dict[key] else 0  for i in range(20)]
#     results[key] = result
#     #print(key, result)
results['time'] = list(range(SIMULATION_RANGE))
keys = [key for key in my_dict]
for i, key in enumerate(my_dict):
    results[citys[i]] = my_dict[key]

data = pd.DataFrame.from_dict(results)
fig = px.line(data, x="time", y=citys, title='RSS simulation', line_shape='linear')
fig.write_image("RSS_simulation.png", format='png')

for city in citys:
    fig = px.line(data, x="time", y=city, title='RSS simulation '+ city, line_shape='linear')
    fig.write_image("RSS_simulation_" + city + ".png", format='png')
import math
for city in citys:
    key = city + 'angle'
    data[key] = [ math.acos(550.0 / value1) /math.pi * 360 for value1 in getattr(data, city)]
    fig = px.line(data, x="time", y=key, title='Angle simulation '+ key, line_shape='linear')
    fig.write_image("Angle_simulation_" + city + ".png", format='png')
    
for i, city_start in enumerate(citys):
    for j, city_end in enumerate(citys):
        if j <= i:
            continue
        key = 'path_from_' + city_start + '_to_' + city_end
        data[key] = [min(value1, value2) for value1, value2 in zip(getattr(data, city_start),getattr(data, city_end))]
        fig = px.line(data, x="time", y=key, title='RSS simulation '+ key, line_shape='linear')
        fig.write_image("THOUGHPUT_simulation_" + key + ".png", format='png')

for i, city_start in enumerate(citys):
    for j, city_end in enumerate(citys):
        if j <= i:
            continue
        key = 'path_from_' + city_start + '_to_' + city_end
        file_name = key
        file_name = file_name.replace(' ', '_')
        with open(file_name, 'w') as f:
            for k in range(SIMULATION_RANGE):
                -43 - 40 * math.log( distance[0] / 1000, 10)
                value = 100 / (10 ** ((data[key][k] + 43) / -40) * 1000 / 550)
                f.write(str(float(k+1)) + ' ' + str(int(value)) + '\n')

import json
json.dump(results, open('log.log', 'w'))
# # ori_graph = sat.create_spaceX_graph(OBSERVATION_DATE)
# # spaceX_constellation = sat.constellationFromSaVi()
# # spaceX_positions = sat.positionsAtTime(spaceX_constellation, OBSERVATION_DATE)
# # distances = sat.distances_pathLabel(spaceX_positions)
# # graph = ori_graph

# landmark0 = NUMBER_OF_NODES
# landmark1 = NUMBER_OF_NODES + 1
# path = nx.algorithms.shortest_path_length(ori_graph, landmark0, landmark1, weight='weight')
# print(path)
