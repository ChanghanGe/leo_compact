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
SIMULATION_RANGE = 100

citys = ['London', 'Boston', 'Shanghai', 'Hong Kong', 'Los Angeles']
results = {}
for i in range(SIMULATION_RANGE):
    DATA = str(ephem.date(ephem.date(OBSERVATION_DATE) + (i+ 450)/24/60/60))
    
    ori_graph = sat.create_spaceX_graph_with_ground_station(DATA)
    
    for j, city_start in enumerate(citys):
        for k, city_end in enumerate(citys):
            if city_start != 'Hong Kong' or city_end != 'Los Angeles':
                continue
            if k <= j:
                continue
            key = 'path_from_' + city_start + '_to_' + city_end
            if key not in results:
                results[key] = []
            landmark0 = NUMBER_OF_NODES + j
            landmark1 = NUMBER_OF_NODES + k
            path_length = nx.algorithms.shortest_path_length(ori_graph, landmark0, landmark1, weight='weight')
            path_length = path_length / 300 + 0.44
            results[key].append(path_length)
            print(i, key, path_length)
            path = nx.algorithms.shortest_path(ori_graph, landmark0, landmark1, weight='weight')
            print(path)
            print(ori_graph[landmark0])

    # landmark0 = NUMBER_OF_NODES
    # landmark1 = NUMBER_OF_NODES + 1
    # landmark2 = NUMBER_OF_NODES + 2
    # landmark3 = NUMBER_OF_NODES + 3
    # path_length = nx.algorithms.shortest_path_length(ori_graph, landmark0, landmark1, weight='weight')
    # path_length2 = nx.algorithms.shortest_path_length(ori_graph, landmark0, landmark2, weight='weight')
    # path_length3 = nx.algorithms.shortest_path_length(ori_graph, landmark2, landmark3, weight='weight')
    # print(i,path_length, path_length2, path_length3)
    # lengths.append(path_length)

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import pandas as pd
import seaborn as sns
sns.set_style('darkgrid')

import matplotlib.pyplot as plt

results['time'] = list(range(SIMULATION_RANGE))

data = pd.DataFrame.from_dict(results)

keys = []
for i, city_start in enumerate(citys):
    for j, city_end in enumerate(citys):
        if j <= i:
            continue
        key = 'path_from_' + city_start + '_to_' + city_end
        keys.append(key)
        fig = px.line(data, x="time", y=key, title='Delay simulation '+ key, line_shape='linear')
        fig.write_image("Delay_" + key + ".png", format='png')

fig = px.line(data, x="time", y=keys, title='Delay simulation', line_shape='linear')
fig.write_image("Delay.png", format='png')

# data = pd.DataFrame.from_dict(results)
# fig = px.line(data, x="time", y='length', title='RSS simulation', line_shape='linear')
# fig.write_image("length_simulation.png", format='png')
