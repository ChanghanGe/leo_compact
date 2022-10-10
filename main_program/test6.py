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
        max_rss = 10000
        for satellite, second_level in enumerate(first_level):    
            distance = ori_graph[groundstation][satellite]
            if distance[1] == False:
                continue
            d = distance[0]
            # result.append((groundstation, satellite, distance[0]))
            # if satellite in my_dict:
            #     my_dict[satellite][i] = distance[0]
            # else:
            #     my_dict[satellite] = {i: distance[0]}
            max_rss = min(max_rss, d)
        if groundstation not in my_dict:
            my_dict[groundstation] = [max_rss]
        else:
            my_dict[groundstation].append(max_rss)

print(len(my_dict))

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import pandas as pd
import seaborn as sns
sns.set_style('darkgrid')

import matplotlib.pyplot as plt

results = {}
results['time'] = list(range(SIMULATION_RANGE))
keys = [key for key in my_dict]
for i, key in enumerate(my_dict):
    results[citys[i]] = my_dict[key]
    
print(results)

data = pd.DataFrame.from_dict(results)
import math
for city in citys:
    key = city + 'angle'
    print(getattr(data, city))
    data[key] = [ 90 - math.asin(550.0 / value1) /math.pi * 180 for value1 in getattr(data, city)]
    fig = px.line(data, x="time", y=key, title='Angle simulation '+ key, line_shape='linear')
    fig.write_image("Angle_simulation_" + city + ".png", format='png')
