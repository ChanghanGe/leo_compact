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

SIMULATION_RANGE = 1200

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
            d = 100 * 550 / distance[0] 
            result.append((groundstation, satellite, distance[0]))
            if satellite in my_dict:
                my_dict[satellite][i] = d
            else:
                my_dict[satellite] = {i: d}
        break

print(len(my_dict))

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import pandas as pd
import seaborn as sns
sns.set_style('darkgrid')

import matplotlib.pyplot as plt

SIMULATION_RANGE = 600

results = {}
results['time'] = list(range(SIMULATION_RANGE))
for key in my_dict:
    value = [  my_dict[key][i] if i in my_dict[key] else 0  for i in range(SIMULATION_RANGE)]
    if sum(value) == 0:
        continue
    result = []
    for i in range(SIMULATION_RANGE):
        if value[i] == 0:
            result.append(0)
        else:
            count = 0
            while my_dict[key][i+count] != 0:
                count+=1
            result.append(float(count))
    results[key] = result
    print(key, result)
    

data = pd.DataFrame.from_dict(results)
data.to_csv('out2.csv')  