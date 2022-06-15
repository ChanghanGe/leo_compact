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

results = []
for i in range(3*60):
    
    DATA = str(ephem.date(ephem.date(OBSERVATION_DATE) + i/24/60/60))
    # print(DATA)
    ori_graph = sat.create_spaceX_graph_with_ground_station(DATA)
    landmark0 = NUMBER_OF_NODES
    landmark1 = NUMBER_OF_NODES + 1
    path = nx.algorithms.shortest_path(ori_graph, landmark0, landmark1, weight='weight')
    path_length = nx.algorithms.shortest_path_length(ori_graph, landmark0, landmark1, weight='weight')
    path_length = path_length / 300 + 0.44
    # print(path_length)
    # print(path)
    tem = []
    for i in range(len(path) - 1):
        x = path[i]
        y = path[i+1]
        tem.append(ori_graph[x][y]['weight'])
    print(path_length, path, tem)
    results.append(path_length)
    
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
