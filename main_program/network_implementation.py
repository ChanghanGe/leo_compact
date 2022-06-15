import networkx as nx
import matplotlib.pyplot as plt
import math
import random
import sys
import numpy as np
import copy
import leo_satellites as sat
import pickle
from bokeh.plotting import figure, output_file, show
from heapq import heappush, heappop
from operator import itemgetter
import cProfile
import time
import scipy.io as sio


# This is the way to put labels in them
NUMBER_OF_NODES = sat.SATELLITES_PER_ORBIT * sat.NUMBER_OF_ORBITS
NUMBER_OF_LEVELS = 3
MAX_DISTANCE_BETWEEN_SATELLITES = sys.float_info.max
MAX_HOPS = sys.float_info.max

def create_full_mesh(number_of_nodes):
    created_graph = nx.Graph()
    created_graph.add_nodes_from([i for i in range(0, number_of_nodes)])
    for i in range (0, number_of_nodes):
        for j in range (0, number_of_nodes):
            if (i != j):
                #random_int = np.random.randint(low = 1, high = 2)
                created_graph.add_edge(i, j, weight=1)
    return created_graph



def create_max_links(number_of_nodes, number_of_links):
    created_graph = nx.Graph()
    created_graph.add_nodes_from([i for i in range(0, number_of_nodes)])
    number_of_neighbours = [0 for _ in range (0, number_of_nodes)]

    for node1 in range (0, number_of_nodes):
        possible_neighbours = list(created_graph.nodes())
        for link in range (0, number_of_links):
            if number_of_neighbours[node1] < number_of_links:
                node2 = random.choice(possible_neighbours)
                possible_neighbours.remove(node2)
                random_int = np.random.randint(low=1, high=5)
                created_graph.add_edge(node1, node2, weight=random_int)
                number_of_neighbours[node1] += 1
                number_of_neighbours[node2] += 1

    return created_graph


def create_a_levels(graph, k):
    nodes = graph.nodes
    a_list = [list(nodes)]
    n = graph.number_of_nodes()

    # Creates the A set of different k levels
    prob = 1.0 / math.pow(float(n), 1 / float(k))
    for i in range(1, k):
        a_i = list()
        for node in a_list[i - 1]:
            if random.uniform(0, 1) < prob:
                a_i.append(node)
        a_list.append(a_i)

    a_list.append([])
    print(len(a_list[0]))
    print(len(a_list[1]))
    print(len(a_list[2]))

    return a_list

# USED FOR SATELLITE GRAPHS ------------------------------------------------------------------------------------------------------


def add_level_to_nodes(graph, a_list):
    for level, set in enumerate(a_list):
        for node in set:
            graph.nodes[node]['level'] = level


#Creates for each node a bunch on the following form:
#dictionnary{node_in_bunch: distance_to_node_in_bunch, previous_node_before_reaching_this_node}
#Keeps track of the path to the next level following the same format
def bunches_and_clusters(graph, a_list):
    level = len(a_list)-1
    nodes = graph.nodes

    #Calculate the level where a_list[level] is empty. One should never reach this level
    max_level = 0
    while len(a_list[max_level]) > 0:
        max_level+=1

    LAST_LEVEL_NODES = set(a_list[max_level-1])
    for node in nodes:
        nodes[node]['bunch'] = {}
        nodes[node]['cluster'] = {}

    for node in nodes:
        #print(node)
        last_level_nodes = copy.deepcopy(LAST_LEVEL_NODES)

        #used to extract the closest node in heap
        heap_distances = []

        #keeps track of already explored nodes
        explored_nodes = set()

        # this will contain : the distance to the other elements, the routing node.
        # the last node before arrival and the length of the path to arrive there
        distances_to_others = {node: (MAX_DISTANCE_BETWEEN_SATELLITES, None, None, MAX_HOPS) for node in nodes}
        distances_to_others[node] = 0, node, node, 0

        for neighbor in graph.neighbors(node):
            if neighbor != node:
                weight = graph[node][neighbor]['weight']
                distances_to_others[neighbor] = weight, neighbor, node, 1
                heappush(heap_distances, (weight, neighbor))
        explored_nodes.add(node)

        #this list constains : the distance to other levels and the last node to arrive there
        distance_to_levels = {lvl: (MAX_DISTANCE_BETWEEN_SATELLITES, None) for lvl in range(0, level+1)}

        next_level = 0
        while next_level < max_level and graph.nodes[node]['level'] >= next_level:
            distance_to_levels[next_level] = 0, node
            next_level += 1

        bunch = set()
        bunch.add(node)
        last_level_nodes.discard(node)
        while (len(heap_distances) > 0 and next_level < max_level) or len(last_level_nodes) > 0:
            dist_u, u = heappop(heap_distances)
            explored_nodes.add(u)

            while graph.nodes[u]['level'] >= next_level:
                distance_to_levels[next_level] = dist_u, u
                next_level += 1

            if graph.nodes[u]['level'] == next_level - 1:
                assert(distance_to_levels[next_level][0]>=MAX_DISTANCE_BETWEEN_SATELLITES)

            if graph.nodes[u]['level'] >= next_level - 1:
                bunch.add(u)

            if u in last_level_nodes:
                bunch.add(u)
                last_level_nodes.discard(u)

            neighbors = list(graph.neighbors(u))
            u_neighbors = [node for node in neighbors if node not in explored_nodes]

            for v in u_neighbors:
                alt = dist_u + graph[u][v]['weight']
                current_distance = distances_to_others[v][0]
                if alt < current_distance:
                    _, routing_node, _, hops_to_u, = distances_to_others[u]
                    distances_to_others[v] = alt, routing_node, u, hops_to_u+1
                    heappush(heap_distances, (alt, v))

        #keeps track of all the useful distances
        for node2 in bunch:
            dist, next_node, last_node, hops = distances_to_others[node2]
            graph.nodes[node]['bunch'][node2] = dist, next_node, hops
            graph.nodes[node2]['cluster'][node] = dist, last_node, hops

        graph.nodes[node]['next_levels'] = distance_to_levels
        #print('for node', node, 'the bunch is ', bunch, '\ndistance_to_levels is ', distance_to_levels)

    return


# Gives the distance between two nodes using compact routing
def dist_u_v(graph, u, v):
    w = u
    next_level = 0
    while w not in set(graph.nodes[v]['bunch'].keys()):
        next_level += 1
        u, v = v, u
        w = graph.nodes[u]['next_levels'][next_level][1]

    return graph.nodes[u]['bunch'][w][0] + graph.nodes[v]['bunch'][w][0]


# Gives the distance between two nodes using compact routing and the number of hops
def dist_with_hops(graph, u, v):

    w = u
    next_level = 0
    while w not in set(graph.nodes[v]['bunch'].keys()):
        next_level += 1
        u, v = v, u
        w = graph.nodes[u]['next_levels'][next_level][1]

    dist_u_w, _, hops_u_w = graph.nodes[u]['bunch'][w]
    dist_v_w, _, hops_v_w = graph.nodes[v]['bunch'][w]

    return dist_u_w + dist_v_w, _, hops_u_w + hops_v_w


# NB : Constructs all the required parameters for the graph EXCEPT CLUSTERS AND ROUTING TABLES
def initialize(graph):
    a_list = create_a_levels(graph, NUMBER_OF_LEVELS)
    add_level_to_nodes(graph, a_list)
    bunches_and_clusters(graph, a_list)

#-----------------------------------------------------------------------------------------------------------------------

# HERE ARE DIFFERENT BASIC GRAPH FUNCTIONS FOR TESTS
# Node_0---(1)---Node_1----(3)---Node_2-----(6)----Node_3
def test_path_graph():
    #level == 4
    links_graph = nx.Graph()
    links_graph.add_nodes_from([0, 1, 2 ,3])
    links_graph.add_edge(0, 1, weight=1)
    links_graph.add_edge(1, 2, weight=3)
    links_graph.add_edge(2, 3, weight=6)
    print(links_graph.edges())

    # levels
    a_list = [[0, 1, 2, 3], [1, 2, 3], [2, 3], [3], []]
    add_level_to_nodes(links_graph, a_list)

    for node in links_graph.nodes:
        print(node, 'has level ', links_graph.nodes[node]['level'])

    bunches_and_clusters(links_graph, a_list)
    for node in links_graph.nodes:
        print(node, ' has bunch ', links_graph.nodes[node]['bunch'])
    for node in links_graph.nodes:
        print(node, ' has next levels ', links_graph.nodes[node]['next_levels'])
    for node in links_graph.nodes:
        print(node, ' has cluster', links_graph.nodes[node]['cluster'])

# Creates a triangle graph
def test_triangle_graph():
    #level == 3
    triangle_graph = nx.Graph()
    triangle_graph.add_nodes_from([0, 1, 2])
    triangle_graph.add_edge(0, 1, weight = 1)
    triangle_graph.add_edge(1, 2, weight = 6)
    triangle_graph.add_edge(2, 0, weight=3)

    a_list = [[0, 1, 2], [1, 2], [2], []]
    add_level_to_nodes(triangle_graph, a_list)
    bunches_and_clusters(triangle_graph, a_list)
    nodes = triangle_graph.nodes
    for node in triangle_graph:
        print(node, ' has bunch ', nodes[node]['bunch'])
    for node in triangle_graph:
        print(node, 'has next level', nodes[node]['next_level'])
    for node1 in triangle_graph.nodes:
        for node2 in triangle_graph.nodes:
            print('node ', node1, ' to node ', node2, ' : ', dist_u_v(triangle_graph, node1, node2))

# This was used to compare our bunches with the ones already implemented in Go by Cristina
def test_rnd_graph():
    triangle_graph = nx.Graph()
    triangle_graph.add_nodes_from([0, 1, 2, 3, 4, 5])
    triangle_graph.add_edge(0, 1, weight=146.86047800548656)
    triangle_graph.add_edge(0, 2, weight=122.18837915284743)
    triangle_graph.add_edge(0, 3, weight=17)
    triangle_graph.add_edge(0, 4, weight=10.63014581273465)
    triangle_graph.add_edge(0, 5, weight=124.90796611905904)
    triangle_graph.add_edge(1, 2, weight=77.3692445355388)
    triangle_graph.add_edge(1, 3, weight=161.8919392681427)
    triangle_graph.add_edge(1, 4, weight=136.4734406395618)
    triangle_graph.add_edge(1, 5, weight=41.593268686170845)
    triangle_graph.add_edge(2, 3, weight=131.24404748406687)
    triangle_graph.add_edge(2, 4, weight=112.25417586887359)
    triangle_graph.add_edge(2, 5, weight=98.08159868191383)
    triangle_graph.add_edge(3, 4, weight=25.96150997149434)
    triangle_graph.add_edge(3, 5, weight=141.4390328021229)
    triangle_graph.add_edge(4, 5, weight=115.52056094046635)
    a_list = [[0, 3, 5], [4], [1, 2], []]
    add_level_to_nodes(triangle_graph, a_list)
    bunches_and_clusters(triangle_graph, a_list)
    nodes = triangle_graph.nodes
    for node in triangle_graph:
        print(node, ' has bunch ', nodes[node]['bunch'])
    for node in triangle_graph:
        print(node, ' has next levels ', nodes[node]['next_levels'])

# This is a simulation of a small constellation. It has been constructed in leo_satellites.py
def test_small_sat_graph():
    # level == 3
    small_graph, _ = sat.create_small_sat_graph()
    print(small_graph.edges(data=True))
    a_list = create_a_levels(small_graph, 3)
    print(a_list)
    add_level_to_nodes(small_graph, a_list)
    bunches_and_clusters(small_graph, a_list)
    nodes = small_graph.nodes
    for node in small_graph:
        print(node, ' has bunch ', nodes[node]['bunch'])
    for node in small_graph:
        print(node, ' has cluster', nodes[node]['cluster'])
    print(dist_u_v(small_graph, 0, 2))
    for node in small_graph:
        print(small_graph.nodes[node]['routing_table'])

# This is a test function with the real graph
def test_big_graph():
    graph, _ = sat.create_spaceX_graph()
    print(graph.number_of_nodes())
    #a_list = [list(range(0, 1600)), list(range(0, 800)), list(range(0, 400)), []]
    a_list = create_a_levels(graph, NUMBER_OF_LEVELS)
    print(a_list)
    add_level_to_nodes(graph, a_list)
    bunches_and_clusters(graph, a_list)
    for node in graph.nodes:
        for _, dist_and_node in graph.nodes[node]['bunch'].items():
            dist1, _ = dist_and_node

# -----------------------------------------------------------------------------------------------------------------------
# True computations happens here

# Saves the graph into a pickle
def save_graph():
    graph, distances_matrix = sat.create_spaceX_graph()
    pickle_out1 = open('big_graph.pickle', 'wb')
    pickle.dump(graph, pickle_out1)
    pickle_out2 = open('distance_big_graph.pickle', 'wb')
    pickle.dump(distances_matrix, pickle_out2)


# Loads the graph from a pickle
def load_graph():
    pickle_in1 = open('big_graph.pickle', "rb")
    graph = pickle.load(pickle_in1)
    pickle_in2 = open('distance_big_graph.pickle', "rb")
    distances_matrix = pickle.load(pickle_in2)
    return graph, distances_matrix

#
def compact_rd(graph):
    a_list = create_a_levels(graph, NUMBER_OF_LEVELS)
    add_level_to_nodes(graph, a_list)

    bunches_and_clusters(graph, a_list)
    nodes = graph.nodes
    distances = []
    for node1 in nodes:
        for node2 in nodes:
            distances.append(dist_u_v(graph, node1, node2))

    return distances


def test_load():
    graph, distance_matrix = load_graph()
    print(distance_matrix)


#-----------------------------------------------------------------------------------------------------------------------
# Computes the difference between the size of routing table of compact routing and the one of dijkstra

def compact_table_size():
    graph, distances_matrix = sat.create_spaceX_graph()
    a_list = create_a_levels(graph, NUMBER_OF_LEVELS)
    add_level_to_nodes(graph, a_list)
    bunches_and_clusters(graph, a_list)
    routing_table_sizes = []

    nodes = graph.nodes

    for node in nodes:
        bunch = set(nodes[node]['bunch'].keys())
        cluster = set(nodes[node]['cluster'].keys())
        routing_table_sizes.append(len(bunch.union(cluster)))


    return list(nodes), routing_table_sizes

# Plots compact_table_size
def plot_routing_table_sizes(number_of_nodes):
    nodes_list, compact_routing_list = compact_table_size()
    dij_size = number_of_nodes
    dij_size_list = [dij_size for _ in nodes_list]

    # output to static HTML file
    output_file("./figures/routing_tables_size.html")

    p = figure(plot_width=700, plot_height=700, title="Routing table size (dijkstra vs compact)")

    # add a circle renderer with a size, color, and alpha
    p.circle(nodes_list, compact_routing_list, size=3, color="blue", alpha=0.5)
    p.circle(nodes_list, dij_size_list, size=1, color="orange", alpha=0.5)

    # labels
    p.yaxis.axis_label = "Size of routing tables"
    p.xaxis.axis_label = "Node number"

    # show plot
    show(p)

#-----------------------------------------------------------------------------------------------------------------------
# Computes the difference between the AVERAGE size of routing table of compact routing and the one of Dijkstra

# Computes the sum of bunches sizes in a graph
def total_compact_table_size():
    graph, distances_matrix = sat.create_spaceX_graph()
    a_list = create_a_levels(graph, NUMBER_OF_LEVELS)
    add_level_to_nodes(graph, a_list)
    bunches_and_clusters(graph, a_list)
    routing_table_size = 0

    nodes = graph.nodes
    for node in nodes:
        bunch = set(nodes[node]['bunch'].keys())
        cluster = set(nodes[node]['cluster'].keys())
        routing_table_size += len(bunch.union(cluster))

    return routing_table_size

# Loops on total_compact_table_size() and stores the result in a pickle file
def multiple_routing_table_comparison():
    number_of_comparisons = 20
    compact_routing_table = []
    for i in range(0, number_of_comparisons):
        compact_routing_table.append(total_compact_table_size())

    pickle_out = open('./spacex_graph_data/routing_table_size.pickle', 'wb')
    pickle.dump(compact_routing_table, pickle_out)

# Plots a comparison between average size of compact routing tables and those of dijkstra
def plot_average_routing_tables(number_of_nodes):
    pickle_in = open('./spacex_graph_data/routing_table_size.pickle', "rb")
    compact_routing_table = pickle.load(pickle_in)
    compact_routing_table = [total/number_of_nodes for total in compact_routing_table]
    indexes_list = [i+1 for i, _ in enumerate(compact_routing_table)]
    dij_size_list = [number_of_nodes for _ in indexes_list]

    print(compact_routing_table)

    # output to static HTML file
    output_file("./figures/routing_tables_average.html")

    p = figure(plot_width=700, plot_height=700,  y_axis_type="log", title="Routing tables average size (mult. runs)")

    # add a circle renderer with a size, color, and alpha
    p.circle(indexes_list, compact_routing_table, size=10, color="blue", alpha=0.5)
    p.line(indexes_list, dij_size_list, color="orange", line_width=3)

    # labels
    p.yaxis.axis_label = "Average size of routing tables (log)"
    p.xaxis.axis_label = "Run number"
    
    #show plot
    show(p)


#-----------------------------------------------------------------------------------------------------------------------

# Computes all the distances using dijkstra, compact_routing, the upper stretch as well as the line of sight and
# stores them into a pickle file
def dij_vs_compact():
    sat_graph, distances_matrix = sat.create_spaceX_graph()
    number_of_nodes = sat_graph.number_of_nodes()

    #Puts all the direct sight values into a single dimension list
    direct_list = []
    for dim1 in distances_matrix:
        for dim2 in dim1:
            for dim3 in dim2:
                for dim4 in dim3:
                    direct_list.append(dim4)


    dij_list = []

    t = time.time()
    dij_distance = dict(nx.all_pairs_dijkstra_path_length(sat_graph))
    print('Dijkstra Running Time: ' + str(time.time()-t))

    t = time.time()
    compact_list = compact_rd(sat_graph)
    print('Compact Routing Running Time: ' + str(time.time()-t))

    bound_list = []

    for node1 in range(0, number_of_nodes):
        dist_dict = dij_distance[node1]
        for node2 in range(0, number_of_nodes):
            dij_list.append(dist_dict[node2])
            bound_list.append((2*NUMBER_OF_LEVELS-1)*dist_dict[node2])

    pickle_compact = open('./spacex_graph_data/compact_list2.pickle', 'wb')
    pickle.dump((direct_list, dij_list, compact_list, bound_list), pickle_compact)

    r_list = []
    for i in range(len(compact_list)):
        if dij_list[i] != 0 and compact_list[i] != 0:
            r = compact_list[i]/dij_list[i]
        elif dij_list[i] == 0 and compact_list[i] == 0:
            r = 1
        elif dij_list[i] == 0 and compact_list[i] != 0:
            r = np.inf
        elif dij_list[i] != 0 and compact_list[i] == 0:
            r = 0
        r_list.append(r)
    r_list = np.array(r_list)

    print('ratio_average_over_elemet = ' + str(np.mean(r_list[r_list != np.inf])))
    print('ratio_total = ' + str(np.sum(compact_list)/np.sum(dij_list)))
    sio.savemat('distance_ratio.mat',{'ratio':np.array(r_list)})

# Plots dij_vs_compact() for a small number of values
def plot_lines():
    pickle_in = open('./spacex_graph_data/compact_list2.pickle', "rb")
    direct_list, dij_list, compact_list, bound_list = pickle.load(pickle_in)
    # Draw point based on above x, y axis values.
    plt.plot(direct_list, dij_list, 'bs', direct_list, compact_list, 'g^', direct_list, bound_list, 'r--')

    # Set chart title.
    plt.title("Distances for different routing algorithms")

    # Set x, y label text.
    plt.xlabel('Line of sight')
    plt.ylabel('Bound - red \n Dijkstra dist - blue \n Compact dist - green')
    plt.show()

# Plots dij_vs_compact() for a big number of values
def plot_dij_vs_compact_large():
    pickle_in = open('./spacex_graph_data/compact_list2.pickle', "rb")
    direct_list, dij_list, compact_list, bound_list = pickle.load(pickle_in)

    #selects the same 1000 distances sample from each list, at random
    list_rnd = random.sample(range(0, NUMBER_OF_NODES*NUMBER_OF_NODES-1), 1000)
    direct_list = itemgetter(*list_rnd)(direct_list)
    dij_list = itemgetter(*list_rnd)(dij_list)
    compact_list = itemgetter(*list_rnd)(compact_list)
    bound_list = itemgetter(*list_rnd)(bound_list)

    # output to static HTML file
    output_file("./figures/plot_comparison.html")

    p = figure(plot_width=700, plot_height=700, title="Algorithms depending on the line of sight")

    # add a circle renderer with a size, color, and alpha
    p.circle(direct_list, dij_list, size=5, color="orange", alpha=0.5)
    p.circle(direct_list, compact_list, size=5, color="navy", alpha=0.5)
    p.circle(direct_list, bound_list, size=5, color="red", alpha=0.5)

    # labels
    p.yaxis.axis_label = "Total Routing Distance (km)"
    p.xaxis.axis_label = "Line of sight between pairs of satellites (km)"

    # show the results
    show(p)

#plot_routing_table_sizes(NUMBER_OF_NODES)
#plot_average_routing_tables(NUMBER_OF_NODES)
#plot_dij_vs_compact_large()