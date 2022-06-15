import ephem
import csv
import math
from astropy import constants as const
import networkx as nx
import copy
import geometry_functions as geom


#NB: In the following code, all satellites are labeled by [orbitNumber][satelliteNumber]

#Modify the constants here in case we use a constellation with different parameters.
#---------------------------------------------------------

TCL_FILE_NAME = 'SpaceX_constellation1_1150km.tcl'
ALTITUDE = 1150
SATELLITES_PER_ORBIT = 50
NUMBER_OF_ORBITS = 32

#---------------------------------------------------------
EPOCH = '2019/03/18'
OBSERVATION_DATE = '2019/03/19 08:00:00'
#---------------------------------------------------------

EARTH_RADIUS = 6371
MAX_DISTANCE_BETWEEN_SATS = 2*(ALTITUDE+EARTH_RADIUS)
MEAN_MOTION = 7528.14 #this is always constant from a circular orbit for a given altitude
ECCENTRICITY_ADJUSTMENT = math.pow(10, -20)
INCLINATION = 53.


NUMBER_OF_LEVELS = 3
links_number = 5

#Creates a Pyephem satellite constellation with the data retrieved from SaVi
def constellationFromSaVi():
    DELETE_FROM_BEGINING = 3
    DELETE_FROM_END = 2
    constellation = list()
    with open(TCL_FILE_NAME, 'r') as tclfile:
        #we read from the TCL file and obtain a list containing SATELLITES_PER_ORBIT*NUMBER_OF_ORBITS
        #Attention : these elements are each one list(string) of size 1
        spaceX_SaVi = list(csv.reader(tclfile, delimiter='\n'))

        #we delete the lines that do not give any infomation for the satellites
        for i in range (0, DELETE_FROM_BEGINING):
            spaceX_SaVi.pop(0)
        for i in range (0, DELETE_FROM_END):
            spaceX_SaVi.pop(len(spaceX_SaVi)-1)
        if len(spaceX_SaVi) != SATELLITES_PER_ORBIT*NUMBER_OF_ORBITS:
            raise Exception('The total number of satellites is not correct')



    for i in range (0, NUMBER_OF_ORBITS):
        orbit_i = list()
        for j in range (0, SATELLITES_PER_ORBIT):
            SaVi_line = spaceX_SaVi.pop(0)[0].split()
            for i in range (0, 2):
                SaVi_line.pop(0)
            to_add_sat = ephem.EarthSatellite()
            to_add_sat._epoch = EPOCH #ok
            to_add_sat._n = geom.semi_major_to_mean_motion(float(SaVi_line[0])) #ok
            to_add_sat._e = float(SaVi_line[1]) + ECCENTRICITY_ADJUSTMENT #ok
            to_add_sat._inc = float(SaVi_line[2]) #ok
            to_add_sat._raan = float(SaVi_line[3]) #ok
            to_add_sat._ap = float(SaVi_line[4]) #ok
            to_add_sat._M = geom.time_to_periapsis_to_mean_anomaly(float(SaVi_line[5]),to_add_sat._n) #ok
            to_add_sat._drag = 0.0
            orbit_i.append(to_add_sat)
        constellation.append(orbit_i)
    return constellation


#Simulates the position of these satellites at a given "time" and returns it as a matrix
def positionsAtTime(input_constellation, time):
    all_positions = list()
    for orbit in input_constellation:
        positions_satellites_orbit = list()
        for satellite in orbit:
            satellite.compute(ephem.date(time))
            positions_satellites_orbit.append((satellite.sublong, satellite.sublat))
        all_positions.append(positions_satellites_orbit)
    return all_positions



#This computes the distances between each satellite and each other satellite. They are labeled by its [orbitNumber][satelliteNumber]
#Example : the distance btw Satellite[3][42] and Satellite[12][8] is given by distances(positions)[3][42][12][8]
def distances(positions):
    all_distances = list()
    for orbit1 in positions:
        all_distances_for_orbit = list()
        for satellite1 in orbit1:
            tempList2 = list()
            for orbit2 in positions:
                tempList = list()
                for satellite2 in orbit2:
                    distance = geom.haversine(satellite1, satellite2, ALTITUDE)
                    #TODO : Compute how many LEOs in a 1800km radius around a satellite/Add a flag if two satellites are in this radius
                    tempList.append(distance)
                tempList2.append(tempList)
            all_distances_for_orbit.append(tempList2)
        all_distances.append(all_distances_for_orbit)
    return all_distances


#input : orbit_number into the mathematical interval [0, NUMBER_OF_ORBITS-1]
#input : sat_number into the mathematical interval [0, SATELLITES_PER_ORBIT-1]
#output : node number in range [0, NUMBER_OF_ORBITS*SATELLITES_PER_ORBIT-1]
def get_node_number(orbit_number, sat_number):
    if orbit_number >= NUMBER_OF_ORBITS or orbit_number < 0 :
        raise AttributeError('The orbit number is out of the interval [0, NUMBER_OF_ORBITS-1]')
    if sat_number >= SATELLITES_PER_ORBIT or sat_number < 0:
        raise AttributeError('The sat number is out of the interval [0, SATELLITES_PER_ORBIT-1]')
    return SATELLITES_PER_ORBIT*orbit_number+sat_number


#input : node number in range [0, NUMBER_OF_ORBITS*SATELLITES_PER_ORBIT-1]
#output : orbit_number into the mathematical interval [0, NUMBER_OF_ORBITS-1]
#output : sat_number into the mathematical interval [0, SATELLITES_PER_ORBIT-1]
def get_orbit_and_sat_number(node_number):
    if node_number > NUMBER_OF_ORBITS*SATELLITES_PER_ORBIT-1 or node_number<0:
        raise AttributeError('The orbit number is out of the interval [0, NUMBER_OF_ORBITS*SATELLITES_PER_ORBIT-1]')
    sat_number = node_number%SATELLITES_PER_ORBIT
    orbit_number = (node_number-sat_number)/SATELLITES_PER_ORBIT
    return (int(orbit_number), sat_number)



#input : distances between satellites in a constellation
#output : nx graph with each vertex corresponding to a satellite and each edge corresponding to the distance between two satellites
#Note : this does not select the five closest links
def graph_full_mesh_from_constellation(constellation):
    n_orbits = len(constellation)
    n_satellites = len(constellation[0])
    new_graph = nx.path_graph(n_orbits*n_satellites)
    for orbit1, first_level in enumerate(constellation):
        for satellite1, second_level in enumerate(first_level):
            for orbit2, third_level in enumerate(second_level):
                for satellite2, distance in enumerate(third_level):
                    new_graph.add_edge(get_node_number(orbit1, satellite1), get_node_number(orbit2, satellite2),
                                    weight=distance)
    return new_graph


# input : distances between satellites in a constellation
# output : nx graph with each vertex corresponding to a satellite and each edge corresponding to the distance between two satellites
# Note : here each satellite has exactly five links
def graph_five_links_from_constellation(input_constellation, links_number=5):
    constellation = copy.deepcopy(input_constellation)  # we use this deepcopy to allow modification inside constellation

    def add_if_not_present(orbit1, orbit2, satellite1, satellite2, n_links, in_graph):
        node1 = get_node_number(orbit1, satellite1)
        node2 = get_node_number(orbit2, satellite2)
        if node1 == node2:
            raise Exception("The two nodes shouldn't be equal")
        if not(in_graph.has_edge(node1, node2)) and n_links[node1] < links_number and n_links[node2] < links_number:
            distance = constellation[orbit1][satellite1][orbit2][satellite2]
            constellation[orbit1][satellite1][orbit2][satellite2] = MAX_DISTANCE_BETWEEN_SATS
            # will be useful for the find_add_min step
            in_graph.add_edge(node1, node2, weight=distance)
            n_links[node1] += 1
            n_links[node2] += 1

    def add_min(target_orbit, target_sat, distances_from_node, n_links, in_graph):
        closest_orbit, closest_sat , min_dist = None, None, MAX_DISTANCE_BETWEEN_SATS
        for i, orbits in enumerate(distances_from_node):
            for j, distance in enumerate(orbits):
                if i != target_orbit and j != target_sat and distance < min_dist:
                    closest_orbit = i
                    closest_sat = j
                    min_dist = distance
        node1 = get_node_number(target_orbit, target_sat)
        node2 = get_node_number(closest_orbit, closest_sat)
        if n_links[node1] < links_number and n_links[node2] < links_number:
            distance = constellation[target_orbit][target_sat][closest_orbit][closest_sat]
            in_graph.add_edge(node1, node2, weight=distance)
            n_links[node1] += 1
            n_links[node2] += 1

    n_orbits = len(constellation)
    n_satellites = len(constellation[0])
    new_graph = nx.Graph()
    new_graph.add_nodes_from([i for i in range (0, n_orbits* n_satellites)])
    number_of_links = [0 for _ in range (0, n_orbits*n_satellites)]
    for orbit, first_level in enumerate(constellation):
        for satellite, second_level in enumerate(first_level):
            next_sat_on_orbit = (satellite+1) % n_satellites
            add_if_not_present(orbit, orbit, satellite, next_sat_on_orbit, number_of_links, new_graph)
            prev_sat_on_orbit = (satellite-1) % n_satellites
            add_if_not_present(orbit, orbit, satellite, prev_sat_on_orbit, number_of_links, new_graph)
            next_orbit = (orbit+1) % n_orbits
            add_if_not_present(orbit, next_orbit, satellite, satellite, number_of_links, new_graph)
            prev_orbit = (orbit-1) % n_orbits
            add_if_not_present(orbit, prev_orbit, satellite, satellite, number_of_links, new_graph)
            add_min(orbit, satellite, second_level, number_of_links, new_graph)

    return new_graph


# Here are the functions to call in order to create the desired graphs ------------------------------------------------

#Creates the nx_graph for our constellation
def create_spaceX_graph():
    spaceX_constellation = constellationFromSaVi()
    spaceX_positions = positionsAtTime(spaceX_constellation, OBSERVATION_DATE)
    all_distances = distances(spaceX_positions)
    graph = graph_five_links_from_constellation(all_distances)
    return graph, all_distances

#Creates a nx_graph corresponding to our constellation, with edges between all the nodes
#The weights correspond to the line of sight between two satellites
def create_full_mesh():
    spaceX_constellation = constellationFromSaVi()
    spaceX_positions = positionsAtTime(spaceX_constellation, OBSERVATION_DATE)
    all_distances = distances(spaceX_positions)
    full_mesh_graph = graph_full_mesh_from_constellation(all_distances)
    return full_mesh_graph


# 3 orbits and 3 satellites per orbit, distributed in a symetric way
def create_small_sat_graph():
    first_angle = math.pi #lattitude and longitude varies between 0 and pi
    second_angle = math.pi
    positions = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    positions[0][0] = 0, 0
    positions[0][1] = first_angle/2, 0
    positions[0][2] = first_angle, 0
    positions[1][0] = 0, second_angle/2
    positions[1][1] = first_angle / 2, second_angle/2
    positions[1][2] = first_angle, second_angle/2
    positions[2][0] = 0, second_angle
    positions[2][1] = first_angle / 2, second_angle
    positions[2][2] = first_angle, second_angle
    all_distances = distances(positions)
    graph = nx.Graph()
    graph.add_nodes_from(range(0, 8))
    for orbit1 in range(0 , 3):
        for sat1 in range(0, 3):
            for orbit2 in range(0, 3):
                for sat2 in range(0, 3):
                    first_node = 3*orbit1+sat1
                    second_node = 3*orbit2+sat2
                    if not graph.has_edge(first_node, second_node):
                        graph.add_edge(first_node, second_node, weight=all_distances[orbit1][sat1][orbit2][sat2])
    return graph, all_distances









