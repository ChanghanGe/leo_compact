import ephem
import csv
import math
from astropy import constants as const
import networkx as nx
import copy
import geometry_functions as geom
import numpy as np
import sys
from multiprocessing_functions import kernel_function_visible_sat

#Creates a Pyephem Fsatellite constellation with the data retrieved from SaVi
def constellationFromSaVi(file_name = TCL_FILE_NAME, OBSERVATION_DATE = EPOCH):
    DELETE_FROM_BEGINING = 0
    DELETE_FROM_END = 0
    constellation = list()
    with open(file_name, 'r') as tclfile:
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

    for i in range(NUMBER_OF_ORBITS):
        orbit_i = list()
        for j in range(SATELLITES_PER_ORBIT):
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
            to_add_sat._drag = 0
            orbit_i.append(to_add_sat)
            # print(float(SaVi_line[5]), geom.time_to_periapsis_to_mean_anomaly(float(SaVi_line[5]),to_add_sat._n))
        constellation.append(orbit_i)
    return constellation

def groundstationFromTable_single_gs(city, EPOCH = EPOCH, OBSERVATION_DATE = EPOCH):
    groundstation = ephem.city(city)
    groundstation.epoch = EPOCH #ok
    groundstation.date = ephem.date(OBSERVATION_DATE) #ok
    return groundstation

def kernel_function_visible_sat(city, OBSERVATION_DATE, FoV = '40'):
    cur_constellation = constellationFromSaVi(OBSERVATION_DATE=OBSERVATION_DATE)
    cur_groundstation = groundstationFromTable_single_gs(city, OBSERVATION_DATE=OBSERVATION_DATE)

    cur_visible_sats = []
    for orbit_id, orbit in enumerate(cur_constellation):
        for sat_id, satellite in enumerate(orbit):
            satellite.compute(cur_groundstation)
            if satellite.alt >= ephem.degrees(FoV):
                cur_visible_sats.append([orbit_id, sat_id])

    return cur_constellation, cur_groundstation, cur_visible_sats