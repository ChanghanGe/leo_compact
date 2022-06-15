import math
from astropy import constants as const
import numpy as np

# haversine formula found on https://www.geeksforgeeks.org/haversine-formula-to-find-distance-between-two-points-on-a-sphere/
# slightly modified function to stick with our data
# INPUT : positions of satellites in radians, altitude in km
# OUTPUT : Approximate distance between two satellites in kilometers, as a straight line
def haversine(positions1, positions2, altitude1, altitude2 = None):
    lat1 = positions1[0]
    lon1 = positions1[1]
    lat2 = positions2[0]
    lon2 = positions2[1]

    # distance between latitudes
    # and longitudes
    dLat = (lat2 - lat1)
    dLon = (lon2 - lon1)


    def hav(angle):
        return (1 - math.cos(angle))/2

    # apply formulae
    haver_angle = 2 * math.asin(math.sqrt(hav(dLat)+(math.cos(lat1)*math.cos(lat2)*hav(dLon))))

    if altitude2 == None:
        r = (const.R_earth.value / 1000) + altitude1
        distance = r * haver_angle
        chord = 2*r*math.sin((0.5*distance/r))
        return chord

    else:
        r1 = (const.R_earth.value / 1000) + altitude1
        r2 = (const.R_earth.value / 1000) + altitude2

        a = np.min([r1,r2]) * math.sin(haver_angle)
        b = np.abs(r1-r2)

        return math.sqrt(a**2 + b**2)


def haversine_pathLabel(positions1, positions2, altitude):
    lat1 = positions1[0]
    lon1 = positions1[1]
    lat2 = positions2[0]
    lon2 = positions2[1]

    # distance between latitudes
    # and longitudes
    dLat = (lat2 - lat1)
    dLon = (lon2 - lon1)


    def hav(angle):
        return (1 - math.cos(angle))/2

    #Earth radius
    R = const.R_earth.value/1000

    # apply formulae
    r = R + altitude
    distance = 2*r*math.asin(math.sqrt(hav(dLat)+(math.cos(lat1)*math.cos(lat2)*hav(dLon))))
    chord = 2*r*math.sin((0.5*distance/r))

    vert = r*math.cos((0.5*distance/r))
    pathLabel = True
    if vert < R:
        pathLabel = False
    return chord, pathLabel


SECONDS_IN_DAY = 86400
#Input : semi_major_axis_length in km
#Output : mean_motion in revolutions/day (float)
def semi_major_to_mean_motion(axis_length):
    return SECONDS_IN_DAY*math.sqrt(const.G.value*(const.M_earth.value)/math.pow(axis_length*1000,3))/(2*math.pi)


#Input : mean motion (in rotations/day) and time to periapsis
#Output : mean anomaly in degrees (string format)
# (Mean Anomaly / 360) = (time to periapsis / duration of orbit)
def time_to_periapsis_to_mean_anomaly(time_to_periapsis, mean_motion):
    time_to_periapsis = -time_to_periapsis
    return 360*time_to_periapsis/(SECONDS_IN_DAY/mean_motion)


def fibonacci_sphere(samples=1000):

    points = []
    points_rad = []
    phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append((x, y, z))
        if x == 0:
            if y > 0:
                temp = math.atan(math.inf)
            else:
                temp = math.atan(-math.inf)
        else:
            temp = math.atan(y/x)
            
        points_rad.append((temp,np.angle(np.exp(1j*theta))))

    return points, points_rad