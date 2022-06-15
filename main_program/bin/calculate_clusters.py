import os
import sys
import datetime
import networkx as nx
import math

#import landmark set

dataFile = open('results/stage2/kshells-' + sys.argv[1] + '.dat', 'r')
line = dataFile.readline()
if int(sys.argv[2]) == 0:
    line = dataFile.readline()
elif int(sys.argv[2]) == 1:
    line = dataFile.readline()
    line = dataFile.readline()

line = line.replace('[', '')
line = line.replace(']', '')
line = line.replace(',', '')
line = line.replace("'", '')
line = line.replace('\n', '')
landmarkSet = line.split()
dataFile.close()

#import graph from the data file and set up the correct landmark nodes
graph = nx.Graph()
dataFile = open('./results/stage1/data-' + sys.argv[1] + '.dat', 'r')

line = dataFile.readline()
line = dataFile.readline()
line = dataFile.readline()

line = dataFile.readline()

while line != '':
    line = line.replace('(', '')
    line = line.replace(')', '')
    line = line.replace(',', '')
    line = line.replace("'", '')
    line.strip()
    line = line.split()

    node1 = line[0]
    node2 = line[1]

    graph.add_edge(node1, node2)

    line = dataFile.readline()

dataFile.close()

pathSet = {}
for landmark in landmarkSet:
    paths = nx.algorithms.shortest_path(graph, landmark)
    del paths[landmark]
    pathSet[landmark] = paths

for node in graph.nodes():
    graph.nodes[node]['paths'] = {}
    for landmark in landmarkSet:
        if node != landmark:
            graph.nodes()[node]['paths'][landmark] = pathSet[landmark][node]


pathLenTotal = 0
numPaths = 0
pathLenFreq = {}

#add a destination node to a node's cluster if the path to the destination is shorter than to any landmark node
for nodeId in graph.nodes:
    node = graph.nodes[nodeId]
    cluster = {}
    paths = nx.algorithms.shortest_paths.shortest_path(graph, nodeId)

    #find the closest landmark
    if nodeId == landmarkSet[0]:
        minLandmark = landmarkSet[1]
    else:
        minLandmark = landmarkSet[0]
    minDist = len(node['paths'][minLandmark])
    for landmark in landmarkSet:
        if nodeId != landmark:
            if len(node['paths'][landmark]) < minDist:
                minLandmark = landmark
                minDist = len(node['paths'][landmark])

    #if the path to the node is shorter through another node than a landmark, add the node to the other node's cluster
    for destinationId in graph.nodes:
        if nodeId != destinationId:

            pathLenTotal = pathLenTotal + len(paths[destinationId])
            numPaths = numPaths + 1
            if len(paths[destinationId]) not in pathLenFreq.keys():
                pathLenFreq[len(paths[destinationId])] = 1
            else:
                pathLenFreq[len(paths[destinationId])] = pathLenFreq[len(paths[destinationId])] + 1

            if len(paths[destinationId]) < minDist:
                destination = graph.nodes[destinationId]
                try:
                    destination['cluster'][nodeId] = paths[destinationId][-2]
                except:
                    destination['cluster'] = {}
                    destination['cluster'][nodeId] = paths[destinationId][-2]
    graph.nodes[nodeId]['cluster'] = cluster

#check that all clusters are under the size limit of 4 * sqrt(n log n)
violations = 0
n = graph.number_of_nodes()
limit = 4*math.sqrt(n * math.log(n))

total_cluster_size = 0
for nodeId in graph.nodes:
    total_cluster_size = total_cluster_size + len(graph.nodes()[nodeId]['cluster'])
    if len(graph.nodes[nodeId]['cluster']) > limit:
        violations = violations + 1

if not os.path.exists('./results/stage3'):
    os.makedirs('./results/stage3')

pathFile = open('results/stage3/path-data-' + sys.argv[1] + '.dat', 'w+')
pathFile.write('average path length: ' + str(round(float(pathLenTotal) / float(numPaths), 2)) + '\n')
averageCSize = round(float(total_cluster_size) / float(graph.number_of_nodes()))
pathFile.write('average cluster size: ' + str(averageCSize) + '\n')
for length in pathLenFreq.keys():
    pathFile.write(str(length) + ': ' + str(pathLenFreq[length]) + '\n')
pathFile.close()

saveFile = open('results/stage3/cluster-' + sys.argv[1] + '.dat', 'w+')
saveFile.write('violations: ' + str(violations) + '\n')

for nodeId in graph.nodes:
    node = graph.nodes[nodeId]
    output = (nodeId + ': ')

    for destId in node['cluster'].keys():
        output = output + '[' + destId + ', ' + node['cluster'][destId] + '], '
    saveFile.write(output + '\n')
saveFile.close()
