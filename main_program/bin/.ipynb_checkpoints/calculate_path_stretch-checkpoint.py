import sys
import networkx as nx
import random
import os

#import landmark set

dataFile = open('./results/stage2/kshells-' + sys.argv[1] + '.dat', 'r')
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

#import graph
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

# import node clusters
dataFile = open('./results/stage3/cluster-' + sys.argv[1] + '.dat', 'r')
line = dataFile.readline()
line = dataFile.readline()
while line != '':
    line = line.replace('\n', '')
    line = line.replace(',', '')
    line = line.replace("'", '')
    line = line.split(':')

    nodeId = line[0]
    node = graph.nodes()[nodeId]
    node['cluster'] = {}

    if len(line) > 1 and line[1] != ' ':
        line = line[1]
        line = line.split()
        for l in range(0, len(line), 2):
            line[l] = line[l].replace('[', '')
            line[l+1] = line[l+1].replace(']', '')
            line[l] = line[l].replace(',', '')
            line[l+1] = line[l+1].replace(',', '')

            node['cluster'][line[l]] = line[l+1]

    line = dataFile.readline()
dataFile.close()

#choose n/100 nodes at random
n = int(graph.number_of_nodes() / 100) + 1
chosenNodes = {}

while len(chosenNodes.keys()) < n:
    node = random.randint(0, graph.number_of_nodes()-1)
    if graph.nodes().keys()[node] not in chosenNodes.keys():
        chosenNodes[graph.nodes().keys()[node]] = {}

#for each chosen node calculate landmark paths and set closest landmark
for nodeId in graph.nodes().keys():
    paths = nx.algorithms.shortest_path(graph, nodeId)

    minLandmark = ''
    if nodeId == landmarkSet[0]:
        minLandmark = landmarkSet[1]
    else:
        minLandmark = landmarkSet[0]
    landmarkPath = paths[minLandmark]

    for landmark in landmarkSet:
        if len(paths[landmark]) < len(landmarkPath):
            minLandmark = landmark
            landmarkPath = paths[landmark]

    graph.nodes()[nodeId]['closestLandmark'] = [minLandmark, landmarkPath]

totalPathStretch = 0.0
totalPaths = 0

#for each chosen node calculate shortest paths
for nodeId in chosenNodes.keys():
    paths = nx.algorithms.shortest_path(graph, nodeId)
    #start by finding the path to the destination's closest landmark

    for destinationId in graph.nodes():
        if nodeId != destinationId:
            destination = graph.nodes()[destinationId]
            destinationLandmark = graph.nodes()[destinationId]['closestLandmark'][0]
            pathToLandmark = paths[destinationLandmark]
            clusterFound = False
            pathLength = 0
            #check if any nodes on the path to the landmark have the destination in their cluster
            for hop in pathToLandmark:
                if destinationId in graph.nodes()[hop]['cluster'].keys():

                    clusterFound = True
                    pathLength = pathToLandmark.index(hop)
                    pathFromClusterNode = nx.shortest_path(graph, hop, destinationId)
                    pathLength = pathLength + len(pathFromClusterNode)

                    pathStretch = 100/float(len(paths[destinationId])) * float(pathLength)
                    shortestPath = paths[destinationId]
                    averagePathStretch = totalPathStretch + pathStretch
                    totalPaths = totalPaths + 1

                    break

            #if no cluster shortcut is found, use the shortest path
            if not clusterFound:
                pathLength = len(pathToLandmark) + len(graph.nodes()[destinationId]['closestLandmark'][1]) -1

                pathStretch = (100/float(len(paths[destinationId]))) * float(pathLength)
                totalPathStretch = totalPathStretch + pathStretch
                totalPaths = totalPaths + 1

print round((float(totalPathStretch) / float(totalPaths)) /100, 2)

if not os.path.exists('./results/stage4'):
    os.makedirs('./results/stage4')

saveFile = open('results/stage4/path_stretch-' + sys.argv[1] + '.dat', 'w+')
saveFile.write(str(round((float(totalPathStretch) / float(totalPaths)) /100, 2)))