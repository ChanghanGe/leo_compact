import networkx as nx
import sys
import os
import math

graph = nx.Graph()
filename = sys.argv[1]
dataFile = open('./results/stage1/data-' + filename + '.dat', 'r')

line = dataFile.readline()
line = dataFile.readline()
line = dataFile.readline()
line = dataFile.readline()

while line != '':

    line = line.replace('(', '')
    line = line.replace(')', '')
    line = line.replace("'", '')
    line = line.replace(',', '')
    line = line.split()

    graph.add_edge(int(line[0]), int(line[1]))
    
    line = dataFile.readline()
    
dataFile.close()
#graph imported

#calculate k-max shell
k_max = nx.algorithms.core.k_shell(graph)

shellSizes = []
nodesInShell = k_max.number_of_nodes()
shellSizes.append(nodesInShell)
x = 1
kc = []
#calculate k-core shells until a shell with at least root(n log(n)) nodes is found
while nodesInShell < math.sqrt(graph.number_of_nodes() * math.log(graph.number_of_nodes())):
    kc = nx.algorithms.core.k_shell(graph, x)
    nodesInShell = kc.number_of_nodes()
    shellSizes.append(nodesInShell)
    x = x + 1

#write results to file
if not os.path.exists('./results/stage2'):
    os.makedirs('./results/stage2')

saveFile = open('./results/stage2/kshells-' + sys.argv[1] + '.dat', 'w+')
saveFile.write(str(shellSizes) + '\n')
saveFile.write(str(k_max.nodes()) + '\n')
saveFile.write(str(kc.nodes()) + '\n')