import matplotlib.pyplot as plot
from matplotlib.ticker import MaxNLocator
from matplotlib.pyplot import figure
import os

#read in data


all_edges = []
oldest_edges = []

for year in range(1, 20):
    if year <10:
        dataFile = open('results/stage1/data-all-200' + str(year) +'.dat', 'r')
        line = dataFile.readline()
        line = dataFile.readline()
        line = line.split()
        all_edges.append(int(line[3]))
        dataFile.close()

        dataFile = open('results/stage1/data-oldest-200' + str(year) +'.dat', 'r')
        line = dataFile.readline()
        line = dataFile.readline()
        line = line.split()
        oldest_edges.append(int(line[3]))
        dataFile.close()

    else:
        dataFile = open('results/stage1/data-all-20' + str(year) + '.dat', 'r')
        line = dataFile.readline()
        line = dataFile.readline()
        line = line.split()
        all_edges.append(int(line[3]))
        dataFile.close()

        dataFile = open('results/stage1/data-oldest-20' + str(year) +'.dat', 'r')
        line = dataFile.readline()
        line = dataFile.readline()
        line = line.split()
        oldest_edges.append(int(line[3]))
        dataFile.close()

years = [2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]
graph = figure().gca()
plot.plot(years, all_edges, label='All Collectors')
plot.plot(years, oldest_edges, label='Oldest Collectors')

plot.xlim(2000, 2020)
graph.xaxis.set_major_locator(MaxNLocator(integer=True))
plot.legend()
plot.title('Number of Edges')
plot.xlabel('Year')
plot.ylabel('Number of edges in the snapshot')

if not os.path.exists('./dissertation/images'):
    os.makedirs('./dissertation/images')

plot.savefig('dissertation/images/num_edges.pdf')
