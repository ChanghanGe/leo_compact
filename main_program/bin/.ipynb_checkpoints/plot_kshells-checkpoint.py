import matplotlib.pyplot as plot
from matplotlib.ticker import MaxNLocator
from matplotlib.pyplot import figure
import os
import math

shells = []

for year in range(1, 20):
    if year <10:
        dataFile = open('results/stage2/kshells-all-200' + str(year) +'.dat', 'r')
        totalFile = open('results/stage1/data-all-200' + str(year) +'.dat', 'r')
    else:
        dataFile = open('results/stage2/kshells-all-20' + str(year) + '.dat', 'r')
        totalFile = open('results/stage1/data-all-20' + str(year) + '.dat', 'r')

    line = dataFile.readline()
    line = line.replace('[', '')
    line = line.replace(']', '')
    line = line.replace(',', '')
    line = line.split()

    for x in range(len(line)):
        line[x] = float(line[x])

    totalSize = float(totalFile.readline().split()[3])
    limit = (4.0 * float(math.sqrt(totalSize* math.log(totalSize))))
    line.append(totalSize)
    line.append(limit)
    shells.append(line)

years = [2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]
graph = figure().gca()

kMax = []
for entry in shells:
    kMax.append(entry[0])
plot.plot(years, kMax, label = 'k-max')

k1 = []
for entry in shells:
    k1.append(entry[1])
plot.plot(years, k1, label = 'k-1')

t = []
for entry in shells:
    t.append(entry[2])
plot.plot(years, t, label = 'total')

s = []
for entry in shells:
    s.append(entry[3])
plot.plot(years, s, label = 's')

plot.xlim(2000, 2020)
graph.xaxis.set_major_locator(MaxNLocator(integer=True))
plot.legend()
plot.xlabel('Year')
plot.ylabel('Number of nodes')

if not os.path.exists('./dissertation/images'):
    os.makedirs('./dissertation/images')

plot.savefig('dissertation/images/k-shells.pdf')
