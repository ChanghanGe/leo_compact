import os
import matplotlib.pyplot as plot
import matplotlib.ticker as mtick
from matplotlib.pyplot import figure

frequencies = []
for year in range(1, 20):
    if year < 10:
        dataFile = open('results/stage3/path-data-all-200' + str(year) + '.dat', 'r')
        line = dataFile.readline()
        line = dataFile.readline()
        line = dataFile.readline()
        yearEntry = [0, 0]
        totalPaths = 0
        while line != '':
            line = line.replace("'", '')
            line = line.replace(':', '')
            line = line.replace('\n', '')
            line = line.split()
            yearEntry.append(int(line[1]))
            totalPaths = totalPaths + float(line[1])
            line = dataFile.readline()
        for x in range(len(yearEntry)):
            if yearEntry[x] > 0.0:
                yearEntry[x] = 100.0/float(totalPaths) * float(yearEntry[x])
        frequencies.append(yearEntry)
    else:
        dataFile = open('results/stage3/path-data-all-20' + str(year) + '.dat', 'r')
        line = dataFile.readline()
        line = dataFile.readline()
        line = dataFile.readline()
        yearEntry = [0, 0]
        totalPaths = 0
        while line != '':
            line = line.replace("'", '')
            line = line.replace(':', '')
            line = line.replace('\n', '')
            line = line.split()
            yearEntry.append(float(line[1]))
            totalPaths = totalPaths + float(line[1])
            line = dataFile.readline()
        for x in range(len(yearEntry)):
            if yearEntry[x] > 0.0:
                yearEntry[x] = 100.0/float(totalPaths) * float(yearEntry[x])
        frequencies.append(yearEntry)

xAxis = range(20)
graph = figure().gca()

for year in range(19):
    while len(frequencies[year]) < 20:
        frequencies[year].append(0.0)
    if year < 9:
        plot.plot(xAxis, frequencies[year], label='200' + str(year+1))
    else:
        plot.plot(xAxis, frequencies[year], label='20' + str(year+1))
plot.xlim(0, 20)
graph.xaxis.set_major_locator(mtick.MaxNLocator(integer=True))
graph.yaxis.set_major_locator(mtick.MaxNLocator(integer=True))
plot.legend(prop={'size': 6})
plot.title('Path Length Frequencies')
plot.xlabel('Path Length')
plot.ylabel('Percentage of paths')

if not os.path.exists('./dissertation/images'):
    os.makedirs('./dissertation/images')

plot.savefig('dissertation/images/path-frequencies.pdf')
