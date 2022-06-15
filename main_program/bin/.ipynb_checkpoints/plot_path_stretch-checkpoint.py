import matplotlib.pyplot as plot
from matplotlib.ticker import MaxNLocator
from matplotlib.pyplot import figure
import os

path_stretch = []

for year in range(1, 20):
    if year <10:
        dataFile = open('results/stage4/path_stretch-all-200' + str(year) +'.dat', 'r')
        line = dataFile.readline()
        path_stretch.append(float(line))
        dataFile.close()

    else:
        dataFile = open('results/stage4/path_stretch-all-20' + str(year) + '.dat', 'r')
        line = dataFile.readline()
        path_stretch.append(float(line))
        dataFile.close()


years = [2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]
graph = figure().gca()
plot.plot(years, path_stretch, label='Average Path Stretch')

plot.xlim(2000, 2020)
graph.xaxis.set_major_locator(MaxNLocator(integer=True))
plot.title('Average path stretch for each snapshot using TZ routing')
plot.xlabel('Year')
plot.ylabel('Average path stretch')

if not os.path.exists('./dissertation/images'):
    os.makedirs('./dissertation/images')

plot.savefig('dissertation/images/path-stretch.pdf')
