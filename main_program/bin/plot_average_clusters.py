import os
import matplotlib.pyplot as plot
from matplotlib.ticker import MaxNLocator
from matplotlib.pyplot import figure

# import node clusters
nodesWithClusters = []
averageClusterSize = []

for year in range(1, 20):
    if year < 10:
        dataFile = open('./results/stage3/cluster-all-200' + str(year) + '.dat', 'r')
        line = dataFile.readline()
        line = dataFile.readline()

        clusterCount = 0
        totalClusterSize = 0
        while line != '':
            line = line.replace('\n', '')
            line = line.replace(',', '')
            line = line.replace("'", '')
            line = line.split(':')

            if len(line) > 1 and line[1] != ' ':
                line = line[1]
                line = line.split()
                clusterCount = clusterCount + 1
                totalClusterSize = totalClusterSize + len(line)

            line = dataFile.readline()

        averageClusterSize.append(round(float(totalClusterSize) / float(clusterCount)))
        sizeFile = open('./results/stage1/data-all-200' + str(year) + '.dat', 'r')
        graphSize = sizeFile.readline()
        graphSize = graphSize.split()[3]
        sizeFile.close()
        n = round((100.0 / float(graphSize)) * float(clusterCount), 2)
        nodesWithClusters.append(n)

        dataFile.close()
    else:
        dataFile = open('./results/stage3/cluster-all-20' + str(year) + '.dat', 'r')
        line = dataFile.readline()
        line = dataFile.readline()

        clusterCount = 0
        totalClusterSize = 0
        while line != '':
            line = line.replace('\n', '')
            line = line.replace(',', '')
            line = line.replace("'", '')
            line = line.split(':')

            if len(line) > 1 and line[1] != ' ':
                line = line[1]
                line = line.split()
                clusterCount = clusterCount + 1
                totalClusterSize = totalClusterSize + len(line)

            line = dataFile.readline()

        averageClusterSize.append(round(float(totalClusterSize) / float(clusterCount)))
        sizeFile = open('./results/stage1/data-all-20' + str(year) + '.dat', 'r')
        graphSize = sizeFile.readline()
        graphSize = graphSize.split()[3]
        sizeFile.close()
        n = round((100.0 / float(graphSize)) * float(clusterCount), 2)
        nodesWithClusters.append(n)

        dataFile.close()


years = [2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]
graph = figure().gca()
plot.plot(years, averageClusterSize, label='Average Cluster Size')

plot.xlim(2000, 2020)
graph.xaxis.set_major_locator(MaxNLocator(integer=True))
plot.title('Average Cluster Size in Each Snapshot')
plot.xlabel('Year')
plot.ylabel('Average number of nodes in each cluster')

if not os.path.exists('./dissertation/images'):
    os.makedirs('./dissertation/images')

plot.savefig('dissertation/images/average-cluster-size.pdf')

plot.clf()


years = [2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]
graph = figure().gca()
plot.plot(years, nodesWithClusters, label='Nodes with Clusters')

plot.xlim(2000, 2020)
graph.xaxis.set_major_locator(MaxNLocator(integer=True))
plot.title('Percentage of nodes with clusters in each snapshot')
plot.xlabel('Year')
plot.ylabel('Percentage of nodes that have a cluster')

if not os.path.exists('./dissertation/images'):
    os.makedirs('./dissertation/images')

plot.savefig('dissertation/images/nodes-with-clusters.pdf')
