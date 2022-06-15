
set_percentage = []

totalPercentage = 0.0
for year in range(1, 20):

    if year < 10:
        dataFile = open('./results/stage1/data-all-200' + str(year) + '.dat', 'r')
    else:
        dataFile = open('./results/stage1/data-all-20' + str(year) + '.dat', 'r')

    totalNumber = float(dataFile.readline().split()[3])
    dataFile.readline()
    setNumber = float(dataFile.readline().split()[2])
    set_percentage.append((100/totalNumber) * setNumber)
    totalPercentage = totalPercentage + (100/totalNumber) * setNumber

print min(set_percentage)
print max(set_percentage)
print totalPercentage / 19.0