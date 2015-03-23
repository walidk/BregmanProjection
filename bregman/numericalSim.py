__author__ = 'Walid Krichene walid@eecs.berkeley.edu'

import projections as p
import time
import numpy

# Test for increasing dimensions
f = open('execution_times', 'w')

def generateGradientVector(d):
    return numpy.random.normal(0, 1/d, d)

avgSortProjectTimes = []
avgQuickProjectTimes = []

nb_sim = 30
epsilon = .1
ds = [100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000, 10000000, 30000000]
# ds = [200000, 400000, 600000, 800000, 1000000, 1200000, 1400000, 1600000, 1800000, 2000000, 2200000, 2400000, 2600000, 3000000]
# ds = [100, 300, 1000, 3000, 10000, 30000]
ds = []
for d in ds:
    print("d={}".format(d))
    x = [1/d]*d
    sortTimes = []
    quickTimes = []
    for sim in range(0, nb_sim):
        # generate a random vector to project
        g = generateGradientVector(d)
        # method 1: sort
        startTime = time.time()
        xx = p.expProjectionSort(epsilon, x, g)
        sortTimes.append(time.time() - startTime)
        # method 2: quick project
        startTime = time.time()
        xx = p.expQuickProjection(epsilon, x, g)
        quickTimes.append(time.time() - startTime)
    avgSortProjectTimes.append(sum(sortTimes)/len(sortTimes))
    avgQuickProjectTimes.append(sum(quickTimes)/len(quickTimes))

f.write("d, t1, t2\n")
for (d, t1, t2) in zip(ds, avgSortProjectTimes, avgQuickProjectTimes):
    f.write("{}, {}, {}\n".format(d, t1, t2))

print(avgQuickProjectTimes)
print(avgSortProjectTimes)

