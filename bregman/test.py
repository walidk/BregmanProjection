__author__ = 'Walid Krichene walid@eecs.berkeley.edu'
import math
import projections as p

# Compare the solution to the exponential projection using the different methods.
epsilon = .05

def expPotential(u):
    return math.exp(u - 1) - epsilon

def expInversePotential(u):
    return 1 + math.log(u+epsilon)

x = [1/2, 1/2]
g = [1.5, 1]

x1 = p.potentialProjection(expPotential, expInversePotential)(x, g, .00001)
x2 = p.expProjectionSort(epsilon, x, g)
x3 = p.expQuickProjection(epsilon, x, g)

print(sum(abs(xi - yi) for (xi, yi) in zip(x1, x3)))
print(x1, x2, x3)

