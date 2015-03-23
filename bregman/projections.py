__author__ = 'Walid Krichene walid@eecs.berkeley.edu'
import math
import random


def identity(x):
    """Identity function."""
    return x


def potentialProjection(potential = identity, inversePotential = identity):
    """Creates a method which computes the Bregman projection of a vector x given a gradient vector g, using a
    bisection method.
    One must specify the potential function (default is identity) and its inverse.
    Returns a method projection(x, g, precision) which takes as input
    - the current iterate x
    - the gradient vector (scaled by the step size) g
    - the desired precision (in the ell one norm)
    """
    def projection(x, g, precision):
        d = len(x)
        y = max(inversePotential(xi) - gi for (xi, gi) in zip(x, g))
        nuU = inversePotential(1) - y
        nuL = inversePotential(1/d) - y
        def xBar(nu):
            return (max(0, potential(inversePotential(xi) - gi + nu)) for (xi,gi) in zip(x, g))
        while(sum(a - b for (a, b) in zip(xBar(nuU), xBar(nuL))) > precision):
            nuM = (nuU+nuL)/2
            if(sum(xBar(nuM)) > 1):
                nuU = nuM
            else:
                nuL = nuM
        return list(xBar(nuU))
    return projection

def expProjectionSort(epsilon, x, g):
    """Computes the Bregman projection, with exponential potential, of a vector x given a gradient vector g, using a
    sorting method. The complexity of this method is O(d log d), where d is the size of x.
    Takes as input
    - the parameter epsilon of the exponential potential.
    - the current iterate x
    - the gradient vector (scaled by the step size) g
    """
    d = len(x)
    y = list((xi+epsilon)*math.exp(-gi) for (xi, gi) in zip(x, g))
    yy = sorted(y)
    S = sum(yy)
    j = 0
    while((1+epsilon*(d-j))*yy[j]/S - epsilon <= 0):
        S -= yy[j]
        j += 1
    return list(max(0, -epsilon+(1+epsilon*(d-j))*yi/S) for yi in y)


def expQuickProjection(epsilon, x, g):
    """Computes the Bregman projection, with exponential potential, of a vector x given a gradient vector g, using a
    randomized pivot method. The expected complexity of this method is O(d), where d is the size of x.
    Takes as input
    - the parameter epsilon of the exponential potential.
    - the current iterate x
    - the gradient vector (scaled by the step size) g
    """
    d = len(x)
    y = list((xi+epsilon)*math.exp(-gi) for (xi, gi) in zip(x, g))
    J = range(0,d)
    S = 0
    C = 0
    s = d+1
    while(len(J) > 0):
        j = random.choice(J)
        pivot = y[j]
        JP = list(i for i in J if y[i] >= pivot)
        JM = list(i for i in J if y[i] < pivot)
        CP = len(JP)
        SP = sum(y[i] for i in JP)
        gamma = (1+epsilon*(C+CP)*pivot - epsilon*(S+SP))
        if(gamma > 0):
            J = JM
            S += SP
            C += CP
            s = j
        else:
            J = JP
    Z = S/(1+epsilon*C)
    return list(max(0, -epsilon+yi/Z) for yi in y)
