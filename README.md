# BregmanProjection

This provides a simple python implementation of the Bregman Projection with omega potentials described in [this paper](http://www.eecs.berkeley.edu/~walid/papers/CDC-efficient-projection.pdf), currently under review at the [IEEE Conference on Decision and Control (CDC)](http://www.cdc2015.ctrl.titech.ac.jp/).

Two implementations are provided:

## Approximate projections for omega potentials
Uses a bisection method to compute an approximate solution up to a precision `epsilon`. The complexity is `O(d log 1/epsilon)`, where `d` is the size of `x`.
Here is an example use (for the exponential projection)
```
import projections as p

def expPotential(u):
    return math.exp(u - 1) - epsilon

def expInversePotential(u):
    return 1 + math.log(u+epsilon)

expProjection = p.potentialProjection(expPotential, expInversePotential)

x = [1/2, 1/2]
g = [1.5, 1]
epsilon = 0.0001

x1 = expProjection(x, g, epsilon)
```
We first start by defining two functions: `expPotential` and `expInversePotential`, then we instantiate the potentialProjection with this potential.
Then to project the current iterate `x`, with gradient `g`, up to a precision `epsilon`, we call `expProjection(x, g, epsilon)`.


## Exact projections for exponential potentials

In the case of exponential projections, the exact solution can be computed. We provide two methods:

### Sort and project
This algorithm creates a vector of weights, sorts it, and finds a threshold to compute the projection. The complexity is `O(d log d)` where `d` is the size of `x`.
Example code:
```
import projections as p

x = [1/2, 1/2]
g = [1.5, 1]
epsilon = 0.0001

x1 = expProjectionSort(x, g, epsilon)
```

### Quick project
This algorithm uses a random pivot at each iteration, splits the list into two parts according to the pivot (the criterion is detailed in the paper), then inducts on one of the lists. The expected complexity is `O(d)`, where d is the size of `x`.
Example code
```
import projections as p

x = [1/2, 1/2]
g = [1.5, 1]
epsilon = 0.0001

x1 = expQuickProjection(x, g, epsilon)
```
