import numpy as np
cimport numpy as np

# Interfaces for functions from 'math.h' library

cdef extern from "math.h":
    double log(double x)
    double exp(double x)

# Target distribution to be sampled from

def target_dist(double n):
    return 125*log(2+exp(n)/(1+exp(n))) - 74*log(1+exp(n)) + 35*n

# Random walk Metropolis-Hastings algorithm

def rand_walk_metropolis(int draws):
    cdef int i
    cdef double current = 0.
    cdef double candidate, prob

    cdef np.ndarray[double, ndim=1] proposal = np.random.uniform(-1, 1, draws)
    cdef np.ndarray[double, ndim=1] target = np.empty(draws)
    cdef np.ndarray[double, ndim=1] uniform = np.random.uniform(0, 1, draws) 

    target[0] = current

    for i in xrange(1, draws):
        candidate = current + proposal[i]
        prob = min([1., target_dist(candidate)/target_dist(current)])
        
        if uniform[i] < prob:
            current = candidate
            target[i] = candidate
        else:
            target[i] = current

    return target
