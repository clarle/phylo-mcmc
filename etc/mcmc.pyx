import numpy as np
cimport numpy as np

# Interfaces for functions from 'math.h' library

cdef extern from "math.h":
    double log(double x)
    double exp(double x)
    double sqrt(double x)

# Target distribution to be sampled from

cdef target_dist(double n):
    return 125*log(2+exp(n)/(1+exp(n))) - 74*log(1+exp(n)) + 35*n

# Standard normal pdf for proposal

cdef sdnorm(double z):
    return (-z*z/2.) - log(sqrt(2*np.pi))

# Random walk Metropolis-Hastings algorithm

def rand_walk_metropolis(int draws):
    cdef int i
    cdef double accepted = 0.
    cdef double current = 0.
    cdef double candidate, prob, ratio

    cdef np.ndarray[double, ndim=1] proposal = np.random.laplace(0., 1., draws)
    cdef np.ndarray[double, ndim=1] target = np.empty(draws)
    cdef np.ndarray[double, ndim=1] rand_draw = np.random.rand(draws)

    target[0] = current

    for i in xrange(1, draws):
        candidate = current + proposal[i]
        prob = min([1., exp(target_dist(candidate) - target_dist(current))])
        
        if rand_draw[i] < prob:
            current = candidate
            accepted += 1.
        
        target[i] = current

    ratio = accepted/draws
    
    return target, ratio

def indep_metropolis(int draws):
    cdef int i
    cdef double accepted = 0.
    cdef double current = 0.
    cdef double candidate, ratio, prob

    cdef np.ndarray[double, ndim=1] proposal = np.random.uniform(-1, 1, draws)
    cdef np.ndarray[double, ndim=1] target = np.empty(draws)
    cdef np.ndarray[double, ndim=1] rand_draw = np.random.rand(draws)

    target[0] = current

    for i in xrange(1, draws):
        candidate = proposal[i]
        prob = min(1., exp((target_dist(candidate) + sdnorm(candidate)) - (target_dist(current) + sdnorm(current))))

        if rand_draw[i] < prob:
            current = candidate
            accepted += 1

        target[i] = current

    ratio = accepted/draws

    return target, ratio
