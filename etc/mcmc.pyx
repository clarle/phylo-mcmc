import numpy as np
from scipy.stats import norm, t
from scipy.optimize import fmin
cimport numpy as np

# Interfaces for functions from 'math.h' library

cdef extern from "math.h":
    double log(double x)
    double exp(double x)
    double sqrt(double x)

# Target distribution to be sampled from

cdef double target_dist(double n):
    return 125*log(2+exp(n)/(1+exp(n))) - 74*log(1+exp(n)) + 35*n

# Log of normal PDF

def log_norm(double n):
    return log(norm.pdf(n))

# Modeling data with Cauchy errors

def cauchy_error_post(theta, data):
    mu = theta[0]
    sigma = exp(theta[1])
    neg_cauchy_err = lambda x: -log(t.pdf((x-mu)/sigma, 1)/sigma)
    logf = np.vectorize(neg_cauchy_err)
    return np.sum(logf(data))

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
        prob = min(1., exp(target_dist(candidate) - target_dist(current)))
        
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

    cdef np.ndarray[double, ndim=1] proposal = np.random.norm(-1, 1, draws)
    cdef np.ndarray[double, ndim=1] target = np.empty(draws)
    cdef np.ndarray[double, ndim=1] rand_draw = np.random.rand(draws)

    target[0] = current

    for i in xrange(1, draws):
        candidate = proposal[i]
        prob = min(1., exp((target_dist(candidate) + log_norm(candidate)) - (target_dist(current) + log_norm(current))))

        if rand_draw[i] < prob:
            current = candidate
            accepted += 1

        target[i] = current

    ratio = accepted/draws

    return target, ratio

def laplace(log_post, mode, data=()):
    mode = fmin(log_post, mode, args=data)
    return mode
