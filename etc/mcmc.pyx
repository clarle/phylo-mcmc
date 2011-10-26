import numpy as np
import numdifftools as nd
from numpy import matrix
from numpy.linalg import inv, cholesky
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

cdef double log_norm(double n):
    return log(exp(-n**2/2.)/sqrt(2*np.pi))

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

    cdef np.ndarray[double, ndim=1] proposal = np.random.normal(0., 1, draws)
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

def rw_matrix_metropolis(log_post, var, scale, start, int draws, data=None):
    cdef np.ndarray[double, ndim=2] target = np.zeros([draws, start.size])
    cdef np.ndarray[double, ndim=2] b = matrix.getT(matrix(start))
    cdef np.ndarray[double, ndim=1] rand_draw = np.random.rand(draws)
    cdef np.ndarray[double, ndim=2] lower_tri = cholesky(matrix(var))

    cdef np.ndarray[double, ndim=2] bc
    cdef double cur, ratio
    cdef double can = log_post(start, data)
    cdef double accepted = 0.
  
    for i in xrange(0, draws):
        bc = b + (scale * matrix.getT(lower_tri)) * matrix.getT(matrix(np.random.normal(0., 1, start.size)))
        cur = log_post(bc, data)
        prob = exp(can - cur)

        if rand_draw[i] < prob:
            cur = can
            b = bc
            accepted += 1.

        target[i][0] = b[0]
        target[i][1] = b[1]

    ratio = accepted/draws

    return target, ratio

def laplace(log_post, mode, data=None):
    mode = fmin(log_post, mode, args=data)
    Hfun = nd.Hessian(lambda x: -log_post(x, data))
    var = -inv(Hfun(mode))
    return mode, var
