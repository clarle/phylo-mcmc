from mcmc import *
import numpy as np
from scipy.optimize import fmin

data = np.asarray((36, 13, 23, 6, 20, 12, 23, 93, 98, 91, 89, 100, 90, 95, 90, 87))
start = np.asarray((90, 2))

print cauchy_error_post(start, data)

print laplace(cauchy_error_post, start, data=(data,))
