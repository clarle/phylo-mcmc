from mcmc import *
import numpy as np
from scipy.optimize import fmin
from matplotlib.pylab import hist, ylabel, xlabel, title, suptitle, show

data = np.asarray((36., 13., 23., 6., 20., 12., 23., 93., 98., 91., 89., 100., 90., 95., 90., 87.))
start = np.asarray((10., 4.))

mode, var = laplace(cauchy_error_post, start, data=(data,))

print mode, var

data, ratio = rw_matrix_metropolis(cauchy_error_post, var, 2., start, 1000, data)

mean = np.mean(data, 0)
std = np.std(data, 0)

print data
print ratio
print mean
print std
