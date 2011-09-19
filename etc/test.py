from random_metro import *
from matplotlib.pylab import *

data = rand_walk_metropolis(10000)
hist(data)
show()
