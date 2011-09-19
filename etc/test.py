from random_metro import rand_walk_metropolis
from matplotlib.pylab import hist, ylabel, xlabel, show

data = rand_walk_metropolis(100000)
hist(data, bins=30, normed=1)
xlabel('x')
ylabel('Frequency')
show()
