from mcmc import rand_walk_metropolis, indep_metropolis
from matplotlib.pylab import hist, ylabel, xlabel, title, suptitle, show

data, ratio = indep_metropolis(100000)
hist(data, bins=30, normed=1)
xlabel('x')
ylabel('Frequency')
title('Acceptance rate: %f' % ratio)
show()
