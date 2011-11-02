import numpy as np
import scipy as sp
from random import shuffle
from itertools import product
from scipy.integrate import quad
from models import HKY85
from tree import Tree
from utils import calc_nuc_sites, prod, four_category_dist

class MCMC:
    def __init__(self, freq, k, alpha, tau, filename):
        self.tau = tau
        self.model = HKY85(freq, k, alpha)
        self.proposal = np.random.normal(0., 1., draws)
        self.tree = Tree(data=filename, tau)
        self.old_topo = None
        self.new_topo = None
        self.matrix = self.init_fun_matrix() # 3D dictionary
        self.fun_list = []
        self.nuc_sum = None
        self.density = None

    def step(self):
        nodes = self.tree.nodes()
        
        # make sure we can't choose the root as the target
        rand = np.random.randint(1, len(nodes))
        
        # STEP 1: Choose a non-root node as a target for the MCMC
        target = nodes[rand]
        topo = self.generate_topology(target)

        # STEP 2: Select a possible time from the continuous density
        time = select_time(target, topo)

        # STEP 3: Select a possible nucleotide sequence
        sequence = select_sequence(target, topo, time)

        return target, topo, time, sequence

    def calc_mutation_prob(self, seq1, seq2, t):
        freqs = calc_nuc_sites(seq1, seq2)
        mutate_probs = []
        
        for change in freqs:
            start_nuc = change[0]
            end_nuc = change[1]
            prob = self.model[start_nuc][end_nuc](t)
            mutate_probs.append(prob)

        return prod(mutate_probs)

    def init_fun_matrix(self):
        fun_dict = {}
        generator = product("ATCG", repeat=3):
        
        for seq in generator:
            fun_dict[seq] = None

        return fun_dict

    def generate_topology(self, target):
        self.old_topo = [target.left, target.right, target.sibling]
        self.new_topo = shuffle(self.old_topo)
        return self.new_topo

    def calc_fun_matrix(self, target, topo):
        parent = target.parent
        for seq in self.matrix:
            for nuc in ('A', 'T', 'C', 'G'):
                p = lambda t: self.model.matrix[seq[0]][nuc](parent.time - t)
                c1 = lambda t: self.model.matrix[nuc][seq[1]](t - topo[1].time)
                c2 = lambda t: self.model.matrix[nuc][seq[2]](t - topo[2].time)
                self.fun_list.append(lambda t: p(t) * c1(t) * c2(t))
            self.matrix[seq] = lambda t: \
                self.fun_list[0](t) + self.fun_list[1](t) + self.fun_list[2](t) + self.fun_list[4](t)
                

    def prob_integral(self, target, topo):
        upper = target.time
        lower = max(topo[1].time, topo[2].time) # max of child times
        result = 0
        for prob in self.matrix:
            result += quad(self.matrix[prob], lower, upper)

        return result

    def sum_numer(self):
        return lambda t: sum(prob(t) for prob in self.matrix.values) 
    
    def calc_density(self, target, topo)
        self.matrix = self.calc_fun_matrix(target, topo)
        self.nuc_sum = self.sum_numer()
        density = lambda t: self.nuc_sum(t)/self.prob_integral(target, topo)
        return density

    def select_time(self, target, topo):
        self.density = self.calc_density(target, topo)
        sample = self.rejection_sample(self.density)
        return sample

    def rejection_sample(self, density):
        while True:
            uniform = np.random.uniform()
            beta = np.random.beta(0., 1., 1.5)
            if beta < density(uniform)
                return uniform

    def select_sequence(self, target, topo, t):
        nucs = ('A', 'T', 'C', 'G')
        freqs = {'A': 0., 'T': 0., 'C': 0., 'G': 0.}
		parent = target.parent
		for seq in self.matrix:
        	for nuc in nucs:
            	p = self.model.matrix[seq[0]][nuc](parent.time - t)
            	c1 = self.model.matrix[nuc][seq[1]](t - topo[1].time)
            	c2 = self.model.matrix[nuc][seq[2]](t - topo[2].time)
				freqs[nuc] += p * c1 * c2

        four_category_dist(1, freqs, nucs)[0]

    def transition_prob(self, target, topo, time)
        return self.prob_integral(target, topo) * \
            calc_mutation_prob(target.sequence, topo[0].sequence, time) *\
            len(self.tree.get_leaves())

	def alpha(self):
        return min(self.transition_prob(target, self.old_topo)/\
                   self.transition_prob(target, self.new_topo), 1.)

    def rand_walk_metropolis(draws):
        accepted = 0
        current = 0

        proposal = np.random.laplace(0., 1., draws)
        target = np.empty(draws)
        rand_draw = np.random.rand(draws)

        target[0] = current

        for i in xrange(1, draws):
            candidate = current + proposal[i]
            prob = self.alpha()

            if rand_draw[i] < prob:
                current = candidate
                accepted += 1.

            target[i] = current

        ratio = accepted/draws

        return target, ratio
