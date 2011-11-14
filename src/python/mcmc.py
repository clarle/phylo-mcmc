import numpy as np
import scipy as sp
from copy import deepcopy
from random import shuffle, choice
from itertools import product
from scipy.integrate import quad
from models import HKY85
from tree import Tree
from utils import calc_nuc_sites, prod, four_category_dist, unique_perms

class MCMC:
    def __init__(self, freq, k, alpha, tau, draws, filename):
        self.tau = tau
        self.model = HKY85(freq, k, alpha)
        self.draws = draws
        self.tree = Tree(data=filename, tau=tau)
        self.target = None
        self.old_topo = None
        self.new_topo = None
        self.matrix = self.init_fun_matrix() # 3D dictionary
        self.fun_list = []
        self.nuc_sum = None
        self.density = None

    def metropolis(self):
        accepted = 0
        first = self.tree
        tree_log = [first]
        print [node.sequence for node in tree_log[0].nodes()]
        print [node.time for node in tree_log[0].nodes()]
        rand_draw = np.random.rand(self.draws)

        for i in xrange(self.draws):
            target, topo, time, sequence = self.step()
            prob = self.alpha(target)

            if rand_draw[i] < prob:
                current = self.accept(target, topo, time, sequence)
                accepted += 1
                self.tree = current
                tree_log.append(current)
            else:
                copy = self.tree
                tree_log.append(copy)
            
            print [node.sequence for node in tree_log[i].nodes()]
            print [node.time for node in tree_log[i].nodes()]

        ratio = accepted/self.draws
        return tree_log, ratio

    def accept(self, target, topo, time, sequence):
        new = target

        topo[0].sibling = new
        topo[1].parent = new
        topo[2].parent = new
        topo[1].sibling = topo[2]
        topo[2].sibling = topo[1]

        new.sequence = sequence
        new.time = time
        new.sibling = topo[0]
        new.left = topo[1]
        new.right = topo[2]
       
        print "New Target:" + str(new.sequence)
        print "New Sibling:" + str(new.sibling.sequence)
        print "New Left:" + str(new.left.sequence)
        print "New Right:" + str(new.right.sequence)

        self.tree.update(target, new)
       
        # old_sib = topo[0]
        # new_sib = topo[0]
        # new_sib.sibling = new

        # self.tree.update(old_sib, new_sib)

        return self.tree

    def step(self):
        nodes = self.tree.get_internal()
        
        # make sure we can't choose the root as the target
        rand = np.random.randint(1, len(nodes))
        
        # STEP 1: Choose a non-root node as a target for the MCMC
        self.target = nodes[rand]
        print "Old Target:" + str(self.target.sequence)
        print "Old Sibling:" + str(self.target.sibling.sequence)
        print "Old Left:" + str(self.target.left.sequence)
        print "Old Right:" + str(self.target.right.sequence)
        topo = self.generate_topology(self.target)
        self.calc_fun_matrix(self.target, topo)

        # STEP 2: Select a possible time from the continuous density
        time = self.select_time(self.target, topo)

        # STEP 3: Select a possible nucleotide sequence
        sequence = self.select_sequence(self.target, topo, time)

        return self.target, topo, time, sequence

    def calc_mutation_prob(self, seq1, seq2, t):
        freqs = calc_nuc_sites(seq1, seq2)
        mutate_probs = []
        
        for change in freqs:
            start_nuc = change[0]
            end_nuc = change[1]
            prob = self.model.matrix[start_nuc][end_nuc](t)
            mutate_probs.append(prob)

        return prod(mutate_probs)

    def init_fun_matrix(self):
        fun_dict = {}
        generator = product("ATCG", repeat=3)
        
        for seq in generator:
            fun_dict[seq] = None

        return fun_dict

    def generate_topology(self, target):
        self.old_topo = [target.left, target.right, target.sibling]
        choices = unique_perms(self.old_topo)
        self.new_topo = choice(choices)
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
        upper = target.parent.time
        lower = max(topo[1].time, topo[2].time) # max of child times
        result = 0
        for prob in self.matrix:
            result += quad(self.matrix[prob], lower, upper)[0]
        return result

    def sum_numer(self):
        return lambda t: sum(prob(t) for prob in self.matrix.values()) 
    
    def calc_density(self, target, topo):
        # self.matrix = self.calc_fun_matrix(target, topo)
        self.nuc_sum = self.sum_numer()
        density = lambda t: self.nuc_sum(t)/self.prob_integral(target, topo)
        return density

    def select_time(self, target, topo):
        self.density = self.calc_density(target, topo)
        sample = self.inverse_cdf(self.density, target, topo, 100)
        return sample

    def inverse_cdf(self, density, target, topo, draws):
        lower = max(topo[1].time, topo[2].time)
        upper = target.parent.time
        grid = np.linspace(lower, upper, draws)
        uniform = np.random.uniform(0, 1, draws)
        pdf_values = np.zeros(draws)
        cum_cdf_values = np.zeros(draws)
        theta = np.zeros(draws)

        for i in xrange(draws):
            pdf_values[i] = density(grid[i])
            cum_cdf_values[i] = np.sum(pdf_values)/(upper - lower)
            theta[i] = grid[(cum_cdf_values[0:i] <= uniform[0:i]).sum()]

        return choice(theta) 

    def rejection_sample(self, density, target, topo):
        count = 0
        while True:
            count += 1
            lower = max(topo[1].time, topo[2].time)
            upper = target.parent.time
            uniform = np.random.uniform(lower, upper)
            prop = np.random.normal()
            dense = density(uniform)
            if prop < dense:
                print count
                print prop
                return uniform

    def select_sequence(self, target, topo, t):
        nucs = ('A', 'T', 'C', 'G')
        freqs = {'A': 0., 'T': 0., 'C': 0., 'G': 0.}
        new_seq = []
        parent = target.parent
        for seq in self.target.sequence:
            for nuc in nucs:
                p = self.model.matrix[seq][nuc](parent.time - t)
                c1 = self.model.matrix[nuc][seq](t - topo[1].time)
                c2 = self.model.matrix[nuc][seq](t - topo[2].time)
                freqs[nuc] += (p * c1 * c2)/4

            new_nuc = four_category_dist(1, freqs.values(), nucs)[0]
            new_seq.append(new_nuc)

        return new_seq

    def transition_prob(self, target, topo):
        integral = self.prob_integral(target, topo) 
        mutation_prob = self.calc_mutation_prob(target.parent.sequence, topo[0].sequence, target.parent.time - topo[0].time)
        return integral * mutation_prob

    def alpha(self, target):
        new_prob = self.transition_prob(target, self.new_topo)

        old_matrix = self.calc_fun_matrix(target, self.old_topo)
        old_prob = self.transition_prob(target, self.old_topo)

        prob = new_prob/old_prob

        return min(prob, 1.)

def main():
    freq = {'A': .292, 'C': .213, 'G': .237, 'T': .258}
    k = 22.2
    alpha = .0005
    tau = 1
    draws = 10
    filename = 'infile'
    mcmc = MCMC(freq, k, alpha, tau, draws, filename)
    return mcmc.metropolis()    

if __name__ == '__main__':
    main()

