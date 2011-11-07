from models import HKY85
from mcmc import MCMC
from math import exp
from random import shuffle
import unittest

class ModelTest(unittest.TestCase):
    
    def setUp(self):
        freq = {'A': .292, 'C': .213, 'G': .237, 'T': .258}
        k = 22.2
        alpha = .0005
        self.model = HKY85(freq, k, alpha)
        
    def test_lamda_a(self):
        self.assertEqual(self.model.lamda('A'), .292 + .237)
        
    def test_lamda_c(self):
        self.assertEqual(self.model.lamda('C'), .213 + .258)
    
    def test_gamma_a(self):
        self.assertEqual(self.model.gamma('A'), 1 + (.292+.237)*(22.2-1))

    def test_expo(self):
        self.assertEqual(self.model.expo(100), exp(-0.0005 * 100))

    def test_expo_gamma(self):
        gamma = self.model.gamma('A')
        self.assertEqual(self.model.expo(100, gamma=gamma), exp(-0.0005 * 100 * gamma))

    def test_calc_probs_AA(self): # same nucleotide
        prob_fun = self.model.calc_probs('A', 'A')
        self.assertAlmostEqual(prob_fun(100), 0.7825546) 

    def test_calc_probs_AG(self): # transitional event
        prob_fun = self.model.calc_probs('A', 'G')
        self.assertAlmostEqual(prob_fun(100), 0.19447445) 
    
    def test_calc_probs_AC(self): # transversional event
        prob_fun = self.model.calc_probs('A', 'C')
        self.assertAlmostEqual(prob_fun(100), 0.0103881)

    def test_matrix_AA(self): # test full matrix
        prob = self.model.matrix['A']['A'](100)
        self.assertAlmostEqual(prob, 0.7825546)

    def test_matrix_AG(self):
        prob = self.model.matrix['A']['G'](100)
        self.assertAlmostEqual(prob, 0.19447445)

    def test_matrix_AC(self):
        prob = self.model.matrix['A']['C'](100)
        self.assertAlmostEqual(prob, 0.0103881)

class MCMCTest(unittest.TestCase):
    def setUp(self):
        freq = {'A': .292, 'C': .213, 'G': .237, 'T': .258}
        k = 22.2
        alpha = .0005
        tau = 1
        draws = 100
        filename = 'infile'
        self.mcmc = MCMC(freq, k, alpha, tau, draws, filename)

    def test_constructor(self):
        self.assertEqual(self.mcmc.tau, 1)

    def test_fun_matrix(self):
        self.assertEqual(len(self.mcmc.matrix), 64)

    def test_topology(self):
        nodes = self.mcmc.tree.get_internal()
        self.mcmc.target = nodes[0]
        self.mcmc.generate_topology(self.mcmc.target)
        # print [node.sequence for node in self.mcmc.old_topo]
        # print [node.sequence for node in self.mcmc.new_topo]

    def test_calc_matrix(self):
        nodes = self.mcmc.tree.nodes()
        # print len(nodes)
        # for node in nodes:
        #    print node.sequence
        self.mcmc.target = nodes[1]
        # print self.mcmc.target.sequence
        # print self.mcmc.target.parent.sequence
        # print self.mcmc.target.sibling.sequence
        # print self.mcmc.target.left.sequence
        # print self.mcmc.target.right.sequence
        self.mcmc.calc_fun_matrix(self.mcmc.target, self.mcmc.new_topo)
        self.assertEqual(len(self.mcmc.matrix), 64)

    def test_select_time(self):
        nodes = self.mcmc.tree.get_internal()
        self.mcmc.target = nodes[0]
        topo = self.mcmc.generate_topology(self.mcmc.target)
        self.mcmc.calc_fun_matrix(self.mcmc.target, self.mcmc.new_topo)

    def test_select_seq(self):
        nodes = self.mcmc.tree.get_internal()
        self.mcmc.target = nodes[0]
        target = self.mcmc.target
        topo = self.mcmc.generate_topology(self.mcmc.target)
        self.mcmc.calc_fun_matrix(self.mcmc.target, self.mcmc.new_topo)
        time = self.mcmc.select_time(target, topo)
        result = self.mcmc.select_sequence(target, topo, time)
        # print target.sequence
        # print result
        # print target.time
        # print time

    def test_accept(self):
        nodes = self.mcmc.tree.get_internal()
        self.mcmc.target = nodes[0]
        target = self.mcmc.target
        topo = self.mcmc.generate_topology(self.mcmc.target)
        self.mcmc.calc_fun_matrix(self.mcmc.target, self.mcmc.new_topo)
        # print [node.sequence for node in self.mcmc.old_topo]
        # print [node.sequence for node in self.mcmc.new_topo]
        # print self.mcmc.alpha(target)


if __name__ == '__main__':
    unittest.main()
