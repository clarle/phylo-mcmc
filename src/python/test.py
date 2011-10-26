from models import HKY85
from math import exp
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

if __name__ == '__main__':
    unittest.main()
