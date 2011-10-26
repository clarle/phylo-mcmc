import numpy as np
from math import exp
from collections import defaultdict

class HKY85:

    def __init__(self, freq, k, alpha):
        self.freq = freq
        self.k = k
        self.alpha = alpha
        self.matrix = self.calc_matrix()

    def lamda(self, base):
        if base is 'A' or base is 'G':
            return self.freq['A'] + self.freq['G']
        else:
            return self.freq['C'] + self.freq['T']

    def gamma(self, base):
        return 1. + ((self.k - 1.) * self.lamda(base))

    def expo(self, t, gamma=1):
        return exp(-self.alpha * gamma * t)

    def is_transition(self, start_base, end_base):
        nuc_dict = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}
        if nuc_dict[start_base] is end_base:
            return True
        else:
            return False

    def calc_matrix(self):
        matrix = defaultdict(defaultdict) # 2D dictionaries
        nucs = ['A', 'C', 'G', 'T']

        for start_nuc in nucs:
            for end_nuc in nucs:
                matrix[start_nuc][end_nuc] = self.calc_probs(start_nuc, end_nuc)

        return matrix

    def calc_probs(self, start, end):
        freq = self.freq
        lamda = self.lamda
        expo = self.expo
        gamma = self.gamma
        
        if start is end:
            return lambda t: freq[end] + freq[end] * ((1/lamda(end))-1) * expo(t) + ((lamda(end)-freq[end])/lamda(end)) * expo(t, gamma=gamma(end))    
        elif self.is_transition(start, end):
            return lambda t: freq[end] + freq[end] * ((1/lamda(end))-1) * expo(t) - (freq[end]/lamda(end)) * expo(t, gamma=gamma(end)) 
        else:
            return lambda t: freq[end] * (1 - expo(t))

    def update(self, t):
        self.t = t
        self.matrix = self.calc_matrix()
