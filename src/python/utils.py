import operator
import numpy as np
from collections import Counter
from itertools import permutations
from random import choice, uniform

def prod(terms):
    return reduce(operator.mul, terms)

def parse_data(name='infile'):
    f = open(name, 'r')
    sequences = []
    meta_data = f.readline().split()
    num_sequences = int(meta_data[0])
    len_sequences = int(meta_data[1])

    for i in range(len_sequences):
        seq = f.readline().strip('\n')
        sequences.append(seq)

    return sequences, num_sequences, len_sequences

def calc_nuc_sites(seq1, seq2):
    return Counter(((seq1[nuc], seq2[nuc]) for nuc in range(len(seq1))))

def four_category_dist(draws, probs, alphabet):
    a = np.random.multinomial(draws, probs)
    b = np.empty(draws, dtype='|S1')

    upper = np.cumsum(a)
    lower = upper - a

    for value in range(len(a)):
        b[lower[value]:upper[value]] = alphabet[value]

    np.random.shuffle(b)
    return b

def merge_seqs(seq1, seq2):
    parent_seq = []
    
    for i in range(len(seq1)):
        nuc = choice((seq1[i], seq2[i]))
        parent_seq.append(nuc)

    return parent_seq

def merge_time(time1, time2, tau):
    return uniform(min(time1, time2), tau)

def unique_perms(a):
    b=permutations(a)
    f=lambda per:all(x==y for x,y in zip(a,per)) or all(x!=y for x,y in zip(a,per))
    return filter(f,b)
