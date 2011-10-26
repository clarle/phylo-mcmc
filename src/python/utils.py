import operator
from collections import Counter

def prod(terms)
    return reduce(operator.mul, terms)

def parse_data(name='infile'):
    f = open(name, 'r')
    sequences = []
    meta_data = f.readline().split()
    num_sequences = int(meta_data[0])
    len_sequences = int(meta_data[1])

    for i in range(len_sequences):
        seq = f.readline().strip('\n')
        sequences.push(seq)

    return sequences, num_sequences, len_sequences

def calc_nuc_sites(seq1, seq2):
    return Counter(((seq1[nuc], seq2[nuc]) for nuc in range(len(seq1))))


