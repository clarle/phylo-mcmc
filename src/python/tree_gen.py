from models import HKY85
from tree import Tree, Node
from utils import four_category_dist

def make_nodes(node, model, height, time):
    if height == 0:
        print ''.join(node.sequence)
        return
    else:
        time += (1./21.) * height
        # print time
        model.update(time)
        
        left_seq = mutate(node.sequence, model, time)
        right_seq = mutate(node.sequence, model, time)

        node.left = Node(time=time, sequence=left_seq)
        node.right = Node(time=time, sequence=right_seq)
        # print str(height) + ":" + str(node.left.sequence)
        # print str(height) + ":" + str(node.right.sequence)

        make_nodes(node.left, model, height - 1, time)
        make_nodes(node.right, model, height - 1, time)

def mutate(parent_seq, model, time):
    nucs = ('A','C','G','T')
    new_seq = []
    for parent_nuc in parent_seq:
        probs = []
        
        for child_nuc in nucs:
            probs.append(model.matrix[parent_nuc][child_nuc](time))
        # print probs
        new_nuc = four_category_dist(1, probs, nucs)[0]
        new_seq.append(new_nuc)

    return new_seq

alpha = 0.0005
k = 22.2
freq = {'A': 0.292, 'C': 0.213, 'G': 0.237, 'T': 0.258}

model = HKY85(freq, k, alpha)
root = Node(time=0, sequence=['A','T','G','G','G','C'])

make_nodes(root, model, 6, 0)
print "6 64"
