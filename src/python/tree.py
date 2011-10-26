from numpy import random
from collections import Counter

class Tree:
    def __init__(self, **kwargs):
        self.root = kwargs.get("root", None)

    def nodes(self, filters=None, sorter=None):
        nodes = [node for node in self.root.pre_order_traverse(filter_by=filters)]
        if sorter:
            nodes.sort(sorter)
        return nodes
        
    def get_internal(self):
        return self.nodes(filters=lambda x: not x.is_leaf())

    def get_leaves(self):
        return self.nodes(filters=lambda x: x.is_leaf())

    def find(self, filters):
        return self.nodes(filters=filters)

    def length(self):
        total = 0
        for edge in self.postorder_edge():
            if edge.length is not None:
                total += edge.length
        return total

    def local_rearrange(self, rand_fn):
        rearrange = rand_fn(1, 3) # random integer between 1 and 3
        if rearrange is 1:
            if self.root is self.root.parent.left:
                self.root.right = self.root.parent.right
            else:
                self.root.right = self.root.parent.left
        elif rearrange is 2:
            if self.root is self.root.parent.left:
                self.root.left = self.root.parent.right
            else:
                self.root.left = self.root.parent.left
        else:
            pass # keep the same arrangement

    def write_newick(self):
        # TO DO LATER
        pass

class Node:
    def __init__(self, **kwargs):
        self.time = kwargs.get("time", None)
        self.sequence = kwargs.get("sequence", None)
        self.edge = kwargs.get("edge", None)
        self.parent = kwargs.get("parent", None)
        self.left = kwargs.get("left", None)
        self.right = kwargs.get("right", None)

    def is_leaf(self):
        return not bool(self.left or self.right)

    def is_internal(self):
        return bool(self.left or self.right)

    def level(self):
        if self.parent:
            return self.parent.level() + 1
        else:
            return 0

    def height(node):
        if node == null:
            return -1
        else:
            return max(height(node.left), height(node.right)) + 1

    def pre_order_traverse(node, filter_by=None):
        yield node
        pre_order_traverse(node.left, filter_by=filter_by)
        pre_order_traverse(node.right, filter_by=filter_by)

    def post_order_traverse(self, filter_by=None):
        post_order_traverse(node.left, filter_by=filter_by)
        post_order_traverse(node.right, filter_by=filter_by)
        yield node

    def write_newick(self):
        # TO DO SOON
        pass

class Edge:
    def __init__(self, **kwargs):
        self.tail = kwargs.get("tail", None)
        self.head = kwargs.get("head", None)
        self.rootedge = kwargs.get("rootedge", False)
        self.length = kwargs.get("length", None)

    def collapse(self):
        pass

    def invert(self):
        pass

class Sequence:

    def __init__(self, **kwargs):
        self.nucs = kwargs.get("nucs", [])
        self.length = len(self.nucs)
        self.frequency = calculate_freq(self.nucs)

    def __len__(self):
        return self.length

    def calculate_freq(nucs):
        c = Counter(self.nucs)
        return list(c.elements())

    def compare_sequence(self, seq):
        mutated_index = []

        for i in self.length:
            if self.nucs[i] is seq.nucs[i]:
                mutated_index.push(i)

        return mutated_index


