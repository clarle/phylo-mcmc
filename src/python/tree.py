from numpy import random
from collections import Counter, deque
from utils import parse_data, merge_seqs, merge_time

class Tree:
    def __init__(self, **kwargs):
        tau = kwargs.get("tau", None)
        filename = kwargs.get("filename", "infile")
        seqs, num_seqs, len_seqs = parse_data(filename)
        leaf_nodes = self.generate_leaves(seqs, tau)
        self.root = self.generate_tree(leaf_nodes, tau)
        self.curr = self.root

    def generate_leaves(self, seqs, tau):
        leaf_nodes = []
        prev = None
        
        for seq in seqs:
            if prev is None:
                node = Node(time=0, sequence=list(seq))
                leaf_nodes.append(node) 
                prev = node
            else:
                node = Node(time=0, sequence=list(seq), sibling=prev)
                leaf_nodes[-1].sibling = node
                leaf_nodes.append(node)
                prev = None
        
        return deque(leaf_nodes)

    def generate_tree(self, leaves, tau):
        prev = None

        while True:
            left = leaves.popleft()
            right = leaves.popleft()
            p_seq = merge_seqs(left.sequence, right.sequence)
            p_time = merge_time(left.time, right.time, tau)
            
            if prev is None or len(leaves) == 0:
                node = Node(time=p_time, sequence=p_seq, left=left, right=right)
                left.parent = node
                right.parent = node
                leaves.append(node)
                prev = node
                if len(leaves) == 1:
                    break
            else:
                node = Node(time=p_time, sequence=p_seq, left=left, right=right, sibling=prev)
                left.parent = node
                right.parent = node
                leaves[-1].sibling = node
                leaves.append(node)
                prev = node

        prev.time = tau # root node must have time zero
        return prev # last node is root
        

    def nodes(self, filters=None, sorter=None):
        nodes = list(self.root.pre_order_iter(filter_fn=filters))
        if sorter:
            nodes.sort(sorter)
        return nodes
        
    def get_internal(self):
        return self.nodes(filters=lambda x: not x.is_leaf() and not x.is_root())

    def get_leaves(self):
        return self.nodes(filters=lambda x: x.is_leaf())

    def find(self, filters):
        return self.nodes(filters=filters)

    def update(self, old, new):
        stack = [self.root]
        node = None
        while stack:
            node = stack.pop(0)
            if node is not None and node is not old:
                child_nodes = [node.left, node.right]
                child_nodes.extend(stack)
                stack = child_nodes
            elif node is old:
                node = new
                break

        while node.parent != None:
            if node.parent.left == old:
                node.parent.left = node
                node = node.parent
            else:
                node.parent.right = node
                node = node.parent
        
        self.root = node

    def length(self):
        total = 0
        for edge in self.postorder_edge():
            if edge.length is not None:
                total += edge.length
        return total

    def local_rearrange(self, rand_fn=random.randint):
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
        self.parent = kwargs.get("parent", None)
        self.sibling = kwargs.get("sibling", None)
        self.left = kwargs.get("left", None)
        self.right = kwargs.get("right", None)

    def is_root(self):
        return not bool(self.parent)
    
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

    def pre_order_iter(self, filter_fn=None):
        stack = [self]
        while stack:
            node = stack.pop(0)
            if node is not None and (filter_fn is None or filter_fn(node)):
                yield node
            if node is not None:
                child_nodes = [node.left, node.right]
                child_nodes.extend(stack)
                stack = child_nodes

    def pre_order_traverse(self, node, filter_by=None):
        yield node
        self.pre_order_traverse(node.left, filter_by=filter_by)
        self.pre_order_traverse(node.right, filter_by=filter_by)

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


