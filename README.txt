The provided code can be compiled with g++.

Input:

   graph.txt:
      m n - m: the node size of the tree; n: the number of relationships
      x_i y_i - x_i is the parent of y_i

   freq1.txt:
      m - the weighted node size of the tree T_1
      x_i y_i - x_i is the node ID with y_i is the weight of node x_i

   freq2.txt:
      m - the weighted node size of the tree T_2
      x_i y_i - x_i is the node ID with y_i is the weight of node x_i

Output:

   number of nodes in similarity set S_1/difference set S_2
   the representative nodes for similarity set S_1
   the representative nodes for difference set S_2