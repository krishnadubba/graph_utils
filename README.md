graph_utils
===========

Some graph utilities that I use for my work.

graph_hash:
-----------
Generates hash for entire graph. See Figure 4 in 'kLog: A Language for Logical and Relational Learning with Kernels' for the algorithm. Requires igraph package as the graph is expected to be generated using igraph. 

Input: An igraph graph, node_name attribute and edge_name_attribute. Note that edge_name_attribute is optional i.e. for graphs without edge labels or to ignore edge labels, edge_name_attribute is None.

Output: graph hash value.
