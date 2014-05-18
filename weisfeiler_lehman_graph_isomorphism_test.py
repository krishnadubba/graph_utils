from igraph import Graph

def weisfeiler_lehman_graph_isomorphism_test(G1,G2):
    """ Performs the Weisfeiler-Lehman test of Isomorphism:
    Weisfeiler and Lehman: A reduction of a graph to a canonical form 
    and an algebra arising during this reduction, 
    Nauchno-Technicheskaya Informatsiya, Ser. 2, no. 9 (1968), 12-16 (in Russian).
    I used this paper to implement the algorithm:
    Nino Shervashidze et.al.: Weisfeiler-Lehman Graph Kernels, Journal of Machine Learning Research (2011)
    """
    name_codes = {}
    MAX_ITRES  = 2
    node_count = 0
    compressed_labels = {}
    for h in range(MAX_ITRES):
        multi_set_labels = {}
        for g in [G1, G2]:
            multi_set_labels[g['name']] = {}
            for node in g.vs:
                neighbours = g.neighbors(node)
                lables = g.vs[neighbours]['name']
                lables.sort()
                new_node_name = node['name'] + '_'.join(lables)
                multi_set_labels[g['name']][node.index] = new_node_name
                if new_node_name not in compressed_labels:
                    
                    compressed_labels[new_node_name] = 'c' + repr(node_count)
                    print new_node_name, compressed_labels[new_node_name]
                    node_count += 1
                
        for g in [G1, G2]:
            for node in g.vs:
                node['name'] = compressed_labels[multi_set_labels[g['name']][node.index]]
        
if __name__ == '__main__':
    G1 = Graph([(0,4), (1,4), (4,5), (3,5), (3,4), (2,5), (2,3)])
    G2 = Graph([(0,4), (1,3), (4,5), (3,5), (3,4), (2,5), (2,4)])
    
    G1.vs["name"] = ["1", "1", "2", "3", "4", "5"]
    G2.vs["name"] = ["1", "2", "2", "3", "4", "5"]
    
    G1["name"] = 'G1'
    G2["name"] = 'G2'
    
    weisfeiler_lehman_graph_isomorphism_test(G1, G2)
    
    print G1
    print G2