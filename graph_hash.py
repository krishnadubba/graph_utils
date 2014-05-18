from igraph import Graph

def graph_hash(G, node_name_attribute='name', edge_name_attribute=None):
    """
    See Figure 4 in 'kLog: A Language for Logical and Relational Learning with Kernels'
    for the algorithm.
    
    Takes an igraph graph, node_name attribute and edge_name_attribute. Note that
    edge_name_attribute is optional i.e. for graphs without edge labels or to ignore edge labels, 
    edge_name_attribute is None.
    """
    # First get a hash for each node
    for node in G.vs:
        paths = G.get_shortest_paths(node)
        node_hashes = []
        for path in paths:
            if len(path) != 0:
                node_name = G.vs[path[-1]][node_name_attribute]
                if node_name == None:
                    node_name = repr(None)
                node_hashes.append((len(path), node_name))
        node_hashes.sort()
        node_hashes_string = ':'.join([repr(i) for i in node_hashes])
        node['hash_name'] = hash(node_hashes_string)
        
    # Use node hashes and generate a hash for each edge    
    edge_hashes = []
    if edge_name_attribute:
        edge_hashes = [(G.vs[edge.source]['hash_name'], G.vs[edge.target]['hash_name'], \
                                   edge[edge_name_attribute]) for edge in G.es]
    else:
        edge_hashes = [(G.vs[edge.source]['hash_name'], G.vs[edge.target]['hash_name']) \
                       for edge in G.es]
        
    # Combine these hashes and get a hash for the whole graph   
    edge_hashes.sort()
    edge_hashes_string = ':'.join([repr(i) for i in edge_hashes])
    return (hash(edge_hashes_string), G)

if __name__ == '__main__':
    g = Graph([(0,1), (0,2), (2,3), (3,4), (4,2), (2,5), (5,0), (6,3), (5,6)])
    # Give some names to vertices
    g.vs["name"] = ["Alice", "Bob", "Claire", "Dennis", "Esther", "Frank", "George"]
    # Give some names to edges. Note edge names are optional
    g.es["name"] = ["A", "B", "C", "D", "E", "F", "G", "H", "K"]
    ghash = graph_hash(g, "name", "name")
    print ghash
    