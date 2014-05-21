import numpy as np
    
def compute_similarity_score_from_node_label_hists(G1, G2):
    score_matrix   = []
    sim = np.zeros((len(G2.vs()),len(G1.vs())),float)
    common_spatial_rels = {}
    for G2_node in G2.vs():
        current_node_score = []
        row = G2_node.index
        for G1_node in G1.vs():
            col = G1_node.index
            if G1_node['type'] != G2_node['type']:
                sim[row,col] = 0.0
            else:  
                sim[row,col] = float(np.sum(np.minimum(G1_node['hist'],G2_node['hist'])))/np.sum(np.maximum(G1_node['hist'],G2_node['hist']))            
            current_node_score.append(sim[row,col])
        score_matrix.append(current_node_score)
    return sim, score_matrix
   
def compute_node_label_hist(G, node_labels):
    num_node_labels = len(node_labels)
    for node in G.vs():
        predecessor_hist = np.zeros(num_node_labels)
        successor_hist   = np.zeros(num_node_labels)        
        predecessors     = G.predecessors(node)
        predecessors_labels = G.vs()[predecessors]['name']
        for label in predecessors_labels:
            predecessor_hist[node_labels.index(label)] += 1
        successors = G.successors(node)
        successors_labels = G.vs()[successors]['name']
        for label in successors_labels:
            successor_hist[node_labels.index(label)] += 1
        G.vs[node.index]['pred_hist'] = predecessor_hist
        G.vs[node.index]['succ_hist'] = successor_hist
    return G        

def compute_graph_similarity_score(G1, G2):
    """
    Deriving phylogenetic trees from the similarity analysis of metabolic pathways
    Maureen Heymans, Ambuj K. Singh
    """
    from hungarian_assignment import Munkres, make_cost_matrix
        
    # Initializing the similarity matrix. Compute the score by matching the histograms
    sim, score_matrix = compute_similarity_score_from_node_label_hists(G1,G2)
    # normalize similarity matrix    
    sim = sim/np.linalg.norm(sim)
    
    simm = sim
    sim_old = sim        
    norm = 10000
    iters = 0
    MAX_ITERS = 10
    while (norm > 0.0001 and iters < MAX_ITERS):
        (total_rows, total_cols) = np.shape(sim_old)
        for vi in range(total_rows):
            for tj in range(total_cols):
                simi = 0.0
                for vk in range(total_rows):
                    if vi == vk:
                        continue
                    for tl in range(total_cols):
                        if tj == tl:
                            continue
                        for node_label in node_labels:
                            node_label_index = node_labels.index(node_label)
                            vi_label_indeg  = G2.vs[vi]['pred_hist'][node_label_index]
                            tj_label_indeg  = G1.vs[tj]['pred_hist'][node_label_index]
                            vi_label_outdeg = G2.vs[vi]['succ_hist'][node_label_index]
                            tj_label_outdeg = G1.vs[tj]['succ_hist'][node_label_index]
                            v_incoming_label = v_outgoing_label = t_incoming_label = t_outgoing_label = None
                            if G2.are_connected(vk,vi):
                                v_incoming_label = G2.vs[vk]['name']
                            if G2.are_connected(vi,vk):
                                v_outgoing_label = G2.vs[vk]['name']
                            if G1.are_connected(tl,tj):
                                t_incoming_label = G1.vs[tl]['name']
                            if G1.are_connected(tj,tl):
                                t_outgoing_label = G1.vs[tl]['name']
                                
                            if v_incoming_label == node_label and t_incoming_label == node_label:
                                # A1(i,j)                                    
                                if vi_label_indeg != 0 and tj_label_indeg != 0:
                                    simi += sim_old[vk][tl]/(vi_label_indeg*tj_label_indeg)
                                elif vi_label_indeg == 0 and tj_label_indeg == 0:
                                    simi += sim_old[vk][tl]/(total_rows*total_cols)
                                else:
                                    simi += 0.0
                                
                                # A2(i,j)
                                if vi_label_outdeg != 0 and tj_label_outdeg != 0:
                                    simi += sim_old[vk][tl]/(vi_label_outdeg*tj_label_outdeg)                                    
                                elif vi_label_outdeg == 0 and tj_label_outdeg == 0:
                                    simi += sim_old[vk][tl]/(total_rows*total_cols)
                                else:
                                    simi += 0.0
                            
                            if v_incoming_label != node_label and t_incoming_label != node_label:
                                # A3(i,j)                                    
                                if vi_label_indeg != total_rows and tj_label_indeg != total_cols:
                                    simi += sim_old[vk][tl]/((total_rows - vi_label_indeg) * (total_cols - tj_label_indeg))
                                elif vi_label_indeg == 0 and tj_label_indeg == 0:
                                    simi += sim_old[vk][tl]/(total_rows*total_cols)
                                else:
                                    simi += 0.0
                                
                                # A4(i,j)
                                if vi_label_outdeg != total_rows and tj_label_outdeg != total_cols:
                                    simi += sim_old[vk][tl]/((total_rows - vi_label_outdeg) * (total_cols - tj_label_outdeg))                                    
                                elif vi_label_outdeg == 0 and tj_label_outdeg == 0:
                                    simi += sim_old[vk][tl]/(total_rows*total_cols)
                                else:
                                    simi += 0.0
                                    
                            if v_incoming_label == node_label and t_incoming_label != node_label:
                                # D1(i,j)
                                if vi_label_indeg != 0 and tj_label_indeg != total_cols:
                                    simi += sim_old[vk][tl]/(vi_label_indeg * (total_cols - tj_label_indeg))
                                elif vi_label_indeg == 0 and (total_cols - tj_label_indeg) == 0:  
                                    simi += sim_old[vk][tl]/(total_rows*total_cols)
                                else:    
                                    simi += 0.0
                                    
                            if v_incoming_label != node_label and t_incoming_label == node_label:
                                # D2(i,j)
                                if vi_label_indeg != total_rows and tj_label_indeg != 0:
                                    simi += sim_old[vk][tl]/((total_rows - vi_label_indeg) * tj_label_indeg)
                                elif (total_rows - vi_label_indeg) == 0 and tj_label_indeg == 0:
                                    simi += sim_old[vk][tl]/(total_rows*total_cols)
                                else:    
                                    simi += 0.0
                                    
                            if v_outgoing_label == node_label and t_outgoing_label != node_label:
                                # D3(i,j)
                                if vi_label_outdeg != 0 and tj_label_outdeg != total_cols:
                                    simi += sim_old[vk][tl]/(vi_label_outdeg * (total_cols - tj_label_outdeg))
                                elif vi_label_outdeg == 0 and (total_cols - tj_label_outdeg) == 0:
                                    simi += sim_old[vk][tl]/(total_rows*total_cols)
                                else:    
                                    simi += 0.0
                                    
                            if v_outgoing_label != node_label and t_outgoing_label == node_label:
                                # D4(i,j)
                                if vi_label_outdeg != total_rows and tj_label_outdeg != 0:
                                    simi += sim_old[vk][tl]/((total_rows - vi_label_outdeg) * tj_label_outdeg)
                                elif (total_rows - vi_label_outdeg) == 0 and tj_label_outdeg == 0:  
                                    simi += sim_old[vk][tl]/(total_rows*total_cols)
                                else:    
                                    simi += 0.0
    
        norm = np.linalg.norm(simm - sim_old)
        sim_old = simm
        iters += 1
        
    hungarian_alg = Munkres()    
    score_matrix = simm
    max_score = np.max(score_matrix) * 2
    # Apply the Hungarian algorithm to the score matrix
    # First compute the cost matrix from score matrix
    # Cost is opposite of score, so we use cost = (max_score * 2 - score) to compute cost.
    # Scores are sometimes very less, so using sys.maxint instead of max_score sometimes makes all scores equal
    cost_matrix = make_cost_matrix(score_matrix,
                                   lambda cost: max_score - cost)
    indexes = hungarian_alg.compute(cost_matrix)
    return score_matrix, indexes


if __name__  == "__main__":
    from igraph import Graph
    
    node_labels = ["before", "after", "s", "d", "c"]

    g1 = Graph([(0,1), (0,2), (2,3), (3,4), (4,2), (2,5), (5,0), (6,3), (5,6)],directed=True)
    g1.vs["name"] = ["before", "after", "s", "d", "c", "s", "d"]
    g1.vs["type"] = ["temporal", "temporal", "spatial", "spatial", "spatial", "spatial", "spatial"]
    compute_node_label_hist(g1,node_labels)
    g1.vs["hist"] = [[1,2,4],[1,2,4],[1,2,4],[1,2,4],[1,2,4],[1,2,4],[1,2,4]]
    g2 = Graph([(0,1), (0,2), (2,3), (3,4), (4,2), (2,5), (5,0), (6,3), (5,6)],directed=True)
    g2.vs["name"] = ["before", "after", "s", "d", "s", "s", "d"]
    g2.vs["type"] = ["temporal", "temporal", "spatial", "spatial", "spatial", "spatial", "spatial"]
    g2.vs["hist"] = [[1,2,4],[1,2,4],[1,2,4],[1,2,4],[1,2,4],[1,2,4],[1,2,4]]
    
    g1 = compute_node_label_hist(g1, node_labels)
    g2 = compute_node_label_hist(g2, node_labels)    
    print compute_similarity_score_from_node_label_hists(g1, g2)
    print compute_graph_similarity_score(g1, g2)