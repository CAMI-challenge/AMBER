# This will return the EMDUnifrac distance only
def EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q):
    '''
    (Z, diffab) = EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
    This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
    and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
    Returns the weighted Unifrac distance Z and the flow F. The flow F is a dictionary with keys of the form (i,j) where
    F[(i,j)] == num means that in the calculation of the Unifrac distance, a total mass of num was moved from the node
    nodes_in_order[i] to the node nodes_in_order[j].
    '''
    num_nodes = len(nodes_in_order)
    Z = 0
    diffab = dict()
    partial_sums = P - Q
    for i in range(num_nodes - 1):
        val = partial_sums[i]
        partial_sums[Tint[i]] += val
        if val != 0:
            diffab[(i, Tint[i])] = lint[i, Tint[i]] * val  # Captures diffab
        Z += lint[i, Tint[i]] * abs(val)
    return (Z, diffab)


# This will return the EMDUnifrac distance only
def EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q):
    '''
    (Z, diffab) = EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q)
    This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
    and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
    Returns the unweighted Unifrac distance Z and the flow F. The flow F is a dictionary with keys of the form (i,j) where
    F[(i,j)] == 1 means that in the calculation of the Unifrac distance, a total mass of 1 was moved from the node
    nodes_in_order[i] to the node nodes_in_order[j].
    '''
    num_nodes = len(nodes_in_order)
    Z = 0
    diffab = dict()
    for i in range(num_nodes):
        if P[i] > 0:
            P[i] = 1
        if Q[i] > 0:
            Q[i] = 1
    partial_sums = P - Q
    for i in range(num_nodes - 1):
        val = partial_sums[i]
        partial_sums[Tint[i]] += val
        if val != 0:
            diffab[(i, Tint[i])] = lint[i, Tint[i]] * val  # Captures diffab
        Z += lint[i, Tint[i]] * abs(val)
    return Z, diffab
