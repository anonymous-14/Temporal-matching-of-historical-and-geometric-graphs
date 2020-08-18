def all_gamma_edges(G, Tmax, gamma):
    (n, L) = G

    res = [[] for i in range(Tmax+1)]

    # if j<i, last_seen[i][j] stores the interval
    # where the edge i-j was seen most recently
    last_seen = [[(Tmax+2, Tmax+2) for j in range(i)] for i in range(n)]
    
    for e in reversed(L):
        (t, u, v) = e
        if u < v:
            v, u = u, v
        # now v <= u

        if last_seen[u][v][0] == t+1:
            # extend the interval
            last_seen[u][v] = (t, last_seen[u][v][1])
            if last_seen[u][v][1] - last_seen[u][v][0] + 1 >= gamma:
                res[t].append((u, v))
        else:
            last_seen[u][v] = (t, t)
    
    return res

def all_gamma_edges_time(G, Tmax, gamma):
    (n, L) = G

    res = [[] for i in range(Tmax+1)]

    # if j<i, last_seen[i][j] stores the interval
    # where the edge i-j was seen most recently
    last_seen = [[(Tmax+2, Tmax+2) for j in range(i)] for i in range(n)]
    
    for e in reversed(L):
        (t, u, v) = e
        if u < v:
            v, u = u, v
        # now v <= u

        if last_seen[u][v][0] == t+1:
            # extend the interval
            last_seen[u][v] = (t, last_seen[u][v][1])
            if last_seen[u][v][1] - last_seen[u][v][0] + 1 >= gamma:
                res[t].append((t, u, v))
        else:
            last_seen[u][v] = (t, t)
    
    return res