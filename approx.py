from all_gamma_edges import all_gamma_edges_time, all_gamma_edges
from itertools import product, chain, combinations
from math import log2
from bisect import bisect_left

def powerset(l):
    s = list(l)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

# M = list of gamma edges
def is_matching(M, n, Tmax, gamma):
    seen = set()
    for e in M:
        (t, u, v) = e
        for t_ in range(t, t+gamma):
            if (t_, u) in seen or (t_, v) in seen:
                return False
            seen.add((t_, u))
            seen.add((t_, v))
    return True

def thickness(A, Tmax):
    nbA = [0 for t in range(Tmax+1)]
    for a in A:
        (t, _, _) = a
        nbA[t] += 1
    return max(nbA)

def density_of_LS(G, positions, Tmax, gamma, dim, speed):
    (n, L) = G
    flatten = lambda l: [item for sublist in l for item in sublist]
    gamma_edges = flatten(all_gamma_edges_time(G, Tmax, gamma))
    
    X = [tuple((positions[u][t][i] + positions[v][t][i])/(2 + 2*(gamma-1)*speed) for i in range(dim)) for (t, u, v) in gamma_edges]
    return density(dim, gamma_edges, X, Tmax)

def density(n, A, X, Tmax):
    """
    A is the set of gamma-edges, and X = [c[a][n-p:] for a in A]
    We can shift each hypercube C of length 1 defined by (x1, x2, ..., xn) and (x1+1, x2+1, ..., xn+1)
    such that the hypercube still contains the same centers of gamma-edges and such that
    each face of coordinate x = x1, y = x2, ... contains a center of a gamma-edge contained in C 
    """
    def in_hypercube(pos, C):
        for i in range(n):
            if not (C[i][i] <= pos[i] <= C[i][i]+1):
                return False
        return True

    if n == 0:
        return thickness(A, Tmax)

    # for each C_list in X^n
    maxi = 0
    for C_list in product(*[X]*n):
        # C = [C_list[i][i] for i in range(n)]
        # C defines an hypercube of length 1
        A_ = [A[i] for i in range(len(A)) if in_hypercube(X[i], C_list)]
        """for i in range(len(A)):
            if in_hypercube(X[i], C_list):
                A_.append(A[i])"""
        maxi = max(maxi, thickness(A_, Tmax))

    return maxi

def get_MM(G, Tmax, gamma = 2):
    flatten = lambda l: [item for sublist in l for item in sublist]

    (n, L) = G
    
    # get all gamma-edges
    gammas = all_gamma_edges_time(G, Tmax, gamma)
    return base_case(flatten(gammas), n, Tmax, gamma)

# G = (nb of vertices, edges/link-stream) = (n, L)
def base_case(gammas_flattened, n, Tmax, gamma):
    flatten = lambda l: [item for sublist in l for item in sublist]
    
    if Tmax + 1 < gamma:
        return 0

    gammas = [[] for t in range(Tmax+1)]
    for (t, u, v) in gammas_flattened:
        gammas[t].append((t, u, v))

    powerset_gammas = [tuple(A for A in powerset(gammas[t]) if is_matching(A, n, Tmax, gamma)) for t in range(len(gammas))]

    # vertices are numbered 0 to n-1
    M = dict()
    for S_list in product(*powerset_gammas[:gamma - 1]):
        union = flatten(S_list)
        if is_matching(union, n, Tmax, gamma):
            M[0, S_list] = len(union)

    # while t + gamma  - 1 <= Tmax - gamma + 1
    for t in range(0, Tmax - 2*gamma + 3):
        for S_list in product(*powerset_gammas[t+1:t+gamma]):
            T = flatten(S_list)
            for S_ in powerset(gammas[t]):
                union = T + list(S_)
                if is_matching(union, n, Tmax, gamma):
                    # print((t+1, S_list))
                    N = len(S_list[-1]) + M[t, (S_, *(S_list[:-1]))]
                    M[t+1, S_list] = max(M[t+1, S_list], N) if (t+1, S_list) in M else N

    mm = 0

    t = max(Tmax - 2*gamma + 2, -1)

    for S_list in product(*[tuple(powerset(gammas[t_])) for t_ in range(t+1, t+gamma)]):
        if is_matching(flatten(S_list), n, Tmax, gamma):
            mm = max(mm, M[t+1, S_list])

    return mm

def approx(G, n, q, positions, Tmax, gamma, speed):
    c = dict()
    A = []
    gammas = all_gamma_edges_time(G, Tmax, gamma)
    for gammas_t in gammas:
        for (t, u, v) in gammas_t:
            # vector calculus
            c[t, u, v] = tuple((positions[u][t][i] + positions[v][t][i])/(2 + 2*(gamma-1)*speed) for i in range(n))
            A.append((t, u, v))
    m = len(A)

    def aux(p, A):
        if len(A) == 0:
            return 0
        f = q**(n-p+1) * 2**(n-p-1) * log2(m) / gamma
        A = sorted(A, key = lambda a: c[a][n-p])
        X = [c[a][n-p:] for a in A]
        d = density(p, A, X, Tmax)
        k = int(f/d)

        xs = [pos[0] for pos in X]
        unique_xs = sorted(set(xs))
        boundaries = [unique_xs[0]-0.51]

        for i in range(len(unique_xs) - 2):
            # next candidate for x_i
            candidate = (unique_xs[i] + unique_xs[i+1])/2
            delim1 = boundaries[0] if len(boundaries) == 1 else boundaries[-1]+0.5
            delim2 = candidate - 0.5
            i1 = bisect_left(xs, delim1)
            i2 = bisect_left(xs, delim2)
            P = A[i1:i2]
            
            dH = density(p-1, P, [c[a][n-p+1:] for a in P], Tmax)
            if dH >= f:
                boundaries.append(candidate)
        
        boundaries.append(xs[-1] + 2*max(0.5, f/d) + 0.001)


        def P(s):
            res = 0
            for i in range(len(boundaries) - 1):
                # next candidate for x_i
                delim1 = boundaries[0]+s if i == 0 else boundaries[i]+s+0.5
                delim2 = boundaries[-1]+s if i == len(boundaries)-2 else boundaries[i+1]+s-0.5
                i1 = bisect_left(xs, delim1)
                i2 = bisect_left(xs, delim2)
                P = A[i1:i2]
                if p != 1:
                    res += aux(p-1, P)
                else:
                    res += base_case(P, p-1, Tmax, gamma)
            return res
        
        maxi = 0
        for s in range(k):
            maxi = max(maxi, P(s))

        return maxi
    
    return aux(n, A)
