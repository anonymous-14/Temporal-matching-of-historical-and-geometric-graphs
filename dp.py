from all_gamma_edges import all_gamma_edges
from functools import reduce
from operator import or_ 
# G = (nb of vertices, edges/link-stream) = (n, L)

def is_in_set(x, A):
    return ((A >> x) & 1) == 1

def add_to_set(x, A):
    return A | (1<<x)

# time stamps go from 0 to Tmax included
def get_MM(G, Tmax, gamma = 2):
    (n, L) = G

    # get all gamma-edges
    gammas = all_gamma_edges(G, Tmax, gamma)

    # vertices are numbered 0 to n-1
    MM_t = []

    origin = 0

    calls = 0

    def MM(t, A_list):
        nonlocal calls
        if t == -1:
            return 0

        # if not already calculated
        if A_list not in MM_t[t-origin]:
            # variable for the max
            mm = MM(t-1, (*A_list[1:], 0))
            calls += 1

            for e in gammas[t]:
                # i-j is a valid gamma edge
                (i, j) = e
                # check if i and j are not in union = A_1 ∪ A_2 ∪ ...
                # where A_list = (A_1, A_2, ...)
                # we calculate union by folding bitwise or
                # on the elements of A_list
                union = reduce(or_, A_list, 0)
                if not (is_in_set(i, union) or is_in_set(j, union)):
                    # if i and j are not in union, add it to A_list[-1], we select this edge
                    new_A_list = (*A_list[:-1], add_to_set(i, add_to_set(j, A_list[-1])))
                    mm = max(mm, 1 + MM(t, new_A_list))
                    calls += 1
            
            MM_t[t-origin][A_list] = mm

        return MM_t[t-origin][A_list]

    for t in range(Tmax):
        # print(t)
        MM_t.append(dict())
        MM(t, tuple(0 for _ in range(gamma)))
        if len(MM_t) >= gamma+1:
            MM_t.pop(0)
            origin += 1

    mm = MM(Tmax - 1, tuple(0 for _ in range(gamma)))
    # print("Number of calls:", calls)
    return mm
