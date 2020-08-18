import numpy as np
import pickle

def open_dyn(file):
    n = 0
    L = []
    with open(file, "r") as file:
        for line in file:
            line = line.strip()
            [u, v, t] = list(map(int, line.split(" ")))
            L.append((t, u, v))
            n = max(n, u, v)    

    Tmax = L[0][0]
    L.reverse()

    return (L, n+1, Tmax)

def write_dyn(file, L):
    with open(file, 'w') as f:
        for e in reversed(L):
            (t, i, j) = e
            f.write(str(i)+" "+str(j)+" "+str(t)+"\n")

def LS_from_positions(positions):
    L = []
    n = len(positions)
    Tmax = len(positions[0]) - 1
    for t in range(0, Tmax+1):
        for i in range(n):
            for j in range(i):
                pos1 = np.array(positions[i][t])
                pos2 = np.array(positions[j][t])
                distsq = np.sum((pos1-pos2)**2)
                if distsq <= 1:
                    L.append((t, i, j))
    
    return L, n, Tmax

def write_pos(file, positions):
    with open(file, "wb") as f:
        pickle.dump(positions, f)

def open_pos(file):
    with open(file, "rb") as f:
        return pickle.load(f)