from all_gamma_edges import all_gamma_edges
import numpy as np
from random import random


def calculate_vel(n, positions, Tmax):
    vel = 0
    for position_i in positions:
        for t in range(Tmax):
            pos1 = np.array(position_i[t])
            pos2 = np.array(position_i[t+1])
            dist = np.sum((pos1-pos2)**2)**(0.5)
            vel = max(vel, dist)
    return vel

# radius is the distance required for the points to be connected by an edge
def random_link_stream_with_velocity(n, Tmax, radius, vel, dim, boxlength = 1):
    L = []

    # copied from http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
    def new_pos(old_pos):
        u = np.random.normal(0, 1, dim)  # an array of d normally distributed random variables
        norm = np.sum(u**2)**(0.5)
        r = random()**(1.0/dim)
        # x is chosen uniformly in the d-ball of radius vel
        x = vel*r*u/norm
        return tuple(np.array(old_pos) + x)

    positions = [[tuple(boxlength*random() for _ in range(dim))] for _ in range(n)]

    for t in range(1, Tmax+1):
        for i in range(n):
            positions[i].append(new_pos(positions[i][-1]))

    for t in range(0, Tmax+1):
        for i in range(n):
            for j in range(i):
                pos1 = np.array(positions[i][t])
                pos2 = np.array(positions[j][t])
                distsq = np.sum((pos1-pos2)**2)
                if distsq <= radius**2:
                    L.append((t, i, j))
    
    return L, positions

# r is the distance required for the points to be connected by an edge
def random_link_stream(n, Tmax, radius, dim, boxlength = 1):
    L = []

    positions = [[tuple(boxlength*random() for _ in range(dim)) for _ in range(Tmax+1)] for _ in range(n)]

    for t in range(0, Tmax+1):
        for i in range(n):
            for j in range(i):
                pos1 = np.array(positions[i][t])
                pos2 = np.array(positions[j][t])
                distsq = np.sum((pos1-pos2)**2)
                if distsq <= r**2:
                    L.append((t, i, j))
    
    return L, positions, calculate_vel(n, positions, Tmax)