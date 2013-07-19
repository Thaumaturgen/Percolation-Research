from numpy import *
import numpy as np
cc =2 
num = len (cc)
NUM_SITES = len(bb)
NUM_NEIGH = 6

def id_to_vector(site_id, L):
    s = site_id % NUM_SITES
    u = (site_id / NUM_SITES) % L
    v = (site_id / (NUM_SITES * L)) % L
    w = (site_id / (NUM_SITES * L**2)) % L
    return np.array([s, u, v, w])


def generate_vectors(s_displace, u_displace, sites):
    """Generates the list of displacement vectors of neighbors for every site."""

    L = 7 
    vec = zeros((NUM_SITES, NUM_NEIGH, 3))

    # for each site
    for site_id in xrange(NUM_SITES):

        # Pick the center of the lattice
        center_unit_id = site_id + NUM_SITES*int((L**3)/2)
        center_coord = id_to_vector(center_unit_id, L)

        # get the neighbors
        neighbors = array(map(lambda i : id_to_vector(i, L), sites(L, center_unit_id)))

        for i in xrange(NUM_NEIGH):

            # the current neighbor
            neighbor_coord = neighbors[i]

            # basis transformation (u,v,w)->(x,y,z)
            neighbor_disp = array(matrix(neighbor_coord[1:]) * matrix(u_displace))[0]
            origin_disp = array(matrix(center_coord[1:]) * matrix(u_displace))[0]

            # calculate the displacement vector
            vec[site_id, i, :] =  s_displace[neighbor_coord[0]] + neighbor_disp
            vec[site_id, i, :] -= s_displace[center_coord[0]] + origin_disp

    return vec

def diamond_site(L, site_id):
    """calculates the nearest neighbors for each bond on the lattice"""
    s = site_id % NUM_SITES
    u = (site_id / NUM_SITES) % L
    v = (site_id / (NUM_SITES * L)) % L
    w = (site_id / (NUM_SITES * L**2)) % L
    if s == 0:
        V = [[ 5, u, v, w-1  ], #
            [  4, u, v, w  ], #
            [  3, u, v, w  ], #
            [  2, u, v, w ]] #
    elif s == 1:
        V = [[ 4, u+1, v, w-1  ], #
            [  2, u+1, v-1, w  ], #
            [  5, u, v, w  ], #
            [  3, u, v, w ]] #
    elif s == 2: 
        V = [[ 1, u-1, v+1, w ], #
            [  0, u, v, w ]] #
    elif s == 3:
        v = [[0 ,u, v, w ], #
            [1, u, v, w ]]  #       
    elif s == 4:
        v = [[0, u, v, w], #
            [ 1, u-1,v, w+1]]  #       
    elif s == 5: 
        V = [[ 1, u, v, w ], #
            [  0, u, v, w+1 ]] # 
    
    else:
        # this should not occur
        return None
        
    neighbors = []
        
    # adjust for wrapping
    for coord in V:
        coord[1] %= L
        coord[2] %= L
        coord[3] %= L
        
        # add index to list
        neighbors.append(coord[0] + NUM_SITES*(coord[1] + L*(coord[2] + coord[3]*L)))
        
    return neighbors

def vecdiamond_site():
    """vector displacements for the (6,4) lattice for each site on the unit lattice"""

    bb = zeros((NUM_SITES,3))
    tv = zeros((3,3))
    cc = zeros(( num, 3))
    bb[0,:] = array([0.,0.,0.]) #0
    bb[1,:] = array([1.,0.,0.]) #1 
    cc[2,:] = array([0.,0.5,0.]) #2 
    cc[3,:] = array([0.5,0.,0.]) #3
    cc[4,:] = array ([0.,0.,0.5])
    cc[5,:] = array ([1.,0.,0.5])
    
    tv[0,:] = array([2,0,0]) #u
    tv[1,:] = array([1,1,0])  #v
    tv[2,:] = array([1,0,1]) #w

    return generate_vectors(bb, tv, diamond_site)

# debugging information
vectors = vecdiamond_site()
for site in vectors:
    print site

# the maximum length
maxlen = sqrt(np.max(array([[array(v).dot(array(v)) for v in s] for s in vectors])))
print maxlen
