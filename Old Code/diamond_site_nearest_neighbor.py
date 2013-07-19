import numpy as np

NUM_SITES = 2
NAME = "diamond_site"

def id_to_vector(site_id, L):
    s = site_id % NUM_SITES
    u = (site_id / NUM_SITES) % L
    v = (site_id / (NUM_SITES * L)) % L
    w = (site_id / (NUM_SITES * L**2)) % L
    return np.array([s, u, v, w], dtype=np.int32)

def generate_vectors(s_displace, u_displace, site_function):
    """Generates the list of displacement vectors of neighbors for every site."""

    # get the list of neighbors
    null, L = site_function(0)

    disp_sites = []

    # for each site
    for site_id in xrange(NUM_SITES):

        # Pick the center of the lattice
        center_unit_id = site_id + NUM_SITES*int((L**3)/2)
        center_coord = id_to_vector(center_unit_id, L)

        # get the neighbors
        neighboring_sites, null = site_function(center_unit_id)
        neighbors = map(lambda ns : id_to_vector(ns, L), neighboring_sites)

        # allocate space for the neighbors
        neighbor_temp = []

        for i in xrange(len(neighbors)):

            # the current neighbor
            neighbor_coord = neighbors[i]

            # basis transformation (u,v,w)->(x,y,z)
            neighbor_disp = np.array(np.matrix(neighbor_coord[1:]) * np.matrix(u_displace))[0]
            origin_disp = np.array(np.matrix(center_coord[1:]) * np.matrix(u_displace))[0]

            # calculate the displacement vector
            displacement_temp = s_displace[neighbor_coord[0]] + neighbor_disp \
                                - s_displace[center_coord[0]] - origin_disp

            # add to temporary array
            neighbor_temp.append(displacement_temp)

        # add to list of site displacements
        disp_sites.append(neighbor_temp)

    return disp_sites

def lattice_sites(site_id, L=3):
    """Calculates the nearest neighbors for each bond on the lattice.  If there
    is no value for L specified, it uses the minimum lattice size necessary to
    calculate the displacement vectors.

    Lattice name: Diamond
    Lattice type: Bond
    Last Updated: July, 18 2013"""

    s, u, v, w = id_to_vector(site_id, L)

    assert s < NUM_SITES, "Error: s=%d >= NUM_SITES=%d" % (s, NUM_SITES)

    if s == 0:
        V = [[ 1, u-1, v, w+1  ],
            [  1, u, v-1 , w  ],
            [  1, u, v, w-1  ],
            [  1, u, v, w ]]
            
    elif s == 1:
        V = [[ 0, u, v, w+1  ],
            [  0, u, v+1, w  ],
            [  0, u+1, v, w-1  ], 
            [  0, u, v, w]]
   
    else:
        print "Error: Site %d does not exist in the unit cell!" % s
        exit() 
        
    neighbors = []
        
    # adjust for wrapping
    V = np.array(V)
    for coord in V:
        # wrap
        coord %= L
        
        # add index to list
        neighbors.append(coord[0] + NUM_SITES*(coord[1] + L*(coord[2] + coord[3]*L)))
        
    return np.array(neighbors), L

def lattice_displacement():
    """Vector displacements for the lattices; for each site on the unit lattice.
    
    Lattice name: Diamond
    Lattice type: Site
    Last Updated: July, 18 2013"""

    site_displacements = []   # Do not edit
    unit_displacements = []   # Do not edit

    # for each site
    site_displacements.append(np.array([0.,0.,0.]))
    site_displacements.append(np.array([0,0.,1]))
    site_displacements.append(np.array([0.,1,1]))
    site_displacements.append(np.array([1,0.,1]))
     
    
    
    # for each direction (u, v, w)
    unit_displacements.append(np.array([2,0,0]))
    unit_displacements.append(np.array([1,1,0]))
    unit_displacements.append(np.array([1,0,1]))

    return generate_vectors(site_displacements, unit_displacements, lattice_sites)

# ! Do not edit the following! #############

# debugging information
vectors = lattice_displacement()
for site in vectors:
    print site

print [[v for v in s] for s in vectors]

# the maximum length
maxlen = np.sqrt(np.max(np.max( \
            np.array([[np.array(v).dot(np.array(v)) 
                    for v in s] 
                for s in vectors]))))
print maxlen
