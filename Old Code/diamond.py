from numpy import *

def vgen(total_sites, total_neighbors, bb, tv):
    vec = [ [[0,0,0] for i in xrange(total_neighbors)] for j in xrange(total_sites)]

    BB = len(bb)
    for i in range(total_sites):
        cs = i + total_sites*13
        res_list = diamond_site(3, cs) #change this line
        j = 0
        for ts in res_list:
            s = cs%BB
            u = cs/BB
            a = u%3
            c = u/(9)
            b = (u-a-9*c)/3
            vtemp = bb[s] + a*tv[0,:] + b*tv[1,:] + c*tv[2,:]
            s = ts%BB
            u = ts/BB
            a = u%3
            c = u/(9)
            b = (u-a-9*c)/3
            vtemp = ( bb[s] + a*tv[0,:] + b*tv[1,:] + c*tv[2,:] ) - vtemp
            for zoidberg in range(3):
                vec[i][j][zoidberg] = vtemp[zoidberg]

            j+=1

    return vec

def vecdiamond_site():
    """vector displacements for the (6,4) lattice for each site on the unit lattice"""
    nsites = 2
    nres = 4

    bb = zeros((2,3))
    tv = zeros((3,3))

    bb[0,:] = array([0.,0.,0.]) #
    bb[1,:] = array([1.,0.,0.]) # 
    bb[2,:] = array([0.,0.5,0.]) # 
    bb[3,:] = array([0.5,0.,0.]) #
    bb[4,:] = array ([0.,0.,0.5])#
    bb[5,:] = array ([1.,0.,0.5])#
    
    tv[0,:] = array([2,0,0]) #u
    tv[1,:] = array([1,1,0])  #v
    tv[2,:] = array([1,0,1]) #w

    
    
    return vgen(nsites, nres, bb, tv)


def diamond_site(L, i):
    """calculates the nearest neighbors for each site on the lattice"""
    s = i % 2
    u = (i/2) % L
    v = (i/(2*L)) % L
    w = (i/(2*L*L)) % L

    if s == 0:
        V = [[ 1, u-1, v,   w+1  ],
            [  1, u,   v-1, w  ],
            [  1, u,   v,   w-1  ],
            [  1, u,   v,   w ]]
    elif s == 1:
        V = [[ 0, u,   v,   w  ],
            [  0, u+1, v,   w-1  ],
            [  0, u,   v+1, w  ],
            [  0, u,   v,   w+1 ]]
    else:
        # this should not occur
        return None
        
    nn = []
        
    # adjust for wrapping
    for coord in V:
        coord[1] %= L
        coord[2] %= L
        coord[3] %= L
        
        # add index to list
        nn.append(coord[0] + 2*(coord[1] + L*(coord[2] + coord[3]*L)))
        
    return nn




    return vgen(nsites, nres, bb, tv)

print vecdiamond_site()