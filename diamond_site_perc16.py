# Jonathan Tran August 16, 2012
# Documentation by Shane Stahlheber August 10, 2012
# Percolation wrapping detection on the (10,3)-b lattice

# update: displacements between sites now written as an array which is called whenever
# needed further simplifying the functions. these displacements are generated in the imported 
# nearest neighbor function.

import time
import math
import random as rand
from numpy import *
from scipy import *

# import .py file containing all the nearest neighbor functions for the k4 lattice
import diamond_site_neighbor as site_definition

lattice = site_definition.lattice_sites
vecb = site_definition.lattice_displacement

# global constants
# check debug mode
DEBUGMODE = False
if DEBUGMODE:
    null, L = lattice(0)
    NUM_RUN = 1
else:
    L = 16                      # Number of units cells on a lattice axes.
    NUM_CPU = 1
    NUM_RUN = 1000/NUM_CPU + 1  # The number of lattices to test.

S = site_definition.NUM_SITES   # Do not edit! The number of sites per unit.
N = (L**3)*S                    # Do not edit! Total number of sites.
EMPTY = -N-1                    # Do not edit!

EPSILON = 0.00000001            # Do not edit!


ptr = zeros(N)                      # array of N pointers  
order = [ i for i in xrange(N) ]    # occupation order
data = []                           # stores all the calculated percolation thresholds

def initialize():
    """Initializes all of the arrays used."""
    global nn, du, dv, dw
    
    # initialize displacement arrays
    du = array([ 0. for i in xrange(N)])
    dv = array([ 0. for i in xrange(N)])
    dw = array([ 0. for i in xrange(N)])

def boundaries():
    """Fills the nearest-neighbors array with the neighbors of each site in the lattice."""

    global vec
    vec = vecb()


def permutation():
    """Determines the random order sites will be occupied (Shuffle)."""
    for i in xrange(N):
        j = int(i + (N - i) * random.rand())
        order[i], order[j] = order[j], order[i]

def findroot(i):
    """Recursive function used to traverse through parent sites until a root site is found."""
    
    # Check if we have found the root
    if ptr[i] < 0: 
        du[i] = dv[i] = dw[i] = 0
        return i, 0., 0., 0.
        
    # Continue searching for the root if not found yet
    root, dispu, dispv, dispw = findroot(ptr[i])
    
    # We have found the root; commence path compression
    ptr[i] = root
    
    # update the vector displacements for parent sites on way to root
    du[i] += float(dispu)
    dv[i] += float(dispv)
    dw[i] += float(dispw)
    
    # return the root index and the displacement from thus far
    return root, du[i], dv[i], dw[i]

def delta3(r2, s1, s2, u, v, w):
    """."""
    du[r2] = du[s1] - du[s2] - u
    dv[r2] = dv[s1] - dv[s2] - v
    dw[r2] = dw[s1] - dw[s2] - w


def displacement(j,s1,s2): 
    """Finds the displacement of a newly placed site to root with respect to the neighbor."""
    i = s1%S
    u = vec[i][j][0]
    v = vec[i][j][1]
    w = vec[i][j][2]
    du[s1] = du[s2] + u
    dv[s1] = dv[s2] + v
    dw[s1] = dw[s2] + w

def detection(s1,s2,maxlen=1.0):
    """Checks if percolation has occurred."""
    
    # get the difference in displacements
    delu = float(du[s1] - du[s2])
    delv = float(dv[s1] - dv[s2])
    delw = float(dw[s1] - dw[s2])
    
    # calculate the magnitude of displacement difference
    netdisplacement = sqrt(delu**2 + delv**2 + delw**2)
    
    # greater than one with an epsilon value because we are dealing with floating point displacements
    return (netdisplacement > EPSILON+site_definition.maxlen)

def percolate(iteration):
    """Run percolation testing."""  
    hasPercolated = False       # initially no percolation
    
    # reset the pointer array
    for i in xrange(N):
        ptr[i] = EMPTY
    
    # for all sites
    for i in xrange(N):
        
        # set the current site as its own root
        r1 = order[i]
        s1 = r1
        ptr[s1] = -1
        du[s1] = 0
        dv[s1] = 0
        dw[s1] = 0

        nn, null = lattice(s1, L)
        
        # for all neighbors
        for j in xrange(len(nn)):
        
            # define the second site as the first's neighbor
            s2 = nn[j]
            
            # check if the neighbor is already occupied
            if ptr[s2] != EMPTY:
            
                # find the root of the second site (and do path compression)
                r2, un1, un2, un3 = findroot(s2)   # un1,un2,un3 are unused variables
                
                # check if they have different roots
                if r2 != r1:
                    if ptr[r1] >= ptr[r2]: # damn, remember those ptr[r1] and ptr[r2] are negative
                        
                        # root2 is larger than root1; add root1 size to root2
                        ptr[r2] += ptr[r1]
                        
                        # set the roots of the initial sites to the second root
                        ptr[r1] = r2
                        ptr[s1] = r2
                        
                        # update the vector from one cluster root to another
                        s = s1%S
                        u = vec[s][j][0]
                        v = vec[s][j][1]
                        w = vec[s][j][2]

                        du[r1] = du[s2] - du[s1] + u
                        dv[r1] = dv[s2] - dv[s1] + v
                        dw[r1] = dw[s2] - dw[s1] + w


                        # add the displacement of site2 (currently pointing directly to root2) to the displacement of site1
                        displacement(j, s1, s2)
                        
                        # set the current root to root2
                        r1 = r2
                        
                    else:
                    
                        # root1 is larger than root2; add root2 size to root1
                        ptr[r1] += ptr[r2]
                        
                        # set the root of the neighbor's root to the current root
                        ptr[r2] = r1
                        
                        # update the vector from one cluster root to another
                        s = s1%S
                        u = vec[s][j][0]
                        v = vec[s][j][1]
                        w = vec[s][j][2]
                        du[r2] = du[s1] - du[s2] - u
                        dv[r2] = dv[s1] - dv[s2] - v
                        dw[r2] = dw[s1] - dw[s2] - w

                else:
                    # They have the same root; check for percolation
                    hasPercolated = detection(s1, s2)
                    
                    # if percolation has occurred, return
                    if hasPercolated:
                        data[iteration,1] = i
                        data[iteration,2] = order[i]
                        return i

def main():
    global data
    
    # initialize array of results
    data = zeros((NUM_RUN, 3))
    for i in xrange(NUM_RUN):
        data[i,0] = i                        

    initialize()
    boundaries()
    for iteration in xrange(NUM_RUN):
        print iteration
        permutation()
        percolate(iteration)
    
    TIMESTAMP = time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())
    fout = open("%s_data_%d_%s%s.txt" % (site_definition.NAME, L, TIMESTAMP, "DEBUG" if DEBUGMODE else ""), "w")
    fout.write("Percolation threshold data\n")
    fout.write('Total number of sites: %d\n\n' % (N))
    fout.write('Run\tNum\tKey Site\n')
    fout.write('---\t---\t--------\n')
    for j in xrange(NUM_RUN):
        fout.write('%d\t%d\t%d\n' % (data[j,0], data[j,1], data[j,2]))
        #fout.write(str(data[j,0]) + '\t' + str(data[j,1]) + '\n')
    
    if DEBUGMODE:
        fout.write("\nOccupied Sites\n")
        fout.write("--------------\n")
        for n in xrange(len(ptr)):
            if ptr[n] != EMPTY:
                v = site_definition.id_to_vector(n, L)
                if ptr[n] < 0:
                    fout.write("Site %d (%d, %d, %d, %d): Root of size %d\r\n" % (n, v[0], v[1], v[2], v[3], ptr[n]+1+N))
                else:
                    fout.write("Site %d (%d, %d, %d, %d): %d\r\n" % (n, v[0], v[1], v[2], v[3], ptr[n]))

    fout.close()

    print L

# run percolation test
main()
