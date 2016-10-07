import re
from sys import stdin
from os import listdir
from os.path import isfile, join
import numpy as np

def load_data(f):
    """Reads and parses an open file for information on the lattice tested."""

    # load the input
    filedata = f.read()

    # look at the header for the number of sites
    count_ = re.search(r"[a-zA-Z\s]*: (?P<num_sites>[0-9]+)", filedata)

    num_sites = int(count_.groupdict()['num_sites'])

    # initialize data array
    num_occupied = []

    # gather data (number of sites occupied when percolation occurs)
    it = re.finditer(r"(?P<run>[0-9]+)\t(?P<occupation>[0-9]+)\t(?P<key_site>[0-9]+)\s*", filedata)
    
    for s in it:
        num_occupied.append(int(s.groupdict()['occupation']))

    return num_sites, np.array(num_occupied, dtype=np.int64)

def get_file_list(f):
    """Retrieve a list of files to analyze from a text file."""

    # load the input
    filedata = f.read()

    # load header information
    name_dict = re.search(r"(?i)Lattice:\s*(?P<name>[^\r\n]+)", filedata)
    name = name_dict.groupdict()['name']

    side_dict = re.search(r"(?i)Side Length:\s*(?P<length>[0-9]+)\s*[\r\n]+", filedata)
    side = int(side_dict.groupdict()['length'])

    sside_dict = re.search(r"(?i)Smallest Length:\s*(?P<small>[0-9]+)\s*[\r\n]+", filedata)
    smallest_side = int(sside_dict.groupdict()['small'])

    size_dict = re.search(r"(?i)Unit Size:\s*(?P<size>[0-9]+)\s*[\r\n]+", filedata)
    size = int(size_dict.groupdict()['size'])

    R_dict = re.search(r"(?i)R:\s*(?P<R>[0-9]*[.][0-9]+)\s*[\r\n]+", filedata)
    
    if R_dict is None:
        R = None
    else:
        if 'R' in R_dict.groupdict():
            R = float(R_dict.groupdict()['R'])
        else:
            R = None

    slope_dict = re.search(r"(?i)Slope:\s*(?P<slope>[0-9]*[.][0-9]+)\s*[\r\n]+", filedata)
    
    if slope_dict is None:
        slope = None
    else:
        if 'slope' in slope_dict.groupdict():
            slope = float(slope_dict.groupdict()['slope'])
        else:
            slope = None
        
    sp_dict = re.search(r"(?i)Starting P:\s*(?P<sp>[01]*[.][0-9]+)\s*[\r\n]+", filedata)
    sp = float(sp_dict.groupdict()['sp'])

    ep_dict = re.search(r"(?i)Ending P:\s*(?P<ep>[01]*[.][0-9]+)\s*[\r\n]+", filedata)
    ep = float(ep_dict.groupdict()['ep'])

    assert ep > sp

    # initialize data array
    files = []

    # gather data (number of sites occupied when percolation occurs)
    it = re.finditer(r"[\"'](?P<path>[^\0\n\r\"']*)[\"']", filedata)

    for s in it:
        files.append(s.groupdict()['path'])

    return name, side, smallest_side, R, slope, size, sp, ep, files

def combine_data(database, size, length):
    """Adds data points until 1000 samples have been collected."""
    combined = []
    while True:
        for entry in database:
            assert entry['num_sites'] == size*length**3
            for n in entry['occupied']:
                combined.append(n)
                if len(combined) >= 1000:
                    return np.array(combined)
    assert len(combined) == 1000

def get_pdf(combined, num_sites):
    """Computes a histogram."""
    pdf = np.zeros(num_sites)
    for n in combined:
        pdf[n] += 1.
    return pdf

def get_cdf(pdf, num_sites):
    """Computes the cummulative distribution function."""
    cdf = np.zeros(num_sites)
    count = 0
    for i in xrange(num_sites):
        count += pdf[i]
        cdf[i] = count
    return cdf / np.max(cdf)

def log_factorial(n):
    """Stirling's Approximation"""
    return n*np.log(n) - n + 0.5*np.log(2.0*np.pi*n)

def binomial(n, k, p):
    """Calculates the binomial distribution for a given probability $p$,
       the number of occupied sites $k$, and the number of sites in the
       lattice $n$."""
    return np.exp(log_factorial(n) - log_factorial(k) - log_factorial(n-k) \
                + k*np.log(p) + (n-k)*np.log(1.0-p))

def convolve(cdf, P, num_sites):
    """Convolves the cdf with a binomial distribution to get probabilities."""
    K = np.array([ i+1 for i in xrange(num_sites-1) ], dtype=float)
    R = np.zeros(len(P))
    for i in xrange(len(P)):
        R[i] = np.sum(binomial(num_sites, K, P[i])*cdf[:-1], axis=0)
    return R

def main():
    """The entry point of the program"""
    # read the job list
    name, length, smallest_length, R, slope, size, Pi, Pf, filelist = get_file_list(stdin)

    # read in the data from the files
    database = []
    for f in filelist:
        num_sites, occupied_array = load_data(open(f, 'r'))
        database.append({'num_sites': num_sites, 'occupied': occupied_array})

    # collect 1000 data points
    combined = combine_data(database, size, length)

    # generate normalized cdf (range 0.0 to 1.0)
    cdf = get_cdf(get_pdf(combined, num_sites), num_sites)

    if R is None and slope is None:
        smallest_size = size*smallest_length**3
    else:
        assert R is not None
        assert slope is not None
        smallest_size = int(1.0/(np.sqrt(R*(1.-R)/1000.)/(2.*slope)))

    P = np.arange(Pi, Pf, 1.0/float(smallest_size))
    R = convolve(cdf, P, num_sites)

    for p, n in zip(P, R):
        print p, ',', n

if __name__ == "__main__":
    main()