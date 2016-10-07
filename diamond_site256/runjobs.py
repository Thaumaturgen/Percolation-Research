import subprocess

for i in xrange(7):
    p = subprocess.Popen(['python', './process%d/diamond_site_perc256.py'%(i+1)])
