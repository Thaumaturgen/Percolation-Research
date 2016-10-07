import subprocess

for i in xrange(5):
    p = subprocess.Popen(['python', './process%d/diamond_site_perc128.py'%(i+1)])