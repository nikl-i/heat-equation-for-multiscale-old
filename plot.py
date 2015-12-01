#!/usr/bin/python3
from subprocess import Popen, PIPE, call
import platform
import time
import sys

PLOTSCRIPT = """
load 'style.gp'
set terminal 'pngcairo'
set output '{}'
set key outside
set grid

set xrange [0.0:1.0]
set yrange [0.0:1.0] 
set xlabel 'x'
set ylabel 'U'
set title '{}'

plot '{}' using 1:2 w lp ls 2 notitle
quit
"""

class gnuplot:
    """ Gnuplot calling interface """
    def run(script):
        if platform.system() == 'Windows' or platform.system() == 'Linux':
            if platform.system() == 'Windows':
                gnuplot = r'C:\Program Files (x86)\gnuplot\bin\pgnuplot'
            if platform.system() == 'Linux':
                gnuplot = 'gnuplot'
        plot = Popen([gnuplot, '-persist'],
                     stdin = PIPE, stdout = PIPE, stderr = PIPE)      
        plot.stdin.write(bytes(script, 'UTF-8'))  
        plot.stdin.flush()    

print("Plotting...")
startTime = time.time()

for i in range(1,len(sys.argv)):
    gnuplot.run(PLOTSCRIPT.format(sys.argv[i].replace(".out",".png"), 
                                  sys.argv[i], sys.argv[i]))
    if (i % 1000 == 0): print("Step: {} / {}".format(i, len(sys.argv)-1)) 
print('Done! ({}s)'.format(time.time() - startTime))





#PLOTHIST2DSCRIPT = """
#set cbrange [0.0:0.1]
#count and plot
#set pm3d map 
#set pm3d interpolate 0,0
#splot '-' matrix
#{}
#EOF
#quit
#"""
