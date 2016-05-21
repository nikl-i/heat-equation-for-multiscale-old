#!/usr/bin/python3
from subprocess import Popen, PIPE, call
import platform
import time
import sys
import os
import numpy as np

WORKDIR = "/tmp/heat/"
RUN_FILE      = WORKDIR + "bin/heat"
SETTINGS_FILE = "src/settings.h"
OUTPUTDIR = WORKDIR

PLOTSCRIPT1D = """
set encoding utf8
unset key

# Color
set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red
set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2 # --- green
set style line 3 lc rgb '#800000' lt 1 lw 2
set style line 4 lc rgb '#ff0000' lt 1 lw 2
set style line 5 lc rgb '#ff4500' lt 1 lw 2
set style line 6 lc rgb '#ffa500' lt 1 lw 2
set style line 7 lc rgb '#006400' lt 1 lw 2
set style line 8 lc rgb '#0000ff' lt 1 lw 2
set style line 9 lc rgb '#9400d3' lt 1 lw 2

# Grey
set style line 11 lc rgb 'gray30' lt 1 lw 2
set style line 12 lc rgb 'gray40' lt 1 lw 2
set style line 13 lc rgb 'gray70' lt 1 lw 2
set style line 14 lc rgb 'gray90' lt 1 lw 2
set style line 15 lc rgb 'black' lt 1 lw 1.5
set style line 16 lc rgb "black" lt 4 lw 2

# Borders, etc.
set style line 21 lc rgb 'black' lt 1 lw 1.5
set border 3 back ls 21
set tics nomirror
set style line 22 lc rgb 'grey20' lt 0 lw 1
set grid back ls 22

# Palette
set palette maxcolors 3
set palette defined ( 0 '#8b1a0e',\
                      1 '#5e9c36',\
                      2 '#800000' )

set terminal '{}'
set output '{}'
set key outside
set grid

set title '{}'
set xlabel '{}'
set ylabel '{}'
set xrange [{}:{}]
set yrange [{}:{}]

plot '-' using 0:1 w lp ls 2 notitle
{}

EOF
quit
"""

PLOTSCRIPT2D  = """
set terminal '{}'
set output '{}'

set title '{}'
set xlabel '{}'
set ylabel '{}'
set xrange [{}:{}]
set yrange [{}:{}]
set cbrange [{}:{}]

set pm3d map
splot '-' matrix
{}

EOF
quit
"""

PLOTSCRIPT3D = """
set terminal '{}' size 1024,1024
set output '{}'

set multiplot layout 3,3 title '{}'

set xrange [{}:{}]
set yrange [{}:{}]
set cbrange [{}:{}]
set pm3d map

set xlabel 'x'
set ylabel 'y'

set title '{}'
splot '-' matrix
{}

EOF

set title '{}'
splot '-' matrix
{}

EOF

set title '{}'
splot '-' matrix
{}

EOF

set xlabel 'x'
set ylabel 'z'

set title '{}'
splot '-' matrix
{}

EOF

set title '{}'
splot '-' matrix
{}

EOF

set title '{}'
splot '-' matrix
{}

EOF

set xlabel 'y'
set ylabel 'z'

set title '{}'
splot '-' matrix
{}

EOF

set title '{}'
splot '-' matrix
{}

EOF

set title '{}'
splot '-' matrix
{}

EOF

#unset multiplot
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

def plot1d(step_first, step_last):
    u.shape = (maxn+1, Nx+1)
    for n in range(step_first, step_last):
        SCRIPT = PLOTSCRIPT1D.format("png", OUTPUTDIR+"plot-t={:06.3f}.png".format(n*tau),
                                     "t = {:6.3f}".format(n*tau), "x", "U",
                                     0.0, Nx, u.min(), u.max(),
                                     np.array2string(u[n,:]).translate(trnstbl).replace("  ", " ").replace(" ", "\n"))
        gnuplot.run(SCRIPT)
    return

def plot2d(step_first, step_last):
    u.shape = (maxn+1, Ny+1, Nx+1)
    for n in range(step_first, step_last):
        SCRIPT = PLOTSCRIPT2D.format("png", OUTPUTDIR+"plot-t={:06.3f}.png".format(n*tau),
                                    "t = {:6.3f} ps".format(n*tau), "x", "y",
                                     0.0, Nx, 0.0, Ny, u.min(), u.max(),
                                     np.array2string(u[n,:,:]).translate(trnstbl))
        gnuplot.run(SCRIPT)
    return

def plot3d(step_first, step_last):
    u.shape = (maxn+1, Nz+1, Ny+1, Nx+1)
    startTime = time.time()
    for n in range(step_first, step_last):
        SCRIPT = PLOTSCRIPT3D.format(
            "png", OUTPUTDIR+"plot-t={:06.3f}.png".format(n*tau),
            "t = {:6.3f} ps".format(n*tau), 0, Nx, 0, Ny, u.min(), u.max(),
            "z = 0.0",
            np.array2string(u[n,0,:,:]).translate(trnstbl),
            "z = {}".format(int(Nz/2)/float(Nz)),
            np.array2string(u[n,int(Nz/2),:,:]).translate(trnstbl),
            "z = {}".format(1.0),
            np.array2string(u[n,Nz,:,:]).translate(trnstbl),
            "y = 0.0",
            np.array2string(u[n,:,0,:]).translate(trnstbl),
            "y = {}".format(int(Ny/2)/float(Ny)),
            np.array2string(u[n,:,int(Ny/2),:]).translate(trnstbl),
            "y = {}".format(1.0),
            np.array2string(u[n,:,Ny,:]).translate(trnstbl),
            "x = 0.0",
            np.array2string(u[n,:,:,0]).translate(trnstbl),
            "x = {}".format(int(Nx/2)/float(Nx)),
            np.array2string(u[n,:,:,int(Nx/2)]).translate(trnstbl),
            "x = {}".format(1.0),
            np.array2string(u[n,:,:,Nx]).translate(trnstbl))
        gnuplot.run(SCRIPT)
        print("Step: {} ({} ps / {} ps) {} seconds lasts".format(
            n, n*tau,  step_last*tau, time.time() - startTime))

    return

def main():
    step_first = 0
    step_last = maxn
    print(maxn)

    print("Plotting...")
    startTime = time.time()
    if (dim == 1):
        plot1d(step_first, step_last)
    elif (dim == 2):
        plot2d(step_first, step_last)
    elif (dim == 3):
        plot3d(step_first, step_last)
    print('Done! ({}s)'.format(time.time() - startTime))
    return




if __name__ == "__main__":
    if (len(sys.argv) != 2):
        trnstbl = {ord(c): None for c in "[]\""}
        problemType = None
        SETTINGS_NUM_OF_ITERATIONS_TO_SKIP = 1;
        with open(SETTINGS_FILE, 'r') as f:
            for line in f:
                tokens = line.split()
                if (len(tokens) == 3):
                    if (tokens[0] == "#define"):
                        if (tokens[1] == "SETTINGS_PROBLEM_TYPE"):
                            if (tokens[2] == "Problem1d"):
                                dim = 1
                            elif (tokens[2] == "Problem2d"):
                                dim = 2
                            elif (tokens[2] == "Problem3d"):
                                dim = 3
                        elif (tokens[1] == "SETTINGS_NUM_OF_ITERATIONS_TO_SKIP"):
                            SETTINGS_NUM_OF_ITERATIONS_TO_SKIP = int(tokens[2])
                        elif (tokens[1] == "SETTINGS_MAXT"):
                            SETTINGS_MAXT = float(tokens[2])
                        elif (tokens[1] == "SETTINGS_TAU"):
                            SETTINGS_TAU = float(tokens[2])
                        elif (tokens[1] == "SETTINGS_NX"):
                            Nx = int(tokens[2])
                        elif (tokens[1] == "SETTINGS_NY"):
                            Ny = int(tokens[2])
                        elif (tokens[1] == "SETTINGS_NZ"):
                            Nz = int(tokens[2])
                        elif (tokens[1] == "SETTINGS_SOLUTION_FILE"):
                            SETTINGS_SOLUTION_FILE = tokens[2].translate(trnstbl)

        maxn = int(SETTINGS_MAXT/(SETTINGS_TAU*SETTINGS_NUM_OF_ITERATIONS_TO_SKIP))
        tau = SETTINGS_NUM_OF_ITERATIONS_TO_SKIP*SETTINGS_TAU

        print(tau, maxn, Nx, SETTINGS_SOLUTION_FILE)
        u = np.fromfile(SETTINGS_SOLUTION_FILE, dtype=np.double)
        np.set_printoptions(threshold=np.nan, linewidth = np.nan, formatter={'float': '{: 012.8f}'.format})
        main()
    elif (sys.argv[1] == "clean"):
        print("***************** Removing plot files: ****************************")
        cmd = "rm " + OUTPUTDIR + "*.png"
        print(cmd)
        os.system(cmd)

#plot text
#plot all
#print one
#print all
#print with step
#video
