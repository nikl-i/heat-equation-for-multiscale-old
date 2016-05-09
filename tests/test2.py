#!/usr/bin/python

import numpy as np
import os
import sys

WORKDIR = "/tmp/heat/"
RUN_FILE      = WORKDIR + "bin/heat"
SETTINGS_FILE = WORKDIR + "settings"
INITCOND_FILE = WORKDIR + "u0.bin"
SOLUTION_FILE = WORKDIR + "u.bin"

SETTINGS_TEMPLATE = ur"""dim = {}
tau  = {}
maxt = {}
Nx   = {}
Dx   = {}
Qx   = {}

Ny   = {}
Dy   = {}
Qy   = {}
"""

dim = 2

tau  = 0.005
maxt = 1.0

Nx   = 100
Dx   = 1.0
Qx   = 0.0

Ny   = 100
Dy   = 1.0
Qy   = 0.0

hx = 1.0 / Nx
hy = 1.0 / Ny
maxn = int(maxt / tau)

def initcond(x,y):
    return np.sin(np.pi*x)*np.sin(np.pi*y)

def solution(x, y, t):
    return np.sin(np.pi*x)*np.sin(np.pi*y)*np.exp(-2.0*np.pi**2 * t)

def gen():
    ################ Generate settings ##############
    print("\n******************************************************************")
    print("******** Generating settings")
    print("******************************************************************\n")
    SETTINGS = SETTINGS_TEMPLATE.format(dim, tau, maxt, Nx, Dx, Qx, Ny, Dy, Qy)
    print(SETTINGS)
    with open(SETTINGS_FILE, 'w+') as out:
        out.write(SETTINGS)
    #################################################

    ################ Generate IC ####################
    print("\n******************************************************************")
    print("******** Generating initial conditions ")
    print("******************************************************************\n")
    totalSize = (Nx+1)*(Ny+1)
    u = np.zeros(totalSize)
    for j in xrange(Ny+0):
    	for i in xrange(Nx+1):
            u[j*(Nx+1) + i] = initcond(i*hx,j*hy)

    with open(INITCOND_FILE, 'wb+') as out:
        out.write(u)
    #################################################
    print("Done.")

def run():
    totalSize = (Nx+1)*(Ny+1)
    cmd = RUN_FILE + " -p CPU -s " + SETTINGS_FILE + " -i " + \
    INITCOND_FILE + " -o " + SOLUTION_FILE
    print("\n******************************************************************")
    print("******** Running the program:\n"+cmd)
    print("******************************************************************\n")
    os.system("rm " + SOLUTION_FILE)
    os.system(cmd)
    u = np.fromfile(SOLUTION_FILE, dtype=np.double)

    print("\n******************************************************************")
    print("******** Comparing with the analytic solution")
    print("******************************************************************\n")

    maxerr = 0.0
    maxerr_t = 0.0
    maxerr_x = 0.0
    maxerr_y = 0.0

    for n in xrange(maxn+1):
        t = n*tau
        locerr = 0.0
        locerr_t = 0.0
        locerr_x = 0.0
        locerr_y = 0.0

        for j in xrange(Ny+1):
            for i in xrange(Nx+1):
               x = i*hx
               y = j*hy
               err = abs(u[totalSize*n + j*(Nx+1) + i] - solution(x, y, t) )
               # if (i == 2 and j == 1):
               #     print("Interation t={:6.4} x={:6.4} y={:6.4} : n={}, i={}, j={}".format(t,x,y,n,i,j))
               #     print("err = {} = abs( res={} - sol={} )\n".format(err, u[totalSize*n + j*(Ny+1) + i], solution(x, y, t)))

               if (err > locerr):
                   locerr = err
                   locerr_x = x
                   locerr_y = y
                   locerr_t = t

            if (locerr > maxerr):
                maxerr = locerr
                maxerr_x = x
                maxerr_y = y
                maxerr_t = t
	#       if (n % 100 == 0):
        print("Interation t={:12.8f}: locerr = {} (x={:10.9f}, y={:10.9f})".format(t, locerr, locerr_x, locerr_y))
    print("maxerr = {} | ({},{},{})".format(maxerr, maxerr_t, maxerr_x, maxerr_y))
    print("Done.")
    return

def main():
    if (len(sys.argv) != 2):
        gen()
        run()
    elif (sys.argv[1] == "gen"):
        gen()
    elif (sys.argv[1] == "run"):
        run()

if __name__ == "__main__":
    main()
