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
bctype = {}

tau  = {}
maxt = {}
iter2skip = {}

Nx   = {}
Dx   = {}
Qx   = {}

Ny   = {}
Dy   = {}
Qy   = {}

Nz   = {}
Dz   = {}
Qz   = {}
"""

dim = 3
tau  = 0.005
maxt = 1.0
inter2skip = 10
bctype = 1
Nx   = 100
Dx   = 1.0
Qx   = 0.0
Ny   = 100
Dy   = 1.0
Qy   = 0.0
Nz   = 100
Dz   = 1.0
Qz   = 0.0

hx = 1.0 / Nx
hy = 1.0 / Ny
hz = 1.0 / Nz
maxn = int(maxt / tau)

def initcond(x,y,z):
    return np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z)

def boundcond(x,y,z):
    return 0.0

def solution(x, y, z, t):
    return np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z)*np.exp(-3.0*np.pi**2*t)

def gen():
    print("\n******************************************************************")
    print("******** Generating settings *************************************")
    print("******************************************************************\n")
    SETTINGS = SETTINGS_TEMPLATE.format(dim, tau, maxt, Nx, Dx, Qx,
                                        Ny, Dy, Qy, Nz, Dz, Qz)
    print(SETTINGS)
    with open(SETTINGS_FILE, 'w+') as out:
        out.write(SETTINGS)

    print("\n******************************************************************")
    print("******** Generating initial conditions ***************************")
    print("******************************************************************\n")
    u = np.fromfunction(lambda k, j, i: initcond(i*hx, j*hy, k*hz),
                        (Nz+1, Ny+1, Nx+1), dtype=int)
    with open(INITCOND_FILE, 'wb+') as out:
        out.write(u)
    print("Done.")

def run():
    cmd = RUN_FILE + " -p CPU -s " + SETTINGS_FILE + " -i " + \
    INITCOND_FILE + " -o " + SOLUTION_FILE
    print("\n******************************************************************")
    print("******** Running the program:\n"+cmd)
    print("******************************************************************\n")
    os.system("rm " + SOLUTION_FILE)
    os.system(cmd)

def compare():
    print("\n******************************************************************")
    print("******** Comparing with the analytic solution ********************")
    print("******************************************************************\n")
    u = np.fromfile(SOLUTION_FILE, dtype=np.double)
    u.shape = (maxn+1, Nz+1, Ny+1, Nx+1)
    errors = np.zeros(maxn+1)
    for n in xrange(maxn+1):
        err = np.abs(u[n] - np.fromfunction(lambda k, j, i: solution(i*hx, j*hy, k*hz, n*tau),
                        (Nz+1, Ny+1, Nx+1), dtype=int))
        k,j,i = np.unravel_index(err.argmax(), err.shape)
        errors[n] = err[k][j][i]
        if (n % 10 == 0):
            print("Step: {} | Max error is {} on ({},{},{})".format(n, err[k][j][i], k,j,i))
    print("Maximum error is {} on step {}".format(errors.max(), errors.argmax()))
    print("Done.")
    return

def main():
    if (len(sys.argv) != 2):
        gen()
        run()
        compare()
    elif (sys.argv[1] == "gen"):
        gen()
    elif (sys.argv[1] == "run"):
        run()
    elif (sys.argv[1] == "cmp"):
        compare()

if __name__ == "__main__":
    main()
