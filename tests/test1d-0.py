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
iter2skip = {}

Nx   = {}
Dx   = {}
Qx   = {}
bcx0 = {}
bcx1 = {}
"""

dim = 1
bcx0 = 0
bcx1 = 0

tau  = 0.005
maxt = 1.0
iter2skip = 10

Nx   = 100
Dx   = 1.0
Qx   = 0.0

hx = 1.0 / Nx
maxm = int(maxt / tau) / iter2skip #size of the array with solution

def initcond(x):
    return np.sin(np.pi*x)

def solution(x, t):
    return np.sin(np.pi*x)*np.exp(- np.pi**2 * t)

def gen():
    print("\n******************************************************************")
    print("******** Generating settings *************************************")
    print("******************************************************************\n")
    SETTINGS = SETTINGS_TEMPLATE.format(dim, tau, maxt, iter2skip, Nx, Dx, Qx, bcx0, bcx1)
    print(SETTINGS)
    with open(SETTINGS_FILE, 'w+') as out:
        out.write(SETTINGS)

    print("\n******************************************************************")
    print("******** Generating initial conditions ***************************")
    print("******************************************************************\n")
    u = np.fromfunction(lambda i: initcond(i*hx), (Nx+1,))
    with open(INITCOND_FILE, 'wb+') as out:
        out.write(u)
    print("Done.")
    return

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
    u.shape = (maxm+1, Nx+1)
    errors = np.zeros(maxm+1)
    for m in xrange(maxm+1):
        err = np.abs(u[m] - np.fromfunction(lambda i: solution(i*hx,m*iter2skip*tau), (Nx+1,)))
        i = np.unravel_index(err.argmax(), err.shape)
        errors[m] = err[i]
        print("Step: {} | Max error is {} on {}".format(m, err[i], i))
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
