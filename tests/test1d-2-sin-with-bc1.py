#!/usr/bin/python

import numpy as np
import os
import sys

WORKDIR = "/tmp/heat/"
RUN_FILE      = WORKDIR + "bin/heat"
SETTINGS_FILE = "src/settings.h"

SETTINGS_INITCOND_FILE = WORKDIR + "u0.bin"
SETTINGS_SOLUTION_FILE = WORKDIR + "u.bin"

SETTINGS_DIM = 1
SETTINGS_NUM_OF_ITERATIONS_TO_SKIP = 10
SETTINGS_MAXT = 1.0
SETTINGS_TAU  = 0.005

SETTINGS_DX = 1.0
SETTINGS_NX  = 100
SETTINGS_BCX0_TYPE = 1
SETTINGS_BCX1_TYPE = 1

SETTINGS_BCQX0_FUNCTION = "0.0"
SETTINGS_BCQX1_FUNCTION = "0.0"
SETTINGS_BCPX0_FUNCTION = "0.0"
SETTINGS_BCPX1_FUNCTION = "0.0"


SETTINGS_TEMPLATE = ur"""#define SETTINGS_INITCOND_FILE "{}"
#define	SETTINGS_SOLUTION_FILE "{}"

#define SETTINGS_PROBLEM_TYPE Problem{}d

#define SETTINGS_NUM_OF_ITERATIONS_TO_SKIP {}
#define SETTINGS_MAXT {}
#define SETTINGS_TAU {}

#define SETTINGS_DX {}
#define SETTINGS_NX {}

#define SETTINGS_BCX0_TYPE {}
#define SETTINGS_BCQX0_FUNCTION(n,j,k) {}
#define SETTINGS_BCPX0_FUNCTION(n,j,k) {}

#define SETTINGS_BCX1_TYPE {}
#define SETTINGS_BCQX1_FUNCTION(n,j,k) {}
#define SETTINGS_BCPX1_FUNCTION(n,j,k) {}
"""

def initcond(x):
    return np.sin(np.pi*x)

def solution(x, t):
    return np.sin(np.pi*x)*np.exp(- np.pi**2 * t)

def gen():
    print("***************** Generating settings ****************************")
    SETTINGS = SETTINGS_TEMPLATE.format(SETTINGS_INITCOND_FILE,
                                        SETTINGS_SOLUTION_FILE,
                                        SETTINGS_DIM,
                                        SETTINGS_NUM_OF_ITERATIONS_TO_SKIP,
                                        SETTINGS_MAXT, SETTINGS_TAU,
                                        SETTINGS_DX, SETTINGS_NX,
                                        SETTINGS_BCX0_TYPE,
                                        SETTINGS_BCQX0_FUNCTION,
                                        SETTINGS_BCPX0_FUNCTION,
                                        SETTINGS_BCX1_TYPE,
                                        SETTINGS_BCQX1_FUNCTION,
                                        SETTINGS_BCPX1_FUNCTION)
    print(SETTINGS)
    with open(SETTINGS_FILE, 'w+') as out:
        out.write(SETTINGS)

    print(SETTINGS_FILE)
    print("***************** Generating initial conditions ******************")
    u = np.fromfunction(lambda i: initcond(i*hx), (SETTINGS_NX+1,))
    with open(SETTINGS_INITCOND_FILE, 'wb+') as out:
        out.write(u)
    return


def build():
    print("***************** Building the program: ****************************")
    cmd = "mkdir " + WORKDIR
    print(cmd)
    os.system(cmd)

    gen()

    cmd = "make -B"
    print(cmd)
    os.system(cmd)

def run():
    cmd = RUN_FILE
    print("***************** Running the program: ****************************")
    print(cmd)
    os.system("rm " + SETTINGS_SOLUTION_FILE)
    os.system(cmd)

def compare():
    print("********** Comparing with the analytic solution ******************")
    u = np.fromfile(SETTINGS_SOLUTION_FILE, dtype=np.double)
    u.shape = (maxm+1, SETTINGS_NX+1)
    errors = np.zeros(maxm+1)
    for m in xrange(maxm+1):
        err = np.abs(u[m] - np.fromfunction(lambda i: solution(i*hx,m*taum), (SETTINGS_NX+1,)))
        i = np.unravel_index(err.argmax(), err.shape)
        errors[m] = err[i]
        print("Step: {} | Max error is {} on {}".format(m, err[i], i))
    print("Maximum error is {} on step {}".format(errors.max(), errors.argmax()))
    print("Done.")
    return

def main():
    if (len(sys.argv) != 2):
        build()
        run()
        compare()
    elif (sys.argv[1] == "gen"):
        gen()
    elif (sys.argv[1] == "run"):
        run()
    elif (sys.argv[1] == "cmp"):
        compare()

if __name__ == "__main__":

    hx = 1.0 / SETTINGS_NX
    maxstep = int(SETTINGS_MAXT / SETTINGS_TAU)
    maxm = maxstep / SETTINGS_NUM_OF_ITERATIONS_TO_SKIP #size of the array with solution
    taum = SETTINGS_NUM_OF_ITERATIONS_TO_SKIP*SETTINGS_TAU
    main()
