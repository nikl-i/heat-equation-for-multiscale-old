#!/usr/bin/python

import numpy as np
import os
import sys

WORKDIR = "/tmp/heat/"
RUN_FILE      = WORKDIR + "bin/heat"
SETTINGS_FILE = WORKDIR + "src/settings"
INITCOND_FILE = WORKDIR + "u0.bin"
SOLUTION_FILE = WORKDIR + "u.bin"
BOUNDCOND_QX0_FILE = WORKDIR + "qx0.bin"
BOUNDCOND_QX1_FILE = WORKDIR + "qx1.bin"
BOUNDCOND_PX0_FILE = WORKDIR + "px0.bin"
BOUNDCOND_PX0_FILE = WORKDIR + "px1.bin"

SETTINGS_DIM = 1
SETTINGS_NUM_OF_ITERATIONS_TO_SKIP = 10
SETTINGS_MAXT = 1.0
SETTINGS_TAU  = 0.005

SETTINGS_DX = 1.0
SETTINGS_NX  = 100
SETTINGS_BCX0_TYPE = 1
SETTINGS_BCX1_TYPE = 1

SETTINGS_DY = 1.0
SETTINGS_NY = 100
SETTINGS_BCY0_TYPE = 1
SETTINGS_BCY1_TYPE = 1

SETTINGS_DZ = 1.0
SETTINGS_NZ = 100
SETTINGS_BCZ0_TYPE = 1
SETTINGS_BCZ1_TYPE = 1

SETTINGS_BCQX0_FUNCTION = 0.0
SETTINGS_BCQX1_FUNCTION = 0.0
SETTINGS_BCPX0_FUNCTION = 0.0
SETTINGS_BCPX1_FUNCTION = 0.0

SETTINGS_BCQY0_FUNCTION = 0.0
SETTINGS_BCQY1_FUNCTION = 0.0
SETTINGS_BCPY0_FUNCTION = 0.0
SETTINGS_BCPY1_FUNCTION = 0.0

SETTINGS_BCQZ0_FUNCTION = 0.0
SETTINGS_BCQZ1_FUNCTION = 0.0
SETTINGS_BCPZ0_FUNCTION = 0.0
SETTINGS_BCPZ1_FUNCTION = 0.0

SETTINGS_TEMPLATE = ur"""
#define SETTINGS_PROBLEM_TYPE Problem{}d
#define SETTINGS_DIM {}
#define SETTINGS_NUM_OF_ITERATIONS_TO_SKIP {}
#define SETTINGS_MAXT {}
#define SETTINGS_TAU {}

#define SETTINGS_DX {}
#define SETTINGS_NX {}
#define SETTINGS_BCX0_TYPE {}
#define SETTINGS_BCX1_TYPE {}

#define SETTINGS_DY {}
#define SETTINGS_NY {}
#define SETTINGS_BCY0_TYPE {}
#define SETTINGS_BCY1_TYPE {}

#define SETTINGS_DZ {}
#define SETTINGS_NZ {}
#define SETTINGS_BCZ0_TYPE {}
#define SETTINGS_BCZ1_TYPE {}

#define SETTINGS_BCQX0_FUNCTION(n,j,k) {}
#define SETTINGS_BCQX1_FUNCTION(n,j,k) {}
#define SETTINGS_BCPX0_FUNCTION(n,j,k) {}
#define SETTINGS_BCPX1_FUNCTION(n,j,k) {}

#define SETTINGS_BCQY0_FUNCTION(n,i,k) {}
#define SETTINGS_BCQY1_FUNCTION(n,i,k) {}
#define SETTINGS_BCPY0_FUNCTION(n,i,k) {}
#define SETTINGS_BCPY1_FUNCTION(n,i,k) {}

#define SETTINGS_BCQZ0_FUNCTION(n,i,j) {}
#define SETTINGS_BCQZ1_FUNCTION(n,i,j) {}
#define SETTINGS_BCPZ0_FUNCTION(n,i,j) {}
#define SETTINGS_BCPZ1_FUNCTION(n,i,j) {}
"""

hx = 1.0 / SETTINGS_NX
maxstep = int(SETTINGS_MAXT / SETTINGS_TAU)
maxm = maxstep / SETTINGS_NUM_OF_ITERATIONS_TO_SKIP #size of the array with solution

def initcond(x):
    return np.sin(np.pi*x)

def solution(x, t):
    return np.sin(np.pi*x)*np.exp(- np.pi**2 * t)

def gen():
    print("\n******************************************************************")
    print("***************** Generating settings ****************************")
    print("******************************************************************\n")

    SETTINGS = SETTINGS_TEMPLATE.format(SETTINGS_DIM, SETTINGS_DIM,
                                    SETTINGS_NUM_OF_ITERATIONS_TO_SKIP,
                                    SETTINGS_MAXT, SETTINGS_TAU,
                                    SETTINGS_DX, SETTINGS_NX,
                                    SETTINGS_BCX0_TYPE, SETTINGS_BCX1_TYPE,
                                    SETTINGS_DY, SETTINGS_NY,
                                    SETTINGS_BCY0_TYPE, SETTINGS_BCY1_TYPE,
                                    SETTINGS_DZ, SETTINGS_NZ,
                                    SETTINGS_BCZ0_TYPE, SETTINGS_BCZ1_TYPE,
                                    SETTINGS_BCQX0_FUNCTION,
                                    SETTINGS_BCQX1_FUNCTION,
                                    SETTINGS_BCPX0_FUNCTION,
                                    SETTINGS_BCPX1_FUNCTION,
                                    SETTINGS_BCQY0_FUNCTION,
                                    SETTINGS_BCQY1_FUNCTION,
                                    SETTINGS_BCPY0_FUNCTION,
                                    SETTINGS_BCPY1_FUNCTION,
                                    SETTINGS_BCQZ0_FUNCTION,
                                    SETTINGS_BCQZ1_FUNCTION,
                                    SETTINGS_BCPZ0_FUNCTION,
                                    SETTINGS_BCPZ1_FUNCTION)

    print(SETTINGS)
    with open(SETTINGS_FILE, 'w+') as out:
        out.write(SETTINGS)

    print("\n******************************************************************")
    print("***************** Generating initial conditions ******************")
    print("******************************************************************\n")
    u = np.fromfunction(lambda i: initcond(i*hx), (Nx+1,))
    with open(INITCOND_FILE, 'wb+') as out:
        out.write(u)
    return


def build():
    print("***************** Building the program: ****************************")
    cmd = "cp -r ./src " + WORKDIR + "/src"
    print(cmd)
    os.system(cmd)
    gen()
    cmd = "cd " + WORKDIR + "; make -B debug; cd -"
    print(cmd)
    os.system(cmd)

def run():
    cmd = RUN_FILE
    print("\n*****************************************************************")
    print("***************** Running the program: ****************************")
    print(cmd)
    print("******************************************************************\n")
    os.system("rm " + SOLUTION_FILE)
    os.system(cmd)

def compare():
    print("\n******************************************************************")
    print("********** Comparing with the analytic solution ******************")
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
    main()
