#!/usr/bin/python

import numpy as np
import os
import stat
import sys

if (len(sys.argv) != 2):
    print("No output dir provided")
    sys.exit(0)

if not os.path.exists(sys.argv[1]):
    os.mkdir(sys.argv[1])


################ Generate IC ####################
Nx = 100
u = np.zeros(Nx+1)
for i in xrange(Nx+1):
    u[i] = np.sin(np.pi*i/Nx)

path = sys.argv[1]+"/u0.bin"

with open(path, 'wb+') as out:
    out.write(u)
#################################################

########### Generate info #######################
INFO = ur"""#+STARTUP: latexpreview

\begin{equation}
x=\sqrt{b}
\end{equation}

If $a^2=b$ and \( b=2 \), then the solution must be
either $$ a=+\sqrt{2} $$ or \[ a=-\sqrt{2} \].
"""

path = sys.argv[1]+"/info.org"

with open(path, 'w+') as out:
    out.write(INFO)
#################################################

########## Generate run script ##################
RUN_SH = ur"""#!/bin/bash

../../bin/heat settings u0.bin
"""

path = sys.argv[1]+"/run.sh"

with open(path, 'w+') as out:
    out.write(RUN_SH)

os.chmod(path, stat.S_IRWXU);
#################################################
