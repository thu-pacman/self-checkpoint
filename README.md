This is the repo for Self-Checkpoint project.
Self-checkpoint is an in-memory checkpoint method using less memory space (about 1/2), so more can be left for applications.
There is a plan to build some APIs, but currently this method is hard coded inside applications.

## SKT-HPL
SKT-HPL is a fault-tolerant HPL implementation based on self-checkpoint.

It is based on HPL benchmark (version 2.2). Modification can be found in 
* include/hpl\_grid.h
* src/pgesv/HPL\_pdgesv0.c
* src/grid/HPL\_grid\_init.c
* testing/ptest/HPL\_pdtest.c

## Step-by-step tutorial
Scripts are provided in SKT-HPL/bin. A README.txt is also provided there.
