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

## Install and Run
The installation of SKT-HPL is the same with original HPL.
Running scripts are provided in skt-hpl/bin/scripts. 
A step-by-step tutorial in PDF format is provided. There is also a README.txt in skt-hpl/bin/scripts, but not so deatiled.
