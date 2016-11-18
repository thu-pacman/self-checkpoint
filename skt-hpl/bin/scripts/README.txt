##### REQUIREMENTS #####
SKT-HPL has been tested with IntelCompiler (icc) + IntelMPI/MPICH + SLURM.

##### INSTALL #####
We provide a precompiled executable file xhpl in SKT-HPL/bin/nft.
But it may be incompatible with your system.

If you need to change compiler or BLAS path, modify Make.nft
$> vim Make.nft

In the root directory of SKT-HPL, type
$> make arch=nft

After that, you can the executable file in SKT-HPL/bin/nft, and this is the run directory.

##### NODELIST #####
Prepare two nodelist files: worklist and sparelist.
The *worklist* is a valid MPI ranklist. It specifies where to run SKT-HPL at first. 
If there is a failure and some nodes are no longer available, then SKT-HPL will try to run on new nodes, which are listed in *sparelist*.

##### CONFIG HPL.dat #####
HPL.dat is the HPL benchmark input file. HPL.dat for SKT-HPL is the same with standard HPL but with some constraints.
The matrix size N must be multiple of NB*P and multiple of NB*Q, i.e.
          N | NB*GCD(P, Q)

##### RUN SKT-HPL #####
Instead of directly submitting SKT-HPL, we use *hpl-daemon.sh* as a automatic driver. It tries to run SKT-HPL in a loop until it completes successfully.
The parameters for SKT-HPL are listed as below (and also described in hpl-daemon.sh):

* RPN: How many MPI processes on a node.
* TSIZE: Group size described in the paper, 8 and 16 are good choices.
* SNAPSHOT: Checkpoint interval. It means *how many iterations*, NOT *how many seconds*.
* REUSEMAT: If set to 0, generate a new problem, otherwise, read from checkpoint. Should be 0 for a fresh run and set to 1 for recovery.
* ONLY: For debugging purpose. If set to 1, only one checkpoint will be made. If set to 0, will do checkpoint periodically.
* RECOVER: For debugging purpose. If set to 0, checkpoint will not be read after restart.
* RUN_TIME_LIMIT: For debugging purpose. A complete SKT-HPL could last for hours, but it will exit after RUN_TIME_LIMIT seconds. 

Example for normal run:
RPN=xxx TSIZE=8 SNAPSHOT=30 REUSEMAT=1 ONLY=0 RECOVER=1 RUN_TIME_LIMIT=600

##### SHM ISSUES #####
SKT-HPL uses SysV* shared memory. Although we can be allocated and freed inside application code, but we may be unable to do so under failure (The code for free will not be executed after a fault).
Addtional scripts are used to clean environment.
*clean.sh* and *clr.sh* are for this purpose. They are invoked by *hpl-daemon.sh*

##### FAULT INJECTION #####
Currently, faults are inject manually.
There are two kinds of fault in our experiment, failure in computing and failure in checkpoint updating.

While SKT-HPL is computing, it prints something like
$> Column=000896 Fraction=0.015 Mflops=607545.31

One can power-off a single node, or kill MPI processes on a node then clean their shm area as if the node is restart and no data in memory.
$> ssh to a compute node
$> cd to the run dir
$> killall xhpl
$> ./clr.sh

After a while, the job will exit with an error code, then hpl-daemon.sh will try to restart SKT-HPL.
$> srun: tasks 16-127: running      
$> srun: tasks 0-15: exited abnormally
$> srun: Job step aborted: Waiting up to 32 seconds for job step to finish.
$> MPI FAILED Tue Nov 15 11:49:20 CST 2016
$> No xhpl running, going to (re)start

When SKT-HPL starts to update checkpoint, it prints 
$> partially overwrite chkpt data

A fault should be inject here (the injection operation is same as above), before the updating finish. 
When SKT-HPL finish updating, it prints something like
$> SNAPSHOT    131 MB/rank, MEMSET  0.02 sec, ENCODE  0.56 sec ......

