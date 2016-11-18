#!/bin/bash
################################################
# Run this script instead of HPL.
# SKT-HPL will be run in a loop until finish.
################################################

REUSE_MAT_N=0
while true; do
xcnt=`ps -ef | grep 'xhpl' | wc -l`
if [ $xcnt -le 1 ]; then
   echo 'No xhpl running, going to (re)start'
# Generate a nodelist
# Remove bad nodes from worklist, and replace them by good nodes in sparelist.
   rm -f newlist
   sline=1
   for node in `cat worklist`; do
       ping -c 1 $node 2>&1 1>/dev/null
       if [ $? -eq 0 ]; then
           #echo "$node is fine"
           echo $node >> newlist
       else
           rnode=`sed -n ${sline}p sparelist`
           echo "Replace $node --> $rnode"
           echo $rnode >> newlist
           ((sline=$sline+1))
       fi
   done
   tail -n +$sline sparelist > tmpsparelist
   mv tmpsparelist sparelist
   mv newlist worklist

# Run SKT-HPL with new nodelist
# RPN: How many MPI processes on a node.
# TSIZE: Group size described in the paper, 8 and 16 are good choices.
# SNAPSHOT: Checkpoint interval. It means *how many iterations*, NOT *how many seconds*.
# REUSEMAT: If set to 0, generate a new problem, otherwise, read from checkpoint. Should be 0 for a fresh run and set to 1 for recovery.
# ONLY: For debugging purpose. If set to 1, only one checkpoint will be made. If set to 0, will do checkpoint periodically.
# RECOVER: For debugging purpose. If set to 0, checkpoint will not be read after restart.
# RUN_TIME_LIMIT: For debugging purpose. A complete SKT-HPL could last for hours, but it will exit after RUN_TIME_LIMIT seconds. 

   RPN=16 TSIZE=8 SNAPSHOT=30 REUSEMAT=$REUSE_MAT_N ONLY=0 RECOVER=1 RUN_TIME_LIMIT=3000 srun -p work -n 128 --nodelist=./worklist --ntasks-per-node=16 ../skt/xhpl

# if a fail run return 0, we need to repeat the running for several times.
   if [ $? -eq 0 ]; then
       echo ""
       echo 'MPI JOB DONE '`date`

       ncnt=`cat worklist | wc -l`
       echo "Cleanning SysV* shm on $ncnt nodes"
       srun -n $ncnt --nodelist=./worklist --ntasks-per-node=1 ./clr.sh
       echo "Environment Cleaned"
       RSM=0
       exit
   fi
   REUSE_MAT_N=1
fi
echo 'MPI FAILED '`date`
done
