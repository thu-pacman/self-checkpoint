#!/bin/bash
for node in `cat worklist`; do
   ping -c 1 $node 2>&1 1>/dev/null
   if [ $? -eq 0 ]; then
       echo $node >> newlist
   else
       echo "Replace $node with a new node"
   fi
done
