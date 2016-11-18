#!/bin/bash
# used to clean the SysV* shm
for key in `ipcs -m | tail -n +4 | awk '{print $1}'`; do ipcrm -M $key &2>/dev/null; done
sleep 1
for key in `ipcs -m | tail -n +4 | awk '{print $1}'`; do ipcrm -M $key &2>/dev/null; done
sleep 1
for key in `ipcs -m | tail -n +4 | awk '{print $1}'`; do ipcrm -M $key &2>/dev/null; done
