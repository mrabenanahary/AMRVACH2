#!/bin/bash
for i in {0..5}
do
   python3 iteramrvac.py $AMRVAC_DIR Parfiles/Article/ test_grackle $i
   mpiexec 
done

