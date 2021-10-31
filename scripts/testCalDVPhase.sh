# !/bin/bash

for i in {4..5};
do
    bash ./gpumine.sh 2>&1 | tee gpuminesinaweibo${i}.log
    bash ./gpumineSingleV.sh 2>&1 | tee gpumineSVsinaweibo${i}.log
done

