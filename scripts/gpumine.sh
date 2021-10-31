# !bin/bash
# gpumine for all query

# datasets=(DD Enron FIRSTMM-DB Gowalla LiveJournal Orkut Patents REDDIT-MULTI-12K Wiki-Vote mico)

dataset=sinaweiboquery

#for file in ${datasets[@]};

for file in `ls ../datasets/${dataset}/`
do
    if [ "${file##*.}"x = "query"x ];then
        queryName=`echo ${file%.query}`
        datasetName=`echo ${queryName%%_*}`
        echo "GPUmine for query ${queryName} on dataset ${datasetName}"
        ./gpumine ${datasetName} ./${dataset}/${queryName}
        echo -e "GPUmine for query ${queryName} on dataset ${datasetName} processed.\n"
        sleep 1
    fi
done
