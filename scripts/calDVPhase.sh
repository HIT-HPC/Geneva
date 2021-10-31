# !bin/bash
# gpumine for all query

# datasets=(DD Enron FIRSTMM-DB Gowalla LiveJournal Orkut Patents REDDIT Wiki-Vote mico)

dataset=Patents

#for file in ${datasets[@]};

for file in `ls ../datasets/${dataset}/`
do
    if [ "${file##*.}"x = "query"x ];then
        queryName=`echo ${file%.query}`
        datasetName=`echo ${queryName%%_*}`
        echo "CalDVPhase for query ${queryName} on dataset ${datasetName}"
        ./calDVPhase ${datasetName} ./${dataset}/${queryName}
#        echo -e "CalDVPhase for query ${queryName} on dataset ${datasetName} processed.\n"
        sleep 1
    fi
done
