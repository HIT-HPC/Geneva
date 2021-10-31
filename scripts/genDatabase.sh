# !/bin/bash
# generate database for all datasets in folder ../datasets/ , some
# dataset may download as the .g file originally which don't need to
# convert to our format.

for file in `ls ../datasets/`
do
    if [ "${file##*.}"x = "g"x ];then
        datasetName=`echo ${file%.g}`
        echo -e "generate database for dataset ${datasetName} \n"
        ./genDatabase ${datasetName}
        echo -e "dataset ${datasetName} processed \n"
        sleep 1
    fi
done
