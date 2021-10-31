# !/bin/bash
# generate query for all datasets which has generated .mygraph file,
# input query node number, edge number, query number for each graph

# datasets=`ls ../datasets/`
# datasets=(DD Enron FIRSTMM-DB Gowalla LiveJournal Orkut Patents REDDIT-MULTI-12K Wiki-Vote mico)
# datasets=(Patents.mygraph Wiki-Vote.mygraph LiveJournal.mygraph)
datasets=(Enron.mygraph)

# for file in `ls ./`
for file in ${datasets[@]};
do
    if [ "${file##*.}"x = "mygraph"x ];then
        graphName=`echo ${file%.mygraph}`
        for x in {4..12}; 
        do
            echo -e "generate query with node ${x}, edge ${x}, dataset ${file}\n"
            ./genQuery ${graphName} $x $x 1
            sleep 2
            y=`echo $((${x}+${x}/2))`
            echo -e "generate query with node ${x}, edge ${y}, dataset ${file}\n"
            ./genQuery ${graphName} $x $y 1
            sleep 2
        done
        echo -e "dataset ${graphName} generate query processed.\n"
        sleep 1
    fi
done
