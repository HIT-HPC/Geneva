#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <bitset>
#include <list>
#include "common.h"
using namespace std;

class labelInfo{
public:
    string label;
    uint posInGlobal;
    uint instanceCount;
};
class subPartVid{
public:
    list<uint> subpart;
    uint vid;
};
class subPartInfo{
public:
    uint len;
    uint svid;
    uint evid;
    uint prelen;
    uint index;
};
uint32_t MurmurHash2(const void * key, int len, uint32_t seed)
{
    // 'm' and 'r' are mixing constants generated offline.
    // They're not really 'magic', they just happen to work well.
    const uint32_t m = 0x5bd1e995;
    const int r = 24;
    // Initialize the hash to a 'random' value
    uint32_t h = seed ^ len;
    // Mix 4 bytes at a time into the hash
    const unsigned char * data = (const unsigned char *) key;
    while (len >= 4) 
    {
        uint32_t k = *(uint32_t*) data;
        k *= m;
        k ^= k >> r;
        k *= m;
        h *= m;
        h ^= k;
        data += 4;
        len -= 4;
    }
    // Handle the last few bytes of the input array
    switch (len) 
    {
        case 3:
            h ^= data[2] << 16;
        case 2:
            h ^= data[1] << 8;
        case 1:
          h ^= data[0];
          h *= m;
    };
    // Do a few final mixes of the hash to ensure the last few
    // bytes are well-incorporated.
    h ^= h >> 13;
    h *= m;
    h ^= h >> 15;
    return h;
}
void extractEUINT(string &line, uint *e_info){
    int start = 0,len=0,i=0;
    int count=0;
    while(i<line.length()){
        if(line[i]>='0' && line[i]<='9'){
            start = i;
            len=0;
            while(i<line.length()){
                if(line[i]>='0' && line[i]<='9'){
                    len++;
                }else{
                    break;
                }
                i++;
            }
            string temp = line.substr(start,len);
            e_info[count] = stoul(temp);
            count++;
            if(count==3)
                break;
        }
        i++;
    }
}
void extractUINT(string &line, map<uint,uint>&v_label,uint &v_num){
    int start = 0,len=0,i=0,vid,lid;
    while(i<line.length()){
        if(line[i]>='0' && line[i]<='9'){
            start = i;
            while(i<line.length()){
                if(line[i]>='0' && line[i]<='9'){
                    len++;
                }else{
                    break;
                }
                i++;
            }
            break;
        }
        i++;
    }
    string temp = line.substr(start,len);
    vid = stoul(temp);
    len=0;
    while(i<line.length()){
        if(line[i]>='0' && line[i]<='9'){
            start = i;
            while(i<line.length()){
                if(line[i]>='0' && line[i]<='9'){
                    len++;
                }else{
                    break;
                }
                i++;
            }
            break;
        }
        i++;
    }
    temp = line.substr(start,len);
    lid = stoul(temp);
    if(vid+1>v_num){
        v_num=vid+1;
    }
    v_label.insert(pair<uint,uint>(vid,lid));
}
bool sortFuncReorder(map<uint,map<uint,vector<uint> > >& m1,map<uint,map<uint,vector<uint> > >& m2){
    return m1.size() > m2.size();
}
bool sortFuncVec(vector<uint> &m1, vector<uint> &m2){
    return m1.size()>m2.size();
}
bool sortFunc(map<uint,bitset<1> >& m1,map<uint,bitset<1> >& m2){
    return m1.size() > m2.size();
}
bool sortFuncUINT(uint m1, uint m2){
    return m1 < m2;
}

//basepart should be new created map, because we need to modify basepart
void intersect(map<uint,bitset<1> > &basePartition, map<uint,bitset<1> > &curPartition,
                vector<vector<uint> >&subPartitions){

    map<uint,bitset<1> >::iterator it;
    map<uint,bitset<1> >::iterator itt;
    vector<uint> commonV;
    //find common vertices in basePartition and curPartition
    for(it=basePartition.begin();it!=basePartition.end();it++){
        uint key = it->first;
        itt = curPartition.find(key);
        if(itt!=curPartition.end()){
            commonV.push_back(key);
        }
    }

    //find common vetices of curPartition and previous subPartitions
    uint totsize = subPartitions.size();
    for(uint i=0;i<totsize;++i){
        //divide the sub partition into new and old. new contains vertices that
        //belong to both basePartition and curPartition while old contains vertices
        //that only exist in base partition.
        vector<uint> &tmpPart = subPartitions[i];
        vector<uint> tmpNewPart;
        vector<uint> tmpOldPart;
        for(uint j=0;j<tmpPart.size();++j){
            itt = curPartition.find(tmpPart[j]);
            if(itt!=curPartition.end()){ tmpNewPart.push_back(tmpPart[j]); }
            else{ tmpOldPart.push_back(tmpPart[j]); }
        }
        //if all vertices in this subpartition are existed in curPartition
        //we do not need to generate new sub partition, add label partition id to this sub partition,
        if(tmpPart.size()==tmpNewPart.size() || tmpPart.size()==tmpOldPart.size()){ /*do nothing*/ }
        else if(tmpPart.size()>tmpNewPart.size()){subPartitions[i] = tmpOldPart; subPartitions.push_back(tmpNewPart); }
        else{ cout<<"wrong"<<endl; }
    }

    if(commonV.size()>0){
        subPartitions.push_back(commonV);
    }
    for(uint i=0;i<commonV.size();++i){
        uint key = commonV[i];
        basePartition.erase(key);
    }
}

//first, we need to evaluate whether there are
//common subpart ids between a and b.
//If this is true, we only keep common subpart ids.
//If not, we keep a as it is.
void overlapIn(list<uint> &a, map<uint,bool>&b){
    list<uint> c;
    list<uint>::iterator ite;
    for(ite=a.begin();ite!=a.end();++ite){
        uint base = (*ite);
        if(b.find(base)!=b.end()){
            c.push_back(base);
        }
    }
    if(c.size()>0){
        a.assign(c.begin(),c.end());
    }
}
void deleteSubpart(list<uint> &finalsubpart,list<uint> &tmp){
    list<uint>::iterator iter = finalsubpart.begin(),itera;
    while(iter!=finalsubpart.end()){
        uint base = (*iter);
        uint flag=0;
        for(itera=tmp.begin();itera!=tmp.end();++itera){
            if((*iter)==(*itera)){
                iter = finalsubpart.erase(iter);
                flag = 1;
                break;
            }
        }
        if(flag==0){
            ++iter;
        }
    }
}
void performReorder(vector<subPartVid>&reorderVec,vector<vector<uint> >&result,map<uint,uint>&reordermap){
    list<uint>::iterator iter;
    for(uint i=0;i<reorderVec.size();++i){
        list<uint> &tmp = reorderVec[i].subpart;
        uint vid = reorderVec[i].vid-1;
        for(iter=tmp.begin();iter!=tmp.end();++iter){
            vector<uint> &vertices = result[(*iter)];
            for(uint k=0;k<vertices.size();++k){
                reordermap.insert(pair<uint,uint>(vertices[k],vid+1));
                --vid;
            }
        }
    }
}

void performReorder(map<uint,bitset<1> >&basePartition,uint newVid,map<uint,uint> &reordermap){
    map<uint,bitset<1> >::iterator ite;
    uint startid = newVid-basePartition.size()+1;
    for(ite=basePartition.begin();ite!=basePartition.end();++ite){
        reordermap.insert(pair<uint,uint>(ite->first,startid));
        ++startid;
    }
}

//graph[i] is the ith label partition, each partition is a map
//the key is the id of source vertex, and the value is a vector that contains
//neighbors of source vertex that comply the label.
/*TODO: we need to verify that, in a edge label partition,
 * if the length of contiguous part a is greater than the length of contiguous part b,
 * then the biggest vertex id in a is smaller than the smallest vertex id in b.
 **/
void genReorderMap(vector<map<uint,map<uint,vector<uint> > > > &graph,
    map<uint,uint> &reordermap){

    //in reorder, we do not need neighbors, therefore we use bitset<1> instead
    vector<map<uint,bitset<1> > > graphCopy;
    cout<<"generating simplified data for label partitions"<<endl;
    for(uint i=0;i<graph.size();++i){
        map<uint,map<uint,vector<uint> > > &tmpa = graph[i];
        map<uint,map<uint,vector<uint> > >::iterator iteraa;
        map<uint,bitset<1> > tmpmap;
        for(iteraa=tmpa.begin();iteraa !=tmpa.end();iteraa++){
            bitset<1> a;
            tmpmap.insert(pair<uint,bitset<1> >(iteraa->first,a));
        }
        graphCopy.push_back(tmpmap);
    }
    cout<<"generate simplified data done"<<endl;
    //sort edge label partitions by the number of source vertices the partition has in decreasing order
    sort(graphCopy.begin(),graphCopy.end(),sortFunc);
    //each label labelPartition contains many vertices, our idea is
    //the more vertices in the labelPartition, the less contiguous subparts it contains
    //therefore, we sort labelPartitions based on the number of vertices it contains.
    //in each iteration, we reorder all vertices in this label partition
    uint totLabelPartition = graphCopy.size();
    uint newVid;
    uint reorderVNum = 0;
    for(uint ii=0;ii<totLabelPartition;++ii){
        cout<<"take label partition "<<ii<<" as the base label partition"<<endl;
        //basePartition is the labelPartition contains the most number of vertices
        //in each iteration, we reorder vertex ids for basePartition
        map<uint,bitset<1> > &basePartition = graphCopy[0];
        newVid = reorderVNum+basePartition.size();
        reorderVNum += basePartition.size();
        uint testBasePartSize = basePartition.size();
        //basePartition will be divided into sub-partitions that overlaps with other label partitions
        //each sub-partition contains contiguous vertices
        //the index of a sub-partition in subPartitions is the id of this sub-partition
        vector<vector<uint> > subPartitions;
        map<uint,map<uint,bool> > labelPartToSubPartId;
        //iterate over all labelPartitions, except the first one, to divide the base partition into sub partitions
        for(uint i=1;i<graphCopy.size();++i){
            map<uint,bitset<1> > &tmpPartition = graphCopy[i];
            //utilize intersect to divide the base labelPartition into subpartions, each subpartition overlaps with at least
            //one label partitions. Vertices in subpartitions are deleted from base partition.
            //Remaining vertices in base partition are not overlapped with other label partitions.
            //Subpartions are sotred in subPartitions, each entry in subPartitions represents a sub partition that contains
            //vertices can be assigned with contiguous ids. The entry index is the id of this sub partition.
            intersect(basePartition, tmpPartition, subPartitions);
        }
        //now construct labelpartidtosubid, remaining vertices in basepartition are added into subpartitions after this loop
        for(uint i=0;i<subPartitions.size();++i){
            vector<uint> &subpartv = subPartitions[i];
            uint tmpvertex = subpartv[0];
            bool found = false;
            for(uint j=1;j<graphCopy.size();++j){
                map<uint,bitset<1>> &tmppart = graphCopy[j];
                if(tmppart.find(tmpvertex)!=tmppart.end()){
                    found = true;
                    map<uint,map<uint,bool>>::iterator ite = labelPartToSubPartId.find(j);
                    if(ite!=labelPartToSubPartId.end()){
                        ite->second[i] = true;
                    }else{
                        map<uint,bool> tmpmap;
                        tmpmap[i] = true;
                        labelPartToSubPartId[j] = tmpmap;
                    }
                }
            }
            if(!found){
                cout << "do not find subpartition "<<i<<endl;
            }
            for(uint j=0;j<subpartv.size();++j){
                tmpvertex = subpartv[j];
                for(uint k=1;k<graphCopy.size();++k){
                    map<uint,bitset<1>> &tmppart = graphCopy[k];
                    if(tmppart.find(tmpvertex)!=tmppart.end()){
                        tmppart.erase(tmpvertex);
                    }
                }
            }
        }
        cout<<endl;


        //cout<<"generate intersection done. tmptestsize="<<tmptestsize<<" tot="<<tmptestsize+basePartition.size()<<" ori="<<testBasePartSize<<endl;
        cout<<"tot subpartitions="<<subPartitions.size()<<endl;
        vector<bool> reorderedSubParts(subPartitions.size(),false);
        map<uint,map<uint,bool> >::iterator iter;
        //after above loop, the base label partition has been divided into several sub partitions (each contains contiguous vertices)
        //and remaining unoverlapped vertices.
        //labelPartToSubPartId contains sub partitions' ids that label partition i contains
        //in the below loop, we reorder sub partitions based on the size of label partitions that contain it
        for(iter=labelPartToSubPartId.begin();iter != labelPartToSubPartId.end();iter++){
            vector<subPartVid> reorderVec;
            map<uint,map<uint,bool> >::iterator iter1;
            map<uint,bool> &basesubpart = iter->second;
            list<uint> finalsubpart;
            for(map<uint,bool>::iterator i=basesubpart.begin();i!=basesubpart.end();++i){
                uint subidtmp = i->first;
                if(!reorderedSubParts[subidtmp]){
                    finalsubpart.push_back(subidtmp);
                    reorderedSubParts[subidtmp] = true;
                }
            }
            uint tmptotvs=0;
            for(list<uint>::iterator i=finalsubpart.begin();i!=finalsubpart.end();++i){
                uint subid = (*i);
                tmptotvs += subPartitions[subid].size();
            }
            cout<<"reorder sub partition, there are "<<finalsubpart.size() <<" subpartitions, totv="<<tmptotvs<<" newvid="<<newVid<<endl;
            while(finalsubpart.size()>0){
                iter1 = iter;
                iter1++;
                list<uint> tmp(finalsubpart);
                for(;iter1!=labelPartToSubPartId.end();iter1++) {
                    map<uint,bool> &compsubpart = iter1->second;
                    overlapIn(tmp, compsubpart);
                    if (tmp.size() == 1) { break; }
                }
                subPartVid tmpsub;
                tmpsub.subpart = tmp;
                tmpsub.vid = newVid;
                reorderVec.push_back(tmpsub);
                uint totvertices=0;
                for(list<uint>::iterator i=tmp.begin();i!=tmp.end();++i){
                    totvertices += subPartitions[(*i)].size();
                }
                if(newVid<totvertices){
                    cout << "hahah, i got it. newVid="<<newVid<<" tot="<<totvertices<<endl;
                }
                newVid -= totvertices;
                deleteSubpart(finalsubpart,tmp);
            }//while(finalsubpart.size()>0)
            uint recordedparts=0;
            tmptotvs = 0;
            for(uint i=0;i<reorderVec.size();++i){
                list<uint> &tmplist = reorderVec[i].subpart;
                recordedparts += tmplist.size();
                for(list<uint>::iterator itc=tmplist.begin();itc!=tmplist.end();++itc){
                    uint subid = (*itc);
                    tmptotvs += subPartitions[subid].size();
                }
            }
            cout<<"recorded partitions="<<recordedparts<<" totvs="<<tmptotvs<<" newvid="<<newVid<<endl;
            performReorder(reorderVec,subPartitions,reordermap);
        }
        cout<<"reorde unoverlapped sub partition, totvs="<<basePartition.size()<<" newvid="<<newVid<<endl;
        performReorder(basePartition,newVid,reordermap);
        //reorder vertices for remaining vertices in base label partition
        graphCopy.erase(graphCopy.begin());
        //when performing intersect, we delete common vertices in other label partitions,
        //therefore, we need resort partitions in each iteration
        sort(graphCopy.begin(),graphCopy.end(),sortFunc);
    }
}

void genNeigIndex(vector<uint> &neig){
    vector<uint> newneig;
    if(neig.size()<VBLOCKSIZE){
        return;
    }else if(neig.size()<=32*VBLOCKSIZE){
        uint blocknum = neig.size()/VBLOCKSIZE;
        for(uint i=1;i<=blocknum;i++){
            newneig.push_back(neig[i*VBLOCKSIZE-1]);
        }
        if(neig.size()%VBLOCKSIZE>0){
            newneig.push_back(neig[neig.size()-1]);
        }
        neig.insert(neig.begin(),newneig.begin(),newneig.end());
    }else{
        uint blocknum = neig.size()/32;
        for(uint i=1;i<=31;++i){
            newneig.push_back(neig[i*blocknum-1]);
        }
        newneig.push_back(neig[neig.size()-1]);
        neig.insert(neig.begin(),newneig.begin(),newneig.end());
    }
}
//data format for our interval index format
//uint0 uint1 uint2 uint3 uint4 array1 array2 array3 array4
//uint0: the number of uints in this label partition.
//uint1: ther number of source vertices (including blank vertices) in this partition
//uint2: position of interval vertices neighbor addr (array2) and the lenght of this interval and previous intervals,
//uint3: position of hash table index.
//uint4: position of neighbor data
//array1: boundries of interval index, at most has 256 entries,
//  each entry contains 2 uints, including start vid, prelen(the length of this interval and previous intervals),
//  loaded into shared memory.
//array2: interval vertices neighbor addr and one extra uint at the end of this section to record the end position of the last one
//array3: hash index,
//array4: neighbors
//subpartsinfo stores all contiguous subparts of this label partition
void writeHybridDataStructure(vector<subPartInfo>&subpartsinfo, map<uint,map<uint,vector<uint> > >& newLabelPartition,
            ofstream &graph, ofstream &readableGraph,uint totNeigNum,labelInfo &label, uint &totWriteNum,map<uint,uint>&reorderVLabel){
    uint intervalNum = subpartsinfo.size();
    if(intervalNum>256){
	std::cout<<"intervalNum > 256, this label partition can not be used"<<std::endl;
        //sort(subpartsinfo.begin(),subpartsinfo.end(),[ ](subPartInfo &a, subPartInfo &b) { return a.len>b.len; });
        //sort(subpartsinfo.begin(),subpartsinfo.begin()+256,[ ](subPartInfo &a, subPartInfo &b) { return a.svid < b.svid; });
        //sort(subpartsinfo.begin(),subpartsinfo.end(),[ ](subPartInfo &a, subPartInfo &b) { return a.svid<b.svid; });
        intervalNum = 256;
    }
    uint totIntervalIndexNum=0,totHashIndexNum=0;
    for(uint i=0;i<intervalNum;++i)                   { totIntervalIndexNum += subpartsinfo[i].len; }
    //for(uint i=intervalNum;i<subpartsinfo.size();++i) { totHashIndexNum     += subpartsinfo[i].len; }
    uint totStoreNum = 5+intervalNum*2+1+1+totIntervalIndexNum+1+totHashIndexNum*32+totNeigNum;
    uint *writePos = (uint*)calloc(totStoreNum,sizeof(uint));
    writePos[0] = totStoreNum;
    writePos[1] = totIntervalIndexNum+totHashIndexNum*32/2;
    writePos[2] = intervalNum;
    writePos[3] = intervalNum*2+1+1+totIntervalIndexNum+1;
    writePos[4] = writePos[3]+totHashIndexNum*32;
    readableGraph<<"five indexs: "<<writePos[0]<<" "<<writePos[1]<<" "<<writePos[2]<<" "<<writePos[3]<<" "<<writePos[4]<<endl;
    uint intervalAddr = 5;
    uint intervalIndexAddr = 5+intervalNum*2+1+1;
    uint neigAddr = intervalIndexAddr+totIntervalIndexNum+1+totHashIndexNum*32;
    uint baseNeigAddr = neigAddr;
    uint preLen = 0;
    writePos[intervalIndexAddr+intervalNum] = 0;
    for(uint i=0;i<intervalNum;++i){
        writePos[intervalAddr+i] = subpartsinfo[i].svid;
        writePos[intervalAddr+i+intervalNum+1] = preLen+subpartsinfo[i].len;
        preLen = preLen+subpartsinfo[i].len;
        for(uint j=subpartsinfo[i].svid;j<=subpartsinfo[i].evid;++j){
            map<uint,vector<uint> > &neig = newLabelPartition[j];
            writePos[intervalIndexAddr] = neigAddr-baseNeigAddr;
            map<uint,vector<uint> >::iterator ite;
            uint k=0;
            for(ite=neig.begin();ite!=neig.end();++ite){
                vector<uint> &vlabelneig = ite->second;
                for(uint l=0;l<vlabelneig.size();++l){
                    writePos[neigAddr+k] = vlabelneig[l];
                    ++k;
                }
            }
            neigAddr += k;
            intervalIndexAddr++;
        }
    }
	std::cout<<"write hash index"<<std::endl;
    writePos[intervalIndexAddr] = neigAddr-baseNeigAddr;
    vector<uint> v_count_in_group(totHashIndexNum,0);
    uint graph_group_num = 32*totHashIndexNum;
    uint hashAddr = writePos[3];
    /*for(uint i=intervalNum;i<subpartsinfo.size();++i){
        uint startvid = subpartsinfo[i].svid;
        uint endvid = subpartsinfo[i].evid;
        for(uint j=startvid;j<=endvid;++j){
            uint key = j;
            uint group_index = MurmurHash2(&key,4,HASHSEED)%totHashIndexNum;
            if(v_count_in_group[group_index]==16){
                cout<<"we need another group"<<endl;
            }else{
                uint stored = v_count_in_group[group_index];
                writePos[hashAddr+group_index*32+stored*2+0] = key;
                stored++;
                if(stored==15){
                    writePos[hashAddr+group_index*32+stored*2+0] =  0xffffffff;
                    stored++;
                }
                v_count_in_group[group_index] = stored;
            }
        }
    }*/
    for(uint j=0;j<graph_group_num;j=j+2){
        if(writePos[hashAddr+j]==0){
        }else if(writePos[hashAddr+j]==0xffffffff){
            writePos[hashAddr+j+1] = neigAddr;
        }else{
            uint key = writePos[hashAddr+j];
            map<uint,vector<uint> > &neig = newLabelPartition[key];
            writePos[hashAddr+j+1] = neigAddr;
            map<uint,vector<uint> >::iterator ite;
            uint k=0;
            for(ite=neig.begin();ite!=neig.end();++ite){
                vector<uint> &vlabelneig=ite->second;
                for(uint l=0;l<neig.size();++l){
                    writePos[neigAddr+k] = vlabelneig[l];
                    ++k;
                }
            }
            neigAddr += k;
        }
    }
    graph.write((char*)writePos,totStoreNum*sizeof(uint));
    readableGraph<<"interval index"<<endl<<"svid: ";
    for(uint i=0;i<intervalNum;++i){
        readableGraph << writePos[intervalAddr+i]<<" ";
    }
    readableGraph<<"interval prelen: ";
    for(uint i=0;i<intervalNum+1+1;++i){
        readableGraph<<writePos[intervalAddr+intervalNum+i]<<" ";
    }
    readableGraph<<endl;
    readableGraph<<"interval addr"<<endl;
    intervalIndexAddr = 5+intervalNum*2+1+1;
    for(uint i=0;i<totIntervalIndexNum+1;++i){
        readableGraph<<writePos[intervalIndexAddr+i]<<" ";
    }
    readableGraph<<endl;
    readableGraph<<"neigbors"<<endl;
    for(uint i=0;i<intervalNum;++i){
        for(uint j=subpartsinfo[i].svid;j<=subpartsinfo[i].evid;++j){
            readableGraph<<"label:"<<reorderVLabel[j]<<" svid:"<<j<<" ";
            map<uint,vector<uint> > &neig = newLabelPartition[j];
            map<uint,vector<uint> >::iterator ite;
            for(ite=neig.begin();ite!=neig.end();++ite){
                uint evidlabel = ite->first;
                readableGraph<<"evidlabel:"<<evidlabel<<" ";
                readableGraph<<"neigs:{";
                vector<uint> &vlabelneig = ite->second;
                for(uint l=0;l<vlabelneig.size();++l){
                    readableGraph<<vlabelneig[l]<<" ";
                }
                readableGraph<<"} ";
            }
            readableGraph<<endl;
        }
    }
    readableGraph<<endl;
    free(writePos);
}

void addEdgeToGraphDebug(uint *e_info, map<uint,uint>&v_label){

    uint e_label = e_info[2];
    uint svid = e_info[0];
    uint evid = e_info[1];
    if(svid==evid){
        cout<<"wrong, equal vids"<<endl;
    }
    auto ite=v_label.find(svid);
    if(ite!=v_label.end()){
        v_label.erase(svid);
    }
    auto ite1=v_label.find(evid);
    if(ite1!=v_label.end()){
        v_label.erase(evid);
    }
}

void addEdgeToGraph(uint *e_info, vector<map<uint,map<uint,vector<uint> > > > &graph, 
    map<uint,uint>&v_label){

    uint e_label = e_info[2];
    uint svid = e_info[0];
    uint evid = e_info[1];
    uint evlabel = v_label[evid];
    map<uint,map<uint,vector<uint> > > &edgeLabelPart = graph[e_label];
    map<uint,map<uint,vector<uint> > >::iterator ite = edgeLabelPart.find(svid);
    if(ite==edgeLabelPart.end()){
        vector<uint> tmp;
        tmp.push_back(evid);
        map<uint, vector<uint> > tmpmap;
        tmpmap.insert(pair<uint,vector<uint> >(evlabel,tmp));
        edgeLabelPart.insert(pair<uint,map<uint,vector<uint> > >(svid,tmpmap));
    }else{
        map<uint,vector<uint> >&tmpmap = ite->second;
        map<uint,vector<uint> >::iterator itea = tmpmap.find(evlabel);
        if(itea==tmpmap.end()){
            vector<uint> tmp;
            tmp.push_back(evid);
            tmpmap.insert(pair<uint,vector<uint> >(evlabel,tmp));
        }else{
            vector<uint> &tmp = itea->second;
            bool found = false;
            for(uint i=0;i<tmp.size();++i){
                if(tmp[i]==evid){
                    found = true;
                    break;
                }
            }
            if(!found){
                tmp.push_back(evid);
            }
        }
    }
    svid = e_info[1];
    evid = e_info[0];
    evlabel = v_label[evid];
    ite = edgeLabelPart.find(svid);
    if(ite==edgeLabelPart.end()){
        vector<uint> tmp;
        tmp.push_back(evid);
        map<uint, vector<uint> > tmpmap;
        tmpmap.insert(pair<uint,vector<uint> >(evlabel,tmp));
        edgeLabelPart.insert(pair<uint,map<uint,vector<uint> > >(svid,tmpmap));
    }else{
        map<uint,vector<uint> >&tmpmap = ite->second;
        map<uint,vector<uint> >::iterator itea = tmpmap.find(evlabel);
        if(itea==tmpmap.end()){
            vector<uint> tmp;
            tmp.push_back(evid);
            tmpmap.insert(pair<uint,vector<uint> >(evlabel,tmp));
        }else{
            vector<uint> &tmp = itea->second;
            bool found = false;
            for(uint i=0;i<tmp.size();++i){
                if(tmp[i]==evid){
                    found = true;
                    break;
                }
            }
            if(!found){
                tmp.push_back(evid);
            }
        }
    }
}

uint genSubParts(vector<subPartInfo> &subpartsinfo,map<uint,map<uint,vector<uint> > > &labelPartition){
    if(labelPartition.size()==0){
        return 0;
    }
    map<uint,map<uint,vector<uint> > >::iterator itera=labelPartition.begin();
    uint totv = labelPartition.size();
    uint curvid = itera->first;
    uint previd = itera->first;
    uint start = previd;
    uint tmpPreLen=0;
    uint totNeigNum=0;
    map<uint,vector<uint> >::iterator ite;
    map<uint,vector<uint> >&vneig = itera->second;
    for(ite=vneig.begin();ite!=vneig.end();++ite){
        totNeigNum += ite->second.size();
    }
    itera++;
    while(itera != labelPartition.end()){
        curvid = itera->first;
        map<uint,vector<uint> >&vneig = itera->second;
        for(ite=vneig.begin();ite!=vneig.end();++ite){
            totNeigNum += ite->second.size();
        }
        if(curvid==previd+1){
            previd = curvid;
        }else if(curvid>previd+1){
            subPartInfo tmpinfo;
            tmpinfo.svid = start;
            tmpinfo.evid = previd;
            tmpinfo.prelen = tmpPreLen;
            tmpinfo.len = previd-start+1;
            //tmp.index = sortsubparts.size();
            subpartsinfo.push_back(tmpinfo);
            tmpPreLen= tmpPreLen+previd-start+1;
            start=curvid;
            previd=curvid;
        }else{
            cout << "wrong in "<< __FUNCTION__ << " line: "<< __LINE__ <<endl;
        }
        itera++;
    }
    subPartInfo tmpinfo;
    tmpinfo.svid = start;
    tmpinfo.evid = previd;
    tmpinfo.prelen = tmpPreLen;
    tmpinfo.len = previd-start+1;
    subpartsinfo.push_back(tmpinfo);
    return totNeigNum;
}

void genReorderGraph(vector<map<uint,map<uint,vector<uint> > > > &graph,vector<map<uint,map<uint,vector<uint> > > > &reorderGraph,
                        map<uint,uint>&reordermap){
    for(uint i=0;i<graph.size();++i){
        map<uint,map<uint,vector<uint> > >::iterator itera;
        map<uint,map<uint,vector<uint> > > &tmpPartition = graph[i];
        map<uint,map<uint,vector<uint> > > newmap;
        for(itera=tmpPartition.begin();itera != tmpPartition.end();itera++){
            uint orisvid = itera->first;
            uint newsvid = reordermap[orisvid];
            map<uint,vector<uint> > &orimap = itera->second;
            map<uint,vector<uint> >::iterator ite;
            map<uint,vector<uint> > tmpmap;
            for(ite=orimap.begin();ite!=orimap.end();ite++) {
                uint orievidlabel = ite->first;
                uint newevidlabel = orievidlabel;
                vector<uint> &orivec = ite->second;
                vector<uint> newvec;
                for (uint j = 0; j < orivec.size(); ++j) {
                    uint newvid = reordermap[orivec[j]];
                    newvec.push_back(newvid);
                }
                sort(newvec.begin(), newvec.end(), sortFuncUINT);
                if (newvec.size() > VBLOCKSIZE) {
                    genNeigIndex(newvec);
                }
                tmpmap.insert(pair<uint, vector<uint>>(newevidlabel, newvec));
            }
            newmap.insert(pair<uint,map<uint,vector<uint> > >(newsvid,tmpmap));
        }
        reorderGraph.push_back(newmap);
    }
}

//for now, we consider labeld graph
//and then stores only one copy of this part
int main(int argc, char* argv[]){
    uint v_num=0, e_num=0;
    string baseFileName = string(argv[1]);
    //string path = "/home/lgz/projects/clion/";
    string path = "../datasets/";
    ifstream input(path+baseFileName+".g",std::ios::in);
    ofstream output(path+baseFileName+".mygraph",std::ios::binary|std::ios::out);
    if(!input.is_open()){
        cout<<"open "<<baseFileName+".g"<<" is wrong"<<endl;
        return 1;
    }
    if(!output.is_open()){
        cout<<"open "<<baseFileName+".mygraph"<<" is wrong"<<endl;
        return 1;
    }
    map<uint,uint> v_label;
    string line;
    while(getline(input,line)){
        if(line.length()==0){ continue; }
        if(line[0]=='t'){ continue; }
        if(line[0]=='v'){ 
            extractUINT(line,v_label,v_num);
        }
        else if(line[0]=='e'){ break; }
    }
    cout<<"reading vertices="<<v_num<<" "<<v_label.size()<<endl;
    //v_num = v_label.size();
    //map<v-e-v label, map<source vertex, neighbors>>
    vector<map<uint,map<uint,vector<uint> > > > graph;
    map<string,uint> labeltoid;
    vector<labelInfo> labels;
    //map<string,uint> idtolabel;
    uint label_count=0;
    while(true){
        if(line[0]=='e'){
            uint e_info[3];
            extractEUINT(line,e_info);
            if(e_info[2]+1>graph.size()){
                graph.resize(e_info[2]+1);
            }
            addEdgeToGraph(e_info, graph, v_label);
            //addEdgeToGraphDebug(e_info, v_label);
        }
        if(!getline(input,line)){
            break;
        }
    }
    input.close();
    /*cout<<"remaining vertices num="<<v_label.size()<<endl;
    for(auto ite=v_label.begin();ite!=v_label.end();++ite){
        cout<<" "<<ite->first;
    }
    cout << endl;
    return 0;*/
    cout<<"reading data done: "<< graph.size()<<endl;
    //first old vid, second new vid
    map<uint,uint> reordermap;
    //construct reorder map
    genReorderMap(graph,reordermap);
    map<uint,uint> reorderVLabel;
    for(auto ite=reordermap.begin();ite!=reordermap.end();++ite){
        uint oldvid = ite->first;
        uint newvid = ite->second;
        uint vidlabel = v_label[oldvid];
        reorderVLabel[newvid] = vidlabel;
    }
    cout<<"generate reorder map done: "<<reordermap.size()<<endl;
    //use reorder map to build new label graph
    vector<map<uint,map<uint,vector<uint> > > > reorderGraph;
    genReorderGraph(graph,reorderGraph,reordermap);
    cout<<"generate reorder graph done"<<endl;

    ofstream reordergraphfile(baseFileName+".reorder",std::ios::out);
    reordergraphfile<<reorderGraph.size()<<std::endl;
    for(uint i=0;i<reorderGraph.size();++i){
        map<uint,map<uint,vector<uint>>> &tmplabelpart = reorderGraph[i];
        reordergraphfile<<tmplabelpart.size()<<std::endl;
        for(auto ite = tmplabelpart.begin();ite!=tmplabelpart.end();++ite){
            uint svid = ite->first;
            map<uint,vector<uint>> &neigs = ite->second;
            reordergraphfile<<svid<<" "<<neigs.size()<<std::endl;
            for(auto ite1 = neigs.begin();ite1!=neigs.end();++ite1){
                uint evidlabel = ite1->first;
                vector<uint> &labelneigs = ite1->second;
                reordergraphfile<<evidlabel<<" "<<labelneigs.size()<<std::endl;
                for(uint j=0;j<labelneigs.size();++j){
                    reordergraphfile<<labelneigs[j]<<" ";
                }
                reordergraphfile<<std::endl;
            }
        }
    }
    reordergraphfile.close();
	std::cout<<"write readableGraph"<<std::endl;
    ofstream readableGraph(baseFileName+"_readableGraph.txt",std::ios::out);
    //output readable graph

    uint totWriteNum = 0;
    //construct index section for each graph and write all data to binary file
    //ofstream test("genraw.txt",std::ios::out);
    //test<<reorderGraph.size()<<endl;
    //sort(reorderGraph.begin(),reorderGraph.end(),sortFuncReorder);
    uint totpart = reorderGraph.size()-1;
	std::cout<<"totpart="<<totpart<<std::endl;
    output.write((char*)(&totpart),4);
    for(uint i=0;i<reorderGraph.size();++i){
        vector<subPartInfo> subpartsinfo;
        uint totNeigNum = genSubParts(subpartsinfo,reorderGraph[i]);
        if(subpartsinfo.size()==0){
            continue;
        }
        readableGraph<<"edge label: "<< i <<endl;
	std::cout<<"write edge label: "<< i <<std::endl;
        //write data to binary file
        writeHybridDataStructure(subpartsinfo,reorderGraph[i], output,readableGraph,totNeigNum,labels[i],totWriteNum,reorderVLabel);
    }
	readableGraph.close();
	std::cout<<"write metadata file"<<std::endl;
    //test.close();
    ofstream meta(path+baseFileName+".metadata",ios::out);
    map<uint,uint>::iterator iter;
    for(iter=reordermap.begin();iter!=reordermap.end();iter++){
        //new id, label, old id
        meta<<iter->second<<" "<<v_label[iter->first] <<" "<< iter->first<<endl;
    }
    meta.close();
	std::cout<<"write metadata file done"<<std::endl;
    input.close();
    output.close();
    return 0;
}
