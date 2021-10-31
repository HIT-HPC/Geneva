#include <iostream>
#include <fstream>
#include <algorithm>
#include <hash_map>
#include "common.h"
#include "genMatchOrderCPU.h"
#include "loadDataBase.h"

using namespace __gnu_cxx;
class recordPosInfo{
public:
    uint size;
    uint* addr;
};

std::map<uint,uint> posOfVisitedV;
std::map<uint64_t,recordPosInfo> recordPosData;
std::map<uint,std::vector<uint>> equalPos;
std::vector<Phase> matchOrder;
std::map<uint64_t,uint> repeatFlag;
std::map<uint,std::map<uint,bool>> edgeLabelStat;
std::map<uint,uint> restricts;
std::map<uint,uint*> allEdgeLabelPartitions;
std::map<uint,uint> queryLabel;
std::vector<hash_map<uint,hash_map<uint,std::vector<uint> > > > reorderGraph;
uint *vLabel;

void constructRecordPos(){
    for(auto ite=edgeLabelStat.begin();ite!=edgeLabelStat.end();++ite){
        uint edgelabel = ite->first;
        uint *tmplabelpart = allEdgeLabelPartitions[edgelabel];
        uint vertexNum = tmplabelpart[1];
        uint typeNum = ite->second.size();
        uint *base = new uint[vertexNum*typeNum*2]();
        uint64_t key = edgelabel;
        key = (key<<32)|0xffffffff;
        recordPosInfo tmp;
        tmp.size = vertexNum*typeNum*2;
        tmp.addr = base;
        recordPosData.insert(std::pair<uint64_t,recordPosInfo>(key,tmp));
        for(auto itea=ite->second.begin();itea!=ite->second.end();++itea){
            uint vlabel = itea->first;
            uint64_t keya = edgelabel;
            keya=(keya<<32)|vlabel;
            recordPosInfo tmpa;
            tmpa.size = vertexNum*2;
            tmpa.addr = base;
            recordPosData.insert(std::pair<uint64_t,recordPosInfo>(keya,tmpa));
            base = base + vertexNum*2;
        }
    }
}

uint genInit(std::ofstream &outfile, uint *partialRes, uint uintsNum){

    Phase &firstPhase = matchOrder[0];
    uint edgeLabel = firstPhase.edgeLabel;
    std::vector<uint> &pairs = firstPhase.pairs;
    uint totRowNum=0;
    if(pairs.size()==2){
        uint maxRowNum = uintsNum/2;
        uint svid = pairs[0], evid = pairs[1];
        posOfVisitedV.insert(std::pair<uint,uint>(svid,0));
        posOfVisitedV.insert(std::pair<uint,uint>(evid,1));
        uint svidlabel = queryLabel[svid], evidlabel = queryLabel[evid];
        uint64_t key = edgeLabel;
        key = (key<<32)|evidlabel;
        bool isRecord = (repeatFlag[key]&0x0000ffff)>1;
        bool isExtRestrict = restricts[evid]==svid;
        if(evidlabel==svidlabel){
            std::vector<uint> tmpv;
            tmpv.push_back(0);
            tmpv.push_back(1);
            equalPos.insert(std::pair<uint,std::vector<uint>>(svidlabel,tmpv));
        }else{
            std::vector<uint> tmpv1;
            std::vector<uint> tmpv2;
            tmpv1.push_back(0);
            tmpv2.push_back(1);
            equalPos.insert(std::pair<uint,std::vector<uint>>(svidlabel,tmpv1));
            equalPos.insert(std::pair<uint,std::vector<uint>>(evidlabel,tmpv2));
        }
        uint rowIndex=0;
        hash_map<uint,hash_map<uint,std::vector<uint>>> &edgelabepart = reorderGraph[edgeLabel];
        for(auto ite=edgelabepart.begin();ite!=edgelabepart.end();++ite){
            uint tmpsvid = ite->first;
            if(vLabel[tmpsvid]==svidlabel){
                hash_map<uint,std::vector<uint>> &labelneigs = ite->second;
                auto ite1 = labelneigs.find(evidlabel);
                if(ite1!=labelneigs.end()){
                    std::vector<uint> &neigs = ite1->second;
                    for(uint i=0;i<neigs.size();++i){
                        std::vector<uint> tmpresult;
                        if(isExtRestrict){
                            if(tmpsvid>neigs[i]){
                                partialRes[rowIndex*2+0] = tmpsvid;
                                partialRes[rowIndex*2+1] = neigs[i];
                                ++rowIndex;
                                ++totRowNum;
                                if(rowIndex==maxRowNum){
                                    outfile.write((char*)partialRes,maxRowNum*2*sizeof(uint));
                                    rowIndex=0;
                                }
                            }
                        }else{
                            partialRes[rowIndex*2+0] = tmpsvid;
                            partialRes[rowIndex*2+1] = neigs[i];
                            ++rowIndex;
                            ++totRowNum;
                            if(rowIndex==maxRowNum){
                                outfile.write((char*)partialRes,maxRowNum*2*sizeof(uint));
                                rowIndex=0;
                            }
                        }
                    }
                }
            }
        }
        outfile.write((char*)partialRes,rowIndex*2*sizeof(uint));
    }else{
        uint maxRowNum = uintsNum/3;
        uint svid1 = pairs[0],evid1 = pairs[1],svid2 = svid1,evid2 = pairs[3];
        uint svid1label = queryLabel[svid1], evid1label = queryLabel[evid1], svid2label = svid1label, evid2label = queryLabel[evid2];
        uint64_t key1 = edgeLabel, key2=edgeLabel;
        key1 = (key1<<32)|evid1label; key2 = (key2<<32)|evid2label;
        posOfVisitedV.insert(std::pair<uint,uint>(svid1,0)); posOfVisitedV.insert(std::pair<uint,uint>(evid1,1)); posOfVisitedV.insert(std::pair<uint,uint>(evid2,2));
        bool isExt1Restrict = restricts[evid1]!=0xffffffff;
        bool isExt1and2Restrict = restricts[evid2]!=0xffffffff?restricts[evid2]==evid1:false;
        bool isExt2Restrict = isExt1and2Restrict?false:restricts[evid2]!=0xffffffff;
        bool isExt1and2SameLabel = evid1label==evid2label;
        bool isRecord1 = (repeatFlag[key1]&0x0000ffff)>1;
        bool isRecord2 = (repeatFlag[key2]&0x0000ffff)>1;
        if(svid1label==evid1label && svid1label==evid2label){
            std::vector<uint> tmpv;
            tmpv.push_back(0); tmpv.push_back(1); tmpv.push_back(2);
            equalPos.insert(std::pair<uint,std::vector<uint>>(svid1label,tmpv));
        }else if(svid1label==evid1label){
            std::vector<uint> tmpv1; std::vector<uint> tmpv2;
            tmpv1.push_back(0); tmpv1.push_back(1); tmpv2.push_back(2);
            equalPos.insert(std::pair<uint,std::vector<uint>>(svid1label,tmpv1)); equalPos.insert(std::pair<uint,std::vector<uint>>(evid2label,tmpv2));
        }else if(svid1label==evid2label){
            std::vector<uint> tmpv1; std::vector<uint> tmpv2;
            tmpv1.push_back(0); tmpv1.push_back(2); tmpv2.push_back(1);
            equalPos.insert(std::pair<uint,std::vector<uint>>(svid1label,tmpv1)); equalPos.insert(std::pair<uint,std::vector<uint>>(evid1label,tmpv2));
        }else if(evid1label==evid2label){
            std::vector<uint> tmpv1; std::vector<uint> tmpv2;
            tmpv1.push_back(0); tmpv2.push_back(1); tmpv2.push_back(2);
            equalPos.insert(std::pair<uint,std::vector<uint>>(svid1label,tmpv1)); equalPos.insert(std::pair<uint,std::vector<uint>>(evid1label,tmpv2));
        }else{
            std::vector<uint> tmpv1; std::vector<uint> tmpv2; std::vector<uint> tmpv3;
            tmpv1.push_back(0); tmpv2.push_back(1); tmpv3.push_back(2);
            equalPos.insert(std::pair<uint,std::vector<uint>>(svid1label,tmpv1)); equalPos.insert(std::pair<uint,std::vector<uint>>(evid1label,tmpv2)); equalPos.insert(std::pair<uint,std::vector<uint>>(evid2label,tmpv3));
        }
        uint rowIndex=0;
        hash_map<uint,hash_map<uint,std::vector<uint>>> &edgelabepart = reorderGraph[edgeLabel];
        for(auto ite=edgelabepart.begin();ite!=edgelabepart.end();++ite){
            uint tmpsvid = ite->first;
            if(vLabel[tmpsvid]==svid1label){
                hash_map<uint,std::vector<uint>> &labelneigs = ite->second;
                auto ite1 = labelneigs.find(evid1label);
                auto ite2 = labelneigs.find(evid2label);
                if(ite1!=labelneigs.end() && ite2!=labelneigs.end()){
                    std::vector<uint> &neigs1 = ite1->second;
                    std::vector<uint> &neigs2 = ite2->second;
                    for(uint i=0;i<neigs1.size();++i){
                        for(uint k=0;k<neigs2.size();++k){
                            if(neigs2[k]!=neigs1[i]){
                                if(isExt1Restrict && tmpsvid<neigs1[i]){
                                    continue;
                                }
                                if(isExt2Restrict && tmpsvid<neigs2[k]){
                                    continue;
                                }
                                if(isExt1and2Restrict && neigs1[i]<neigs2[k]){
                                    continue;
                                }
                                partialRes[rowIndex*3+0] = tmpsvid;
                                partialRes[rowIndex*3+1] = neigs1[i];
                                partialRes[rowIndex*3+2] = neigs2[k];
                                ++rowIndex;
                                ++totRowNum;
                                if(rowIndex==maxRowNum){
					std::cout<<"excedd in init"<<std::endl;
                                    outfile.write((char*)partialRes,maxRowNum*3*sizeof(uint));
                                    rowIndex=0;
                                }
                            }
                        }
                    }
                }
            }
        }
        outfile.write((char*)partialRes,rowIndex*3*sizeof(uint));
    }
    return totRowNum;
}

uint genExtEmb(Phase &curPhase, bool isLastPhase, std::string inputfilename,std::string outputfilename,
    uint *partialRes,uint uintsNum, uint partialRowNum){

    std::ifstream inputfile(inputfilename,std::ios::in|std::ios::binary);
    std::ofstream outputfile(outputfilename,std::ios::out|std::ios::binary);
    uint edgeLabel = curPhase.edgeLabel;
    std::vector<uint> &pairs = curPhase.pairs;
    uint totRowNum=0;
    uint preEmbLen = posOfVisitedV.size();
    if(pairs.size()==2){
        uint maxRowNum = uintsNum/(preEmbLen+1);
        uint svid = pairs[0], evid = pairs[1];
        uint index = posOfVisitedV.size();
        posOfVisitedV.insert(std::pair<uint,uint>(evid,index));
        uint svidlabel = queryLabel[svid], evidlabel = queryLabel[evid];
        uint64_t key = edgeLabel;
        key = (key<<32)|evidlabel;
        bool isRecord = (repeatFlag[key]&0x0000ffff)>1;
        bool isExtRestrict = restricts[evid]!=0xffffffff;
        std::cout<<"restrict:"<<evid<<"<"<<restricts[evid]<<" full 1:"<<0xffffffff<<" isExtRestrict="<<isExtRestrict<<std::endl;
        bool useRecord = (repeatFlag[key]&0x0000ffff)<((repeatFlag[key]>>16)&0x0000ffff);

        std::vector<uint> auxArray(8);
        auxArray[0]=evidlabel;auxArray[1]=evidlabel;
        auxArray[4]=posOfVisitedV[svid];auxArray[5]=posOfVisitedV[svid];
        auxArray[2]=isExtRestrict?posOfVisitedV[restricts[evid]]:0;
        auxArray[3]=auxArray[2];
        if(equalPos.find(evidlabel)!=equalPos.end()){
            std::vector<uint> &tmpv = equalPos[evidlabel];
            auxArray[6]=tmpv.size();
            auxArray[7]=tmpv.size();
            auxArray.insert(auxArray.end(),tmpv.begin(),tmpv.end());
            tmpv.push_back(posOfVisitedV[evid]);
        }else{
            auxArray[5]=0;
            auxArray[6]=0;
            std::vector<uint> tmpv;
            tmpv.push_back(posOfVisitedV[evid]);
            equalPos.insert(std::pair<uint,std::vector<uint>>(evidlabel,tmpv));
        }
        uint lesspos = isExtRestrict?posOfVisitedV[restricts[evid]]:0;
        hash_map<uint,hash_map<uint,std::vector<uint>>> &edgelabepart = reorderGraph[edgeLabel];
        uint rowIndex=0, embLen=posOfVisitedV.size();
        uint *emb = (uint *)malloc(32*sizeof(uint));
        for(uint i=0;i<partialRowNum;++i){
            if(i%20000==0){
                std::cout<<"extension, tot="<<partialRowNum<<" iter="<<i<<std::endl;
            }
            inputfile.read((char*)emb,preEmbLen*sizeof(uint));
            uint svidpos = posOfVisitedV[svid];
            uint tmpsvid = emb[svidpos];
            auto ite = edgelabepart.find(tmpsvid);
            if(ite==edgelabepart.end()) { continue; }
            hash_map<uint,std::vector<uint>> &labelneigs = ite->second;
            auto ite1 = labelneigs.find(evidlabel);
            if(ite1==labelneigs.end()) { continue; }
            std::vector<uint> &neigs = ite1->second;
            for(uint j=0;j<neigs.size();++j){
                bool found = false;
                for(uint k=0;k<preEmbLen;++k){
                    if(neigs[j]==emb[k]){
                        found = true;
                        break;
                    }
                }
                if(!found){
                    if(isExtRestrict){
                        uint tmppos = posOfVisitedV[restricts[evid]];
                        uint tmpvid = emb[tmppos];
                        if(tmpvid>neigs[j]){
                            for(uint l=0;l<preEmbLen;++l){
                                partialRes[rowIndex*embLen+l] = emb[l];
                            }
                            partialRes[rowIndex*embLen+preEmbLen] = neigs[j];
                            ++rowIndex;
                            ++totRowNum;
                            if(rowIndex==maxRowNum){
                                outputfile.write((char*)partialRes,maxRowNum*embLen*sizeof(uint));
                                rowIndex=0;
                                std::cout<<"write to file done,embLen="<<embLen<<" tot="<<totRowNum<<std::endl;
                            }
                        }
                    }else{
                        for(uint l=0;l<preEmbLen;++l){
                            partialRes[rowIndex*embLen+l] = emb[l];
                        }
                        partialRes[rowIndex*embLen+preEmbLen] = neigs[j];
                        ++rowIndex;
                        ++totRowNum;
                        if(rowIndex==maxRowNum){
                            outputfile.write((char*)partialRes,maxRowNum*embLen*sizeof(uint));
                            rowIndex=0;
                            std::cout<<"write to file done,embLen="<<embLen<<" tot="<<totRowNum<<std::endl;
                        }
                    }
                }
            }
        }
        outputfile.write((char*)partialRes,rowIndex*embLen*sizeof(uint));
        free(emb);
    }else{
        uint maxRowNum = uintsNum/(preEmbLen+2);
        uint svid1 = pairs[0],evid1 = pairs[1],svid2 = pairs[2],evid2 = pairs[3];
        uint index = posOfVisitedV.size();
        posOfVisitedV.insert(std::pair<uint,uint>(evid1,index));
        posOfVisitedV.insert(std::pair<uint,uint>(evid2,index+1));
        uint svid1label = queryLabel[svid1],evid1label = queryLabel[evid1];
        uint svid2label = queryLabel[svid2],evid2label = queryLabel[evid2];
        uint64_t key1 = edgeLabel, key2 = edgeLabel;
        key1 = (key1<<32)|evid1label; 
        key2 = (key2<<32)|evid2label;

        bool isExt1Restrict = restricts[evid1]!=0xffffffff;
        bool isExt1and2Restrict = restricts[evid2]!=0xffffffff?restricts[evid2]==evid1:false;
        bool isExt2Restrict = isExt1and2Restrict?false:restricts[evid2]!=0xffffffff;
        bool isRecord1 = (repeatFlag[key1]&0x0000ffff)>1;
        bool isRecord2 = (repeatFlag[key2]&0x0000ffff)>1;
        bool useRecord1 = (repeatFlag[key1]&0x0000ffff)<((repeatFlag[key1]>>16)&0x0000ffff);
        bool useRecord2 = (repeatFlag[key2]&0x0000ffff)<((repeatFlag[key2]>>16)&0x0000ffff);
        bool isExt1and2SameLabel = evid1label==evid2label;
        bool isExt1and2SameSrc = svid1==svid2;

        std::vector<uint> auxArray(8);
        auxArray[0]=evid1label;auxArray[1]=evid2label;auxArray[4]=posOfVisitedV[svid1];auxArray[5]=posOfVisitedV[svid2];
        auxArray[2] = isExt1Restrict?posOfVisitedV[restricts[evid1]]:0;
        auxArray[3] = isExt2Restrict?posOfVisitedV[restricts[evid2]]:0;
        if(equalPos.find(evid1label)!=equalPos.end()){
            std::vector<uint> &tmpv = equalPos[evid1label];
            auxArray[6]=tmpv.size();
            auxArray.insert(auxArray.end(),tmpv.begin(),tmpv.end());
            uint index = posOfVisitedV[evid1];
            tmpv.push_back(index);
        }else{
            std::vector<uint> tmpv;
            uint index = posOfVisitedV[evid1];
            tmpv.push_back(index);
            equalPos.insert(std::pair<uint,std::vector<uint>>(evid1label,tmpv));
            auxArray[6]=0;
        }
        if(!isExt1and2SameLabel){
            if(equalPos.find(evid2label)!=equalPos.end()){
                std::vector<uint> &tmpv = equalPos[evid2label];
                uint tmpsize = tmpv.size();
                auxArray[7]=tmpsize+auxArray[6];
                auxArray.insert(auxArray.end(),tmpv.begin(),tmpv.end());
                uint index = posOfVisitedV[evid2];
                tmpv.push_back(index);
            }else{
                auxArray[7]=0+auxArray[6];
                std::vector<uint> tmpv;
                uint index = posOfVisitedV[evid2];
                tmpv.push_back(index);
                equalPos.insert(std::pair<uint,std::vector<uint>>(evid2label,tmpv));
            }
        }else{
            auxArray[7] = auxArray[6];
            std::vector<uint> &tmpv = equalPos[evid1label];
            uint index = posOfVisitedV[evid2];
            tmpv.push_back(index);
        }
        uint embLen = posOfVisitedV.size(),rowIndex=0;
        uint *emb = (uint *)malloc(32*sizeof(uint));
        std::cout<<"begin cpu search "<<partialRowNum<<std::endl;
        hash_map<uint,hash_map<uint,std::vector<uint>>> &edgelabepart = reorderGraph[edgeLabel];
        for(uint i=0;i<partialRowNum;++i){
            if(i%20000==0){
                std::cout<<"extension, tot="<<partialRowNum<<" iter="<<i<<std::endl;
            }

            inputfile.read((char*)emb,preEmbLen*sizeof(uint));
            uint svid1pos = posOfVisitedV[svid1];
            uint svid2pos = posOfVisitedV[svid2];
            uint tmpsvid1 = emb[svid1pos];
            uint tmpsvid2 = emb[svid2pos];
            auto itea = edgelabepart.find(tmpsvid1);
            auto iteb = edgelabepart.find(tmpsvid2);
            if(itea==edgelabepart.end() || iteb==edgelabepart.end()) { continue; }
            hash_map<uint,std::vector<uint>> &labelneigs1 = itea->second;
            hash_map<uint,std::vector<uint>> &labelneigs2 = iteb->second;
            auto ite1 = labelneigs1.find(evid1label);
            auto ite2 = labelneigs2.find(evid2label);
            if(ite1==labelneigs1.end() || ite2==labelneigs2.end()) { continue; }
            std::vector<uint> &neigs1 = ite1->second;
            std::vector<uint> &neigs2 = ite2->second;
            for(uint j=0;j<neigs1.size();++j){
                bool found = false;
                for(uint k=0;k<preEmbLen;++k){
                    if(neigs1[j]==emb[k]){
                        found = true;
                        break;
                    }
                }
                if(found){
                    continue;
                }
                if(isExt1Restrict){
                    uint tmppos = posOfVisitedV[restricts[evid1]];
                    if(neigs1[j]>emb[tmppos]){
                        continue;
                    }
                }
                for(uint l=0;l<neigs2.size();++l){
                    found = false;
                    for(uint k=0;k<preEmbLen;++k){
                        if(neigs2[l]==emb[k]){
                            found = true;
                            break;
                        }
                    }
                    if(!found && neigs2[l]!=neigs1[j]){
                        if(isExt2Restrict){
                            uint tmppos = posOfVisitedV[restricts[evid2]];
                            if(neigs2[l]>emb[tmppos]){
                                continue;
                            }
                        }
                        if(isExt1and2Restrict && neigs2[l]>neigs1[j]){
                            continue;
                        }
                        for(uint h=0;h<preEmbLen;++h){
                            partialRes[rowIndex*embLen+h] = emb[h];
                        }
                        partialRes[rowIndex*embLen+preEmbLen] = neigs1[j];
                        partialRes[rowIndex*embLen+preEmbLen+1] = neigs2[l];
                        ++rowIndex;
                        ++totRowNum;
                        if(rowIndex==maxRowNum){
                            outputfile.write((char*)partialRes,maxRowNum*embLen*sizeof(uint));
                            rowIndex=0;
                        }
                    }
                }
            }
        }
        outputfile.write((char*)partialRes,rowIndex*embLen*sizeof(uint));
        free(emb);
        std::cout<<"cpu search done"<<std::endl;
    }
    outputfile.close();
    inputfile.close();
    return totRowNum;
}

void sortReducePairs(std::vector<uint>& pairs){
    uint size = pairs.size()/2;
    for(int i=0;i<size;++i){
        uint tmpsize = size-i;
        uint minsvid = pairs[0];
        uint minevid = pairs[1];
        uint minindex=0;
        for(int j=1;j<tmpsize;++j){
            uint cursvid = pairs[j*2];
            uint curevid = pairs[j*2+1];
            if(cursvid<minsvid){
                minsvid = cursvid;
                minevid = curevid;
                minindex = j;
            }else if(cursvid==minsvid){
                if(queryLabel[cursvid]<queryLabel[minsvid]){
                    minsvid = cursvid;
                    minevid = curevid;
                    minindex = j;
                }else if(queryLabel[curevid]==queryLabel[minevid]){
                    if(curevid<minevid){
                        minsvid = cursvid;
                        minevid = curevid;
                        minindex = j;
                    }
                }
            }
        }
        uint lastindex = tmpsize-1;
        uint tmp = pairs[minindex*2];
        pairs[minindex*2] = pairs[lastindex*2];
        pairs[lastindex*2] = tmp;
        tmp = pairs[minindex*2+1];
        pairs[minindex*2+1] = pairs[lastindex*2+1];
        pairs[lastindex*2+1] = tmp;
    }
}

//the format for indexForIndex[0:2*256] is vs1,vs2,vs3,len1,len2,len3. len1 is the numberr of vertices in interval 1, len2
//is the number of vertices in all previsous intervals (including interval 2).
//the format for indexForIndex[2*256:] is totUintNum(uint)(only include pairs),(start1,end1),(start2,end2),evid labels(16),recordPos(16).
// first, all same start are grouped together, in each group
//all end with same label are grouped. each number uses 5 bits, each time we at most process 16 pairs. we use one uint to represent use
//use record and isrecord.

uint reduction(Phase &curPhase,bool isLastPhase, std::string inputfilename, std::string outputfilename,
    uint *partialRes,uint uintsNum, uint partialRowNum){
    std::ifstream inputfile(inputfilename,std::ios::in|std::ios::binary);
    std::ofstream outputfile(outputfilename,std::ios::out|std::ios::binary);
    uint edgelabel = curPhase.edgeLabel;
    uint64_t key = edgelabel;
    key = (key<<32) | 0xffffffff;
    uint embLen = posOfVisitedV.size();
    uint maxRowNum = uintsNum/embLen;
    uint *recordpos = recordPosData[key].addr;
    uint recordsize = recordPosData[key].size;
    std::vector<uint> &pairs = curPhase.pairs;
    sortReducePairs(pairs);
    if(pairs.size()>32){
        std::cout<<"not process this situation now"<<std::endl;
        return 0;
    }
    uint totUintNum = (pairs.size()+5)/6;
    std::vector<uint> auxArray(1+6+16+16+1,0);
    auxArray[0]=totUintNum;
    uint groupsize = 0,tmp=0,index=1;
    for(uint i=0;i<pairs.size()/2;++i){
        uint svid = pairs[i*2];
        uint evid = pairs[i*2+1];
        std::cout<<"after sort: svid="<<svid<<" evid="<<evid<<std::endl;
        uint evidlabel = queryLabel[evid];
        auxArray[1+6+i] = evidlabel;
        key = edgelabel;
        key = (key<<32)|0xffffffff;
        uint *tmprecordpos = recordPosData[key].addr;
        auxArray[1+6+16+i] = tmprecordpos-recordpos;
        uint spos = posOfVisitedV[svid];
        uint epos = posOfVisitedV[evid];
        tmp = tmp | (epos<<groupsize);
        tmp = tmp | (spos<<(groupsize+5));
        groupsize = groupsize+10;
        if(groupsize==30){
            auxArray[index]=tmp;
            tmp=0;
            groupsize=0;
            ++index;
        }
    }
    if(tmp>0){
        auxArray[index] = tmp;
    }
    uint rowIndex=0, totRowNum=0;
    hash_map<uint,hash_map<uint,std::vector<uint>>> &edgelabepart = reorderGraph[edgelabel];
    uint *emb = (uint*)malloc(32*sizeof(uint));
    for(uint i=0;i<partialRowNum;++i){
        if(i%20000==0){
            std::cout<<"reduction, tot="<<partialRowNum<<" iter="<<i<<std::endl;
        }
        inputfile.read((char*)emb,embLen*sizeof(uint));
        bool allfound = true;
        for(uint j=0;j<pairs.size()/2;++j){
            uint svid = pairs[j*2];
            uint evid = pairs[j*2+1];
            uint evidlabel = queryLabel[evid];
            uint spos = posOfVisitedV[svid];
            uint epos = posOfVisitedV[evid];
            svid = emb[spos];
            evid = emb[epos];
            auto ite = edgelabepart.find(svid);
            if(ite==edgelabepart.end()) { allfound = false; break; }
            hash_map<uint,std::vector<uint>> &labelneigs = ite->second;
            auto itea = labelneigs.find(evidlabel);
            if(itea==labelneigs.end()) { allfound = false; break; }
            std::vector<uint> &neigs = itea->second;
            bool found = false;
            for(uint k=0;k<neigs.size();++k){
                if(neigs[k]==evid){
                    found = true;
                    break;
                }
            }
            if(!found){
                allfound = false;
                break;
            }
        }
        if(allfound){
            for(uint l=0;l<embLen;++l){
                partialRes[rowIndex*embLen+l] = emb[l]; 
            }
            ++rowIndex;
            ++totRowNum;
            if(rowIndex==maxRowNum){
                outputfile.write((char*)partialRes,maxRowNum*embLen*sizeof(uint));
                rowIndex = 0;
            }
        }
    }
    outputfile.write((char*)partialRes,rowIndex*embLen*sizeof(uint));
    inputfile.close();
    outputfile.close();
    free(emb);
    return totRowNum;
}


void printPartialEmbCPU(std::string inputfilename, std::string outputfilename, uint partialRowNum, uint embLen){
    std::ifstream inputfile(inputfilename,std::ios::in|std::ios::binary);
    std::ofstream outputfile(outputfilename,std::ios::out);
    uint *emb=(uint*)malloc(32*sizeof(uint));
    for(uint i=0;i<partialRowNum;++i){
        inputfile.read((char*)emb,embLen*sizeof(uint));
        for(uint j=0;j<embLen;++j){
            outputfile<<emb[j]<<" ";
        }
        outputfile<<std::endl;
    }
    outputfile.close();
    inputfile.close();
    free(emb);
}

int main(int argc, char *argv[]){
	bool check;
	if(argc==3){
		check=false;
	}else if(argc==5){
		check=true;
	}else{
		std::cout<<"wrong parameters"<<std::endl;
		return -1;
	}
	//suffix .reorder, .query
    std::string reorderFileName = std::string(argv[1]);
    std::string queryFileName = std::string(argv[2]);
    //std::string checkFileName = path+baseFileName+".g";

	if(check){
		//suffix .mygraph, .metadata
    		std::string inputFileName = std::string(argv[3]);
    		std::string metaFileName = std::string(argv[4]);
    		allEdgeLabelPartitions = loadDataBase(inputFileName, metaFileName);
    		vLabel = allEdgeLabelPartitions[0xffffffff];
	}
    std::cout<<"generating graph"<<std::endl;
    std::ifstream reorderfile(reorderFileName,std::ios::in);
    uint totreorderlabelpart;
    reorderfile>> totreorderlabelpart;
    std::cout<<"total label partitions: "<<totreorderlabelpart<<std::endl;
    for(uint i=0;i<totreorderlabelpart;++i){
        uint svidnum;
        reorderfile >> svidnum;
        std::cout<<"svidnum: "<<svidnum<<std::endl;
        hash_map<uint,hash_map<uint,std::vector<uint> > > labelpart(svidnum);
        for(uint j=0;j<svidnum;++j){
            uint svid, evidlabelnum;
            reorderfile>>svid>>evidlabelnum;
            //std::cout<<"svid:"<<svid<<" evidlabelnum:"<<evidlabelnum<<std::endl;
            hash_map<uint,std::vector<uint>> &evidlabelneigs=labelpart[svid];
            for(uint k=0;k<evidlabelnum;++k){
                uint evidlabel, neigsnum;
                reorderfile>>evidlabel>>neigsnum;
                //std::cout<<"evidlabel:"<<evidlabel<<" neigsum:"<<neigsnum<<std::endl;
                std::vector<uint> neigs(neigsnum);
                for(uint l=0;l<neigsnum;++l){
                    uint neigvid;
                    reorderfile>>neigvid;
                    neigs[l] = neigvid;
                }
                evidlabelneigs[evidlabel] = neigs;
            }
            //std::cout<<"insert "<<svid<<" done"<<std::endl;
        }
        reorderGraph.push_back(labelpart);
    }
    reorderfile.close();
	std::cout<<"Generte graph done. Next, generate query order."<<std::endl;

    genMatchOrder(queryFileName,matchOrder,repeatFlag,restricts,edgeLabelStat,queryLabel);
    Phase &tmpPhase = matchOrder[0];
    for(auto ite=restricts.begin();ite!=restricts.end();++ite){
        uint key = ite->first;
        uint value = ite->second;
        std::cout<<"restrict: "<<key<<"<"<<value<<std::endl;
    }
    for(uint i=0;i<matchOrder.size();++i){
        Phase &phase = matchOrder[i];
        std::cout << "phase "<<i<<" type "<<phase.phaseType<<" edge label:"<<phase.edgeLabel<<std::endl;
        std::vector<uint> &pairs = phase.pairs;
        for(uint j=0;j<pairs.size()/2;++j){
            std::cout<<"svid:"<<pairs[j*2]<<" evid:"<<pairs[j*2+1]<<" ";
        }
        std::cout<<std::endl;
    }
    uint label = tmpPhase.edgeLabel;
    uint embLen;
    if(matchOrder[0].pairs.size()==2){
        embLen=2;
    }else{
        embLen=3;
    }
    std::ofstream cpuresult("tmpcpuresult.txt",std::ios::out|std::ios::binary);
    uint *partialResult = (uint*)malloc(sizeof(uint)*128000);
    uint totwriteNum = genInit(cpuresult,partialResult,128000);
    std::cout<<"after init, tot="<<totwriteNum<<std::endl;
    cpuresult.close();
    //printPartialEmbCPU("tmpcpuresult.txt","partialEmbCPU_0.txt",totwriteNum,embLen);

    uint i;
    std::string inputfilename = "tmpcpuresult.txt";
    std::string outputfilename = "tmpcpuresult1.txt";
    for(i=1;i<matchOrder.size();++i){
        bool isLastPhase = i==(matchOrder.size()-1);
        Phase &curPhase = matchOrder[i];
        if(curPhase.phaseType=='E'){
            if(curPhase.pairs.size()==2){ embLen += 1; }
            else if(curPhase.pairs.size()==4){ embLen += 2; }
            std::cout<<"in Extension"<<std::endl;
            totwriteNum = genExtEmb(curPhase,isLastPhase,inputfilename,outputfilename,partialResult,128000,totwriteNum);
            std::string tmp = inputfilename;
            inputfilename = outputfilename;
            outputfilename = tmp;
            //std::string filename = "partialEmbCPU_ext_"+std::to_string(i)+".txt";
            //std::cout<<"write to "<<filename<<std::endl;
            //printPartialEmbCPU(inputfilename,filename,totwriteNum,embLen);
            //std::cout<<"write to "<<filename<<" done"<<std::endl;
        }else{
            totwriteNum = reduction(curPhase,isLastPhase,inputfilename,outputfilename,partialResult,128000,totwriteNum);
            std::cout<<"in Reduction:"<<totwriteNum<<std::endl;
            std::string tmp = inputfilename;
            inputfilename = outputfilename;
            outputfilename = tmp;
            //std::string filename = "partialEmbCPU_red_"+std::to_string(i)+".txt";
            //std::cout<<"write to "<<filename<<std::endl;
            //printPartialEmbCPU(inputfilename,filename,totwriteNum,embLen);
            //std::cout<<"write to "<<filename<<" done"<<std::endl;
        }
    }
    std::cout<<"write to final result cpu, tot="<<totwriteNum<<std::endl;
    printPartialEmbCPU(inputfilename,"partialResultsCPU.txt",totwriteNum,embLen);
    return 0;
}
