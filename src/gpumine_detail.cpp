#include <iostream>
#include <fstream>
#include <cuda_runtime.h>
#include <algorithm>
#include <hash_map>
#include <ctime>
#include "common.h"
#include "loadDataBase.h"
#include "genMatchOrder.h"
#include "miningPattern.h"
#include "reduction.h"

using namespace __gnu_cxx;
class recordPosInfo{
public:
    uint size;
    uint* addr;
};
size_t gpuAvailSpace;
uint maxRowNum;
std::vector<uint*> partResults;
std::vector<uint> partResultRowNum;
int GPU_SM_NUM;
int debugNum=0;
int debugNumend=0;
std::map<uint,uint> posOfVisitedV;
std::map<uint64_t,recordPosInfo> recordPosData;
std::map<uint,std::vector<uint>> equalPos;
std::vector<Phase> matchOrder;
std::map<uint64_t,uint> repeatFlag;
std::map<uint,std::map<uint,bool>> edgeLabelStat;
std::map<uint,uint> restricts;
std::map<uint,uint*> allEdgeLabelPartitions;
std::map<uint,uint> queryLabel;
uint *vLabel;


void printPartialEmbGPUDebug(uint *partialEmb_dev, uint &partialEmbRowNum, uint embLen, std::string filename){
    //partialEmbRowNum = debugNumend-debugNum;
    uint *partialEmb_host = (uint*)malloc(sizeof(uint)*partialEmbRowNum*embLen);
    std::ifstream inputfile(filename,std::ios::in);
    for(uint i=0;i<partialEmbRowNum;++i){
        for(uint j=0;j<embLen;++j){
            inputfile>>partialEmb_host[i*embLen+j];
        }
    }
    cudaMemcpy(partialEmb_dev,partialEmb_host,sizeof(uint)*partialEmbRowNum*embLen,cudaMemcpyHostToDevice);
    inputfile.close();
    free(partialEmb_host);
}

void printPartialEmbGPU(uint *partialEmb_dev, uint partialEmbRowNum, uint embLen, std::string filename){
    uint *partialEmb_host = (uint*)malloc(sizeof(uint)*partialEmbRowNum*embLen);
    cudaMemcpy(partialEmb_host,partialEmb_dev,sizeof(uint)*partialEmbRowNum*embLen,cudaMemcpyDeviceToHost);
    std::ofstream outputfile;
    outputfile.open(filename, std::ios::out);
    if (!outputfile.is_open()) return;
    for(uint i=0;i<partialEmbRowNum;++i){
        outputfile<<partialEmb_host[i*embLen+0];
        for(uint j=1;j<embLen;++j){
            outputfile<<" "<<partialEmb_host[i*embLen+j];
        }
        outputfile<<std::endl;
    }
    outputfile.close();
    free(partialEmb_host);
}
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

uint genInit(uint *vLabel_dev, uint *writeRowNum_dev, uint *edgeLabelPartition_dev, uint *neighborsData_dev, 
    uint *baseAddr_dev, uint stride, uint intervalNum, uint totSrcNum,bool isContinue){

    Phase &firstPhase = matchOrder[0];
    uint edgeLabel = firstPhase.edgeLabel;
    std::vector<uint> &pairs = firstPhase.pairs;
    if(pairs.size()==2){
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

        uint *recordPos_dev=NULL;
        if(isRecord){
            recordPos_dev = baseAddr_dev+stride;
            stride = stride+totSrcNum*2;
        }
        //the generated rows can not occupy half of the available sapec
        maxRowNum = ((gpuAvailSpace-stride)/2-38*1500)/2;
        initEmb_2V_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,recordPos_dev,baseAddr_dev+stride,intervalNum,totSrcNum,svidlabel,evidlabel,isExtRestrict,isRecord,isContinue,maxRowNum);
        if(!isContinue){
            constructRecordPos();
        }
        if(isRecord){
            repeatFlag[key] = repeatFlag[key]-1;
            uint *recordPos_host = recordPosData[key].addr;
            cudaMemcpyAsync(recordPos_host,recordPos_dev,sizeof(uint)*totSrcNum*2,cudaMemcpyDeviceToHost);
        }
    }else{
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
        uint *recordPos1_dev = NULL, *recordPos2_dev = NULL;
        if(isRecord1){ recordPos1_dev = baseAddr_dev+stride; stride = stride+totSrcNum*2; }
        if(!isExt1and2SameLabel && isRecord2){ recordPos2_dev = baseAddr_dev+stride; stride = stride+totSrcNum*2; }
        maxRowNum = ((gpuAvailSpace-stride)/2-38*1500)/3;
        initEmb_3V_1and2_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,recordPos1_dev,recordPos2_dev,
            baseAddr_dev+stride,intervalNum,totSrcNum,svid1label,evid1label,evid2label,isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isExt1and2SameLabel,isRecord1,isRecord2,isContinue,maxRowNum);
        //cudaDeviceSynchronize();
        //cudaError_t err = cudaGetLastError();
        //std::cout<<"after initemb_3v "<<cudaGetErrorString(err)<<std::endl;
        if(!isContinue){
            constructRecordPos();
        }
        uint *recordPos_host;
        if(isRecord1){
            repeatFlag[key1] = repeatFlag[key1]-1;
            recordPos_host = recordPosData[key1].addr;
            std::cout<< "size is:"<<recordPosData[key1].size<<" key1="<<key1<<std::endl;
            cudaMemcpyAsync(recordPos_host,recordPos1_dev,sizeof(uint)*totSrcNum*2,cudaMemcpyDeviceToHost);
        }
        //cudaDeviceSynchronize();
        //err = cudaGetLastError();
        //std::cout<<"after copyasnc "<<cudaGetErrorString(err)<<" totsrcnum="<<totSrcNum<<std::endl;
        if(isRecord2){
            repeatFlag[key1] = repeatFlag[key1]-1;
            if(!isExt1and2SameLabel){
                uint *recordPos_host = recordPosData[key2].addr;
                cudaMemcpyAsync(recordPos_host,recordPos2_dev,sizeof(uint)*totSrcNum*2,cudaMemcpyDeviceToHost);
            }
        }
        //cudaDeviceSynchronize();
        //err = cudaGetLastError();
        //std::cout<<"after copyasnc2 "<<cudaGetErrorString(err)<<std::endl;
    }
    return stride;
}

uint* genExtEmb(uint *vLabel_dev, uint *writeRowNum_dev, uint *writeRowNum, uint *edgeLabelPartition_dev, uint *neighborsData_dev,
    uint *baseAddr_dev, uint &stride, uint *partialEmb_dev, uint intervalNum, uint totSrcNum, uint partialRowNum, 
    Phase &curPhase, bool isLastPhase, bool isContinue){
	cudaError_t err;
    uint edgeLabel = curPhase.edgeLabel;
    std::vector<uint> &pairs = curPhase.pairs;
    if(pairs.size()==2){
        uint svid = pairs[0], evid = pairs[1];
        uint index = posOfVisitedV.size();
        uint preEmbLen = posOfVisitedV.size();
        posOfVisitedV.insert(std::pair<uint,uint>(evid,index));
        uint embLen = posOfVisitedV.size();
        uint svidlabel = queryLabel[svid], evidlabel = queryLabel[evid];
        uint64_t key = edgeLabel;
        key = (key<<32)|evidlabel;
        bool isRecord = (repeatFlag[key]&0x0000ffff)>1;
        bool isExtRestrict = restricts[evid]!=0xffffffff;
        //std::cout<<"restrict:"<<evid<<"<"<<restricts[evid]<<" full 1:"<<0xffffffff<<" isExtRestrict="<<isExtRestrict<<std::endl;
        bool useRecord = (repeatFlag[key]&0x0000ffff)<((repeatFlag[key]>>16)&0x0000ffff);

        uint *recordPos_dev;
        if(isRecord || useRecord){
            recordPos_dev = baseAddr_dev+stride;
            stride = stride+totSrcNum*2;
            cudaMemcpyAsync(recordPos_dev,recordPosData[key].addr,sizeof(uint)*totSrcNum*2,cudaMemcpyHostToDevice);
        }

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
        uint *auxArray_dev = baseAddr_dev+stride;
        stride=stride+auxArray.size();
        cudaMemcpyAsync(auxArray_dev,auxArray.data(),sizeof(uint)*auxArray.size(),cudaMemcpyHostToDevice);
        uint lesspos = isExtRestrict?posOfVisitedV[restricts[evid]]:0;
        //std::cout << "intervalNum="<<intervalNum << std::endl;
        maxRowNum = (gpuAvailSpace-stride-partialRowNum*preEmbLen-1500*38)/embLen;
        extEmb_1V_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,recordPos_dev,auxArray_dev,
            baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,embLen,isExtRestrict,isRecord,useRecord,isLastPhase,false,maxRowNum);
        //cudaDeviceSynchronize();
        //std::cout<<"first extEmb_1V_NoHash is done"<<std::endl;
        cudaMemcpy(writeRowNum,writeRowNum_dev,sizeof(uint)*2,cudaMemcpyDeviceToHost);
        uint tmppartialRowNum = writeRowNum[0];
        //std::cout<<"first extEmb_1V_NoHash is done:"<<tmppartialRowNum<<" max:"<<maxRowNum<<std::endl;
        //printPartialEmbGPU(baseAddr_dev+stride, tmppartialRowNum,embLen,"testext1.txt");
	/*if(tmppartialRowNum>=maxRowNum){
		std::cout<<"extEmb_1V_NoHash exceeds the maximum: genrownum="
			 <<tmppartialRowNum<<" maxrownum="<<maxRowNum<<" remove="
			 <<writeRowNum[1]<<std::endl;
	}*/
        std::vector<uint*> tmppartResults;
        std::vector<uint> tmppartResultRowNum;
        while(tmppartialRowNum>=maxRowNum){
            tmppartialRowNum = tmppartialRowNum-writeRowNum[1];
            //std::cout<<"over size detec in extension 2, write row num:"<<tmppartialRowNum<<" maxRowNum:"<<maxRowNum<<std::endl;
            uint allocnum = ((tmppartialRowNum*embLen+31)>>5)<<5;
            uint *tmpaddr = (uint*)malloc(sizeof(uint)*allocnum);
            cudaMemcpyAsync(tmpaddr,baseAddr_dev+stride,sizeof(uint)*tmppartialRowNum*embLen,cudaMemcpyDeviceToHost);
            writeRowNum[0] = 0;
            writeRowNum[1] = 0;
            cudaMemcpyAsync(writeRowNum_dev,writeRowNum,sizeof(uint)*2,cudaMemcpyHostToDevice);
            tmppartResults.push_back(tmpaddr);
            tmppartResultRowNum.push_back(tmppartialRowNum);
            //cudaDeviceSynchronize();
            //err = cudaGetLastError();
            //std::cout<<"in extension 2 "<<cudaGetErrorString(err)<<" before extEmb_1V"<<std::endl;
	//	std::cout<<"isExtRestrict="<<isExtRestrict<<" isRecord="<<isRecord<<" useRecord="<<useRecord<<" isLastPhase="<<isLastPhase<<std::endl;
            extEmb_1V_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,recordPos_dev,auxArray_dev,
                baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,embLen,isExtRestrict,isRecord,useRecord,isLastPhase,true,maxRowNum);           
            cudaMemcpy(writeRowNum,writeRowNum_dev,sizeof(uint)*2,cudaMemcpyDeviceToHost);
            tmppartialRowNum = writeRowNum[0];
            //err = cudaGetLastError();
            //std::cout<<"in extension 2 "<<cudaGetErrorString(err)<<" after extEmb_1V:"<<tmppartialRowNum<<std::endl;
        }
        if(partResults.size()>0){
            //std::cout<<"there are other parts to be processed"<<std::endl;
            uint allocnum = ((tmppartialRowNum*embLen+31)>>5)<<5;
            uint *tmpaddr = (uint*)malloc(sizeof(uint)*allocnum);
            cudaMemcpyAsync(tmpaddr,baseAddr_dev+stride,sizeof(uint)*tmppartialRowNum*embLen,cudaMemcpyDeviceToHost);
            writeRowNum[0] = 0;
            writeRowNum[1] = 0;
            cudaMemcpyAsync(writeRowNum_dev,writeRowNum,sizeof(uint)*2,cudaMemcpyHostToDevice);
            tmppartResults.push_back(tmpaddr);
            tmppartResultRowNum.push_back(tmppartialRowNum);
            for(uint j=0;j<partResults.size();++j){
                tmpaddr = partResults[j];
                partialRowNum = partResultRowNum[j];
                partialEmb_dev = baseAddr_dev+gpuAvailSpace-partialRowNum*preEmbLen-32;
                cudaMemcpyAsync(partialEmb_dev,tmpaddr,sizeof(uint)*partialRowNum*preEmbLen,cudaMemcpyHostToDevice);
                maxRowNum = (gpuAvailSpace-stride-partialRowNum*preEmbLen-1500*38)/embLen;
                extEmb_1V_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,recordPos_dev,auxArray_dev,
                    baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,embLen,isExtRestrict,isRecord,useRecord,isLastPhase,false,maxRowNum);
                cudaMemcpy(writeRowNum,writeRowNum_dev,sizeof(uint)*2,cudaMemcpyDeviceToHost);
                free(tmpaddr);
                tmppartialRowNum = writeRowNum[0];
                while(tmppartialRowNum>=maxRowNum){
                    tmppartialRowNum = tmppartialRowNum-writeRowNum[1];
                    uint allocnum = ((tmppartialRowNum*embLen+31)>>5)<<5;
                    uint *tmpaddr = (uint*)malloc(sizeof(uint)*allocnum);
                    cudaMemcpyAsync(tmpaddr,baseAddr_dev+stride,sizeof(uint)*tmppartialRowNum*embLen,cudaMemcpyDeviceToHost);
                    writeRowNum[0] = 0;
                    writeRowNum[1] = 0;
                    cudaMemcpyAsync(writeRowNum_dev,writeRowNum,sizeof(uint)*2,cudaMemcpyHostToDevice);
                    tmppartResults.push_back(tmpaddr);
                    tmppartResultRowNum.push_back(tmppartialRowNum);
                    extEmb_1V_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,recordPos_dev,auxArray_dev,
                        baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,embLen,isExtRestrict,isRecord,useRecord,isLastPhase,true,maxRowNum);           
                    cudaMemcpy(writeRowNum,writeRowNum_dev,sizeof(uint)*2,cudaMemcpyDeviceToHost);
                    tmppartialRowNum = writeRowNum[0];
                }
            }
        }
        partResults = tmppartResults;
        partResultRowNum = tmppartResultRowNum;
        if(isRecord){
            repeatFlag[key] = repeatFlag[key]-1;
            uint *recordPos_host = recordPosData[key].addr;
            cudaMemcpyAsync(recordPos_host,recordPos_dev,sizeof(uint)*totSrcNum*2,cudaMemcpyDeviceToHost);
        }
    }else{
        uint svid1 = pairs[0],evid1 = pairs[1],svid2 = pairs[2],evid2 = pairs[3];
        uint index = posOfVisitedV.size();
        uint preEmbLen = posOfVisitedV.size();
        posOfVisitedV.insert(std::pair<uint,uint>(evid1,index));
        posOfVisitedV.insert(std::pair<uint,uint>(evid2,index+1));
        uint embLen = posOfVisitedV.size();
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

        uint *recordPos1_dev = NULL, *recordPos2_dev=NULL;
        if(isRecord1 || useRecord1){
            recordPos1_dev = baseAddr_dev+stride;
            stride = stride+totSrcNum*2;
            cudaMemcpyAsync(recordPos1_dev,recordPosData[key1].addr,sizeof(uint)*totSrcNum*2,cudaMemcpyHostToDevice);
        }
        if(!isExt1and2SameLabel && (isRecord2 || useRecord2)){
            recordPos2_dev = baseAddr_dev+stride;
            stride = stride+totSrcNum*2;
            cudaMemcpyAsync(recordPos2_dev,recordPosData[key2].addr,sizeof(uint)*totSrcNum*2,cudaMemcpyHostToDevice);
        }

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

        uint *auxArray_dev = baseAddr_dev+stride;
        stride=stride+auxArray.size();
        cudaMemcpyAsync(auxArray_dev,auxArray.data(),sizeof(uint)*auxArray.size(),cudaMemcpyHostToDevice);
        //cudaDeviceSynchronize();
        //cudaError_t err = cudaGetLastError();
        //std::cout<<"200 "<<cudaGetErrorString(err)<<std::endl;
        maxRowNum = (gpuAvailSpace-stride-partialRowNum*preEmbLen-1500*38)/embLen;
        if(isExt1and2SameSrc && isExt1and2SameLabel) {
            extEmb_2V_1src_1and2_sameLabel_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,
                recordPos1_dev,auxArray_dev,baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,embLen,
                isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isRecord1,useRecord1,false,maxRowNum);
        }else if(isExt1and2SameSrc && !isExt1and2SameLabel){
            extEmb_2V_1src_1and2_notSameLabel_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,
                recordPos1_dev,recordPos2_dev,auxArray_dev,baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,
                embLen,isExt1Restrict,isExt2Restrict,isRecord1,isRecord2,useRecord1,useRecord2,false,maxRowNum);
        }else if(!isExt1and2SameSrc && isExt1and2SameLabel) {
            extEmb_2V_2src_1and2_sameLabel_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,
                recordPos1_dev,auxArray_dev,baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,embLen,
                isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isRecord1,useRecord1,false,maxRowNum);
        }else{
            extEmb_2V_2src_1and2_notSameLabel_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,
                recordPos1_dev,recordPos2_dev,auxArray_dev,baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,
                embLen,isExt1Restrict,isExt2Restrict,isRecord1,isRecord2,useRecord1,useRecord2,false,maxRowNum);
        }
        //cudaDeviceSynchronize();
        //err = cudaGetLastError();
        //std::cout<<"300 "<<cudaGetErrorString(err)<<std::endl;
        cudaMemcpy(writeRowNum,writeRowNum_dev,sizeof(uint)*2,cudaMemcpyDeviceToHost);
        uint tmppartialRowNum = writeRowNum[0];
        std::vector<uint*> tmppartResults;
        std::vector<uint> tmppartResultRowNum;
        while(tmppartialRowNum>=maxRowNum){
            //std::cout<<"over size detec in extension 4"<<std::endl;
            tmppartialRowNum = tmppartialRowNum-writeRowNum[1];
            uint allocnum = ((tmppartialRowNum*embLen+31)>>5)<<5;
            uint *tmpaddr = (uint*)malloc(sizeof(uint)*allocnum);
            cudaMemcpyAsync(tmpaddr,baseAddr_dev+stride,sizeof(uint)*tmppartialRowNum*embLen,cudaMemcpyDeviceToHost);
            writeRowNum[0] = 0;
            writeRowNum[1] = 0;
            cudaMemcpyAsync(writeRowNum_dev,writeRowNum,sizeof(uint)*2,cudaMemcpyHostToDevice);
            tmppartResults.push_back(tmpaddr);
            tmppartResultRowNum.push_back(tmppartialRowNum);
            if(isExt1and2SameSrc && isExt1and2SameLabel) {
                extEmb_2V_1src_1and2_sameLabel_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,
                    recordPos1_dev,auxArray_dev,baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,embLen,
                    isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isRecord1,useRecord1,true,maxRowNum);
            }else if(isExt1and2SameSrc && !isExt1and2SameLabel){
                extEmb_2V_1src_1and2_notSameLabel_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,
                    recordPos1_dev,recordPos2_dev,auxArray_dev,baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,
                    embLen,isExt1Restrict,isExt2Restrict,isRecord1,isRecord2,useRecord1,useRecord2,true,maxRowNum);
            }else if(!isExt1and2SameSrc && isExt1and2SameLabel) {
                extEmb_2V_2src_1and2_sameLabel_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,
                    recordPos1_dev,auxArray_dev,baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,embLen,
                    isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isRecord1,useRecord1,true,maxRowNum);
            }else{
                extEmb_2V_2src_1and2_notSameLabel_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,
                    recordPos1_dev,recordPos2_dev,auxArray_dev,baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,
                    embLen,isExt1Restrict,isExt2Restrict,isRecord1,isRecord2,useRecord1,useRecord2,true,maxRowNum);
            }
            cudaMemcpy(writeRowNum,writeRowNum_dev,sizeof(uint)*2,cudaMemcpyDeviceToHost);
            tmppartialRowNum = writeRowNum[0];
        }
        if(partResults.size()>0){
            uint allocnum = ((tmppartialRowNum*embLen+31)>>5)<<5;
            uint *tmpaddr = (uint*)malloc(sizeof(uint)*allocnum);
            cudaMemcpyAsync(tmpaddr,baseAddr_dev+stride,sizeof(uint)*tmppartialRowNum*embLen,cudaMemcpyDeviceToHost);
            writeRowNum[0] = 0;
            writeRowNum[1] = 0;
            cudaMemcpyAsync(writeRowNum_dev,writeRowNum,sizeof(uint)*2,cudaMemcpyHostToDevice);
            tmppartResults.push_back(tmpaddr);
            tmppartResultRowNum.push_back(tmppartialRowNum);
            for(uint j=0;j<partResults.size();++j){
                tmpaddr = partResults[j];
                partialRowNum = partResultRowNum[j];
                partialEmb_dev = baseAddr_dev+gpuAvailSpace-partialRowNum*preEmbLen-32;
                cudaMemcpyAsync(partialEmb_dev,tmpaddr,sizeof(uint)*partialRowNum*preEmbLen,cudaMemcpyHostToDevice);
                maxRowNum = (gpuAvailSpace-stride-partialRowNum*preEmbLen-1500*38)/embLen;
                if(isExt1and2SameSrc && isExt1and2SameLabel) {
                    extEmb_2V_1src_1and2_sameLabel_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,
                        recordPos1_dev,auxArray_dev,baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,embLen,
                        isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isRecord1,useRecord1,false,maxRowNum);
                }else if(isExt1and2SameSrc && !isExt1and2SameLabel){
                    extEmb_2V_1src_1and2_notSameLabel_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,
                        recordPos1_dev,recordPos2_dev,auxArray_dev,baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,
                        embLen,isExt1Restrict,isExt2Restrict,isRecord1,isRecord2,useRecord1,useRecord2,false,maxRowNum);
                }else if(!isExt1and2SameSrc && isExt1and2SameLabel) {
                    extEmb_2V_2src_1and2_sameLabel_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,
                        recordPos1_dev,auxArray_dev,baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,embLen,
                        isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isRecord1,useRecord1,false,maxRowNum);
                }else{
                    extEmb_2V_2src_1and2_notSameLabel_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,
                        recordPos1_dev,recordPos2_dev,auxArray_dev,baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,
                        embLen,isExt1Restrict,isExt2Restrict,isRecord1,isRecord2,useRecord1,useRecord2,false,maxRowNum);
                }
                cudaMemcpy(writeRowNum,writeRowNum_dev,sizeof(uint)*2,cudaMemcpyDeviceToHost);
                free(tmpaddr);
                tmppartialRowNum = writeRowNum[0];
                while(tmppartialRowNum>=maxRowNum){
                    tmppartialRowNum = tmppartialRowNum-writeRowNum[1];
                    uint allocnum = ((tmppartialRowNum*embLen+31)>>5)<<5;
                    uint *tmpaddr = (uint*)malloc(sizeof(uint)*allocnum);
                    cudaMemcpyAsync(tmpaddr,baseAddr_dev+stride,sizeof(uint)*tmppartialRowNum*embLen,cudaMemcpyDeviceToHost);
                    writeRowNum[0] = 0;
                    writeRowNum[1] = 0;
                    cudaMemcpyAsync(writeRowNum_dev,writeRowNum,sizeof(uint)*2,cudaMemcpyHostToDevice);
                    tmppartResults.push_back(tmpaddr);
                    tmppartResultRowNum.push_back(tmppartialRowNum);
                    if(isExt1and2SameSrc && isExt1and2SameLabel) {
                        extEmb_2V_1src_1and2_sameLabel_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,
                            recordPos1_dev,auxArray_dev,baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,embLen,
                            isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isRecord1,useRecord1,true,maxRowNum);
                    }else if(isExt1and2SameSrc && !isExt1and2SameLabel){
                        extEmb_2V_1src_1and2_notSameLabel_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,
                            recordPos1_dev,recordPos2_dev,auxArray_dev,baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,
                            embLen,isExt1Restrict,isExt2Restrict,isRecord1,isRecord2,useRecord1,useRecord2,true,maxRowNum);
                    }else if(!isExt1and2SameSrc && isExt1and2SameLabel) {
                        extEmb_2V_2src_1and2_sameLabel_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,
                            recordPos1_dev,auxArray_dev,baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,embLen,
                            isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isRecord1,useRecord1,true,maxRowNum);
                    }else{
                        extEmb_2V_2src_1and2_notSameLabel_NoHash(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,
                            recordPos1_dev,recordPos2_dev,auxArray_dev,baseAddr_dev+stride,partialEmb_dev,intervalNum,partialRowNum,
                            embLen,isExt1Restrict,isExt2Restrict,isRecord1,isRecord2,useRecord1,useRecord2,true,maxRowNum);
                    }
                    cudaMemcpy(writeRowNum,writeRowNum_dev,sizeof(uint)*2,cudaMemcpyDeviceToHost);
                    tmppartialRowNum = writeRowNum[0];
                }
            }
        }
        partResults = tmppartResults;
        partResultRowNum = tmppartResultRowNum;


        if(isRecord1){
            repeatFlag[key1] = repeatFlag[key1]-1;
            uint *recordPos_host = recordPosData[key1].addr;
            cudaMemcpyAsync(recordPos_host,recordPos1_dev,sizeof(uint)*totSrcNum*2,cudaMemcpyDeviceToHost);
        }
        if(isRecord2){
            repeatFlag[key2] = repeatFlag[key2]-1;
            if(!isExt1and2SameLabel){
                uint *recordPos_host = recordPosData[key2].addr;
                cudaMemcpyAsync(recordPos_host,recordPos2_dev,sizeof(uint)*totSrcNum*2,cudaMemcpyDeviceToHost);
            }
        }
    }
    return partialEmb_dev;
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

uint* reduction(uint *edgelabelPartition,uint *edgeLabelPartition_dev, uint *vLabel_dev, uint*writeRowNum_dev, uint *writeRowNum,
               uint*baseAddr_dev,uint &stride,uint totSrcNum,uint intervalNum,Phase &curPhase,uint *partialEmb_dev,
               uint partialRowNum,bool isLastPhase){

    uint *labelPart_dev = edgeLabelPartition_dev+edgelabelPartition[2];
    uint *neighborsData_dev = edgeLabelPartition_dev+edgelabelPartition[4];
    uint edgelabel = curPhase.edgeLabel;
    uint embLen = posOfVisitedV.size();
    uint64_t key = edgelabel;
    key = (key<<32) | 0xffffffff;

    uint *recordpos = recordPosData[key].addr;
    uint recordsize = recordPosData[key].size;
    uint *recordpos_dev = baseAddr_dev+stride;
    stride = stride+recordsize;
    cudaMemcpyAsync(recordpos_dev,recordpos,sizeof(uint)*recordsize,cudaMemcpyHostToDevice);
    std::vector<uint> &pairs = curPhase.pairs;
    sortReducePairs(pairs);
    if(pairs.size()>32){
        std::cout<<"not process this situation now"<<std::endl;
        return NULL;
    }
    uint totUintNum = (pairs.size()+5)/6;
    std::vector<uint> auxArray(1+6+16+16+1,0);
    auxArray[0]=totUintNum;
    uint groupsize = 0,tmp=0,index=1;
    for(uint i=0;i<pairs.size()/2;++i){
        uint svid = pairs[i*2];
        uint evid = pairs[i*2+1];
        //std::cout<<"after sort: svid="<<svid<<" evid="<<evid<<std::endl;
        uint evidlabel = queryLabel[evid];
        auxArray[1+6+i] = evidlabel;
        key = edgelabel;
        key = (key<<32)|0xffffffff;
        uint *tmprecordpos = recordPosData[key].addr;
        uint *tmprecordpos_dev = recordpos_dev+(tmprecordpos-recordpos);
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

    uint *aux_dev = baseAddr_dev+stride;
    stride = stride+1+6+16+16+1;
    cudaMemcpyAsync(aux_dev,auxArray.data(),sizeof(uint)*(1+6+16+16+1),cudaMemcpyHostToDevice);
    //cudaDeviceSynchronize();
    //cudaError_t err = cudaGetLastError();
    //std::cout<<"0#### "<<cudaGetErrorString(err)<<"#####partialRowNum="<<partialRowNum<<std::endl;

    reductionPhase(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,recordpos_dev,aux_dev,partialEmb_dev,intervalNum,partialRowNum,embLen,isLastPhase);
    //cudaDeviceSynchronize();
    //err = cudaGetLastError();
    //std::cout<<"1#### "<<cudaGetErrorString(err)<<std::endl;
    std::vector<uint*> tmppartResults;
    std::vector<uint> tmppartResultRowNum;
    cudaMemcpy(writeRowNum,writeRowNum_dev,sizeof(uint)*2,cudaMemcpyDeviceToHost);
    partialRowNum = writeRowNum[0];
    if(partResults.size()>0){
        uint allocnum = ((partialRowNum*embLen+31)>>5)<<5;
        uint *tmpaddr = (uint*)malloc(sizeof(uint)*allocnum);
        cudaMemcpyAsync(tmpaddr,partialEmb_dev,sizeof(uint)*partialRowNum*embLen,cudaMemcpyDeviceToHost);
        tmppartResults.push_back(tmpaddr);
        tmppartResultRowNum.push_back(partialRowNum);
        writeRowNum[0] = 0;
        writeRowNum[1] = 0;
        cudaMemcpyAsync(writeRowNum_dev,writeRowNum,sizeof(uint)*2,cudaMemcpyHostToDevice);
        uint totpartrownum=0;
        maxRowNum = (gpuAvailSpace-stride)/embLen;
        for(uint j=0;j<partResults.size();++j){
            tmpaddr = partResults[j];
            partialRowNum = partResultRowNum[j];
            totpartrownum+=partialRowNum;
            if(totpartrownum>=maxRowNum){
                totpartrownum = totpartrownum-partialRowNum;
                partialEmb_dev = baseAddr_dev+gpuAvailSpace-totpartrownum*embLen;
                reductionPhase(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,recordpos_dev,aux_dev,partialEmb_dev,intervalNum,totpartrownum,embLen,isLastPhase);
                totpartrownum = partialRowNum;
                cudaMemcpy(writeRowNum,writeRowNum_dev,sizeof(uint)*2,cudaMemcpyDeviceToHost);
                partialRowNum = writeRowNum[0];
                allocnum = ((partialRowNum*embLen+31)>>5)<<5;
                tmpaddr = (uint*)malloc(sizeof(uint)*allocnum);
                tmppartResults.push_back(tmpaddr);
                tmppartResultRowNum.push_back(partialRowNum);
                cudaMemcpyAsync(tmpaddr,partialEmb_dev,sizeof(uint)*partialRowNum*embLen,cudaMemcpyDeviceToHost);
                partialEmb_dev = baseAddr_dev+gpuAvailSpace-totpartrownum*embLen;
                tmpaddr = partResults[j];
                cudaMemcpyAsync(partialEmb_dev,tmpaddr,sizeof(uint)*totpartrownum*embLen,cudaMemcpyHostToDevice);
                writeRowNum[0] = 0;
                writeRowNum[1] = 0;
                cudaMemcpyAsync(writeRowNum_dev,writeRowNum,sizeof(uint)*2,cudaMemcpyHostToDevice);
            }else{
                partialEmb_dev = baseAddr_dev+gpuAvailSpace-totpartrownum*embLen;
                cudaMemcpyAsync(partialEmb_dev,tmpaddr,sizeof(uint)*partialRowNum*embLen,cudaMemcpyHostToDevice);
            }
        }
        partialEmb_dev = baseAddr_dev+gpuAvailSpace-totpartrownum*embLen;
        reductionPhase(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,neighborsData_dev,recordpos_dev,aux_dev,partialEmb_dev,intervalNum,totpartrownum,embLen,isLastPhase);
        cudaMemcpy(writeRowNum,writeRowNum_dev,sizeof(uint)*2,cudaMemcpyDeviceToHost);
    }
    for(uint j=0;j<partResults.size();++j){
        uint *tmpaddr = partResults[j];
        free(tmpaddr);
    }
    partResultRowNum = tmppartResultRowNum;
    partResults = tmppartResults;
    cudaMemcpyAsync(recordpos,recordpos_dev,sizeof(uint)*recordsize,cudaMemcpyDeviceToHost);
    //cudaDeviceSynchronize();
    //err = cudaGetLastError();
    //std::cout<<"2#### "<<cudaGetErrorString(err)<<std::endl;
    return partialEmb_dev;
}



int main(int argc, char *argv[]){
    cudaError_t err;
    cudaDeviceGetAttribute(&GPU_SM_NUM, cudaDevAttrMultiProcessorCount, 0);
    std::string inputFileName = std::string(argv[1]);//suffix .mygraph
    std::string metaFileName = std::string(argv[2]);//suffix .metadata
    std::string queryFileName = std::string(argv[3]);//suffix .query

	if(argc==5){
    		std::string checkFileName = std::string(argv[4]);//suffix .g;
    		loadDataBaseAndCheck(inputFileName,metaFileName,checkFileName);
    		return 0;
	}

    std::cout<<"loading database"<<std::endl;
    allEdgeLabelPartitions = loadDataBase(inputFileName, metaFileName);

    //alloc gpu global memory, the unit is byte
    size_t total;
    uint *baseAddr_dev, stride=0, *partialEmb_dev;
    cudaMemGetInfo(&gpuAvailSpace,&total);
    //it seems that cudamemgetinfo triggers the bug
    gpuAvailSpace = ((gpuAvailSpace>>30)<<30);
    std::cout<<"gpu avail mem:"<<gpuAvailSpace<<" gpu tot mem:"<<total<<std::endl;
    cudaMalloc((void**)&baseAddr_dev,gpuAvailSpace);
    //now gpuAvailSpace is how many uints available
    maxRowNum=0;
    gpuAvailSpace = gpuAvailSpace>>2;
    //err = cudaGetLastError();
    //std::cout<<"0 "<<cudaGetErrorString(err)<<std::endl;
    //alloc gpu memory for vLabel
    vLabel = allEdgeLabelPartitions[0xffffffff];
    uint totVNum = vLabel[0];
    uint *vLabel_dev = baseAddr_dev;
    stride = stride + totVNum+1;
    cudaMemcpyAsync(vLabel_dev,vLabel,sizeof(uint)*(totVNum+1),cudaMemcpyHostToDevice);

    //alloc gpu memory for totWriteRowNum
    //writeRowNum[0] is totWriteRowNum, writeRowNum[1] is extraWriteRowNum
    uint writeRowNum[2];
    uint *writeRowNum_dev = baseAddr_dev+stride;
    stride = stride + 2;

    //vLabel_dev and totWriteRowNum_dev are fixed
    baseAddr_dev = baseAddr_dev + stride;
    gpuAvailSpace = gpuAvailSpace-stride;
    stride = 0;
    //generate match order, and transfer the first edge label partition to device once we find the edge label of first phase
    //stride = stride + len;
    clock_t start, end;
    start = clock();
    genMatchOrder(queryFileName,matchOrder,repeatFlag,restricts,edgeLabelStat,queryLabel,allEdgeLabelPartitions,baseAddr_dev,stride);
    end = clock();
    double endtimetmp=(double)(end-start)/CLOCKS_PER_SEC;
    std::cout<<"genmatchorder time: "<<endtimetmp*1000<<std::endl; 
    Phase &tmpPhase = matchOrder[0];
    /*for(auto ite=restricts.begin();ite!=restricts.end();++ite){
        uint key = ite->first;
        uint value = ite->second;
        std::cout<<"restrict: "<<key<<"<"<<value<<std::endl;
    }*/
    for(uint i=0;i<matchOrder.size();++i){
        Phase &phase = matchOrder[i];
        std::cout << "phase "<<i<<" type "<<phase.phaseType<<" edge label:"<<phase.edgeLabel<<std::endl;
        std::vector<uint> &pairs = phase.pairs;
        for(uint j=0;j<pairs.size()/2;++j){
            std::cout<<"svid:"<<pairs[j*2]<<" evid:"<<pairs[j*2+1]<<" ";
        }
        std::cout<<std::endl;
    }
    uint totWriteRowNum=0;
    uint label = tmpPhase.edgeLabel;
    uint *edgeLabelPartition = allEdgeLabelPartitions[label];
    uint len = edgeLabelPartition[0]-5;
    uint totSrcNum = edgeLabelPartition[1];
    uint intervalNum = edgeLabelPartition[2];
    uint *edgeLabelPartition_dev = baseAddr_dev+stride-len;
    uint embLen;
    if(matchOrder[0].pairs.size()==2){
        embLen=2;
    }else{
        embLen=3;
    }
    writeRowNum[0] = 0;
    writeRowNum[1] = 0;
    cudaMemcpyAsync(writeRowNum_dev,writeRowNum,sizeof(uint)*2,cudaMemcpyHostToDevice);
    //generate init embing
    start = clock();
    uint newstride = genInit(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,edgeLabelPartition_dev+edgeLabelPartition[4], 
        baseAddr_dev,stride,intervalNum,totSrcNum,false);
    //get tot write row number
    cudaMemcpy(writeRowNum,writeRowNum_dev,sizeof(uint)*2,cudaMemcpyDeviceToHost);
    uint partialRowNum = writeRowNum[0];
    end = clock();
    endtimetmp=(double)(end-start)/CLOCKS_PER_SEC;
    std::cout<<"init time and data copy: "<<endtimetmp*1000<<std::endl;
	std::cout<<"###partialRowNum="<<partialRowNum<<std::endl;
    //err = cudaGetLastError();
    //std::cout<<"2#### "<<cudaGetErrorString(err)<<std::endl;
    //std::cout<<"after first init, partialRowNum="<<partialRowNum<<" maxRowNum="<<maxRowNum<<std::endl;
    while(partialRowNum>=maxRowNum){
        //std::cout << "over size detect"<<std::endl;
        partialRowNum = partialRowNum-writeRowNum[1];
        totWriteRowNum += partialRowNum;
        uint allocnum = ((partialRowNum*embLen+31)>>5)<<5;
        uint *tmpaddr = (uint*)malloc(sizeof(uint)*allocnum);
        cudaMemcpyAsync(tmpaddr,baseAddr_dev+newstride,sizeof(uint)*partialRowNum*embLen,cudaMemcpyDeviceToHost);
        writeRowNum[0] = 0;
        writeRowNum[1] = 0;
        cudaMemcpyAsync(writeRowNum_dev,writeRowNum,sizeof(uint)*2,cudaMemcpyHostToDevice);
        partResults.push_back(tmpaddr);
        partResultRowNum.push_back(partialRowNum);
        //std::cout<<"begin oversize init"<<std::endl;
        newstride = genInit(vLabel_dev,writeRowNum_dev,edgeLabelPartition_dev,edgeLabelPartition_dev+edgeLabelPartition[4], 
                            baseAddr_dev,stride,intervalNum,totSrcNum,true);
        cudaMemcpy(writeRowNum,writeRowNum_dev,sizeof(uint)*2,cudaMemcpyDeviceToHost);
        partialRowNum = writeRowNum[0];
        //std::cout<<"after oversize init partialRowNUm="<<partialRowNum<<std::endl;
    }
    stride = newstride;
    
    //printPartialEmbGPU(baseAddr_dev+stride,partialRowNum,embLen,"partialEmbGPU_0.txt");
    writeRowNum[0] = writeRowNum[0]*embLen;
	start = clock();
    partialEmb_dev = baseAddr_dev+gpuAvailSpace-writeRowNum[0];
    copyPartialEmb(baseAddr_dev+stride,partialEmb_dev,writeRowNum[0]);
    //err = cudaGetLastError();
    //std::cout<<"2 "<<cudaGetErrorString(err)<<std::endl;
    //printPartialEmbGPU(partialEmb_dev,partialRowNum,embLen,"partialEmbGPU_tmp.txt");
    //initialize totWriteRowNum_dev to 0
    writeRowNum[0] = 0;
    writeRowNum[1] = 0;
    cudaMemcpyAsync(writeRowNum_dev,writeRowNum,sizeof(uint)*2,cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
	end = clock();
    endtimetmp=(double)(end-start)/CLOCKS_PER_SEC;
    std::cout<<"after init, transfer data time: "<<endtimetmp*1000<<std::endl;
    //err = cudaGetLastError();
    //std::cout<<"3 "<<cudaGetErrorString(err)<<std::endl;
    uint preedgelabel = label,i;
    for(i=1;i<matchOrder.size();++i){
	//std::cout<<"phase "<<i<<std::endl;
        bool isLastPhase = i==(matchOrder.size()-1);
        stride = 0;
        Phase &curPhase = matchOrder[i];
        if(curPhase.edgeLabel!=preedgelabel){
            preedgelabel = curPhase.edgeLabel;
            edgeLabelPartition = allEdgeLabelPartitions[preedgelabel];
            len = edgeLabelPartition[0]-5;
            totSrcNum = edgeLabelPartition[1];
            intervalNum = edgeLabelPartition[2];
		start = clock();
            cudaMemcpyAsync(edgeLabelPartition_dev,edgeLabelPartition+5,sizeof(uint)*len,cudaMemcpyHostToDevice);
	    cudaDeviceSynchronize();
		end = clock();
    	endtimetmp=(double)(end-start)/CLOCKS_PER_SEC;
    	std::cout<<"transfer label partition data time: "<<endtimetmp*1000<<std::endl;
        }
    	//cudaDeviceSynchronize();
        //err = cudaGetLastError();
        //std::cout<<"4 "<<i<<" "<<cudaGetErrorString(err)<<std::endl;
        stride = stride+len;
        if(curPhase.phaseType=='E'){
            if(curPhase.pairs.size()==2) { embLen += 1; }
            else { embLen += 2; }
            //uint distance = (partialEmb_dev-baseAddr_dev)-stride;
            //std::cout<<"in Extension row num:"<<partialRowNum<<" it can fit in: "<<distance/4<<" rows"<<std::endl;
       		start = clock(); 
            genExtEmb(vLabel_dev,writeRowNum_dev,writeRowNum,edgeLabelPartition_dev,edgeLabelPartition_dev+edgeLabelPartition[4],
                baseAddr_dev,stride,partialEmb_dev,intervalNum,totSrcNum,partialRowNum,curPhase,isLastPhase,false);
            partialRowNum = writeRowNum[0];
		end = clock();
     		endtimetmp=(double)(end-start)/CLOCKS_PER_SEC;
    		std::cout<<"query phase "<<i<<", genext time: "<<endtimetmp*1000<<std::endl;
		std::cout<<"###partialRowNum="<<partialRowNum<<std::endl;
		if(partialRowNum==0){ break; }
            //cudaDeviceSynchronize();
            //err = cudaGetLastError();
            //std::cout<<"7 "<<i<<" "<<cudaGetErrorString(err)<<" partialnum="<<partialRowNum<<std::endl;
            partialEmb_dev = baseAddr_dev+gpuAvailSpace-partialRowNum*embLen;
            //std::cout<<"test position:"<<partialRowNum<<" ="<<(partialEmb_dev-baseAddr_dev-stride)/embLen<<std::endl;
		start = clock();
            copyPartialEmb(baseAddr_dev+stride,partialEmb_dev,partialRowNum*embLen);
            //std::string filename = "partialEmbGPU_ext_"+std::to_string(i)+".txt";
            //printPartialEmbGPU(partialEmb_dev,partialRowNum,embLen,filename);
            //printPartialEmbGPUDebug(partialEmb_dev,partialRowNum,embLen,filename);
            //err = cudaGetLastError();
            //std::cout<<"9 "<<i<<" "<<cudaGetErrorString(err)<<" partialnum="<<partialRowNum<<std::endl;
            writeRowNum[0] = 0;
            writeRowNum[1] = 0;
            cudaMemcpyAsync(writeRowNum_dev,writeRowNum,sizeof(uint)*2,cudaMemcpyHostToDevice);
            cudaDeviceSynchronize();
		end=clock();
     		endtimetmp=(double)(end-start)/CLOCKS_PER_SEC;
    		std::cout<<"query phase "<<i<<", after genext, transfer time: "<<endtimetmp*1000<<std::endl;
            //err = cudaGetLastError();
            //std::cout<<"10 "<<i<<" "<<cudaGetErrorString(err)<<std::endl;
        }else{
            //std::cout<<"in Reduction"<<std::endl;
       		start = clock(); 
            partialEmb_dev = reduction(edgeLabelPartition,edgeLabelPartition_dev,vLabel_dev,writeRowNum_dev,writeRowNum,baseAddr_dev,stride,totSrcNum,intervalNum,curPhase,partialEmb_dev,partialRowNum,isLastPhase);
            partialRowNum = writeRowNum[0];
		end = clock();
     		endtimetmp=(double)(end-start)/CLOCKS_PER_SEC;
    		std::cout<<"query phase "<<i<<", reduction time: "<<endtimetmp*1000<<std::endl;
		std::cout<<"###partialRowNum="<<partialRowNum<<std::endl;
		if(partialRowNum==0){ break; }
            //cudaDeviceSynchronize();
            //err = cudaGetLastError();
            //std::cout<<"10 "<<i<<" "<<cudaGetErrorString(err)<<" partialnum="<<partialRowNum<<std::endl;
		//std::cout<<"reduction done:"<<partialRowNum<<std::endl;
       		start = clock(); 
            if(i<matchOrder.size()-1 && matchOrder[i+1].phaseType=='E'){
                uint *newpartialEmb_dev = baseAddr_dev+gpuAvailSpace-partialRowNum*embLen;
       		start = clock(); 
                copyPartialEmb(partialEmb_dev,newpartialEmb_dev,partialRowNum*embLen);
                partialEmb_dev = newpartialEmb_dev;
            	cudaDeviceSynchronize();
		end = clock();
     		endtimetmp=(double)(end-start)/CLOCKS_PER_SEC;
    		std::cout<<"transfer data  "<<endtimetmp*1000<<std::endl;
            }
            //std::cout<<"500 "<<i<<" "<<cudaGetErrorString(err)<<std::endl;
            //std::string filename = "partialEmbGPU_red_"+std::to_string(i)+".txt";
            //std::cout<<"rownum="<<partialRowNum<<std::endl;
            //printPartialEmbGPU(partialEmb_dev,partialRowNum,embLen,filename);
            //printPartialEmbGPUDebug(partialEmb_dev,partialRowNum,embLen,filename);
            writeRowNum[0] = 0;
            writeRowNum[1] = 0;
		start = clock();
            cudaMemcpyAsync(writeRowNum_dev,writeRowNum,sizeof(uint)*2,cudaMemcpyHostToDevice);
            cudaDeviceSynchronize();
		end = clock();
     		endtimetmp=(double)(end-start)/CLOCKS_PER_SEC;
    		std::cout<<"transfer data "<<endtimetmp*1000<<std::endl;
            //err = cudaGetLastError();
            //std::cout<<"11 "<<i<<" "<<cudaGetErrorString(err)<<" partialnum="<<partialRowNum<<std::endl;
        }
    }
    cudaDeviceSynchronize();
    end = clock();
    double endtime=(double)(end-start)/CLOCKS_PER_SEC;
    std::cout<<"time: "<<endtime*1000<<std::endl;
    err = cudaGetLastError();
    std::cout<<cudaGetErrorString(err)<<std::endl;
    std::cout<<"Result: "<<partialRowNum<<" extra size"<<partResults.size()<<std::endl;
    std::cout<<"PartialRowNum: "<<partialRowNum<<" embLen:  "<<embLen<<std::endl;
    printPartialEmbGPU(partialEmb_dev,partialRowNum,embLen,"partialResultsGPU.txt");
    return 0;
}
