#include <iostream>
#include "common.h"
extern int GPU_SM_NUM;
#define GETBASEADDR                                                                         \
    lowerLimit=0;                                                                           \
    predicate=0;                                                                            \
    len = (intervalNum+31)&0xffffffe0;                                                      \
    for(uint l=laneId;l<len;l=l+32){                                                             \
        lowerLimit = indexForIndex[l];                                                      \
        upperLimit = indexForIndex[l+intervalNum+1];                                        \
        upperLimit = lowerLimit+upperLimit-indexForIndex[l+intervalNum];                    \
        predicate = l<intervalNum?(svid>=lowerLimit && svid<upperLimit):0;                  \
        predicate = __ballot_sync(0xffffffff,predicate);                                    \
        if(predicate>0){                                                                    \
            index = svid-lowerLimit+indexForIndex[l+intervalNum];                           \
            uint tmpIndex = __ffs(predicate)-1;                                               \
            index = __shfl_sync(0xffffffff,index,tmpIndex);                                 \
            break;                                                                          \
        }                                                                                   \
    }                                                                                       \
    if(predicate==0){                                                                       \
        goto not_valid_nextEmb;                                                                             \
    }                                                                                       \
    baseNum = edgeLabelPartition[index];                                                    \
    baseAddr = neighborsData+baseNum;                                                       \
    baseLen = edgeLabelPartition[index+1]-baseNum;


#define GETLABEL_ADDR_LEN(pairIndex)                                                        \
    useRecord = (recordFlag>>pairIndex)&0x00000001;                                         \
    isRecord = (recordFlag>>pairIndex)&00010000;                                            \
    indexNum = 0;                                                                           \
    if(useRecord==1){                                                                       \
        uint distance = recordPos[index*2];                                                 \
        tmpNeigLen = recordPos[index*2+1];                                                  \
        if(distance==0 && tmpNeigLen==0){                                                   \
            lowerLimit = 0;                                                                 \
            upperLimit = baseLen;                                                       \
            FIND_LABEL_LIMIT(baseAddr, lowerLimit, upperLimit, evidlabel)                   \
            tmpNeigLen = upperLimit - lowerLimit;                                       \
            if (tmpNeigLen == 0) {                                                          \
                if(isRecord){                                                               \
                    if(laneId==0){ recordPos[index*2]=1; recordPos[index*2+1]=0; }          \
                }                                                                           \
                goto not_valid_nextEmb;                                                     \
            }                                                                               \
            if (tmpNeigLen > VBLOCKSIZE+2 && tmpNeigLen <= 32 * VBLOCKSIZE + 32){             \
                indexNum = (tmpNeigLen+VBLOCKSIZE) / (VBLOCKSIZE + 1);                            \
            }                                                                               \
            else if (tmpNeigLen > 32 * VBLOCKSIZE + 32) { indexNum= 32; }                   \
            tmpNeigAddr=baseAddr+lowerLimit;                                                \
            tmpNeigLen=tmpNeigLen-indexNum;                                                 \
            if (isRecord) {                                                                 \
                if (laneId == 0) {                                                          \
                    recordPos[index * 2] = lowerLimit;                                      \
                    recordPos[index * 2 + 1] = tmpNeigLen;                                  \
                }                                                                           \
            }                                                                               \
        }else if(distance==1&&tmpNeigLen==0){ goto not_valid_nextEmb; }                     \
    }else {                                                                                 \
        lowerLimit = 0; upperLimit = baseLen;                                           \
        FIND_LABEL_LIMIT(baseAddr, lowerLimit, upperLimit, evidlabel)                       \
        tmpNeigLen = upperLimit - lowerLimit;                                           \
        if(tmpNeigLen==0){                                                                  \
            if(isRecord){                                                                   \
                if(laneId==0){                                                              \
                    recordPos[index*2]=1;                                                   \
                    recordPos[index*2+1]=0;                                                 \
                }                                                                           \
            }                                                                               \
            goto not_valid_nextEmb;                                                         \
        }                                                                                   \
        if (tmpNeigLen > VBLOCKSIZE+2 && tmpNeigLen <= 32 * VBLOCKSIZE + 32) {                \
            indexNum = (tmpNeigLen+VBLOCKSIZE) / (VBLOCKSIZE + 1);                                \
        }                                                                                   \
        else if (tmpNeigLen > 32 * VBLOCKSIZE + 32) { indexNum = 32; }                      \
        tmpNeigAddr = baseAddr + lowerLimit;                                                \
        tmpNeigLen = tmpNeigLen - indexNum;                                                 \
        if (isRecord) {                                                                     \
            if (laneId == 0) {                                                              \
                recordPos[index * 2] = lowerLimit;                                          \
                recordPos[index * 2 + 1] = tmpNeigLen;                                      \
            }                                                                               \
        }                                                                                   \
    }


#define FINDEVID                                                                            \
    if(indexNum>0){                                                                         \
        upperLimit = laneId<indexNum?tmpNeigAddr[laneId]:0;                                 \
        lowerLimit = __shfl_up_sync(0xffffffff,upperLimit,1);                               \
        if(laneId==0) { lowerLimit = tmpNeigAddr[indexNum]-1; }                             \
        predicate = (evid>lowerLimit && evid<=upperLimit)?1:0;                              \
        predicate = __ballot_sync(0xffffffff,predicate);                                    \
        if(predicate==0){ goto not_valid_nextEmb; }                                         \
        predicate = __ffs(predicate)-1;                                                       \
        lowerLimit = predicate*blockSize;                                                   \
        upperLimit = predicate==indexNum-1?tmpNeigLen:lowerLimit+blockSize;                 \
    }                                                                                       \
    for(uint l=lowerLimit+laneId;l<((upperLimit+31)&0xffffffe0);l=l+32){                         \
        uint tmpV = l<upperLimit?tmpNeigAddr[l]:0;                                          \
        predicate = tmpV==evid;                                                             \
        predicate = __ballot_sync(0xffffffff,predicate);                                    \
        if(predicate>0){ break; }                                                           \
    }                                                                                       \
    if(predicate==0) { goto not_valid_nextEmb; }


//the format for indexForIndex[0:2*256] is vs1,vs2,vs3,len1,len2,len3. len1 is the numberr of vertices in interval 1, len2
//is the number of vertices in all previsous intervals (including interval 2).
//the format for indexForIndex[2*256:] is totUintNum(uint),(start1,end1),(start2,end2),labels(16,we do not need svid label),recordPos(16).
// first, all same start are grouped together, in each group
//all end with same label are grouped. each number uses 5 bits, each time we at most process 16 pairs. we use one uint to represent use
//use record and isrecord.
template<bool isLastPhase>
__global__ void reductionPhase_kernel(uint *vLabel, uint *totWriteRowNum, uint *edgeLabelPartition, uint *neighborsData,
    uint *baseRecordPos, uint *auxArray, uint *partialEmb, uint intervalNum, uint partialRowNum, uint embLen){

    uint laneId, gridWarpNum, warpIdInBlock, totUintNum, i, j;
    laneId = threadIdx.x & 31;
    gridWarpNum = (gridDim.x * blockDim.x)>>5;
    warpIdInBlock = threadIdx.x >> 5;
    __shared__ uint neigV[WARPPERBLOCK][33];
    __shared__ uint indexForIndex[256*2+1+1+64];
    for(i=threadIdx.x;i<intervalNum*2+1+1;i=i+blockDim.x){
        indexForIndex[i] = edgeLabelPartition[i];
    }
    for(i=threadIdx.x;i<64;i=i+blockDim.x){
        indexForIndex[256*2+1+1+i] = auxArray[i];
    }
    i = (blockIdx.x * blockDim.x + threadIdx.x)>>5;
    edgeLabelPartition = edgeLabelPartition + intervalNum*2+1+1;
    if(laneId==0) { neigV[warpIdInBlock][32] = 0; }
    __syncthreads();
    totUintNum = indexForIndex[256*2+2];
    uint recordFlag = indexForIndex[256*2+2+1+6+16+16];
    uint writePos = i;
    while (i < partialRowNum){
        if(laneId==0) { neigV[warpIdInBlock][32] = 0; }
        uint predicate,pairIndex,tmp,evidlabel,indexNum=0,svid,evid,presvid,preevid,preevidlabel;
        neigV[warpIdInBlock][laneId] = laneId<embLen?partialEmb[i * embLen + laneId]:1;
        predicate = neigV[warpIdInBlock][laneId]==0;
        predicate = __ballot_sync(0xffffffff,predicate);
        if(predicate>0){ goto zero_nextEmb; }
        pairIndex = 0;
        tmp = indexForIndex[256*2+2+1];
        evidlabel = indexForIndex[256*2+2+1+6];
        svid = (tmp>>5)&0x0000001f;
        evid = tmp&0x0000001f;
        svid = neigV[warpIdInBlock][svid];
        evid = neigV[warpIdInBlock][evid];
        presvid = svid;
        preevid = evid;
        preevidlabel = evidlabel;
        uint lowerLimit, upperLimit, index,baseLen,baseNum,len,useRecord,isRecord,found,lessThan,greatThan;
        uint *baseAddr,*tmpNeigAddr,*recordPos, tmpNeigLen,blockSize;
        recordPos = baseRecordPos+indexForIndex[256*2+2+1+6+16];
        GETBASEADDR
        GETLABEL_ADDR_LEN(0)

        blockSize = tmpNeigLen<=32*VBLOCKSIZE?VBLOCKSIZE:(tmpNeigLen>>5);
        lowerLimit = 0; upperLimit = tmpNeigLen;
        FINDEVID

        for(j=1;j<3;++j){
            pairIndex++;
            evidlabel = indexForIndex[256*2+2+1+6+pairIndex];
            evid = (tmp>>(j*2*5))&0x0000001f;
            svid = (tmp>>((j*2+1)*5))&0x0000001f;
            if(svid==0 || evid==0) { goto normal_nextEmb;}
            svid = neigV[warpIdInBlock][svid];
            evid = neigV[warpIdInBlock][evid];
            if(svid==presvid){
                if(evidlabel==preevidlabel){
                    FINDEVID
                }else{
                    preevidlabel = evidlabel;
                    recordPos = baseRecordPos+indexForIndex[256*2+2+1+6+16+pairIndex];
                    GETLABEL_ADDR_LEN(pairIndex)
                    uint blockSize = tmpNeigLen<=32*VBLOCKSIZE?VBLOCKSIZE:(tmpNeigLen>>5);
                    lowerLimit = 0; upperLimit = tmpNeigLen;
                    FINDEVID
                }
            }else{
                presvid = svid;
                preevidlabel = evidlabel;
                GETBASEADDR
                recordPos = baseRecordPos+indexForIndex[256*2+2+1+6+16+pairIndex];
                GETLABEL_ADDR_LEN(pairIndex)
                uint blockSize = tmpNeigLen<=32*VBLOCKSIZE?VBLOCKSIZE:(tmpNeigLen>>5);
                lowerLimit = 0; upperLimit = tmpNeigLen;
                FINDEVID
            }
        }

        for(j=1;j<totUintNum;++j){
            uint tmp = indexForIndex[256*2+2+1+j];
            for(uint k=0;k<3;++k){
                pairIndex++;
                evidlabel = indexForIndex[256*2+2+1+6+pairIndex];
                evid = (tmp>>(j*2*5))&0x0000001f;
                svid = (tmp>>((j*2+1)*5))&0x0000001f;
                if(svid==0 || evid==0) { goto normal_nextEmb;}
                svid = neigV[warpIdInBlock][svid];
                evid = neigV[warpIdInBlock][evid];
                if(svid==presvid){
                    if(evidlabel==preevidlabel){
                        FINDEVID
                    }else{
                        preevidlabel = evidlabel;
                        recordPos = baseRecordPos+indexForIndex[256*2+2+1+6+16+pairIndex];
                        GETLABEL_ADDR_LEN(pairIndex)
                        uint blockSize = tmpNeigLen<=32*VBLOCKSIZE?VBLOCKSIZE:(tmpNeigLen>>5);
                        lowerLimit = 0; upperLimit = tmpNeigLen;
                        FINDEVID
                    }
                }else{
                    presvid = svid;
                    preevidlabel = evidlabel;
                    GETBASEADDR
                    recordPos = baseRecordPos+indexForIndex[256*2+2+1+6+16+pairIndex];
                    GETLABEL_ADDR_LEN(pairIndex)
                    uint blockSize = tmpNeigLen<=32*VBLOCKSIZE?VBLOCKSIZE:(tmpNeigLen>>5);
                    lowerLimit = 0; upperLimit = tmpNeigLen;
                    FINDEVID
                }
            }
        }
        if(j==totUintNum) { goto normal_nextEmb; }
not_valid_nextEmb:
        if(!isLastPhase){
            if(laneId==0) { partialEmb[i*embLen] = 0; }
        }
zero_nextEmb:
        i = i+gridWarpNum;
        continue;
normal_nextEmb:
        if(!isLastPhase){
            if(writePos<i){
                if(laneId<embLen) { 
                    partialEmb[writePos*embLen+laneId] = neigV[warpIdInBlock][laneId]; 
                    partialEmb[i*embLen+laneId] = 0; 
                }
            }
            writePos = writePos+gridWarpNum;
            if(laneId==0) { neigV[warpIdInBlock][32] = 1; }
        }else{
            if(laneId==0) { writePos = atomicAdd(totWriteRowNum,1); }
            writePos = __shfl_sync(0xffffffff,writePos,0);
            if(laneId<embLen) { partialEmb[writePos*embLen+laneId] = neigV[warpIdInBlock][laneId]; }
        }
        i = i+gridWarpNum;
    }
    if(!isLastPhase){
        if(i>=gridWarpNum) { i=i-gridWarpNum; }
        if(laneId==0 && i<partialRowNum) {
            if(neigV[warpIdInBlock][32]==1){
                writePos = writePos-gridWarpNum;
            }
            writePos = writePos+1;
            atomicMax(totWriteRowNum,writePos); 
        }
    }
}

void reductionPhase(uint *vLabel, uint *totWriteRowNum, uint *edgeLabelPartition, uint *neighborsData, uint *baseRecordPos, uint *auxArray, 
    uint *partialEmb, uint intervalNum, uint partialRowNum, uint embLen, bool isLastPhase){

    int numBlocks;
    if(isLastPhase){
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocks,reductionPhase_kernel<true>,WARPPERBLOCK*32,0);
        reductionPhase_kernel<true><<<numBlocks*GPU_SM_NUM,WARPPERBLOCK*32>>>(vLabel,totWriteRowNum,edgeLabelPartition,neighborsData,baseRecordPos,auxArray,partialEmb,intervalNum,partialRowNum,embLen);
    }else{
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocks,reductionPhase_kernel<false>,WARPPERBLOCK*32,0);
        reductionPhase_kernel<false><<<numBlocks*GPU_SM_NUM,WARPPERBLOCK*32>>>(vLabel,totWriteRowNum,edgeLabelPartition,neighborsData,baseRecordPos,auxArray,partialEmb,intervalNum,partialRowNum,embLen);
    }

}

__global__ void copyPartialEmb_kernel(uint *src, uint *dst, uint numUints){
    uint index = blockIdx.x*blockDim.x+threadIdx.x;
    if(index<numUints){
        dst[numUints-1-index] = src[numUints-1-index];
    }
}

void copyPartialEmb(uint *src, uint *dst, uint numUints){
    uint blocksize = 8*32;
    uint blockNum = (numUints+blocksize-1)/blocksize;
    copyPartialEmb_kernel<<<blockNum,blocksize>>>(src, dst, numUints);

    /*uint blocksize = 8*32;
    uint blockNum;
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(&blockNum,copyPartialEmb_kernel,blocksize,0);
    copyPartialEmb_kernel<<<blockNum,blocksize>>>(src, dst, numUints);*/
}
