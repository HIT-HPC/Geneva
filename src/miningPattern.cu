#include <iostream>
#include "miningPattern.h"
#include "common.h"
//the number of uints
#define INITSHAREDMEMSIZE 128
//#define NBYNMINUSONE (2*256+31)

#define GENEMB_32BY32_1TRANS_RESTRICT(VID,EMBLEN)                                  \
    uint l,tmp1,tmp2,tmp3,tmp4;                                                    \
    uint writeVertex;                                                              \
    for(l=0;l<15*EMBLEN;++l){                                                       \
        uint index = (l<<5)+laneId;                                                \
        tmp1 = index/EMBLEN;                                                       \
        tmp2 = index-tmp1*EMBLEN;                                                  \
        tmp3 = tmp1>>5;                                                            \
        tmp4 = tmp1&31;                                                            \
        if(tmp4<=tmp3){ tmp3 = 30 - tmp3; tmp4 = 31 - tmp4; }                      \
        indexs[0] = __shfl_sync(0xffffffff,VID,tmp3);                              \
        indexs[1] = __shfl_sync(0xffffffff,VID,tmp4);                              \
        writeVertex = tmp2<EMBLEN-2?tmpmem_s[tmp2]:indexs[pos1Index];              \
        writeVertex = tmp2<EMBLEN-1?writeVertex:indexs[pos2Index];                 \
        newEmb_g[index] = writeVertex;                                             \
    }                                                                              \
    newEmb_g = newEmb_g + 15*32*EMBLEN;\
    for(l=laneId;l<((16*EMBLEN+31)&0xffffffe0);l=l+32){                                                       \
        tmp1 = l/EMBLEN;                                                       \
        tmp2 = l-tmp1*EMBLEN;                                                  \
        indexs[0] = __shfl_sync(0xffffffff,VID,15);                              \
        indexs[1] = __shfl_sync(0xffffffff,VID,tmp1+16);                              \
        writeVertex = tmp2<EMBLEN-2?tmpmem_s[tmp2]:indexs[pos1Index];              \
        writeVertex = tmp2<EMBLEN-1?writeVertex:indexs[pos2Index];                 \
        if(l<16*EMBLEN) { newEmb_g[l] = writeVertex; }                                             \
    }                                                                              \
    newEmb_g = newEmb_g + 16*EMBLEN;


#define GENEMB_NBYN_1TRANS_RESTRICT(N,VID,EMBLEN)                                  \
    uint l,tmp1,tmp2,tmp3,tmp4,writeVertex;                                        \
    uint tot = ((N*(N-1))>>1)*EMBLEN;                                          \
    uint tot_32 = (tot+31)&0xffffffe0;                                         \
    if((N&1)==1){                                                                  \
        for(uint l=laneId;l<tot_32;l=l+32){                                        \
            tmp1 = l/EMBLEN;                                                       \
            tmp2 = l-tmp1*EMBLEN;                                                  \
            tmp3 = tmp1/N;                                                         \
            tmp4 = tmp1-tmp3*N;                                                    \
            if(tmp3>=tmp4){ tmp3 = N-2-tmp3; tmp4 = N-1-tmp4; }                         \
            indexs[0] = __shfl_sync(0xffffffff,VID,tmp3);                          \
            indexs[1] = __shfl_sync(0xffffffff,VID,tmp4);                          \
            writeVertex = tmp2<EMBLEN-2?tmpmem_s[tmp2]:indexs[pos1Index];          \
            writeVertex = tmp2<EMBLEN-1?writeVertex:indexs[pos2Index];             \
            if(l<tot) { newEmb_g[l] = writeVertex; }                            \
        }                                                                          \
        newEmb_g = newEmb_g + tot;                                                 \
    }else{                                                                         \
        for(uint l=laneId;l<tot_32;l=l+32){                                        \
            tmp1 = l/EMBLEN;                                                       \
            tmp2 = l-tmp1*EMBLEN;                                                  \
            tmp3 = tmp1/(N-1);                                                         \
            tmp4 = tmp1-tmp3*(N-1)+1;                                                    \
            if(tmp4<=tmp3){ tmp3 = N-tmp3-1; tmp4 = N-tmp4; }                    \
            indexs[0] = __shfl_sync(0xffffffff,VID,tmp3);                          \
            indexs[1] = __shfl_sync(0xffffffff,VID,tmp4);                          \
            writeVertex = tmp2<EMBLEN-2?tmpmem_s[tmp2]:indexs[pos1Index];          \
            writeVertex = tmp2<EMBLEN-1?writeVertex:indexs[pos2Index];             \
            if(l<tot) { newEmb_g[l] = writeVertex; }                               \
        }                                                                          \
        newEmb_g = newEmb_g + tot;                                                 \
    }


#define GENEMB_NBYM_1OR2TRANS_NORESTRICT_EVALEQ(OUTERLEN,INNERLEN,OUTERVID,INNERVID,EMBLEN)     \
    uint tot = (OUTERLEN*INNERLEN*EMBLEN);                                                      \
    uint tot_32 = (tot+31)&0xffffffe0;                                                          \
    for(uint l=laneId;l<tot_32;l=l+32){                                                         \
        uint tmp1 = l/EMBLEN;                                                                   \
        uint tmp2 = l-tmp1*EMBLEN;                                                              \
        uint tmp3 = tmp1/INNERLEN;                                                              \
        uint tmp4 = tmp1-tmp3*INNERLEN;                                                         \
        indexs[0] = __shfl_sync(0xffffffff,OUTERVID,tmp3);                                      \
        indexs[1] = __shfl_sync(0xffffffff,INNERVID,tmp4);                                      \
        indexs[0] = indexs[0]==indexs[1]?0:indexs[0];                                           \
        uint writeVertex = tmp2<EMBLEN-2?tmpmem_s[tmp2]:indexs[pos1Index];                      \
        writeVertex = tmp2<EMBLEN-1?writeVertex:indexs[pos2Index];                              \
        if(l<tot) { newEmb_g[l] = writeVertex; }                                                \
    }                                                                                           \
    newEmb_g = newEmb_g + tot;


#define GENEMB_NBYM_2TRANS_NOEVALEQ(OUTERLEN,INNERLEN,OUTERVID,INNERVID,EMBLEN)    \
    uint tot = (OUTERLEN*INNERLEN*EMBLEN);                                         \
    uint tot_32 = (tot+31)&0xffffffe0;                                             \
    for(uint l=laneId;l<tot_32;l=l+32){                                            \
        uint tmp1 = l/EMBLEN;                                                      \
        uint tmp2 = l-tmp1*EMBLEN;                                                 \
        uint tmp3 = tmp1/INNERLEN;                                                 \
        uint tmp4 = tmp1-tmp3*INNERLEN;                                            \
        indexs[0] = __shfl_sync(0xffffffff,OUTERVID,tmp3);                         \
        indexs[1] = __shfl_sync(0xffffffff,INNERVID,tmp4);                         \
        uint writeVertex = tmp2<EMBLEN-2?tmpmem_s[tmp2]:indexs[pos1Index];         \
        writeVertex = tmp2<EMBLEN-1?writeVertex:indexs[pos2Index];                 \
        if(l<tot) { newEmb_g[l] = writeVertex; }                                   \
    }                                                                              \
    newEmb_g = newEmb_g + tot;


#define GENEMB_NBY32_2TRANS_NOEVALEQ(N,VID1,VID2,EMBLEN)                           \
    for(uint l=0;l<N;++l) {                                                        \
        indexs[0] = __shfl_sync(0xffffffff,VID1,l);                                \
        for(uint m=0;m<EMBLEN;++m) {                                               \
            uint index = (m<<5)+laneId;                                            \
            uint tmp1 = index/EMBLEN;                                              \
            uint tmp2 = index-tmp1*EMBLEN;                                         \
            indexs[1] = __shfl_sync(0xffffffff,VID2,tmp1);                         \
            uint writeVertex = tmp2<EMBLEN-2?tmpmem_s[tmp2]:indexs[pos1Index];     \
            writeVertex = tmp2<EMBLEN-1?writeVertex:indexs[pos2Index];             \
            newEmb_g[index] = writeVertex;                                         \
        }                                                                          \
        newEmb_g = newEmb_g+32*EMBLEN;                                             \
    }


#define GENEMB_32BYN_2TRANS_NOEVALEQ(N,VID1,VID2,EMBLEN)                           \
    for(uint l=0;l<N;++l) {                                                        \
        indexs[1] = __shfl_sync(0xffffffff,VID2,l);                                \
        for(uint m=0;m<EMBLEN;++m) {                                               \
            uint index = (m<<5)+laneId;                                            \
            uint tmp1 = index/EMBLEN;                                              \
            uint tmp2 = index-tmp1*EMBLEN;                                         \
            indexs[0] = __shfl_sync(0xffffffff,VID1,tmp1);                         \
            uint writeVertex = tmp2<EMBLEN-2?tmpmem_s[tmp2]:indexs[pos1Index];     \
            writeVertex = tmp2<EMBLEN-1?writeVertex:indexs[pos2Index];             \
            newEmb_g[index] = writeVertex;                                         \
        }                                                                          \
        newEmb_g = newEmb_g+(EMBLEN<<5);                                           \
    }


#define GENEMB_32BY32_2TRANS_NOEVALEQ(VID1,VID2,EMBLEN)                            \
    for(uint l=0;l<32;++l) {                                                       \
        indexs[0] = __shfl_sync(0xffffffff,VID1,l);                                \
        for(uint m=0;m<EMBLEN;++m) {                                               \
            uint index = (m<<5)+laneId;                                            \
            uint tmp1 = index/EMBLEN;                                              \
            uint tmp2 = index-tmp1*EMBLEN;                                         \
            indexs[1] = __shfl_sync(0xffffffff,VID2,tmp1);                         \
            uint writeVertex = tmp2<EMBLEN-2?tmpmem_s[tmp2]:indexs[pos1Index];     \
            writeVertex = tmp2<EMBLEN-1?writeVertex:indexs[pos2Index];             \
            newEmb_g[index] = writeVertex;                                         \
        }                                                                          \
        newEmb_g = newEmb_g+(EMBLEN<<5);                                           \
    }


#define GENEMB_32BY32_1TRANS_NORESTRICT_EVALEQ(VID,EMBLEN)                         \
    for(uint l=0;l<32;++l) {                                                       \
        indexs[0] = __shfl_sync(0xffffffff,VID,l);                                 \
        for(uint m=0;m<EMBLEN;++m) {                                               \
            uint index = (m<<5)+laneId;                                            \
            uint tmp1 = index/EMBLEN;                                              \
            uint tmp2 = index-tmp1*EMBLEN;                                         \
            indexs[1] = __shfl_sync(0xffffffff,VID,tmp1);                          \
            indexs[1] = indexs[1]==indexs[0]?0:indexs[1];                          \
            uint writeVertex = tmp2<EMBLEN-2?tmpmem_s[tmp2]:indexs[pos1Index];     \
            writeVertex = tmp2<EMBLEN-1?writeVertex:indexs[pos2Index];             \
            newEmb_g[index] = writeVertex;                                         \
        }                                                                          \
        newEmb_g = newEmb_g+(EMBLEN<<5);                                           \
    }


#define REARRANGE(EXTVID,BITFLAG,TOTCOUNT,EQUALS,EQUALE,SHAREDADDR)                           \
    BITFLAG = EXTVID==0?0:1;                                                                   \
    for(uint ii = EQUALS; ii < EQUALE; ++ii) {                                     \
        uint tmp = __shfl_sync(0xffffffff, equalVid, ii);                          \
        if(tmp==EXTVID){ BITFLAG = 0; EXTVID = 0; }                                \
    }                                                                              \
    BITFLAG = __ballot_sync(0xffffffff, BITFLAG );                                 \
    if(BITFLAG==0xffffffff){ TOTCOUNT = 32; }                                      \
    else{                                                                          \
        TOTCOUNT = __popc(BITFLAG);                                                \
        uint mask = 0xffffffff >> (31 - laneId);                                   \
        BITFLAG = mask & BITFLAG;                                                  \
        uint index = __popc(BITFLAG)-1;                                            \
        if(EXTVID>0){ tmpmem_s[SHAREDADDR+index] = EXTVID; }                               \
        EXTVID = tmpmem_s[SHAREDADDR+laneId];                                              \
    }


#define CAL_EQUAL(EXTVID,EQS,EQE)                                                  \
    for(uint ii=EQS;ii<EQE;++ii) {                                                 \
        uint tmp = __shfl_sync(0xffffffff,equalVid,ii);                            \
        EXTVID = tmp==EXTVID?0:EXTVID;                                             \
    }


#define CAL_EQUAL_NUM(EXTVID,COUNT,EQS,EQE)                                        \
    COUNT = 1;                                                                     \
    for(uint ii=EQS;ii<EQE;++ii) {                                                 \
        uint tmp = __shfl_sync(0xffffffff,equalVid,ii);                            \
        if(tmp==EXTVID) { COUNT = 0; EXTVID = 0; }                                 \
    }                                                                              \
    COUNT = __ballot_sync(0xffffffff, COUNT );


#define DETECT_PHASE(EXTVID,EQUALS,EQUALE,PHASEPOS,COUNTINDEX,LABELNUM)              \
    for(uint k=EQUALS;k<EQUALE;++k) {                                              \
        uint tmp = __shfl_sync(0xffffffff,equalVid,k);                             \
        if(tmp==EXTVID) {                                                          \
            uint tmppos = tmpmem_s[COUNTINDEX];                                    \
            tmpmem_s[COUNTINDEX] = tmppos+1;                                      \
            tmpmem_s[32+tmppos+LABELNUM]=PHASEPOS;                                 \
            EXTVID=0;                                                              \
        }                                                                          \
    }

#define GEN_1BY32(INNERVID)                                                         \
    for(uint l=laneId;l<32*embLen;l=l+32) {                                       \
        uint tmp1 = l/embLen;                                                      \
        uint tmp2 = l-tmp1*embLen;                                                 \
        indexs[1] = __shfl_sync(0xffffffff,INNERVID,tmp1);                           \
        uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[pos1Index];         \
        writeVertex = tmp2<embLen-1?writeVertex:indexs[pos2Index];                 \
        newEmb_g[l] = writeVertex;                                                 \
    }                                                                              \
    newEmb_g = newEmb_g + 32*embLen;



#define GEN_1BYN(N,INNERVID)                                                        \
    uint tmptot_32 = (N*embLen+31)&0xffffffe0;                                      \
    for(uint l=laneId;l<tmptot_32;l=l+32) {                                         \
        uint tmp1 = l/embLen;                                                       \
        uint tmp2 = l-tmp1*embLen;                                                  \
        indexs[1] = __shfl_sync(0xffffffff,INNERVID,tmp1);                          \
        uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[pos1Index];          \
        writeVertex = tmp2<embLen-1?writeVertex:indexs[pos2Index];                  \
        if(l<N*embLen) { newEmb_g[l] = writeVertex; }                               \
    }                                                                               \
    newEmb_g = newEmb_g + N*embLen;


#define ALLOCSPACE(N,EMBLEN)                                                       \
    if(laneId==0){ writePos = atomicAdd(totWriteRowNum,N); }                      \
    writePos = __shfl_sync(0xffffffff,writePos,0);                             \
    newEmb_g = basenewEmb_g+writePos*EMBLEN;                                   



__device__ uint* genInitEmbCore_3V_1and2_noOverlap(uint *neigData1_g, uint *neigData2_g, uint len1, uint len2,
           uint *newEmb_g, uint laneId, uint *tmpmem_s){
    const int embLen = 3;
    uint i,j, len1_32 = len1&0xffffffe0, len2_32 = len2&0xffffffe0, indexs[2];
    if(len1<len2){
        const uint pos1Index = 0, pos2Index = 1;
        for(i=laneId;i<len1_32;i=i+32){
            uint outerVid = neigData1_g[i];
            for(j=laneId;j<len2_32;j=j+32){ uint innerVid = neigData2_g[j]; GENEMB_32BY32_2TRANS_NOEVALEQ(outerVid,innerVid,3) }
            if(len2_32<len2){               uint innerVid = neigData2_g[j]; GENEMB_32BYN_2TRANS_NOEVALEQ((len2-len2_32),outerVid,innerVid,3) }
        }
        if(len1_32<len1){
            uint outerVid = neigData1_g[i];
            for(j=laneId;j<len2_32;j=j+32){ uint innerVid = neigData2_g[j]; GENEMB_NBY32_2TRANS_NOEVALEQ((len1-len1_32),outerVid,innerVid,3) }
            if(len2_32<len2){               uint innerVid = neigData2_g[j]; GENEMB_NBYM_2TRANS_NOEVALEQ((len1-len1_32),(len2-len2_32),outerVid,innerVid,embLen) }
        }
    }else{
        const uint pos1Index = 1, pos2Index = 0;
        for(i=laneId;i<len2_32;i=i+32){
            uint outerVid = neigData2_g[i];
            for(j=laneId;j<len1_32;j=j+32){ uint innerVid = neigData1_g[j]; GENEMB_32BY32_2TRANS_NOEVALEQ(outerVid,innerVid,3) }
            if(len1_32<len1){               uint innerVid = neigData1_g[j]; GENEMB_32BYN_2TRANS_NOEVALEQ((len2-len2_32),outerVid,innerVid,3) }
        }
        if(len2_32<len2){
            uint outerVid = neigData2_g[i];
            for(j=laneId;j<len1_32;j=j+32){ uint innerVid = neigData1_g[j]; GENEMB_NBY32_2TRANS_NOEVALEQ((len1-len1_32),outerVid,innerVid,3) }
            if(len1_32<len1){               
                uint innerVid = neigData1_g[j]; 
                GENEMB_NBYM_2TRANS_NOEVALEQ((len2-len2_32),(len1-len1_32),outerVid,innerVid,embLen)
            }
        }
    }
    return newEmb_g;
}

__device__ uint* genInitEmbCore_3V_1and2_fullOverlap_norestrict(uint *neigData_g, uint len, uint *newEmb_g, uint laneId, uint *tmpmem_s){
    const uint pos1Index = 0, pos2Index = 1;
    uint i,j, len_32 = len&0xffffffe0, indexs[2];
    for(i=laneId;i<len_32;i=i+32) {
        uint outerVid = neigData_g[i];
        GENEMB_32BY32_1TRANS_NORESTRICT_EVALEQ(outerVid,3)
        for(j=laneId;j<i;j=j+32) { uint innerVid = neigData_g[j]; GENEMB_32BY32_2TRANS_NOEVALEQ(outerVid,innerVid,3) }
        for(j=i+32;j<len_32;j=j+32) { uint innerVid = neigData_g[j]; GENEMB_32BY32_2TRANS_NOEVALEQ(outerVid,innerVid,3) }
        if(len_32<len) { uint innerVid = neigData_g[j]; GENEMB_32BYN_2TRANS_NOEVALEQ((len-len_32),outerVid,innerVid,3) }
    }
    if(len_32<len) {
        uint outerVid = neigData_g[i];
        GENEMB_NBYM_1OR2TRANS_NORESTRICT_EVALEQ((len-len_32),(len-len_32),outerVid,outerVid,3)
        for(j=laneId;j<len_32;j=j+32){ uint innerVid = neigData_g[j]; GENEMB_NBY32_2TRANS_NOEVALEQ((len-len_32),outerVid,innerVid,3) }
    }
    return newEmb_g;
}

__device__ uint* genInitEmbCore_3V_1and2_fullOverlap_restrict(uint *neigData_g, uint len, uint *newEmb_g, uint laneId, uint *tmpmem_s){
    uint i,j, len_32 = len&0xffffffe0, indexs[2];

    const uint pos1Index = 1, pos2Index = 0;
    for(i=laneId;i<len_32;i=i+32){
        uint outerVid = neigData_g[i];
        GENEMB_32BY32_1TRANS_RESTRICT(outerVid,3)
        for(j=i+32;j<len_32;j=j+32){ uint innerVid = neigData_g[j]; GENEMB_32BY32_2TRANS_NOEVALEQ(outerVid,innerVid,3) }
        if(len_32<len){ uint innerVid = neigData_g[j]; GENEMB_32BYN_2TRANS_NOEVALEQ((len-len_32),outerVid,innerVid,3) }
    }
    if(len_32<len){ 
        uint outerVid = neigData_g[i]; 
        GENEMB_NBYN_1TRANS_RESTRICT((len-len_32),outerVid,3)
    }
    return newEmb_g;
}

//this is modified
__device__ void genInitEmb_2V(uint *neigData_g,uint len,uint *basenewEmb_g, uint *totWriteRowNum, uint laneId, uint svid,
    uint *tmpmem_s,uint maxRowNum) {

    uint *newEmb_g, writePos;
    ALLOCSPACE(len,2)
    if(writePos+len>=maxRowNum){
        if(laneId==0) {tmpmem_s[0] = maxRowNum;}
        return;
    }
    uint i,j,index=0;
    uint len_32 = len&0xffffffe0;
    uint isodd = laneId&1, fetchIndex = laneId>>1;
    for(i=laneId;i<len_32;i=i+32){
        uint extvid = neigData_g[i];
        uint tmpvid1 = __shfl_sync(0xffffffff,extvid,fetchIndex);
        uint tmpvid2 = __shfl_sync(0xffffffff,extvid,fetchIndex+16);
        uint outerVid = isodd==1?tmpvid1:svid;
        uint innerVid = isodd==1?tmpvid2:svid;
        newEmb_g[index*64+laneId] = outerVid;
        newEmb_g[index*64+32+laneId] = innerVid;
        newEmb_g = newEmb_g + 64;
    }
    if(len_32<len){
        uint extvid = neigData_g[i];
        uint tot = (len-len_32)*2;
        uint tot_32 = (tot+31)&0xffffffe0;
        for(j=laneId;j<tot_32;j=j+32){
            uint tmp1 = j>>1;
            uint tmp2 = j&1;
            uint tmpvid = __shfl_sync(0xffffff,extvid,tmp1);
            uint extvid = tmp2==0?svid:tmpvid;
            if(j<tot){
                newEmb_g[j] = extvid;
            }
        }
    }
    if(laneId==0) { tmpmem_s[0] = 0; }
}

//this is modified
__device__ void genInitEmb_3V_1and2_restric(uint *neigData1_g, uint *neigData2_g, uint len1, uint len2,
           uint *tmpmem_s, uint *basenewEmb_g, uint *totWriteRowNum, uint laneId, uint svid,uint maxRowNum){
    uint *newEmb_g, totWrite=0,writePos;
    if(neigData1_g+len1<=neigData2_g) { if(laneId==0){ tmpmem_s[0] = 0; }return; }
    if(neigData2_g+len2<=neigData1_g){
        ALLOCSPACE((len1*len2),3)
        if(writePos+len1*len2>=maxRowNum){
            if(laneId==0){
                tmpmem_s[0] = maxRowNum;
                tmpmem_s[1] = len1*len2;
                return;
            }
        }
        genInitEmbCore_3V_1and2_noOverlap(neigData1_g,neigData2_g,len1,len2,newEmb_g,laneId,tmpmem_s);
        if(laneId==0){
            tmpmem_s[0]=0;
        }
        return;
    }
    if(neigData1_g < neigData2_g){ len1 = len1 - (neigData2_g-neigData1_g); neigData1_g = neigData2_g; }
    if(neigData1_g+len1<neigData2_g+len2){ len2 = len2-((neigData2_g+len2)-(neigData1_g+len1)); }
    uint tmpn = neigData1_g-neigData2_g;
    uint tmpm = (neigData1_g+len1)-(neigData2_g+len2);
    uint tmph = len1-tmpm;
    totWrite = (tmph*(2*tmpn+tmph-1))/2+tmpm*len2;
    ALLOCSPACE(totWrite,3)
    if(writePos+totWrite>=maxRowNum){
        if(laneId==0){
            tmpmem_s[0] = maxRowNum;
            tmpmem_s[1] = totWrite;
            return;
        }
    }
    if(tmpn>0){
        newEmb_g = genInitEmbCore_3V_1and2_noOverlap(neigData1_g,neigData2_g,len1,tmpn,newEmb_g,laneId,tmpmem_s);
    }
    if(tmpm>0){
        newEmb_g = genInitEmbCore_3V_1and2_noOverlap(neigData2_g+len2,neigData2_g,tmpm,len2,newEmb_g,laneId,tmpmem_s);
    }
    genInitEmbCore_3V_1and2_fullOverlap_restrict(neigData1_g,len1-tmpm,newEmb_g,laneId,tmpmem_s);
    if(laneId==0){
        tmpmem_s[0]=0;
    }
}

//this is modified
__device__ void genInitEmb_3V_1and2_noRestrict_sameLabel(uint *neigData1_g, uint *neigData2_g, uint len1, uint len2,
           uint *tmpmem_s, uint *basenewEmb_g, uint *totWriteRowNum, uint laneId,uint maxRowNum) {
    uint *newEmb_g, writePos;
    ALLOCSPACE((len1*len2),3)
    if(writePos+len1*len2>=maxRowNum){
        if(laneId==0){
            tmpmem_s[0] = maxRowNum;
            tmpmem_s[1] = len1*len2;
            return;
        }
    }
    if(neigData1_g+len1<=neigData2_g || neigData2_g+len2<=neigData1_g){
        genInitEmbCore_3V_1and2_noOverlap(neigData1_g,neigData2_g,len1,len2,newEmb_g,laneId,tmpmem_s);
        if(laneId==0){
            tmpmem_s[0]=0;
        }
        return;
    }
    if(neigData1_g<neigData2_g){
        uint tmplen = neigData2_g-neigData1_g;
        newEmb_g = genInitEmbCore_3V_1and2_noOverlap(neigData1_g,neigData2_g,tmplen,len2,newEmb_g,laneId,tmpmem_s);
        neigData1_g = neigData2_g;
        len1 = len1-tmplen;
    }
    if(neigData2_g<neigData1_g){
        uint tmplen = neigData1_g-neigData2_g;
        newEmb_g = genInitEmbCore_3V_1and2_noOverlap(neigData1_g,neigData2_g,len1,tmplen,newEmb_g,laneId,tmpmem_s);
        neigData2_g = neigData1_g;
        len2 = len2-tmplen;
    }
    if(len1<len2){
        uint tmplen = len2-len1;
        len2 = len2 - tmplen;
        uint *tmpaddr = neigData2_g+len2;
        newEmb_g = genInitEmbCore_3V_1and2_noOverlap(neigData1_g,tmpaddr,len1,tmplen,newEmb_g,laneId,tmpmem_s);
    }
    if(len2<len1){
        uint tmplen = len1-len2;
        len1 = len1 - tmplen;
        uint *tmpaddr = neigData1_g+len1;
        newEmb_g = genInitEmbCore_3V_1and2_noOverlap(tmpaddr,neigData2_g,tmplen,len2,newEmb_g,laneId,tmpmem_s);
    }
    genInitEmbCore_3V_1and2_fullOverlap_norestrict(neigData1_g,len1,newEmb_g,laneId,tmpmem_s);
    if(laneId==0){
        tmpmem_s[0] = 0;
    }
}

__device__ void genInitEmb_3V_1and2_notsamelabel(uint *neigData1_g, uint *neigData2_g, uint len1, uint len2,
           uint *tmpmem_s, uint *basenewEmb_g, uint *totWriteRowNum, uint laneId, uint maxRowNum) {
    uint *newEmb_g, writePos;
    ALLOCSPACE((len1*len2),3)
    if(writePos+len1*len2>=maxRowNum){
        if(laneId==0) { tmpmem_s[0] = maxRowNum; }
        return;
    }
    genInitEmbCore_3V_1and2_noOverlap(neigData1_g,neigData2_g,len1,len2,newEmb_g,laneId,tmpmem_s);
    if(laneId==0) { tmpmem_s[0] = 0; }
}

__device__ uint* genExtEmbCore_2V_1and2_sameLabel_noRestrict_fulloverlap_detectPhase(uint outerVidnew,uint *neigData_g,uint len,
           uint* tmpmem_s,uint equalVid,uint equalS, uint equalE,uint *newEmb_g,uint laneId,uint embLen) {
    uint j,indexs[2];
    indexs[0] = __shfl_sync(0xffffffff,outerVidnew,0);
    uint len_32 = len&0xffffffe0;
    for(j=laneId;j<len_32;j=j+32) {
        uint innerVidnew = neigData_g[j];
        DETECT_PHASE(innerVidnew,equalS,equalE,j,32,1)
        for(uint l=laneId;l<32*embLen;l=l+32) {                                       
            uint tmp1 = l/embLen;                                                    
            uint tmp2 = l-tmp1*embLen;                                               
            indexs[1] = __shfl_sync(0xffffffff,innerVidnew,tmp1);                       
            uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[0];
            writeVertex = tmp2<embLen-1?writeVertex:indexs[1]; 
            newEmb_g[l] = writeVertex;
            writeVertex = tmp2<embLen-2?writeVertex:indexs[1];
            writeVertex = tmp2<embLen-1?writeVertex:indexs[0]; 
            newEmb_g[l+32*embLen] = writeVertex;
        }                                                                          
        newEmb_g = newEmb_g + 32*embLen*2;
    }
    if(len_32<len){
        uint innerVidnew = laneId<len-len_32?neigData_g[j]:0;
        DETECT_PHASE(innerVidnew,equalS,equalE,j,32,1)
        uint tmptot = (len-len_32)*embLen;
        uint tmptot_32 = (tmptot+31)&0xffffffe0;                            
        for(uint l=laneId;l<tmptot_32;l=l+32) {                               
            uint tmp1 = l/embLen;                                             
            uint tmp2 = l-tmp1*embLen;                                        
            indexs[1] = __shfl_sync(0xffffffff,innerVidnew,tmp1);                
            uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[0];
            writeVertex = tmp2<embLen-1?writeVertex:indexs[1];
            if(l<tmptot) { newEmb_g[l] = writeVertex; }                     
            writeVertex = tmp2<embLen-2?writeVertex:indexs[1];
            writeVertex = tmp2<embLen-1?writeVertex:indexs[0];
            if(l<tmptot) { newEmb_g[l+tmptot] = writeVertex; } 
        }                                                                     
        newEmb_g = newEmb_g + tmptot*2;
    }
    return newEmb_g;
}

template<bool len1great2, bool evaluateEqual>
__device__ uint* genExtEmbCore_2V_1and2_sameLabel_detectPhase(uint outerVidnew,uint *neigData_g,uint len,
           uint* tmpmem_s,uint equalVid,uint equalS, uint equalE,uint *newEmb_g,uint laneId,uint embLen) {
    uint j,indexs[2];
    if(len1great2) { 
        const uint pos1Index = 1, pos2Index = 0;
        indexs[0] = __shfl_sync(0xffffffff,outerVidnew,0);
        uint len_32 = len&0xffffffe0;
        for(j=laneId;j<len_32;j=j+32) {
            uint innerVidnew = neigData_g[j];
            if(evaluateEqual) { innerVidnew = innerVidnew==indexs[0]?0:innerVidnew; }
            DETECT_PHASE(innerVidnew,equalS,equalE,j,32,1)
            //GEN_1BY32(innerVid)
            for(uint l=laneId;l<32*embLen;l=l+32) {                                       
                uint tmp1 = l/embLen;                                                    
                uint tmp2 = l-tmp1*embLen;                                               
                indexs[1] = __shfl_sync(0xffffffff,innerVidnew,tmp1);                       
                uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[pos1Index];       
                writeVertex = tmp2<embLen-1?writeVertex:indexs[pos2Index]; 

                newEmb_g[l] = writeVertex;                                              
            }                                                                          
            newEmb_g = newEmb_g + 32*embLen;
        }
        if(len_32<len){
            uint innerVidnew = laneId<len-len_32?neigData_g[j]:0;
            if(evaluateEqual) { innerVidnew = innerVidnew==indexs[0]?0:innerVidnew; }
            DETECT_PHASE(innerVidnew,equalS,equalE,j,32,1)
            //GEN_1BYN((len-len_32),innerVid)
            uint tmptot = (len-len_32)*embLen;
            uint tmptot_32 = (tmptot+31)&0xffffffe0;                            
            for(uint l=laneId;l<tmptot_32;l=l+32) {                               
                uint tmp1 = l/embLen;                                             
                uint tmp2 = l-tmp1*embLen;                                        
                indexs[1] = __shfl_sync(0xffffffff,innerVidnew,tmp1);                
                uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[pos1Index];
                writeVertex = tmp2<embLen-1?writeVertex:indexs[pos2Index];
     
                if(l<tmptot) { newEmb_g[l] = writeVertex; }                     
            }                                                                     
            newEmb_g = newEmb_g + tmptot;
        }
    }else { 
        const uint pos1Index = 0, pos2Index = 1; 
        indexs[0] = __shfl_sync(0xffffffff,outerVidnew,0);
        uint len_32 = len&0xffffffe0;
        for(j=laneId;j<len_32;j=j+32) {
            uint innerVidnew = neigData_g[j];
            if(evaluateEqual) { innerVidnew = innerVidnew==indexs[0]?0:innerVidnew; }
            DETECT_PHASE(innerVidnew,equalS,equalE,j,32,1)
            //GEN_1BY32(innerVid)
            for(uint l=laneId;l<32*embLen;l=l+32) {                                       
                uint tmp1 = l/embLen;                                                    
                uint tmp2 = l-tmp1*embLen;                                               
                indexs[1] = __shfl_sync(0xffffffff,innerVidnew,tmp1);                       
                uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[pos1Index];       
                writeVertex = tmp2<embLen-1?writeVertex:indexs[pos2Index]; 

                newEmb_g[l] = writeVertex;                                              
            }                                                                          
            newEmb_g = newEmb_g + 32*embLen;

        }
        if(len_32<len){
            uint innerVidnew = laneId<len-len_32?neigData_g[j]:0;
            if(evaluateEqual) { innerVidnew = innerVidnew==indexs[0]?0:innerVidnew; }
            DETECT_PHASE(innerVidnew,equalS,equalE,j,32,1)

            //GEN_1BYN((len-len_32),innerVid)
            uint tmptot = (len-len_32)*embLen;
            uint tmptot_32 = (tmptot+31)&0xffffffe0;                                      
            for(uint l=laneId;l<tmptot_32;l=l+32) {                                       
                uint tmp1 = l/embLen;                                                     
                uint tmp2 = l-tmp1*embLen;                                               
                indexs[1] = __shfl_sync(0xffffffff,innerVidnew,tmp1);                       
                uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[pos1Index];      
                writeVertex = tmp2<embLen-1?writeVertex:indexs[pos2Index];              

                if(l<tmptot) { newEmb_g[l] = writeVertex; }         
            }                                                       
            newEmb_g = newEmb_g + tmptot;
        }
    }
    return newEmb_g;
}

template<bool len1great2, bool evaluateEqual>
__device__ uint* genExtEmbCore_2V_1and2_noRestrict_sameLabel_usePhase(uint outerVid,uint *neigData_g,uint totcount, uint len,
           uint* tmpmem_s,uint phaseNum,uint *newEmb_g,uint laneId,uint embLen) {
    uint i,j,indexs[2];
    if(len1great2){
        const uint pos1Index = 1, pos2Index = 0;
        uint lowerLimit = 0;
        for(j=0;j<phaseNum;++j) {
            uint upperLimit = tmpmem_s[33+j];
            uint tmpLen = upperLimit-lowerLimit;
            uint tmpLen_32 = tmpLen&0xffffffe0;
            uint k;
            for(k=laneId;k<tmpLen_32;k=k+32) {
                uint innerVid = neigData_g[lowerLimit+k];
                for(i=0;i<totcount;++i) {
                    indexs[0] = __shfl_sync(0xffffffff,outerVid,i);
                    if(evaluateEqual) { innerVid = innerVid==indexs[0]?0:innerVid; }
                    //GEN_1BY32(innerVid)
                    for(uint l=laneId;l<32*embLen;l=l+32) {                                       
                        uint tmp1 = l/embLen;                                                    
                        uint tmp2 = l-tmp1*embLen;                                               
                        indexs[1] = __shfl_sync(0xffffffff,innerVid,tmp1);                      
                        uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[pos1Index];     
                        writeVertex = tmp2<embLen-1?writeVertex:indexs[pos2Index];             

                        newEmb_g[l] = writeVertex;                                             
                    }                                                                          
                    newEmb_g = newEmb_g + 32*embLen;

                }
            }
            if(tmpLen_32<tmpLen) {
                uint innerVid = k<tmpLen?neigData_g[lowerLimit+k]:0;
                if(evaluateEqual) {
                    //GENEMB_NBYM_1OR2TRANS_NORESTRICT_EVALEQ(totcount,(tmpLen-tmpLen_32),outerVid,innerVid,embLen)
                    uint tot = (totcount*(tmpLen-tmpLen_32)*embLen);                                             
                    uint tot_32 = (tot+31)&0xffffffe0;                                                         
                    for(uint l=laneId;l<tot_32;l=l+32){                                                      
                        uint tmp1 = l/embLen;                                                               
                        uint tmp2 = l-tmp1*embLen;                                                          
                        uint tmp3 = tmp1/(tmpLen-tmpLen_32);                                               
                        uint tmp4 = tmp1-tmp3*(tmpLen-tmpLen_32);                                          
                        indexs[0] = __shfl_sync(0xffffffff,outerVid,tmp3);                                
                        indexs[1] = __shfl_sync(0xffffffff,innerVid,tmp4);                                
                        indexs[0] = indexs[0]==indexs[1]?0:indexs[0];                                     
                        uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[pos1Index];                
                        writeVertex = tmp2<embLen-1?writeVertex:indexs[pos2Index];                        
                        if(l<tot) { newEmb_g[l] = writeVertex; }                                               
                    }                                                                                         
                    newEmb_g = newEmb_g + tot;
                }else {
                    GENEMB_NBYM_2TRANS_NOEVALEQ(totcount,(tmpLen-tmpLen_32),outerVid,innerVid,embLen)
                }
            }
            lowerLimit = upperLimit + 1;
        }
    }else { 
        const uint pos1Index = 0, pos2Index = 1; 
        uint lowerLimit = 0;
        for(j=0;j<phaseNum;++j) {
            uint upperLimit = tmpmem_s[33+j];
            uint tmpLen = upperLimit-lowerLimit;
            uint tmpLen_32 = tmpLen&0xffffffe0;
            uint k;
            for(k=laneId;k<tmpLen_32;k=k+32) {
                uint innerVid = neigData_g[lowerLimit+k];
                for(i=0;i<totcount;++i) {
                    indexs[0] = __shfl_sync(0xffffffff,outerVid,i);
                    if(evaluateEqual) { innerVid = innerVid==indexs[0]?0:innerVid; }
                    GEN_1BY32(innerVid)
                }
            }
            if(tmpLen_32<tmpLen) {
                uint innerVid = k<tmpLen?neigData_g[lowerLimit+k]:0;
                if(evaluateEqual) {
                    GENEMB_NBYM_1OR2TRANS_NORESTRICT_EVALEQ(totcount,(tmpLen-tmpLen_32),outerVid,innerVid,embLen)
                }else {
                    GENEMB_NBYM_2TRANS_NOEVALEQ(totcount,(tmpLen-tmpLen_32),outerVid,innerVid,embLen)
                }
            }
            lowerLimit = upperLimit + 1;
        }
    }
    return newEmb_g;
}


//invoked by extemb_1src_2L->genExtEmb_2V_1and2_notSameLabel_withPhase
//tmpmem_s needs to be 96
//this is modified
template<int labelNum, bool isContinue>
__device__ void genExtEmbCore_2V_1and2_noOverlap_withPhase(uint *neigData1_g,uint *neigData2_g,uint len1,uint len2,uint *tmpmem_s,
           uint equalVid, uint equalVNum1,uint equalVNum2, uint *totWriteRowNum, uint *basenewEmb_g, uint laneId, uint embLen,uint maxRowNum) {
    uint equalS1,equalE1,equalS2,equalE2, *newEmb_g,i,totcount,predicate,writePos;
    if(labelNum==1){ equalS1 = 0; equalS2 = 0; equalE1 = equalVNum1; equalE2 = equalVNum2; }
    else {equalS1 = 0; equalE1 = equalVNum1; equalS2 = equalVNum1; equalE2 = equalVNum2; }
    if(len1<len2) {
        if(!isContinue){
            i=laneId; totcount=0;
            uint outerVid, tot_32 = (len1+31)&0xffffffe0;
            while(totcount==0 && i<tot_32) {
                outerVid = i<len1?neigData1_g[i]:0;
                REARRANGE(outerVid,predicate,totcount,equalS1,equalE1,64)
                i=i+32;
            }
            if(totcount==0) { if(laneId==0) { tmpmem_s[66] = 0; } return; }
            tmpmem_s[32+laneId] = 0;
            ALLOCSPACE(len2,embLen)
            if(writePos+len2>=maxRowNum){
                if(laneId==0) { tmpmem_s[64] = 1; tmpmem_s[65] = len2; tmpmem_s[66] = maxRowNum; }
                return;
            }
            genExtEmbCore_2V_1and2_sameLabel_detectPhase<false,false>(outerVid,neigData2_g,len2,tmpmem_s,equalVid,equalS2,equalE2,newEmb_g,laneId,embLen);
            uint phaseNum = tmpmem_s[32];
            totcount = totcount - 1;
            if(totcount==0) { if(laneId==0) { tmpmem_s[66] = 0; } return; }
            outerVid = __shfl_down_sync(0xffffffff,outerVid,1);
            ALLOCSPACE(((len2-phaseNum)*totcount),embLen)
            if(tmpmem_s[32+1+phaseNum-1]!=len2-1) {
                if(laneId==0) { tmpmem_s[32+1+phaseNum] = len2; }
                phaseNum = phaseNum+1;
            }
            if(writePos+(len2-phaseNum)*totcount>=maxRowNum){
                if(laneId==0) { tmpmem_s[64] = 2; tmpmem_s[65] = (len2-phaseNum)*totcount; tmpmem_s[66]=maxRowNum; }
                return;
            }
            genExtEmbCore_2V_1and2_noRestrict_sameLabel_usePhase<false,false>(outerVid,neigData2_g,totcount,len2,tmpmem_s,phaseNum,newEmb_g,laneId,embLen);
            if(len1>32) {
                for(;i<tot_32;i=i+32) {
                    outerVid = i<len1?neigData1_g[i]:0;
                    REARRANGE(outerVid,predicate,totcount,equalS1,equalE1,64)
                    if(totcount==0) { continue; }
                    ALLOCSPACE(((len2-phaseNum)*totcount),embLen)
                    if(writePos+(len2-phaseNum)*totcount>=maxRowNum){
                        if(laneId==0) { tmpmem_s[64] = ((i>>5)<<16)|3; tmpmem_s[65] = (len2-phaseNum)*totcount; tmpmem_s[66] = maxRowNum; }
                        return;
                    }
                    genExtEmbCore_2V_1and2_noRestrict_sameLabel_usePhase<false,false>(outerVid,neigData2_g,totcount,len2,tmpmem_s,phaseNum,newEmb_g,laneId,embLen);
                }
            }
        }else{
            uint retval = tmpmem_s[64];
            uint loopindex = ((retval&0xffff0000)>>16)<<5;
            uint stopindex = retval&0x0000ffff;
            i=laneId; totcount=0;
            uint outerVid, tot_32 = (len1+31)&0xffffffe0;
            while(totcount==0 && i<tot_32) {
                outerVid = i<len1?neigData1_g[i]:0;
                REARRANGE(outerVid,predicate,totcount,equalS1,equalE1,64)
                i=i+32;
            }
            if(totcount==0) { if(laneId==0) { tmpmem_s[66] = 0; } return; }
            if(stopindex==1){
                tmpmem_s[32+laneId] = 0;
                ALLOCSPACE(len2,embLen)
                genExtEmbCore_2V_1and2_sameLabel_detectPhase<false,false>(outerVid,neigData2_g,len2,tmpmem_s,equalVid,equalS2,equalE2,newEmb_g,laneId,embLen);
            }
            uint phaseNum = tmpmem_s[32];
            if(stopindex<=2){
                totcount = totcount - 1;
                if(totcount==0) { if(laneId==0) { tmpmem_s[66] = 0; } return; }
                outerVid = __shfl_down_sync(0xffffffff,outerVid,1);
                ALLOCSPACE(((len2-phaseNum)*totcount),embLen)
                if(tmpmem_s[32+1+phaseNum-1]!=len2-1) {
                    if(laneId==0) { tmpmem_s[32+1+phaseNum] = len2; }
                    phaseNum = phaseNum+1;
                }
                genExtEmbCore_2V_1and2_noRestrict_sameLabel_usePhase<false,false>(outerVid,neigData2_g,totcount,len2,tmpmem_s,phaseNum,newEmb_g,laneId,embLen);
            }
            if(stopindex<=3){
                if(len1>32) {
                    if(stopindex==3){
                        i = loopindex+laneId;
                    }
                    for(;i<tot_32;i=i+32) {
                        outerVid = i<len1?neigData1_g[i]:0;
                        REARRANGE(outerVid,predicate,totcount,equalS1,equalE1,64)
                        if(totcount==0) { continue; }
                        ALLOCSPACE(((len2-phaseNum)*totcount),embLen)
                        genExtEmbCore_2V_1and2_noRestrict_sameLabel_usePhase<false,false>(outerVid,neigData2_g,totcount,len2,tmpmem_s,phaseNum,newEmb_g,laneId,embLen);
                    }
                }
            }
        }
    }else {
        if(!isContinue){
            i=laneId; totcount=0;
            uint outerVid, tot_32 = (len2+31)&0xffffffe0;
            while(totcount==0 && i<tot_32) {
                outerVid = i<len2?neigData2_g[i]:0;
                REARRANGE(outerVid,predicate,totcount,equalS2,equalE2,64)
                i=i+32;
            }
            if(totcount==0) { if(laneId==0) { tmpmem_s[66] = 0; } return; }
            tmpmem_s[32+laneId] = 0;
            ALLOCSPACE(len1,embLen)
            if(writePos+len1>=maxRowNum){
                if(laneId==0) { tmpmem_s[64] = 1; tmpmem_s[65] = len1; tmpmem_s[66]=maxRowNum; }
                return;
            }
            genExtEmbCore_2V_1and2_sameLabel_detectPhase<true,false>(outerVid,neigData1_g,len1,tmpmem_s,equalVid,equalS1,equalE1,newEmb_g,laneId,embLen);
            uint phaseNum = tmpmem_s[32];
            totcount = totcount - 1;
            if(totcount==0) { if(laneId==0) { tmpmem_s[66] = 0; } return; }
            outerVid = __shfl_down_sync(0xffffffff,outerVid,1);
            ALLOCSPACE((len1-phaseNum)*totcount,embLen)
            if(tmpmem_s[32+1+phaseNum-1]!=len1-1) {
                if(laneId==0) { tmpmem_s[32+1+phaseNum] = len1; }
                phaseNum = phaseNum+1;
                if(laneId==0) { tmpmem_s[32] = phaseNum; }
            }
            if(writePos+(len1-phaseNum)*totcount>=maxRowNum){
                if(laneId==0) { tmpmem_s[64] = 2; tmpmem_s[65] = (len1-phaseNum)*totcount; tmpmem_s[66] = maxRowNum; }
                return;
            }
            genExtEmbCore_2V_1and2_noRestrict_sameLabel_usePhase<true,false>(outerVid,neigData1_g,totcount,len1,tmpmem_s,phaseNum,newEmb_g,laneId,embLen);
            if(len2>32) {
                for(;i<tot_32;i=i+32) {
                    outerVid = i<len2?neigData2_g[i]:0;
                    REARRANGE(outerVid,predicate,totcount,equalS2,equalE2,64)
                    if(totcount==0) { continue; }
                    ALLOCSPACE(((len1-phaseNum)*totcount),embLen)
                    if(writePos+(len1-phaseNum)*totcount>=maxRowNum){
                        if(laneId==0) { tmpmem_s[64] = ((i>>5)<<16)|3; tmpmem_s[65] = (len1-phaseNum)*totcount; tmpmem_s[66]=maxRowNum; }
                        return;
                    }
                    genExtEmbCore_2V_1and2_noRestrict_sameLabel_usePhase<true,false>(outerVid,neigData1_g,totcount,len1,tmpmem_s,phaseNum,newEmb_g,laneId,embLen);
                }
            }
        }else{
            uint retval = tmpmem_s[64];
            uint loopindex = ((retval&0xffff0000)>>16)<<5;
            uint stopindex = retval&0x0000ffff;
            i=laneId; totcount=0;
            uint outerVid, tot_32 = (len2+31)&0xffffffe0;
            while(totcount==0 && i<tot_32) {
                outerVid = i<len2?neigData2_g[i]:0;
                REARRANGE(outerVid,predicate,totcount,equalS2,equalE2,64)
                i=i+32;
            }
            if(totcount==0) { if(laneId==0) { tmpmem_s[66] = 0; } return; }
            if(stopindex==1){
                tmpmem_s[32+laneId] = 0;
                ALLOCSPACE(len1,embLen)
                genExtEmbCore_2V_1and2_sameLabel_detectPhase<true,false>(outerVid,neigData1_g,len1,tmpmem_s,equalVid,equalS1,equalE1,newEmb_g,laneId,embLen);
            }
            uint phaseNum = tmpmem_s[32];
            if(stopindex<=2){
                totcount = totcount - 1;
                if(totcount==0) { if(laneId==0) { tmpmem_s[66] = 0; } return; }
                outerVid = __shfl_down_sync(0xffffffff,outerVid,1);
                ALLOCSPACE((len1-phaseNum)*totcount,embLen)
                if(tmpmem_s[32+1+phaseNum-1]!=len1-1) {
                    if(laneId==0) { tmpmem_s[32+1+phaseNum] = len1; }
                    phaseNum = phaseNum+1;
                    if(laneId==0) { tmpmem_s[32] = phaseNum; }
                }
                genExtEmbCore_2V_1and2_noRestrict_sameLabel_usePhase<true,false>(outerVid,neigData1_g,totcount,len1,tmpmem_s,phaseNum,newEmb_g,laneId,embLen);
            }
            if(stopindex<=3){
                if(len2>32) {
                    if(stopindex==3){
                        i=loopindex+laneId;
                    }
                    for(;i<tot_32;i=i+32) {
                        outerVid = i<len2?neigData2_g[i]:0;
                        REARRANGE(outerVid,predicate,totcount,equalS2,equalE2,64)
                        if(totcount==0) { continue; }
                        ALLOCSPACE(((len1-phaseNum)*totcount),embLen)
                        genExtEmbCore_2V_1and2_noRestrict_sameLabel_usePhase<true,false>(outerVid,neigData1_g,totcount,len1,tmpmem_s,phaseNum,newEmb_g,laneId,embLen);
                    }
                }
            }
        }
    }
}

//this need to be modified according to
//genExtEmbCore_2V_1and2_noRestrict_fullOverlap_withPhase
//this is modified
template<uint pos1Index, uint pos2Index, bool isContinue>
__device__ void genExtEmbCore_2V_1and2_restrict_fullOverlap_withPhase(uint *neigData_g,uint len,uint *tmpmem_s, uint equalVid,
           uint equalVNum, uint *totWriteRowNum, uint *basenewEmb_g, uint laneId, uint embLen, uint maxRowNum) {
    uint len_32 = len&0xffffffe0,*newEmb_g, tot_32 = (len+31)&0xffffffe0, indexs[2];
    uint i,j,predicate,totcount,outerVid,outerLowerLimit,writePos;

    i=laneId; totcount = 0;
    uint skipNum = 0;
    uint retval = tmpmem_s[64];
    while(totcount==0 && i<tot_32) {
        outerVid = i<len?neigData_g[i]:0;
        predicate = outerVid==0?0:1;                                                               
        for(uint ii = 0; ii < equalVNum; ++ii) {                                     
            uint tmp = __shfl_sync(0xffffffff, equalVid, ii);                        
            if(tmp==outerVid){ predicate = 0; outerVid = 0; }                        
        }                                                                           
        predicate = __ballot_sync(0xffffffff, predicate); 
        if(predicate==0xffffffff){ 
            totcount = 32; 
            break;
        }
        totcount = __popc(predicate);                                            
        if(totcount==0){
            i=i+32;
            continue;
        }
        skipNum = (i&0xffffffe0)+__ffs(predicate)-1;
        uint mask = 0xffffffff >> (31 - laneId);                               
        predicate = mask & predicate;                                              
        uint index = __popc(predicate)-1;                                        
        if(outerVid>0){ tmpmem_s[64+index] = outerVid; }                   
        outerVid = tmpmem_s[64+laneId];
        break;                  
    }
    if(totcount==0){ if(laneId==0){ tmpmem_s[66]=0; }return; }
    neigData_g = neigData_g+skipNum+1;
    len = len-skipNum-1;
    uint phaseNum;
    if(!isContinue){
        tmpmem_s[32+laneId] = 0;
        ALLOCSPACE(len,embLen)
        if(writePos+len>=maxRowNum){
            if(laneId==0) { tmpmem_s[64] = 1; tmpmem_s[65] = len; tmpmem_s[66] = maxRowNum; }
            return;
        }
        if(pos1Index>pos2Index){
            genExtEmbCore_2V_1and2_sameLabel_detectPhase<true,false>(outerVid, neigData_g, len,tmpmem_s, equalVid, 0, equalVNum, newEmb_g,laneId, embLen);
        }else{
            genExtEmbCore_2V_1and2_sameLabel_detectPhase<false,false>(outerVid, neigData_g, len,tmpmem_s, equalVid, 0, equalVNum, newEmb_g,laneId, embLen);
        }
        phaseNum = tmpmem_s[32];
        ALLOCSPACE((((len-phaseNum)*(len-phaseNum-1))>>1),embLen)
        if(writePos+(((len-phaseNum)*(len-phaseNum-1))>>1)>=maxRowNum){
            if(laneId==0) { tmpmem_s[64] = 2; tmpmem_s[65] = (((len-phaseNum)*(len-phaseNum-1))>>1); tmpmem_s[66] = maxRowNum; }
            return;
        }
        if(tmpmem_s[32+1+phaseNum-1]!=len-1) {
            if(laneId==0) { tmpmem_s[32+1+phaseNum] = len; }
            phaseNum = phaseNum+1;
        }
    }else{
        uint stopindex = retval & 0x0000000f;
        if(stopindex==1){
            tmpmem_s[32+laneId] = 0;
            ALLOCSPACE(len,embLen)
            if(pos1Index>pos2Index){
                genExtEmbCore_2V_1and2_sameLabel_detectPhase<true,false>(outerVid, neigData_g, len,tmpmem_s, equalVid, 0, equalVNum, newEmb_g,laneId, embLen);
            }else{
                genExtEmbCore_2V_1and2_sameLabel_detectPhase<false,false>(outerVid, neigData_g, len,tmpmem_s, equalVid, 0, equalVNum, newEmb_g,laneId, embLen);
            }
        }
        if(stopindex<=2){
            phaseNum = tmpmem_s[32];
            ALLOCSPACE((((len-phaseNum)*(len-phaseNum-1))>>1),embLen)
            if(tmpmem_s[32+1+phaseNum-1]!=len-1) {
                if(laneId==0) { tmpmem_s[32+1+phaseNum] = len; }
                phaseNum = phaseNum+1;
            }
        }
    }
    outerLowerLimit = 0;
    for(i=0;i<phaseNum;++i) {
        uint outerUpperLimit = tmpmem_s[33+i];
        uint outerLen = outerUpperLimit-outerLowerLimit;
        uint outerLen_32 = outerLen&0xffffffe0;
        for(j=laneId;j<outerLen_32;j=j+32) {
            outerVid = neigData_g[outerLowerLimit+j];
            GENEMB_32BY32_1TRANS_RESTRICT(outerVid,embLen)
            uint k;
            for(k=j+32;k<outerLen_32;k=k+32) {
                uint innerVid = neigData_g[outerLowerLimit+k];
                GENEMB_32BY32_2TRANS_NOEVALEQ(outerVid, innerVid, embLen)
            }
            if(outerLen_32<outerLen) {
                uint innerVid = neigData_g[outerLowerLimit+k];
                GENEMB_32BYN_2TRANS_NOEVALEQ((outerLen - outerLen_32), outerVid, innerVid, embLen)
            }
            uint innerLowerLimit = outerUpperLimit+1, innerUpperLimit;
            for(k=i+1;k<phaseNum;++k) {
                innerUpperLimit = tmpmem_s[33+k];
                uint innerLen = innerUpperLimit-innerLowerLimit;
                uint innerLen_32 = innerLen&0xffffffe0;
                uint g;
                for(g=laneId;g<innerLen_32;g=g+32) {
                    uint innerVid = neigData_g[innerLowerLimit+g];
                    GENEMB_32BY32_2TRANS_NOEVALEQ(outerVid, innerVid, embLen)
                }
                if(innerLen_32<innerLen) {
                    uint innerVid = neigData_g[innerLowerLimit+g];
                    GENEMB_32BYN_2TRANS_NOEVALEQ((innerLen - innerLen_32), outerVid, innerVid, embLen)
                }
                innerLowerLimit = innerUpperLimit + 1;
            }
        }
        if(outerLen_32<outerLen) {
            outerVid = neigData_g[outerLowerLimit+j];
            GENEMB_NBYN_1TRANS_RESTRICT((outerLen - outerLen_32), outerVid, embLen)
            uint innerLowerLimit = outerUpperLimit+1, innerUpperLimit;
            uint k;
            for(k=i+1;k<phaseNum;++k) {
                innerUpperLimit = tmpmem_s[33+k];
                uint innerLen = innerUpperLimit-innerLowerLimit;
                uint innerLen_32 = innerLen&0xffffffe0;
                uint g;
                for(g=laneId;j<innerLen_32;j=j+32) {
                    uint innerVid = neigData_g[innerLowerLimit+g];
                    GENEMB_NBY32_2TRANS_NOEVALEQ((outerLen - outerLen_32), outerVid, innerVid, embLen)
                }
                if(innerLen_32<innerLen) {
                    uint innerVid = neigData_g[innerLowerLimit+g];
                    GENEMB_NBYM_2TRANS_NOEVALEQ((outerLen-outerLen_32),(innerLen-innerLen_32),outerVid,innerVid,embLen)
                }
                innerLowerLimit = innerUpperLimit+1;
            }
        }
        outerLowerLimit = outerUpperLimit+1;
    }
    if(laneId==0){
        tmpmem_s[66] = 0;
    }
}

//this is modified
template<bool isContinue>
__device__ void genExtEmbCore_2V_1and2_noRestrict_fullOverlap_withPhase(uint *neigData_g,uint len,uint *tmpmem_s, uint equalVid,
           uint equalS,uint equalE, uint *totWriteRowNum, uint *basenewEmb_g, uint laneId, uint embLen,uint maxRowNum) {
    uint *newEmb_g, indexs[2];
    uint i,j,predicate, totcount,outerVid,outerLowerLimit,tot_32,writePos;
    const uint pos1Index = 0, pos2Index = 1;
    i=laneId; totcount = 0, tot_32 = (len+31)&0xffffffffe0;
    uint skipNum=0;
    uint retval = tmpmem_s[64];
    while(totcount==0 && i<tot_32) {
        outerVid = i<len?neigData_g[i]:0;
        predicate = outerVid==0?0:1;                                                         
        for(uint ii = equalS; ii < equalE; ++ii) {                                     
            uint tmp = __shfl_sync(0xffffffff, equalVid, ii);                        
            if(tmp==outerVid){ predicate = 0; outerVid = 0; }                        
        }                                                                           
        predicate = __ballot_sync(0xffffffff, predicate); 
        if(predicate==0xffffffff){ 
            totcount = 32; 
            break;
        }
        totcount = __popc(predicate);                                            
        if(totcount==0){
            i=i+32;
            continue;
        }
        skipNum = (i&0xffffffe0)+__ffs(predicate)-1;
        uint mask = 0xffffffff >> (31 - laneId);                               
        predicate = mask & predicate;                                              
        uint index = __popc(predicate)-1;                                        
        if(outerVid>0){ tmpmem_s[64+index] = outerVid; }                   
        outerVid = tmpmem_s[64+laneId];
        break;                  
    }
    if(totcount==0){ if(laneId==0) { tmpmem_s[66]=0; } return; }
    neigData_g = neigData_g + skipNum+1;
    len = len-skipNum-1;
    if(len==0) { if(laneId==0) { tmpmem_s[66]=0; } return; }
    uint phaseNum;
    if(!isContinue){
        tmpmem_s[32+laneId] = 0;
        ALLOCSPACE((len*2),embLen)
        if(writePos+len*2>=maxRowNum){
            if(laneId==0) { tmpmem_s[64]=1; tmpmem_s[65]=len*2; tmpmem_s[66]=maxRowNum; }
            return;
        }
        genExtEmbCore_2V_1and2_sameLabel_noRestrict_fulloverlap_detectPhase(outerVid,neigData_g,len,tmpmem_s,equalVid,equalS,equalE,newEmb_g,laneId,embLen);
        phaseNum = tmpmem_s[32];
        if(len-phaseNum-1==0) { if(laneId==0) { tmpmem_s[66]=0; } return; }
        ALLOCSPACE(((len-phaseNum-1)*(len-phaseNum)),embLen)
        if(writePos+((len-phaseNum-1)*(len-phaseNum))>=maxRowNum){
            if(laneId==0) { tmpmem_s[64]=2; tmpmem_s[65]=((len-phaseNum-1)*(len-phaseNum)); tmpmem_s[66]=maxRowNum; }
            return;
        }
        if(tmpmem_s[32+1+phaseNum-1]!=len-1) {
            if(laneId==0) { tmpmem_s[32+1+phaseNum] = len; }
            phaseNum = phaseNum+1;
        }
    }else{
        uint stopindex = retval&0x0000000f;
        if(stopindex==1){
            tmpmem_s[32+laneId] = 0;
            ALLOCSPACE((len*2),embLen)
            genExtEmbCore_2V_1and2_sameLabel_noRestrict_fulloverlap_detectPhase(outerVid,neigData_g,len,tmpmem_s,equalVid,equalS,equalE,newEmb_g,laneId,embLen);
        }
        if(stopindex<=2){
            phaseNum = tmpmem_s[32];
            ALLOCSPACE(((len-phaseNum-1)*(len-phaseNum)),embLen)
            if(tmpmem_s[32+1+phaseNum-1]!=len-1) {
                if(laneId==0) { tmpmem_s[32+1+phaseNum] = len; }
                phaseNum = phaseNum+1;
            }
        }
    }
    outerLowerLimit = 0;
    for(i=0;i<phaseNum;++i) {
        uint outerUpperLimit = tmpmem_s[33+i];
        uint outerLen = outerUpperLimit-outerLowerLimit;
        uint outerLen_32 = outerLen&0xffffffe0;
        for(j=laneId;j<outerLen_32;j=j+32) {
            outerVid = neigData_g[outerLowerLimit+j];
            for(uint l=laneId;l<32*31*embLen;l=l+32){
                uint tmp1 = l/embLen;
                uint tmp2 = l-tmp1*embLen;
                uint col = (tmp1>>5);
                uint row = tmp1-col*32;
                col = row<=col?col+1:col;
                row = __shfl_sync(0xffffffff,outerVid,row);                       
                col = __shfl_sync(0xffffffff,outerVid,col);
                uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:row;
                writeVertex = tmp2<embLen-1?writeVertex:col;
                newEmb_g[l] = writeVertex;
            }
            newEmb_g = newEmb_g + 32*31*embLen;
            for(uint k=j+32;k<outerLen_32;k=k+32){
                uint innerVid = neigData_g[outerLowerLimit+k];
                for(uint l=laneId;l<32*32*embLen;l=l+32){
                    uint tmp1 = l/embLen;
                    uint tmp2 = l-tmp1*embLen;
                    uint row = tmp1>>5;
                    uint col = tmp1-row*32;
                    row = __shfl_sync(0xffffffff,outerVid,row);                       
                    col = __shfl_sync(0xffffffff,innerVid,col);
                    uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:row;
                    writeVertex = tmp2<embLen-1?writeVertex:col;
                    newEmb_g[l] = writeVertex;
                    writeVertex = tmp2<embLen-2?writeVertex:col;
                    writeVertex = tmp2<embLen-1?writeVertex:row;
                    newEmb_g[32*32*embLen+l] = writeVertex;
                }
                newEmb_g = newEmb_g + 32*32*embLen*2;
            }
            if(outerLen_32<outerLen){
                uint innerVid = laneId<outerLen-outerLen_32?neigData_g[outerLowerLimit+outerLen_32+laneId]:0;
                uint tmptot = 32*(outerLen-outerLen_32)*embLen;
                for(uint l=laneId;l<tmptot;l=l+32){
                    uint tmp1 = l/embLen;
                    uint tmp2 = l-tmp1*embLen;
                    uint row = tmp1>>5;
                    uint col = tmp1-row*32;
                    row = __shfl_sync(0xffffffff,innerVid,row);                       
                    col = __shfl_sync(0xffffffff,outerVid,col);
                    uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:row;
                    writeVertex = tmp2<embLen-1?writeVertex:col;
                    newEmb_g[l] = writeVertex;
                    writeVertex = tmp2<embLen-2?writeVertex:col;
                    writeVertex = tmp2<embLen-1?writeVertex:row;
                    newEmb_g[tmptot+l] = writeVertex;
                }
                newEmb_g = newEmb_g + tmptot*2;
            }
            uint k, innerLowerLimit, innerUpperLimit;
            innerLowerLimit = outerUpperLimit+1;
            for(k=i+1;k<phaseNum;++k) {
                innerUpperLimit = tmpmem_s[33+k];
                uint innerLen = innerUpperLimit-innerLowerLimit;
                uint innerLen_32 = innerLen&0xffffffe0;
                for(uint g=laneId;g<innerLen_32;g=g+32){ 
                    uint innerVid = neigData_g[innerLowerLimit+g];
                    for(uint l=laneId;l<32*32*embLen;l=l+32){
                        uint tmp1 = l/embLen;
                        uint tmp2 = l-tmp1*embLen;
                        uint row = tmp1>>5;
                        uint col = tmp1-row*32;
                        row = __shfl_sync(0xffffffff,outerVid,row);                       
                        col = __shfl_sync(0xffffffff,innerVid,col);
                        uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:row;
                        writeVertex = tmp2<embLen-1?writeVertex:col;
                        newEmb_g[l] = writeVertex;
                        writeVertex = tmp2<embLen-2?writeVertex:col;
                        writeVertex = tmp2<embLen-1?writeVertex:row;
                        newEmb_g[32*32*embLen+l] = writeVertex;
                    }
                    newEmb_g = newEmb_g + 32*32*embLen*2;
                }
                if(innerLen_32<innerLen) { 
                    uint innerVid = laneId<innerLen-innerLen_32?neigData_g[innerLowerLimit+innerLen_32+laneId]:0; 
                    uint tmptot = 32*(innerLen-innerLen_32)*embLen;
                    for(uint l=laneId;l<tmptot;l=l+32){
                        uint tmp1 = l/embLen;
                        uint tmp2 = l-tmp1*embLen;
                        uint row = tmp1>>5;
                        uint col = tmp1-row*32;
                        row = __shfl_sync(0xffffffff,innerVid,row);                       
                        col = __shfl_sync(0xffffffff,outerVid,col);
                        uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:row;
                        writeVertex = tmp2<embLen-1?writeVertex:col;
                        newEmb_g[l] = writeVertex;
                        writeVertex = tmp2<embLen-2?writeVertex:col;
                        writeVertex = tmp2<embLen-1?writeVertex:row;
                        newEmb_g[tmptot+l] = writeVertex;
                    }
                    newEmb_g = newEmb_g + tmptot*2;
                }
                innerLowerLimit = innerUpperLimit + 1;
            }
        }
        if(outerLen_32<outerLen) {
            outerVid = laneId<outerLen-outerLen_32?neigData_g[outerLowerLimit+outerLen_32+laneId]:0;
            uint tmptot = (outerLen-outerLen_32)*(outerLen-outerLen_32-1)*embLen;
            uint tmptot_32 = (tmptot+31)&0xffffffe0;
            for(uint l=laneId;l<tmptot_32;l=l+32){
                uint tmp1 = l/embLen;
                uint tmp2 = l-tmp1*embLen;
                uint col = tmp1/(outerLen-outerLen_32);
                uint row = tmp1-col*(outerLen-outerLen_32);
                col = row<=col?col+1:col;
                row = __shfl_sync(0xffffffff,outerVid,row);                       
                col = __shfl_sync(0xffffffff,outerVid,col);
                uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:row;
                writeVertex = tmp2<embLen-1?writeVertex:col;
                if(l<tmptot) { newEmb_g[l] = writeVertex; }
            }
            newEmb_g = newEmb_g + tmptot;   
            uint k,innerLowerLimit = outerUpperLimit+1,innerUpperLimit;
            for(k=i+1;k<phaseNum;++k) {
                innerUpperLimit = tmpmem_s[33+k];
                uint innerLen = innerUpperLimit-innerLowerLimit;
                uint innerLen_32 = innerLen&0xffffffe0;
                for(uint g=laneId;g<innerLen_32;g=g+32) { 
                    uint innerVid = neigData_g[innerLowerLimit+g];
                    uint tmptot = 32*(outerLen-outerLen_32)*embLen;
                    for(uint l=laneId;l<tmptot;l=l+32){
                        uint tmp1 = l/embLen;
                        uint tmp2 = l-tmp1*embLen;
                        uint row = tmp1>>5;
                        uint col = tmp1-row*32;
                        row = __shfl_sync(0xffffffff,outerVid,row);                       
                        col = __shfl_sync(0xffffffff,innerVid,col);
                        uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:row;
                        writeVertex = tmp2<embLen-1?writeVertex:col;
                        newEmb_g[l] = writeVertex;
                        writeVertex = tmp2<embLen-2?writeVertex:col;
                        writeVertex = tmp2<embLen-1?writeVertex:row;
                        newEmb_g[tmptot+l] = writeVertex;
                    }
                    newEmb_g = newEmb_g + tmptot*2;
                }
                if(innerLen_32<innerLen) { 
                    uint innerVid = laneId<innerLen-innerLen_32?neigData_g[innerLowerLimit+innerLen_32+laneId]:0;
                    uint tmptot = (outerLen-outerLen_32)*(innerLen-innerLen_32)*embLen;
                    uint tmptot_32 = (tmptot+31)&0xffffffe0;
                    for(uint l=laneId;l<tmptot_32;l=l+32){
                        uint tmp1=l/embLen;
                        uint tmp2 = l-tmp1*embLen;
                        uint row = tmp1/(outerLen-outerLen_32);
                        uint col = tmp1-row*(outerLen-outerLen_32);
                        row = __shfl_sync(0xffffffff,innerVid,row);                       
                        col = __shfl_sync(0xffffffff,outerVid,col);
                        uint writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:row;
                        writeVertex = tmp2<embLen-1?writeVertex:col;
                        newEmb_g[l] = writeVertex;
                        writeVertex = tmp2<embLen-2?writeVertex:col;
                        writeVertex = tmp2<embLen-1?writeVertex:row;
                        newEmb_g[tmptot+l] = writeVertex;
                    }
                    newEmb_g = newEmb_g + tmptot*2;
                }
                innerLowerLimit = innerUpperLimit + 1;
            }
        }
        outerLowerLimit = outerUpperLimit+1;
    }
    if(laneId==0){ tmpmem_s[66] = 0; }
}

template<bool isLastPhase, bool isContinue>
__device__ void genExtEmb_1V(uint *neigData_g, uint len, uint *tmpmem_s, uint *equalVertices_s,
    uint *totWriteRowNum, uint *basenewEmb_g, uint laneId, uint embLen, uint maxRowNum) {

    uint *newEmb_g,i,writePos;
    uint equalVNum = equalVertices_s[0];
    uint equalVPos = equalVertices_s[2+laneId];
    uint equalVid = laneId<equalVNum?tmpmem_s[equalVPos]:0;

    if(!isLastPhase){
        ALLOCSPACE(len,embLen)
        if(writePos+len>=maxRowNum){
            if(laneId==0) { tmpmem_s[32] = 0xffff0000; tmpmem_s[33] = len; tmpmem_s[34] = maxRowNum; }
            return;
        }
        uint len_32 = len&0xffffffe0;
        for(i=laneId;i<len_32;i=i+32) {
            uint extvid = neigData_g[i];
            CAL_EQUAL(extvid,0,equalVNum)
            for(uint j=laneId;j<32*embLen;j=j+32) {
                uint tmp1 = j/embLen;
                uint tmp2 = j-tmp1*embLen;
                uint tmpvid = __shfl_sync(0xffffffff,extvid,tmp1);
                uint writeVertex = tmp2<embLen-1?tmpmem_s[tmp2]:tmpvid;
                newEmb_g[j] = writeVertex;
            }
            newEmb_g = newEmb_g+32*embLen;
        }
        if(len_32<len) {
            uint extvid = neigData_g[i];
            CAL_EQUAL(extvid,0,equalVNum)
            uint tmp = (len-len_32)*embLen;
            uint tmp_32 = (tmp+31)&0xffffffe0;
            for(uint j=laneId;j<tmp_32;j=j+32) {
                uint tmp1 = j/embLen;
                uint tmp2 = j-tmp1*embLen;
                uint tmpvid = __shfl_sync(0xffffffff,extvid,tmp1);
                uint writeVertex = tmp2<embLen-1?tmpmem_s[tmp2]:tmpvid;
                if(j<tmp) { newEmb_g[j] = writeVertex; }
            }
        }
    }else{
        uint len_32 = (len+31)&0xffffffe0, predicate,totcount;
        if(isContinue){
            uint retval = tmpmem_s[3];
            i=(retval>>16)+laneId;
        }else{
            i=laneId;
        }
        while(i<len_32) {
            totcount=0;
            uint extvid = i<len?neigData_g[i]:0;
            REARRANGE(extvid,predicate,totcount,0,equalVNum,32)
            ALLOCSPACE(totcount,embLen)
            if(writePos+totcount>=maxRowNum){
                if(laneId==0) { tmpmem_s[32] = ((i>>5)<<16)|0x00000001; tmpmem_s[33] = totcount; tmpmem_s[34] = maxRowNum;}
                return;
            }
            uint tmp = totcount*embLen;
            uint tmp_32 = (tmp+31)&0xffffffe0;
            for(uint j=laneId;j<tmp_32;j=j+32) {
                uint tmp1 = j/embLen;
                uint tmp2 = j-tmp1*embLen;
                uint tmpvid = __shfl_sync(0xffffffff,extvid,tmp1);
                uint writeVertex = tmp2<embLen-1?tmpmem_s[tmp2]:tmpvid;
                if(j<tmp) { newEmb_g[j] = writeVertex; }
            }
            newEmb_g = newEmb_g+tmp;
            i=i + 32;
        }
    }
    if(laneId==0){ tmpmem_s[34]=0; }
}

//invoked by extEmb_2V_1src_1and2_sameLabel_NoHash_kernel
//tmpmem_s 96
//this is modified
template<bool isContinue>
__device__ void genExtEmb_2V_1src_1and2_restrict_withPhase(uint *neigData1_g, uint *neigData2_g,
           uint len1, uint len2, uint *tmpmem_s, uint *equalVertices_s, uint *totWriteRowNum, 
           uint *basenewEmb_g, uint laneId, uint embLen, uint maxRowNum) {

    uint i, predicate,tmp;
    uint equalVNum = equalVertices_s[1];
    uint equalVPos = equalVertices_s[2+laneId];
    uint equalVid = laneId<equalVNum?tmpmem_s[equalVPos]:0;
    if(isContinue) { tmp = tmpmem_s[32+laneId]; }
    for(i=0;i<equalVNum;++i){
        equalVPos = __shfl_sync(0xffffffff,equalVid,i);
        predicate = equalVPos>equalVid;
        predicate = __ballot_sync(0xffffffff, predicate);
        predicate = __popc(predicate)-(32-equalVNum);
        if(laneId==0){
            tmpmem_s[32+predicate] = equalVPos;
        }
    }
    equalVid = laneId<equalVNum?tmpmem_s[32+laneId]:0;
    if(isContinue) { tmpmem_s[32+laneId] = tmp; }
    if(neigData1_g+len1<=neigData2_g) { return; }
    if(neigData2_g+len2<=neigData1_g) {
        genExtEmbCore_2V_1and2_noOverlap_withPhase<1,isContinue>(neigData1_g,neigData2_g,len1,len2,tmpmem_s,equalVid,equalVNum,equalVNum,totWriteRowNum,basenewEmb_g,laneId,embLen,maxRowNum);
        return;
    }
    if(!isContinue){
        if(neigData1_g<neigData2_g) { neigData1_g = neigData2_g; len1 = len1- (neigData2_g-neigData1_g); }
        else if(neigData1_g>neigData2_g){
            uint tmpLen = neigData1_g-neigData2_g;
            genExtEmbCore_2V_1and2_noOverlap_withPhase<1,false>(neigData1_g,neigData2_g,len1,tmpLen,tmpmem_s,equalVid,equalVNum,equalVNum,totWriteRowNum,basenewEmb_g,laneId,embLen,maxRowNum);
            if(tmpmem_s[66]==maxRowNum){
                tmp = tmpmem_s[64];
                tmp = (tmp<<4) | 0x00000001;
                tmp = tmp&0x0000ffff;
                if(laneId==0) { tmpmem_s[64] = (tmpmem_s[64]&0xffff0000) | tmp; }
                return;
            }
            neigData2_g = neigData1_g;
            len2 = len2 - tmpLen;
        }
        if(len1<len2) { len2 = len1; }
        else if(len1>len2){
            uint tmpLen = len1-len2;
            genExtEmbCore_2V_1and2_noOverlap_withPhase<1,false>(neigData2_g+len2,neigData2_g,tmpLen,len2,tmpmem_s,equalVid,equalVNum,equalVNum,totWriteRowNum,basenewEmb_g,laneId,embLen,maxRowNum);
            if(tmpmem_s[66]==maxRowNum){
                tmp = tmpmem_s[64];
                tmp = (tmp<<4) | 0x00000002;
                tmp = tmp&0x0000ffff;
                if(laneId==0) { tmpmem_s[64] = (tmpmem_s[64]&0xffff0000) | tmp; }
                return;
            }
            len1 = len1-tmpLen;
        }
        genExtEmbCore_2V_1and2_restrict_fullOverlap_withPhase<1,0,false>(neigData1_g,len1,tmpmem_s,equalVid,equalVNum,totWriteRowNum,basenewEmb_g,laneId,embLen,maxRowNum);
        if(tmpmem_s[66]==maxRowNum){
            tmp = tmpmem_s[64];
            tmp = (tmp<<4) | 0x00000003;
            tmp = tmp&0x0000ffff;
            if(laneId==0) { tmpmem_s[64] = (tmpmem_s[64]&0xffff0000) | tmp; }
            return;
        }
    }else{
        uint tmp = tmpmem_s[64+5];
        uint stopindex = tmp&0x0000000f;
        tmpmem_s[64+5] = ((tmp&0x0000ffff)>>4)|(tmp&0xffff0000);
        if(neigData1_g<neigData2_g) { neigData1_g = neigData2_g; len1 = len1- (neigData2_g-neigData1_g); }
        else if(neigData1_g>neigData2_g){
            uint tmpLen = neigData1_g-neigData2_g;
            if(stopindex==1){
                genExtEmbCore_2V_1and2_noOverlap_withPhase<1,true>(neigData1_g,neigData2_g,len1,tmpLen,tmpmem_s,equalVid,equalVNum,equalVNum,totWriteRowNum,basenewEmb_g,laneId,embLen,maxRowNum);
            }
            neigData2_g = neigData1_g;
            len2 = len2 - tmpLen;
        }
        if(len1<len2) { len2 = len1; }
        else if(len1>len2){
            uint tmpLen = len1-len2;
            if(stopindex<=2){
                genExtEmbCore_2V_1and2_noOverlap_withPhase<1,true>(neigData2_g+len2,neigData2_g,tmpLen,len2,tmpmem_s,equalVid,equalVNum,equalVNum,totWriteRowNum,basenewEmb_g,laneId,embLen,maxRowNum);
            }
            len1 = len1-tmpLen;
        }
        genExtEmbCore_2V_1and2_restrict_fullOverlap_withPhase<1,0,true>(neigData1_g,len1,tmpmem_s,equalVid,equalVNum,totWriteRowNum,basenewEmb_g,laneId,embLen,maxRowNum);
    }
}

//invoked by extEmb_2V_2src_1L, this uses 96 shared memory
//this is modified
template<bool len1SmallLen2,bool isContinue>
__device__ void genExtEmb_2V_2src_1and2_restrict(uint *neigData1_g, uint *neigData2_g,uint len1, uint len2, uint *tmpmem_s,
           uint *equalVertices_s, uint *totWriteRowNum, uint *newEmb_g,uint *basenewEmb_g, uint laneId, uint embLen,uint maxRowNum) {
    
    uint i,j,indexs[2],outerLen, innerLen, outerVid,predicate,writePos;
    uint equalVNum = equalVertices_s[0];
    uint equalVPos = equalVertices_s[1+laneId];
    uint equalVid = laneId<equalVNum?tmpmem_s[equalVPos]:0;
    for(i=0;i<equalVNum;++i){
        equalVPos = __shfl_sync(0xffffffff,equalVid,i);
        predicate = equalVPos>equalVid;
        predicate = __ballot_sync(0xffffffff, predicate);
        predicate = __popc(predicate)-(32-equalVNum);
        if(laneId==0){
            tmpmem_s[64+predicate] = equalVPos;
        }
    }
    equalVid = laneId<equalVNum?tmpmem_s[64+laneId]:0;
    //this is used to store rearranged outer vertices
    if(len1SmallLen2){ outerLen = len1; innerLen = len2; }
    else { outerLen = len2; innerLen = len1; }
    bool isfirst;
    if(isContinue){
        isfirst=true;
        i = tmpmem_s[32+5]+laneId;
    }else{
        i = laneId;
        tmpmem_s[32+laneId] = 0;
    }
    for(;i<((outerLen+31)&0xffffffe0);i=i+32) {
        if (len1SmallLen2) { outerVid = i < outerLen ? neigData1_g[outerLen - 1 - i] : 0; }
        else { outerVid = i < outerLen ? neigData2_g[i] : 0; }
        uint predicate,totcount;
        REARRANGE(outerVid,predicate,totcount,0,equalVNum,64)
        if(totcount == 0) { continue; }
        int boundryOfNoNeedEval = -1,newInnerLen = 0;
        bool notOverMaxRowNum=true;
        if(isContinue){
            if(isfirst){
                isfirst=false;
                outerVid = __shfl_down_sync(0xffffffff,outerVid,tmpmem_s[32+8]);
                totcount = totcount-tmpmem_s[32+8];
                j = tmpmem_s[32+6]+laneId;
                innerLen = tmpmem_s[32+7]; 
            }else{
                j=laneId;
            }
        }else{
            j = laneId;
        }
        for(;j<((innerLen+31)&0xffffffe0);j=j+32) {
            uint innerVid;
            if (len1SmallLen2) { innerVid = j<innerLen?neigData2_g[innerLen-1-j]:0; }
            else { innerVid = j<innerLen?neigData1_g[j]:0; }
            for(uint k=0;k<equalVNum;++k) { 
                uint tmp=__shfl_sync(0xffffffff,equalVid,k);
                innerVid = tmp==innerVid?0:innerVid; 
            }
            uint transsize = j<(innerLen&0xffffffe0)?32:innerLen-(innerLen&0xffffffe0);
            for(int k=0;k<=boundryOfNoNeedEval;++k) {
                indexs[0] = __shfl_sync(0xffffffff,outerVid,k);
                //tmpmem_s[64+k] is the distance between basenewEmb_g and the current pos of writed uints
                uint tmp = tmpmem_s[64+k];
                newEmb_g = basenewEmb_g + tmp;
                uint len = transsize*embLen;
                uint len_32 = len&0xffffffe0;
                uint l;
                for(l=laneId;l<len_32;l=l+32) {
                    uint tmp1 = l/embLen;
                    uint tmp2 = l-tmp1*embLen;
                    indexs[1] = __shfl_sync(0xffffffff,innerVid,tmp1);
                    uint writeVertex;
                    if(len1SmallLen2){
                        writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[0];
                        writeVertex = tmp2<embLen-1?writeVertex:indexs[1];
                    }else {
                        writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[1];
                        writeVertex = tmp2<embLen-1?writeVertex:indexs[0];
                    }
                    newEmb_g[l] = writeVertex;
                }
                if(len_32<len){
                    uint tmp1 = (len_32+laneId)/embLen;
                    uint tmp2 = (len_32+laneId)-tmp1*embLen;
                    indexs[1] = __shfl_sync(0xffffffff,innerVid,tmp1);
                    uint writeVertex;
                    if(len1SmallLen2){
                        writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[0];
                        writeVertex = tmp2<embLen-1?writeVertex:indexs[1];
                    }else {
                        writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[1];
                        writeVertex = tmp2<embLen-1?writeVertex:indexs[0];
                    }
                    if(laneId<len-len_32) { newEmb_g[l] = writeVertex; }
                }
                if(laneId==0) { tmp += transsize*embLen; tmpmem_s[64+k] = tmp; }
            }
            uint tmptot=0;
            newInnerLen = 0;
            int newboundryOfNoNeedEval = -1;
            if(notOverMaxRowNum){
                for(int k=boundryOfNoNeedEval+1;k<totcount;++k) {
                    indexs[0] = __shfl_sync(0xffffffff,outerVid,k);

                    if (len1SmallLen2) {
                        predicate = innerVid>0&&innerVid<indexs[0]?1:0;
                        predicate = __ballot_sync(0xffffffff,predicate);
                        if(predicate==0) { newInnerLen = innerLen-((j&0xffffffe0)+transsize); break; }
                        else {
                            predicate = __ffs(predicate)-1;
                            newboundryOfNoNeedEval = k;
                            tmptot += innerLen-(j&0xffffffe0)-predicate;
                            if(laneId==0) { tmpmem_s[32+k] = predicate; }
                        }
                    }else {
                        predicate = innerVid>0&&innerVid>indexs[0]?1:0;
                        predicate = __ballot_sync(0xffffffff,predicate);
                        if(predicate==0) { newInnerLen = innerLen-((j&0xffffffe0)+transsize); neigData1_g += ((j&0xffffffe0)+transsize);  break; }
                        else {
                            predicate = __ffs(predicate)-1;
                            newboundryOfNoNeedEval = k;
                            tmptot += innerLen-(j&0xffffffe0)-predicate;
                            if(laneId==0) { tmpmem_s[32+k] = predicate; }
                        }
                    }
                }
                if(newboundryOfNoNeedEval>=0) {
                    ALLOCSPACE(tmptot,embLen)
                    if(writePos+tmptot>=maxRowNum){
                        notOverMaxRowNum = false;
                        if(laneId==0){
                            tmpmem_s[32] = i>>5;
                            tmpmem_s[33] = j>>5;
                            tmpmem_s[34] = maxRowNum;
                            tmpmem_s[35] = innerLen;
                            tmpmem_s[36] = (boundryOfNoNeedEval+1);
                            tmpmem_s[37] = tmptot;
                        }
                        if(boundryOfNoNeedEval==-1){
                            tmpmem_s[64+laneId] = tmpmem_s[32+laneId];
                            return;
                        }
                    }else{
                        tmptot = newEmb_g-basenewEmb_g;
                        for(int k=boundryOfNoNeedEval+1;k<=newboundryOfNoNeedEval;++k) {
                            indexs[0] = __shfl_sync(0xffffffff,outerVid,k);
                            predicate = tmpmem_s[32+k];
                            uint len = (transsize-predicate)*embLen;
                            uint len_32 = len&0xffffffe0;
                            uint l;
                            for(l=laneId;l<len_32;l=l+32) {
                                uint tmp1 = l/embLen;
                                uint tmp2 = l-tmp1*embLen;
                                indexs[1] = __shfl_sync(0xffffffff,innerVid,tmp1);
                                uint writeVertex;
                                if(len1SmallLen2){
                                    writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[0];
                                    writeVertex = tmp2<embLen-1?writeVertex:indexs[1];
                                }else {
                                    writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[1];
                                    writeVertex = tmp2<embLen-1?writeVertex:indexs[0];
                                }
                                newEmb_g[l] = writeVertex;
                            }
                            if(len_32<len){
                                uint tmp1 = (len_32+laneId)/embLen;
                                uint tmp2 = (len_32+laneId)-tmp1*embLen;
                                indexs[1] = __shfl_sync(0xffffffff,innerVid,tmp1);
                                uint writeVertex;
                                if(len1SmallLen2){
                                    writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[0];
                                    writeVertex = tmp2<embLen-1?writeVertex:indexs[1];
                                }else {
                                    writeVertex = tmp2<embLen-2?tmpmem_s[tmp2]:indexs[1];
                                    writeVertex = tmp2<embLen-1?writeVertex:indexs[0];
                                }
                                if(laneId<len-len_32) { newEmb_g[l] = writeVertex; }
                            }
                            if(laneId==0) { tmpmem_s[64+k] = tmptot+(transsize-predicate)*embLen; }
                            tmptot = tmptot+(transsize-predicate+innerLen-(j&0xffffffe0))*embLen;
                            newEmb_g = basenewEmb_g + tmptot;
                        }
                        boundryOfNoNeedEval = newboundryOfNoNeedEval;
                    }
                }
            }
        }
        if(!notOverMaxRowNum){
            tmpmem_s[64+laneId] = tmpmem_s[32+laneId];
            return;
        }
        innerLen = newInnerLen>0?newInnerLen:innerLen;
    }
    if(laneId==0) { tmpmem_s[66] = 0; }
}

//invoked by extEmb_2V_1src_1and2_sameLabel_NoHash_kernel
//tmpmem_s 96
//this is modified
template <bool isContinue>
__device__ void genExtEmb_2V_1src_1and2_noRestrict_sameLabel_evalEq_withPhase(uint *neigData1_g, uint *neigData2_g,
    uint len1, uint len2, uint *tmpmem_s, uint *equalVertices_s, uint *totWriteRowNum, uint *basenewEmb_g, 
    uint laneId, uint embLen, uint maxRowNum) {

    uint *newEmb_g, i, predicate,tmp;
    uint equalVNum = equalVertices_s[1];
    uint equalVPos = equalVertices_s[2+laneId];
    uint equalVid = laneId<equalVNum?tmpmem_s[equalVPos]:0;
    if(isContinue){ tmp = tmpmem_s[32+laneId]; }
    for(i=0;i<equalVNum;++i){
        equalVPos = __shfl_sync(0xffffffff,equalVid,i);
        predicate = equalVPos>equalVid;
        predicate = __ballot_sync(0xffffffff, predicate);
        predicate = __popc(predicate)-(32-equalVNum);
        if(laneId==0){
            tmpmem_s[32+predicate] = equalVPos;
        }
    }
    equalVid = laneId<equalVNum?tmpmem_s[32+laneId]:0;
    if(isContinue){ tmpmem_s[32+laneId] = tmp; }
    if(neigData1_g+len1<=neigData2_g || neigData2_g+len2<=neigData1_g) {
        genExtEmbCore_2V_1and2_noOverlap_withPhase<1,isContinue>(neigData1_g, neigData2_g, len1, len2, tmpmem_s, equalVid, equalVNum,equalVNum, totWriteRowNum, basenewEmb_g, laneId, embLen,maxRowNum);
        return;
    }

    if(!isContinue){
        if(neigData1_g<neigData2_g) {
            uint tmpLen = neigData2_g-neigData1_g;
            genExtEmbCore_2V_1and2_noOverlap_withPhase<1,false>(neigData1_g, neigData2_g, tmpLen, len2, tmpmem_s, equalVid, equalVNum,equalVNum, totWriteRowNum, basenewEmb_g, laneId, embLen,maxRowNum);
            if(tmpmem_s[66] == maxRowNum){
                tmp = tmpmem_s[64];
                tmp = ((tmp<<4) | 1)&0x0000ffff;
                if(laneId==0) { tmpmem_s[64] = (tmpmem_s[64]&0xffff0000)|tmp; }
                return;
            }
            neigData1_g = neigData2_g;
            len1 = len1 - tmpLen;
        }else if(neigData2_g<neigData1_g) {
            uint tmpLen = neigData1_g-neigData2_g;
            genExtEmbCore_2V_1and2_noOverlap_withPhase<1,false>(neigData1_g, neigData2_g, len1, tmpLen, tmpmem_s, equalVid, equalVNum,equalVNum, totWriteRowNum, basenewEmb_g, laneId, embLen,maxRowNum);
            if(tmpmem_s[66] == maxRowNum){
                tmp = tmpmem_s[64];
                tmp = ((tmp<<4) | 2)&0x0000ffff;
                if(laneId==0) { tmpmem_s[64] = (tmpmem_s[64]&0xffff0000)|tmp; }
                return;
            }
            neigData2_g = neigData1_g;
            len2 = len2 - tmpLen;
        }
        if(len1<len2) {
            uint tmpLen = len2-len1;
            genExtEmbCore_2V_1and2_noOverlap_withPhase<1,false>(neigData1_g, neigData2_g+len1, len1, tmpLen, tmpmem_s, equalVid, equalVNum, equalVNum, totWriteRowNum, basenewEmb_g, laneId, embLen,maxRowNum);
            if(tmpmem_s[66] == maxRowNum){
                tmp = tmpmem_s[64];
                tmp = ((tmp<<4) | 3)&0x0000ffff;
                if(laneId==0) { tmpmem_s[64] = (tmpmem_s[64]&0xffff0000)|tmp; }
                return;
            }
            len2 = len1;
        }else if(len1>len2){
            uint tmpLen = len1-len2;
            genExtEmbCore_2V_1and2_noOverlap_withPhase<1,false>(neigData1_g+len2, neigData2_g, tmpLen, len2, tmpmem_s, equalVid, equalVNum,equalVNum, totWriteRowNum, basenewEmb_g, laneId, embLen,maxRowNum);
            if(tmpmem_s[66] == maxRowNum){
                tmp = tmpmem_s[64];
                tmp = ((tmp<<4) | 4)&0x0000ffff;
                if(laneId==0) { tmpmem_s[64] = (tmpmem_s[64]&0xffff0000)|tmp; }
                return;
            }
            len1 = len2;
        }
        genExtEmbCore_2V_1and2_noRestrict_fullOverlap_withPhase<false>(neigData1_g,len1,tmpmem_s,equalVid,0,equalVNum,totWriteRowNum,basenewEmb_g,laneId,embLen,maxRowNum);
        if(tmpmem_s[66] == maxRowNum){
            tmp = tmpmem_s[64];
            tmp = ((tmp<<4) | 5)&0x0000ffff;
            if(laneId==0) { tmpmem_s[64] = (tmpmem_s[64]&0xffff0000)|tmp; }
            return;
        }
    }else{
        uint stopindex = tmpmem_s[64]&0x0000000f;
        if(laneId==0) { tmpmem_s[64] = ((tmpmem_s[64]&0x0000ffff)>>4)| (tmpmem_s[64]&0xffff0000); }
        if(neigData1_g<neigData2_g) {
            uint tmpLen = neigData2_g-neigData1_g;
            if(stopindex==1){
                genExtEmbCore_2V_1and2_noOverlap_withPhase<1,true>(neigData1_g, neigData2_g, tmpLen, len2, tmpmem_s, equalVid, equalVNum,equalVNum, totWriteRowNum, basenewEmb_g, laneId, embLen,maxRowNum);
            }
            neigData1_g = neigData2_g;
            len1 = len1 - tmpLen;
        }else if(neigData2_g<neigData1_g) {
            uint tmpLen = neigData1_g-neigData2_g;
            if(stopindex<=2){
                genExtEmbCore_2V_1and2_noOverlap_withPhase<1,true>(neigData1_g, neigData2_g, len1, tmpLen, tmpmem_s, equalVid, equalVNum,equalVNum, totWriteRowNum, basenewEmb_g, laneId, embLen,maxRowNum);
            }
            neigData2_g = neigData1_g;
            len2 = len2 - tmpLen;
        }
        if(len1<len2) {
            uint tmpLen = len2-len1;
            if(stopindex<=3){
                genExtEmbCore_2V_1and2_noOverlap_withPhase<1,true>(neigData1_g, neigData2_g+len1, len1, tmpLen, tmpmem_s, equalVid, equalVNum, equalVNum, totWriteRowNum, basenewEmb_g, laneId, embLen,maxRowNum);
            }
            len2 = len1;
        }else if(len1>len2){
            uint tmpLen = len1-len2;
            if(stopindex<=4){
                genExtEmbCore_2V_1and2_noOverlap_withPhase<1,true>(neigData1_g+len2, neigData2_g, tmpLen, len2, tmpmem_s, equalVid, equalVNum,equalVNum, totWriteRowNum, basenewEmb_g, laneId, embLen,maxRowNum);
            }
            len1 = len2;
        }
        genExtEmbCore_2V_1and2_noRestrict_fullOverlap_withPhase<true>(neigData1_g,len1,tmpmem_s,equalVid,0,equalVNum,totWriteRowNum,basenewEmb_g,laneId,embLen,maxRowNum);
    }
}

//incoked by extEmb_2V_2src_1L
//needs 96 shared mem
//this is modified
template<bool isContinue>
__device__ void genExtEmb_2V_2src_1and2_noRestrict_sameLabel_evalEq_withPhase(uint *neigData1_g, uint *neigData2_g,
           uint len1, uint len2, uint *tmpmem_s, uint *equalVertices_s, uint *totWriteRowNum,uint *basenewEmb_g,
           uint laneId, uint embLen, uint maxRowNum) {

    uint i,predicate,totcount, *newEmb_g,writePos,tmp,stopindex,loopindex;
    uint equalVNum = equalVertices_s[1];
    uint equalVPos = equalVertices_s[2+laneId];
    uint equalVid = laneId<equalVNum?tmpmem_s[equalVPos]:0;
    if(isContinue) { tmp = tmpmem_s[32+laneId]; }
    for(i=0;i<equalVNum;++i){
        equalVPos = __shfl_sync(0xffffffff,equalVid,i);
        predicate = equalVPos>equalVid;
        predicate = __ballot_sync(0xffffffff, predicate);
        predicate = __popc(predicate)-(32-equalVNum);
        if(laneId==0){
            tmpmem_s[32+predicate] = equalVPos;
        }
    }
    equalVid = laneId<equalVNum?tmpmem_s[32+laneId]:0;
    if(isContinue) { tmpmem_s[32+laneId] = tmp; }
    if(isContinue){
        stopindex = tmpmem_s[64+5]&0x0000000f;
        loopindex = (tmpmem_s[64+5]&0xffff0000)>>16;
    }
    if(len1<len2) {
        i=laneId; totcount=0;
        uint outerVid, tot_32 = (len1+31)&0xffffffe0;
        while(totcount==0 && i<tot_32) {
            outerVid = i<len1?neigData1_g[i]:0;
            REARRANGE(outerVid,predicate,totcount,0,equalVNum,64)
            i=i+32;
        }
        if(totcount==0) { if(laneId==0) { tmpmem_s[66]=0;} return; }
        if(!isContinue){
            tmpmem_s[32+laneId] = 0;
            ALLOCSPACE(len2,embLen)
            if(writePos+len2>=maxRowNum){
                if(laneId==0){
                    tmpmem_s[64] = 1;
                    tmpmem_s[65] = len2;
                    tmpmem_s[66] = maxRowNum;
                }
                return;
            }
            genExtEmbCore_2V_1and2_sameLabel_detectPhase<false,true>(outerVid,neigData2_g,len2,tmpmem_s,equalVid,0,equalVNum,newEmb_g,laneId,embLen);
        }else{
            if(stopindex==1){
                tmpmem_s[32+laneId] = 0;
                ALLOCSPACE(len2,embLen)
                genExtEmbCore_2V_1and2_sameLabel_detectPhase<false,true>(outerVid,neigData2_g,len2,tmpmem_s,equalVid,0,equalVNum,newEmb_g,laneId,embLen);
            }
        }
        uint phaseNum = tmpmem_s[32];
        totcount = totcount - 1;
        if(totcount==0) { if(laneId==0) { tmpmem_s[66]=0;}return; }
        outerVid = __shfl_down_sync(0xffffffff,outerVid,1);
        if(!isContinue){
            ALLOCSPACE(((len2-phaseNum)*totcount),embLen)
            if(writePos+((len2-phaseNum)*totcount)>=maxRowNum){
                if(laneId==0){
                    tmpmem_s[64] = 2;
                    tmpmem_s[65] = ((len2-phaseNum)*totcount);
                    tmpmem_s[66] = maxRowNum;
                }
                return;
            }
        }else{
            if(stopindex<=2){
                ALLOCSPACE(((len2-phaseNum)*totcount),embLen)
            }
        }       
        if(tmpmem_s[32+1+phaseNum-1]!=len2-1) {
            if(laneId==0) { tmpmem_s[32+1+phaseNum] = len2; }
            phaseNum = phaseNum+1;
            //if(laneId==0) { tmpmem_s[32] = phaseNum; }
        }
        if(!isContinue){
            genExtEmbCore_2V_1and2_noRestrict_sameLabel_usePhase<false,true>(outerVid,neigData2_g,totcount,len2,tmpmem_s,phaseNum,newEmb_g,laneId,embLen);
        }else{
            if(stopindex<=2){
                genExtEmbCore_2V_1and2_noRestrict_sameLabel_usePhase<false,true>(outerVid,neigData2_g,totcount,len2,tmpmem_s,phaseNum,newEmb_g,laneId,embLen);
            }
        }
        if(len1>32) {
            if(isContinue){
                if(stopindex==3){
                    i=loopindex+laneId;
                }
            }
            for(;i<tot_32;i=i+32) {
                outerVid = i<len1?neigData1_g[i]:0;
                REARRANGE(outerVid,predicate,totcount,0,equalVNum,64)
                if(totcount==0) { continue; }
                ALLOCSPACE(((len2-phaseNum)*totcount),embLen)
                if(!isContinue){
                    if(writePos+((len2-phaseNum)*totcount)>=maxRowNum){
                        if(laneId==0){
                            tmpmem_s[64] = ((i>>5)<<16)|3;
                            tmpmem_s[65] = ((len2-phaseNum)*totcount);
                            tmpmem_s[66] = maxRowNum;
                        }
                        return;
                    }
                }
                genExtEmbCore_2V_1and2_noRestrict_sameLabel_usePhase<false,true>(outerVid,neigData2_g,totcount,len2,tmpmem_s,phaseNum,newEmb_g,laneId,embLen);
            }
        }
    }else {
        i=laneId; totcount=0;
        uint innerVid, tot_32 = (len2+31)&0xffffffe0;
        while(totcount==0 && i<tot_32) {
            innerVid = i<len2?neigData2_g[i]:0;
            REARRANGE(innerVid,predicate,totcount,0,equalVNum,64)
            i=i+32;
        }
        if(totcount==0) { if(laneId==0) {tmpmem_s[66] = 0;} return; }
        if(!isContinue){
            ALLOCSPACE(len1,embLen)
            tmpmem_s[32+laneId] = 0;
            if(writePos+len1>=maxRowNum){
                if(laneId==0){
                    tmpmem_s[64] = 1;
                    tmpmem_s[65] = len1;
                    tmpmem_s[66] = maxRowNum;
                }
                return;
            }
            genExtEmbCore_2V_1and2_sameLabel_detectPhase<true,true>(innerVid,neigData1_g,len1,tmpmem_s,equalVid,0,equalVNum,newEmb_g,laneId,embLen);
        }else{
            if(stopindex==1){
                ALLOCSPACE(len1,embLen)
                tmpmem_s[32+laneId] = 0;
                genExtEmbCore_2V_1and2_sameLabel_detectPhase<true,true>(innerVid,neigData1_g,len1,tmpmem_s,equalVid,0,equalVNum,newEmb_g,laneId,embLen);
            }
        }
        uint phaseNum = tmpmem_s[32];
        totcount = totcount - 1;
        if(totcount==0){ if(laneId==0) {tmpmem_s[66] = 0;} return; }
        innerVid = __shfl_down_sync(0xffffffff,innerVid,1);
        if(!isContinue){
            ALLOCSPACE(((len1-phaseNum)*totcount),embLen)
            if(writePos+((len1-phaseNum)*totcount)>=maxRowNum){
                if(laneId==0){
                    tmpmem_s[64] = 2;
                    tmpmem_s[65] = ((len1-phaseNum)*totcount);
                    tmpmem_s[66] = maxRowNum;
                }
                return;
            }
        }else{
            if(stopindex<=2){
                ALLOCSPACE(((len1-phaseNum)*totcount),embLen)
            }
        }
        if(tmpmem_s[32+1+phaseNum-1]!=len1-1) {
            if(laneId==0) { tmpmem_s[32+1+phaseNum] = len1; }
            phaseNum = phaseNum+1;
            //if(laneId==0) { tmpmem_s[32] = phaseNum; }
        }
        if(!isContinue){
            genExtEmbCore_2V_1and2_noRestrict_sameLabel_usePhase<true,true>(innerVid,neigData1_g,totcount,len1,tmpmem_s,phaseNum,newEmb_g,laneId,embLen);
        }else{
            if(stopindex<=2){
                genExtEmbCore_2V_1and2_noRestrict_sameLabel_usePhase<true,true>(innerVid,neigData1_g,totcount,len1,tmpmem_s,phaseNum,newEmb_g,laneId,embLen);
            }
        }
        if(len2>32) {
            if(isContinue){
                if(stopindex==3){
                    i=loopindex+laneId;
                }
            }
            for(;i<tot_32;i=i+32) {
                innerVid = i<len2?neigData2_g[i]:0;
                REARRANGE(innerVid,predicate,totcount,0,equalVNum,64)
                if(totcount==0) { continue; }
                ALLOCSPACE(((len1-phaseNum)*totcount),embLen)
                if(!isContinue){
                    if(writePos+((len1-phaseNum)*totcount)>=maxRowNum){
                        if(laneId==0){
                            tmpmem_s[64] = ((i>>5)<<16)|3;
                            tmpmem_s[65] = ((len1-phaseNum)*totcount);
                            tmpmem_s[66] = maxRowNum;
                        }
                        return;
                    }
                }
                genExtEmbCore_2V_1and2_noRestrict_sameLabel_usePhase<true,true>(innerVid,neigData1_g,totcount,len1,tmpmem_s,phaseNum,newEmb_g,laneId,embLen);
            }
        }
    }
    if(laneId==0){
        tmpmem_s[66] = 0;
    }
}

//invoked by extemb_2V_1src_2L, extEmb_2V_2src_2L
//equalVNum1 is the number for extvid1label, equalVNum2 is the number of both extvid1label and extvid2label
//tmpmem_s 96
//this is modified
template<bool isContinue>
__device__ void genExtEmb_2V_1and2_notSameLabel_withPhase(uint *neigData1_g, uint *neigData2_g, uint len1, uint len2, uint *tmpmem_s,
           uint *equalVertices_s, uint *totWriteRowNum, uint *basenewEmb_g,uint laneId, uint embLen,uint maxRowNum) {
    
    uint i, predicate,tmp;
    uint equalVNum1 = equalVertices_s[0];
    uint equalVNum2 = equalVertices_s[1];
    uint equalVPos = equalVertices_s[2+laneId];
    uint equalVid = laneId<equalVNum2?tmpmem_s[equalVPos]:0;
    if(isContinue){
        tmp = tmpmem_s[32+laneId];
    }
    for(i=0;i<equalVNum1;++i){
        equalVPos = __shfl_sync(0xffffffff,equalVid,i);
        predicate = laneId<equalVNum1?equalVPos>equalVid:0;
        predicate = __ballot_sync(0xffffffff, predicate);
        predicate = __popc(predicate);
        if(laneId==0){
            tmpmem_s[32+predicate] = equalVPos;
        }
    }
    for(i=equalVNum1;i<equalVNum2;++i){
        equalVPos = __shfl_sync(0xffffffff,equalVid,i);
        predicate = laneId>=equalVNum1 && laneId<equalVNum2?equalVPos>equalVid:0;
        predicate = __ballot_sync(0xffffffff, predicate);
        predicate = __popc(predicate);
        if(laneId==0){
            tmpmem_s[32+equalVNum1+predicate] = equalVPos;
        }
    }
    equalVid = laneId<equalVNum2?tmpmem_s[32+laneId]:0;
    if(isContinue){
        tmpmem_s[32+laneId] = tmp;
    }
    genExtEmbCore_2V_1and2_noOverlap_withPhase<2,isContinue>(neigData1_g,neigData2_g,len1,len2,tmpmem_s,equalVid,equalVNum1,equalVNum2,totWriteRowNum,basenewEmb_g,laneId,embLen,maxRowNum);
}




//return the number of vertices that lower than or equal to greatV
__device__ uint findLessVPosNoIndex(uint *neigData_g, uint neigLen, uint lessV, uint laneId) {
    uint tmp2, i, predicate = 0;
    tmp2 = ((neigLen + 31) & 0xffffffe0);
    for (i = laneId; i < tmp2; i = i + 32) {
        predicate = i < neigLen ? (neigData_g[i] < lessV) : 0;
        predicate = __ballot_sync(0xffffffff, predicate);
        if (predicate < 0xffffffff) { predicate = __popc(predicate); return (i & 0xffffffe0) + predicate; }
    }
    return neigLen;
}

__device__ uint findLessVPosWithIndex(uint *neigData_g, uint neigLen, uint lessV, uint laneId, uint indexNum){
    uint tmp2, i, predicate=0, upperLimit,lowerLimit;
    upperLimit = laneId<indexNum?neigData_g[laneId]:0;
    lowerLimit = __shfl_up_sync(0xffffffff,upperLimit,1);
    if(laneId==0){ lowerLimit = neigData_g[indexNum]-1; } //neigData_g[indexNum] is the first and smallest vertex of neighbors
    predicate = (lowerLimit<lessV && lessV<=upperLimit);
    predicate = __ballot_sync(0xffffffff,predicate);
    if(predicate==0){ return 0; }
    tmp2 = __ffs(predicate);
    uint blockSize = neigLen<=32*VBLOCKSIZE?VBLOCKSIZE:(neigLen>>5);
    lowerLimit = tmp2*blockSize;
    blockSize = tmp2==(indexNum-1)?neigLen-blockSize*tmp2:blockSize;
    upperLimit = (blockSize+31)&0xffffffe0;
    for(i=laneId;i<upperLimit;i=i+32){
        predicate = i<blockSize?(neigData_g[indexNum+lowerLimit+i]<lessV):0;
        predicate = __ballot_sync(0xffffffff,predicate);
        if(predicate<0xffffffff){ predicate = __popc(predicate); return lowerLimit+(i&0xffffffe0)+predicate; }
    }
    return lowerLimit+blockSize;
}

#define INIT_VARS_SHAREDMEM_NOAUX(SIZE)                                         \
    laneId = threadIdx.x & 31;                                                  \
    gridWarpNum = (gridDim.x * blockDim.x)>>5;                                       \
    warpIdInBlock = threadIdx.x >> 5;                                           \
    __shared__ uint neigV[WARPPERBLOCK][SIZE];                                  \
    __shared__ uint indexForIndex[256*2+1+1+64+1];                                    \
    for(i=threadIdx.x;i<intervalNum*2+1+1;i=i+blockDim.x){                        \
        indexForIndex[i] = edgeLabelPartition[i];                               \
    }                                                                           \
    i = blockIdx.x * blockDim.x + threadIdx.x;                                  \
    i >>= 5;                                                                    \
    edgeLabelPartition = edgeLabelPartition+intervalNum*2+1+1;\
    if(threadIdx.x==0){ indexForIndex[256*2+66] = 0; }

//before this macro, edgeLabelPartition points to array0, after this macro
//edgeLabelPartition points to array1
#define INIT_VARS_SHAREDMEM(SIZE)                                               \
    laneId = threadIdx.x & 31;                                                  \
    gridWarpNum = (gridDim.x * blockDim.x)>>5;                                  \
    warpIdInBlock = threadIdx.x >> 5;                                           \
    __shared__ uint neigV[WARPPERBLOCK][SIZE];                                  \
    __shared__ uint indexForIndex[256*2+1+1+64+1];                                \
    for(i=threadIdx.x;i<intervalNum*2+1+1;i=i+blockDim.x){                      \
        indexForIndex[i] = edgeLabelPartition[i];                               \
    }                                                                           \
    for(i=threadIdx.x;i<64;i=i+blockDim.x){                                     \
        indexForIndex[256*2+1+1+i] = auxArray[i];                               \
    }                                                                           \
    i = blockIdx.x * blockDim.x + threadIdx.x;                                  \
    i >>= 5;                                                                    \
    edgeLabelPartition = edgeLabelPartition+intervalNum*2+1+1;\
    if(threadIdx.x==0){ indexForIndex[256*2+66] = 0; }

//init: generate embeddings of the first partial pattern
//ext: extension
//genRec: generate positions of edge label--end vertex label combination
//edgeLabelPart: the position of array2
//vLabel: label id of each vertex
//totWriteRowNum is a global number, all lane 0 accumulate their number of generated embeddings onto this variable
//recordPos: store adresses of neighbor into recordPos. recordPos[i*2]=0,recordPos[i*2+1]=0 means the label doesn't fit, we do not
//have the change to find its neighbor address for this label. recordPos[i*2]=1,recordPos[i*2+1]=0 means this svid does not have
//neighbors of this label
//restrictFlag:0-7 bits, extvid1<v[pos]; 8-15 bits, extvid2<v[pos]; 16-23 bits, is extvid1<extvid2, 0, 1; 0xff means no restrict
//indexForIndex[0:2*256] are interval indexs, indexForIndex[2*256] and indexForIndex[2*256+1] are labels for extvid1 and extvid2
//indexForIndex[2*256+2] and indexForIndex[2*256+3] are positions of lessvid for extvid1 and extvid2
//the format for indexForIndex[0:2*256] is vs1,vs2,vs3,len1,len2,len3. len1 is the numberr of vertices in interval 1, len2
//is the number of vertices in all previsous intervals (including interval 2).
//neighborsData is the address of array4
//this is modified
template<bool isExt1Restrict, bool isExt2Restrict, bool isExt1and2Restrict, bool isExt1and2SameLabel, bool isRecord1, bool isRecord2>
__global__ void initEmb_3V_1and2_NoHash_kernel(uint *vLabel, uint *totWriteRowNum, uint *edgeLabelPartition, 
    uint *neighborsData, uint *recordPos1, uint *recordPos2, uint *newEmb, uint intervalNum, uint partialRowNum, 
    uint svidlabel,uint evid1label,uint evid2label, bool isContinue, uint maxRowNum){

    uint laneId, gridWarpNum, warpIdInBlock, i, j;
    INIT_VARS_SHAREDMEM_NOAUX(32)
    __syncthreads();
    uint intervalStart = laneId,predicate;
    if(isContinue){
        uint writeEle = laneId<7?newEmb[(maxRowNum+i)*7+laneId]:0;
        neigV[warpIdInBlock][laneId] = writeEle; 
        uint *neigAddr1 = neighborsData+neigV[warpIdInBlock][0];
        uint neigLen1 = neigV[warpIdInBlock][1];
        uint *neigAddr2 = neighborsData+neigV[warpIdInBlock][2];
        uint neigLen2 = neigV[warpIdInBlock][3];
        i = neigV[warpIdInBlock][4];
        intervalStart = neigV[warpIdInBlock][5]+laneId;
        uint svid = neigV[warpIdInBlock][6];
        if(laneId==0) { neigV[warpIdInBlock][0] = svid;}
        if(isExt1and2SameLabel){
            if (isExt1and2Restrict) {
                genInitEmb_3V_1and2_restric(neigAddr1,neigAddr2,neigLen1,neigLen2,neigV[warpIdInBlock],newEmb,totWriteRowNum,laneId,svid,maxRowNum);
            } else {
                genInitEmb_3V_1and2_noRestrict_sameLabel(neigAddr1,neigAddr2,neigLen1,neigLen2,neigV[warpIdInBlock],newEmb,totWriteRowNum,laneId,maxRowNum);
            }
        }else{
            genInitEmb_3V_1and2_notsamelabel(neigAddr1,neigAddr2,neigLen1,neigLen2,neigV[warpIdInBlock],newEmb,totWriteRowNum,laneId,maxRowNum);
        }
        i=i+gridWarpNum;
    }
    while (i < partialRowNum){
        if(isRecord1){
            if(laneId<2) { recordPos1[i*2+laneId] = 0; }
        }
        if(isRecord2){
            if(laneId<2) { recordPos2[i*2+laneId] = 0; }
        }
        uint lowerLimit, upperLimit, svid;
        uint len_32 = (intervalNum+31)&0xffffffe0;
        for(j=intervalStart;j<len_32;j=j+32){
            upperLimit = indexForIndex[j+intervalNum+1];
            lowerLimit = indexForIndex[j+intervalNum];
            //i+1 because vid starts from 1
            predicate = j<intervalNum?(i+1>lowerLimit && i+1<=upperLimit):0;
            predicate = __ballot_sync(0xffffffff,predicate);
            if(predicate>0) {
                upperLimit=i-lowerLimit;
                lowerLimit=indexForIndex[j];
                svid=lowerLimit+upperLimit;
                intervalStart = j;
                break;
            }
        }
        if(predicate>0){
            uint tmpIndex = __ffs(predicate);
            svid = __shfl_sync(0xffffffff,svid,tmpIndex-1);
        }
        else{ /*std::cout<<"wrong"<<std::endl;*/ return; }

        if(vLabel[svid]!=svidlabel){ 
            i=i+gridWarpNum; 
            continue; 
        }

        if(laneId==0) { neigV[warpIdInBlock][0] = svid;}
        uint *neigAddr = neighborsData + edgeLabelPartition[i];
        uint neigLen = edgeLabelPartition[i+1]-edgeLabelPartition[i];
        lowerLimit = 0; 
        upperLimit = neigLen;
        uint lessThan, greatThan;
        if(isExt1and2SameLabel) {
            uint found;
            FIND_LABEL_LIMIT(neigAddr,lowerLimit,upperLimit,evid1label)
            neigLen = upperLimit - lowerLimit;
            if (neigLen == 0) { 
                if(isRecord1){ 
                    if (laneId==0) {
                        recordPos1[i*2]=1;recordPos1[i*2+1]=0;
                    } 
                } 
                i=i+gridWarpNum; 
                continue; 
            }
            uint indexNum = 0,neigNum;
            if(neigLen>VBLOCKSIZE+2 && neigLen<=32*VBLOCKSIZE+32) {indexNum=(neigLen+VBLOCKSIZE)/(VBLOCKSIZE+1);}
            else if (neigLen > 32 * VBLOCKSIZE + 32) { indexNum = 32; }
            neigAddr = neigAddr + lowerLimit;//address starts from index terms
            neigLen = neigLen - indexNum;
            //recorded positon starts from index terms, recorded length is the length of neghbors (not inlcude index terms)
            if (isRecord1) { if(laneId==0) {recordPos1[i*2]=lowerLimit;recordPos1[i*2+1]=neigLen;} }
            if (isExt1Restrict && isExt2Restrict) {
                neigNum = neigLen;
                if (indexNum == 0) { neigLen = findLessVPosNoIndex(neigAddr, neigNum, svid, laneId); }
                else { neigLen = findLessVPosWithIndex(neigAddr, neigNum, svid, laneId, indexNum); }
                if (neigLen == 0) { i=i+gridWarpNum;continue; }
                neigNum = neigLen;
                if (isExt1and2Restrict) {
                    genInitEmb_3V_1and2_restric(neigAddr+indexNum,neigAddr+indexNum,neigLen,neigNum,neigV[warpIdInBlock],newEmb,totWriteRowNum,laneId,svid,maxRowNum);
                } else {
                    genInitEmb_3V_1and2_noRestrict_sameLabel(neigAddr+indexNum,neigAddr+indexNum,neigLen,neigNum,neigV[warpIdInBlock],newEmb,totWriteRowNum,laneId,maxRowNum);
                }
            } else if (isExt1Restrict) {
                neigNum = neigLen;
                if (indexNum == 0) { neigLen = findLessVPosNoIndex(neigAddr, neigNum, svid, laneId); }
                else { neigLen = findLessVPosWithIndex(neigAddr, neigNum, svid, laneId, indexNum); }
                if (neigLen == 0) { i=i+gridWarpNum;continue; }
                if (isExt1and2Restrict) {
                    genInitEmb_3V_1and2_restric(neigAddr+indexNum,neigAddr+indexNum,neigLen,neigNum,neigV[warpIdInBlock],newEmb,totWriteRowNum,laneId,svid,maxRowNum);
                } else {
                    genInitEmb_3V_1and2_noRestrict_sameLabel(neigAddr+indexNum,neigAddr+indexNum,neigLen,neigNum,neigV[warpIdInBlock],newEmb,totWriteRowNum,laneId,maxRowNum);
                }
            } else if (isExt2Restrict) {
                if (indexNum == 0) { neigNum = findLessVPosNoIndex(neigAddr, neigLen, svid, laneId); }
                else { neigNum = findLessVPosWithIndex(neigAddr, neigLen, svid, laneId, indexNum); }
                if (neigNum == 0) { i=i+gridWarpNum;continue; }
                if (isExt1and2Restrict) {
                    genInitEmb_3V_1and2_restric(neigAddr+indexNum,neigAddr+indexNum,neigLen,neigNum,neigV[warpIdInBlock],newEmb,totWriteRowNum,laneId,svid,maxRowNum);
                } else {
                    genInitEmb_3V_1and2_noRestrict_sameLabel(neigAddr+indexNum,neigAddr+indexNum,neigLen,neigNum,neigV[warpIdInBlock],newEmb,totWriteRowNum,laneId,maxRowNum);
                }
            } else {
                neigNum = neigLen;
                if (isExt1and2Restrict) {
                    genInitEmb_3V_1and2_restric(neigAddr+indexNum,neigAddr+indexNum,neigLen,neigNum,neigV[warpIdInBlock],newEmb,totWriteRowNum,laneId,svid,maxRowNum);
                } else {
                    genInitEmb_3V_1and2_noRestrict_sameLabel(neigAddr+indexNum,neigAddr+indexNum,neigLen,neigNum,neigV[warpIdInBlock],newEmb,totWriteRowNum,laneId,maxRowNum);
                }
            }
            if(neigV[warpIdInBlock][0]==maxRowNum){
                uint index = blockIdx.x*blockDim.x+threadIdx.x;
                uint writeEle = neigAddr+indexNum-neighborsData;
                if(laneId==1){
                    writeEle = neigLen;
                }else if(laneId==2){
                    writeEle = neigAddr+indexNum-neighborsData;
                }else if(laneId==3){
                    writeEle = neigNum;
                }else if(laneId==4){
                    writeEle = i;
                }else if(laneId==5){
                    writeEle = intervalStart>>5;
                }else if(laneId==6){
                    writeEle = svid;
                }
                if(laneId<7){
                    newEmb[(maxRowNum+index)*7+laneId] = writeEle;
                }
                writeEle = neigV[warpIdInBlock][1];
                if(laneId==0){
                    atomicAdd(indexForIndex+256*2+66,writeEle);
                }
                __syncthreads();
                if(threadIdx.x==0){
                    atomicAdd(totWriteRowNum+1,indexForIndex[256*2+66]);
                }
                break;
            }
        }else {
            uint upperLimit2 = upperLimit, lowerLimit2 = lowerLimit, found;
            FIND_LABEL_LIMIT(neigAddr,lowerLimit,upperLimit,evid1label)
            neigLen = upperLimit - lowerLimit;
            if (neigLen == 0) { if(isRecord1){ if (laneId==0) {recordPos1[i*2]=1;recordPos1[i*2+1]=0;} } i=i+gridWarpNum; continue; }
            FIND_LABEL_LIMIT(neigAddr, lowerLimit2, upperLimit2, evid2label)
            uint neigLen2 = upperLimit2 - lowerLimit2;
            if (neigLen2 == 0) { if(isRecord2){ if (laneId==0) {recordPos2[i*2]=1;recordPos2[i*2+1]=0;} } i=i+gridWarpNum; continue; }
            uint indexNum = 0,neigNum, indexNum2=0,neigNum2;
            if(neigLen>VBLOCKSIZE+2 && neigLen<=32*VBLOCKSIZE+32) {indexNum=(neigLen+VBLOCKSIZE)/(VBLOCKSIZE+1);}
            else if (neigLen > 32 * VBLOCKSIZE + 32) { indexNum = 32; }
            if(neigLen2>VBLOCKSIZE+2 && neigLen2<=32*VBLOCKSIZE+32) {indexNum2=(neigLen2+VBLOCKSIZE)/(VBLOCKSIZE+1);}
            else if (neigLen2 > 32 * VBLOCKSIZE + 32) { indexNum2 = 32; }

            uint *neigAddr2 = neigAddr + lowerLimit2;
            neigLen2 = neigLen2 - indexNum2;
            neigAddr = neigAddr + lowerLimit;
            neigLen = neigLen - indexNum;
            if (isRecord1) { if(laneId==0) {recordPos1[i*2]=lowerLimit;recordPos1[i*2+1]=neigLen;} }
            if (isRecord2) { if(laneId==0) {recordPos2[i*2]=lowerLimit2;recordPos2[i*2+1]=neigLen2;} }
            if (isExt1Restrict) {
                if (indexNum == 0) { neigLen = findLessVPosNoIndex(neigAddr, neigLen, svid, laneId); }
                else { neigLen = findLessVPosWithIndex(neigAddr, neigLen, svid, laneId, indexNum); }
                if (neigLen == 0) { i=i+gridWarpNum;continue; }
                genInitEmb_3V_1and2_notsamelabel(neigAddr+indexNum,neigAddr2+indexNum2,neigLen,neigLen2,neigV[warpIdInBlock],newEmb,totWriteRowNum,laneId,maxRowNum);
            } else if (isExt2Restrict) {
                if (indexNum2 == 0) { neigLen2 = findLessVPosNoIndex(neigAddr2, neigLen2, svid, laneId); }
                else { neigLen2 = findLessVPosWithIndex(neigAddr2, neigLen2, svid, laneId, indexNum2); }
                if (neigLen2 == 0) { i=i+gridWarpNum;continue; }
                genInitEmb_3V_1and2_notsamelabel(neigAddr+indexNum,neigAddr2+indexNum2,neigLen,neigLen2,neigV[warpIdInBlock],newEmb,totWriteRowNum,laneId,maxRowNum);
            } else {
                genInitEmb_3V_1and2_notsamelabel(neigAddr+indexNum,neigAddr2+indexNum2,neigLen,neigLen2,neigV[warpIdInBlock],newEmb,totWriteRowNum,laneId,maxRowNum);
            }
            if(neigV[warpIdInBlock][0]==maxRowNum){
                uint index = blockIdx.x*blockDim.x+threadIdx.x;
                uint writeEle = neigAddr+indexNum-neighborsData;
                if(laneId==1){
                    writeEle = neigLen;
                }else if(laneId==2){
                    writeEle = neigAddr2+indexNum2-neighborsData;
                }else if(laneId==3){
                    writeEle = neigLen2;
                }else if(laneId==4){
                    writeEle = i;
                }else if(laneId==5){
                    writeEle = intervalStart>>5;
                }else if(laneId==6){
                    writeEle = svid;
                }
                if(laneId<7){
                    newEmb[(maxRowNum+index)*7+laneId] = writeEle;
                }
                if(laneId==0){
                    atomicAdd(indexForIndex+256*2+66,neigLen*neigLen2);
                }
                __syncthreads();
                if(threadIdx.x==0){
                    atomicAdd(totWriteRowNum+1,indexForIndex[256*2+66]);
                }
                break;
            }
        }
        i = i+gridWarpNum;
    }
}

//this is modified
template<bool isExtRestrict, bool isRecord>
__global__ void initEmb_2V_NoHash_kernel(uint *vLabel, uint *totWriteRowNum, uint *edgeLabelPartition, uint *neighborsData,uint *recordPos,
    uint *newEmb, uint intervalNum, uint partialRowNum, uint svidlabel,uint evidlabel, bool isContinue,uint maxRowNum){

    uint laneId, gridWarpNum, warpIdInBlock, i, j;
    INIT_VARS_SHAREDMEM_NOAUX(32)
    __syncthreads();
    uint intervalStart = laneId;
    if(isContinue){
        uint writeEle = laneId<5?newEmb[(maxRowNum+i)*5+laneId]:0;
        neigV[warpIdInBlock][laneId] = writeEle; 
        i = neigV[warpIdInBlock][2];
        uint svid = neigV[warpIdInBlock][3];
        intervalStart = neigV[warpIdInBlock][4]+laneId;
        uint *neigAddr = neighborsData+neigV[warpIdInBlock][0];
        uint neigLen = neigV[warpIdInBlock][1];
        genInitEmb_2V(neigAddr,neigLen,newEmb,totWriteRowNum,laneId,svid,neigV[warpIdInBlock],maxRowNum);
        i=i+gridWarpNum;
    }
    while (i < partialRowNum){
        uint lowerLimit, upperLimit, svid,predicate;
        uint len_32 = (intervalNum+31)&0xffffffe0;
        for(j=intervalStart;j<len_32;j=j+32){
            upperLimit = indexForIndex[j+intervalNum+1];
            lowerLimit = indexForIndex[j+intervalNum];
            //i+1 because vid starts from 1
            predicate = j<intervalNum?(i+1>lowerLimit && i+1<=upperLimit):0;
            predicate = __ballot_sync(0xffffffff,predicate);
            if(predicate>0) {
                upperLimit=i-lowerLimit;
                lowerLimit=indexForIndex[j];
                svid=lowerLimit+upperLimit;
                intervalStart = j;
                break;
            }
        }
        if(predicate>0){
            uint tmpIndex = __ffs(predicate);
            svid = __shfl_sync(0xffffffff,svid,tmpIndex-1);
        }else{ return; }
        if(vLabel[svid]!=svidlabel){ i=i+gridWarpNum; continue; }
        if(laneId==0) { neigV[warpIdInBlock][0] = svid;}
        uint *neigAddr = neighborsData + edgeLabelPartition[i];
        uint neigLen = edgeLabelPartition[i+1]-edgeLabelPartition[i],found;

        lowerLimit = 0;
        upperLimit = neigLen;
        uint lessThan,greatThan;
        FIND_LABEL_LIMIT(neigAddr,lowerLimit, upperLimit, evidlabel)
        neigLen = upperLimit - lowerLimit;
        if (neigLen == 0) { if(isRecord){ if (laneId==0) {recordPos[i*2]=1;recordPos[i*2+1]=0;} } i=i+gridWarpNum; continue; }
        uint indexNum = 0,neigNum;
        if(neigLen>VBLOCKSIZE+2 && neigLen<=32*VBLOCKSIZE+32) {indexNum=(neigLen+VBLOCKSIZE)/(VBLOCKSIZE+1);}
        else if (neigLen > 32 * VBLOCKSIZE + 32) { indexNum = 32; }
        neigAddr = neigAddr + lowerLimit;
        neigLen = neigLen - indexNum;
        if (isRecord) { if(laneId==0) {recordPos[i*2]=lowerLimit;recordPos[i*2+1]=neigLen;} }
        if(isExtRestrict) {
            if (indexNum == 0) { neigNum = findLessVPosNoIndex(neigAddr, neigLen, svid, laneId); }
            else { neigNum = findLessVPosWithIndex(neigAddr, neigLen, svid, laneId, indexNum); }
            if (neigNum == 0) { i=i+gridWarpNum; continue; }
            genInitEmb_2V(neigAddr+indexNum,neigNum,newEmb,totWriteRowNum,laneId,svid,neigV[warpIdInBlock],maxRowNum);
            if(neigV[warpIdInBlock][0]==maxRowNum){
                uint index = (blockIdx.x*blockDim.x+threadIdx.x)>>5;
                uint writeEle = neigAddr+indexNum-neighborsData;
                if(laneId==1){
                    writeEle = neigNum;
                }else if(laneId==2){
                    writeEle = i;
                }else if(laneId==3){
                    writeEle = svid;
                }else if(laneId==4){
                    writeEle = intervalStart>>5;
                }
                if(laneId<5){
                    newEmb[(maxRowNum+index)*5+laneId]=writeEle;
                }
                if(laneId==0){
                    atomicAdd(indexForIndex+256*2+66,neigNum);
                }
                __syncthreads();
                if(threadIdx.x==0){
                    atomicAdd(totWriteRowNum+1,indexForIndex[256*2+66]);
                }
                break;
            }
        }else{
            genInitEmb_2V(neigAddr+indexNum,neigLen,newEmb,totWriteRowNum,laneId,svid,neigV[warpIdInBlock],maxRowNum);
            if(neigV[warpIdInBlock][0]==maxRowNum){
                uint index = (blockIdx.x*blockDim.x+threadIdx.x)>>5;
                uint writeEle = neigAddr+indexNum-neighborsData;
                if(laneId==1){
                    writeEle = neigLen;
                }else if(laneId==2){
                    writeEle = i;
                }else if(laneId==3){
                    writeEle = svid;
                }else if(laneId==4){
                    writeEle = intervalStart>>5;
                }
                if(laneId<5){
                    newEmb[(maxRowNum+index)*5+laneId]= writeEle;
                }
                if(laneId==0){
                    atomicAdd(indexForIndex+256*2+66,neigLen);
                }
                __syncthreads();
                if(threadIdx.x==0){
                    atomicAdd(totWriteRowNum+1,indexForIndex[256*2+66]);
                }
                break;
            }
        }
        i = i+gridWarpNum;
    }
}


//this is modified
template<bool isExt1Restrict, bool isExt2Restrict, bool isExt1and2Restrict,bool isRecord, bool useRecord>
__global__ void extEmb_2V_2src_1and2_sameLabel_NoHash_kernel(uint *vLabel, uint *totWriteRowNum, uint *edgeLabelPartition, 
    uint *neighborsData, uint *recordPos, uint *auxArray, uint *newEmb, uint *partialEmb, uint intervalNum, 
    uint partialRowNum, uint embLen, bool isContinue, uint maxRowNum) {

    uint laneId, gridWarpNum, warpIdInBlock, i, j;
    INIT_VARS_SHAREDMEM(96)
    __syncthreads();
    if(isContinue){
        if(isExt1and2Restrict){
            uint writeEle = laneId<9?newEmb[(maxRowNum+i)*9+laneId]:0;
            neigV[warpIdInBlock][32+laneId] = writeEle; 
            i = neigV[warpIdInBlock][32+4];
            j = laneId<(embLen-2)?partialEmb[i*(embLen-2)+laneId]:1;
            neigV[warpIdInBlock][laneId] = j;
            uint *neigAddr1 = neighborsData+neigV[warpIdInBlock][32+0];
            uint neigLen1 = neighborsData[32+1];
            uint *neigAddr2 = neighborsData+neigV[warpIdInBlock][32+2];
            uint neigLen2 = neighborsData[32+3];
            if(neigLen1<neigLen2){
                genExtEmb_2V_2src_1and2_restrict<true,true>(neigAddr1,neigAddr2,neigLen1,neigLen2,neigV[warpIdInBlock],indexForIndex+2*256+1+1+6,totWriteRowNum,newEmb,newEmb,laneId,embLen,maxRowNum);
            }else{
                genExtEmb_2V_2src_1and2_restrict<false,true>(neigAddr1,neigAddr2,neigLen1,neigLen2,neigV[warpIdInBlock],indexForIndex+2*256+1+1+6,totWriteRowNum,newEmb,newEmb,laneId,embLen,maxRowNum);
            }
        }else{
            uint writeEle = newEmb[(maxRowNum+i)*38+laneId];
            neigV[warpIdInBlock][32+laneId] = writeEle;
            writeEle = laneId<6?newEmb[(maxRowNum+i)*38+32+laneId]:0;
            neigV[warpIdInBlock][64+laneId] = writeEle;
            i = neigV[warpIdInBlock][64+4];
            j = laneId<(embLen-2)?partialEmb[i*(embLen-2)+laneId]:1;
            neigV[warpIdInBlock][laneId] = j;
            uint *neigAddr1 = neighborsData+neigV[warpIdInBlock][64+0];
            uint neigLen1 = neighborsData[64+1];
            uint *neigAddr2 = neighborsData+neigV[warpIdInBlock][64+2];
            uint neigLen2 = neighborsData[64+3];
            genExtEmb_2V_2src_1and2_noRestrict_sameLabel_evalEq_withPhase<true>(neigAddr1,neigAddr2,neigLen1,neigLen2,neigV[warpIdInBlock],indexForIndex+2*256+1+1+6,totWriteRowNum,newEmb,laneId,embLen,maxRowNum);
        }
        i = i + gridWarpNum;
    }
    while (i < partialRowNum) {
        j = laneId<(embLen-2)?partialEmb[i*(embLen-2)+laneId]:1;
        neigV[warpIdInBlock][laneId] = j;
        j = __ballot_sync(0xffffffff,j);
        if(j<0xffffffff) { i=i+gridWarpNum; continue; }

        uint svid1 = indexForIndex[256*2+1+1+4];
        uint svid2 = indexForIndex[256*2+1+1+5];
        svid1 = neigV[warpIdInBlock][svid1];
        svid2 = neigV[warpIdInBlock][svid2];
        uint lowerLimit = 0, upperLimit, predicate = 0, index;
        uint len = (intervalNum+31)&0xffffffe0;
        for(j=laneId;j<len;j=j+32){
            lowerLimit = indexForIndex[j];
            upperLimit = indexForIndex[j+intervalNum+1];
            upperLimit = lowerLimit+upperLimit-indexForIndex[j+intervalNum];
            predicate = j<intervalNum?(svid1>=lowerLimit && svid1<upperLimit):0;
            predicate = __ballot_sync(0xffffffff,predicate);
            if(predicate>0){
                index = svid1-lowerLimit+indexForIndex[j+intervalNum];
                uint tmpIndex = __ffs(predicate)-1;
                index = __shfl_sync(0xffffffff,index,tmpIndex);
                break;
            }
        }
        if(predicate==0){ 
            i=i+gridWarpNum;
            continue; 
        }
        uint baseNum = edgeLabelPartition[index];
        uint *neigAddr1 = neighborsData+baseNum;
        uint neigLen1 = edgeLabelPartition[index+1]-baseNum;
        uint indexNum1 = 0;
        if(useRecord){
            uint distance = recordPos[index*2];
            uint tmpNeigLen1 = recordPos[index*2+1];
            if(distance==0 && tmpNeigLen1==0){
                uint evidlabel = indexForIndex[2*256+1+1+0],found,lessThan,greatThan;
                lowerLimit = 0; upperLimit = neigLen1;
                FIND_LABEL_LIMIT(neigAddr1, lowerLimit, upperLimit, evidlabel)
                neigLen1 = upperLimit - lowerLimit;
                if (neigLen1 == 0) {
                    if (isRecord) { if (laneId == 0) { recordPos[index * 2] = 1; recordPos[index * 2 + 1] = 0; } }
                    i = i + gridWarpNum;
                    continue;
                }
                if (neigLen1 > VBLOCKSIZE+2 && neigLen1 <= 32 * VBLOCKSIZE + 32) { indexNum1 = (neigLen1+VBLOCKSIZE) / (VBLOCKSIZE + 1); }
                else if (neigLen1 > 32 * VBLOCKSIZE + 32) { indexNum1 = 32; }
                neigAddr1 = neigAddr1 + lowerLimit;
                neigLen1 = neigLen1 - indexNum1;
                if (isRecord) { if (laneId == 0) { recordPos[index * 2] = lowerLimit; recordPos[index * 2 + 1] = neigLen1; } }
            }else if(distance==1 && tmpNeigLen1==0){
                i = i+gridWarpNum;
                continue;
            }else{
                neigAddr1 = neigAddr1+distance;
                neigLen1 = tmpNeigLen1;
            }
        }else {
            uint evidlabel = indexForIndex[2*256+1+1+0],found,lessThan,greatThan;
            lowerLimit = 0; upperLimit = neigLen1;
            FIND_LABEL_LIMIT(neigAddr1, lowerLimit, upperLimit, evidlabel)
            neigLen1 = upperLimit - lowerLimit;
            if (neigLen1 == 0) {
                if (isRecord) { if (laneId == 0) { recordPos[index * 2] = 1; recordPos[index * 2 + 1] = 0; } }
                i = i + gridWarpNum;
                continue;
            }
            if (neigLen1 > VBLOCKSIZE+2 && neigLen1 <= 32 * VBLOCKSIZE + 32) { indexNum1 = (neigLen1+VBLOCKSIZE) / (VBLOCKSIZE + 1); }
            else if (neigLen1 > 32 * VBLOCKSIZE + 32) { indexNum1 = 32; }
            neigAddr1 = neigAddr1 + lowerLimit;
            neigLen1 = neigLen1 - indexNum1;
            if (isRecord) { if (laneId == 0) { recordPos[index * 2] = lowerLimit; recordPos[index * 2 + 1] = neigLen1; } }
        }
        if(isExt1Restrict){
            uint lessvid = indexForIndex[2*256+1+1+2];
            lessvid = neigV[warpIdInBlock][lessvid];
            if (indexNum1 == 0) { neigLen1 = findLessVPosNoIndex(neigAddr1, neigLen1, lessvid, laneId); }
            else { neigLen1 = findLessVPosWithIndex(neigAddr1, neigLen1, lessvid, laneId, indexNum1); }
            if (neigLen1 == 0) { i=i+gridWarpNum; continue; }
        }



        for(j=laneId;j<len;j=j+32){
            lowerLimit = indexForIndex[j];
            upperLimit = indexForIndex[j+intervalNum+1];
            upperLimit = lowerLimit+upperLimit-indexForIndex[j+intervalNum];
            predicate = j<intervalNum?(svid2>=lowerLimit && svid2<upperLimit):0;
            predicate = __ballot_sync(0xffffffff,predicate);
            if(predicate>0){
                index = svid2-lowerLimit+indexForIndex[j+intervalNum];
                uint tmpIndex = __ffs(predicate)-1;
                index = __shfl_sync(0xffffffff,index,tmpIndex);
                break;
            }
        }
        if(predicate==0){ 
            i=i+gridWarpNum;
            continue; 
        }
        baseNum = edgeLabelPartition[index];
        uint *neigAddr2 = neighborsData+baseNum;
        uint neigLen2 = edgeLabelPartition[index+1]-baseNum;
        uint indexNum2 = 0;
        if(useRecord){
            uint distance = recordPos[index*2];
            uint tmpNeigLen2 = recordPos[index*2+1];
            if(distance==0 && tmpNeigLen2==0){
                uint evidlabel = indexForIndex[2*256+1+1+1],found,lessThan,greatThan;
                lowerLimit = 0; upperLimit = neigLen2;
                FIND_LABEL_LIMIT(neigAddr2, lowerLimit, upperLimit, evidlabel)
                neigLen2 = upperLimit - lowerLimit;
                if (neigLen2 == 0) {
                    if (isRecord) { if (laneId == 0) { recordPos[index * 2] = 1; recordPos[index * 2 + 1] = 0; } }
                    i = i + gridWarpNum;
                    continue;
                }
                if (neigLen2 > VBLOCKSIZE+2 && neigLen2 <= 32 * VBLOCKSIZE + 32) { indexNum2 = (neigLen2+VBLOCKSIZE) / (VBLOCKSIZE + 1); }
                else if (neigLen2 > 32 * VBLOCKSIZE + 32) { indexNum2 = 32; }
                neigAddr2 = neigAddr2 + lowerLimit;
                neigLen2 = neigLen2 - indexNum2;
                if (isRecord) { if (laneId == 0) { recordPos[index * 2] = lowerLimit; recordPos[index * 2 + 1] = neigLen2; } }
            }else if(distance==1 && tmpNeigLen2==0){
                i = i+gridWarpNum;
                continue;
            }else{
                neigAddr2 = neigAddr2+distance;
                neigLen2 = tmpNeigLen2;
            }
        }else {
            uint evidlabel = indexForIndex[2*256+1+1+1],found,lessThan,greatThan;
            lowerLimit = 0; upperLimit = neigLen2;
            FIND_LABEL_LIMIT(neigAddr2, lowerLimit, upperLimit, evidlabel)
            neigLen2 = upperLimit - lowerLimit;
            if (neigLen2 == 0) {
                if (isRecord) { if (laneId == 0) { recordPos[index * 2] = 1; recordPos[index * 2 + 1] = 0; } }
                i = i + gridWarpNum;
                continue;
            }
            if (neigLen2 > VBLOCKSIZE+2 && neigLen2 <= 32 * VBLOCKSIZE + 32) { indexNum2 = (neigLen2+VBLOCKSIZE) / (VBLOCKSIZE + 1); }
            else if (neigLen2 > 32 * VBLOCKSIZE + 32) { indexNum2 = 32; }
            neigAddr2 = neigAddr2 + lowerLimit;
            neigLen2 = neigLen2 - indexNum2;
            if (isRecord) { if (laneId == 0) { recordPos[index * 2] = lowerLimit; recordPos[index * 2 + 1] = neigLen2; } }
        }
        if(isExt2Restrict){
            uint lessvid = indexForIndex[2*256+1+1+3];
            lessvid = neigV[warpIdInBlock][lessvid];
            if (indexNum2 == 0) { neigLen2 = findLessVPosNoIndex(neigAddr2, neigLen2, lessvid, laneId); }
            else { neigLen2 = findLessVPosWithIndex(neigAddr2, neigLen2, lessvid, laneId, indexNum2); }
            if (neigLen2 == 0) { i=i+gridWarpNum; continue; }
        }
        if(isExt1and2Restrict){
            if(neigLen1<neigLen2){
                genExtEmb_2V_2src_1and2_restrict<true,false>(neigAddr1+indexNum1,neigAddr2+indexNum2,neigLen1,neigLen2,neigV[warpIdInBlock],indexForIndex+2*256+1+1+6,totWriteRowNum,newEmb,newEmb,laneId,embLen,maxRowNum);
            }else{
                genExtEmb_2V_2src_1and2_restrict<false,false>(neigAddr1+indexNum1,neigAddr2+indexNum2,neigLen1,neigLen2,neigV[warpIdInBlock],indexForIndex+2*256+1+1+6,totWriteRowNum,newEmb,newEmb,laneId,embLen,maxRowNum);
            }
            if(neigV[warpIdInBlock][66]==maxRowNum){
                index = (blockIdx.x*blockDim.x+threadIdx.x)>>5;
                uint writeEle = neigAddr1+indexNum1-neighborsData;
                if(laneId==1){ writeEle = neigLen1; }
                else if(laneId==2){ writeEle = neigAddr2+indexNum2-neighborsData; }
                else if(laneId==3){ writeEle = neigLen2; }
                else if(laneId==4){ writeEle = i; }
                else if(laneId==5){ writeEle = neigV[warpIdInBlock][32]; }
                else if(laneId==6){ writeEle = neigV[warpIdInBlock][33]; }
                else if(laneId==7){ writeEle = neigV[warpIdInBlock][35]; }
                else if(laneId==8){ writeEle = neigV[warpIdInBlock][36]; }
                if(laneId<9){
                    newEmb[(maxRowNum+index)*9+laneId] = writeEle;
                }
                writeEle = neigV[warpIdInBlock][37];
                if(laneId==0){
                    atomicAdd(indexForIndex+256*2+66,writeEle);
                }
                __syncthreads();
                if(threadIdx.x==0){
                    atomicAdd(totWriteRowNum+1,indexForIndex[256*2+66]);
                }
                break;
            }
        }else{
            genExtEmb_2V_2src_1and2_noRestrict_sameLabel_evalEq_withPhase<false>(neigAddr1+indexNum1,neigAddr2+indexNum2,neigLen1,neigLen2,neigV[warpIdInBlock],indexForIndex+2*256+1+1+6,totWriteRowNum,newEmb,laneId,embLen,maxRowNum);
            if(neigV[warpIdInBlock][66]==maxRowNum){
                index = (blockIdx.x*blockDim.x+threadIdx.x)>>5;
                uint writeEle = neigAddr1+indexNum1-neighborsData;
                if(laneId==1){ writeEle = neigLen1; }
                else if(laneId==2){ writeEle = neigAddr2+indexNum2-neighborsData; }
                else if(laneId==3){ writeEle = neigLen2; }
                else if(laneId==4){ writeEle = i; }
                else if(laneId==5){ writeEle = neigV[warpIdInBlock][64];}
                newEmb[(maxRowNum+index)*38+laneId] = neigV[warpIdInBlock][32+laneId];
                if(laneId<6){
                    newEmb[(maxRowNum+index)*38+32+laneId] = writeEle;
                }
                writeEle = neigV[warpIdInBlock][65];
                if(laneId==0){
                    atomicAdd(indexForIndex+256*2+66,writeEle);
                }
                __syncthreads();
                if(threadIdx.x==0){
                    atomicAdd(totWriteRowNum+1,indexForIndex[256*2+66]);
                }
                break;
            }
        }

        i = i + gridWarpNum;
    }
}

//this is modified
template<bool isExt1Restrict,bool isExt2Restrict,bool isRecord1,bool isRecord2,bool useRecord1,bool useRecord2>
__global__ void extEmb_2V_2src_1and2_notSameLabel_NoHash_kernel(uint *vLabel, uint *totWriteRowNum, uint *edgeLabelPartition, 
    uint *neighborsData, uint *recordPos1, uint *recordPos2, uint *auxArray, uint *newEmb, uint *partialEmb, uint intervalNum, 
    uint partialRowNum, uint embLen,bool isContinue, uint maxRowNum) {

    uint laneId, gridWarpNum, warpIdInBlock, i, j;
    INIT_VARS_SHAREDMEM(96)
    __syncthreads();
    if(isContinue){
        uint writeEle = laneId<6?newEmb[(maxRowNum+i)*38+laneId]:0;
        neigV[warpIdInBlock][laneId+64] = writeEle; 
        neigV[warpIdInBlock][32+laneId] = newEmb[(maxRowNum+i)*38+6+laneId];
        i = neigV[warpIdInBlock][64+4];
        j = laneId<(embLen-2)?partialEmb[i*(embLen-2)+laneId]:1;
        neigV[warpIdInBlock][laneId] = j;
        uint *neigAddr1 = neighborsData+neigV[warpIdInBlock][64+0];
        uint neigLen1 = neighborsData[64+1];
        uint *neigAddr2 = neighborsData+neigV[warpIdInBlock][64+2];
        uint neigLen2 = neighborsData[64+3];
        genExtEmb_2V_1and2_notSameLabel_withPhase<true>(neigAddr1,neigAddr2,neigLen1,neigLen2,neigV[warpIdInBlock],indexForIndex+2*256+1+1+6,totWriteRowNum,newEmb,laneId,embLen,maxRowNum);
        i = i + gridWarpNum;
    }
    while (i < partialRowNum) {
        j = laneId<(embLen-2)?partialEmb[i*(embLen-2)+laneId]:1;
        neigV[warpIdInBlock][laneId] = j;
        j = __ballot_sync(0xffffffff,j);
        if(j<0xffffffff) { i=i+gridWarpNum; continue; }

        uint svid1 = indexForIndex[256*2+1+1+4];
        uint svid2 = indexForIndex[256*2+1+1+5];
        svid1 = neigV[warpIdInBlock][svid1];
        svid2 = neigV[warpIdInBlock][svid2];
        uint lowerLimit = 0, upperLimit, predicate = 0, index;
        uint len = (intervalNum+31)&0xffffffe0;
        for(j=laneId;j<len;j=j+32){
            lowerLimit = indexForIndex[j];
            upperLimit = indexForIndex[j+intervalNum+1];
            upperLimit = lowerLimit+upperLimit-indexForIndex[j+intervalNum];
            predicate = j<intervalNum?(svid1>=lowerLimit && svid1<upperLimit):0;
            predicate = __ballot_sync(0xffffffff,predicate);
            if(predicate>0){
                index = svid1-lowerLimit+indexForIndex[j+intervalNum];
                uint tmpIndex = __ffs(predicate)-1;
                index = __shfl_sync(0xffffffff,index,tmpIndex);
                break;
            }
        }
        if(predicate==0){ 
            i=i+gridWarpNum;
            continue; 
        }
        uint baseNum = edgeLabelPartition[index];
        uint *neigAddr1 = neighborsData+baseNum;
        uint neigLen1 = edgeLabelPartition[index+1]-baseNum;
        uint indexNum1 = 0;
        if(useRecord1){
            uint distance = recordPos1[index*2];
            uint tmpNeigLen1 = recordPos1[index*2+1];
            if(distance==0 && tmpNeigLen1==0){
                uint evid1label = indexForIndex[2*256+1+1+0],found,lessThan,greatThan;
                lowerLimit = 0; upperLimit = neigLen1;
                FIND_LABEL_LIMIT(neigAddr1, lowerLimit, upperLimit, evid1label)
                neigLen1 = upperLimit - lowerLimit;
                if (neigLen1 == 0) {
                    if (isRecord1) { if (laneId == 0) { recordPos1[index * 2] = 1; recordPos1[index * 2 + 1] = 0; } }
                    i = i + gridWarpNum;
                    continue;
                }
                if (neigLen1 > VBLOCKSIZE+2 && neigLen1 <= 32 * VBLOCKSIZE + 32) { indexNum1 = (neigLen1+VBLOCKSIZE) / (VBLOCKSIZE + 1); }
                else if (neigLen1 > 32 * VBLOCKSIZE + 32) { indexNum1 = 32; }
                neigAddr1 = neigAddr1 + lowerLimit;
                neigLen1 = neigLen1 - indexNum1;
                if (isRecord1) { if (laneId == 0) { recordPos1[index * 2] = lowerLimit; recordPos1[index * 2 + 1] = neigLen1; } }
            }else if(distance==1 && tmpNeigLen1==0){
                i = i+gridWarpNum;
                continue;
            }else{
                neigAddr1 = neigAddr1+distance;
                neigLen1 = tmpNeigLen1;
            }
        }else {
            uint evid1label = indexForIndex[2 * 256+1+1+0],found,lessThan,greatThan;
            lowerLimit = 0; upperLimit = neigLen1;
            FIND_LABEL_LIMIT(neigAddr1, lowerLimit, upperLimit, evid1label)
            neigLen1 = upperLimit - lowerLimit;
            if (neigLen1 == 0) {
                if (isRecord1) { if (laneId == 0) { recordPos1[index * 2] = 1; recordPos1[index * 2 + 1] = 0; } }
                i = i + gridWarpNum;
                continue;
            }
            if (neigLen1 > VBLOCKSIZE+2 && neigLen1 <= 32 * VBLOCKSIZE + 32) { indexNum1 = (neigLen1+VBLOCKSIZE) / (VBLOCKSIZE + 1); }
            else if (neigLen1 > 32 * VBLOCKSIZE + 32) { indexNum1 = 32; }
            neigAddr1 = neigAddr1 + lowerLimit;
            neigLen1 = neigLen1 - indexNum1;
            if (isRecord1) { if (laneId == 0) { recordPos1[index * 2] = lowerLimit; recordPos1[index * 2 + 1] = neigLen1; } }
        }
        if(isExt1Restrict){
            uint lessvid = indexForIndex[2*256+1+1+2];
            lessvid = neigV[warpIdInBlock][lessvid];
            if (indexNum1 == 0) { neigLen1 = findLessVPosNoIndex(neigAddr1, neigLen1, lessvid, laneId); }
            else { neigLen1 = findLessVPosWithIndex(neigAddr1, neigLen1, lessvid, laneId, indexNum1); }
            if (neigLen1 == 0) { i=i+gridWarpNum; continue; }
        }

        for(j=laneId;j<len;j=j+32){
            lowerLimit = indexForIndex[j];
            upperLimit = indexForIndex[j+intervalNum+1];
            upperLimit = lowerLimit+upperLimit-indexForIndex[j+intervalNum];
            predicate = j<intervalNum?(svid2>=lowerLimit && svid2<upperLimit):0;
            predicate = __ballot_sync(0xffffffff,predicate);
            if(predicate>0){
                index = svid2-lowerLimit+indexForIndex[j+intervalNum];
                uint tmpIndex = __ffs(predicate)-1;
                index = __shfl_sync(0xffffffff,index,tmpIndex);
                break;
            }
        }
        if(predicate==0){ 
            i=i+gridWarpNum;
            continue; 
        }
        baseNum = edgeLabelPartition[index];
        uint *neigAddr2 = neighborsData+baseNum;
        uint neigLen2 = edgeLabelPartition[index+1]-baseNum;
        uint indexNum2 = 0;
        if(useRecord2){
            uint distance = recordPos2[index*2];
            uint tmpNeigLen2 = recordPos2[index*2+1];
            if(distance==0 && tmpNeigLen2==0){
                uint evid2label = indexForIndex[2 * 256+1+1+1],found,lessThan,greatThan;
                lowerLimit = 0; upperLimit = neigLen2;
                FIND_LABEL_LIMIT(neigAddr2, lowerLimit, upperLimit, evid2label)
                neigLen2 = upperLimit - lowerLimit;
                if (neigLen2 == 0) {
                    if (isRecord2) { if (laneId == 0) { recordPos2[index * 2] = 1; recordPos2[index * 2 + 1] = 0; } }
                    i = i + gridWarpNum;
                    continue;
                }
                if (neigLen2 > VBLOCKSIZE+2 && neigLen2 <= 32 * VBLOCKSIZE + 32) { indexNum2 = (neigLen2+VBLOCKSIZE) / (VBLOCKSIZE + 1); }
                else if (neigLen2 > 32 * VBLOCKSIZE + 32) { indexNum2 = 32; }
                neigAddr2 = neigAddr2 + lowerLimit;
                neigLen2 = neigLen2 - indexNum2;
                if (isRecord2) { if (laneId == 0) { recordPos2[index * 2] = lowerLimit; recordPos2[index * 2 + 1] = neigLen2; } }
            }else if(distance==1 && tmpNeigLen2==0){
                i = i+gridWarpNum;
                continue;
            }else{
                neigAddr2 = neigAddr2+distance;
                neigLen2 = tmpNeigLen2;
            }
        }else {
            uint evid2label = indexForIndex[2*256+1+1+1],found,lessThan,greatThan;
            lowerLimit = 0; upperLimit = neigLen2;
            FIND_LABEL_LIMIT(neigAddr2, lowerLimit, upperLimit, evid2label)
            neigLen2 = upperLimit - lowerLimit;
            if (neigLen2 == 0) {
                if (isRecord2) { if (laneId == 0) { recordPos2[index * 2] = 1; recordPos2[index * 2 + 1] = 0; } }
                i = i + gridWarpNum;
                continue;
            }
            if (neigLen2 > VBLOCKSIZE+2 && neigLen2 <= 32 * VBLOCKSIZE + 32) { indexNum2 = (neigLen2+VBLOCKSIZE) / (VBLOCKSIZE + 1); }
            else if (neigLen2 > 32 * VBLOCKSIZE + 32) { indexNum2 = 32; }
            neigAddr2 = neigAddr2 + lowerLimit;
            neigLen2 = neigLen2 - indexNum2;
            if (isRecord2) { if (laneId == 0) { recordPos2[index * 2] = lowerLimit; recordPos2[index * 2 + 1] = neigLen2; } }
        }
        if(isExt2Restrict){
            uint lessvid = indexForIndex[2*256+1+1+3];
            lessvid = neigV[warpIdInBlock][lessvid];
            if (indexNum2 == 0) { neigLen2 = findLessVPosNoIndex(neigAddr2, neigLen2, lessvid, laneId); }
            else { neigLen2 = findLessVPosWithIndex(neigAddr2, neigLen2, lessvid, laneId, indexNum2); }
            if (neigLen2 == 0) { i=i+gridWarpNum; continue; }
        }
        genExtEmb_2V_1and2_notSameLabel_withPhase<false>(neigAddr1+indexNum1,neigAddr2+indexNum2,neigLen1,neigLen2,neigV[warpIdInBlock],indexForIndex+2*256+1+1+6,totWriteRowNum,newEmb,laneId,embLen,maxRowNum);
        if(neigV[warpIdInBlock][66]==maxRowNum){
            index = (blockIdx.x*blockDim.x+threadIdx.x)>>5;
            uint writeEle = neigAddr1+indexNum1-neighborsData;
            if(laneId==1){ writeEle = neigLen1; }
            else if(laneId==2){ writeEle = neigAddr2+indexNum2-neighborsData; }
            else if(laneId==3){ writeEle = neigLen2; }
            else if(laneId==4){ writeEle = i; }
            else if(laneId==5){ writeEle = neigV[warpIdInBlock][64]; }
            if(laneId<6){
                newEmb[(maxRowNum+index)*38+laneId] = writeEle;
            }
            newEmb[(maxRowNum+index)*38+6+laneId] = neigV[warpIdInBlock][32+laneId];
            writeEle = neigV[warpIdInBlock][65];
            if(laneId==0){
                atomicAdd(indexForIndex+256*2+66,writeEle);
            }
            __syncthreads();
            if(threadIdx.x==0){
                atomicAdd(totWriteRowNum+1,indexForIndex[256*2+66]);
            }
            break;
        }
        i = i + gridWarpNum;
    }
}

//this is modified
template<bool isExt1Restrict, bool isExt2Restrict, bool isExt1and2Restrict,bool isRecord, bool useRecord>
__global__ void extEmb_2V_1src_1and2_sameLabel_NoHash_kernel(uint *vLabel, uint *totWriteRowNum, 
    uint *edgeLabelPartition, uint *neighborsData, uint *recordPos, uint *auxArray, uint *newEmb, uint *partialEmb, uint intervalNum, 
    uint partialRowNum, uint embLen, bool isContinue, uint maxRowNum) {

    uint laneId, gridWarpNum, warpIdInBlock, i, j;
    INIT_VARS_SHAREDMEM(96)
    __syncthreads();
    if(isContinue){
        uint writeEle = laneId<6?newEmb[(maxRowNum+i)*38+laneId]:0;
        neigV[warpIdInBlock][laneId+64] = writeEle; 
        neigV[warpIdInBlock][32+laneId] = newEmb[(maxRowNum+i)*38+6+laneId];
        i = neigV[warpIdInBlock][64+4];
        j = laneId<(embLen-2)?partialEmb[i*(embLen-2)+laneId]:1;
        neigV[warpIdInBlock][laneId] = j;
        uint *neigAddr1 = neighborsData+neigV[warpIdInBlock][64+0];
        uint neigLen1 = neighborsData[64+1];
        uint *neigAddr2 = neighborsData+neigV[warpIdInBlock][64+2];
        uint neigLen2 = neighborsData[64+3];
        if(isExt1and2Restrict){
            genExtEmb_2V_1src_1and2_restrict_withPhase<true>(neigAddr1,neigAddr2,neigLen1,neigLen2,neigV[warpIdInBlock],indexForIndex+2*256+1+1+6,totWriteRowNum,newEmb,laneId,embLen,maxRowNum);
        }else{
            genExtEmb_2V_1src_1and2_noRestrict_sameLabel_evalEq_withPhase<true>(neigAddr1,neigAddr2,neigLen1,neigLen2,neigV[warpIdInBlock],indexForIndex+2*256+1+1+6,totWriteRowNum,newEmb,laneId,embLen,maxRowNum);
        }
        i = i + gridWarpNum;
    }
    while (i < partialRowNum) {
        j = laneId<(embLen-2)?partialEmb[i*(embLen-2)+laneId]:1;
        neigV[warpIdInBlock][laneId] = j;
        j = __ballot_sync(0xffffffff,j);
        if(j<0xffffffff) { i=i+gridWarpNum; continue; }

        uint svid = indexForIndex[256*2+1+1+4];
        svid = neigV[warpIdInBlock][svid];
        uint lowerLimit = 0, upperLimit, predicate = 0, index;
        uint len = (intervalNum+31)&0xffffffe0;
        for(j=laneId;j<len;j=j+32){
            lowerLimit = indexForIndex[j];
            upperLimit = indexForIndex[j+intervalNum+1];
            upperLimit = lowerLimit+upperLimit-indexForIndex[j+intervalNum];
            predicate = j<intervalNum?(svid>=lowerLimit && svid<upperLimit):0;
            predicate = __ballot_sync(0xffffffff,predicate);
            if(predicate>0){
                index = svid-lowerLimit+indexForIndex[j+intervalNum];
                uint tmpIndex = __ffs(predicate)-1;
                index = __shfl_sync(0xffffffff,index,tmpIndex);
                break;
            }
        }
        if(predicate==0){ 
            i=i+gridWarpNum;
            continue; 
        }
        uint baseNum = edgeLabelPartition[index];
        uint *neigAddr = neighborsData+baseNum;
        uint neigLen = edgeLabelPartition[index+1]-baseNum;
        uint indexNum = 0;
        if(useRecord){
            uint distance = recordPos[index*2];
            uint tmpNeigLen = recordPos[index*2+1];
            if(distance==0 && tmpNeigLen==0){
                uint evidlabel = indexForIndex[2*256+1+1+0],found,lessThan,greatThan;
                lowerLimit = 0; upperLimit = neigLen;
                FIND_LABEL_LIMIT(neigAddr, lowerLimit, upperLimit, evidlabel)
                neigLen = upperLimit - lowerLimit;
                if (neigLen == 0) {
                    if (isRecord) { if (laneId == 0) { recordPos[index*2] = 1; recordPos[index*2+1] = 0; } }
                    i = i + gridWarpNum;
                    continue;
                }
                if (neigLen > VBLOCKSIZE+2 && neigLen <= 32 * VBLOCKSIZE + 32) { indexNum = (neigLen+VBLOCKSIZE) / (VBLOCKSIZE + 1); }
                else if (neigLen > 32 * VBLOCKSIZE + 32) { indexNum = 32; }
                neigAddr = neigAddr + lowerLimit;
                neigLen = neigLen - indexNum;
                if (isRecord) { if (laneId == 0) { recordPos[index * 2] = lowerLimit; recordPos[index * 2 + 1] = neigLen; } }
            }else if(distance==1 && tmpNeigLen==0){
                i = i+gridWarpNum;
                continue;
            }else{
                neigAddr = neigAddr+distance;
                neigLen = tmpNeigLen;
            }
        }else {
            uint evidlabel = indexForIndex[2*256+1+1+0],found,lessThan,greatThan;
            lowerLimit = 0; upperLimit = neigLen;
            FIND_LABEL_LIMIT(neigAddr, lowerLimit, upperLimit, evidlabel)
            neigLen = upperLimit - lowerLimit;
            if (neigLen == 0) {
                if (isRecord) { if (laneId == 0) { recordPos[index*2] = 1; recordPos[index*2+1] = 0; } }
                i = i + gridWarpNum;
                continue;
            }
            if (neigLen > VBLOCKSIZE+2 && neigLen <= 32 * VBLOCKSIZE + 32) { indexNum = (neigLen+VBLOCKSIZE) / (VBLOCKSIZE + 1); }
            else if (neigLen > 32 * VBLOCKSIZE + 32) { indexNum = 32; }
            neigAddr = neigAddr + lowerLimit;
            neigLen = neigLen - indexNum;
            if (isRecord) { if (laneId == 0) { recordPos[index * 2] = lowerLimit; recordPos[index * 2 + 1] = neigLen; } }
        }
        uint neigLen1 = neigLen, neigLen2 = neigLen;
        if(isExt1Restrict){
            uint lessvid = indexForIndex[2*256+1+1+2];
            lessvid = neigV[warpIdInBlock][lessvid];
            if (indexNum == 0) { neigLen1 = findLessVPosNoIndex(neigAddr, neigLen, lessvid, laneId); }
            else { neigLen1 = findLessVPosWithIndex(neigAddr, neigLen, lessvid, laneId, indexNum); }
            if (neigLen1 == 0) { 
                i=i+gridWarpNum;
                continue; 
            }
        }
        if(isExt2Restrict){
            uint lessvid = indexForIndex[2*256+1+1+3];
            lessvid = neigV[warpIdInBlock][lessvid];
            if (indexNum == 0) { neigLen2 = findLessVPosNoIndex(neigAddr, neigLen, lessvid, laneId); }
            else { neigLen2 = findLessVPosWithIndex(neigAddr, neigLen, lessvid, laneId, indexNum); }
            if (neigLen2 == 0) { 
                i=i+gridWarpNum;
                continue; 
            }
        }

        if(isExt1and2Restrict){
            if(isContinue){
                genExtEmb_2V_1src_1and2_restrict_withPhase<true>(neigAddr+indexNum,neigAddr+indexNum,neigLen1,neigLen2,neigV[warpIdInBlock],indexForIndex+2*256+1+1+6,totWriteRowNum,newEmb,laneId,embLen,maxRowNum);
            }else{
                genExtEmb_2V_1src_1and2_restrict_withPhase<false>(neigAddr+indexNum,neigAddr+indexNum,neigLen1,neigLen2,neigV[warpIdInBlock],indexForIndex+2*256+1+1+6,totWriteRowNum,newEmb,laneId,embLen,maxRowNum);
            }
        }else{
            if(isContinue){
                genExtEmb_2V_1src_1and2_noRestrict_sameLabel_evalEq_withPhase<true>(neigAddr+indexNum,neigAddr+indexNum,neigLen1,neigLen2,neigV[warpIdInBlock],indexForIndex+2*256+1+1+6,totWriteRowNum,newEmb,laneId,embLen,maxRowNum);
            }else{
                genExtEmb_2V_1src_1and2_noRestrict_sameLabel_evalEq_withPhase<false>(neigAddr+indexNum,neigAddr+indexNum,neigLen1,neigLen2,neigV[warpIdInBlock],indexForIndex+2*256+1+1+6,totWriteRowNum,newEmb,laneId,embLen,maxRowNum);
            }
        }
        if(neigV[warpIdInBlock][66]==maxRowNum){
            index = (blockIdx.x*blockDim.x+threadIdx.x)>>5;
            uint writeEle = neigAddr+indexNum-neighborsData;
            if(laneId==1){ writeEle = neigLen1; }
            else if(laneId==2){ writeEle = neigAddr+indexNum-neighborsData; }
            else if(laneId==3){ writeEle = neigLen2; }
            else if(laneId==4){ writeEle = i; }
            else if(laneId==5){ writeEle = neigV[warpIdInBlock][64]; }
            if(laneId<6){
                newEmb[(maxRowNum+index)*38+laneId] = writeEle;
            }
            newEmb[(maxRowNum+index)*38+6+laneId] = neigV[warpIdInBlock][32+laneId];
            writeEle = neigV[warpIdInBlock][65];
            if(laneId==0){
                atomicAdd(indexForIndex+256*2+66,writeEle);
            }
            __syncthreads();
            if(threadIdx.x==0){
                atomicAdd(totWriteRowNum+1,indexForIndex[256*2+66]);
            }
            break;
        }
        i = i + gridWarpNum;
    }
}

//this is modified
template<bool isExt1Restrict,bool isExt2Restrict,bool isRecord1,bool isRecord2,bool useRecord1,bool useRecord2>
__global__ void extEmb_2V_1src_1and2_notSameLabel_NoHash_kernel(uint *vLabel, uint *totWriteRowNum, uint *edgeLabelPartition, 
    uint *neighborsData, uint *recordPos1, uint *recordPos2, uint *auxArray, uint *newEmb, uint *partialEmb, uint intervalNum, 
    uint partialRowNum, uint embLen,bool isContinue, uint maxRowNum){

	uint laneId, gridWarpNum, warpIdInBlock, i, j;
    INIT_VARS_SHAREDMEM(96)
    __syncthreads();
    if(isContinue){
        uint writeEle = laneId<6?newEmb[(maxRowNum+i)*38+laneId]:0;
        neigV[warpIdInBlock][laneId+64] = writeEle; 
        neigV[warpIdInBlock][32+laneId] = newEmb[(maxRowNum+i)*38+6+laneId];
        i = neigV[warpIdInBlock][64+4];
        j = laneId<(embLen-2)?partialEmb[i*(embLen-2)+laneId]:1;
        neigV[warpIdInBlock][laneId] = j;
        uint *neigAddr1 = neighborsData+neigV[warpIdInBlock][64+0];
        uint neigLen1 = neighborsData[64+1];
        uint *neigAddr2 = neighborsData+neigV[warpIdInBlock][64+2];
        uint neigLen2 = neighborsData[64+3];
        genExtEmb_2V_1and2_notSameLabel_withPhase<true>(neigAddr1,neigAddr2,neigLen1,neigLen2,neigV[warpIdInBlock],indexForIndex+2*256+1+1+6,totWriteRowNum,newEmb,laneId,embLen,maxRowNum);
        i = i + gridWarpNum;
    }

    while (i < partialRowNum) {
        j = laneId<(embLen-2)?partialEmb[i*(embLen-2)+laneId]:1;
        neigV[warpIdInBlock][laneId] = j;
        j = __ballot_sync(0xffffffff,j);
        if(j<0xffffffff){ i=i+gridWarpNum; continue; }

        uint svid = indexForIndex[256*2+1+1+4];
        svid = neigV[warpIdInBlock][svid];
        uint lowerLimit = 0, upperLimit, predicate = 0, index;
        uint len = (intervalNum+31)&0xffffffe0;
        for(j=laneId;j<len;j=j+32){
            lowerLimit = indexForIndex[j];
            upperLimit = indexForIndex[j+intervalNum+1];
            upperLimit = lowerLimit+upperLimit-indexForIndex[j+intervalNum];
            predicate = j<intervalNum?(svid>=lowerLimit && svid<upperLimit):0;
            predicate = __ballot_sync(0xffffffff,predicate);
            if(predicate>0){
                index = svid-lowerLimit+indexForIndex[j+intervalNum];
                uint tmpIndex = __ffs(predicate)-1;
                index = __shfl_sync(0xffffffff,index,tmpIndex);
                break;
            }
        }
        if(predicate==0){ 
            i=i+gridWarpNum;
            continue; 
        }
        uint baseNum = edgeLabelPartition[index];
        uint *neigAddr = neighborsData+baseNum;
        uint neigLen = edgeLabelPartition[index+1]-baseNum;
        uint *neigAddr1 = neigAddr, neigLen1 = neigLen, indexNum1 = 0;
        if(useRecord1){
            uint distance = recordPos1[index*2];
            //uint *tmpNeigAddr = neigAddr + recordPos1[index*2];
            uint tmpNeigLen = recordPos1[index*2+1];
            if(distance==0 && tmpNeigLen==0){
                uint evid1label = indexForIndex[2*256+1+1+0],found,lessThan,greatThan;
                lowerLimit = 0; upperLimit = neigLen;
                FIND_LABEL_LIMIT(neigAddr, lowerLimit, upperLimit, evid1label)
                neigLen1 = upperLimit - lowerLimit;
                if (neigLen1 == 0) {
                    if (isRecord1) { if (laneId == 0) { recordPos1[index*2] = 1; recordPos1[index * 2 + 1] = 0; } }
                    i = i + gridWarpNum;
                    continue;
                }
                if (neigLen1 > VBLOCKSIZE+2 && neigLen1 <= 32 * VBLOCKSIZE + 32) { indexNum1 = (neigLen1+VBLOCKSIZE) / (VBLOCKSIZE + 1); }
                else if (neigLen1 > 32 * VBLOCKSIZE + 32) { indexNum1 = 32; }
                neigAddr1 = neigAddr + lowerLimit;
                neigLen1 = neigLen1 - indexNum1;
                if (isRecord1) { if (laneId == 0) { recordPos1[index * 2] = lowerLimit; recordPos1[index * 2 + 1] = neigLen1; } }
            }else if(distance==1 && tmpNeigLen==0){
                i = i+gridWarpNum;
                continue;
            }else{
                neigAddr1 = neigAddr+distance;
                neigLen1 = tmpNeigLen;
            }
        }else {
            uint evid1label = indexForIndex[2*256+1+1+0],found,lessThan,greatThan;
            lowerLimit = 0; upperLimit = neigLen;
            FIND_LABEL_LIMIT(neigAddr, lowerLimit, upperLimit, evid1label)
            neigLen1 = upperLimit - lowerLimit;
            if (neigLen1 == 0) {
                if (isRecord1) { if (laneId == 0) { recordPos1[index * 2] = 1; recordPos1[index * 2 + 1] = 0; } }
                i = i + gridWarpNum;
                continue;
            }
            if (neigLen1 > VBLOCKSIZE+2 && neigLen1 <= 32 * VBLOCKSIZE + 32) { indexNum1 = (neigLen1+VBLOCKSIZE) / (VBLOCKSIZE + 1); }
            else if (neigLen1 > 32 * VBLOCKSIZE + 32) { indexNum1 = 32; }
            neigAddr1 = neigAddr + lowerLimit;
            neigLen1 = neigLen1 - indexNum1;
            if (isRecord1) { if (laneId == 0) { recordPos1[index * 2] = lowerLimit; recordPos1[index * 2 + 1] = neigLen1; } }
        }
        if(isExt1Restrict){
            uint lessvid = indexForIndex[2*256+1+1+2];
            lessvid = neigV[warpIdInBlock][lessvid];
            if (indexNum1 == 0) { neigLen1 = findLessVPosNoIndex(neigAddr1, neigLen1, lessvid, laneId); }
            else { neigLen1 = findLessVPosWithIndex(neigAddr1, neigLen1, lessvid, laneId, indexNum1); }
            if (neigLen1 == 0) { i=i+gridWarpNum; continue; }
        }

        uint *neigAddr2 = neigAddr, neigLen2 = neigLen, indexNum2 = 0;
        if(useRecord2){
            uint distance = recordPos2[index*2];
            //uint *tmpNeigAddr2 = neigAddr2 + recordPos2[index*2];
            uint tmpNeigLen2 = recordPos2[index*2+1];
            if(distance==0 && tmpNeigLen2==0){
                uint evid2label = indexForIndex[2*256+1+1+1],found,lessThan,greatThan;
                //evid2 is always > evid1label
                lowerLimit = neigAddr1+neigLen1-neigAddr; 
                upperLimit = neigLen2;
                FIND_LABEL_LIMIT(neigAddr, lowerLimit, upperLimit, evid2label)
                neigLen2 = upperLimit - lowerLimit;
                if (neigLen2 == 0) {
                    if (isRecord2) { if (laneId == 0) { recordPos2[index * 2] = 1; recordPos2[index * 2 + 1] = 0; } }
                    i = i + gridWarpNum;
                    continue;
                }
                if (neigLen2 > VBLOCKSIZE+2 && neigLen2 <= 32 * VBLOCKSIZE + 32) { indexNum2 = (neigLen2+VBLOCKSIZE) / (VBLOCKSIZE + 1); }
                else if (neigLen2 > 32 * VBLOCKSIZE + 32) { indexNum2 = 32; }
                neigAddr2 = neigAddr2 + lowerLimit;
                neigLen2 = neigLen2 - indexNum2;
                if (isRecord2) { if (laneId == 0) { recordPos2[index * 2] = lowerLimit; recordPos2[index * 2 + 1] = neigLen2; } }
            }else if(distance==1 && tmpNeigLen2==0){
                i = i+gridWarpNum;
                continue;
            }else{
                neigAddr2 = neigAddr+distance;
                neigLen2 = tmpNeigLen2;
            }
        }else {
            uint evid2label = indexForIndex[2 * 256+1+1+1],found,lessThan,greatThan;
            lowerLimit = neigAddr1+neigLen1-neigAddr; 
            upperLimit = neigLen2;
            FIND_LABEL_LIMIT(neigAddr2, lowerLimit, upperLimit, evid2label)
            neigLen2 = upperLimit - lowerLimit;
            if (neigLen2 == 0) {
                if (isRecord2) { if (laneId == 0) { recordPos2[index * 2] = 1; recordPos2[index * 2 + 1] = 0; } }
                i = i + gridWarpNum;
                continue;
            }
            if (neigLen2 > VBLOCKSIZE+2 && neigLen2 <= 32 * VBLOCKSIZE + 32) { indexNum2 = (neigLen2+VBLOCKSIZE) / (VBLOCKSIZE + 1); }
            else if (neigLen2 > 32 * VBLOCKSIZE + 32) { indexNum2 = 32; }
            neigAddr2 = neigAddr2 + lowerLimit;
            neigLen2 = neigLen2 - indexNum2;
            if (isRecord2) { if (laneId == 0) { recordPos2[index * 2] = lowerLimit; recordPos2[index * 2 + 1] = neigLen2; } }
        }
        if(isExt2Restrict){
            uint lessvid = indexForIndex[2*256+1+1+3];
            lessvid = neigV[warpIdInBlock][lessvid];
            if (indexNum2 == 0) { neigLen2 = findLessVPosNoIndex(neigAddr2, neigLen2, lessvid, laneId); }
            else { neigLen2 = findLessVPosWithIndex(neigAddr2, neigLen2, lessvid, laneId, indexNum2); }
            if (neigLen2 == 0) { i = i + gridWarpNum;continue; }
        }
        genExtEmb_2V_1and2_notSameLabel_withPhase<false>(neigAddr1+indexNum1,neigAddr2+indexNum2,neigLen1,neigLen2,neigV[warpIdInBlock],indexForIndex+2*256+1+1+6,totWriteRowNum,newEmb,laneId,embLen,maxRowNum);
        if(neigV[warpIdInBlock][66]==maxRowNum){
            index = (blockIdx.x*blockDim.x+threadIdx.x)>>5;
            uint writeEle = neigAddr1+indexNum1-neighborsData;
            if(laneId==1){ writeEle = neigLen1; }
            else if(laneId==2){ writeEle = neigAddr2+indexNum2-neighborsData; }
            else if(laneId==3){ writeEle = neigLen2; }
            else if(laneId==4){ writeEle = i; }
            else if(laneId==5){ writeEle = neigV[warpIdInBlock][64]; }
            if(laneId<6){
                newEmb[(maxRowNum+index)*38+laneId] = writeEle;
            }
            newEmb[(maxRowNum+index)*38+6+laneId] = neigV[warpIdInBlock][32+laneId];
            writeEle = neigV[warpIdInBlock][65];
            if(laneId==0){
                atomicAdd(indexForIndex+256*2+66,writeEle);
            }
            __syncthreads();
            if(threadIdx.x==0){
                atomicAdd(totWriteRowNum+1,indexForIndex[256*2+66]);
            }
            break;
        }
        i = i + gridWarpNum;
    }
}

//this is modified
template<bool isExtRestrict,bool isRecord,bool useRecord, bool isLastPhase>
__global__ void extEmb_1V_NoHash_kernel(uint *vLabel, uint *totWriteRowNum, uint *edgeLabelPartition, 
    uint *neighborsData, uint *recordPos, uint *auxArray, uint *newEmb, uint *partialEmb, uint intervalNum, uint partialRowNum, 
    uint embLen, bool isContinue, uint maxRowNum) {

	uint laneId, gridWarpNum, warpIdInBlock, i, j;
    INIT_VARS_SHAREDMEM(64)
	__syncthreads();

    if(isContinue){
        uint writeEle = laneId<4?newEmb[(maxRowNum+i)*4+laneId]:0;
        neigV[warpIdInBlock][laneId+32] = writeEle; 
        i = neigV[warpIdInBlock][32+2];
        j = laneId<(embLen-1)?partialEmb[i*(embLen-1)+laneId]:1;
        neigV[warpIdInBlock][laneId] = j;
        uint *neigAddr = neighborsData+neigV[warpIdInBlock][32];
        uint neigLen = neigV[warpIdInBlock][32+1];
        genExtEmb_1V<isLastPhase,true>(neigAddr,neigLen,neigV[warpIdInBlock],indexForIndex+256*2+1+1+6,totWriteRowNum,newEmb,laneId,embLen,maxRowNum);
        /*if(retval>0){
            index = (blockIDx.x*blockDim.x+threadIdx.x)>>5;
            uint writeEle = neigAddr-neighborsData;
            if(laneId==1){ writeEle = neigLen; }
            else if(laneId==2){ writeEle = i; }
            else if(laneId==3){ writeEle = retval; }
            if(laneId<4){
                newEmb[(maxRowNum+index)*4+laneId] = writeEle;
            }
            break;
        }*/
        i = i + gridWarpNum;
    }
    while (i < partialRowNum) {
        j = laneId<(embLen-1)?partialEmb[i*(embLen-1)+laneId]:1;
        neigV[warpIdInBlock][laneId] = j;
        j = __ballot_sync(0xffffffff,j);
        if(j<0xffffffff){ i = i + gridWarpNum; continue; }

        uint svid = indexForIndex[256*2+1+1+4];
        svid = neigV[warpIdInBlock][svid];

        uint lowerLimit = 0, upperLimit, predicate = 0, index;
        uint len = (intervalNum+31)&0xffffffe0;
        for(j=laneId;j<len;j=j+32){
            lowerLimit = indexForIndex[j];
            upperLimit = indexForIndex[j+intervalNum+1];
            upperLimit = lowerLimit+upperLimit-indexForIndex[j+intervalNum];
            predicate = j<intervalNum?(svid>=lowerLimit && svid<upperLimit):0;
            predicate = __ballot_sync(0xffffffff,predicate);
            if(predicate>0){
                index = svid-lowerLimit+indexForIndex[j+intervalNum];
                uint tmpIndex = __ffs(predicate)-1;
                index = __shfl_sync(0xffffffff,index,tmpIndex);
                break;
            }
        }
        if(predicate==0){ 
            i=i+gridWarpNum;
            continue; 
        }
        uint baseNum = edgeLabelPartition[index];
        uint *neigAddr = neighborsData+baseNum;
        uint neigLen = edgeLabelPartition[index+1]-baseNum;
        uint *neigAddr1 = neigAddr, neigLen1 = neigLen, indexNum1 = 0;
        if(useRecord){
            uint distance = recordPos[index*2];
            //uint *tmpNeigAddr = neigAddr + recordPos[index*2];
            uint tmpNeigLen = recordPos[index*2+1];
            if(distance==0 && tmpNeigLen==0){
                uint evidlabel = indexForIndex[2 * 256+1+1+0],found,lessThan,greatThan;
                lowerLimit = 0; upperLimit = neigLen;
                FIND_LABEL_LIMIT(neigAddr, lowerLimit, upperLimit, evidlabel)
                neigLen1 = upperLimit - lowerLimit;
                if (neigLen1 == 0) {
                    if (isRecord) { if (laneId == 0) { recordPos[index * 2] = 1; recordPos[index * 2 + 1] = 0; } }
                    i = i + gridWarpNum;
                    continue;
                }
                if (neigLen1 > VBLOCKSIZE+2 && neigLen1 <= 32 * VBLOCKSIZE + 32) { indexNum1 = (neigLen1+VBLOCKSIZE) / (VBLOCKSIZE + 1); }
                else if (neigLen1 > 32 * VBLOCKSIZE + 32) { indexNum1 = 32; }
                neigAddr1 = neigAddr + lowerLimit;
                neigLen1 = neigLen1 - indexNum1;
                if (isRecord) { if (laneId == 0) { recordPos[index * 2] = lowerLimit; recordPos[index * 2 + 1] = neigLen1; } }
            }else if(distance==1 && tmpNeigLen==0){
                i = i+gridWarpNum;
                continue;
            }else{
                neigAddr1 = neigAddr+distance;
                neigLen1 = tmpNeigLen;
            }
        }else {
            uint evidlabel = indexForIndex[2 * 256+1+1+0],found,lessThan,greatThan;
            lowerLimit = 0; upperLimit = neigLen;
			FIND_LABEL_LIMIT(neigAddr, lowerLimit, upperLimit, evidlabel)
            neigLen1 = upperLimit- lowerLimit;
            if (neigLen1 == 0) {
                if (isRecord) { if (laneId == 0) { recordPos[index * 2] = 1; recordPos[index * 2 + 1] = 0; } }
                i = i + gridWarpNum;
                continue;
            }
            if (neigLen1 > VBLOCKSIZE+2 && neigLen1 <= 32 * VBLOCKSIZE + 32) { indexNum1 = (neigLen1+VBLOCKSIZE) / (VBLOCKSIZE + 1); }
            else if (neigLen1 > 32 * VBLOCKSIZE + 32) { indexNum1 = 32; }
            neigAddr1 = neigAddr + lowerLimit;
            neigLen1 = neigLen1 - indexNum1;
            if (isRecord) { if (laneId == 0) { recordPos[index * 2] = lowerLimit; recordPos[index * 2 + 1] = neigLen1; } }
        }
        if(isExtRestrict){
            uint lesspos = indexForIndex[256*2+1+1+2];
            uint lessvid = neigV[warpIdInBlock][lesspos];
            if (indexNum1 == 0) { neigLen1 = findLessVPosNoIndex(neigAddr1, neigLen1, lessvid, laneId); }
            else { neigLen1 = findLessVPosWithIndex(neigAddr1, neigLen1, lessvid, laneId, indexNum1); }
            if (neigLen1 == 0) { i = i + gridWarpNum; continue; }
        }
        genExtEmb_1V<isLastPhase,false>(neigAddr1+indexNum1,neigLen1,neigV[warpIdInBlock],indexForIndex+256*2+1+1+6,totWriteRowNum,newEmb,laneId,embLen,maxRowNum);
        if(neigV[warpIdInBlock][34]==maxRowNum){
            index = (blockIdx.x*blockDim.x+threadIdx.x)>>5;
            uint writeEle = neigAddr1+indexNum1-neighborsData;
            if(laneId==1){ writeEle = neigLen1; }
            else if(laneId==2){ writeEle = i; }
            else if(laneId==3){ writeEle = neigV[warpIdInBlock][32]; }
            if(laneId<4){
                newEmb[(maxRowNum+index)*4+laneId] = writeEle;
            }
            writeEle = neigV[warpIdInBlock][33];
            if(laneId==0){
                atomicAdd(indexForIndex+256*2+66,writeEle);
            }
            __syncthreads();
            if(threadIdx.x==0){
                atomicAdd(totWriteRowNum+1,indexForIndex[256*2+66]);
            }
            break;
        }
        i = i + gridWarpNum;
    }
}

extern int GPU_SM_NUM;

#define INITEMB_2V(isExtRestrict,isRecord) \
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocks,initEmb_2V_NoHash_kernel<isExtRestrict,isRecord>,WARPPERBLOCK * 32, 0);\
    initEmb_2V_NoHash_kernel<isExtRestrict,isRecord><<<numBlocks*GPU_SM_NUM,WARPPERBLOCK*32>>>(vLabel,writeRowNum,edgeLabelPartition,neighborsData,recordPos,newEmb,intervalNum,totSrcNum,svidlabel,evidlabel,isContinue,maxRowNum);


#define INITEMB_3V(isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isExt1and2SameLabel,isRecord1,isRecord2) \
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocks,initEmb_3V_1and2_NoHash_kernel<isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isExt1and2SameLabel,isRecord1,isRecord2>,WARPPERBLOCK * 32, 0);\
    initEmb_3V_1and2_NoHash_kernel<isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isExt1and2SameLabel,isRecord1,isRecord2><<<numBlocks*GPU_SM_NUM,WARPPERBLOCK*32>>>(vLabel,writeRowNum,edgeLabelPartition,neighborsData,recordPos1,recordPos2,newEmb,intervalNum,totSrcNum,svidlabel,evid1label,evid2label,isContinue,maxRowNum);

#define EXTEMB_1V(isExtRestrict,isRecord,useRecord,isLastPhase) \
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocks,extEmb_1V_NoHash_kernel<isExtRestrict,isRecord,useRecord,isLastPhase>,WARPPERBLOCK*32,0);\
    extEmb_1V_NoHash_kernel<isExtRestrict,isRecord,useRecord,isLastPhase><<<numBlocks*GPU_SM_NUM,WARPPERBLOCK*32>>>(vLabel,writeRowNum,edgeLabelPartition,neighborsData,recordPos,auxArray,newEmb,partialEmb,intervalNum,partialRowNum,embLen,isContinue,maxRowNum);\

#define EXTEMB_2V_2S_1L(isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isRecord,useRecord) \
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocks,extEmb_2V_2src_1and2_sameLabel_NoHash_kernel<isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isRecord,useRecord>,WARPPERBLOCK*32,0);\
    extEmb_2V_2src_1and2_sameLabel_NoHash_kernel<isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isRecord,useRecord><<<numBlocks*GPU_SM_NUM,WARPPERBLOCK*32>>>(vLabel,writeRowNum,edgeLabelPartition,neighborsData,recordPos,auxArray,newEmb,partialEmb,intervalNum,partialRowNum,embLen,isContinue,maxRowNum);

#define EXTEMB_2V_1S_2L(isExt1Restrict,isExt2Restrict,isRecord1,isRecord2,useRecord1,useRecord2) \
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocks,extEmb_2V_1src_1and2_notSameLabel_NoHash_kernel<isExt1Restrict,isExt2Restrict,isRecord1,isRecord2,useRecord1,useRecord2>,WARPPERBLOCK*32,0);\
    extEmb_2V_1src_1and2_notSameLabel_NoHash_kernel<isExt1Restrict,isExt2Restrict,isRecord1,isRecord2,useRecord1,useRecord2><<<numBlocks*GPU_SM_NUM,WARPPERBLOCK*32>>>(vLabel,writeRowNum,edgeLabelPartition,neighborsData,recordPos1,recordPos2,auxArray,newEmb,partialEmb,intervalNum,partialRowNum,embLen,isContinue,maxRowNum);

#define EXTEMB_2V_1S_1L(isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isRecord,useRecord) \
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocks,extEmb_2V_1src_1and2_sameLabel_NoHash_kernel<isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isRecord,useRecord>,WARPPERBLOCK*32,0);\
    extEmb_2V_1src_1and2_sameLabel_NoHash_kernel<isExt1Restrict,isExt2Restrict,isExt1and2Restrict,isRecord,useRecord><<<numBlocks*GPU_SM_NUM,WARPPERBLOCK*32>>>(vLabel,writeRowNum,edgeLabelPartition,neighborsData,recordPos,auxArray,newEmb,partialEmb,intervalNum,partialRowNum,embLen,isContinue,maxRowNum);

#define EXTEMB_2V_2S_2L(isExt1Restrict,isExt2Restrict,isRecord1,isRecord2,useRecord1,useRecord2) \
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocks,extEmb_2V_2src_1and2_notSameLabel_NoHash_kernel<isExt1Restrict,isExt2Restrict,isRecord1,isRecord2,useRecord1,useRecord2>,WARPPERBLOCK*32,0);\
    extEmb_2V_2src_1and2_notSameLabel_NoHash_kernel<isExt1Restrict,isExt2Restrict,isRecord1,isRecord2,useRecord1,useRecord2><<<numBlocks*GPU_SM_NUM,WARPPERBLOCK*32>>>(vLabel,writeRowNum,edgeLabelPartition,neighborsData,recordPos1,recordPos2,auxArray,newEmb,partialEmb,intervalNum,partialRowNum,embLen,isContinue,maxRowNum);


void initEmb_3V_1and2_NoHash(uint *vLabel,  uint *writeRowNum, uint *edgeLabelPartition, uint *neighborsData, uint *recordPos1,
    uint *recordPos2, uint *newEmb, uint intervalNum, uint totSrcNum, uint svidlabel,uint evid1label,uint evid2label,
    bool isExt1Restrict, bool isExt2Restrict, bool isExt1and2Restrict, bool isExt1and2SameLabel, bool isRecord1, bool isRecord2,
    bool isContinue, uint maxRowNum){
    int numBlocks;
    if(isExt1and2SameLabel && isRecord1){
        if(isExt1Restrict && isExt1and2Restrict){   INITEMB_3V(true,false,true,true,true,false)    }
        else if(isExt1Restrict){                    INITEMB_3V(true,false,false,true,true,false)   }
        else if(isExt2Restrict){                    INITEMB_3V(false,true,false,true,true,false)   }
        else if(isExt1and2Restrict){                INITEMB_3V(false,false,true,true,true,false)   }
        else{                                       INITEMB_3V(false,false,false,true,true,false)  }
    }else if(isExt1and2SameLabel && !isRecord1){
        if(isExt1Restrict && isExt1and2Restrict){   INITEMB_3V(true,false,true,true,false,false)   }
        else if(isExt1Restrict){                    INITEMB_3V(true,false,false,true,false,false)  }
        else if(isExt2Restrict){                    INITEMB_3V(false,true,false,true,false,false)  }
        else if(isExt1and2Restrict){                INITEMB_3V(false,false,true,true,false,false)  }
        else{                                       INITEMB_3V(false,false,false,true,false,false) }
    }else if(!isExt1and2SameLabel && isRecord1 && isRecord2){
        if(isExt1Restrict){                         INITEMB_3V(true,false,false,false,true,true)   }
        else if(isExt2Restrict){                    INITEMB_3V(false,true,false,false,true,true)   }
        else{                                       INITEMB_3V(false,false,false,false,true,true)  }
    }else if(!isExt1and2SameLabel && isRecord1){
        if(isExt1Restrict){                         INITEMB_3V(true,false,false,false,true,false)  }
        else if(isExt2Restrict){                    INITEMB_3V(false,true,false,false,true,false)  }
        else{                                       INITEMB_3V(false,false,false,false,true,false) }
    }else if(!isExt1and2SameLabel && isRecord2){
        if(isExt1Restrict){                         INITEMB_3V(true,false,false,false,false,true)  }
        else if(isExt2Restrict){                    INITEMB_3V(false,true,false,false,false,true)  }
        else{                                       INITEMB_3V(false,false,false,false,false,true) }
    }else if(!isExt1and2SameLabel){
        if(isExt1Restrict){                         INITEMB_3V(true,false,false,false,false,false) }
        else if(isExt2Restrict){                    INITEMB_3V(false,true,false,false,false,false) }
        else{                                       INITEMB_3V(false,false,false,false,false,false)}
    }
}


void initEmb_2V_NoHash(uint *vLabel, uint *writeRowNum, uint *edgeLabelPartition, uint *neighborsData, uint *recordPos,
    uint *newEmb, uint intervalNum, uint totSrcNum, uint svidlabel, uint evidlabel, uint isExtRestrict, uint isRecord,
    bool isContinue, uint maxRowNum){
    
    int numBlocks;
    if(isRecord && isExtRestrict) {       INITEMB_2V(true,true)  }
    else if(isRecord && !isExtRestrict){  INITEMB_2V(false,true) }
    else if(!isRecord && isExtRestrict){  INITEMB_2V(true,false) }
    else{                                 INITEMB_2V(false,false)}                             
}

void extEmb_2V_2src_1and2_sameLabel_NoHash(uint *vLabel, uint *writeRowNum, uint *edgeLabelPartition, uint *neighborsData,
    uint *recordPos, uint *auxArray, uint *newEmb, uint *partialEmb, uint intervalNum, uint partialRowNum, uint embLen,
    bool isExt1Restrict, bool isExt2Restrict, bool isExt1and2Restrict,bool isRecord, bool useRecord, bool isContinue, uint maxRowNum){

    int numBlocks;
    if (isExt1Restrict) {
        if (isExt1and2Restrict) {
            if(isRecord && useRecord){        EXTEMB_2V_2S_1L(true, false, true, true,true) }
            else if(isRecord && !useRecord) { EXTEMB_2V_2S_1L(true, false, true, true, false) }
            else if(!isRecord && useRecord) { EXTEMB_2V_2S_1L(true, false, true, false, true) }
            else {                              EXTEMB_2V_2S_1L(true, false, true, false, false) }
        } else if (isExt2Restrict) {
            if(isRecord && useRecord){        EXTEMB_2V_2S_1L(true, true, false, true, true) }
            else if(isRecord && !useRecord) { EXTEMB_2V_2S_1L(true, true, false, true, false) }
            else if(!isRecord && useRecord) { EXTEMB_2V_2S_1L(true, true, false, false, true) }
            else {                              EXTEMB_2V_2S_1L(true, true, false, false, false) }
        } else {
            if(isRecord && useRecord){        EXTEMB_2V_2S_1L(true, false, false, true, true) }
            else if(isRecord && !useRecord) { EXTEMB_2V_2S_1L(true, false, false, true, false) }
            else if(!isRecord && useRecord) { EXTEMB_2V_2S_1L(true, false, false, false, true) }
            else {                              EXTEMB_2V_2S_1L(true, false, false, false, false) }
        }
    } else {
        if (isExt1and2Restrict) {
            if(isRecord && useRecord){        EXTEMB_2V_2S_1L(false, false, true, true,true) }
            else if(isRecord && !useRecord) { EXTEMB_2V_2S_1L(false, false, true, true, false) }
            else if(!isRecord && useRecord) { EXTEMB_2V_2S_1L(false, false, true, false, true) }
            else {                              EXTEMB_2V_2S_1L(false, false, true, false, false) }
        } else if (isExt2Restrict) {
            if(isRecord && useRecord){        EXTEMB_2V_2S_1L(false, true, false, true, true) }
            else if(isRecord && !useRecord) { EXTEMB_2V_2S_1L(false, true, false, true, false) }
            else if(!isRecord && useRecord) { EXTEMB_2V_2S_1L(false, true, false, false, true) }
            else {                              EXTEMB_2V_2S_1L(false, true, false, false, false) }
        } else {
            if(isRecord && useRecord){        EXTEMB_2V_2S_1L(false, false, false, true, true) }
            else if(isRecord && !useRecord) { EXTEMB_2V_2S_1L(false, false, false, true, false) }
            else if(!isRecord && useRecord) { EXTEMB_2V_2S_1L(false, false, false, false, true) }
            else {                              EXTEMB_2V_2S_1L(false, false, false, false, false) }
        }
    }
}

void extEmb_2V_2src_1and2_notSameLabel_NoHash(uint *vLabel, uint *writeRowNum, uint *edgeLabelPartition, uint*neighborsData, 
    uint *recordPos1, uint *recordPos2, uint *auxArray, uint *newEmb, uint *partialEmb, uint intervalNum, uint partialRowNum, 
    uint embLen, bool isExt1Restrict, bool isExt2Restrict,bool isRecord1,bool isRecord2,bool useRecord1,bool useRecord2,
    bool isContinue, uint maxRowNum){

    int numBlocks;
    if(isExt1Restrict && isExt2Restrict){
        if(isRecord1 && isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_2S_2L(true,true,true,true,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_2S_2L(true,true,true,true,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_2S_2L(true,true,true,true,false,true) }
            else {                               EXTEMB_2V_2S_2L(true,true,true,true,false,false) }
        } else if(isRecord1 && !isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_2S_2L(true,true,true,false,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_2S_2L(true,true,true,false,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_2S_2L(true,true,true,false,false,true) }
            else {                               EXTEMB_2V_2S_2L(true,true,true,false,false,false) }
        }else if(!isRecord1 && isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_2S_2L(true,true,false,true,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_2S_2L(true,true,false,true,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_2S_2L(true,true,false,true,false,true) }
            else {                               EXTEMB_2V_2S_2L(true,true,false,true,false,false) }
        }else{
            if (useRecord1 && useRecord2) {      EXTEMB_2V_2S_2L(true,true,false,false,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_2S_2L(true,true,false,false,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_2S_2L(true,true,false,false,false,true) }
            else {                               EXTEMB_2V_2S_2L(true,true,false,false,false,false) }
        }
    } else if(isExt1Restrict && !isExt2Restrict){
        if(isRecord1 && isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_2S_2L(true,false,true,true,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_2S_2L(true,false,true,true,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_2S_2L(true,false,true,true,false,true) }
            else {                               EXTEMB_2V_2S_2L(true,false,true,true,false,false) }
        } else if(isRecord1 && !isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_2S_2L(true,false,true,false,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_2S_2L(true,false,true,false,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_2S_2L(true,false,true,false,false,true) }
            else {                               EXTEMB_2V_2S_2L(true,false,true,false,false,false) }
        }else if(!isRecord1 && isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_2S_2L(true,false,false,true,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_2S_2L(true,false,false,true,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_2S_2L(true,false,false,true,false,true) }
            else {                               EXTEMB_2V_2S_2L(true,false,false,true,false,false) }
        }else{
            if (useRecord1 && useRecord2) {      EXTEMB_2V_2S_2L(true,false,false,false,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_2S_2L(true,false,false,false,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_2S_2L(true,false,false,false,false,true) }
            else {                               EXTEMB_2V_2S_2L(true,false,false,false,false,false) }
        }
    }else if(!isExt1Restrict && isExt2Restrict){
        if(isRecord1 && isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_2S_2L(false,true,true,true,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_2S_2L(false,true,true,true,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_2S_2L(false,true,true,true,false,true) }
            else {                               EXTEMB_2V_2S_2L(false,true,true,true,false,false) }
        } else if(isRecord1 && !isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_2S_2L(false,true,true,false,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_2S_2L(false,true,true,false,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_2S_2L(false,true,true,false,false,true) }
            else {                               EXTEMB_2V_2S_2L(false,true,true,false,false,false) }
        }else if(!isRecord1 && isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_2S_2L(false,true,false,true,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_2S_2L(false,true,false,true,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_2S_2L(false,true,false,true,false,true) }
            else {                               EXTEMB_2V_2S_2L(false,true,false,true,false,false) }
        }else{
            if (useRecord1 && useRecord2) {      EXTEMB_2V_2S_2L(false,true,false,false,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_2S_2L(false,true,false,false,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_2S_2L(false,true,false,false,false,true) }
            else {                               EXTEMB_2V_2S_2L(false,true,false,false,false,false) }
        }
    }else{
        if(isRecord1 && isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_2S_2L(false,false,true,true,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_2S_2L(false,false,true,true,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_2S_2L(false,false,true,true,false,true) }
            else {                               EXTEMB_2V_2S_2L(false,false,true,true,false,false) }
        } else if(isRecord1 && !isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_2S_2L(false,false,true,false,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_2S_2L(false,false,true,false,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_2S_2L(false,false,true,false,false,true) }
            else {                               EXTEMB_2V_2S_2L(false,false,true,false,false,false) }
        }else if(!isRecord1 && isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_2S_2L(false,false,false,true,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_2S_2L(false,false,false,true,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_2S_2L(false,false,false,true,false,true) }
            else {                               EXTEMB_2V_2S_2L(false,false,false,true,false,false) }
        }else{
            if (useRecord1 && useRecord2) {      EXTEMB_2V_2S_2L(false,false,false,false,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_2S_2L(false,false,false,false,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_2S_2L(false,false,false,false,false,true) }
            else {                               EXTEMB_2V_2S_2L(false,false,false,false,false,false) }
        }
    }
}


void extEmb_2V_1src_1and2_sameLabel_NoHash(uint *vLabel, uint *writeRowNum, uint *edgeLabelPartition, uint*neighborsData, 
    uint *recordPos, uint *auxArray, uint *newEmb, uint *partialEmb, uint intervalNum, uint partialRowNum, uint embLen, 
    bool isExt1Restrict, bool isExt2Restrict, bool isExt1and2Restrict,bool isRecord, bool useRecord, bool isContinue, uint maxRowNum){

    int numBlocks;
    if (isExt1Restrict) {
        if (isExt1and2Restrict) {
            if(isRecord && useRecord){        EXTEMB_2V_1S_1L(true, false, true, true,true) }
            else if(isRecord && !useRecord) { EXTEMB_2V_1S_1L(true, false, true, true, false) }
            else if(!isRecord && useRecord) { EXTEMB_2V_1S_1L(true, false, true, false, true) }
            else {                            EXTEMB_2V_1S_1L(true, false, true, false, false) }
        } else if (isExt2Restrict) {
            if(isRecord && useRecord){        EXTEMB_2V_1S_1L(true, true, false, true, true) }
            else if(isRecord && !useRecord) { EXTEMB_2V_1S_1L(true, true, false, true, false) }
            else if(!isRecord && useRecord) { EXTEMB_2V_1S_1L(true, true, false, false, true) }
            else {                            EXTEMB_2V_1S_1L(true, true, false, false, false) }
        } else {
            if(isRecord && useRecord){        EXTEMB_2V_1S_1L(true, false, false, true, true) }
            else if(isRecord && !useRecord) { EXTEMB_2V_1S_1L(true, false, false, true, false) }
            else if(!isRecord && useRecord) { EXTEMB_2V_1S_1L(true, false, false, false, true) }
            else {                            EXTEMB_2V_1S_1L(true, false, false, false, false) }
        }
    } else {
        if (isExt1and2Restrict) {
            if(isRecord && useRecord){        EXTEMB_2V_1S_1L(false, false, true, true,true) }
            else if(isRecord && !useRecord) { EXTEMB_2V_1S_1L(false, false, true, true, false) }
            else if(!isRecord && useRecord) { EXTEMB_2V_1S_1L(false, false, true, false, true) }
            else {                            EXTEMB_2V_1S_1L(false, false, true, false, false) }
        } else if (isExt2Restrict) {
            if(isRecord && useRecord){        EXTEMB_2V_1S_1L(false, true, false, true, true) }
            else if(isRecord && !useRecord) { EXTEMB_2V_1S_1L(false, true, false, true, false) }
            else if(!isRecord && useRecord) { EXTEMB_2V_1S_1L(false, true, false, false, true) }
            else {                            EXTEMB_2V_1S_1L(false, true, false, false, false) }
        } else {
            if(isRecord && useRecord){        EXTEMB_2V_1S_1L(false, false, false, true, true) }
            else if(isRecord && !useRecord) { EXTEMB_2V_1S_1L(false, false, false, true, false) }
            else if(!isRecord && useRecord) { EXTEMB_2V_1S_1L(false, false, false, false, true) }
            else {                            EXTEMB_2V_1S_1L(false, false, false, false, false) }
        }
    }
}


void extEmb_2V_1src_1and2_notSameLabel_NoHash(uint *vLabel, uint *writeRowNum, uint *edgeLabelPartition, uint *neighborsData, 
    uint *recordPos1, uint *recordPos2, uint *auxArray, uint *newEmb, uint *partialEmb, uint intervalNum, uint partialRowNum, 
    uint embLen, bool isExt1Restrict, bool isExt2Restrict, bool isRecord1, bool isRecord2, bool useRecord1, bool useRecord2,
    bool isContinue, uint maxRowNum){

    int numBlocks;
    if (isExt1Restrict && isExt2Restrict) {
        if(isRecord1 && isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_1S_2L(true,true,true,true,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_1S_2L(true,true,true,true,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_1S_2L(true,true,true,true,false,true) }
            else {                               EXTEMB_2V_1S_2L(true,true,true,true,false,false) }
        } else if(isRecord1 && !isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_1S_2L(true,true,true,false,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_1S_2L(true,true,true,false,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_1S_2L(true,true,true,false,false,true) }
            else {                               EXTEMB_2V_1S_2L(true,true,true,false,false,false) }
        }else if(!isRecord1 && isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_1S_2L(true,true,false,true,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_1S_2L(true,true,false,true,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_1S_2L(true,true,false,true,false,true) }
            else {                               EXTEMB_2V_1S_2L(true,true,false,true,false,false) }
        }else{
            if (useRecord1 && useRecord2) {      EXTEMB_2V_1S_2L(true,true,false,false,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_1S_2L(true,true,false,false,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_1S_2L(true,true,false,false,false,true) }
            else {                               EXTEMB_2V_1S_2L(true,true,false,false,false,false) }
        }
    } else if(isExt1Restrict && !isExt2Restrict){
        if(isRecord1 && isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_1S_2L(true,false,true,true,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_1S_2L(true,false,true,true,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_1S_2L(true,false,true,true,false,true) }
            else {                               EXTEMB_2V_1S_2L(true,false,true,true,false,false) }
        } else if(isRecord1 && !isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_1S_2L(true,false,true,false,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_1S_2L(true,false,true,false,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_1S_2L(true,false,true,false,false,true) }
            else {                               EXTEMB_2V_1S_2L(true,false,true,false,false,false) }
        }else if(!isRecord1 && isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_1S_2L(true,false,false,true,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_1S_2L(true,false,false,true,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_1S_2L(true,false,false,true,false,true) }
            else {                               EXTEMB_2V_1S_2L(true,false,false,true,false,false) }
        }else{
            if (useRecord1 && useRecord2) {      EXTEMB_2V_1S_2L(true,false,false,false,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_1S_2L(true,false,false,false,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_1S_2L(true,false,false,false,false,true) }
            else {                               EXTEMB_2V_1S_2L(true,false,false,false,false,false) }
        }
    }else if(!isExt1Restrict && isExt2Restrict){
        if(isRecord1 && isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_1S_2L(false,true,true,true,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_1S_2L(false,true,true,true,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_1S_2L(false,true,true,true,false,true) }
            else {                               EXTEMB_2V_1S_2L(false,true,true,true,false,false) }
        } else if(isRecord1 && !isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_1S_2L(false,true,true,false,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_1S_2L(false,true,true,false,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_1S_2L(false,true,true,false,false,true) }
            else {                               EXTEMB_2V_1S_2L(false,true,true,false,false,false) }
        }else if(!isRecord1 && isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_1S_2L(false,true,false,true,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_1S_2L(false,true,false,true,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_1S_2L(false,true,false,true,false,true) }
            else {                               EXTEMB_2V_1S_2L(false,true,false,true,false,false) }
        }else{
            if (useRecord1 && useRecord2) {      EXTEMB_2V_1S_2L(false,true,false,false,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_1S_2L(false,true,false,false,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_1S_2L(false,true,false,false,false,true) }
            else {                               EXTEMB_2V_1S_2L(false,true,false,false,false,false) }
        }
    }else{
        if(isRecord1 && isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_1S_2L(false,false,true,true,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_1S_2L(false,false,true,true,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_1S_2L(false,false,true,true,false,true) }
            else {                               EXTEMB_2V_1S_2L(false,false,true,true,false,false) }
        } else if(isRecord1 && !isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_1S_2L(false,false,true,false,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_1S_2L(false,false,true,false,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_1S_2L(false,false,true,false,false,true) }
            else {                               EXTEMB_2V_1S_2L(false,false,true,false,false,false) }
        }else if(!isRecord1 && isRecord2){
            if (useRecord1 && useRecord2) {      EXTEMB_2V_1S_2L(false,false,false,true,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_1S_2L(false,false,false,true,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_1S_2L(false,false,false,true,false,true) }
            else {                               EXTEMB_2V_1S_2L(false,false,false,true,false,false) }
        }else{
            if (useRecord1 && useRecord2) {      EXTEMB_2V_1S_2L(false,false,false,false,true,true) }
            else if(useRecord1 && !useRecord2) { EXTEMB_2V_1S_2L(false,false,false,false,true,false) }
            else if(!useRecord1 && useRecord2) { EXTEMB_2V_1S_2L(false,false,false,false,false,true) }
            else {                               EXTEMB_2V_1S_2L(false,false,false,false,false,false) }
        }
    }
}


void extEmb_1V_NoHash(uint *vLabel, uint *writeRowNum, uint *edgeLabelPartition, uint *neighborsData, uint *recordPos, 
    uint *auxArray, uint *newEmb, uint *partialEmb, uint intervalNum, uint partialRowNum, uint embLen, bool isExtRestrict, 
    bool isRecord, bool useRecord, bool isLastPhase, bool isContinue, uint maxRowNum){

    int numBlocks;
    if(isLastPhase){
        if(isExtRestrict){
            if(isRecord){
                if(useRecord){ EXTEMB_1V(true,true,true,true) }
                else{          EXTEMB_1V(true,true,false,true)}
            }else{
                if(useRecord){ EXTEMB_1V(true,false,true,true) }
                else{          EXTEMB_1V(true,false,false,true)}
            }
        }else{
            if(isRecord){
                if(useRecord){ EXTEMB_1V(false,true,true,true) }
                else{          EXTEMB_1V(false,true,false,true)}
            }else{
                if(useRecord){ EXTEMB_1V(false,false,true,true)}
                else{          EXTEMB_1V(false,false,false,true)}
            }
        }
    }else{
        if(isExtRestrict){
            if(isRecord){
                if(useRecord){ EXTEMB_1V(true,true,true,false) }
                else{          EXTEMB_1V(true,true,false,false)}
            }else{
                if(useRecord){ EXTEMB_1V(true,false,true,false) }
                else{          EXTEMB_1V(true,false,false,false)}
            }
        }else{
            if(isRecord){
                if(useRecord){ EXTEMB_1V(false,true,true,false) }
                else{          EXTEMB_1V(false,true,false,false)}
            }else{
                if(useRecord){ EXTEMB_1V(false,false,true,false)}
                else{          EXTEMB_1V(false,false,false,false)}
            }
        }
    }

}
