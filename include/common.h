//#include <sm_30_intrinsics.h>
//#include <cuda_runtime.h>
//#include <device_launch_parameters.h>
#include <map>
#include <string>

#define nullV 0xffffffff
#define WARPPERBLOCK 8
//how many vertices in each block
#define VBLOCKSIZE 128
//2^n=32*VBLOCKSIZE
#define VL2BLOCKPOW 11
//2^n=VBLOCKSIZE
#define VL1BLOCKPOW 6
#define HASHSEED 17
#define DENO (1.0/(VBLOCKSIZE+1))

#define FIND_LABEL_LIMIT(neigAddr,lowerLimit,upperLimit,label)                                          \
    found = 0;                                                                                          \
    while(upperLimit-lowerLimit>1024){                                                                  \
        uint lessThan, greatThan;                                                                       \
        uint blockSize = (upperLimit-lowerLimit)/31;                                                    \
        uint tmpPos = laneId==31?upperLimit-1:lowerLimit+blockSize*laneId;                               \
        uint tmpV = neigAddr[tmpPos];                                                                   \
        uint firstLevelLabel = vLabel[tmpV];                                                            \
        lessThan = firstLevelLabel<label;                                                               \
        greatThan = firstLevelLabel>label;                                                              \
        lessThan = __ballot_sync(0xffffffff,lessThan);                                                  \
        greatThan = __ballot_sync(0xffffffff,greatThan);                                                \
        if(lessThan==0xffffffff || greatThan==0xffffffff){ found=1; lowerLimit=upperLimit; break; }     \
        if(lessThan==0 && greatThan==0) { found = 1; break; }                                           \
        lessThan = __popc(lessThan);                                                                    \
        greatThan = __popc(greatThan);                                                                  \
        lowerLimit = lessThan==0?lowerLimit:lowerLimit+blockSize*(lessThan-1)+1;                        \
        if(greatThan>0){                                                                                \
            upperLimit = greatThan == 1 ? upperLimit-1 : lowerLimit + blockSize * (32-greatThan)+1;     \
        }                                                                                               \
    }                                                                                                   \
    if(found==0){                                                                                       \
        int len = upperLimit-lowerLimit;                                                                \
        int tmplen_32 = (len+31)&0xffffffe0;                                                            \
        for(int k=laneId;k<tmplen_32;k=k+32){                                                           \
            uint tmpV = neigAddr[lowerLimit+k];                                                         \
            uint firstLevelLabel = k<len?vLabel[tmpV]:0xffffffff;                                       \
            lessThan = firstLevelLabel==label;                                                          \
            lessThan = __ballot_sync(0xffffffff,lessThan);                                              \
            if(lessThan>0){                                                                             \
                found = 1;                                                                              \
                for(int l=tmplen_32-1-laneId;l>=(k&0xffffffe0);l=l-32){                                 \
                    uint tmpV1 = neigAddr[lowerLimit+l];                                                \
                    firstLevelLabel = l<len?vLabel[tmpV1]:0xffffffff;                                   \
                    greatThan = firstLevelLabel==label;                                                 \
                    greatThan = __ballot_sync(0xffffffff,greatThan);                                    \
                    if(greatThan>0){                                                                    \
                        greatThan = __ffs(greatThan);                                                   \
                        upperLimit = lowerLimit + (l&0xffffffe0)+32-greatThan+1;                        \
                        break;                                                                          \
                    }                                                                                   \
                }                                                                                       \
                lessThan = __ffs(lessThan)-1;                                                           \
                lowerLimit = lowerLimit + (k&0xffffffe0)+lessThan;                                      \
                break;                                                                                  \
            }                                                                                           \
        }                                                                                               \
        if(found==0) { lowerLimit = upperLimit; }                                                       \
    }


/*uint MurmurHash2(const void * key, int len, uint seed)
{
    // 'm' and 'r' are mixing constants generated offline.
    // They're not really 'magic', they just happen to work well.
    const uint m = 0x5bd1e995;
    const int r = 24;
    // Initialize the hash to a 'random' value
    uint h = seed ^ len;
    // Mix 4 bytes at a time into the hash
    const unsigned char * data = (const unsigned char *) key;
    while (len >= 4)
    {
        uint k = *((uint*) data);
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
}*/


std::map<uint,uint*> loadDataBase(std::string filename, std::string metafilename);
