#include <iostream>
#include <fstream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <ctime>
#include <map>
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

#define INVALID 0xffffffff
#define HASHSEED 17
typedef int LABEL;
typedef int VID;
typedef int EID;
typedef int GID;
typedef long PID;
typedef long LENGTH;
class Neighbor
{
public:
	VID vid;
	LABEL elb;
	Neighbor()
	{
		vid = -1;
		elb = -1;
	}
	Neighbor(int _vid, int _elb)
	{
		vid = _vid;
		elb = _elb;
	}
	bool operator<(const Neighbor& _nb) const
	{
		if(this->elb == _nb.elb)
		{
			return this->vid < _nb.vid;
		}
		else
		{
			return this->elb < _nb.elb;
		}
	}
};

class Element
{
public:
	int label;
	int id;
	bool operator<(const Element& _ele) const
	{
		if(this->label == _ele.label)
		{
			return this->id <_ele.id;
		}
		else
		{
			return this->label < _ele.label;
		}
	}
};

class Vertex
{
public:
	//VID id;
	LABEL label;
	//NOTICE:VID and EID is just used in this single graph
	vector<Neighbor> in;
	vector<Neighbor> out;
	Vertex()
	{
		label = -1;
	}
	Vertex(LABEL lb):label(lb)
	{
	}
};

class PCSR
{
public:
    unsigned* row_offset;  //the size is 32*key_num
    unsigned* column_index;
    unsigned key_num;  //also the group number
    unsigned edge_num;
    PCSR()
    {
        row_offset = NULL;
        column_index = NULL;
        key_num = 0;
        edge_num = 0;
    }
    ~PCSR()
    {
        delete[] row_offset;
        delete[] column_index;
    }
    inline unsigned getEdgeNum() const
    {
        return this->edge_num;
    }
};
uint32_t MurmurHash2CPU(const void * key, int len, uint32_t seed) 
{
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
class Graph
{
public:
	std::vector<Vertex> vertices;
	void addVertex(LABEL _vlb){
		this->vertices.push_back(Vertex(_vlb));
	}
	void addEdge(VID _from, VID _to, LABEL _elb){
		this->vertices[_from].out.push_back(Neighbor(_to, _elb));
		this->vertices[_to].in.push_back(Neighbor(_from, _elb));
	}

    	unsigned vertexLabelNum, edgeLabelNum;
	//CSR format: 4 pointers
	unsigned vertex_num;
	unsigned* vertex_value;

    	PCSR* csrs_in;
    	PCSR* csrs_out;

	Graph() 
	{
		vertex_num = 0;
        	csrs_in = csrs_out = NULL;
	}
	~Graph() 
	{ 
		delete[] vertex_value;
        	delete[] csrs_in;
        	delete[] csrs_out;
	}

    	void buildPCSR(PCSR* pcsr, std::vector<unsigned>& keys, int label, bool incoming){
		unsigned key_num = keys.size();
    		unsigned* row_offset = new unsigned[key_num * 32];
    		unsigned edge_num = pcsr->edge_num;
    		unsigned* column_index = new unsigned[edge_num];
    		pcsr->key_num = key_num;
    		pcsr->row_offset = row_offset;
    		pcsr->column_index = column_index;
    		for(int i = 0; i < key_num*16; ++i)
    		{
    		    row_offset[2*i] = INVALID;
    		    row_offset[2*i+1] = 0;
    		}
    		for(int i = 0; i < edge_num; ++i)
    		{
    		    column_index[i] = INVALID;
    		}
    		vector<unsigned>* buckets = new vector<unsigned>[key_num];
    		for(int i = 0; i < key_num; ++i)
    		{
    		    unsigned id = keys[i];
    		    unsigned pos = MurmurHash2CPU(&id, 4, HASHSEED) % key_num;
    		    buckets[pos].push_back(id);
    		}
    		queue<unsigned> empty_buckets;
    		for(int i = 0; i < key_num; ++i)
    		{
    		    if(buckets[i].empty())
    		    {
    		        empty_buckets.push(i);
    		    }
    		}
    		for(int i = 0; i < key_num; ++i)
    		{
    		    if(buckets[i].empty())
    		    {
    		        continue;
    		    }
    		    int tsize = buckets[i].size(), j;
    		    if(tsize > 15)
    		    {
    		        cout<<"DETECTED: more than 1 buckets are needed!"<<endl;
    		        exit(1);
    		    }
    		    else if(tsize > 30)
    		    {
    		        cout<<"DETECTED: more than 2 buckets are needed!"<<endl;
    		        exit(1);
    		    }
    		    for(j = 0; j < 15 && j < tsize; ++j)
    		    {
    		        row_offset[32*i+2*j] = buckets[i][j];
    		    }
    		    if(j < tsize)
    		    {
    		        int another_bucket = empty_buckets.front(), k = 0;
    		        empty_buckets.pop();
    		        row_offset[32*i+30] = another_bucket;
    		        while(j < tsize)
    		        {
    		            row_offset[32*another_bucket+2*k] = buckets[i][j];
    		            ++j;
    		            ++k;
    		        }
    		    }
    		}
    		delete[] buckets;

    		unsigned pos = 0;
    		for(int i = 0; i < key_num; ++i)
    		{
    		    int j;
    		    for(j = 0; j < 15; ++j)
    		    {
    		        unsigned id = row_offset[32*i+2*j];
    		        if(id == INVALID)
    		        {
    		            break;
    		        }
    		        vector<Neighbor>* adjs = &this->vertices[id].out;
    		        if(incoming)
    		        {
    		            adjs = &this->vertices[id].in;
    		        }
    		        row_offset[32*i+2*j+1] = pos;
    		        for(int k = 0; k < adjs->size(); ++k)
    		        {
    		            if((*adjs)[k].elb == label)
    		            {
    		                column_index[pos++] = (*adjs)[k].vid;
    		            }
    		        }
    		    }
    		    //set final next offset in this group, also the start offset of next valid ID
    		    row_offset[32*i+2*j+1] = pos;
    		    //row_offset[32*i+31] = pos;
    		}
	}
	void transformToCSR(){
		this->vertex_num = this->vertices.size();
		this->vertex_value = new unsigned[this->vertex_num];
		for(int i = 0; i < this->vertex_num; ++i)
		{
			this->vertex_value[i] = this->vertices[i].label;
        		sort(this->vertices[i].in.begin(), this->vertices[i].in.end());
        		sort(this->vertices[i].out.begin(), this->vertices[i].out.end());
    		}

    		//NOTICE: the edge label begins from 1
    		this->csrs_in = new PCSR[this->edgeLabelNum+1];
    		this->csrs_out = new PCSR[this->edgeLabelNum+1];
    		vector<unsigned>* keys_in = new vector<unsigned>[this->edgeLabelNum+1];
    		vector<unsigned>* keys_out = new vector<unsigned>[this->edgeLabelNum+1];
		for(int i = 0; i < this->vertex_num; ++i)
    		{
        		int insize = this->vertices[i].in.size(), outsize = this->vertices[i].out.size();
        		for(int j = 0; j < insize; ++j)
        		{
        		    	int vid = this->vertices[i].in[j].vid;
        		    	int elb = this->vertices[i].in[j].elb;
        		    	int tsize = keys_in[elb].size();
        		    	if(tsize == 0 || keys_in[elb][tsize-1] != i)
        		    	{
        		    	    keys_in[elb].push_back(i);
        		    	}
        		    	//NOTICE: we do not use C++ reference PCSR& here because it can not change(frpm p-->A to p-->B)
        		    	PCSR* tcsr = &this->csrs_in[elb];
        		    	tcsr->edge_num++;
        		}
        		for(int j = 0; j < outsize; ++j)
        		{
        		    int vid = this->vertices[i].out[j].vid;
        		    int elb = this->vertices[i].out[j].elb;
        		    int tsize = keys_out[elb].size();
			    // cout<<tsize<<endl;
        		    // for (const auto &c : keys_out[elb]) cout << c << " ";
        		    if(tsize == 0 || keys_out[elb][tsize-1] != i)
        		    {
        		        keys_out[elb].push_back(i);
        		    }
        		    PCSR* tcsr = &this->csrs_out[elb];
        		    tcsr->edge_num++;
        		}
    		}

    		for(int i = 1; i <= this->edgeLabelNum; ++i)
    		{
    		    	PCSR* tcsr = &this->csrs_in[i];
    		    	this->buildPCSR(tcsr, keys_in[i], i, true);
    		    	tcsr = &this->csrs_out[i];
    		    	this->buildPCSR(tcsr, keys_out[i], i, false);
    		}
    		delete[] keys_in;
    		delete[] keys_out;
	}

	inline unsigned vSize() const
	{
		return vertex_num;
	}
};



Graph* input(FILE* fp)
{
	char c1, c2;
	int id0, id1, id2, lb;
	bool flag = false;
	Graph* ng = NULL;

	while(true)
	{
		fscanf(fp, "%c", &c1);
		if(c1 == 't')
		{
			if(flag)
			{
				fseek(fp, -1, SEEK_CUR);
				return ng;
			}
			flag = true;
			fscanf(fp, " %c %d\n", &c2, &id0);
			if(id0 == -1)
			{
				return NULL;
			}
			else
			{
				ng = new Graph;
			}
			//read vertex num, edge num, vertex label num, edge label num
			int numVertex, numEdge;
			fscanf(fp, " %d %d %d %d\n", &numVertex, &numEdge, &(ng->vertexLabelNum), &(ng->edgeLabelNum));
		}
		else if(c1 == 'v')
		{
			fscanf(fp, " %d %d\n", &id1, &lb);
			ng->addVertex(lb); 
		}
		else if(c1 == 'e')
		{
			fscanf(fp, " %d %d %d\n", &id1, &id2, &lb);
			//NOTICE:we treat this graph as directed, each edge represents two
			//This may cause too many matchings, if to reduce, only add the first one
			//ng->addEdge(id1, id2, lb+1);
			ng->addEdge(id1, id2, lb);
			//ng->addEdge(id2, id1, lb);
		}
		else if (c1 == '#')
		{
			char tmpc[100];
           		fscanf(fp, "%[^\n]%*c", tmpc);
			continue;
		}
		else 
		{
			cerr<<"ERROR in input() -- invalid char"<<endl;
			return NULL;
		}
	}
	return NULL;
}

__device__ uint32_t 
MurmurHash2GPU(const void * key, int len, uint32_t seed) 
{
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

//searchnum should be a multiple of 32
__global__ void
join_kernel(unsigned* row_offset, unsigned *result, unsigned long keynum, unsigned long searchnum)
{
    	__shared__ unsigned s_pool1[256];
    	__shared__ unsigned s_pool3[256];
	unsigned i,laneId=threadIdx.x;
	unsigned valid=0,invalid=0;
	unsigned tmp1 = result[laneId];
	unsigned tmp2 = result[laneId+32];
    	unsigned idx = threadIdx.x & 0x1f;
    	const uint32_t m = 0x5bd1e995;
    	const int r = 24;
	for(i=1;i<=searchnum;i++){
    		unsigned bgroup = threadIdx.x & 0xffffffe0;  //equal to (x/32)*32
		unsigned vid = i;
    		unsigned bucket;// = MurmurHash2GPU(&vid, 4, HASHSEED) % keynum;

    		uint32_t h = HASHSEED ^ 4;
    		const unsigned char * data = (const unsigned char *) (&vid);
        	uint32_t k = *(uint32_t*) data;
        	k *= m;
        	k ^= k >> r;
        	k *= m;
        	h *= m;
        	h ^= k;
        	data += 4;
    		h ^= h >> 13;
    		h *= m;
    		h ^= h >> 15;
    		bucket= (h&4095);
    		s_pool1[bgroup+idx] = row_offset[32*bucket+idx];
    		if(idx == 0)
    		{
    		    s_pool3[bgroup] = INVALID;
    		}
    		if(idx < 30 && (idx&1)==0)
    		{
    		    if(s_pool1[bgroup+idx] == i)
    		    {
    		        s_pool3[bgroup] = s_pool1[bgroup+idx+1];
    		        s_pool3[bgroup+1] = s_pool1[bgroup+idx+3];
    		    }
    		}
    		if(s_pool3[bgroup] == INVALID)  // not found
    		{
    		    invalid++;
    		}
		else{
			valid++;
		}
	}
	result[laneId] = tmp1 + valid;
	result[laneId+32] = tmp2 + invalid;
}

__global__ void init(unsigned *result){
	unsigned i=blockIdx.x*blockDim.x+threadIdx.x;
	result[i] = 0;
}

void loadDataBase(string filename, string metafilename, map<uint,uint*> &labelPartitions){
    ifstream graph(filename,ios::in|ios::binary);
    if(!graph.is_open()){
        cout<<"open "<<filename <<" wrong"<<endl;
        exit(0);
    }
    ifstream metafile(metafilename,ios::in);
    if(!metafile.is_open()){
        cout<<"open metafile wrong"<<endl;
        exit(0);
    }
    map<uint,uint> v_label;
    uint oldvid,label,newvid,maxid=0;
    while(metafile>>newvid>>label>>oldvid){
        v_label.insert(pair<uint,uint>(newvid,label));
        if(newvid>maxid){
            maxid = newvid;
        }
    }
    metafile.close();
    //cout<<"maxid="<<maxid<<endl;
    //i don't know how many bytes are aligned, but 16 uints should be fine
    uint allocsizev = (((maxid+1+15)>>4)<<4);
    uint *vertex_label = (uint*)malloc(sizeof(uint)*allocsizev);
    vertex_label[0] = maxid;
	//cout<<"in loaddatabase"<<endl;
    map<uint,uint>::iterator ite;
    for(ite=v_label.begin();ite!=v_label.end();++ite){
        uint id = ite->first;
        uint lid = ite->second;
        vertex_label[id] = lid;
    }
	//cout<<"after reading v label"<<endl;
    uint tmp;
    uint j=1,index=0;
    graph.read((char*)&tmp,sizeof(uint));
    uint totpart = tmp;
	//cout<<"in loadDatabase, totpart="<<totpart<<endl;
    while(index<totpart){
        graph.read((char*)&tmp,sizeof(uint));
        uint totLabelPartEleNum = tmp;
	//cout<<"part "<<index+1<<" num="<<totLabelPartEleNum<<endl;
        uint allocsize = (((totLabelPartEleNum+15)>>4)<<4);
        uint *tmpuint = new uint[allocsize];
        graph.read((char*)&tmpuint[1],sizeof(uint)*(totLabelPartEleNum-1));
        tmpuint[0] = totLabelPartEleNum;
        labelPartitions.insert(pair<uint,uint*>(j,tmpuint));
        ++j;
        ++index;
    }
    labelPartitions.insert(pair<uint,uint*>(0xffffffff,vertex_label));
    graph.close();
}


__global__ void accessTimeMy(unsigned* row_offset, unsigned *result, unsigned intervalNum, unsigned long searchnum) {

	uint laneId, i, j;
    	laneId = threadIdx.x & 31;
	unsigned tmp1 = result[laneId];
	unsigned tmp2 = result[laneId+32];
    	__shared__ uint indexForIndex[256*2+1+1+64+1];
    	for(i=threadIdx.x;i<intervalNum*2+1+1;i=i+32){
        	indexForIndex[i] = row_offset[i];
    	}
	uint valid=0, invalid=0;
    	if(laneId==0){ indexForIndex[256*2+66] = 0; }
	row_offset = row_offset+intervalNum*2+1+1;
    	for(i=1;i<searchnum;++i) {
        	uint svid = i;
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
		uint pos = row_offset[index];
        	if(predicate==0){ 
        	  	invalid = pos+invalid;
        	}else{
			valid = pos+valid;
		}
    	}
	result[laneId] = tmp1+valid;
	result[laneId+32] = tmp2+invalid;
}

int main(int argc, const char *argv[])
{
    	cudaError_t err;
	int labelnum;
	unsigned long totnum = 0, totkeynum=0;
	int *pos;
	unsigned *keynum;
	string data = argv[1];
	string inputFileName = "/WORKSPACE/lgz/subgraph_matching/subgraph_matching/datasets/"+data+".mygraph";
	string metaFileName = "/WORKSPACE/lgz/subgraph_matching/subgraph_matching/datasets/"+data+".metadata";
	data = "/WORKSPACE/lgz/subgraph_matching/GSI/data/"+data+"/"+data+".g";
	FILE *fp = fopen(data.c_str(),"r");
	if(fp==NULL){
		cout<<"can not open file "<<data<<endl;
		return 0;
	}
	Graph *datagraph = NULL;
	datagraph = input(fp);
	if(datagraph == NULL){
		cout<<"wrong data graph"<<endl;
		return 0;
	}
	datagraph->transformToCSR();
	labelnum = datagraph->edgeLabelNum;
	cout<<"construct pcsr done, GSI has "<<labelnum<<" edge labels"<<endl;
	pos = new int[labelnum+1];
	keynum = new unsigned[labelnum+1];
	for(int label=1;label<=labelnum;++label){
		PCSR* tcsr;
        	tcsr = &(datagraph->csrs_out[label]);
		pos[label] = totnum;
		totnum = totnum+tcsr->key_num*32;
		keynum[label] = tcsr->key_num;
		totkeynum = totkeynum + tcsr->key_num;
	}

	unsigned *alldata;
    	cudaMalloc(&alldata, sizeof(unsigned)*totnum);
	for(int label=1;label<=labelnum;++label){
		PCSR* tcsr;
        	tcsr = &(datagraph->csrs_out[label]);
		totnum = tcsr->key_num*32;
		unsigned *dev = alldata+pos[label];
    		cudaMemcpy(dev, tcsr->row_offset, totnum*sizeof(unsigned), cudaMemcpyHostToDevice);
	}
	unsigned searchnum = ((datagraph->vertex_num)/64)*64;
	cout<<"searchnum="<<searchnum<<endl;
	unsigned numBlocks=1;
	unsigned *result;
	cudaMalloc(&result,sizeof(unsigned)*numBlocks*64);
	init<<<numBlocks,64>>>(result);
    	cudaDeviceSynchronize();
    
	double start = clock();
	for(int label=1;label<=labelnum;++label){
		unsigned *dev = alldata+pos[label];
		join_kernel<<<numBlocks, 32>>>(dev, result, keynum[label], searchnum);
	}
    	cudaDeviceSynchronize();
	double end = clock();
    	double endtime=(double)(end-start)/CLOCKS_PER_SEC;
    	std::cout<<"GSI time: "<<endtime*1000<<std::endl;
    	err = cudaGetLastError();
    	std::cout<<cudaGetErrorString(err)<<std::endl;
	unsigned *resulthost = (unsigned *)malloc(sizeof(unsigned)*numBlocks*64);
    	cudaMemcpy(resulthost, result, numBlocks*64*sizeof(unsigned), cudaMemcpyDeviceToHost);
	cout<<"valid = "<<resulthost[0]<< " invalid = "<< resulthost[32]<<endl;
	cudaFree(alldata);



    	map<unsigned,unsigned *> allEdgeLabelPartitions;
	loadDataBase(inputFileName, metaFileName,allEdgeLabelPartitions);
	cout<<"my has "<<allEdgeLabelPartitions.size()-1<<" edge labels"<<endl;
	unsigned *intervalNum = new unsigned[labelnum+1];
	unsigned *posmy = new unsigned[labelnum+1];
	unsigned totlen=0, len;
        for(int label=1;label<=labelnum;++label){
		posmy[label] = totlen;
		unsigned *edgeLabelPartition = allEdgeLabelPartitions[label];
        	len = edgeLabelPartition[0]-5;
		totlen = len+totlen;
        	intervalNum[label] = edgeLabelPartition[2];
	}
	unsigned *datadev;
	cudaMalloc(&datadev,sizeof(unsigned)*totlen);
	for(int label=1;label<=labelnum;++label){
		unsigned *edgeLabelPartition = allEdgeLabelPartitions[label];
        	len = edgeLabelPartition[0]-5;
        	cudaMemcpy(datadev+posmy[label],edgeLabelPartition+5,sizeof(unsigned)*len,cudaMemcpyHostToDevice);
	}
	init<<<numBlocks,64>>>(result);
    	cudaDeviceSynchronize();
	start = clock();
	for(int label=1;label<=labelnum;++label){
		unsigned *edgeLabelPartition = allEdgeLabelPartitions[label];
        	len = edgeLabelPartition[0]-5;
		accessTimeMy<<<1,32>>>(datadev+posmy[label],result, intervalNum[label],searchnum);
	}
    	cudaDeviceSynchronize();
	end = clock();
    	endtime=(double)(end-start)/CLOCKS_PER_SEC;
    	std::cout<<"my time: "<<endtime*1000<<std::endl;
    	err = cudaGetLastError();
    	std::cout<<cudaGetErrorString(err)<<std::endl;
    	cudaMemcpy(resulthost, result, numBlocks*64*sizeof(unsigned), cudaMemcpyDeviceToHost);
	cout<<"valid = "<<resulthost[0]<<" invalid="<<resulthost[32]<<endl;
	//cudaFree(datadev);
	cudaFree(result);
	free(resulthost);
}
