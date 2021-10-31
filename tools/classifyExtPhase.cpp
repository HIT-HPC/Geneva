#include <iostream>
#include <fstream>
#include <algorithm>
#include <hash_map>
#include <ctime>
#include "common.h"
#include "calDVPhaseCount.h"

using namespace __gnu_cxx;

std::map<uint64_t,uint> repeatFlag;
std::vector<Phase> matchOrder;
std::map<uint,std::map<uint,bool>> edgeLabelStat;
std::map<uint,uint> restricts;
std::map<uint,uint*> allEdgeLabelPartitions;
std::map<uint,uint> queryLabel;
uint *vLabel;

int main(int argc, char *argv[]){
    std::string path = "../datasets/";
    std::string baseFileName = std::string(argv[1]);
    std::string inputFileName = path+baseFileName+".mygraph";
    std::string checkFileName = path+baseFileName+".g";
    std::string metaFileName = path+baseFileName+".metadata";
    std::string queryFileName = path+std::string(argv[2])+".query";

    genMatchOrder(queryFileName,matchOrder,repeatFlag,restricts,edgeLabelStat,queryLabel);
    Phase &tmpPhase = matchOrder[0];
	int DVPhasecount=0;
    for(uint i=0;i<matchOrder.size();++i){
        Phase &phase = matchOrder[i];
        std::vector<uint> &pairs = phase.pairs;
	if(phase.phaseType=='E' && pairs.size()==4){
		std::cout<<"svid:"<<pairs[0]<<" evid:"<<pairs[1]<<" svid:"<<pairs[2]<<" evid:"<<pairs[3]<<std::endl;
		DVPhasecount++;
	}
    }
	std::cout<<"DV phase count: "<<DVPhasecount<<std::endl;

    return 0;
}
