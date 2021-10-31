#include <iostream>
#include <vector>
class PhaseSingleV{
public:
    char phaseType;
    uint edgeLabel;
    std::vector<uint> pairs;
    std::vector<uint8_t> reducePos;
};
void genMatchOrderSingleV(std::string queryFileName, std::vector<PhaseSingleV>&matchOrder, std::map<uint64_t,uint> &repeatFlag,std::map<uint,uint>&restricts,std::map<uint,std::map<uint,bool>>&recordLabelStat, std::map<uint,uint> &queryLabel,std::map<uint,uint*>&allEdgeLabelPartitions,uint *baseAddr_dev, uint &stride);