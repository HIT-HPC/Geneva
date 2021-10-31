void reductionPhase(uint *vLabel, uint *totWriteRowNum, uint *edgeLabelPartition, uint *neighborsData, uint *baseRecordPos, uint *auxArray, 
    uint *partialEmb, uint intervalNum, uint partialRowNum,uint embLen,bool isLastPhase);

void copyPartialEmb(uint *src_dev, uint *dst_dev, uint numUints);