void initEmb_3V_1and2_NoHash(uint *vLabel_dev, uint *totWriteRowNum_dev, uint *edgeLabelPartition_dev, uint *neighborsData_dev,   
    uint *recordPos1_dev, uint *recordPos2_dev, uint *newEmb_dev, uint intervalNum, uint totSrcNum,  uint svidlabel, uint evid1label,
    uint evid2label,bool isExt1Restrict, bool isExt2Restrict,  bool isExt1and2Restrict, bool isExt1and2SameLabel, bool isRecord1, 
    bool isRecord2, bool isContinue, uint maxRowNum);

void initEmb_2V_NoHash(uint *vLabel_dev, uint *totWriteRowNum_dev, uint *edgeLabelPartition_dev, uint *neighborsData_dev, uint *recordPos_dev,
    uint *newEmb_dev, uint intervalNum, uint totSrcNum, uint svidlabel, uint evidlabel, uint isRestrict, uint isRecord, bool isContinue, uint maxRowNum);

void extEmb_2V_2src_1and2_sameLabel_NoHash(uint *vLabel_dev, uint *totWriteRowNum_dev, uint *edgeLabelPartition_dev,
    uint *neighborsData_dev, uint *recordPos_dev, uint *auxArray_dev, uint *newEmb_dev, uint *partialEmb_dev, uint intervalNum, 
    uint partialRowNum, uint embLen, bool isExt1Restrict, bool isExt2Restrict, bool isExt1and2Restrict,bool isRecord, bool useRecord, 
    bool isContinue, uint maxRowNum);

void extEmb_2V_2src_1and2_notSameLabel_NoHash(uint *vLabel_dev, uint *totWriteRowNum_dev, uint *edgeLabelPartition_dev, 
    uint *neighborsData_dev, uint *recordPos1_dev, uint *recordPos2_dev, uint *auxArray_dev, uint *newEmb_dev, uint *partialEmb_dev, 
    uint intervalNum, uint partialRowNum, uint embLen, bool isExt1Restrict, bool isExt2Restrict, bool isRecord1, bool isRecord2, 
    bool useRecord1, bool useRecord2, bool isContinue, uint maxRowNum);

void extEmb_2V_1src_1and2_sameLabel_NoHash(uint *vLabel_dev, uint *totWriteRowNum_dev, uint *edgeLabelPartition_dev, 
    uint *neighborsData_dev, uint *recordPos_dev, uint *auxArray_dev, uint *newEmb_dev, uint *partialEmb_dev, uint intervalNum, 
    uint partialRowNum, uint embLen, bool isExt1Restrict, bool isExt2Restrict, bool isExt1and2Restrict,bool isRecord, bool useRecord,
    bool isContinue, uint maxRowNum);

void extEmb_2V_1src_1and2_notSameLabel_NoHash(uint *vLabel_dev, uint *totWriteRowNum_dev, uint *edgeLabelPartition_dev, 
    uint*neighborsData_dev, uint *recordPos1_dev, uint *recordPos2_dev, uint *auxArray_dev, uint *newEmb_dev, uint *partialEmb_dev, 
    uint intervalNum, uint partialRowNum, uint embLen, bool isExt1Restrict,bool isExt2Restrict, bool isRecord1, bool isRecord2, 
    bool useRecord1, bool useRecord2, bool isContinue, uint maxRowNum);

void extEmb_1V_NoHash(uint *vLabel_dev, uint *totWriteRowNum_dev, uint *edgeLabelPartition_dev, uint*neighborsData_dev,
    uint *recordPos_dev, uint *auxArray_dev, uint *newEmb_dev, uint *partialEmb_dev, uint intervalNum, uint partialRowNum, 
    uint embLen, bool isExtRestrict, bool isRecord, bool useRecord, bool isLastPhase, bool isContinue, uint maxRowNum);