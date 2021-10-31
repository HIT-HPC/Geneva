#include <fstream>
#include <map>
#include <iostream>
#include "common.h"
#include "loadDataBase.h"

using namespace std;

void extractEUINT(string &line, uint *e_info){
    int start = 0,len=0,i=0;
    int count=0;
    while(i<line.length()){
        if(line[i]>='0' && line[i]<='9'){
            start = i;
            len=0;
            while(i<line.length()){
                if(line[i]>='0' && line[i]<='9'){
                    len++;
                }else{
                    break;
                }
                i++;
            }
            string temp = line.substr(start,len);
            e_info[count] = stoul(temp);
            count++;
            if(count==3)
                break;
        }
        i++;
    }
}
void findLabelBoundry(uint *neigSAddr, uint neigNum, uint evlabel,map<uint,uint>&v_label,uint &leftBoundry,uint &rightBoundry){
    uint i;
    leftBoundry=0xffffffff;
    for(i=0;i<neigNum;++i){
        uint neigId = neigSAddr[i];
        if(v_label[neigId]==evlabel){
            leftBoundry = i;
            break;
        }
    }
    if(leftBoundry==0xffffffff){
        //cout<<"do not find leftBoundry"<<endl;
        /*for(i=0;i<neigNum;++i){
            cout<<neigSAddr[i]<<":"<<v_label[neigSAddr[i]]<<" "<<endl;
        }*/
        rightBoundry=0xffffffff;
        return;
    }
    rightBoundry = neigNum;
    for(;i<neigNum;++i){
        uint neigId = neigSAddr[i];
        if(v_label[neigId]>evlabel){
            rightBoundry = i;
            break;
        }else if(v_label[neigId]<evlabel){
            cout<<"this is the wrong label order;neigNum="<<neigNum<<endl;
        }
    }
}
uint isInNeig(uint *labelPart,uint *neigSAddr, uint neigNum, uint evid,uint evlabel,map<uint,uint>&v_label){
    uint findNum=0;
    uint leftBoundry, rightBoundry;
    findLabelBoundry(neigSAddr,neigNum,evlabel,v_label,leftBoundry,rightBoundry);
    neigNum = rightBoundry-leftBoundry;
    if(neigNum==0){
        //cout<<"do not find evid label interval"<<endl;
        return findNum;
    }
    if(neigNum<=VBLOCKSIZE){
        for(uint i=leftBoundry;i<rightBoundry;++i){
            if(evid==neigSAddr[i]){
                ++findNum;
            }
        }
    }else if(neigNum<=32*VBLOCKSIZE+32){
        uint intervalNum = (neigNum+VBLOCKSIZE) / (VBLOCKSIZE + 1); 
        uint realNeigNum = neigNum-intervalNum;
        neigNum = realNeigNum;
        uint lowerLimit = neigSAddr[leftBoundry+intervalNum]-1, upperLimit;
        for(uint i=0;i<intervalNum;++i){
            upperLimit = neigSAddr[leftBoundry+i];
            if(evid>lowerLimit && evid<=upperLimit){
                uint tmp = (i+1)*VBLOCKSIZE<=neigNum?VBLOCKSIZE:neigNum-i*VBLOCKSIZE;
                for(uint j=0;j<tmp;++j){
                    if(evid==neigSAddr[leftBoundry+intervalNum+i*VBLOCKSIZE+j]){
                        ++findNum;
                    }
                }
            }
            lowerLimit = upperLimit;
        }
    }else{
        neigNum = neigNum- 32;
        uint intervalNum = 32;
        uint intervalSize = neigNum/32;
        uint lowerLimit = neigSAddr[leftBoundry+intervalNum]-1, upperLimit;
        for(uint i=0;i<intervalNum-1;++i){
            upperLimit = neigSAddr[leftBoundry+i];
            if(evid>lowerLimit && evid<=upperLimit){
                for(uint j=0;j<intervalSize;++j){
                    if(evid==neigSAddr[leftBoundry+intervalNum+i*intervalSize+j]){
                        ++findNum;
                    }
                }
            }
            lowerLimit = upperLimit;
        }
        upperLimit = neigSAddr[leftBoundry+intervalNum-1];
        if(evid>lowerLimit && evid<=upperLimit){
            for(uint j=0;j<neigNum-31*intervalSize;++j){
                if(evid==neigSAddr[leftBoundry+intervalNum+31*intervalSize+j]){
                    ++findNum;
                }
            }
        }
    }
    return findNum;
}

uint findPos(uint *labelPart,uint *e_info, map<uint,uint> &v_label){
    uint svid = e_info[0], elabel = e_info[2], evid = e_info[1];
    //cout << "test svid="<<svid<<" edgeLabel="<<elabel<<" evid="<<evid;
    uint evlabel = v_label[evid];
    //cout<<" evid label="<<evlabel<<endl;
    uint totUINT = labelPart[0];
    uint intervalIndexAddr = 5;
    uint intervalNum = labelPart[2];
    uint intervalAddr = 5+intervalNum*2+1+1;
    uint hashIndexAddr = labelPart[3]+5;
    uint neigAddr = labelPart[4]+5;
    uint hashEntryNum = (neigAddr-hashIndexAddr)/32;
    uint *neigDataAddr = labelPart+neigAddr;
    if(intervalNum<256 && hashEntryNum>0){
        cerr << "wrong iterval num and hash num"<<endl;
        exit(0);
    }
    //cout<<"intervalNum="<<intervalNum<<endl;
    uint intervalFoundNum = 0;
    uint exactInterval = 0;
    uint preLen=0;
    for(uint i=0;i<intervalNum;++i){
        uint sbound = labelPart[intervalIndexAddr+i];
        uint totLen = labelPart[intervalIndexAddr+intervalNum+1+i];
        uint ebound = sbound+totLen-preLen-1;
        preLen = totLen;
        if(svid>=sbound && svid<=ebound){
            //cout<<"found svid interval:start boundry:"<<sbound<<" "<<"end boundry:"<<ebound<<endl;
            ++intervalFoundNum;
            uint vNeigAddr = labelPart[intervalAddr+totLen-(ebound-svid)-1];
            uint neigNum = labelPart[intervalAddr+totLen-(ebound-svid)]-vNeigAddr;
            //cout<<"neigNum="<<neigNum<<" vNeigAddr="<<vNeigAddr<<" totLen="<<totLen<<" sbound="<<sbound<<" ebound="<<ebound<<endl;
            uint flag=0;
            flag = isInNeig(labelPart,neigDataAddr+vNeigAddr,neigNum,evid,evlabel,v_label);
            //if(flag==0) { cout<< "find 0 evids in neigbors in interval" << endl; }
            if(flag==1){  ++exactInterval; }
            else if(flag>1){ cout << "find multiple evids in neighbors in interval"<<endl; exit(0);}
        }
    }
    if(intervalFoundNum==0){
        //cerr<<"do not find interval"<<endl;
        return 0;
    }
    if(intervalFoundNum>1){
        cout << "find multiple intervals contain the svid"<<endl;
        exit(0);
    }
    uint hashFoundNum = 0;
    uint exactHash = 0;
    if(hashEntryNum>0) {
        uint groupNum = 0;//MurmurHash2(&svid, 4, HASHSEED) % hashEntryNum;
        for (uint i = 0; i < 32; i = i + 2) {
            uint tmpv = labelPart[hashIndexAddr + groupNum * 32 + i];
            if (tmpv == svid) {
                ++hashFoundNum;
                uint tmpaddr = labelPart[hashIndexAddr + groupNum * 32 + i + 1];
                uint len = labelPart[hashIndexAddr + groupNum * 32 + i + 3] - tmpaddr;
                uint flag = isInNeig(labelPart, neigDataAddr + tmpaddr, len, evid,evlabel,v_label);
                if (flag > 1) { cout << " find multiple evids in neighbors in hash" << endl; exit(0); }
                else if(flag==1){ ++exactHash; }
            }
        }
    }
    return exactInterval+exactHash;
}
void checkData(map<uint,uint*>&allEdgeLabelPartitions, map<uint,uint>&reorderMap, ifstream&checkFile, map<uint,uint>&v_label){
    string line;
    map<uint,uint*>::iterator ite;
    uint forwardfound=0,backwardfound=0;
    while(getline(checkFile,line)){
        if(line.length()==0){ continue; }
        else if(line[0]=='t'){ continue; }
        else if(line[0]=='v'){ continue; }
        else if(line[0]=='e'){
            uint e_info[3];
            extractEUINT(line,e_info);
            for(uint i=0;i<2;++i){
                e_info[i] = reorderMap[e_info[i]];
            }
            uint foundNum=0, foundPart=0,index=0;
            cout<<"searching for svid:"<<e_info[0]<<" and evid:"<<e_info[1]<<" edgeLabel:"<<e_info[2]<<" evidLabel:"<<v_label[e_info[1]]<<endl;
            for(ite=allEdgeLabelPartitions.begin();ite!=allEdgeLabelPartitions.end();ite++){
                //cout<<"############label partition"<<ite->first<<"############"<<endl;
                uint *labelPart = ite->second;
                foundNum += findPos(labelPart,e_info,v_label);
                if(foundNum>0){ foundPart=index; }
                index++;
            }
            if(foundNum==0){
                cout<<"0 matches, wrong for svid:"<<e_info[0]<<" and evid:"<<e_info[1]<<"\tforward"<<endl;
                std::exit(0);
            }else if(foundNum==1) {
                ++forwardfound;
            }else if(foundNum>1){
                cout<<foundNum<<" matches, wrong for svid:"<<e_info[0]<<" and evid:"<<e_info[1]<<"\tforward"<<endl;
                std::exit(0);
            }
            cout<<endl<<endl;
            uint tmp = e_info[0];
            e_info[0] = e_info[1];
            e_info[1] = tmp;
            foundNum=0, foundPart=0,index=0;
            cout<<"searching for svid:"<<e_info[0]<<" and evid:"<<e_info[1]<<" edgeLabel:"<<e_info[2]<<" evidLabel:"<<v_label[e_info[1]]<<endl;
            for(ite=allEdgeLabelPartitions.begin();ite!=allEdgeLabelPartitions.end();ite++){
                //cout<<"############label partition"<<ite->first<<"############"<<endl;
                uint *labelPart = ite->second;
                foundNum += findPos(labelPart,e_info,v_label);
                if(foundNum>0){ foundPart=index; }
                index++;
            }
            if(foundNum==0){
                cout<<"0 matches, wrong for svid:"<<e_info[0]<<" and evid:"<<e_info[1]<<"\tbackward"<<endl;
                std::exit(0);
            }else if(foundNum==1){
                ++backwardfound;
            }else if(foundNum>1){
                cout<<foundNum<<" matches, wrong for svid:"<<e_info[0]<<" and evid:"<<e_info[1]<<"\tbackward"<<endl;
                std::exit(0);
            }

            cout<<endl<<endl;
            //if(foundNum==0){
            //    cerr<<"do not find in any label partition"<<endl;
            //}else if(foundNum==1){
            //    cerr<<"find in label partition:"<<foundPart<<endl;
            //}else{
            //    cerr<<"find in multiple label partitions"<<endl;
            //}
            //cerr<<endl;
        }
    }
    cout << "forward found: "<< forwardfound <<endl;
    cout << "backward found: "<< backwardfound << endl;
}
/*map<uint,uint*> loadDataBase(string filename, string metafilename){
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
    map<uint,uint*> labelPartitions;
    map<uint,uint> v_label;
    uint oldvid,label,newvid,maxid=0;
    while(metafile>>newvid>>label>>oldvid){
        v_label.insert(pair<uint,uint>(newvid,label));
        if(newvid>maxid){
            maxid = newvid;
        }
    }
    metafile.close();
    cout<<"maxid="<<maxid<<endl;
    uint *vertex_label = (uint*)malloc(sizeof(uint)*(maxid+1));
    vertex_label[0] = maxid;
    map<uint,uint>::iterator ite;
    for(ite=v_label.begin();ite!=v_label.end();++ite){
        uint id = ite->first;
        uint lid = ite->second;
        vertex_label[id+1] = lid;
    }
    char tmp[4]={'1','1','1','1'};
    uint i=0,j=0,index=0;
    graph.read(tmp,4);
    cout<<"uint size:"<<sizeof(uint)<<endl;
    uint *tmptotpart = (uint*)tmp;
    uint totpart = (*tmptotpart);
    while(index<totpart){
        graph.read(tmp,4);
        uint *tmpuint = (uint*)tmp;
        uint totLabelPartEleNum = (*tmpuint);
        //std::cout<<totLabelPartEleNum<<" uints"<<std::endl;
        tmpuint = (uint*)malloc(sizeof(uint)*totLabelPartEleNum);
        uint *tmpa = tmpuint+1;
        graph.read((char*)tmpa,sizeof(uint)*(totLabelPartEleNum-1));
        tmpuint[0] = totLabelPartEleNum;
        labelPartitions.insert(pair<uint,uint*>(j,tmpuint));
        ++j;
        ++index;
    }
    labelPartitions.insert(pair<uint,uint*>(0xffffffff,vertex_label));
    graph.close();
    return labelPartitions;
}*/
map<uint,uint*> loadDataBase(string filename, string metafilename){
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
    map<uint,uint*> labelPartitions;
    map<uint,uint> v_label;
    uint oldvid,label,newvid,maxid=0;
    while(metafile>>newvid>>label>>oldvid){
        v_label.insert(pair<uint,uint>(newvid,label));
        if(newvid>maxid){
            maxid = newvid;
        }
    }
    metafile.close();
    cout<<"maxid="<<maxid<<endl;
    //i don't know how many bytes are aligned, but 16 uints should be fine
    uint allocsizev = (((maxid+1+15)>>4)<<4);
    uint *vertex_label = (uint*)malloc(sizeof(uint)*allocsizev);
    vertex_label[0] = maxid;
    map<uint,uint>::iterator ite;
    for(ite=v_label.begin();ite!=v_label.end();++ite){
        uint id = ite->first;
        uint lid = ite->second;
        vertex_label[id] = lid;
    }
    uint tmp;
    uint i=0,j=1,index=0;
    graph.read((char*)&tmp,sizeof(uint));
    uint totpart = tmp;
    while(index<totpart){
        graph.read((char*)&tmp,sizeof(uint));
        uint totLabelPartEleNum = tmp;
        //uint totLabelPartEleNum = 10693;
        //std::cout<<totLabelPartEleNum<<" uints"<<std::endl;
        uint allocsize = (((totLabelPartEleNum+15)>>4)<<4);
        //std::cout<<"alloc size:"<<allocsize<<std::endl;
        uint *tmpuint = new uint[allocsize];
        //std::cout<<"alloc done:"<<allocsize<<std::endl;
        graph.read((char*)&tmpuint[1],sizeof(uint)*(totLabelPartEleNum-1));
        //uint readbyes = graph.gcount();
        //std::cout<<"read: "<<(readbyes*1.0)/sizeof(uint)<<" uints"<<std::endl;
        tmpuint[0] = totLabelPartEleNum;
        labelPartitions.insert(pair<uint,uint*>(j,tmpuint));
        ++j;
        ++index;
    }
    labelPartitions.insert(pair<uint,uint*>(0xffffffff,vertex_label));
    graph.close();
    return labelPartitions;
}

void loadDataBaseAndCheck(string filename, string metafilename, string checkFileName){
    ifstream graph(filename,ios::in|ios::binary);
    if(!graph.is_open()){
        cout<<"open "<<filename <<" wrong"<<endl;
        return;
    }else{
        cout<<"open "<<filename <<" success"<<endl;
    }
    ifstream checkFile(checkFileName,ios::in);
    if(!checkFile.is_open()){
        cout<<"open checkfile wrong"<<endl;
        return;
    }else{
        cout<<"open checkfile:"+checkFileName+" success"<<endl;
    }
    ifstream metafile(metafilename,ios::in);
    if(!metafile.is_open()){
        cout<<"open metafile wrong"<<endl;
        return;
    }else{
        cout<<"open metafile:"+metafilename+" success"<<endl;
    }
    //ofstream rawdata("raw.txt",ios::out);
    map<uint,uint*> labelPartitions;
    map<uint,uint> v_label;
    map<uint,uint> reorderMap;
    uint oldvid,label,newvid;
    while(metafile>>newvid>>label>>oldvid){
        v_label.insert(pair<uint,uint>(newvid,label));
        reorderMap.insert(pair<uint,uint>(oldvid,newvid));
    }
    metafile.close();
    char tmp[4]={'1','1','1','1'};
    uint i=0,j=0,index=0;
    graph.read(tmp,4);
    uint *tmptotpart = (uint*)tmp;
    //rawdata<< *tmptotpart<<endl;;
    uint totpart = (*tmptotpart);
    cout<<"tot edge labes = "<< totpart<<endl;
    while(index<totpart){
        graph.read(tmp,4);
        uint *tmpuint = (uint*)tmp;
        uint totLabelPartEleNum = (*tmpuint);
        tmpuint = (uint*)malloc(sizeof(uint)*totLabelPartEleNum);
        uint *tmpa = tmpuint+1;
        graph.read((char*)tmpa,sizeof(uint)*(totLabelPartEleNum-1));
        tmpuint[0] = totLabelPartEleNum;
        labelPartitions.insert(pair<uint,uint*>(j,tmpuint));
        /*for(uint i=0;i<totLabelPartEleNum;++i){
            rawdata << " "<< tmpuint[i];
        }
        rawdata << endl<<endl;*/
        ++j;
        ++index;
        cout<<tmpuint[0]<<" "<<tmpuint[1]<<" "<<tmpuint[2]<<" "<<tmpuint[3]<<" "<<tmpuint[4]<<endl;
    }
    checkData(labelPartitions,reorderMap,checkFile,v_label);
    checkFile.close();
    graph.close();
    map<uint,uint*>::iterator tmpiter;
    for(tmpiter=labelPartitions.begin();tmpiter!=labelPartitions.end();tmpiter++){
        uint *tmpuint = tmpiter->second;
        free(tmpuint);
    }
}
//int main(){
//    string path = "../";
//    string baseFileName = "soc-Flickr-ASU";
//    string inputFileName = path+baseFileName+".g.graph";
//    string checkFileName = path+baseFileName+".g";
//    string metaFileName = path+"metadata.g";
//    //loadDataBase(path+"BA-2_24_60-L2.g.graph", path+"metadata.g");
//    loadDataBase(inputFileName, checkFileName, metaFileName);
//    return 0;
//}
