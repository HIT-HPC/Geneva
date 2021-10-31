#include <cuda_runtime.h>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <thread>
#include <list>
#include <algorithm>
#include "common.h"
#include "genMatchOrder.h"

using namespace std;
class queryNeigInfo{
public:
    queryNeigInfo(){
        neigs = 0;
    }
    //records all neighbors of a vertex, if the vid of a neighbor is n, then set
    //the nth bit of neigs to 1
    uint neigs;
    //each vertex has many edges with different labels
    //each distinct label is an entry in the map, and the neighbors with this label
    //are stored in the vector[1] as a bit, vector[0] is the total number
    std::map<uint,std::vector<uint> > edgeLabelNeigs;
};
class queryEdge{
public:
    uint svid;
    uint evid;
    uint edgeLabel;
};

class triInfo{
public:
    uint edgeIndex;
    uint triCount;
    uint triScore;
    uint maxScoreVID;
};
class NeigInfo{
public:
    uint commonVNum;
    uint svid;
    uint evid;
    uint edgeLabel;
    uint evLabel;
    uint evNum;
};
class phaseKeyCount{
public:
    phaseKeyCount(){ key = 0; count = 0; }
    uint64_t key;
    uint count;
    std::vector<uint64_t> keyList;
    std::vector<uint> phaseList;//phase ids that contain this label
};
//3 means 3 vertices, ext1 is symmetric with ext2, 2 means 2 vertices, src1 is symmetric with ext1
//1 means no symmetries
uint partialQueryVidMap;
list<uint> partialQueryVid;
vector<int> vidToIndex;
vector<uint> vertexNeig;
vector<uint> adjMatrix;
uint vnum;
bool sortDescend(triInfo& a, triInfo& b){
    return a.triCount>b.triCount;
}

void addToQueryGraph(uint svid, uint evid, uint edgelabel,vector<queryNeigInfo>&query){
    //set the evid bit of neigs to 1
    query[svid].neigs = (query[svid].neigs | (1<<evid));
    map<uint, vector<uint>> &neig = query[svid].edgeLabelNeigs;
    map<uint, vector<uint>>::iterator ite = neig.find(edgelabel);
    if (ite != neig.end()) {
        vector<uint> &elneig = ite->second;
        elneig[0] = elneig[0]+1;
        elneig[1] = (elneig[1] | (1<<evid));
    }else{
        vector<uint> tmpv;
        tmpv.push_back(1);
        uint tmp = (1<<evid);
        tmpv.push_back(tmp);
        neig.insert(pair<uint,vector<uint> >(edgelabel,tmpv));
    }
}

//31 30 29 ... 0
void countTriangle(uint edgeIndex,vector<triInfo> &countResult,vector<queryEdge>&edges, map<uint,uint>&queryLabel){
    uint vnum = vertexNeig.size();
    queryEdge &edge = edges[edgeIndex];
    uint svidNeig = vertexNeig[edge.svid];
    uint evidNeig = vertexNeig[edge.evid];
    uint triangles = svidNeig & evidNeig;
    countResult[edgeIndex].edgeIndex = edgeIndex;
    countResult[edgeIndex].triCount = __builtin_popcount(triangles);
    if(triangles==0){
        return;
    }
    int prei=-1,prescore=-1;
    for(uint i=0;i<32;++i) {
        uint tmp = (1 << i) & triangles;
        if (tmp == 0) { continue; }
        uint e1l = adjMatrix[edge.svid * vnum + i], e2l = adjMatrix[edge.evid * vnum + i], e3l = edge.edgeLabel;
        int edgescore = (e1l == e2l) + (e1l == e3l) + (e2l == e3l);
        edgescore = edgescore * 10;
        uint v1l = queryLabel[edge.svid], v2l = queryLabel[edge.evid], v3l = queryLabel[i];
        int vertexscore = (v1l == v2l) + (v1l == v3l) + (v2l == v3l);
        if (vertexscore + edgescore > prescore) {
            prescore = vertexscore + edgescore;
            prei = i;
        }
    }
    countResult[edgeIndex].triScore = prescore;
    countResult[edgeIndex].maxScoreVID = prei;
}

bool sortByEdgeLabel(NeigInfo &a, NeigInfo &b){
    if(a.edgeLabel<b.edgeLabel){
        return true;
    }else if(a.edgeLabel==b.edgeLabel){
        return a.evLabel<b.evLabel;
    }
    return false;
}
bool multiLevelSort(NeigInfo &a, NeigInfo &b){
    //1st priority, more back edges
    if(a.commonVNum > b.commonVNum){
        return true;
    }else if(a.commonVNum==b.commonVNum){
        //2nd priority, one vertex has many neighbors with same edge labels
        if(a.evNum>b.evNum){
            return true;
        }else if(a.evNum==b.evNum){
            //3rd priority, choose vertices has the same edge label, though they are not the same
            if(a.edgeLabel<b.edgeLabel){
                return true;
            }else if(a.edgeLabel==b.edgeLabel){
                if(a.svid<b.svid){
                    return true;
                }else if(a.svid==b.svid){
                    if(a.evLabel<b.evLabel){
                        return true;
                    }
                }
            }
        }
    }
    return false;
}
class initPartialQueryInfo{
public:
    uint maxSameEdgeLabelNum;
    uint maxSameEdgeEvidLabelNum;
    uint maxNeigofEvidNum;
    uint firstEndVid;
    uint secondEndVid;
    uint svid;
};

//the logic of this func
/**
 * for each svid and evid pair, we build a information and use this information to compare all pairs and
 * choose the best two or one pair with the same svid. and we use a global variable globalMaxSameSvidEdgeLabelNum
 * to indicate if there are two evids of one svid that have the same edge label.
 * 1. each pair has a variable sameSvidAndEdgeLabelNum to record how many pairs that have the same svid
 *      and edge label with this pair. if sameSvidEdgeLabelNum>2, we just set it to 2 because
 *      we can at most extend 2 vertices from one svid.
 * 2. each pair has a variable sameSvidAndEdgeAndEvidLabelNum to record how many pairs that have the same svid
 *      , the same edge label and the same evid label with this pair. if this value>2, then set it to 2.
 * 3. each pair has a variable evidNeigNum to record the number of neighbors of evid.
 *
 * if  globalMaxSameSvidEdgeLabelNum ==2:
 *      we first sort pairs by sameSvidAndEdgeLabelNum and then sameSvidAndEdgeAndEvidLabelNum and then svid
 *      and then edge label and last evidNeigNum.
 **/
void genInitialPartialQueryNoTriangle(vector<queryNeigInfo> &query, map<uint,uint> &queryLabel,
    uint &svid, uint &evid1, uint &evid2){

    uint globalSvid=0xffffffff, globalFirstMax=0, globalSecondMax=0;
    uint globalFirstEndVid=0xffffffff,globalSecondEndVid=0xffffffff;
    map<uint64_t,vector<uint>> sameEdgeEvidLabel;
    map<uint64_t,vector<uint>> sameEdgeLabel;
    map<uint64_t,vector<uint>> sameSvid;
    for(uint i=0;i<vnum;++i){
        uint64_t keysvid = i;
        queryNeigInfo &vNeigInfo = query[i];
        map<uint, vector<uint> > &edgeLabelNeigs = vNeigInfo.edgeLabelNeigs;
        for(auto ite=edgeLabelNeigs.begin();ite!=edgeLabelNeigs.end();++ite){
            uint edgelabel = ite->first;
            uint64_t keysvidedge = edgelabel;
            keysvidedge = (keysvidedge<<20) | keysvid;
            uint sameedgevnum = ite->second[0];
            uint neigs = ite->second[1];
            uint totNeigs = __builtin_popcount(neigs);
            if(totNeigs>=2){
                map<uint,vector<uint>> tmpmap;
                for(uint k=0;k<32;++k){
                    uint evid = ((neigs&(1<<k))>>k);
                    if(evid==0) { continue; }
                    evid = k;
                    auto itetmp = tmpmap.find(queryLabel[evid]);
                    if(itetmp==tmpmap.end()){
                        vector<uint> tmpv;
                        tmpv.push_back(evid);
                        tmpmap.insert(pair<uint,vector<uint>>(queryLabel[evid],tmpv));
                    }else{
                        vector<uint> &tmpv = tmpmap[queryLabel[evid]];
                        tmpv.push_back(evid);
                    }
                }
                for(auto itea=tmpmap.begin();itea!=tmpmap.end();++itea){
                    uint64_t keysvidedgeevid = itea->first;
                    keysvidedgeevid = (keysvidedgeevid<<40) | keysvidedge;
                    vector<uint> &tmpv = itea->second;
                    if(tmpv.size()>=2){
                        sameEdgeEvidLabel.insert(pair<uint64_t,vector<uint>>(keysvidedgeevid,tmpv));
                    }else{
                        auto itesvidedge = sameEdgeLabel.find(keysvidedge);
                        if(itesvidedge==sameEdgeLabel.end()){
                            sameEdgeLabel.insert(pair<uint64_t,vector<uint>>(keysvidedge,tmpv));
                        }else{
                            vector<uint> &tmpv1 = sameEdgeLabel[keysvidedge];
                            tmpv1.insert(tmpv1.end(),tmpv.begin(),tmpv.end());
                        }
                    }
                }
            }else{
                for(uint k=0;k<32;++k){
                    uint evid = ((neigs&(1<<k))>>k);
                    if(evid==0) { continue; }
                    evid = k;
                    auto itesvid = sameSvid.find(keysvid);
                    if(itesvid==sameSvid.end()){
                        vector<uint> tmpv;
                        tmpv.push_back(evid);
                        sameSvid.insert(pair<uint64_t,vector<uint>>(keysvid,tmpv));
                    }else{
                        vector<uint> &tmpv = sameSvid[keysvid];
                        tmpv.push_back(evid);
                    }
                }
            }
        }
    }
    map<uint64_t,vector<uint>> tmpmap;
    if(sameEdgeEvidLabel.size()>0)  { tmpmap = sameEdgeEvidLabel; }
    else if(sameEdgeLabel.size()>0) { tmpmap = sameEdgeLabel; }
    else{                             tmpmap = sameSvid; }
    for(auto ite=tmpmap.begin();ite!=tmpmap.end();++ite){
        uint64_t key = ite->first;
        vector<uint> &tmpv = ite->second;
        uint firstMaximum=0,secondMaximum=0,firstEndVid=0xffffffff,secondEndVid=0xffffffff;
        for(uint i=0;i<tmpv.size();++i){
            uint evid = tmpv[i];
            uint totNeigs = __builtin_popcount(query[evid].neigs);
            if(totNeigs>firstMaximum){
                secondMaximum = firstMaximum;
                secondEndVid = firstEndVid;
                firstMaximum = totNeigs;
                firstEndVid = evid;
            }else if(totNeigs>secondMaximum){
                secondMaximum = totNeigs;
                secondEndVid = evid;
            }
        }
        if(sameEdgeEvidLabel.size()>0 || sameEdgeLabel.size()>0){
            if(firstMaximum+secondMaximum>globalFirstMax+globalSecondMax){
                globalFirstMax = firstMaximum;
                globalFirstEndVid = firstEndVid;
                globalSecondMax = secondMaximum;
                globalSecondEndVid = secondEndVid;
                globalSvid = (uint)(key&0x00000000000fffff);
            }
        }else{
            if(firstMaximum>globalFirstMax){
                globalFirstMax = firstMaximum;
                globalFirstEndVid = firstEndVid;
                globalSvid = (uint)(key&0x00000000000fffff);
            }
        }
    }
    svid = globalSvid;
    evid1 = globalFirstEndVid;
    evid2 = globalSecondEndVid;
}

void addExtPhase(uint svid1,uint extvid1,uint svid2,uint extvid2, uint sameEdgeLabelNum,
                 map<uint64_t,uint> &extendPhaseKeyCount, vector<Phase>&matchOrder,map<uint,uint>&queryLabel){
    if(sameEdgeLabelNum==2 && queryLabel[extvid2]<queryLabel[extvid1]){
        uint tmp = svid1;
        uint tmp1 = extvid1;
        svid1 = svid2;
        extvid1=extvid2;
        svid2 = tmp;
        extvid2 = tmp1; 
    }
    Phase extphase1;
    extphase1.phaseType = 'E';
    extphase1.pairs.push_back(svid1);
    extphase1.pairs.push_back(extvid1);
    extphase1.edgeLabel = adjMatrix[svid1*vnum+extvid1];
    uint64_t key = extphase1.edgeLabel;
    key = (key<<32) | queryLabel[extvid1];
    map<uint64_t,uint>::iterator ite;
    if((ite=extendPhaseKeyCount.find(key))==extendPhaseKeyCount.end()){
        extendPhaseKeyCount.insert(pair<uint64_t,uint>(key,1));
    }else{
        ite->second = ite->second+1;
    }
    partialQueryVidMap = partialQueryVidMap | (1<<extvid1);
    vidToIndex[extvid1] = partialQueryVid.size();
    partialQueryVid.push_back(extvid1);
    vertexNeig[svid1] = (~(1<<extvid1))&vertexNeig[svid1];
    vertexNeig[extvid1] = (~(1<<svid1))&vertexNeig[extvid1];
    if(sameEdgeLabelNum==2){
        extphase1.pairs.push_back(svid2);
        extphase1.pairs.push_back(extvid2);
        key = adjMatrix[svid2*vnum+extvid2];
        key = (key<<32) | queryLabel[extvid2];
        if((ite=extendPhaseKeyCount.find(key))==extendPhaseKeyCount.end()){
            extendPhaseKeyCount.insert(pair<uint64_t,uint>(key,1));
        }else{
            ite->second = ite->second+1;
            //ite->second.phaseList.push_back(matchOrder.size());
        }
        partialQueryVidMap = partialQueryVidMap | (1<<extvid2);
        vidToIndex[extvid2] = partialQueryVid.size();
        partialQueryVid.push_back(extvid2);
        vertexNeig[svid2] = (~(1<<extvid2))&vertexNeig[svid2];
        vertexNeig[extvid2] = (~(1<<svid2))&vertexNeig[extvid2];
    }
    matchOrder.push_back(extphase1);
}
void addReducePhase(vector<uint> &pairs,map<uint64_t,phaseKeyCount> &reducePhaseKeyCount,
                    map<uint64_t,uint> &othersKeyCount,vector<Phase>&matchOrder,map<uint,uint>&queryLabel){
    Phase reducephase;
    reducephase.phaseType = 'R';
    uint8_t totEndVLabel = 1;
    reducephase.reducePos.push_back(totEndVLabel);
    for(uint i=0;i<pairs.size()/2;++i){
        uint svid1 = pairs[i*2], extvid1 = pairs[i*2+1];
        reducephase.pairs.push_back(svid1);
        reducephase.pairs.push_back(extvid1);
        vertexNeig[svid1] = (~(1<<extvid1))&vertexNeig[svid1];
        vertexNeig[extvid1] = (~(1<<svid1))&vertexNeig[extvid1];
        reducephase.edgeLabel = adjMatrix[svid1*vnum+extvid1];
        uint64_t key = reducephase.edgeLabel;
        key = (key<<32) | queryLabel[extvid1];
        uint64_t key1 = adjMatrix[svid1*vnum+extvid1];
        key1 = (key1<<32) | queryLabel[svid1];
        if(key==key1){
            map<uint64_t,uint>::iterator ite = othersKeyCount.find(key);
            if(ite==othersKeyCount.end()){
                othersKeyCount.insert(pair<uint64_t,uint>(key,1));
            }else{
                ite->second = ite->second+1;
            }
        }else{
            map<uint64_t,phaseKeyCount>::iterator ite = reducePhaseKeyCount.find(key);
            if(ite==reducePhaseKeyCount.end()){
                phaseKeyCount tmp;
                tmp.key = key;
                tmp.count = 1;
                tmp.keyList.push_back(key1);
                tmp.phaseList.push_back(matchOrder.size());
                tmp.phaseList.push_back(i);
                reducePhaseKeyCount.insert(pair<uint64_t,phaseKeyCount>(key,tmp));
            }else{
                ite->second.count = ite->second.count+1;
                ite->second.keyList.push_back(key1);
                ite->second.phaseList.push_back(matchOrder.size());
                ite->second.phaseList.push_back(i);
            }
            ite = reducePhaseKeyCount.find(key1);
            if(ite==reducePhaseKeyCount.end()){
                phaseKeyCount tmp;
                tmp.key = key1;
                tmp.count = 1;
                tmp.keyList.push_back(key);
                tmp.phaseList.push_back(matchOrder.size());
                tmp.phaseList.push_back(i);
                reducePhaseKeyCount.insert(pair<uint64_t,phaseKeyCount>(key1,tmp));
            }else{
                ite->second.count = ite->second.count+1;
                ite->second.keyList.push_back(key);
                ite->second.phaseList.push_back(matchOrder.size());
                ite->second.phaseList.push_back(i);
            }
        }
    }
    reducephase.reducePos.push_back(pairs.size()/2);
    reducephase.reducePos[0] = totEndVLabel;
    matchOrder.push_back(reducephase);
}

bool isEquivalent(uint vid1, uint vid2, vector<queryNeigInfo> &query,map<uint,uint>&queryLabel){
    if(vid1==vid2){
        return true;
    }
    if(queryLabel[vid1]!=queryLabel[vid2]){
        return false;
    }
    map<uint,vector<uint> >& vid1Neig = query[vid1].edgeLabelNeigs;
    map<uint,vector<uint> >& vid2Neig = query[vid2].edgeLabelNeigs;
    if(vid1Neig.size()!=vid2Neig.size()){
        return false;
    }
    for(auto vid1Iter=vid1Neig.begin(),vid2Iter=vid2Neig.begin();vid1Iter!=vid1Neig.end();vid1Iter++,vid2Iter++){
        if(vid1Iter->first!=vid2Iter->first || vid1Iter->second[0]!=vid2Iter->second[0]){
            return false;
        }
    }
    return true;
}

//newId has not been added to autCandidate
bool isPossibleAut(vector<uint> &autCandidate, uint newID, vector<queryNeigInfo> &query, map<uint,uint> &queryLabel){
    //new id is the last vid of autCandidate, its position represents the original id
    uint oriIDofNewID = autCandidate.size();
    if(!isEquivalent(oriIDofNewID,newID,query,queryLabel)){
        return false;
    }
    for(int i=0;i<autCandidate.size();++i){
        uint oriEdgeLabel = adjMatrix[i*vnum+oriIDofNewID];
        uint newEdgeLabel = adjMatrix[autCandidate[i]*vnum+newID];
        if(oriEdgeLabel!=newEdgeLabel){
            return false;
        }
    }
    return true;
}

void genAutomorphism(vector<uint> autCandidate, bool isIdentity,vector<queryNeigInfo> &query,
                     vector<vector<uint>>&autResults,map<uint,uint> &queryLabel){
    if(autCandidate.size()==vnum-1){
        for(int i=0;i<vnum;++i){
            bool isExist = false;
            for(int j=0;j<autCandidate.size();++j){ if(autCandidate[j]==i){ isExist = true; break; } }
            if(!isExist){
                if(isIdentity && i==autCandidate.size()){ return; }
                if(isPossibleAut(autCandidate, i, query,queryLabel)){
                    autCandidate.push_back(i);
                    autResults.push_back(autCandidate);
                }
                break;
            }
        }
    }else{
        for(int i=0;i<vnum;++i){
            bool isExist = false;
            for(int j=0;j<autCandidate.size();++j){ if(autCandidate[j]==i){ isExist = true; break; } }
            if(!isExist){
                if(isPossibleAut(autCandidate, i, query,queryLabel)){
                    autCandidate.push_back(i);
                    bool isIdentityTmp = isIdentity&&(i==autCandidate.size()-1);
                    genAutomorphism(autCandidate, isIdentityTmp, query, autResults,queryLabel);
                    autCandidate.pop_back();
                }
            }
        }
    }
}

//the goal of restrictions is to delete all automorphisms until there is no automorphisms
//therefore, once a restriction can eliminate an automorphism, we delete this automorphism
//and do not need to justify other vertices in this automorphism
void genRestriction(vector<vector<uint>> automorphisms, map<uint,uint>& restricts,vector<Phase>&matchOrder){
    vector<bool> isVisited(vnum,false);
    for(int i=0;i<matchOrder.size();++i){
        Phase &tmpPhase = matchOrder[i];
        if(tmpPhase.phaseType=='R') continue;
        vector<uint> &vidPairs = tmpPhase.pairs;
        for(int j=0;j<vidPairs.size();++j){
            uint vid = vidPairs[j];
            if(isVisited[vid]) continue;
            isVisited[vid] = true;
            vector<vector<uint>> stabilized_aut;
            for(int k=0;k<automorphisms.size();++k){
                vector<uint> &aut = automorphisms[k];
                if(aut[vid]==vid){
                    stabilized_aut.push_back(aut);
                }else{
                    uint key = aut[vid];
                    map<uint,uint>::iterator ite;
                    if((ite=restricts.find(key))==restricts.end()){
                        restricts.insert(pair<uint,uint>(key,vid));
                    }else{
                        uint tmp = ite->second;
                        if(tmp<vid){
                            ite->second = vid;
                        }
                    }
                }
            }
            automorphisms = stabilized_aut;
        }
    }
    //for each vertex, we give a restriction. thus we do not need to judge
    //whether the vertex exits in restricts. if the value is 0xffffffff, then
    //no restricts exist for vertex key
    for(int i=0;i<=vnum;++i){
        if(restricts.find(i)==restricts.end()){
            restricts.insert(pair<uint,uint>(i,0xffffffff));
        }
    }
}

void rearrangeReducePhase(map<uint64_t,uint>&extendPhaseKeyCount, map<uint64_t,uint>&othersKeyCount,map<uint64_t,phaseKeyCount>&reducePhaseKeyCount,
                          map<uint64_t,uint>&repeatFlag,vector<Phase>&matchOrder,map<uint,map<uint,bool>>&recordLabelStat,map<uint,uint>&queryLabel){
    map<uint64_t,phaseKeyCount>::iterator maxite;
    uint itercount = reducePhaseKeyCount.size();
    for(uint j=0;j<itercount;++j) {
        uint precount = 0;
        for (auto ite = reducePhaseKeyCount.begin(); ite != reducePhaseKeyCount.end(); ++ite) {
            uint64_t key = ite->first;
            phaseKeyCount &tmp = ite->second;
            uint count = 0;
            if (extendPhaseKeyCount.find(key) != extendPhaseKeyCount.end()) {
                count = count + extendPhaseKeyCount[key];
            }
            if (othersKeyCount.find(key) != othersKeyCount.end()) {
                count = count + othersKeyCount[key];
            }
            count = count + tmp.count;
            if (count > precount) {
                maxite = ite;
                precount = precount;
            }
        }
        uint64_t key = maxite->first;
        phaseKeyCount &keycount = maxite->second;
        vector<uint64_t> &keylist = keycount.keyList;
        for (uint i = 0; i < keylist.size(); ++i) {
            uint64_t tmpkey = keylist[i];
            if ((maxite = reducePhaseKeyCount.find(key)) == reducePhaseKeyCount.end()) {
                cout << "wrong in reducePhaseKeyCount" << endl;
                exit(0);
            }
            if (maxite->second.count > 0) {
                maxite->second.count = maxite->second.count - 1;
            }
        }
        vector<uint> &phaselist = keycount.phaseList;
        for (uint i = 0; i < phaselist.size() / 2; ++i) {
            uint phaseIndex = phaselist[i * 2];
            uint pairIndex = phaselist[i * 2 + 1];
            Phase &phase = matchOrder[phaseIndex];
            vector<uint> &pairs = phase.pairs;
            uint vid1 = pairs[pairIndex * 2];
            uint vid2 = pairs[pairIndex * 2 + 1];
            uint64_t tmp = 0;
            if ((tmp | queryLabel[vid1]) == (key & 0x00000000ffffffff)) {
                pairs[pairIndex * 2] = vid2;
                pairs[pairIndex * 2 + 1] = vid1;
            }
        }
        map<uint64_t,uint>::iterator tmpite;
        if((tmpite=repeatFlag.find(key))==repeatFlag.end()){
            uint tmpkeycount = reducePhaseKeyCount[key].count;
            tmpkeycount = (tmpkeycount<<16)|tmpkeycount;
            repeatFlag.insert(pair<uint64_t,uint>(key,tmpkeycount));
            if(reducePhaseKeyCount[key].count>1){
                uint vlabel = (key&0x00000000ffffffff);
                uint edgeLabel = ((key>>32)&0x00000000ffffffff);
                if(recordLabelStat.find(edgeLabel)==recordLabelStat.end()){
                    map<uint,bool> tmpmap;
                    tmpmap.insert(pair<uint,bool>(vlabel,true));
                    recordLabelStat.insert(pair<uint,map<uint,bool>>(edgeLabel,tmpmap));
                }else{
                    auto ite = recordLabelStat.find(edgeLabel);
                    map<uint,bool> &tmpmap = ite->second;
                    if(tmpmap.find(vlabel)==tmpmap.end()){
                        tmpmap.insert(pair<uint,bool>(vlabel,true));
                    }
                }
            }
        }else{
            cout<<"repeated key"<<endl;
            exit(0);
        }
        reducePhaseKeyCount.erase(key);
    }
    for(auto ite=extendPhaseKeyCount.begin();ite!=extendPhaseKeyCount.end();++ite){
        uint64_t key = ite->first;
        uint count = ite->second;
        map<uint64_t,uint>::iterator tmpite;
        if((tmpite=repeatFlag.find(key))==repeatFlag.end()){
            uint tmpkeycount = count;
            tmpkeycount = (tmpkeycount<<16)|tmpkeycount;
            repeatFlag.insert(pair<uint64_t,uint>(key,tmpkeycount));
            if(count>1){
                uint vlabel = (key&0x00000000ffffffff);
                uint edgeLabel = ((key>>32)&0x00000000ffffffff);
                if(recordLabelStat.find(edgeLabel)==recordLabelStat.end()){
                    map<uint,bool> tmpmap;
                    tmpmap.insert(pair<uint,bool>(vlabel,true));
                    recordLabelStat.insert(pair<uint,map<uint,bool>>(edgeLabel,tmpmap));
                }else{
                    auto ite = recordLabelStat.find(edgeLabel);
                    map<uint,bool> &tmpmap = ite->second;
                    if(tmpmap.find(vlabel)==tmpmap.end()){
                        tmpmap.insert(pair<uint,bool>(vlabel,true));
                    }
                }
            }
        }else {
            tmpite->second = tmpite->second+count;
            if(tmpite->second>1){
                uint vlabel = (key&0x00000000ffffffff);
                uint edgeLabel = ((key>>32)&0x00000000ffffffff);
                if(recordLabelStat.find(edgeLabel)==recordLabelStat.end()){
                    map<uint,bool> tmpmap;
                    tmpmap.insert(pair<uint,bool>(vlabel,true));
                    recordLabelStat.insert(pair<uint,map<uint,bool>>(edgeLabel,tmpmap));
                }else{
                    auto ite = recordLabelStat.find(edgeLabel);
                    map<uint,bool> &tmpmap = ite->second;
                    if(tmpmap.find(vlabel)==tmpmap.end()){
                        tmpmap.insert(pair<uint,bool>(vlabel,true));
                    }
                }
            }
        }
    }
    for(auto ite=othersKeyCount.begin();ite!=othersKeyCount.end();++ite){
        uint64_t key = ite->first;
        uint count = ite->second;
        map<uint64_t,uint>::iterator tmpite;
        if((tmpite=repeatFlag.find(key))==repeatFlag.end()){
            uint tmpkeycount = count;
            tmpkeycount = (tmpkeycount<<16)|tmpkeycount;
            repeatFlag.insert(pair<uint64_t,uint>(key,tmpkeycount));
            if(count>1){
                uint vlabel = (key&0x00000000ffffffff);
                uint edgeLabel = ((key>>32)&0x00000000ffffffff);
                if(recordLabelStat.find(edgeLabel)==recordLabelStat.end()){
                    map<uint,bool> tmpmap;
                    tmpmap.insert(pair<uint,bool>(vlabel,true));
                    recordLabelStat.insert(pair<uint,map<uint,bool>>(edgeLabel,tmpmap));
                }else{
                    auto ite = recordLabelStat.find(edgeLabel);
                    map<uint,bool> &tmpmap = ite->second;
                    if(tmpmap.find(vlabel)==tmpmap.end()){
                        tmpmap.insert(pair<uint,bool>(vlabel,true));
                    }
                }
            }
        }else {
            tmpite->second = tmpite->second+count;
            if(tmpite->second>1){
                uint vlabel = (key&0x00000000ffffffff);
                uint edgeLabel = ((key>>32)&0x00000000ffffffff);
                if(recordLabelStat.find(edgeLabel)==recordLabelStat.end()){
                    map<uint,bool> tmpmap;
                    tmpmap.insert(pair<uint,bool>(vlabel,true));
                    recordLabelStat.insert(pair<uint,map<uint,bool>>(edgeLabel,tmpmap));
                }else{
                    auto ite = recordLabelStat.find(edgeLabel);
                    map<uint,bool> &tmpmap = ite->second;
                    if(tmpmap.find(vlabel)==tmpmap.end()){
                        tmpmap.insert(pair<uint,bool>(vlabel,true));
                    }
                }
            }
        }
    }
}

//we can conver graph mining into fp32 and uint32 mixed approach
void genMatchOrder(string queryFileName, vector<Phase>&matchOrder, map<uint64_t,uint> &repeatFlag,map<uint,uint>&restricts,
                   map<uint,map<uint,bool>>&recordLabelStat, map<uint,uint> &queryLabel,map<uint,uint*>&allEdgeLabelPartitions,
                   uint *baseAddr_dev, uint &stride) {
    char type;
    ifstream queryfile(queryFileName, ios::in);
    if (!queryfile.is_open()) {
        cout << "can not open query file" << endl;
        return;
    }
    //read vertices of the query pattern
    while (queryfile.get(type)) {
        if (type == 'v') {
            uint vid, vlabel;
            queryfile >> vid >> vlabel;
            map<uint, uint>::iterator ite = queryLabel.find(vid);
            if (ite != queryLabel.end()) {
                cout << "duplicate vertex " << vid << endl;
            } else {
                queryLabel.insert(pair<uint, uint>(vid, vlabel));
            }
        } else if (type == 'e') {
            break;
        } else if( (type=='#' || type=='t') || (type>='0' && type<='9')){
            char tmpc[100];
            queryfile.getline(tmpc,100);
        }
    }
    queryfile.unget();
    vnum = queryLabel.size();
    //record neighbors information of each vertex in the query pattern
    vector<queryNeigInfo> query(vnum);
    vector<queryEdge> edges;
    adjMatrix.resize(vnum*vnum,0xffffffff);
    vertexNeig.resize(vnum,0);
    //each vertex label or edge label is represented
    //with 20 bits. svid label+edge label+evid label contitute a uint64_t integer
    //each distinct uint64_t integer has a unique id
    map<uint64_t,uint> vevLabel2Id;
    //read edges of query pattern
    while (queryfile.get(type)) {
        if(type=='e'){
            uint svid, evid, elabel;
            queryfile >> svid >> evid >> elabel;
            //cout << svid <<" "<< evid <<" "<< elabel <<endl;
            addToQueryGraph(svid,evid,elabel,query);
            addToQueryGraph(evid,svid,elabel,query);
            uint64_t key = queryLabel[svid];
            key = (key<<20) | elabel;
            key = (key<<20)|queryLabel[evid];
            //repeatFlag[key] = 0;
            if(vevLabel2Id.find(key)==vevLabel2Id.end()){
                vevLabel2Id[key] = vevLabel2Id.size();
            }
            key = queryLabel[evid];
            key = (key<<20) | elabel;
            key = (key<<20)|queryLabel[svid];
            //repeatFlag[key] = 0;
            if(vevLabel2Id.find(key)==vevLabel2Id.end()){
                vevLabel2Id[key] = vevLabel2Id.size();
            }
            if(adjMatrix[svid*vnum+evid]==0xffffffff){
                    queryEdge tmpedge;
                    if(svid<evid) {tmpedge.svid=svid; tmpedge.evid = evid;}
                    else {tmpedge.svid = evid; tmpedge.evid = svid;}
                    tmpedge.edgeLabel = elabel;
                    edges.push_back(tmpedge);
                    adjMatrix[svid*vnum+evid] = elabel;
                    adjMatrix[evid*vnum+svid] = elabel;
            }
            uint tmp = vertexNeig[svid];
            vertexNeig[svid] = tmp | (1<<evid);
            tmp = vertexNeig[evid];
            vertexNeig[evid] = tmp | (1<<svid);
        }
    }
    queryfile.close();


    map<uint64_t,uint> extendPhaseKeyCount;
    map<uint64_t,phaseKeyCount> reducePhaseKeyCount;
    map<uint64_t,uint> othersKeyCount;
    //in the next, we need to generate access sequence of vertices
    //to generate first two or three vertices, we need to justify the existence of triangles
    vector<triInfo> triangleCount(edges.size());
    for(uint i=0;i<edges.size();++i){
        countTriangle(i,triangleCount,edges,queryLabel);
    }
    //sort results by the number of triangles in descrending order
    sort(triangleCount.begin(),triangleCount.end(),sortDescend);
    vidToIndex.resize(vnum,-1);
    if(triangleCount[0].triCount==0){
        //no triangles
        uint svid,extvid1,extvid2;
        genInitialPartialQueryNoTriangle(query,queryLabel,svid,extvid1,extvid2);
        uint edgelabel = adjMatrix[svid*vnum+extvid1];
        uint *edgeLabelPartition = allEdgeLabelPartitions[edgelabel];
        uint len = edgeLabelPartition[0]-5;
        cudaMemcpyAsync(baseAddr_dev+stride,edgeLabelPartition+5,sizeof(uint)*len,cudaMemcpyHostToDevice);
        stride = stride+len;
        vidToIndex[svid] = 0;
        partialQueryVid.push_back(svid);
        partialQueryVidMap = partialQueryVidMap | (1<<svid);
        addExtPhase(svid,extvid1,svid,extvid2,extvid2==0xffffffff?1:2,extendPhaseKeyCount,matchOrder,queryLabel);
    }else{
        //we need to select the vertex that has at least two vertices with the same label
        uint tricount = triangleCount[0].triCount;
        uint i=1,prei=0;
        while(i<triangleCount.size() && triangleCount[i].triCount==tricount){
            if(triangleCount[i].triScore>triangleCount[prei].triScore){ prei = i; }
            ++i;
        }
        uint score = triangleCount[prei].triScore;
        uint edgeindex = triangleCount[prei].edgeIndex;
        uint svid = edges[edgeindex].svid;
        uint evid = edges[edgeindex].evid;
        uint thirdid = triangleCount[prei].maxScoreVID;
        Phase extphase1, extphase2, reducephase3;
        //at least two edges have the same label
        if(score==33 || score==31 || score==30 || (score>=10 && score<=13)){
            uint sourcvid=svid,extvid1=evid,extvid2=thirdid;
            if(score==31){
                if(queryLabel[svid]==queryLabel[evid]){ sourcvid = thirdid; extvid1=svid; extvid2=evid; }
                else if(queryLabel[svid]==queryLabel[thirdid]){ sourcvid = evid; extvid1=svid; extvid2=thirdid; }
                else if(queryLabel[evid]==queryLabel[thirdid]){ sourcvid = svid; extvid1=evid; extvid2=thirdid; }
            }else if(score>=10 && score<=13){
                uint e1l = adjMatrix[svid * vnum + thirdid], e2l = adjMatrix[evid * vnum + thirdid], e3l = edges[edgeindex].edgeLabel;
                if(e1l==e2l) { sourcvid = thirdid; extvid1 = svid; extvid2=evid; }
                else if(e1l==e3l) { sourcvid = svid; extvid1=evid; extvid2 = thirdid;}
                else if(e2l==e3l) { sourcvid = evid; extvid1=svid; extvid2 = thirdid;}
            }
            vidToIndex[sourcvid] = 0;
            partialQueryVid.push_back(sourcvid);
            partialQueryVidMap = partialQueryVidMap | (1<<sourcvid);
            uint edgelabel = adjMatrix[sourcvid*vnum+extvid1];
            uint *edgeLabelPartition = allEdgeLabelPartitions[edgelabel];
            uint len = edgeLabelPartition[0]-5;
            cudaMemcpyAsync(baseAddr_dev+stride,edgeLabelPartition+5,sizeof(uint)*len,cudaMemcpyHostToDevice);
            stride = stride+len;
            addExtPhase(sourcvid,extvid1,sourcvid,extvid2,2,extendPhaseKeyCount,matchOrder,queryLabel);
            vector<uint> pairs;
            pairs.push_back(extvid1);
            pairs.push_back(extvid2);
            addReducePhase(pairs,reducePhaseKeyCount,othersKeyCount,matchOrder,queryLabel);
        }else if(score<=3){
            uint sourcvid=svid,extvid1=evid,extvid2=thirdid;
            vidToIndex[sourcvid] = 0;
            partialQueryVid.push_back(sourcvid);
            partialQueryVidMap = partialQueryVidMap | (1<<sourcvid);
            uint edgelabel = adjMatrix[sourcvid*vnum+extvid1];
            uint *edgeLabelPartition = allEdgeLabelPartitions[edgelabel];
            uint len = edgeLabelPartition[0]-5;
            cudaMemcpyAsync(baseAddr_dev+stride,edgeLabelPartition+5,sizeof(uint)*len,cudaMemcpyHostToDevice);
            stride = stride+len;
            addExtPhase(sourcvid,extvid1,sourcvid,extvid2,1,extendPhaseKeyCount,matchOrder,queryLabel);
            addExtPhase(sourcvid,extvid2,sourcvid,extvid2,1,extendPhaseKeyCount,matchOrder,queryLabel);
            vector<uint> pairs;
            pairs.push_back(extvid1);
            pairs.push_back(extvid2);
            addReducePhase(pairs,reducePhaseKeyCount,othersKeyCount,matchOrder,queryLabel);
        }else{
            cout<<"wrong"<<endl;
        }
    }
    //now we need to decide how to iterate over partial vertices to see if there are back edges
    list<uint>::iterator itelist;
    //construct match order
    while(partialQueryVid.size()<vnum) {
        vector<NeigInfo> neigPhaseInfo;
        //evids of partial svids, all evids are grouped by edge label
        map<uint,vector<uint>> sameEdgeLabelEvid;
        for (itelist = partialQueryVid.begin(); itelist != partialQueryVid.end(); itelist++) {
            uint svid = (*itelist);
            uint neigs = vertexNeig[svid];
            map<uint, vector<uint> > &totEdgeLabelNeig = query[svid].edgeLabelNeigs;
            //iterate over neighbors (not in partialQueryVid) of svid
            for (uint i = 0; i < 32; ++i) {
                //obtain the vid of the neghbor, the position of non-zero bit is the neighbor vid
                uint neigId = neigs & (1 << i);
                if (neigId == 0) { continue; }
                neigId = i;
                //find out if the neigId has back edges with partialQueryVidMap
                uint tmp = vertexNeig[neigId] & partialQueryVidMap;
                tmp = __builtin_popcount(tmp) - 1;
                NeigInfo tmpneig;
                tmpneig.svid = svid;
                tmpneig.evid = neigId;
                tmpneig.edgeLabel = adjMatrix[svid * vnum + neigId];
                tmpneig.commonVNum = tmp;
                tmpneig.evLabel = queryLabel[neigId];
                uint tmpEdgeLabelNeig = totEdgeLabelNeig[tmpneig.edgeLabel][1];
                tmpneig.evNum = __builtin_popcount(neigs & tmpEdgeLabelNeig);
                neigPhaseInfo.push_back(tmpneig);

                uint isvalidevid = (1<<neigId) & partialQueryVidMap;
                if(isvalidevid==0){
                    auto ite = sameEdgeLabelEvid.find(tmpneig.edgeLabel);
                    if(ite==sameEdgeLabelEvid.end()){
                        vector<uint> tmpv;
                        tmpv.push_back(svid);
                        tmpv.push_back(neigId);
                        sameEdgeLabelEvid.insert(pair<uint,vector<uint>>(tmpneig.edgeLabel,tmpv));
                    }else{
                        vector<uint> &tmpv = ite->second;
                        tmpv.push_back(svid);
                        tmpv.push_back(neigId);
                    }
                }
            }
        }
        sort(neigPhaseInfo.begin(), neigPhaseInfo.end(), multiLevelSort);
        //if there has back edges
        if (neigPhaseInfo[0].commonVNum > 0) {
            uint sourcvid = neigPhaseInfo[0].svid, extvid1 = neigPhaseInfo[0].evid;
            addExtPhase(sourcvid, extvid1, 0, 0, 1,extendPhaseKeyCount,matchOrder,queryLabel);
            neigPhaseInfo.clear();
            uint neigpos = vertexNeig[extvid1] & partialQueryVidMap;
            for (uint i = 0; i < 32; ++i) {
                if ((neigpos & (1 << i))>0) {
                    NeigInfo tmpneig;
                    tmpneig.svid = extvid1;
                    tmpneig.evid = i;
                    tmpneig.edgeLabel = adjMatrix[extvid1 * vnum + i];
                    tmpneig.evLabel = queryLabel[i];
                    neigPhaseInfo.push_back(tmpneig);
                }
            }
            sort(neigPhaseInfo.begin(), neigPhaseInfo.end(), [](NeigInfo &a, NeigInfo &b){
                if(a.edgeLabel<b.edgeLabel){
                    return true;
                }else if(a.edgeLabel==b.edgeLabel){
                    return a.evLabel<b.evLabel;
                }
                return false;
            });
            uint prei = 0;
            vector<uint> prev;
            for (uint i = 0; i < neigPhaseInfo.size(); ++i) {
                if (neigPhaseInfo[i].edgeLabel == neigPhaseInfo[prei].edgeLabel) {
                    prev.push_back(neigPhaseInfo[i].svid);
                    prev.push_back(neigPhaseInfo[i].evid);
                } else {
                    addReducePhase(prev,reducePhaseKeyCount,othersKeyCount,matchOrder,queryLabel);
                    prev.clear();
                    prei = i;
                    --i;
                }
            }
            addReducePhase(prev,reducePhaseKeyCount,othersKeyCount,matchOrder,queryLabel);
        }else{
            //find out if there are two extvids that have the same edge label and there is an edge
            //between two extvids.
            bool found = false;
            for(auto ite=sameEdgeLabelEvid.begin();ite!=sameEdgeLabelEvid.end();++ite){
                vector<uint> &evids = ite->second;
                //reduce pairs, extend pairs
                for(uint i=0;i<evids.size()/2;++i){
                    uint evid1 = evids[i*2+1];
                    for(uint j=i+1;j<evids.size()/2;++j){
                        uint evid2 = evids[j*2+1];
                        if(adjMatrix[evid1*vnum+evid2]!=0xffffffff){
                            found = true;
                            uint sourcvid1 = evids[i*2], extvid1 = evid1;
                            uint sourcvid2 = evids[j*2], extvid2 = evid2;
                            addExtPhase(sourcvid1, extvid1, sourcvid2, extvid2, 2,extendPhaseKeyCount,matchOrder,queryLabel);
                            vector<uint> reducpairs;
                            reducpairs.push_back(extvid1);
                            reducpairs.push_back(extvid2);
                            addReducePhase(reducpairs,reducePhaseKeyCount,othersKeyCount,matchOrder,queryLabel);
                            break;
                        }
                    }
                    if(found){
                        break;
                    }
                }
                if(found){
                    break;
                }
            }
            //TODO: when extend two vertices, we need to judge if there is an edge between two vertices
            if(!found){
                if (neigPhaseInfo[0].evNum > 1) {
                    //no back edges, we choose the vertex that has the most same edge label neighbors,1 src
                    uint sourcvid1 = neigPhaseInfo[0].svid, extvid1 = neigPhaseInfo[0].evid;
                    uint sourcvid2 = neigPhaseInfo[1].svid, extvid2 = neigPhaseInfo[1].evid;
                    addExtPhase(sourcvid1, extvid1, sourcvid2, extvid2, 2,extendPhaseKeyCount,matchOrder,queryLabel);
                } else {
                    bool found=false;
                    for(uint i=0;i<neigPhaseInfo.size()-1;++i){
                        if(neigPhaseInfo[i].edgeLabel==neigPhaseInfo[i+1].edgeLabel){
                            uint sourcvid1 = neigPhaseInfo[i].svid, extvid1 = neigPhaseInfo[i].evid;
                            uint sourcvid2 = neigPhaseInfo[i+1].svid, extvid2 = neigPhaseInfo[i+1].evid;
                            addExtPhase(sourcvid1, extvid1, sourcvid2, extvid2, 2,extendPhaseKeyCount,matchOrder,queryLabel);
                            found=true;
                            break;
                        }
                    }
                    if(!found){
                        uint sourcvid = neigPhaseInfo[0].svid, extvid1 = neigPhaseInfo[0].evid;
                        addExtPhase(sourcvid, extvid1, 0, 0, 1,extendPhaseKeyCount,matchOrder,queryLabel);
                    }
                }
            }
        } 

    }
    //using another thread to detect symmetries
    vector<uint> autCandidate;
    vector<vector<uint>> autResults;
    genAutomorphism(autCandidate, true,query, autResults,queryLabel);
    rearrangeReducePhase(extendPhaseKeyCount,othersKeyCount,reducePhaseKeyCount,repeatFlag,matchOrder,recordLabelStat,queryLabel);
    //above codes are used to generate access sequence
    //next, we need to generate restrictions for vertices based on the access sequence and automorphisms
    //map<uint,uint> restricts;
    /*cout<<"automorphisms count: "<<autResults.size()<<endl;
    for(int i=0;i<autResults.size();++i){
        vector<uint> &tmpV = autResults[i];
        for(int j=0;j<tmpV.size();++j){
            cout<<tmpV[j] <<" ";
        }
        cout<<endl;
    }
    cout<<"-------------------------------"<<endl;*/
    genRestriction(autResults, restricts,matchOrder);
    return;
}
