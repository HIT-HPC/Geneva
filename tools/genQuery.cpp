#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>
#include <map>
#include <ctime>
#include <cstdlib>

std::vector<std::map<uint,std::map<uint,std::vector<uint>>>> reorderGraph;
std::vector<std::vector<uint>> allsvids;
std::map<uint,uint> vLabel;

void genQuery(uint vertexNum, uint maxEdgeNum, uint totreorderlabelpart, std::vector<uint> &queryV,std::vector<uint> &edges, std::map<uint,uint> &reorderV){
    uint seed;
    seed = std::time(0);
    std::srand(seed);
    std::map<uint,int> vidFlag;

    uint edgelabel = std::rand()%totreorderlabelpart;
    std::map<uint,std::map<uint,std::vector<uint>>> &labelpart = reorderGraph[edgelabel];
    std::vector<int> adj(vertexNum*vertexNum,-1);
    std::vector<uint> &svidcollect = allsvids[edgelabel];
    uint svidnum = svidcollect.size();
    uint svidpos = std::rand()%svidnum;
    uint svid = svidcollect[svidpos];
    queryV.push_back(svid);
    vidFlag[svid] = 1;
    reorderV[svid] = 0;
    uint i=1;
    while(i<vertexNum){
        std::map<uint,std::map<uint,std::vector<uint>>> &labelpart = reorderGraph[edgelabel];
        std::map<uint,std::vector<uint>> &labelneigs = labelpart[svid];
        uint evidlabelpos = std::rand()%labelneigs.size();
        uint index=0;
        bool foundNewV=false;
        for(auto ite=labelneigs.begin();ite!=labelneigs.end();++ite){
            if(index==evidlabelpos){
                std::vector<uint> &neigs = ite->second;
                uint evidpos = std::rand()%neigs.size();
                uint evid = neigs[evidpos];
                auto itea = vidFlag.find(evid);
                if(itea!=vidFlag.end()){
                    uint querysvid = reorderV[svid];
                    uint queryevid = reorderV[evid];
                    if(adj[querysvid*vertexNum+queryevid]==-1){
                        edges.push_back(svid);
                        edges.push_back(evid);
                        edges.push_back(edgelabel);
                        adj[querysvid*vertexNum+queryevid] = 1;
                        adj[queryevid*vertexNum+querysvid] = 1;
                    }
                }else{
                    queryV.push_back(evid);
                    vidFlag[evid] = 1;
                    reorderV[evid] = reorderV.size();
                    foundNewV = true;
                    uint querysvid = reorderV[svid];
                    uint queryevid = reorderV[evid];
                    if(adj[querysvid*vertexNum+queryevid]==-1){
                        edges.push_back(svid);
                        edges.push_back(evid);
                        edges.push_back(edgelabel);
                        adj[querysvid*vertexNum+queryevid] = 1;
                        adj[queryevid*vertexNum+querysvid] = 1;
                    }
                }
                break;
            }
        }
        if(foundNewV==true){
            i++;
        }
        svidpos=std::rand()%queryV.size();
        svid = queryV[svidpos];
        std::vector<uint> tmp;
        for(uint j=0;j<reorderGraph.size();++j){
            std::map<uint,std::map<uint,std::vector<uint>>> &tmpmap = reorderGraph[j];
            if(tmpmap.find(svid)!=tmpmap.end()){
                tmp.push_back(j);
            }
        }
        uint edgelabelpos = std::rand()%tmp.size();
        edgelabel = tmp[edgelabelpos];
    }
    if(edges.size()/3<maxEdgeNum){
        for(i=0;i<queryV.size();++i){
            for(uint j=i+1;j<queryV.size();++j){
                uint querysvid = reorderV[queryV[i]];
                uint queryevid = reorderV[queryV[j]];
                uint evidlable = vLabel[queryV[j]];
                if(adj[querysvid*vertexNum+queryevid]==-1){
                    for(uint k=0;k<reorderGraph.size();++k){
                        std::map<uint,std::map<uint,std::vector<uint>>> &tmpmap = reorderGraph[k];
                        if(tmpmap.find(queryV[i])!=tmpmap.end()){
                            std::map<uint,std::vector<uint>> &labelneigs = tmpmap[queryV[i]];
                            if(labelneigs.find(evidlable)!=labelneigs.end()){
                                std::vector<uint> &neigs = labelneigs[evidlable];
                                for(uint l=0;l<neigs.size();++l){
                                    if(neigs[l]==queryV[j]){
                                        if(adj[querysvid*vertexNum+queryevid]==-1){
                                            edges.push_back(queryV[i]);
                                            edges.push_back(queryV[j]);
                                            edges.push_back(k);
                                            adj[querysvid*vertexNum+queryevid]==1;
                                            adj[queryevid*vertexNum+querysvid]==1;
                                        }
                                        if(edges.size()/3>=maxEdgeNum){
                                            break;
                                        }
                                    }
                                }
                                if(edges.size()/3>=maxEdgeNum){
                                    break;
                                }
                            }
                        }
                    }
                    if(edges.size()/3>=maxEdgeNum){
                        break;
                    }
                }
            }
            if(edges.size()/3>=maxEdgeNum){
                break;
            }
        }
    }
}

//gen query_4_5_0 wrong
int main(int argc, char *argv[]){
    std::string path = "../datasets/";
    std::string databaseFileName = std::string(argv[1]);
    int vertexNum = atoi(argv[2]);
    int maxEdgeNum = atoi(argv[3]);
    int queryNum = atoi(argv[4]);
    
    std::ifstream reorderfile(databaseFileName+".reorder",std::ios::in);
    if(!reorderfile.is_open()){
        std::cout << "open " << "reorder" << " is wrong" << std::endl;
        return 1;
    }
    uint totreorderlabelpart;
    reorderfile>> totreorderlabelpart;
    for(uint i=0;i<totreorderlabelpart;++i){
        uint svidnum;
        reorderfile >> svidnum;
        std::map<uint,std::map<uint,std::vector<uint> > > labelpart;
        std::vector<uint> svidcollect;
        for(uint j=0;j<svidnum;++j){
            uint svid, evidlabelnum;
            reorderfile>>svid>>evidlabelnum;
            std::map<uint,std::vector<uint>> evidlabelneigs;
            for(uint k=0;k<evidlabelnum;++k){
                uint evidlabel, neigsnum;
                reorderfile>>evidlabel>>neigsnum;
                std::vector<uint> neigs;
                for(uint l=0;l<neigsnum;++l){
                    uint neigvid;
                    reorderfile>>neigvid;
                    neigs.push_back(neigvid);
                    auto ite=vLabel.find(neigvid);
                    if(ite!=vLabel.end()){
                        uint label = vLabel[neigvid];
                        if(label!=evidlabel){
                            std::cout << "wrong labels"<<std::endl;
                        }
                    }else{
                        vLabel.insert(std::pair<uint,uint>(neigvid,evidlabel));
                    }
                }
                evidlabelneigs[evidlabel] = neigs;
            }
            labelpart[svid] = evidlabelneigs;
            svidcollect.push_back(svid);
        }
        reorderGraph.push_back(labelpart);
        allsvids.push_back(svidcollect);
    }
    reorderfile.close();

    for(int i=0;i<queryNum;++i){
        std::vector<uint> queryV;
        std::vector<uint> edges;
        std::map<uint,uint> reorderV;
        genQuery(vertexNum,maxEdgeNum,totreorderlabelpart,queryV,edges,reorderV);
        std::ofstream queryfile(path+databaseFileName+"_"+std::to_string(vertexNum)+"_"+std::to_string(maxEdgeNum)+"_"+std::to_string(i)+".query",std::ios::out);
        for(int j=0;j<queryV.size();++j){
            queryfile<<"#v "<<queryV[j]<<" "<<vLabel[queryV[j]]<<std::endl;
        }
        for(int j=0;j<edges.size()/3;++j){
            uint querysvid = edges[j*3+0];
            uint queryevid = edges[j*3+1];
            uint edgelabel = edges[j*3+2];
            queryfile<<"#e "<<querysvid<<" "<<queryevid<<" "<<edgelabel<<std::endl;
        }
        for(int j=0;j<queryV.size();++j){
            queryfile<<"v "<<reorderV[queryV[j]]<<" "<<vLabel[queryV[j]]<<std::endl;
        }
        for(int j=0;j<edges.size()/3;++j){
            uint querysvid = edges[j*3+0];
            uint queryevid = edges[j*3+1];
            uint edgelabel = edges[j*3+2];
            queryfile<<"e "<<reorderV[querysvid]<<" "<<reorderV[queryevid]<<" "<<edgelabel<<std::endl;
        }
        queryfile.close();
    }

    return 0;
}
