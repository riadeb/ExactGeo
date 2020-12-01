#include <igl/opengl/glfw/Viewer.h>

#include <list>
#include <unordered_map>
#include <utility>
#include <vector>
#include <queue>
#include <__nullptr>



using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

class Edge;

class Window {
    public :
    Edge* edge;
    double p1,p2;
    double d1,d2,sigma;
    bool tau; //true if in order of edge->v1 , edge->v2
    Window(Edge* e, double v1, double v2, double dd1, double dd2, double sig, bool _tau){
        edge = e; p1 = v2; p2 = v2; d1 = dd1; d2 = dd2; sigma = sig; tau = _tau;
    }

};
class Edge { //edge (oriented from smaller vertex to bigger vertex)

    public :
        int v1,v2; //v1 smaller that v2
        list<Window>  windows;
        Edge(int u1, int u2){
            v1 = min(u1,u2);
            v2 = max(u1,u2);
            windows = list<Window>();
        }
};

class PairHash
{
public:
  size_t operator()(const pair<int,int> &p) const
  {
    return p.first * p.second;
  }
};

class EdgeDS {

    public :

    vector<vector<Edge*> > edges;
    unordered_map<pair<int,int>, vector<pair<int,int> >, PairHash> adjL;
    MatrixXd V; MatrixXi F;

    EdgeDS(MatrixXd _V, MatrixXi _F){
        V = _V;
        F = _F;
        int n = V.rows();
        int f = F.rows();
        edges = vector<vector<Edge*> > (n, vector<Edge*>(n, nullptr));
        for(int i = 0; i < f; i++){
            for(int j = 0; j < 3; j++){
                  int v1 = F(i, j);
                  int v2 = F(i, (j+1)%3);
                  if(edges[v1][v2] == nullptr){
                    Edge* ed = new  Edge(v1,v2);
                     edges[v1][v2] = ed;
                     edges[v2][v1] = ed;
                  }
                 
            }
        }
         for(int i = 0; i < f; i++){
            for(int j = 0; j < 3; j++){
                  int v1 = F(i, j);
                  int v2 = F(i, (j+1)%3);
                  int v3 = F(i, (j+2)%3);
                  adjL[make_pair(v1,v2)].push_back(make_pair(v2,v3));
                    adjL[make_pair(v1,v2)].push_back(make_pair(v3,v1));
                 
            }

        }

    }
    ~EdgeDS(){
        int n = edges.size();
        for(int i = 0; i < n-1; i++){
            for(int j = i+1; j < n ;j++) delete edges[i][j];
        }
    }
    double dist(int v1, int v2){
        return (V.row(v1) - V.row(v2)).norm();
    }
    void exactGeo(int s){
            queue<Window> q;
            int n = edges.size();
            for(int i = 0; i < n; i++){ //initialize queue
                if(edges[s][i] != nullptr){
                    if(s < i){
                         Window w = Window(edges[s][i], 0, dist(s,i), 0, dist(s,i),0,true);
                         edges[s][i]->windows.push_back(w);
                        q.push(w);
                    }
                    else{
                         Window w = Window(edges[s][i], 0, dist(s,i), dist(s,i),0,0,false);
                         edges[s][i]->windows.push_back(w);
                         q.push(w);
                    }
                    for(auto e : adjL[make_pair(s,i)]){
                        if(e.first != s && e.second != s){
                            Edge* edge =  edges[e.first][e.second];
                            if(edge->windows.size() == 0){
                                bool tau = e.first < e.second; //direction of window
                                Window w = Window(edge,0, dist(edge->v1,edge->v2), dist(s,edge->v1), dist(s, edge->v2), 0,tau);
                                q.push(w);
                                edges[s][i]->windows.push_back(w);
                            }
                        }
                    }
                }
            }

        while(!q.empty()){
            Window w = q.front();
            q.pop();
            

        }

    }


    

};


