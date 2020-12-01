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
    Window(Edge* e, int v1, int v2, double dd1, double dd2, double sig, bool _tau){
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


class EdgeDS {

    public :

    vector<vector<Edge*> > edges;
    vector<vector<vector<Edge*> > > adjL;
    MatrixXd V, MatrixXi F;

    EdgeDS(MatrixXd _V, MatrixXi _F){
        V = _V;
        F = _F;
        int n = V.rows();
        int f = F.rows();
        edges = vector<vector<Edge*> > (n, vector<Edge*>(n, nullptr));
        adjL =  vector<vector<vector<Edge*> > > (n ,vector<vector<Edge*> > (n));
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
                  if(v1 < v2){ //the face is in the positive side of the edge
                    adjL[v1][v2][true].pushback(edges[v2][v3]);
                    adjL[v1][v2][true].pushback(edges[v1][v3]);
                  }
                  else{ //the face is in the negative side of the edge
                    adjL[v2][v1][false].pushback(edges[v2][v3]);
                    adjL[v2][v2][false].pushback(edges[v1][v3]);
                  }
                 
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
            for(int i = 0; i < n; i++){
                if(edges[s][i] != nullptr){
                    Window w;
                    if(s < i){
                         w = Window(edges[s][i], 0, dist(s,i), 0, dist(s,i),0,true);
                    }
                    else{
                         w = Window(edges[s][i], 0, dist(s,i), dist(s,i),0,0,false);
                    }
                    edges[s][i]->windows.push_back(w);
                    q.push(w);
                    for(auto e : adjL[edges[s][i]->v1][edges[s][i]->v2][0]){
                        fdd
                    }
                }
            }

    }


    

};


