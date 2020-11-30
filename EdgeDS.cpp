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
    int p1,p2;
    double d1,d2,sigma;
    int face; //true if in order of edge->v1 , edge->v2
    Window(Edge* e, int v1, int v2, double dd1, double dd2, double sig, int face_){
        edge = e; p1 = v2; p2 = v2; d1 = dd1; d2 = dd2; sigma = sig; face =face_;
    }

};
class Edge {

    public :
        int v1,v2; //v1 smaller that v2
        list<Window>  windows;
        vector<int> faces;
        Edge(int u1, int u2){
            v1 = min(u1,u2);
            v2 = max(u1,u2);
            windows = list<Window>();
        }
};


class EdgeDS {

    public :

    vector<vector<Edge*> > edges;

    EdgeDS(MatrixXd V, MatrixXi F){
        int n = V.rows();
        int f = F.rows();
        edges = vector<vector<Edge*> > (n, vector<Edge*>(n, nullptr));
        for(int i = 0; i < f; i++){
            for(int j = 0; j < 3; j++){
                  int v1 = F(i, j);
                  int v2 = F(i, (j+1)%3);
                  if(edges[v1][v2] != nullptr){
                      (edges[v1][v2]->faces). push_back(f);
                      (edges[v2][v1]->faces). push_back(f);
                  }
                  else{
                    Edge* ed = new  Edge(v1,v2);
                     edges[v1][v2] = ed;
                     edges[v2][v1] = ed;
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

    void djikstra(int s){
            queue<Window> q;



    }


    

};


