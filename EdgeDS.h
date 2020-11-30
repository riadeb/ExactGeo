#include <igl/opengl/glfw/Viewer.h>

#include <list>
#include <unordered_map>
#include <utility>
#include <vector>


using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;



class Window {
    public :
    Edge* edge;
    int p1,p2;
    double d1,d2,sigma;
    bool to;
    Window(Edge* e, int v1, int v2, double dd1, double dd2, double sig, bool t);

};
class Edge {
    public :
        int v1,v2; //v1 smaller that v2
        list<Window>  windows;
        Edge(int u1, int u2);
};


class EdgeDS {

    public :

    vector<vector<Edge*> > edges;

    EdgeDS(MatrixXd V, MatrixXi F){
        int n = V.rows();
        int f = F.rows();
        edges = vector<vector<Edge*> > (n, vector<Edge*>(n, 0));
        for(int i = 0; i < f; i++){
            for(int j = 0; j < 3; j++){
                  pair<int,int> p;
                  int v1 = F(i, j);
                  int v2 = F(i, (j+1)%3);
                  Edge* ed = new  Edge(v1,v2);
                  edges[v1][v2] = ed;
                  edges[v2][v1] = ed;
            }

        }

    }
    ~EdgeDS(){
        int n = edges.size();
        for(int i = 0; i < n-1; i++){
            for(int j = i+1; j < n ;j++) delete edges[i][j];
        }
    }


    

};


