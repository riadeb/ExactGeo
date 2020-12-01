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
class WindowInterval {
    public :
        int v1,v2; //directed edge (v1,v2)
        list<Window>  windows;
        WindowInterval(int u1, int u2){
            v1 = u1;
            v2 = u2;
            windows = list<Window>();
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
    WindowInterval propagate(int v1, int v2, int v3, double alphal, double alphar, double p1, double p2, bool flip, double d1, double d2, double sigma){
        Vector3d v21 = V.row(v1) - V.row(v2);
        Vector3d v23 = V.row(v3) - V.row(v2);
        double teta = acos(v21.dot(v23)/(v21.norm()*v23.norm()));
        WindowInterval res(v2,v3);
        if(flip){
            res.v1 = v3;
            res.v2 = v2;
        }
        if(M_PI -  alphal - teta <= 0){
            double wl = v23.norm();
            Window w(edges[v2][v3], 0, wl, p1, wl * sin(teta)/sin(alphal) , d1, false );
            res.windows.push_back(w);
            return res;
        }
        double wl = sin(alphal)*p1/(sin(alphal + teta));
        if(wl >= v23.norm()) {
            Window w(edges[v2][v3], 0, wl, p1, wl * sin(teta)/sin(alphal) , d1, false );
            res.windows.push_back(w);
            return res;
        }
        if(alphar <= teta){
            Window windowl(edges[v2][v3], 0, wl, p1, wl * sin(teta)/sin(alphal) , d1, false);
            res.windows.push_back(windowl);
            Window windowr(edges[v2][v3], wl, v23.norm(), d1 + (wl * sin(teta)/sin(alphal)) , d2 + (v23.norm() *sin(teta)/sin(alphar) ) , sigma, false);
            res.windows.push_back(windowr);
            return res;
        } 
        double wr = sin(alphar)*p2/(sin(alphar - teta));
        Window windowl(edges[v2][v3], 0, wl, p1, wl * sin(teta)/sin(alphal) , d1, false);
        res.windows.push_back(windowl);
        Window windowr(edges[v2][v3], wl, wr, d1 + (wl * sin(teta)/sin(alphal)) , d2 + (wr *sin(teta)/sin(alphar) ) , sigma, false);
        res.windows.push_back(windowr);
        Window windowrr(edges[v2][v3], wr, v23.norm() , wr *sin(teta)/sin(alphar), v23.norm() * sin(teta)/ sin(alphar) , d2, false);
        res.windows.push_back(windowrr);
        return res;
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
            Edge* edge = w.edge;
            int v1,v2;
            if(w.tau) { //v1 -> v2 is opposite edge of "edge"
                v1 = edge->v2; v2 = edge->v1;
            }
            else{
                v1 = edge->v1; v2 = edge->v2;
            }
            double dw = abs(w.p1 - w.p2);
            if(adjL[make_pair(v1,v2)].size() == 0) continue;
            pair<int,int> p = adjL[make_pair(v1,v2)][0];
            int v3 = p.second; //p.firsy = v2
            if(w.tau){
                    double dp1 = w.d1;
                    double dp2 = w.d2;
                    double p1 = w.p1;
                    double p2 = w.p2;
                    double alphal = acos((dw*dw + dp1*dp1 - dp2*dp2) / (2*dw*dp1));
                    double alphar = acos((dw*dw + dp2*dp2 - dp1*dp1) / (2*dw*dp2));
                    propagate(v1,v2,v3,alphal,alphar,p1,p2,false,dp1,dp2,w.sigma);

            }
            else{
                    Vector3d v21 = V.row(v1) - V.row(v2);
                   double dp1 = w.d2;
                    double dp2 = w.d1;
                    double p1 = v21.norm() - w.p2;
                    double p2 = v21.norm() - w.p1;
                    double alphal = acos((dw*dw + dp1*dp1 - dp2*dp2) / (2*dw*dp1));
                    double alphar = acos((dw*dw + dp2*dp2 - dp1*dp1) / (2*dw*dp2));
                    propagate(v1,v2,v3,alphal,alphar,p1,p2,false,dp1,dp2,w.sigma); 
            }

            p = adjL[make_pair(v1,v2)][1];
            int v3 = p.first; //p.second = v1
            int temp = v2;
            v2 = v1; v1 = temp;

            if(!w.tau){
                    double dp1 = w.d1;
                    double dp2 = w.d2;
                    double p1 = w.p1;
                    double p2 = w.p2;
                    double alphal = acos((dw*dw + dp1*dp1 - dp2*dp2) / (2*dw*dp1));
                    double alphar = acos((dw*dw + dp2*dp2 - dp1*dp1) / (2*dw*dp2));
                    propagate(v1,v2,v3,alphal,alphar,p1,p2,true,dp1,dp2,w.sigma);

            }
            else{
                    Vector3d v21 = V.row(v1) - V.row(v2);
                   double dp1 = w.d2;
                    double dp2 = w.d1;
                    double p1 = v21.norm() - w.p2;
                    double p2 = v21.norm() - w.p1;
                    double alphal = acos((dw*dw + dp1*dp1 - dp2*dp2) / (2*dw*dp1));
                    double alphar = acos((dw*dw + dp2*dp2 - dp1*dp1) / (2*dw*dp2));
                    propagate(v1,v2,v3,alphal,alphar,p1,p2,true,dp1,dp2,w.sigma); 
            }
            
            



        }

    }


    

};


