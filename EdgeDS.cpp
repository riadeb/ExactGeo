#include <igl/opengl/glfw/Viewer.h>

#include <list>
#include <unordered_map>
#include <utility>
#include <vector>
#include <queue>
#include <__nullptr>
#include <set>



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
        edge = e; p1 = v1; p2 = v2; d1 = dd1; d2 = dd2; sigma = sig; tau = _tau;
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

class WindowComp
{
public:
  size_t operator()(const Window &w1, const Window &w2) const
  {
      double d1min = min(w1.d1, w1.d2);
    double d2min = min(w2.d1, w2.d2);
    if(d1min != d2min) return d1min < d2min;
    if(w1.edge->v1 != w2.edge->v1) return w1.edge->v1 < w2.edge->v1;
    if(w1.edge->v2 != w2.edge->v2) return w1.edge->v2 < w2.edge->v2;
    return w1.p1 < w2.p1;
  }
};

class WindowInterval {
    public :
        int v1,v2; //directed edge (v1,v2)
        list<Window>  windows;
        WindowInterval(int u1, int u2){
            v1 = u1;
            v2 = u2;
        }
        
    
};
class EdgeDS {
double eps = 0.00001;

    public :

    vector<vector<Edge*> > edges;
    unordered_map<pair<int,int>, int , PairHash> adjL;
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
                  adjL[make_pair(v1,v2)] = v3;
                    adjL[make_pair(v2,v3)] = v1;
                     adjL[make_pair(v3,v1)] = v2;
                 
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
    double dist_to_psource_in_new_window(double x, double dv2s, double angle_sv2_v2x){
        return sqrt(x*x + dv2s*dv2s - 2*dv2s*x*cos(angle_sv2_v2x));
    }
    WindowInterval propagate(int v1, int v2, int v3, double alphal, double alphar, double p1, double p2, double d1, double d2, double sigma, bool tau){

        Vector3d v21 = V.row(v1) - V.row(v2);
        Vector3d v23 = V.row(v3) - V.row(v2);
        double teta = acos(v21.dot(v23)/(v21.norm()*v23.norm()));
        WindowInterval res(v2,v3);

        /*
        if(flip){
            res.v1 = v3;
            res.v2 = v2;
        } 
        */
        double sv2 = sqrt(d1*d1 + p1*p1 + 2* d1*p1*cos(alphal));
        double sinbeta = sin(alphar) * d2/sv2;
        if(sinbeta > 1) sinbeta = 1;
        double beta = asin(sinbeta);
        if(M_PI -  alphal - teta <= -eps){
            double wl = v23.norm();
                   Window w(edges[v2][v3], 0, wl, p1, dist_to_psource_in_new_window(wl, sv2,beta + teta) , d1, tau );
              //  res.windows.push_back(w);
          
           
            return res;
        }
        double wl = sin(alphal)*p1/(sin(alphal + teta)); //alphal + teta > 0 because teta > 0
        if(wl-eps >= v23.norm()) {
            Window w(edges[v2][v3], 0, v23.norm(), p1, dist_to_psource_in_new_window(v23.norm(), sv2,beta + teta) , d1, tau );
         //   res.windows.push_back(w);
            return res;
        }
  
        if(alphar-eps <= teta){
            if(wl > eps){
                    Window windowl(edges[v2][v3], 0, wl, p1, wl * sin(teta)/sin(alphal) , d1, tau);
                //    res.windows.push_back(windowl);
            }

            Window windowr(edges[v2][v3], wl, v23.norm(), d1 + (wl * sin(teta)/sin(alphal)) , dist_to_psource_in_new_window(v23.norm(), sv2,beta + teta) , sigma, tau);
            res.windows.push_back(windowr);
            return res;
        } 
        double wr = sin(alphar)*p2/(sin(alphar - teta));
        if(wl > eps){
            Window windowl(edges[v2][v3], 0, wl, p1, wl * sin(teta)/sin(alphal) , d1, tau);
           // res.windows.push_back(windowl);
        }
        if(wr+eps <= v23.norm()){
             Window windowr(edges[v2][v3], wl, wr, d1 + (wl * sin(teta)/sin(alphal)) , d2 + (wr *sin(teta)/sin(alphar) ) , sigma, tau);
            res.windows.push_back(windowr);
            Window windowrr(edges[v2][v3], wr, v23.norm() , wr *sin(teta)/sin(alphar), v23.norm() * sin(teta)/ sin(alphar) , d2, tau);
            //res.windows.push_back(windowrr);
        }
        else{
             Window windowr(edges[v2][v3], wl, v23.norm(), d1 + (wl * sin(teta)/sin(alphal)) , dist_to_psource_in_new_window(v23.norm(), sv2,beta + teta) , sigma, tau);
             res.windows.push_back(windowr);
        }
       
        return res;
    }
    void flip(WindowInterval &wi){ //flips the edge, reversing the windows as well
                WindowInterval w(wi.v2,wi.v1);
                double de = dist(wi.v1, wi.v2);
                for (auto i = wi.windows.rbegin(); i != wi.windows.rend(); ++i) {
                        Window curr_w(i->edge,de - i->p2, de - i->p1,i->d2,i->d1,i->sigma, !i->tau );
                       w.windows.push_back(curr_w);

                }
                wi = w;

        }
    double dist_to_source(Window& w,double x ){ //returns distance of point x in edge to the source using Window w
                if(w.d1 < eps) return w.sigma +  abs(w.p1 - x);

            double dw = abs(w.p2 - w.p1);
                        if(w.d2 < eps) return w.sigma + abs(dw - x);

            double cosalphal =(dw*dw + w.d1*w.d1 - w.d2*w.d2) / (2*dw*w.d1);
            return  w.sigma + sqrt(w.d1*w.d1 + (x - w.p1)*(x - w.p1) - 2* w.d1*(x - w.p1)*cosalphal);
    }
    double dist_to_psource(Window& w,double x ){ //returns distance of point x in edge to the source using Window w
            if(w.d1 < eps) return abs(w.p1 - x);
            double dw = abs(w.p2 - w.p1);
            if(w.d2 < eps) return abs(dw - x);
            double cosalphal =(dw*dw + w.d1*w.d1 - w.d2*w.d2) / (2*dw*w.d1);
            return  sqrt(w.d1*w.d1 + (x - w.p1)*(x - w.p1) - 2* w.d1*(x - w.p1)*cosalphal);
    }

    void intersect(Edge * edge, WindowInterval newWindows, set<Window,WindowComp> & q){

        if(newWindows.v1 > newWindows.v2) flip(newWindows);
        if(edge->windows.empty()){
            edge->windows = newWindows.windows;
            for(auto w : edge->windows) q.insert(w);
            return;
        }
        for(auto w0 : newWindows.windows){
                 if (w0.p1 >= w0.p2) continue;
                 auto it = edge->windows.begin();
                while( it != edge->windows.end()){
                    Window w1 = *it;


                    double p1_inter = max(w0.p1, w1.p1);
                    double p2_inter = min(w0.p2, w1.p2);
                    if(p1_inter >= p2_inter) {
                        it++;
                        continue;
                    }
                    //cut the existing window into 3 parts
                    Window w1_left = w1;
                    w1_left.p2 = p1_inter;
                    w1_left.d2 = dist_to_psource(w1,w1_left.p2);


                    Window w1_right = w1;
                    w1_right.p1 = p2_inter;
                    w1_right.d1 = dist_to_psource(w1, w1_right.p1);


                    Window w1_intersection(w1.edge,p1_inter,p2_inter,w1_left.d2,w1_right.d1,w1.sigma,w1.tau);
                    Window w0_intersection(w0.edge, p1_inter,p2_inter,dist_to_psource(w0,p1_inter), dist_to_psource(w0, p2_inter), w0.sigma,w0.tau);

                    double p1_w0 = dist_to_source(w0, p1_inter);
                    double p2_w0 = dist_to_source(w0, p2_inter);
                    double p1_w1 = dist_to_source(w1, p1_inter);
                    double p2_w1 = dist_to_source(w1, p2_inter);
                    
                    

                    if(p1_w0 + eps >= p1_w1 && p2_w0 + eps  >= p2_w1){
                        //w1_intersection is smaller, keep it, do nothing
                        it++;
                        
                    }
                    else if (p1_w0 - eps <= p1_w1 && p2_w0- eps <= p2_w1){
                        //w0_intersection is smaller
                      
                        list<Window> newWins;
                        if(w1_left.p1 + eps < w1_left.p2) newWins.push_back(w1_left);
                        newWins.push_back(w0_intersection);
                        if(w1_right.p1 +eps < w1_right.p2) newWins.push_back(w1_right);

                        it = edge->windows.erase(it);
                        edge->windows.insert(it, newWins.begin(), newWins.end());

                         q.insert(w0_intersection);

                         auto todelete = q.find(w1);
                         if(todelete != q.end()){
                             q.erase(w1);
                             if(w1_left.p1 +eps < w1_left.p2) q.insert(w1_left);
                                if(w1_right.p1 +eps < w1_right.p2) q.insert(w1_right);
                         }
                        
                    }
                    else if(p1_w0 >= p1_w1 && p2_w0 <= p2_w1){
                            //w1 is the "left" window
                            auto sol = solveEq(w1_intersection,w0_intersection);
                            double q1 = sol.first; double q2 = sol.second;
                            if(q1 == INFINITY || q1 != q1) {
                                it++;
                                continue;
                            }
                            if((abs(q1 - q2) < eps &&  belongInterval(q1,p1_inter,p2_inter)) || (belongInterval(q1,p1_inter,p2_inter) && !belongInterval(q2,p1_inter,p2_inter)) ||  (belongInterval(q2,p1_inter,p2_inter) && !belongInterval(q1,p1_inter,p2_inter))){
                               //One solution to the equation in the intersection interval
                                if( belongInterval(q2,p1_inter,p2_inter) && !belongInterval(q1,p1_inter,p2_inter) ) q1 = q2;
                                    w1_intersection.p2 = q1;
                                    w1_intersection.d2 = dist_to_psource(w1,q1);

                                    w0_intersection.p1 = q1;
                                    w0_intersection.d1 = dist_to_psource(w0,q1);

                                    list<Window> newWins;
                                    w1_left.p2 = w1_intersection.p2; //stitches the window 1
                                    w1_left.d2 = w1_intersection.d2;
                                    if(w1_left.p1 + eps < w1_left.p2) newWins.push_back(w1_left);
                                    if(w0_intersection.p1 +eps < w0_intersection.p2) newWins.push_back(w0_intersection);
                                    if(w1_right.p1 + eps < w1_right.p2) newWins.push_back(w1_right);

                                    it = edge->windows.erase(it);
                                    edge->windows.insert(it, newWins.begin(), newWins.end());
                                    if(w0_intersection.p1 +eps < w0_intersection.p2)  q.insert(w0_intersection);

                                    auto todelete = q.find(w1);
                                    if(todelete != q.end()){
                                        q.erase(w1);
                                        if(w1_left.p1 + eps < w1_left.p2) q.insert(w1_left);
                                            if(w1_right.p1 + eps< w1_right.p2) q.insert(w1_right);
                                    }


                                
                            }
                            else if(belongInterval(q1,p1_inter,p2_inter) && belongInterval(q2,p1_inter,p2_inter)) //two solutions in interval
                            {
                                    double temp = q1;
                                    q1 = min(q1,q2);
                                    q2 = max(q2,temp);

                                    w0_intersection.p1 = q1;
                                    w0_intersection.d1 = dist_to_psource(w0,q1);
                                    w0_intersection.p2 = q2;
                                    w0_intersection.d2 = dist_to_psource(w0,q2);

                                    Window w1_intersection1 = w1_intersection;
                                    w1_intersection1.p2 = q1;
                                    w1_intersection1.d2 = dist_to_psource(w1,q1);

                                    Window w1_intersection2 = w1_intersection;
                                    w1_intersection1.p1 = q2;
                                    w1_intersection1.d1 = dist_to_psource(w1,q2);

                                    w1_left.p2 = w1_intersection1.p2; //stich windows back togther
                                    w1_left.d2 = w1_intersection1.d2;

                                    w1_right.p1 = w1_intersection2.p1;
                                    w1_right.d1 = w1_intersection2.d1;

                                    list<Window> newWins;
                                    if(w1_left.p1 +eps < w1_left.p2) newWins.push_back(w1_left);
                                   // if(w1_intersection1.p1 + eps< w1_intersection1.p2) newWins.push_back(w1_intersection1);
                                    if(w0_intersection.p1 +eps < w0_intersection.p2)newWins.push_back(w0_intersection);
                                    //if(w1_intersection2.p1 +eps < w1_intersection2.p2)newWins.push_back(w1_intersection2);
                                    if(w1_right.p1 +eps < w1_right.p2)newWins.push_back(w1_right);

                                    it = edge->windows.erase(it);
                                    edge->windows.insert(it, newWins.begin(), newWins.end());

                                    if(w0_intersection.p1 +eps < w0_intersection.p2)  q.insert(w0_intersection);

                                    auto todelete = q.find(w1);
                                    if(todelete != q.end()){
                                        q.erase(w1);
                                        if(w1_left.p1 +eps < w1_left.p2) q.insert(w1_left);
                                            if(w1_right.p1  +eps < w1_right.p2) q.insert(w1_right);
                                    }

                            }
                    }
                    else{
                        //w0 is the "left" window
                            auto sol = solveEq(w0_intersection,w1_intersection);
                            double q1 = sol.first; double q2 = sol.second;
                            if(q1 == INFINITY || q1 != q1) {
                                it++;
                                continue;
                            }
                            if((abs(q1 - q2) < eps &&  belongInterval(q1,p1_inter,p2_inter)) || (belongInterval(q1,p1_inter,p2_inter) && !belongInterval(q2,p1_inter,p2_inter)) ||  (belongInterval(q2,p1_inter,p2_inter) && !belongInterval(q1,p1_inter,p2_inter))){
                               //One solution to the equation in the intersection interval
                                if( belongInterval(q2,p1_inter,p2_inter) && !belongInterval(q1,p1_inter,p2_inter) ) q1 = q2;
                                    w0_intersection.p2 = q1;
                                    w0_intersection.d2 = dist_to_psource(w1,q1);

                                    w1_intersection.p1 = q1;
                                    w1_intersection.d1 = dist_to_psource(w0,q1);

                                     w1_right.p1 = w1_intersection.p1;
                                    w1_right.d1 = w1_intersection.d1;

                                    list<Window> newWins;
                                    

                                    if(w1_left.p1 +eps < w1_left.p2) newWins.push_back(w1_left);
                                    if(w0_intersection.p1 +eps < w0_intersection.p2) newWins.push_back(w0_intersection);
                                    if(w1_intersection.p1 +eps < w1_intersection.p2) newWins.push_back(w1_intersection);
                                    if(w1_right.p1 +eps < w1_right.p2) newWins.push_back(w1_right);

                                    it = edge->windows.erase(it);
                                    edge->windows.insert(it, newWins.begin(), newWins.end());
                                    if(w0_intersection.p1 +eps < w0_intersection.p2)  q.insert(w0_intersection);

                                    auto todelete = q.find(w1);
                                    if(todelete != q.end()){
                                        q.erase(w1);
                                        if(w1_left.p1 +eps < w1_left.p2) q.insert(w1_left);
                                            if(w1_right.p1 +eps < w1_right.p2) q.insert(w1_right);
                                    }

                                
                            }
                            else if(belongInterval(q1,p1_inter,p2_inter) && belongInterval(q2,p1_inter,p2_inter)) //two solutions in interval
                            {
                                    double temp = q1;
                                    q1 = min(q1,q2);
                                    q2 = max(q2,temp);

                                    w1_intersection.p1 = q1;
                                    w1_intersection.d1 = dist_to_psource(w1,q1);
                                    w1_intersection.p2 = q2;
                                    w1_intersection.d2 = dist_to_psource(w1,q2);

                                    Window w0_intersection1 = w0_intersection;
                                    w0_intersection1.p2 = q1;
                                    w0_intersection1.d2 = dist_to_psource(w0,q1);

                                    Window w0_intersection2 = w0_intersection;
                                    w0_intersection1.p1 = q2;
                                    w0_intersection1.d1 = dist_to_psource(w0,q2);

                                    list<Window> newWins;
                            

                                    if(w1_left.p1 +eps < w1_left.p2) newWins.push_back(w1_left);
                                    if(w0_intersection1.p1 +eps < w0_intersection1.p2) newWins.push_back(w0_intersection1);
                                    if(w1_intersection.p1 +eps  < w1_intersection.p2)newWins.push_back(w1_intersection);
                                    if(w0_intersection2.p1 +eps < w0_intersection2.p2)newWins.push_back(w0_intersection2);
                                    if(w1_right.p1 +eps < w1_right.p2)newWins.push_back(w1_right);

                                    it = edge->windows.erase(it);
                                    edge->windows.insert(it, newWins.begin(), newWins.end());

                                    if(w0_intersection1.p1 +eps < w0_intersection1.p2) q.insert(w0_intersection1);
                                    if(w0_intersection2.p1 + eps < w0_intersection2.p2) q.insert(w0_intersection2);

                                     auto todelete = q.find(w1);
                                        if(todelete != q.end()){
                                            q.erase(w1);
                                            if(w1_left.p1 + eps < w1_left.p2) q.insert(w1_left);
                                                if(w1_right.p1 +eps  < w1_right.p2) q.insert(w1_right);
                                                if(w1_intersection.p1 +eps  < w1_intersection.p2) q.insert(w1_intersection);
                                        }

                            }
                    }

                }
        }

    }
    
    pair<double,double> solveEq(Window w0, Window w1) {
        double dw = w0.p2 - w0.p1; //both windows span same interval
        double cos_alphal0 = (dw*dw + w0.d1*w0.d1 - w0.d2*w0.d2)/(2*dw*w0.d1);
        double cos_alphal1 = (dw*dw + w1.d1*w1.d1 - w1.d2*w1.d2)/(2*dw*w1.d1);

        double sin_alphal0 = sin(acos(cos_alphal0));
        double sin_alphal1 = sin(acos(cos_alphal1));

        double s0x = w0.d1*cos_alphal0 + w0.p1;
        double s1x = w1.d1*cos_alphal1 + w0.p1;

        double s0y = w0.d1*sin_alphal0;
        double s1y = w1.d1*sin_alphal1;

        double alpha = w1.d1*cos_alphal1 - w0.d1*cos_alphal0;
        double beta = w1.sigma - w0.sigma;
        double gamma =  (s0x)*s0x + s0y*s0y - s1x*s1x - s1y*s1y - beta*beta;

        double A = alpha*alpha - beta*beta;
        double B = gamma*alpha + 2* s1x*beta*beta;
        double C = 0.25 * gamma*gamma - beta*beta*(s1x*s1x + s1y*s1y);

        if(abs(A) < eps && abs(B) > eps) return make_pair(-C/B, -C/B);
        if(abs(A) < eps && abs(B) < eps) return make_pair(INFINITY, INFINITY);

        double det = B*B - 4*A*C;
        if(abs(det)< eps){
            return make_pair(-B/(2*A), -B/(2*A));
        }
        return make_pair((-B - sqrt(det))/ (2*A),   (-B + sqrt(det))/ (2*A));




    }
    bool belongInterval(double x, double a, double b){
        return (x  >= a + eps) && (x + eps <= b);
    }
    
    void printEdge(Edge* e){
        cout << e->v1 << " --- " << e->v2 << endl;
        for(auto w : e->windows){
                cout << w.p1 << " ( " << w.d1 << " ) --- " <<w.p2 << " ( " << w.d2 << " ) sgima = " << w.sigma << " , tau = " << w.tau << endl; 
        }
    }
    void exactGeo(int s){

            set<Window, WindowComp> q;

            int n = edges.size();
            for(int i = 0; i < n; i++){ //initialize queue
                if(edges[s][i] != nullptr){
                    if(s < i){
                         Window w = Window(edges[s][i], 0, dist(s,i), 0, dist(s,i),0,true);
                         edges[s][i]->windows.push_back(w);
                        //q.push(w);
                    }
                    else{
                         Window w = Window(edges[s][i], 0, dist(s,i), dist(s,i),0,0,false);
                         edges[s][i]->windows.push_back(w);
                         //q.push(w);
                    }
                    auto e = adjL[make_pair(s,i)];
                    Edge* edge =  edges[e][i];
                    if(edge->windows.size() == 0){
                        bool tau = i < e; //direction of window
                        Window w = Window(edge,0, dist(edge->v1,edge->v2), dist(s,edge->v1), dist(s, edge->v2), 0,tau);
                        q.insert(w);
                        
                        edge->windows.push_back(w);
                    }
                        
                    
                }
            }

        while(! q.empty() ){
            Window w = *q.begin();
            q.erase(q.begin());
            Edge* edge = w.edge;
            int v1,v2;
            v1 = edge->v1;
            v2 = edge->v2;
            double dw = abs(w.p1 - w.p2);
            int v3;
            if(w.tau){
                    if(adjL.find(make_pair(v2,v1)) == adjL.end()) continue;
                    v3 = adjL[make_pair(v2,v1)];
            }
            else{
                    if(adjL.find(make_pair(v1,v2)) == adjL.end()) continue;
                    v3 = adjL[make_pair(v1,v2)];
            }
            

            WindowInterval newWindows(v2,v3);
           
            Vector3d v21 = V.row(v1) - V.row(v2);
            double dp1 = w.d2;
            double dp2 = w.d1;
            double p1 = v21.norm() - w.p2;
            double p2 = v21.norm() - w.p1;
            double alphal = acos((dw*dw + dp1*dp1 - dp2*dp2) / (2*dw*dp1));
            double alphar = acos((dw*dw + dp2*dp2 - dp1*dp1) / (2*dw*dp2));
            if(isnan(alphar) || isnan(alphal) ) continue;
            newWindows = propagate(v1,v2,v3,alphal,alphar,p1,p2,dp1,dp2,w.sigma,! w.tau); 


            Edge* edge_to_propagate = edges[v2][v3];
            cout << "edge before " << endl;
            printEdge(edge_to_propagate);
            intersect(edges[v2][v3],newWindows,q);
            cout << "edge after" << endl;
            printEdge(edge_to_propagate);


            int temp = v2;
            v2 = v1; v1 = temp;
             dp1 = w.d1;
             dp2 = w.d2;
             p1 = w.p1;
            p2 = w.p2;
             alphal = acos((dw*dw + dp1*dp1 - dp2*dp2) / (2*dw*dp1));
             alphar = acos((dw*dw + dp2*dp2 - dp1*dp1) / (2*dw*dp2));
            newWindows = propagate(v1,v2,v3,alphal,alphar,p1,p2,dp1,dp2,w.sigma, w.tau);

           
            edge_to_propagate = edges[v2][v3];

            cout << "edge before " << endl;
            printEdge(edge_to_propagate);
            intersect(edges[v2][v3],newWindows,q);
            cout << "edge after" << endl;
            printEdge(edge_to_propagate);




        }

    }


    

};


