#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>
#include <igl/exact_geodesic.h>

#include <igl/gaussian_curvature.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>

#include <utility>

#include <igl/loop.h>

#include "EdgeDS.cpp"



using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V;
MatrixXi F;

VectorXd  test_igl(int s){
                Eigen::VectorXi VS,FS,VT,FT;
            // The selected vertex is the source
            VS.resize(1);
            VS << s;
            // All vertices are the targets
            VT.setLinSpaced(V.rows(),0,V.rows()-1);
            Eigen::VectorXd d;
            igl::exact_geodesic(V,F,VS,FS,VT,FT,d);
            return d;
}

vector<int> subdivide(int no, int q){

            MatrixXd V1;
            MatrixXi F1;
          igl::loop(V,F,V1,F1,no);
          EdgeDS goe = EdgeDS(V1,F1,q);
          int w = goe.exactGeo(0);
          vector<int> res;
          res.push_back(V1.rows());
          res.push_back(F1.rows());
            res.push_back(w);

          return res;


}

// ------------ main program ----------------
int main(int argc, char *argv[])
{   
    if(argc < 2) argv[1] = "../data/cube_tri.off";
    int qtype = 1;
    //qtype = 0 for FIFO, 1 for Min distance in extremeties

    if(argc >= 3){
        qtype =  argv[2][0] - '0';
    }

    std::cout << "reading input file: " << argv[1] << std::endl;

	igl::readOFF(argv[1], V, F);
    EdgeDS ds = EdgeDS(V,F,qtype);
    std::cout << "reading DONE " << argv[1] << std::endl;
    int w = ds.exactGeo(0);
    std::cout << "GEO DONE with NUmber of Windows : " << w << std::endl;
        std::cout << " f = " << F.rows() << std::endl;

    /*
    ofstream out("out.out");
    VectorXd dd = test_igl(0);
    for(int i = 0 ; i < V.rows(); i++){
        double d= ds.distance(i);
        out << d << " , " << dd(i) << endl;
    } */

    for(int no = 1; no < 8; no++){
        vector<int> res = subdivide(no,qtype);
        std::cout << "E =" << 3*res[1] /2 << "F : " << res[1] << " , W = " << res[2] << std::endl;
    }
    
    

    igl::opengl::glfw::Viewer viewer; // create the 3d viewer

	viewer.data().set_mesh(V, F);

	viewer.core(0).align_camera_center(V, F);
	viewer.launch(); // run the viewer


}
