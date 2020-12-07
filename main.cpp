#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>

#include <igl/gaussian_curvature.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>

#include "EdgeDS.cpp"



using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V;
MatrixXi F;


// ------------ main program ----------------
int main(int argc, char *argv[])
{   
    if(argc < 2) argv[1] = "../data/cube_tri.off";
    std::cout << "reading input file: " << argv[1] << std::endl;
	igl::readOFF(argv[1], V, F);
    EdgeDS ds = EdgeDS(V,F);
    std::cout << "reading DONE " << argv[1] << std::endl;
            ds.exactGeo(0);
    std::cout << "GEO DONE " << argv[1] << std::endl;
    igl::opengl::glfw::Viewer viewer; // create the 3d viewer

	viewer.data().set_mesh(V, F);

	viewer.core(0).align_camera_center(V, F);
	viewer.launch(); // run the viewer


}
