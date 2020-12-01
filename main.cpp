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
    cout <<< "hhh";
}
