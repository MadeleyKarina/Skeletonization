#include <iostream>
#include "myDisplay.h"

void myDisplay::push_backMatP(Eigen::MatrixXd& m, Vector3d&& values, std::size_t row)
{
    if(row >= m.rows()) {
        m.conservativeResize(row + 1, Eigen::NoChange);
    }
    m.row(row) = values;
}


Eigen::MatrixXd myDisplay::readEdges2(std::string fileCgaL){
    int cont2;
    ifstream input;
    Eigen::MatrixXd V_box(1,3);
    input.open(fileCgaL);
    while (true) {
        float x1,y1,z1;
        input>>x1>>y1>>z1;
        if( input.eof() ) break;
        push_backMatP(V_box, Vector3d(x1,y1,z1), cont2);
        cont2++;
    }
    input.close();
    return V_box;
}

void myDisplay::showSkel2(igl::viewer::Viewer& viewer, string name){
    Eigen::MatrixXd edges = readEdges2(name);
    for (unsigned i=0; i<edges.rows()-1;i+=2){
        //cout<<"i: "<<i<<endl;
        viewer.data.add_edges
        (
          edges.row(i),
          edges.row(i+1),
          Eigen::RowVector3d(1,0,0)
        );
        viewer.core.show_lines = false;
        viewer.core.show_faces = false;
        viewer.core.show_overlay_depth = true;
        viewer.core.line_width = 1;
    }
}

//void showVertexSkel(string name, igl::viewer::Viewer& viewer){
//    int tmp;
//    Eigen::MatrixXd V_box(1,3);
//    ifstream input;
//    input.open(name);
//    while (true){
//        float x1,y1,z1;
//        input>>x1>>y1>>z1;
//        if(input.eof()) break;
//        push_backMat(V_box, Vector3d(x1,y1,z1), tmp);
//        viewer.data.add_points(V_box.row(0),Eigen::RowVector3d(1,0,1));
//    }
//}
