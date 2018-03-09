#ifndef MYDISPLAY_H
#define MYDISPLAY_H

#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include <iostream>
//#include <GLFW/glfw3.h>
//#include <GL/glew.h>
#include <nanogui/nanogui.h>
#include "meshes_path.h"
#include <vector>
#include <stdio.h>
#include <dirent.h>
#include <sys/stat.h>
#include <cstdlib>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include "remesh.h"
#include "mySkeleton.h"
using namespace Eigen;
using namespace std;

class myDisplay{
    public:
        void push_backMatP(Eigen::MatrixXd& m, Vector3d&& values, std::size_t row);
        Eigen::MatrixXd readEdges2(std::string);
        void showSkel2(igl::viewer::Viewer& viewer, string name);

};


#endif
