#ifndef PROCESS_H
#define PROCESS_H

#include "meshes_path.h"
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include <iostream>
//#include <GLFW/glfw3.h>
//#include <GL/glew.h>
#include <nanogui/nanogui.h>
#include <vector>
#include <fstream>
#include <numeric>


using namespace std;
using namespace Eigen;

class Process{
private:
    int numberEigen = 50;
    int numberShapes = 600;
    string currentShape;
    bool typeSkeleton;
    char typeAdjMatrix;
    bool changeNumberVert;
    string currDiag;

    void saveGraph(Eigen::MatrixXd);
public:
    Process();
    Process(string);
    Process(char,bool);
    Process(string, bool, char);
    Process(string, bool, char, bool);

    double	vectors_distance(const std::vector<float>& a, const std::vector<float>& b);
    string getShapeName(string);
    string getResultFileName(string);
    const char * getEigenFileName(string);
    void push_backMat(Eigen::MatrixXd&, Vector3d&&, std::size_t);
    Eigen::MatrixXd readVertForGraph();
    vector<int> getValency(Eigen::MatrixXd );
    Eigen::MatrixXd getLapGraph(Eigen::MatrixXd, vector<int> );
    Eigen::MatrixXd readEdges(string);
    vector <vector<float>> readEigensLap();
    vector <map<int,int>> readJunct();
    map<int,int> getJunction(Eigen::MatrixXd);
    vector<int> countSurfVert(Eigen::MatrixXd, string );
    vector<float> getActLapEigen();
    map<int,int> getCurJunct();
    vector<map<int,int>> compareJunct(map<int,int>, vector<map<int,int>>);
    vector<pair<float, string>> sumJunct(vector<map<int,int>> dataJunct);
    void compActLapEigen(); //save
    map<string, float> compareLargerDesc(vector<vector<float>>, vector<float>);
    map<string, float> compareSmallerDesc(vector<vector<float>>, vector<float>);
    int startDescript(vector<float> );
    map<string, float> compareSmallerPosDesc(vector<vector<float>> , vector<float> );
    void compEucDist();
    void compJunct();
    void printActualEigenLaplace(vector<float> );
    void printDataBaseEigenLaplace();
    void printCurrentJunct();
    void printDataBaseJunct();
    void saveEigen (string, char);
    void compMatrixDissim();
    vector<pair<int, float>> sortMap(map<string, float> eucDist);
    vector<float> distSurfVert(Eigen::MatrixXd, string );


};

#endif
