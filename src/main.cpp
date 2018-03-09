#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include <iostream>
//#include <GLFW/glfw3.h>
//#include <GL/glew.h>
#include <nanogui/nanogui.h>
#include "funcgui.h"
#include "meshes_path.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/Surface_mesh.h>
#include <boost/foreach.hpp>
#include <vector>
#include <stdio.h>

#include <cstdlib>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include "Mesh.h"
#include "Clock.h"
#include "util.h"
#include "Properties.h"
#include "HarrisDetector.h"
#include "remesh.h"
#include "mySkeleton.h"
#include "myDisplay.h"
#include "Process.h"



using namespace std;
using namespace Eigen;


string path = "/Users/madeley/Documents/examples/project/SampleMeshes/db/newDB/T";
string pathResult = "/mnt/data/Code/KeyPoints/results_query/";
string pathRemesh = "/mnt/data/Code/KeyPoints/project/SampleMeshes/db/Remesh/400/T";
//string path = "/mnt/data/Code/KeyPoints/project/SampleMeshes/db/mine2000/"; //dir
int numShapes = 600;
int numEigen = 30;
string nameSimpSkel = "skel2.cgal";
string nameVertSimpSkel = "vertex.cgal";
string nameMCFSkel = "skel-sm.cgal";
string nameVertMCFSkel = "vertex-sm.cgal";
string nameCorrVert = "corr_skel2.cgal";

string chosenFile = "/m1500.off";

string DataProof = "/mnt/data/Code/KeyPoints/project/SampleMeshes/db/Remesh/400/T";
string pathDBSkel = "/mnt/data/Code/KeyPoints/project/SampleMeshes/db/newDB_Skel/T";

//only needed for the display of the skeleton as maximal polylines

template <typename DynamicEigenMatrix>
void push_backMat(DynamicEigenMatrix& m, Vector3d&& values, std::size_t row)
{
    if(row >= m.rows()) {
        m.conservativeResize(row + 1, Eigen::NoChange);
    }
    m.row(row) = values;
}

void showSkel(igl::viewer::Viewer& viewer, Eigen::MatrixXd edges){
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
void showVertexSkel(string name, igl::viewer::Viewer& viewer){
    int tmp;
    Eigen::MatrixXd V_box(1,3);
    ifstream input;
    input.open(name);
    while (true){
        float x1,y1,z1;
        input>>x1>>y1>>z1;
        if(input.eof()) break;
        push_backMat(V_box, Vector3d(x1,y1,z1), tmp);
        viewer.data.add_points(V_box.row(0),Eigen::RowVector3d(1,0,1));
    }
}



// Generating the file with the key points
string getOutputPath(string filename){
  string outputPath(filename);
  size_t pos;
  string directory = "";

  if((pos = filename.find_last_of('/')) != string::npos){
    directory = filename.substr(0, pos + 1);
    outputPath = filename.substr(pos + 1);
  }

  pos = outputPath.find_last_of('.');

  outputPath = outputPath.substr(0, pos);

  outputPath = directory + outputPath + ".int";

  return outputPath;
}



/*******************************************************************************/



int main(int argc, char *argv[])
{
    string offFilename = "0002.null.0.off";//argv[1];
    string outfilename = getOutputPath(offFilename);
    ifstream inp;
    inp.open(outfilename.c_str(), ifstream::in);
    inp.close();
    if(inp.fail()){
        inp.clear(ios::failbit);
        int i = 0;
        Clock r;
        r.tick();
        Mesh mesh("0002.null.0.off");
        //Mesh mesh(argv[1]);
        r.tick();
        cout<<"->Loaded";
        cout<<"->Loading time:"<<r.getTime()<<"s."<<endl;
        Properties prop;
//		if(argc == 3)
//			prop.load(optFilename);
        HarrisDetector hd(&mesh, &prop);
        vector<int> interestPoints;
        interestPoints.clear();
        Vertex* vertices = mesh.getVertices();
        Clock r1;
        r1.tick();
        cout<<"->Interest points"<<endl;
        hd.detectInterestPoints(interestPoints);
        r1.tick();
        cout<<"->Calculation time:"<<r1.getTime()<<"s."<<endl;
        cout<<"->Saving"<<endl;
        ofstream out(outfilename.c_str());
        //for(int i = 0; i < interestPoints.size(); i++)
        //	cout << interestPoints[i] << endl;
        //cout << "Interest points:" <<interestPoints.size() << endl;
        out<<interestPoints.size()<<endl;
        for(register int i = 0; i < interestPoints.size(); i++) {
            out<<interestPoints[i]<<endl;
        }
        out.close();
        interestPoints.clear();
        cout<<"File saved."<<endl;
    }else{
        cout << "Interest points already exist - Skipped" << endl;
    }


/***********************************Skeleton*******************************************/
//    mySkeleton preCompSkel;
//          for (int i = 0; i < numShapes; i++){
//              string fileSaveSkel = PATH_MESHES + to_string(i) + "/T"+ to_string(i) ;
//              preCompSkel.saveData(fileSaveSkel);
//              cout << i<<endl;
//           }
/***********************************EigenConn*********************************************/
       // Process pDabaBase;
       // for (int i = 0; i < numShapes; i++){
       //   string fileSaveEigen = DataProof + to_string(i)  + "/T" + to_string(i) ;
       //   pDabaBase.saveEigen(fileSaveEigen, 'C');
       //   cout <<fileSaveEigen <<endl;
       // }
/////***********************************EigenLap******************************************/

       // Process pDabaBase;
       // for (int i = 0; i < numShapes; i++){
       //   string fileSaveEigen = path + to_string(i)  + "/T" + to_string(i) ;
       //   pDabaBase.saveEigen(fileSaveEigen, 'L');
       //   cout <<fileSaveEigen <<endl;
       // }

/////***********************************EigenLap******************************************/

      // Process pDabaBase;
      // //for (int i = 234; i < 250; i++){//_bin
      // for (int i = 0; i < numShapes; i++){//_bin250
      //   try {
      //     string fileSaveEigen = path + to_string(i)  + "/T" + to_string(i) ;
      //     pDabaBase.saveEigen(fileSaveEigen, 'A');
      //     cout <<fileSaveEigen <<endl;
      //   }
      //    catch (int e) {
      //       cout <<"Error at i " <<i<<endl;
      //
      //   }
      //
      // }


    /////***********************************EigenCorr******************************************/
//            Process pDB;
//            for (int i = 0; i< numShapes; i++){
//              string fileSaveEigen = PATH_MESHES + to_string(i)  + "/T" + to_string(i) ;
//              pDB.saveEigen(fileSaveEigen, 'N');
//              cout <<fileSaveEigen <<endl;
//            }
    /////***********************************EigenDist******************************************/

//            Process pDabaBase;
//            for (int i = 0; i < numShapes; i++){
//              string fileSaveEigen = PATH_MESHES + to_string(i)  + "/T" + to_string(i) ;
//              pDabaBase.saveEigen(fileSaveEigen, 'N');
//              cout <<fileSaveEigen <<endl;
//            }
/////***********************************EigenWeight******************************************/

//            Process pDabaBase;
//            for (int i = 0; i < numShapes; i++){
//              string fileSaveEigen = PATH_MESHES + to_string(i)  + "/T" + to_string(i) ;
//              pDabaBase.saveEigen(fileSaveEigen, 'A');
//              cout <<fileSaveEigen <<endl;
//            }
/////***********************************EigenWeight_2******************************************/

           // Process pDabaBase;
           // for (int i = 0; i < numShapes; i++){//_bin
           //   string fileSaveEigen = PATH_MESHES + to_string(i)  + "/T" + to_string(i) ;
           //   pDabaBase.saveEigen(fileSaveEigen, 'W');
           //   cout <<fileSaveEigen <<endl;
           // }

    /////***********************************Junc******************************************/

////      for (int i = 0; i < numShapes; i++){
////          string fileSaveJunct = path + to_string(i) + "/T" + to_string(i);
////          saveJunct(fileSaveJunct);
////          cout <<fileSaveJunct <<endl;
////      }
//        for (int i = 0; i < numShapes; i++){
//            string fileSaveJunct = DataProof + to_string(i) + "/T" + to_string(i);
//            saveJunct(fileSaveJunct);
//            cout <<fileSaveJunct <<endl;
//        }
/************************************************************************************/

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    bool skelProperties = true;
    string chosenFile = "/m1500.off";
    igl::readOFF(SAMPLE_MESHES_PATH+chosenFile, V, F);
    // Init the viewer
    igl::viewer::Viewer viewer;
    viewer.data.set_mesh(V, F);
    bool boolVariable = false;

    //objeto to process functions for manipulpara procesar funciones de manipulacion
    Process p;
    // Extend viewer menu
    viewer.callback_init = [&](igl::viewer::Viewer &viewer)
    {
        viewer.ngui->addButton("Load", [&viewer, &V, &F,&chosenFile](){
            chosenFile = nanogui::file_dialog({ { "off", "OFF Triangle Mesh"} }, false);
            cout << chosenFile << endl;
            viewer.data.clear();
            igl::readOFF(chosenFile, V, F);
            viewer.data.set_mesh(V,F);
        });
        viewer.ngui->addGroup("Skeletonization");
        viewer.ngui->addButton("Simple Skeleton",[&viewer,&V, &F,&chosenFile, &skelProperties,&p](){
            viewer.data.clear();
            viewer.data.set_mesh(V, F);
            skelProperties = true;
            mySkeleton skel(chosenFile, skelProperties);
            skel.saveEdges(nameSimpSkel);
            skel.saveVertex(nameVertSimpSkel);
            skel.saveSurfPoints(nameCorrVert);
           // p.countSurfVert(nameVertSimpSkel,"correspondance-poly.cgal");
            string fileSkeleton = nameSimpSkel;
            Eigen::MatrixXd matrixSkel = p.readEdges(fileSkeleton);
            showSkel(viewer,matrixSkel);
        });
        viewer.ngui->addButton("MCF Skeleton",[&viewer,&V, &F,&chosenFile,&skelProperties,&p](){
            cout << "chosen new skeleton: " << chosenFile <<endl;
            viewer.data.clear();
            viewer.data.set_mesh(V, F);
            skelProperties = false;
            mySkeleton skel(chosenFile,skelProperties);
            skel.saveEdges(nameMCFSkel);
            skel.saveVertex(nameVertMCFSkel);
            string fileSkeleton = nameMCFSkel;
            Eigen::MatrixXd matrixSkel = p.readEdges(fileSkeleton);
            showSkel(viewer,matrixSkel);
        });
        viewer.ngui->addVariable<bool>("Show Vertex Skeleton", [&viewer,&V, &F, &skelProperties,&p](bool val){
            if(val){
                if(skelProperties){
                    showVertexSkel(nameSimpSkel, viewer);
                }
                else{
                    showVertexSkel(nameMCFSkel, viewer);
                }
            }
            else{
                viewer.data.clear();
                viewer.data.set_mesh(V, F);
                if(skelProperties){
                    string  fileSkeleton = nameSimpSkel;
                    Eigen::MatrixXd matrixSkel = p.readEdges(fileSkeleton);
                    showSkel(viewer,matrixSkel);
                }
                else{
                    string  fileSkeleton = nameMCFSkel;
                    Eigen::MatrixXd matrixSkel = p.readEdges(fileSkeleton);
                    showSkel(viewer,matrixSkel);
                }
            }
            },[&]() {
            return boolVariable;
        });
        viewer.ngui->addButton("Junction", [&viewer, &V, &F,&chosenFile,&skelProperties](){
            cout <<"Result for Junction"<<endl;
           Process process(chosenFile);
           process.compJunct();
            map<int,int> descJunct = process.getCurJunct();
                    process.printCurrentJunct();
//            vector<map<int,int>> dataJunct= process.readJunct();
//            process.printDataBaseJunct();
         });
        viewer.ngui->addWindow(Eigen::Vector2i(220,10),"Evaluations");
        viewer.ngui->addGroup("Compare - 3000 vertices");
        viewer.ngui->addButton("Compare Eigen - Weight", [&viewer, &V, &F,&chosenFile,&skelProperties](){
            if(skelProperties == true){
                cout <<"Result for Simple Skeleton"<<endl;
            }
            else{
                cout <<"Result for MCF Skeleton"<<endl;
            }
            Process process(chosenFile,skelProperties,'A');
            process.compEucDist();
         });
        viewer.ngui->addButton("Compare Eigen - Lap", [&viewer, &V, &F,&chosenFile,&skelProperties](){
            if(skelProperties == true){
                cout <<"Result for Simple Skeleton"<<endl;
            }
            else{
                cout <<"Result for MCF Skeleton"<<endl;
            }
            Process process(chosenFile,skelProperties,'L');
            process.compEucDist();
         });
        viewer.ngui->addButton("Compare Eigen - Conn", [&viewer, &V, &F,&chosenFile,&skelProperties](){
            if(skelProperties == true){
                cout <<"Result for Simple Skeleton"<<endl;
            }
            else{
                cout <<"Result for MCF Skeleton"<<endl;
            }
            Process process(chosenFile,skelProperties, 'C');
            process.compEucDist();
         });
        viewer.ngui->addButton("Compare Eigen - Corr", [&viewer, &V, &F,&chosenFile,&skelProperties](){
            if(skelProperties == true){
                cout <<"Result for Simple Skeleton"<<endl;
            }
            else{
                cout <<"Result for MCF Skeleton"<<endl;
            }
            Process process(chosenFile,skelProperties,'N');
            process.compEucDist();
         });
        viewer.ngui->addGroup("Compare - 400 vertices");
        viewer.ngui->addButton("Compare Eigen - Lap", [&viewer, &V, &F,&chosenFile,&skelProperties](){
            if(skelProperties == true){
                cout <<"Result for Simple Skeleton"<<endl;
            }
            else{
                cout <<"Result for MCF Skeleton"<<endl;
            }
            Process process(chosenFile,skelProperties,'L', true); // 400 VERTEX
            process.compEucDist();
         });
        viewer.ngui->addButton("Compare Eigen - Conn", [&viewer, &V, &F,&chosenFile,&skelProperties](){
            if(skelProperties == true){
                cout <<"Result for Simple Skeleton"<<endl;
            }
            else{
                cout <<"Result for MCF Skeleton"<<endl;
            }
            Process process(chosenFile,skelProperties, 'C', true);
            process.compEucDist();
         });
        viewer.ngui->addButton("Process Matrix: 3000 Lap", [&viewer, &V, &F,&chosenFile,&skelProperties](){
//            Process process('L',false);
//            process.compMatrixDissim();
//            Process process2('C',false);
//            process2.compMatrixDissim();
//            Process process2('N',false);
            // Process process3('A',false);
            // process3.compMatrixDissim();
            Process process('W',false);
            process.compMatrixDissim();

         });
      viewer.screen->performLayout();
        return false;
    };
    viewer.launch();
}
