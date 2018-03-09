#include "remesh.h"
#include "meshes_path.h"
#include <dirent.h>
#include <sys/stat.h>
using namespace std;

int newRemesh(std::string dataBase,std::string name2, std::string model, int numFaces)
{
    string name_remesh = dataBase + ".off";
    cout << "remeshhhhh "<<name_remesh<<endl;
    Surface_mesh_remesh surface_mesh_remesh;
    std::ifstream is(name_remesh.c_str()) ; is >> surface_mesh_remesh ;
    // This is a stop predicate (defines when the algorithm terminates).
    // In this example, the simplification stops when the number of undirected edges
    // left in the surface mesh drops below the specified number (1000)
    //SMS::Count_stop_predicate<Surface_mesh_remesh> stop(1215);
    SMS::Count_stop_predicate<Surface_mesh_remesh> stop(numFaces);

    // This the actual call to the simplification algorithm.
    // The surface mesh and stop conditions are mandatory arguments.
    // The index maps are needed because the vertices and edges
    // of this surface mesh lack an "id()" field.
    int r = SMS::edge_collapse
              (surface_mesh_remesh
              ,stop
               ,CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index,surface_mesh_remesh))
                                 .halfedge_index_map  (get(CGAL::halfedge_external_index  ,surface_mesh_remesh))
                                 .get_cost (SMS::Edge_length_cost <Surface_mesh_remesh>())
                                 .get_placement(SMS::Midpoint_placement<Surface_mesh_remesh>())
              );

    std::cout << "\nFinished...\n" << r << " edges removed.\n"
              << (surface_mesh_remesh.size_of_halfedges()/2) << " final edges.\n" ;


    std::ofstream os( name2 + "/" + model + ".off" ) ;
    os << surface_mesh_remesh ;
    os.close();

    return 0 ;
  }
void remeshDataBase(int numShapes, int numFaces){
    string nameFile = PATH_REMESHES + to_string(numFaces);
    mkdir(nameFile.c_str(),777);
    string newFile = nameFile + string("/S");
    for (int i = 1; i < numShapes; i++ ){
        string name = newFile + to_string(i);
        mkdir(name.c_str(),777);
        string model = "S" + to_string(i);
        //string pathDB = PATH_MESHES + to_string(i);
        string pathDB = PATH_MESHES_SAMPLE + to_string(i);
        newRemesh(pathDB,name,model,numFaces);
        cout << name << endl;
    }
//    newRemesh("/mnt/data/Code/KeyPoints/project/m178",
//              "/mnt/data/Code/KeyPoints/project/SampleMeshes/db/Remesh/","m178");

}

