#include <iostream>
#include "mySkeleton.h"

using namespace std;

struct Display_polylines{
  const Skeleton& skeleton;
  std::ofstream& out;
  int polyline_size;
  std::stringstream sstr;
  Display_polylines(const Skeleton& skeleton, std::ofstream& out)
    :skeleton(skeleton), out(out)
  {}

  void start_new_polyline(){
    polyline_size=0;
    sstr.str("");
    sstr.clear();
  }

  void add_node(Skeleton_vertex v){
    ++polyline_size;
    sstr << " " << skeleton[v].point << '\n';
    //std::cout << "point: " << skeleton[v].point<< "\n";
  }

  void end_polyline()
  {
    //out << polyline_size << "\n" << sstr.str() << "\n";
      out << sstr.str() << 1000 << "\n" ;
     // std::cout << polyline_size << "\n";
  }
};

mySkeleton::mySkeleton(){
    cout << "For precomputation stage " <<endl;
}
mySkeleton::mySkeleton(string fileModel, bool properties){
    std::ifstream input(fileModel);
    flagProperties = properties;
    if(flagProperties == true){
        input >> tmesh;
        CGAL::extract_mean_curvature_flow_skeleton(tmesh, skeleton);
        std::cout << fileModel <<"Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
    }
    else{
        input >> tmesh_ob;
        if (!CGAL::is_triangle_mesh(tmesh_ob))
        {
          std::cout << "Input geometry is not triangulated." << std::endl;
        }
        Skeletonization_ob mcs(tmesh_ob);
        mcs.set_max_iterations(700);
        mcs.set_max_triangle_angle(110);
        // 1. Contract the mesh by mean curvature flow.
        mcs.contract_geometry();
        // 2. Collapse short edges and split bad triangles.
        mcs.collapse_edges();
        mcs.split_faces();
        // 3. Fix degenerate vertices.
        mcs.detect_degeneracies();
        // Perform the above three steps in one iteration.
        mcs.contract();
        // Iteratively apply step 1 to 3 until convergence.
        mcs.contract_until_convergence();
        // Convert the contracted mesh into a curve skeleton and
        // get the correspondent surface points
        mcs.convert_to_skeleton(skeleton_ob);
    }
}
void mySkeleton::savePolylines(string name){
    std::ofstream outputPol;
    outputPol.open(name,ios::trunc);// name = "skel.cgal"
    Display_polylines display(skeleton,outputPol);
    CGAL::split_graph_into_polylines(skeleton, display);
    outputPol.close();
}
void mySkeleton::saveEdges(string name){
  //  std::cout << chosenFile <<"Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
    std::ofstream outputEdges(name);// name = "skel2.cgal"
    if(flagProperties == true){
        BOOST_FOREACH(Skeleton_edge e, edges(skeleton))
        {
            Point& s = skeleton[source(e, skeleton)].point;
            Point& t = skeleton[target(e, skeleton)].point;
            outputEdges << s << "\n" << t << "\n";
        }
    }
    else{
        BOOST_FOREACH(Skeleton_edge_ob e, edges(skeleton_ob))
        {
          const Point& s = skeleton_ob[source(e, skeleton_ob)].point;
          const Point& t = skeleton_ob[target(e, skeleton_ob)].point;
          outputEdges << s << "\n" << t << "\n";
        }
    }
    outputEdges.close();
}
void mySkeleton::saveVertex(string name){
    std::ofstream outputVertex(name);// name = "vertex.cgal"
    if(flagProperties == true){
        BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton))
     // BOOST_FOREACH(vertex_descriptor vd, skeleton[v].vertices)
            outputVertex << skeleton[v].point << "\n";
    }
    else{
        BOOST_FOREACH(Skeleton_vertex_ob v, vertices(skeleton_ob))
         // BOOST_FOREACH(vertex_descriptor vd, skeleton[v].vertices)
            outputVertex << skeleton_ob[v].point << "\n";
    }
    outputVertex.close();
}

void mySkeleton::saveSurfPoints(string nameCorrVert){
    // Output skeleton points and the corresponding surface points
    std::ofstream output;
    output.open(nameCorrVert);
    BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton))
      BOOST_FOREACH(vertex_descriptor vd, skeleton[v].vertices)
        output <<  skeleton[v].point << "\n "
                       << get(CGAL::vertex_point, tmesh, vd)  << "\n";
    output.close();
}

void mySkeleton::saveData (string chosenFile){
    string nameOFF = chosenFile + ".off";
    std::ifstream input(nameOFF);
    Polyhedron tmesh;
    input >> tmesh;
    Skeleton skeleton;
    CGAL::extract_mean_curvature_flow_skeleton(tmesh, skeleton);
    std::cout << "Number of vertices of the skeleton save: " << boost::num_vertices(skeleton) << "\n";
    string nameSkel = chosenFile + "_skel.cgal";
    std::ofstream output(nameSkel);
    BOOST_FOREACH(Skeleton_edge e, edges(skeleton))
    {
      const Point& s = skeleton[source(e, skeleton)].point;
      const Point& t = skeleton[target(e, skeleton)].point;
      output << s << "\n" << t << "\n";
    }
    output.close();
    // Output all the vertex of the skeleton.
    string nameVertex = chosenFile + "_vertex.cgal";
    std::ofstream output2(nameVertex);
    BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton))
     // BOOST_FOREACH(vertex_descriptor vd, skeleton[v].vertices)
        output2 << skeleton[v].point << "\n";
    output2.close();
    // Output skeleton points and the corresponding surface points
    string nameSurfSkel = chosenFile + "_correspondance-poly.cgal";
    std::ofstream output3(nameSurfSkel);
    BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton))
      BOOST_FOREACH(vertex_descriptor vd, skeleton[v].vertices)
        output3 << skeleton[v].point << "\n"
                       << get(CGAL::vertex_point, tmesh, vd)  << "\n" ;
    output3.close();
}

