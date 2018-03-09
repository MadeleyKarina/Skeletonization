#ifndef MYSKELETON_H
#define MYSKELETON_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Mean_curvature_flow_skeletonization.h>
#include <fstream>
#include <boost/foreach.hpp>
#include <iostream>

typedef CGAL::Simple_cartesian<double>                        Kernel;
typedef Kernel::Point_3                                       Point;
typedef CGAL::Polyhedron_3<Kernel>                            Polyhedron;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron> Skeletonization;
typedef Skeletonization::Skeleton                             Skeleton;
typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
typedef Skeleton::edge_descriptor                             Skeleton_edge;

typedef CGAL::Surface_mesh<Point>                                   Triangle_mesh_ob;
typedef boost::graph_traits<Triangle_mesh_ob>::vertex_descriptor    vertex_descriptor_ob;
typedef CGAL::Mean_curvature_flow_skeletonization<Triangle_mesh_ob> Skeletonization_ob;
typedef Skeletonization_ob::Skeleton                                Skeleton_ob;
typedef Skeleton_ob::vertex_descriptor                           Skeleton_vertex_ob;
typedef Skeleton_ob::edge_descriptor                             Skeleton_edge_ob;

using namespace std;
using namespace Eigen;

class mySkeleton{
  private:
    bool flagProperties;
    Polyhedron tmesh;
    Skeleton skeleton;
    Triangle_mesh_ob tmesh_ob;
    Skeleton_ob skeleton_ob;

  public:
    mySkeleton();
    mySkeleton(string , bool);
    void savePolylines(string);
    void saveEdges(string);
    void saveVertex(string);
    void saveSurfPoints(string);
    void saveData(string );
};

#endif
