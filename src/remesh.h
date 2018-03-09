#ifndef REMESH_H
#define REMESH_H
#include <iostream>
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
typedef CGAL::Simple_cartesian<double> KernelRemesh;
typedef CGAL::Polyhedron_3<KernelRemesh> Surface_mesh_remesh;
namespace SMS = CGAL::Surface_mesh_simplification ;

int newRemesh(std::string ,std::string , std::string, int) ;
void remeshDataBase(int , int );

#endif // REMESH_H
