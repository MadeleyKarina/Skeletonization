#ifndef __MESH_H_
#define __MESH_H_

#include "Vertex.h"
#include "Face.h"
#include <vector>

using namespace std;

class Mesh{

      //Topological information
      Vertex* vertices;
      Face* faces;
	  
	  int numVertices;
      int numFaces;
	  
	  float** D;
      
      //BBox Information
      double xmin, xmax, ymin, ymax, zmin, zmax;
      double diag;


	  public:
	  		void cleanMesh();
      
             Mesh();
             Mesh(const char* nombreArchivo);
             ~Mesh();

             void loadFromFile(const char* filename);
             
             friend ostream& operator<<(ostream& os, Mesh &obj);
			 Vertex* getVertices(){return vertices;}
			 Face* getFaces(){return faces;}
			 
			 int getNumberVertices(){return numVertices;}
			 int getNumFaces(){return numFaces;}
			 
			 double getDiagonal(){return diag;}
			 double getArea();
			 void getPatch(int seed, vector<int> indices, double radius, Vertex center, Mesh* patch, set<int>& vertReturned);
			 void getSpatialPatch(Mesh* patch, Vertex center, double radius);
			 
			 void setNumberVertices(int nv);
			 void setNumberFaces(int nf);
			 
			 void insertVertex(int pos, double x, double y, double z);
			 void insertFace(int pos, int p1, int p2, int p3);
			 			 
};
#endif
