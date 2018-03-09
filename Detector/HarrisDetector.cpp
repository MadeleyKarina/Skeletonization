#include <cassert>
#include <cstdlib>
#include <cfloat>
#include <algorithm>
#include <vector>

#include <Eigen/Dense>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_3.h>

#include "Mesh.h"
#include "HarrisDetector.h"
#include "Clock.h"


typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Search_traits_3<K> Traits;
typedef CGAL::Random_points_in_cube_3<Point_3>       Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Kd_tree<Traits> Tree;
typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;





using namespace std;
using namespace Eigen;

HarrisDetector :: HarrisDetector(){
	typeNeighborhood = ADAPTIVE;
	fractionDiagonal = 0.01;
	numberRingNeighbor = 0;

	k = 0.04;

	numberRingsDetection = 1;

	typeSelection = FRACTION;
	paramSelection = 0.01;

	filteringSteps = 0;
	object = NULL;
	prop = NULL;
}

HarrisDetector :: HarrisDetector(Mesh* obj, Properties* pr){
	object = obj;
	prop = pr;

	processOptions();
}

void HarrisDetector :: processOptions(){
	/*We evaluate the properties*/
	//typeNeighborhood
	string val = prop->getProperty("type-neighborhood");
	if(val.empty()){
		cout << "Option type-neighborhood was not set. ADAPTIVE is assumed." << endl;
		typeNeighborhood = ADAPTIVE;
	}
	else{
		if(!val.compare("spatial"))
			typeNeighborhood = SPATIAL;
		else if(!val.compare("rings"))
			typeNeighborhood = RINGS;
		else if(!val.compare("adaptive"))
			typeNeighborhood = ADAPTIVE;
		else{
			cout << "Option type-neighborhood is not recognized. ADAPTIVE is assumed." << endl;
			typeNeighborhood = ADAPTIVE;
		}
	}

	//typeSelection
	val = prop->getProperty("interest-points-selection");
	if(val.empty()){
		cout << "Option interest-points-detection was not set. FRACTION is assumed." << endl;
		typeSelection = FRACTION;
	}else{
		if(!val.compare("clustering"))
			typeSelection = CLUSTERING;
		else if(!val.compare("fraction"))
			typeSelection = FRACTION;
		else if(!val.compare("number"))
			typeSelection = NUMBER;
		else{
			cout << "Option interest-points-detection is not recognized. FRACTION is assumed." << endl;
		}
	}

	//parameter-neighborhood
	val = prop->getProperty("parameter-neighborhood");
	if(typeNeighborhood == SPATIAL || typeNeighborhood == ADAPTIVE){
		if(val.empty()){
			cout << "Using SPATIAL or ADAPTIVE: Option parameter-neighborhood was not set. Using 0.01 by default" << endl;
			fractionDiagonal = 0.01;
		}
		else
			fractionDiagonal = atof(val.c_str());
	}else if(typeNeighborhood == RINGS){
		if(val.empty()){
			cout << "Using RINGS: Option parameter-neighborhood was not set. Using 2 by default" << endl;
			numberRingNeighbor = 2;
		}
		else
			numberRingNeighbor = atoi(val.c_str());
	}

	//k
	val = prop->getProperty("K");
	if(val.empty()){
		cout << "Parameter K was not set. Using 0.04 by default" << endl;
		k = 0.04;
	}
	else{
		k = atof(val.c_str());
	}

	//numberRingsDetection
	val = prop->getProperty("ring-maxima-detection");
	if(val.empty()){
		cout << "Parameter ring-maxima-detection was not set. Using 1 by default" << endl;
		numberRingsDetection = 1;
	}else{
		numberRingsDetection = atoi(val.c_str());
	}

	//parameter-selection
	val = prop->getProperty("parameter-selection");
	if(typeSelection == FRACTION || typeSelection == CLUSTERING){
		if(val.empty()){
			cout << "Parameter parameter-selection was not set. Using 0.01 by default." << endl;
			paramSelection = 0.01;
		}else{
			paramSelection = atof(val.c_str());
		}
	}else if(typeSelection == NUMBER){
		if(val.empty()){
			cout << "Parameter parameter-selection was not set. Using 100 by default." << endl;
			numberKeypoints = 100;
		}else{
			numberKeypoints = atoi(val.c_str());
		}
	}

	val = prop->getProperty("filtering-steps");
	if(val.empty()){
		cout << "Parameter filtering-steps was not set. Using 0 by default" << endl;
		filteringSteps = 0;
	}else{
		filteringSteps = atoi(val.c_str());
	}
}

void HarrisDetector :: showOptions(){
	cout << "type-neighborhood = ";
	switch(typeNeighborhood){
		case SPATIAL: cout << "SPATIAL" << endl;
					  break;
		case ADAPTIVE: cout << "ADAPTIVE" << endl;
						break;
		case RINGS:	cout << "RINGS" << endl;
					break;
	}

	cout << "parameter-neighborhood = ";
	if(typeNeighborhood == SPATIAL || typeNeighborhood == ADAPTIVE)
		cout << fractionDiagonal << endl;
	else
		cout << numberRingNeighbor << endl;

	cout <<"K = "<< k << endl;
	cout <<"ring-maxima-detection = " << numberRingsDetection << endl;

	cout << "interest-points-selection = ";
	switch(typeSelection){
		case FRACTION: cout << "FRACTION" << endl;
						break;
		case CLUSTERING: cout << "CLUSTERING" << endl;
						break;
	}

	cout << "parameter-selection = " << paramSelection << endl;
}

bool comp(Vertex i, Vertex j){
	return i.getResponse() > j.getResponse();
}

void HarrisDetector :: detectInterestPoints(vector<int>& interestPoints){
	assert(object);

	double max = 0.0;
	int rad;
	int cont = 0;
	double diag = object->getDiagonal();

	Vertex* vertices = object->getVertices();
	Face* faces = object->getFaces();

	cout << "Starting Harris" << endl;
	Tree tree;

	if(typeNeighborhood == SPATIAL){
		vector<Point_3> points;
		for(register int i = 0; i < object->getNumberVertices(); i++){
			points.push_back(Point_3(vertices[i].getX(), vertices[i].getY(), vertices[i].getZ()));
		}
		tree.insert(points.begin(), points.end());
	}


	double* responses = new double[object->getNumberVertices()];
	memset(responses, 0, sizeof(double)*object->getNumberVertices());
	//Process each vertex
	/*This loop computes the Harris response for each vertex*/
	for(register int i = 0; i< object->getNumberVertices(); i++){
		//cout << "Vertex " << i << "->";
		vector<Vertex*> set;

		//Calcular radio adaptativo
		if(typeNeighborhood == ADAPTIVE){
			rad = vertices[i].getRadius(vertices, diag * fractionDiagonal, set);
            cout << "vertices ejemplo: " << endl;
            cout << vertices[1]<< endl;
		}
		else if(typeNeighborhood == SPATIAL){
			Point_3 query(vertices[i].getX(), vertices[i].getY(), vertices[i].getZ());
			Fuzzy_sphere fs(query, diag * fractionDiagonal, 0.0);

			vector<Point_3> puntos;
			tree.search(back_inserter(puntos), fs);

			for(int i = 0; i < puntos.size(); i++){
				Vertex* p = new Vertex(puntos[i].x(), puntos[i].y(), puntos[i].z());
				set.push_back(p);
			}
		}else if(typeNeighborhood == RINGS){
			vertices[i].getNeighborhood(numberRingNeighbor, set, vertices);
		}

		if(set.size() < 6){
			cont++;
			//vertices[i].setResponse(DBL_MIN);
			responses[i] = DBL_MIN;
			continue;
		}
		//Process "set", first, calculate the centroid
		double xc = 0, yc = 0, zc = 0;
		for(register int j = 0; j< set.size();j++){
			xc += set[j]->getX();
			yc += set[j]->getY();
			zc += set[j]->getZ();
		}

		xc /= set.size();
		yc /= set.size();
		zc /= set.size();

		//Translate the vertex, in order the centroid is in [0 0 0]
		for(register int j = 0; j< set.size();j++){
			set[j]->setX(set[j]->getX() - xc);
			set[j]->setY(set[j]->getY() - yc);
			set[j]->setZ(set[j]->getZ() - zc);
		}

		//Aplicar PCA para encontrar una pose de la nube de puntos de manera que la mayor dispersión de los puntos esté en el plano XY
		//La media de las 3 coordenadas ya es (0, 0, 0), así que en realidad no es necesario calcularla, directamente calculamos la
		//matriz de covarianza
		double A[9];
		memset(A, 0, sizeof(double)*9);

		for(register int j = 0; j < set.size(); j++){
			double x = set[j]->getX();
			double y = set[j]->getY();
			double z = set[j]->getZ();

			A[0] += x*x;	A[1] += x*y;	A[2] += x*z;
							A[4] += y*y;	A[5] += y*z;
											A[8] += z*z;
		}
		A[3] = A[1];	A[6] = A[2];	A[7] = A[5];
		for(int j = 0; j < 9; j++)
			A[j] /= (set.size()-1);

		Matrix3d m;
		m(0,0) = A[0];	m(0,1) = A[1];	m(0,2) = A[2];
		m(1,0) = A[3];	m(1,1) = A[4];	m(1,2) = A[5];
		m(2,0) = A[6];	m(2,1) = A[7];	m(2,2) = A[8];

		SelfAdjointEigenSolver<MatrixXd>	eigensolver(m);
		Matrix3d evec = eigensolver.eigenvectors();

		//cout << "The eigenvalues are " << eigensolver.eigenvalues() << endl;
		double x_1 = evec(0,2);	double x_2 = evec(1,2);	double x_3 = evec(2,2);
		double y_1 = evec(0,1);	double y_2 = evec(1,1);	double y_3 = evec(2,1);
		double z_1 = evec(0,0);	double z_2 = evec(1,0);	double z_3 = evec(2,0);

		double x2 = set[0]->getX() - xc;
		double y2 = set[0]->getY() - yc;
		double z2 = set[0]->getZ() - zc;

		if((z_1*x2 + z_2*y2 + z_3*z2) < 0){
			z_1 = -z_1;
			z_2 = -z_2;
			z_3 = -z_3;

			double aux_x1 = x_1;
			double aux_x2 = x_2;
			double aux_x3 = x_3;

			x_1 = y_1;
			x_2 = y_2;
			x_3 = y_3;
			y_1 = aux_x1;
			y_2 = aux_x2;
			y_3 = aux_x3;

		}
		//Realizamos la rotacion, con el nuevo sistema de coordenadas
		for(register int j = 0; j < set.size(); j++){
			double x = set[j]->getX();
			double y = set[j]->getY();
			double z = set[j]->getZ();

			set[j]->setX(x_1*x + x_2*y + x_3*z);
			set[j]->setY(y_1*x + y_2*y + y_3*z);
			set[j]->setZ(z_1*x + z_2*y + z_3*z);
		}

		//Movemos todos los puntos para que el punto de analisis se encuentre en el origen del plano XY. Solo movemos en las coordenadas X e Y
		double x = set[0]->getX();
		double y = set[0]->getY();

		for(register int j = 0; j < set.size(); j++){
			set[j]->setX(set[j]->getX() - x);
			set[j]->setY(set[j]->getY() - y);
		}

		//Fit a quadratic surface
		double C[36];
		double D[6];

		memset(C, 0, sizeof(double)*36);
		memset(D, 0, sizeof(double)*6);

        for(register int j = 0; j < set.size(); j++){
			double x = set[j]->getX();
			double y = set[j]->getY();
			double z = set[j]->getZ();

			double x2 = x*x;
			double y2 = y*y;
			double x3 = x2*x;
			double y3 = y2*y;

			C[0] += x*x3;	C[1] += x3*y;	C[2] += x2*y2;	C[3] += x3;		C[4] += x2*y;	C[5] += x2;
							C[7] += x2*y2;	C[8] += x*y3;	C[9] += x2*y;	C[10] += x*y2;	C[11] += x*y;
											C[14] += y*y3;	C[15] += x*y2;	C[16] += y3;	C[17] += y2;
															C[21] += x2;	C[22] += x*y;	C[23] += x;
																			C[28] += y2;	C[29] += y;
			D[0] += z*x2;	D[1] += z*x*y;	D[2] += z*y2;	D[3] += z*x;	D[4] += z*y;	D[5] += z;
		}

        C[6] = C[1];
		C[12] = C[2];	C[13] = C[8];
		C[18] = C[3];	C[19] = C[9];	C[20] = C[15];
		C[24] = C[4];	C[25] = C[10];	C[26] = C[16];	C[27] = C[22];
		C[30] = C[5];	C[31] = C[11];	C[32] = C[17];	C[33] = C[23];	C[34] = C[29];

		C[35] = (double)set.size();

		MatrixXd m1(6,6);
		VectorXd  b1(6);

		for(int j = 0; j < 6; j++){
			for(int k = 0; k < 6; k++)
				m1(j,k) = C[j*6 + k];
			b1(j) = D[j];
		}

		ColPivHouseholderQR<MatrixXd> dec(m1);
		VectorXd x1 = dec.solve(b1);

		double c20_2 = x1(0);	//p1
		double c11 = x1(1);		//p2
		double c02_2 = x1(2);	//p3
		double c10 = x1(3);		//p4
		double c01 = x1(4);		//p5
		double c0 = x1(5);		//p6

		//double c20 = c20_2*2;
		//double c02 = c02_2*2;
		double c20 = c20_2;
		double c02 = c02_2;

		double fx2 = c10*c10 + 2*c20*c20 + 2*c11*c11; //A
		double fy2 = c01*c01 + 2*c11*c11 + 2*c02*c02; //B
		double fxfy = c10*c01 + 2*c20*c11 + 2*c11*c02; //C

		//double k = 0.04;
		double resp = fx2*fy2 - fxfy*fxfy - k*(fx2 + fy2)*(fx2 + fy2);

		//vertices[i].setResponse(resp);
		responses[i] = resp;

		if(resp > max)
			max = resp;

		for(register int j = 0; j < set.size(); j++)
			delete set[j];
		set.clear();
	}

	double* auxResp = new double[object->getNumberVertices()];

	cout << "Executing filtering in " << filteringSteps << " steps" << endl;
	//Perform the filtering before the use of the responses
	for(int i = 0; i < filteringSteps; i++){
		memcpy(auxResp, responses, sizeof(double)*object->getNumberVertices());

		for(int j = 0; j < object->getNumberVertices(); j++){
			set<int> neighbors = vertices[j].getAdjacentVertices();
			double sum=0.0;
			for(set<int>::iterator it = neighbors.begin(); it!=neighbors.end(); it++){
				sum += auxResp[*it];
			}
			responses[j] = sum/neighbors.size();
		}
	}

	delete[] auxResp;
	for(int i = 0; i < object->getNumberVertices(); i++)
		vertices[i].setResponse(responses[i]);

	delete[] responses;
	vector<Vertex> candidatePoints;

	//Search for local maximum
	for(register int i = 0; i< object->getNumberVertices();i++){
		vertices[i].processMaximum(vertices, numberRingsDetection);
		if(vertices[i].getInterest()){
			candidatePoints.push_back(vertices[i]);
		}
	}

	cout << "Candidates:" << candidatePoints.size() << endl;
	sort(candidatePoints.begin(), candidatePoints.end(), comp);

	if(typeSelection == FRACTION){
		//Seleccionar los puntos de mayor respuesta
		int numPoints = paramSelection * object->getNumberVertices() > candidatePoints.size()? candidatePoints.size():paramSelection*object->getNumberVertices();
		for(register int i = 0; i < numPoints; i++)
			interestPoints.push_back(candidatePoints[i].getIndex());
	}else if(typeSelection == CLUSTERING){
		//Aplicar proceso de Clustering
		for(int i = 0; i < candidatePoints.size(); i++){
			bool found = false;
			int j = 0;
			while(j < interestPoints.size() && !found){
				double distX = vertices[interestPoints[j]].getX() - candidatePoints[i].getX();
				double distY = vertices[interestPoints[j]].getY() - candidatePoints[i].getY();
				double distZ = vertices[interestPoints[j]].getZ() - candidatePoints[i].getZ();

				if(sqrt(distX*distX + distY*distY + distZ*distZ) < (paramSelection * diag))
					found = true;
				j++;
			}
			if(!found)
				interestPoints.push_back(candidatePoints[i].getIndex());
		}
	}else if(typeSelection == NUMBER){
		if(candidatePoints.size() < numberKeypoints)
			numberKeypoints = candidatePoints.size();
		for(int i = 0; i < numberKeypoints; i++)
			interestPoints.push_back(candidatePoints[i].getIndex());
	}

	//for(int i = 0; i < interestPoints.size(); i++)
	//	cout << interestPoints[i].getIndex() << endl;
	//cout << "Interest points:" <<interestPoints.size() << endl;
}
