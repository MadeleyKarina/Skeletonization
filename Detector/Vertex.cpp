#include "Vertex.h"
#include "util.h"
#include <iostream>
#include <cmath>
#include <map>

void Vertex::getNeighborhood(int rad, vector<Vertex*>& V, Vertex* vertices){
	queue<Vertex*> Q;
	vector<Vertex*> marked;
	
	Q.push(this);
	this->setMark(true);
	this->setDepth(0);
	marked.push_back(this);
	
	while(!Q.empty()){
		Vertex* v0 = Q.front();
		Q.pop();
		V.push_back(new Vertex(v0->getX(), v0->getY(), v0->getZ())); //Indeed, copy vertex information rather than return the same vertex
		
		int dep = v0->getDepth();
		if(dep <= rad){
			set<int> listVertices = v0->getAdjacentVertices();
			set<int> :: iterator it;
			for(it = listVertices.begin(); it!=listVertices.end(); it++){
				Vertex* v1 = &vertices[*it];
					if(!v1->isMarked()){
						Q.push(v1);
						v1->setMark(true);
						v1->setDepth(dep + 1);
						marked.push_back(v1);
					}
				}
			}
	}
	
	vector<Vertex*>::iterator ini = marked.begin();
	
	while(ini<marked.end()){
		(*ini)->setMark(false);
		(*ini)->setDepth(0);
		ini++;
	}
}

int Vertex :: getRadius(Vertex*  vertices, double radius, vector<Vertex*>& V){
	vector<Vertex*> marked; //Store the marked vertices
	map<int, double> distances; //Store the distances relatives to the current vertex
	map<int, Vertex*> markedRing; //Elements in a new ring
	queue<Vertex*> Q; 
	double maxDistance = 0.0; 
	int rad = -1; 
	
	Q.push(this);
	this->setMark(true);
	markedRing.insert(pair<int, Vertex*>(this->index, this));
		
	distances[this->index] = 0.0;
	
	while(!Q.empty()){
		Vertex* v0 = Q.front();
		Q.pop();
		
		int dep = v0->getDepth();
		if(dep != rad){ //First vertex in the new ring
			map<int, Vertex*>::iterator it;
			double max = 0.0;
			
			//Mark the previous ring
			for(it = markedRing.begin(); it!=markedRing.end(); it++){
				Vertex* mar = (*it).second;
				mar->setMark(true);
				marked.push_back(mar);
				V.push_back(new Vertex(mar->getX(), mar->getY(), mar->getZ()));
				if(distances[(*it).first] > max)
					max = distances[(*it).first];
			}
			
			rad++;
			markedRing.clear(); 
			maxDistance = max;
			if(maxDistance > radius)
				break;
		}
		
		set<int> listVertices = v0->getAdjacentVertices();
		set<int> :: iterator it;
		
		for(it = listVertices.begin(); it!=listVertices.end(); it++){
				Vertex* v1 = &vertices[*it];
				if(!v1->isMarked()){ 
					if(distances[v1->getIndex()] == 0.0){ //Distance is not set
						Q.push(v1);   	
						v1->setDepth(dep + 1);
					}
					markedRing.insert(pair<int, Vertex*>(v1->getIndex(), v1));
					double dist = distanceL2(v0, v1);
					double newDistance = distances[v0->getIndex()] + dist;
					if(distances[v1->getIndex()] == 0.0){ //First time on this vertex
						distances[v1->getIndex()] = newDistance;
					}else if(newDistance  < distances[v1->getIndex()]){
						distances[v1->getIndex()] = newDistance;
					}
				}
			}
		//}
	}
	
	if(!markedRing.empty()){
			map<int, Vertex*>::iterator it;
			double max = 0.0;
			
		
			for(it = markedRing.begin(); it!=markedRing.end(); it++){
				Vertex* mar = (*it).second;
				mar->setMark(true);
				marked.push_back(mar);
				V.push_back(new Vertex(mar->getX(), mar->getY(), mar->getZ()));
				if(distances[(*it).first] > max)
					max = distances[(*it).first];
			}
			
			rad++;
			markedRing.clear(); 
			maxDistance = max;
			
	}
	
	//Unmark all vertices
	vector<Vertex*>::iterator ini = marked.begin();
	while(ini  < marked.end()){
		(*ini)->setMark(false);
		(*ini)->setDepth(0);
		ini++;
	}
	return rad;
}

void Vertex::processMaximum(Vertex* vertices, int numRings){
		set<int> :: iterator it;
		for(it = adjacentVertices.begin(); it!=adjacentVertices.end(); it++){
			Vertex* v1 = &vertices[*it];
			if(v1!=this){
				if(response < v1->getResponse())
					return;
			}
		}
	isInterest = true;
}

void Vertex :: getPatch(Vertex* vertices, vector<int> indices, set<int>& returned, set<int>& faceR, double radius, Vertex center){
	set<int> waiting;
	queue<int> visited;
	visited.push(this->index); //Este vertice va a la cola
	
	waiting.insert(indices.begin(), indices.end());
	waiting.erase(this->index); // Eliminamos este vertice del conjunto de faltantes
	
	while(!waiting.empty() && !visited.empty()){
		int ind = visited.front();
		visited.pop();
		
		if(!vertices[ind].isMarked()){
			returned.insert(ind);
			vertices[ind].setMark(true);
			waiting.erase(ind);
			
			set<int> listVertices = vertices[ind].getAdjacentVertices();
			set<int> :: iterator it;
			
			for(it = listVertices.begin(); it!=listVertices.end(); it++){
				int ind = *it;
				if(!vertices[ind].isMarked()){
					double distX = vertices[ind].x() - center.x();
					double distY = vertices[ind].y() - center.y();
					double distZ = vertices[ind].z() - center.z();
					double dist = sqrt(distX*distX + distY*distY + distZ*distZ);
					if(dist < radius)
						visited.push(*it);
				}
			}
			
			vector<int> fac = vertices[ind].getFaces();
			vector<int>::iterator it1;
			for(it1 = fac.begin(); it1!=fac.end(); it1++)
				faceR.insert(*it1);
		}
		
	}
	
	set<int>::iterator it2;
	for(it2 = returned.begin(); it2!=returned.end(); it2++)
		vertices[*it2].setMark(false);
	
}

ostream& operator<<(ostream& out, Vertex& point){
	out << point.x() <<" "<<point.y()<<" "<<point.z()<<endl;
	
	return out;
}
