#ifndef __CLOCK_H_
#define __CLOCK_H_
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

class Clock{
      private:
              clock_t current;
              clock_t previous;    
              
      public:
             Clock(){previous = current = clock();}
             virtual ~Clock(){}
             
             void tick(){
                  previous = current;     
                  current = clock();
             }
             
             double getTime(){
                    return (((double)(current - previous))/CLOCKS_PER_SEC);
             }
             
             string getFormattedTime(){
                 stringstream aux;
				 unsigned long t = (unsigned long)getTime();
				 int hour = t / 3600;
				 int min = (t%3600)/60;
				 int sec = (t%3600)%60;
				 
				 aux << hour << " hours - " << min << " min. - " << sec << " sec.";
				 return aux.str();
			 }
};

#endif
