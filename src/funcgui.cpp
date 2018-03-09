#include "funcgui.h"

#include <iostream>

using namespace std;
funcGUI::funcGUI()
{
}

void funcGUI::showDirectory(int a){
    DIR *pdir = NULL; // pointer to a directory
    pdir = opendir ("."); // "." will refer to the current directory
    struct dirent *pent = NULL;

    if (pdir == NULL) // if pdir wasn't initialised correctly
    {
        cout << "\nERROR! pdir could not be initialised correctly";
        exit (3);
    } // end if

    while (pent = readdir(pdir)) // while there is still something in the directory to list
    {
        if (pent == NULL) // if pent has not been initialised correctly
        { // print an error message, and exit the program
            cout << "\nERROR! pent could not be initialised correctly";
            exit (3);
        }
        // otherwise, it was initialised correctly. Let's print it on the console:
        cout << pent->d_name << endl;
       // cout<< "aaaaa" <<endl;
    }
    // finally, let's close the directory
    closedir (pdir);
    //cin.get (); // pause for input
    //return EXIT_SUCCESS; // everything went OK
}
