#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <xdrfile/xdrfile.h>
#include <xdrfile/xdrfile_xtc.h>
#include <gmx_reader.h>

int main( int argc, char* argv[] )
{
    int currentSample = 0;

    if ( argc != 2 ){
        printf("Program expects only one argument, which is the name of \nan input file containing the details of the analysis.\nAborting...\n");
        exit(EXIT_FAILURE);
    }

    // get filename
    string inpf(argv[1]);
    // attempt to initialize reading of gmx file
    gmx_reader reader( inpf );

    for ( currentSample = 0; currentSample < reader.nsamples; currentSample ++ ){
        if ( currentSample != 0 )
        {
            // if at sample 0, the first frame was already read by the class initialization
            reader.read_next();
        }
        cout << "Current time: " << reader.gmxtime << endl;
    }
}
