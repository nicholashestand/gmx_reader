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
    int frameno;

    if ( argc != 2 ){
        printf("Program expects only one argument, which is the name of \nan input file containing the details of the analysis.\nAborting...\n");
        exit(EXIT_FAILURE);
    }

    // get filename
    string inpf(argv[1]);
    // attempt to initialize reading of gmx file
    gmx_reader reader( inpf );

    // do some checks of random access reading
    frameno = reader.get_frame_number( 10 );
    cout <<  frameno << " " << reader.nframes << endl;
    if ( frameno != -1 ){
        reader.find_frame(frameno);
        reader.read_next_frame();
        cout << "frame: " << frameno << " time: " << reader.gmxtime << " (ps)" << endl;
    }
    for (int i = frameno; i < frameno + 10; i ++ )
    {
        reader.read_next_frame();
        cout << "time of next frame: " << reader.gmxtime << " (ps)" << endl;
    }

    frameno = reader.get_frame_number( 100.2 );
    if ( frameno != -1 ){
        reader.read_frame(frameno);
        cout << "frame: " << frameno << " time: " << reader.gmxtime << " (ps)" << endl;
    }

    frameno = reader.get_frame_number( 1. );
    if ( frameno != -1 ){
        reader.read_frame(frameno);
        cout << "frame: " << frameno << " time: " << reader.gmxtime << " (ps)" << endl;
    }
}
