/*
 * This is a generic class to read a gromacs trajectory
 *
 */

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <xdrfile.h>
#include <xdrfile_xtc.h>
#include <xtc_seek.h>
#include "gmx_reader.h"

using namespace std;

// Define class functions
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

gmx_reader::gmx_reader( string _inpf_ )
// Default constructor
{
    gmx_reader::read_in_param( _inpf_ );
    gmx_reader::xtcf_init( );
}

gmx_reader::~gmx_reader()
// Default Destructor
{
    // close the xdrfile
    xdrfile_close( trj );

    // free memory
    delete[] x;
    delete[] frame_offset;
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


// *******************************************************
// Read parameter file for trajectory
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

void gmx_reader::read_in_param( string _inpf_ )
// read the input file to get the parameters for the analysis
{
    string line, para, value;

    ifstream inpf( _inpf_.c_str() );
    if ( !inpf.is_open() )
    {
        cerr << "ERROR: Could not open " << _inpf_ << \
            ". The first argument should contain  a  vaild\n"<<
            "file name that points to a file containing the simulation parameters.";
        exit(EXIT_FAILURE);
    }
    cout << ">>> Reading parameters from input file: " << _inpf_ << endl;

    // initialize user parameters
    nuParams = 0;

    // Parse input file
    while ( getline( inpf, line ) )
    {
        // create a string stream and parse the parameter name and parameter value
        istringstream line2(line);
        line2 >> para >> value;

        // skip comment lines
        if (para[0] == '#') continue;

        if      ( para.find("xtcf") != string::npos )        xtcf        = value;
        else if ( para.find("offsetf") != string::npos )     offsetf     = value;
        else if ( para.find("natoms_mol")  != string::npos ) natoms_mol  = stoi(value);

        // Save unknown parameters for access later
        else {
            if ( nuParams == nuParamsMax - 1 ) {
                cout << "To many uParams. Max is: " << nuParams << endl;
            }
            cout << "\tWARNING: Parameter " << para << " not recognized. Saving as uParam["\
                << nuParams << "] with uValue " << value << "." << endl ;
            uParams[nuParams] = para;
            uValues[nuParams] = value;
            nuParams += 1;
        }
    }

    // check that offset file has a name, and if not, set it to a default
    if ( offsetf.empty() ) offsetf = xtcf + ".fnx";
    
    inpf.close();
    cout << "\tSetting xtcf file to: " << xtcf.c_str() << endl;
    cout << "\tSetting offset file to: " << offsetf.c_str() << endl;
    cout << "\tSetting natoms_mol to: " << natoms_mol << endl;
    printf(">>> Done reading input file and setting parameters\n");
}

void gmx_reader::xtcf_init()
// Initialize the trajectory file for reading
{
    int est_nframes;
    float time;

    // Get the number of atoms, molecules, and allocate space for the positions
    cout << "Will open and read trajectory from: " << xtcf << endl;
    trj = xdrfile_open( xtcf.c_str(), "r");
    read_xtc_natoms( (char *) xtcf.c_str(), &natoms );
    cout << "Found " << natoms << " atoms." << endl; 
    nmol = natoms/natoms_mol; 
    x    = new rvec[ natoms ];

    // find frame time offsets -- assume they are regular;
    read_xtc( trj, natoms, &step, &gmxtime, box, x, &prec );
    time = gmxtime;
    read_xtc( trj, natoms, &step, &gmxtime, box, x, &prec );
    dt = round((gmxtime - time)*(prec/10))/((prec/10)); // make slightly less than full precision to avoid rounding problems
    cout << "Frame time offset is: " << dt << " (ps)" << endl;

    // close the xdr file
    xdrfile_close( trj );


    // Create a frame index for random frame access of the trajectory file
    // check if the frame index exists
    ifstream fin(offsetf.c_str());
    if ( fin ){
        // if it exists, read the offsets
        cout << "Reading indexes for " << xtcf << ". ";
        read_offsets();
        cout << "Done." <<  endl;

    }
    else{
        // if it doesnt exist, we need to index the frames for random access
        cout << "Now indexing the xtc file: " << xtcf \
            << ". This may take some time, be patient." << endl;
        read_xtc_n_frames( (char *) xtcf.c_str(), &nframes, &est_nframes, &frame_offset );
        // write them to a file so you dont have to index it again
        write_offsets();
        cout << "Done indexing the xtc file. Next time will be faster." << endl;
    }

    // Open the xdr file for random access reading
    trj = xdrfile_open( xtcf.c_str(), "r");
}

void gmx_reader::write_offsets()
// write frame offsets to file
{
    FILE *file = fopen( offsetf.c_str(), "wb" );
    fwrite( &nframes, sizeof(int), 1, file );
    fwrite( frame_offset, sizeof(int64_t), nframes, file );
    fclose(file);
}

void gmx_reader::read_offsets()
// read frame offsets from file
{
    FILE *file = fopen( offsetf.c_str(), "rb" );
    fread( &nframes, sizeof(int), 1, file );
    // allocate frame_offset array
    frame_offset = new int64_t[ nframes ];
    fread( frame_offset, sizeof(int64_t), nframes, file);
    fclose(file);
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// Reading Trajectory
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void gmx_reader::read_next_frame()
// read next frame
{
    int xdrinfo;

    xdrinfo = read_xtc( trj, natoms, &step, &gmxtime, box, x, &prec );
    if ( xdrinfo != exdrOK ){
        cout << "WARNING:: read_xtc returned error " << xdrinfo << "." << endl;
        xdrfile_close( trj );
        exit(EXIT_FAILURE);
    }
}

void gmx_reader::find_frame(int frame)
// set file pointer to the beginning of a certain frame
{
    int xdrinfo;
    int64_t offset;

    offset = frame_offset[frame];
    // set pointer to desired frame
    xdrinfo = xdr_seek( trj, offset, SEEK_SET );
    if ( xdrinfo != exdrOK ){
        cout << "WARNING:: xdr_seek returned error " << xdrinfo << "." << endl;
        xdrfile_close( trj );
        exit(EXIT_FAILURE);
    }
}

void gmx_reader::read_frame(int frame)
// read a specific frame
{
    find_frame(frame);
    read_next_frame();
}

bool gmx_reader::checktime(double time)
// make sure the time we want is the same time as in the gmxtime
{
    if ( fabs( time - gmxtime ) > 1E-2 ) return false;
    else return true;
}

int gmx_reader::get_frame_number(double time)
// returns the frame number for a given time
{
    int frame;

    // check if time is valid -- if it isnt evenly divisible by dt, then that frame is not availible
    if ( fabs(remainder( time, dt )) > 1E-2 ){
        cout << "Warning get_frame_number failed. Aborting. " << endl;
        exit(0);
        return -1;
    }
    // get the frame number from the time
    frame = (int) round(time/dt);
    // check to make sure we have that many frames
    if (frame > nframes) return -1;
    return frame;
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


// Useful vector operations
// *****************************************************
void gmx_reader::minImage( float dx[3] )
// min image of a vector
{
    int i;
    for ( i = 0; i < 3; i++ )
    {
        dx[i] -= box[i][i]*round(dx[i]/box[i][i]);
    }
}

float gmx_reader::mag3( float dx[3] )
// The magnitude of a 3 dimensional vector
{
    return sqrt( dot3( dx, dx ) );
}

float gmx_reader::dot3( float x[3], float y[3] )
// The dot product of a 3 dimensional vector

{
    return  x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

void gmx_reader::cross3( float in1[3], float in2[3], float out[3] )
// The cross product of two three dimensional vectors
{
    float i, j, k;
    out[0] = in1[1] * in2[2] - in1[2] * in2[1];
    out[1] = in1[2] * in2[0] - in1[0] * in2[2];
    out[2] = in1[0] * in2[1] - in1[1] * in2[0];
}
// *******************************************************


