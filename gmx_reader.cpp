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
#include <xdrfile/xdrfile.h>
#include <xdrfile/xdrfile_xtc.h>
#include "gmx_reader.h"

using namespace std;

// Define class functions
// Default constructor
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
gmx_reader::gmx_reader( string _inpf_ )
{
    gmx_reader::read_in_param( _inpf_ );
    gmx_reader::xtcf_init( startTime );
}

// Default Destructor
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
gmx_reader::~gmx_reader()
{
    delete[] x;
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// Reading Trajectory
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// read next frame
void gmx_reader::read_next()
{
    int xdrinfo = read_xtc( trj, natoms, &step, &gmxtime, box, x, &prec );
    if ( xdrinfo != 0 )
    { 
        cout << "WARNING:: read_xtc returned error " << xdrinfo << "." << " Is the trajectory long enough?" << endl;
        cout << "gmxtime: " << gmxtime;
        exit(EXIT_FAILURE);
    }
}

// search for frame at given time
void gmx_reader::search_for_time(float time)
{
    do 
    {
        read_next();
    }
    while ( fabs(time - gmxtime) > 0.001 ); // allow some tolerance for numerical precision
}

// advance to next sample
void gmx_reader::search_for_sample(int sample)
{
    do
    {
        read_next();
        //printf("%.8f %.8f %.8f\n", gmxtime, startTime + sampleEvery*sample, gmxtime - (startTime + sampleEvery*sample));
    }
    while (fabs(gmxtime - (startTime + sampleEvery*sample))>0.001 ); // allow some tolerance for numerical precision
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// Useful vector operations
// *****************************************************
void gmx_reader::minImage( float dx[3] )
{
    int i;
    for ( i = 0; i < 3; i++ )
    {
        dx[i] -= box[i][i]*round(dx[i]/box[i][i]);
    }
}

// The magnitude of a 3 dimensional vector
float gmx_reader::mag3( float dx[3] )
{
    return sqrt( dot3( dx, dx ) );
}

// The dot product of a 3 dimensional vector
float gmx_reader::dot3( float x[3], float y[3] )
{
    return  x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

// The cross product of two three dimensional vectors
void gmx_reader::cross3( float in1[3], float in2[3], float out[3] )
{
    float i, j, k;
    out[0] = in1[1] * in2[2] - in1[2] * in2[1];
    out[1] = in1[2] * in2[0] - in1[0] * in2[2];
    out[2] = in1[0] * in2[1] - in1[1] * in2[0];
}
// *******************************************************

// Read parameter file for trajectory
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void gmx_reader::read_in_param( string _inpf_ )
{
    string line, para, value;

    ifstream inpf( _inpf_.c_str() );
    if ( !inpf.is_open() )
    {
        cerr << "ERROR: Could not open " << _inpf_ << ". The first argument should contain  a  vaild\nfile name that points to a file containing the simulation parameters.";
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

        if      ( para.find("xtcf") != string::npos )       xtcf = value;
        else if ( para.find("natoms_mol")  != string::npos ) natoms_mol  = stoi(value);
        else if ( para.find("nsamples")    != string::npos ) nsamples    = stoi(value);
        else if ( para.find("sampleEvery") != string::npos ) sampleEvery = stof(value);
        else if ( para.find("startTime")   != string::npos ) startTime   = stof(value);
        // We could just save these to access later if necessary
        else {
            if ( nuParams == nuParamsMax - 1 ) {
                cout << "To many uParams. Max is: " << nuParams << endl;
            }
            cout << "\tWARNING: Parameter " << para << " in input file " << _inpf_ << " not recognized by gmx_reader. Saving as uParam[" << nuParams << "] with uValue " << value << "." << endl ;
            uParams[nuParams] = para;
            uValues[nuParams] = value;
            nuParams += 1;
        }
    }

    inpf.close();
    cout << "\tSetting xtcf file to: " << xtcf.c_str() << endl;
    cout << "\tSetting natoms_mol to: " << natoms_mol << endl;
    cout << "\tSetting nsamples to: " << nsamples << endl;
    cout << "\tSetting sampleEvery to: " << sampleEvery << endl;
    cout << "\tSetting startTime to: " << startTime << endl;
    printf(">>> Done reading input file and setting parameters\n");
    }
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


// Initialize file for reading and set relevant variables and allocate memory
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void gmx_reader::xtcf_init(float startTime)
// Initialize the trajectory file for reading
{
    cout << "Will open and read trajectory from: " << xtcf << endl;
    trj = xdrfile_open( xtcf.c_str(), "r");
    read_xtc_natoms( (char *) xtcf.c_str(), &natoms );
    cout << "Found " << natoms << " atoms." << endl; 
            
    // Set some variables and allocate space for molecular positions
    nmol = natoms/natoms_mol;
    x      = new rvec[ natoms ];

    // Fast forward to the desired frame
    gmx_reader::search_for_time( startTime );
    cout << "Fast forward trajectory to: " << gmxtime << " ps." << endl;
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
