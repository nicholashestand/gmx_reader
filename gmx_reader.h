#include <string.h>

using namespace std;

#ifndef gmx_H 
#define gmx_H
class gmx_reader
{
    protected:
        // Class variables
        XDRFILE *trj;
        static const int nuParamsMax = 100;
        int64_t *frame_offset;

    public:
        // Class variables
        string  xtcf, offsetf;
        rvec    *x;
        int     natoms, nmol, natoms_mol, nframes, nuParams, step;
        matrix  box;
        float   gmxtime, prec;
        double  dt;
        string  uParams[nuParamsMax], uValues[nuParamsMax];

        // Default constructor and destructor
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        gmx_reader( string _inpf_ );
        ~gmx_reader();

        // Reading Trajectory
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        void read_next_frame();
        void find_frame(int frame);
        void read_frame(int frame);
        void index_frames();
        int  get_frame_number(double time);
        bool checktime(float time);
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        // Useful vector operations
        // *****************************************************
        void  minImage( float dx[3] );
        float mag3( float dx[3] );
        float dot3( float x[3], float y[3] );
        void  cross3( float in1[3], float in2[3], float out[3] );
        // *******************************************************

    private:
        // Read parameter file for trajectory
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        void read_in_param( string _inpf_ );
        void xtcf_init();
        void write_offsets();
        void read_offsets();
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
};
#endif
