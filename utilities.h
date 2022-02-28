#include <cstddef>
#include <cmath>

static const double PI = 3.14159265358979323846;
static const double rootPI = sqrt(PI);
static const double kappa = 1.0;
static const int P = 3; // fixed for tricubic approximation
static const int Pflat = (P + 1)*(P + 1)*(P + 1);
const bool UseSleep = false; // for testing memory usage purpose

extern int N_cube; // number of particles

extern int fdiff, mgrid;
extern double hgrid, rdr, twoh, rtwoh, hsq, rfhsq, r8hcb;
extern double* fval; //Pointer to double

extern double sweight; //weight of source particle in cluster-particle
extern double sfweight; //scaled weight of source particle in cluster-particle
extern int sID; //index of source particle in cluster-particle
extern double *cpvelo; // cp treecode approximation of the potential
extern double *cpdvx; // cp treecode approximation of the x-component of the field
extern double *cpdvy; // cp treecode approximation of the y-component of the field
extern double *cpdvz; // cp treecode approximation of the z-component of the field

extern double **Binv; //B-inverse matrix

//**********************************************************//

struct vec_4d
{
    double val[4];
};


//**********************************************************//

struct xyz // particle coordinates (physical)
{
	double* x;
	double* y;
	double* z;
	size_t* index;
	size_t* old_index;
	size_t size;
	xyz(size_t N_cube_in)
	{
		size = N_cube_in;
		x = new double[size];
		y = new double[size];
		z = new double[size];
		index = new size_t[size];
		old_index = new size_t[size];
	}
	~xyz()
	{
		delete[] x;
		delete[] y;
		delete[] z;
		delete[] index;
		delete[] old_index;
	}
};

//**************//

long getTickCount();

void init_vec(double vec[]);

void read_particle_data(int N_cube,
			struct xyz &particles,
			double *lambda);

double fEval(double rrr);
double f1derivOrd2(double r2h, double t);

