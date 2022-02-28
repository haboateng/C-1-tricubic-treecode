#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/times.h>
#include <cmath>

#include "utilities.h"

using namespace std;

int N_cube; // number of particles

int fdiff, mgrid;
double hgrid, rdr, twoh, rtwoh, hsq, rfhsq, r8hcb;
double* fval = NULL; //Pointer to double
double sweight, sfweight;
int sID;
double* cpvelo = NULL;
double* cpdvx  = NULL;
double* cpdvy  = NULL;
double* cpdvz  = NULL;
double** Binv = NULL;

//************************************************************************//

long getTickCount()
{
  tms tm;
  return times(&tm);
}

//************************************************************************//

// initialize a vector --- assumes size is 64

void init_vec(double vec[])
{
  for (int i=0; i<64; i++)
    vec[i]=0.0;
}

//************************************************************************//

// read particle coordinates and weights from files

void read_particle_data(int N_cube,
			struct xyz &particles,
			double *lambda)
{
  FILE * fp;
  char S_data_file[64] = {0};
  sprintf(S_data_file, "./rand_%d.txt", N_cube);
  fp = fopen(S_data_file, "r");
  
  double x1, x2, x3;
  int count = -1;
  if (fp == NULL)
    {
      cout << "Cannot open random points file" << endl;
      exit(-1);
    }
  
  while (true)
    {
      count++;
      fscanf(fp,"%lf%lf%lf", &x1, &x2, &x3);
      if (feof(fp))
	break;
      if (count >= N_cube)
	{
	  cout << "Out of range" << endl;
	  exit(-1);
	}
      
      particles.x[count] = x1;
      particles.y[count] = x2;
      particles.z[count] = x3;
      particles.index[count] = -1;
      particles.old_index[count] = count;
    }
  
  
  // ********* read particle weights from a file *********
  
  char lambda_Str_data_file[64] = {0};
  sprintf(lambda_Str_data_file, "./lambda_%d.txt", N_cube);
  fp = fopen(lambda_Str_data_file, "r");
  
  count = -1;
  if (fp == NULL)
    {
      cout << "Cannot open lambda file" << endl;
      exit(-1);
    }
  
  while (true)
    {
      count++;
      fscanf(fp,"%lf", &x1);
      if (feof(fp))
	break;
      if (count >= N_cube)
	{
	  cout << "Out of range" << endl;
	  exit(-1);
	}
      
      lambda[count] = x1;
    }
}

//*********************************************************************//

// Use 3-point interpolation to compute function value

double fEval(double rrr) 
{
  int ll     = static_cast<int>(rrr*rdr);
  int l1     = ll+1;
  int l2     = ll+2;
  double ppp = rrr*rdr-static_cast<double>(ll); 

  double vk0 = fval[ll];
  double vk1 = fval[l1];
  double vk2 = fval[l2];
  double t1  = vk0+(vk1-vk0)*ppp;
  double t2  = vk1+(vk2-vk1)*(ppp-1.0);

  double fvalue = (t1+(t2-t1)*ppp*0.5);
  return (t1+(t2-t1)*ppp*0.5);

}
//********************************************************************//

// Second order finite difference

double f1derivOrd2(double r2h, double t)
{
  double tth = t*twoh;
  return rtwoh*(fEval(sqrt(r2h - tth)) - fEval(sqrt(r2h + tth)));
  // the input is (x_b - x) and we want the derivative wrt to x,
  // so it's the opposite:
  // i.e. (f(x-h) - f(x+h))/2h
}
