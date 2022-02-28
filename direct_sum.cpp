
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

#include "kernel_utils.h"

using namespace std;

//************************************************************************//

void Ds_FULL(int N_cube,
	     struct xyz &particles,
	     double *lambda,
	     double v_true[],
             double fx_true[],
             double fy_true[],
             double fz_true[],
	     long &ds_cpu_time,
	     const Kernel& kernel)
{
  long Start_ds, End_ds;
  Start_ds = getTickCount(); // Get currenct CPU time
  
  for (int i = 0; i < N_cube-1; i++)
    {
      double temp_x = particles.x[i];
      double temp_y = particles.y[i];
      double temp_z = particles.z[i];
      double sum  = 0.0;
      double sumx = 0.0;
      double sumy = 0.0;
      double sumz = 0.0;
      double wi   = lambda[i];
      
      for (int j = i+1; j < N_cube; j++)
	{
	  double xx = temp_x - particles.x[j];
	  double yy = temp_y - particles.y[j];
	  double zz = temp_z - particles.z[j];

          double wj = lambda[j];
	  double kerVal  = kernel.formula(xx, yy, zz);
	  double kerValF = kernel.Fformula(xx, yy, zz);
	  
          double fx = xx*kerValF; 
          double fy = yy*kerValF;
          double fz = zz*kerValF;

	  sum += wj*kerVal;
	  v_true[j] += wi*kerVal;

          sumx += wj*fx;
          sumy += wj*fy;
          sumz += wj*fz;

          fx_true[j] += wi*fx;
          fy_true[j] += wi*fy;
          fz_true[j] += wi*fz;  
	}
      
      v_true[i] += sum;
      fx_true[i] -= sumx;
      fy_true[i] -= sumy;
      fz_true[i] -= sumz;
    }
  
  End_ds = getTickCount(); // Get currenct CPU time
  ds_cpu_time = End_ds - Start_ds;
} 

//*************************************************************************//

int main()
{
  ifstream ifile("input_params.txt");
  string treeMethod;
  string kernelName;
  if (ifile.is_open())
    { 
      string line;
      getline(ifile,line);
      getline(ifile,line);
      stringstream stream;
      stream << line << endl;
      stream >> treeMethod;
      getline(ifile,line);

      getline(ifile,line);
      getline(ifile,line);
      stream << line << endl;
      stream >> kernelName;
      getline(ifile,line);

      getline(ifile,line);
      getline(ifile,line);
      stream << line << endl;
      stream >> N_cube;
      getline(ifile,line);
    }
  ifile.close();
  Kernel* kernel;
  if (kernelName == "Coulomb") {
    kernel = new CoulombKernel();
  } else if (kernelName == "ScreenedCoulomb") {
    kernel = new ScreenedCoulombKernel();
  } else if (kernelName == "RSEwaldSum") {
    kernel = new RSEwaldSumKernel();
  } else {
    cerr << "Invalid kernel name: " << kernelName << endl;
    exit(-1);
  }
  cout << "Kernel = " << kernelName << endl;
  cout << "N_cube = " << N_cube << endl;
  struct xyz particles(N_cube);
  double *lambda = new double[N_cube];
  
  read_particle_data(N_cube, particles, lambda);

  
  //***************** Direct summation *****************
  
  double *v_true = new double[N_cube];
  double *fx_true = new double[N_cube];
  double *fy_true = new double[N_cube];
  double *fz_true = new double[N_cube];

  for (int i = 0; i < N_cube; i++)
    {
      v_true[i]  = 0.0;
      fx_true[i] = 0.0;
      fy_true[i] = 0.0;
      fz_true[i] = 0.0;
    }
  long ds_cpu_time;
  
  Ds_FULL(N_cube,
	  particles,
	  lambda,
	  v_true,
          fx_true,
          fy_true,
          fz_true,
	  ds_cpu_time,
	  *kernel);
  
  cout << "ds time in seconds*100 = " << ds_cpu_time << endl;
  
  
  // write results to a file
  
  char file_name[256];
  sprintf(file_name, "exact_sum_%s_N%d", kernelName.c_str(), N_cube);
  ofstream output_file(file_name);
  
  for (int i = 0; i < N_cube; i++)
    {
      output_file << i << "   " << setprecision(16) << v_true[i] <<"   "
                  <<fx_true[i]<<"   "<<fy_true[i]<<"   "<<fz_true[i]<<endl;
    }
  
  output_file << setprecision(16) << ds_cpu_time << endl;
  output_file.close();	   

  delete [] lambda;
  delete [] v_true;
  delete [] fx_true;
  delete [] fy_true;
  delete [] fz_true;

  delete kernel;
  
  return 0;
}
