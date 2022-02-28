#include <vector>
#include <string>

#include "kernel_utils.h"

using namespace std;

extern size_t node_count;

extern int N0; // leaf size
extern double mtheta;
extern double sq_theta; // theta^2
extern int max_level;

extern int MACflag;
extern double maxCradius; 

//**************//

struct panel
{
  size_t members[2];
  double xinterval[2];
  double yinterval[2];
  double zinterval[2];
  double xc; // panel center x coordinate
  double yc; // panel center y coordinate
  double zc; // panel center z coordinate
  double rxl; // inverse of panel length parallel to x axis
  double ryl; // inverse of panel length parallel to y axis
  double rzl; // inverse of panel length parallel to z axis
  vector<size_t> children;
  double radius;
  double MAC; // r^2 / theta^2
  double moments[Pflat];
  double dxmoments[Pflat];
  double dymoments[Pflat];
  double dzmoments[Pflat];

  int moment_flag;
  panel() // initialization
  {
    moment_flag = 0;
    members[0] = 0;
    members[1] = -1;
    for (size_t kk = 0; kk < Pflat + 1; kk++){
      moments[kk] = 0.0;
      dxmoments[kk] = 0.0;
      dymoments[kk] = 0.0;
      dzmoments[kk] = 0.0;
   }
  }
};

extern vector<panel> tree;
extern vector<size_t> leaf;
//********************************************************************//

void read_tree_params(string &treeMethod,
                      string &kernelName,
		      int &N_cube,
		      int &N0,
		      double &theta);

void read_direct_sum(const string& kernelName,
		     int N_cube,
		     double v_true[],
                     double fx_true[],
                     double fy_true[],
                     double fz_true[],
		     long &ds_cpu_time);

void build_tree_init(int newMAC, struct xyz &particles);

void Swap(size_t i, size_t j, double *lambda, struct xyz &s);

void split_tree_node(int newMAC,
		     size_t panel_index,
		     double *lambda,
		     struct xyz &particles);

void build_tree_3D_Recursive(int newMAC,
			     size_t panel_index,
			     double *lambda,
			     struct xyz &particles,
			     int level);

vec_4d Call_Ds_PC(int limit_1, int limit_2, 
	       double p_x, double p_y, double p_z,
	       struct xyz &particles,
	       double *lambda, const Kernel& kernel);

void Call_Ds_CP(int limit_1, int limit_2, 
               double p_x, double p_y, double p_z,
               struct xyz &particles,
               const Kernel& kernel);

void write_treecode_sum(const string& kernelName,
			const string &type,
			int N_cube,
			long ds_cpu_time,
			long treecode_cpu_time,
			const double v_true[],
                        const double fx_true[],
                        const double fy_true[],
                        const double fz_true[],
			const double velo[],
                        const double dvx[],
                        const double dvy[],
                        const double dvz[]);

void compute_error(const double v_true[],
                   const double fx_true[],
                   const double fy_true[],
                   const double fz_true[],
		   const double velo[],
                   const double dvx[],
                   const double dvy[],
                   const double dvz[],
		   double &e_n,
		   double &e_d,
		   double &e_d_ex,
		   double &E,
                   double &fd,
                   double &fe,
                   double &FE,
                   double &eL1,
                   double &fL1);
