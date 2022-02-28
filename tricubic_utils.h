
#include "tree_utils.h"

//*******************************************************************//

void Panel_Moment_Tricubic(size_t panel_index,
			   double *lambda,
			   struct xyz &particles,
			   double m[],
                           double mx[], double my[], double mz[]);

void getClusterCorners_B2(size_t panel_index,
			  double &dx, double &dy, double &dz,
			  double ccx[], double ccy[], double ccz[]);

void getClusterCorners_LM(size_t panel_index,
			  double &dx, double &dy, double &dz,
			  double ccx[], double ccy[], double ccz[]);

void getClusterPoints_B6(size_t panel_index,
                          double &dx, double &dy, double &dz,
                          double ccx[], double ccy[], double ccz[]);

void scaleDerivatives_B2(double b[], double dx, double dy, double dz);

void scaleDerivatives_LM(double b[], double dx, double dy, double dz);

void Compute_bvec_B2(double bvec[],
		     double x, double y, double z,
		     size_t panel_index,
		     const Kernel& kernel);

void Compute_bvec_LM(double bvec[],
		     double x, double y, double z,
		     size_t panel_index);

void Compute_bvec_B6(double bvec[],
                     double x, double y, double z,
                     size_t panel_index);

void Compute_Coeff_B2(double x, double y, double z,
                     size_t panel_index,
                     const Kernel& kernel);

void Compute_Coeff_LM(double x, double y, double z,
                     size_t panel_index,
                     const Kernel& kernel);

void Compute_Coeff_B6(double x, double y, double z,
                     size_t panel_index,
                     const Kernel& kernel);

vec_4d assemble_velocity(const double vec[], int panel_index);

vec_4d Call_Tricubic_B2(double px, double py, double pz,
			int panel_index,
			const Kernel& kernel);

vec_4d Call_Tricubic_LM(double px, double py, double pz,
			int panel_index,
                        const Kernel& kernel);

vec_4d Call_Tricubic_B6(double px, double py, double pz,
                        int panel_index,
                        const Kernel& kernel);

void BinvMultiply_Transpose_B2(const double b[], double a[]);

void BinvMultiply_B2(const double b[], double a[]);

void BinvMultiply_Transpose_LM(const double b[], double a[]);

void BinvMultiply_LM(const double b[], double a[]);

void BinvMultiply_Transpose_B6(const double b[], double a[]);

void BinvMultiply_B6(const double b[], double a[]);

