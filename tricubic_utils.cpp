#include <fstream>

#include "tricubic_utils.h"

//***********************************************************************//

void Panel_Moment_Tricubic(size_t panel_index,
			   double *lambda,
			   struct xyz &particles,
			   double m[], 
                           double mx[],
                           double my[],
                           double mz[])
{
  double xmin = tree[panel_index].xinterval[0];
  double ymin = tree[panel_index].yinterval[0];
  double zmin = tree[panel_index].zinterval[0];
 
  double rdx = tree[panel_index].rxl; 
  double rdy = tree[panel_index].ryl; 
  double rdz = tree[panel_index].rzl;

  double tp0 = tree[panel_index].members[0];
  double tp1 = tree[panel_index].members[1];
  
  for (int i = 0; i < Pflat; i++){
      m[i] = 0.0; 
      mx[i]=0.0; my[i]=0.0; mz[i]=0.0;
  }

  for (size_t tp_j = tp0; tp_j <= tp1; tp_j++)
    {
      double x = (particles.x[tp_j] - xmin) * rdx;
      double y = (particles.y[tp_j] - ymin) * rdy;
      double z = (particles.z[tp_j] - zmin) * rdz;
      
      double si = lambda[tp_j]; 
      for (int i = 0; i < P + 1; i++)
	{
	  double sij = si; 
	  for (int j = 0; j < P + 1; j++)
	    {
              int i4j = i + 4*j;
	      
              double s = sij;
	      for (int k = 0; k < P + 1; k++)
		{
		  int kk = i4j + 16*k;
		  
		  m[kk] += s;

                  s *= z;
		}
                sij *= y;
	    }
            si *= x;
        }
    }

      for (int i = 1; i < P + 1; i++)
        {
          double di = static_cast<double>(i);

          int fourI = 4*i; int sixtI = 16*i;
          int fII = 4*(i-1); int sII = 16*(i-1);

          for (int j = 0; j < P + 1; j++)
            { 
              int i4j = i + 4*j;
              int j4i = j + fourI;

              int fourJ = 4*j;
              int fJsI  = fourJ + sixtI;
              int jfII  = j     + fII;
              int fJsII = fourJ + sII;

              for (int k = 0; k < P + 1; k++)
                { 
                  int sk = 16*k;
                  int ii = i4j + sk;
                  int jj = j4i + sk;
                  int kk = k   + fJsI;

                  int njj = jfII + sk;
                  int nkk = k    + fJsII;
 
                  mx[ii] += di*m[ii-1];
                  my[jj] += di*m[njj];
                  mz[kk] += di*m[nkk];       
                  
                }
            }
        }

}

//***********************************************************************//

void getClusterCorners_B2(size_t panel_index,
		       double &dx, double &dy, double &dz,
		       double ccx[], double ccy[], double ccz[])
{
  double xmin = tree[panel_index].xinterval[0];   
  double xmax = tree[panel_index].xinterval[1];   
  double ymin = tree[panel_index].yinterval[0];   
  double ymax = tree[panel_index].yinterval[1];   
  double zmin = tree[panel_index].zinterval[0];   
  double zmax = tree[panel_index].zinterval[1];   

  dx = xmax-xmin;
  dy = ymax-ymin;
  dz = zmax-zmin;

  ccx[0]=xmin; ccy[0]=ymin; ccz[0]=zmin;
  ccx[1]=xmax; ccy[1]=ymin; ccz[1]=zmin; 
  ccx[2]=xmin; ccy[2]=ymax; ccz[2]=zmin;
  ccx[3]=xmax; ccy[3]=ymax; ccz[3]=zmin;
  ccx[4]=xmin; ccy[4]=ymin; ccz[4]=zmax;
  ccx[5]=xmax; ccy[5]=ymin; ccz[5]=zmax; 
  ccx[6]=xmin; ccy[6]=ymax; ccz[6]=zmax;
  ccx[7]=xmax; ccy[7]=ymax; ccz[7]=zmax;


  ccx[8]=xmin+dx/4.0;  ccy[8]=ymin+dy/4.0; ccz[8]=zmin+dz/4.0;
  
  ccx[9] =ccx[8]+dx/2.0; ccy[9] =ccy[8];        ccz[9] =ccz[8];
  ccx[10]=ccx[8];        ccy[10]=ccy[8]+dy/2.0; ccz[10]=ccz[8];
  ccx[11]=ccx[9];        ccy[11]=ccy[10];       ccz[11]=ccz[8];
  ccx[12]=ccx[8];        ccy[12]=ccy[8];        ccz[12]=ccz[8]+dz/2.0;
  ccx[13]=ccx[9];        ccy[13]=ccy[8];        ccz[13]=ccz[12];
  ccx[14]=ccx[8];        ccy[14]=ccy[10];       ccz[14]=ccz[12];
  ccx[15]=ccx[9];        ccy[15]=ccy[10];       ccz[15]=ccz[12];

}

//************************************************************************//

void getClusterCorners_LM(size_t panel_index,
			  double &dx, double &dy, double &dz,
			  double ccx[], double ccy[], double ccz[])
{
  double xmin = tree[panel_index].xinterval[0];   
  double xmax = tree[panel_index].xinterval[1];   
  double ymin = tree[panel_index].yinterval[0];   
  double ymax = tree[panel_index].yinterval[1];   
  double zmin = tree[panel_index].zinterval[0];   
  double zmax = tree[panel_index].zinterval[1];   

  dx = xmax-xmin;
  dy = ymax-ymin;
  dz = zmax-zmin;

  ccx[0]=xmin; ccy[0]=ymin; ccz[0]=zmin;
  ccx[1]=xmax; ccy[1]=ymin; ccz[1]=zmin; 
  ccx[2]=xmin; ccy[2]=ymax; ccz[2]=zmin;
  ccx[3]=xmax; ccy[3]=ymax; ccz[3]=zmin;
  ccx[4]=xmin; ccy[4]=ymin; ccz[4]=zmax;
  ccx[5]=xmax; ccy[5]=ymin; ccz[5]=zmax; 
  ccx[6]=xmin; ccy[6]=ymax; ccz[6]=zmax;
  ccx[7]=xmax; ccy[7]=ymax; ccz[7]=zmax;
}

//************************************************************************//

void getClusterPoints_B6(size_t panel_index,
                          double &dx, double &dy, double &dz,
                          double ccx[], double ccy[], double ccz[])
{
  double xmin = tree[panel_index].xinterval[0];
  double xmax = tree[panel_index].xinterval[1];
  double ymin = tree[panel_index].yinterval[0];
  double ymax = tree[panel_index].yinterval[1];
  double zmin = tree[panel_index].zinterval[0];
  double zmax = tree[panel_index].zinterval[1];

  dx = xmax-xmin;
  dy = ymax-ymin;
  dz = zmax-zmin;

  ccx[0]=xmin; ccy[0]=ymin; ccz[0]=zmin;
  ccx[1]=xmax; ccy[1]=ymin; ccz[1]=zmin;
  ccx[2]=xmin; ccy[2]=ymax; ccz[2]=zmin;
  ccx[3]=xmax; ccy[3]=ymax; ccz[3]=zmin;
  ccx[4]=xmin; ccy[4]=ymin; ccz[4]=zmax;
  ccx[5]=xmax; ccy[5]=ymin; ccz[5]=zmax;
  ccx[6]=xmin; ccy[6]=ymax; ccz[6]=zmax;
  ccx[7]=xmax; ccy[7]=ymax; ccz[7]=zmax;

  ccx[8]=xmin + 0.5*dx; ccy[8]=ymin + 0.5*dy; ccz[8]=zmin + 0.5*dz; 
}

//*******************************************************************//

void scaleDerivatives_B2(double b[], double dx, double dy, double dz)
{
  b[1] *= dx;
  b[2] *= dy;
  b[3] *= dz;
}

//*******************************************************************//

void scaleDerivatives_LM(double b[], double dx, double dy, double dz)
{
  b[1] *= dx;
  b[2] *= dy;
  b[3] *= dz;
  b[4] *= dx*dy;
  b[5] *= dx*dz;
  b[6] *= dy*dz;
  b[7] *= dx*dy*dz;
}

//*******************************************************************//

void Compute_bvec_B2(double bvec[],
		     double x, double y, double z,
		     size_t panel_index,
		     const Kernel& kernel)
{
  double dx, dy, dz;
  double ccx[16], ccy[16], ccz[16]; // Corners of cluster
  getClusterCorners_B2(panel_index, dx, dy, dz, ccx, ccy, ccz);
  
  // loop through corners of cluster
  for (int i=0; i<16; i++) 
    {
      double btemp[8];
      // evaluate kernel and derivatives at corner i
      //kerEval_B2(x-ccx[i], y-ccy[i], z-ccz[i], btemp);
      kernel.eval_bvec(x-ccx[i], y-ccy[i], z-ccz[i], btemp);
      scaleDerivatives_B2(btemp, dx, dy, dz);

      // stack computed values into b vector
      // btemp = [f, dfdx, dfdy, dfdz]
      for (int j=0; j<4; j++)
	{
	  int jj = 16*j + i;
	  bvec[jj] = btemp[j]; ///27.0;
	}
    }
}

//***********************************************************************//

void Compute_bvec_LM(double bvec[],
		     double x, double y, double z,
		     size_t panel_index,
                     const Kernel& kernel)
{
  double dx, dy, dz;
  double ccx[8], ccy[8], ccz[8]; // Corners of cluster
  getClusterCorners_LM(panel_index, dx, dy, dz, ccx, ccy, ccz);
  
  // loop through corners of cluster
  for (int i=0; i<8; i++) 
    {
      double btemp[8];
      // evaluate kernel and derivatives at corner i
     // kerEval_LM(x-ccx[i], y-ccy[i], z-ccz[i], btemp);
      kernel.eval_bvec(x-ccx[i], y-ccy[i], z-ccz[i], btemp);
      scaleDerivatives_LM(btemp, dx, dy, dz);

      // stack computed values into b vector
      // btemp = [f, dfdx, dfdy, dfdz, d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz]
      for (int j=0; j<8; j++)
	{
	  int jj = 8*j + i;
	  bvec[jj] = btemp[j];
	}
    }
}

//***********************************************************************//

void Compute_bvec_B6(double bvec[],
                     double x, double y, double z,
                     size_t panel_index,
                     const Kernel& kernel)
{
  double dx, dy, dz;
  double ccx[9], ccy[9], ccz[9]; // Corners of cluster
  getClusterPoints_B6(panel_index, dx, dy, dz, ccx, ccy, ccz);

  // loop through corners of cluster
  double tempb;
  for (int i=0; i<9; i++)
    {
      double btemp[8];
      // evaluate kernel and derivatives at corner i
     // kerEval_LM(x-ccx[i], y-ccy[i], z-ccz[i], btemp);
      kernel.eval_bvec(x-ccx[i], y-ccy[i], z-ccz[i], btemp);
      scaleDerivatives_LM(btemp, dx, dy, dz);

      // stack computed values into b vector
      // btemp = [f, dfdx, dfdy, dfdz, d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz]
      for (int j=0; j<7; j++)
        {
          int jj = 9*j + i;
          bvec[jj] = btemp[j];
        }
     tempb = btemp[7];
    }
    bvec[63] = tempb;
}

//*******************************************************************//

void Compute_Coeff_B2(double x, double y, double z,
                     size_t panel_index,
                     const Kernel& kernel)
{
  double dx, dy, dz;
  double ccx[16], ccy[16], ccz[16]; // Corners of cluster
  getClusterCorners_B2(panel_index, dx, dy, dz, ccx, ccy, ccz);
  
  // loop through corners of cluster
  for (int i=0; i<16; i++) 
    {
      double btemp[8];
      // evaluate kernel and derivatives at corner i
      kernel.eval_bvec(x-ccx[i], y-ccy[i], z-ccz[i], btemp);
      scaleDerivatives_B2(btemp, dx, dy, dz);
     
      // sum up b-vectors
      for (int j=0; j<4; j++)
        {
          int jj = 16*j + i;
          tree[panel_index].moments[jj] += sfweight * btemp[j]; ///27.0;
        }
    }
}

//***********************************************************************//

void Compute_Coeff_LM(double x, double y, double z,
                     size_t panel_index,
                     const Kernel& kernel)
{
  double dx, dy, dz;
  double ccx[8], ccy[8], ccz[8]; // Corners of cluster
  getClusterCorners_LM(panel_index, dx, dy, dz, ccx, ccy, ccz);
  
  // loop through corners of cluster
  for (int i=0; i<8; i++) 
    {
      double btemp[8];
      // evaluate kernel and derivatives at corner i
      //kerEval_LM(x-ccx[i], y-ccy[i], z-ccz[i], btemp);
      kernel.eval_bvec(x-ccx[i], y-ccy[i], z-ccz[i], btemp);
      scaleDerivatives_LM(btemp, dx, dy, dz);

      // stack computed values into b vector
      // btemp = [f, dfdx, dfdy, dfdz, d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz]
      for (int j=0; j<8; j++)
        {
          int jj = 8*j + i;
          tree[panel_index].moments[jj] += sweight * btemp[j];
        }
    }
}

//***********************************************************************//

void Compute_Coeff_B6(double x, double y, double z,
                     size_t panel_index,
                     const Kernel& kernel)
{
  double dx, dy, dz;
  double ccx[9], ccy[9], ccz[9]; // Corners of cluster
  getClusterPoints_B6(panel_index, dx, dy, dz, ccx, ccy, ccz);

  // loop through corners of cluster
  double tempb;
  for (int i=0; i<9; i++)
    {
      double btemp[8];
      // evaluate kernel and derivatives at corner i
      //kerEval_LM(x-ccx[i], y-ccy[i], z-ccz[i], btemp);
      kernel.eval_bvec(x-ccx[i], y-ccy[i], z-ccz[i], btemp);
      scaleDerivatives_LM(btemp, dx, dy, dz);

      // stack computed values into b vector
      // btemp = [f, dfdx, dfdy, dfdz, d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz]
      for (int j=0; j<7; j++)
        {
          int jj = 9*j + i;
          tree[panel_index].moments[jj] += sweight * btemp[j];
        }
       tempb = btemp[7];
    }
    tree[panel_index].moments[63] += sweight * tempb;
}

//************************************************************************//

vec_4d assemble_velocity(const double vec[], int panel_index)
{
  vec_4d velocity;
  velocity.val[0] = 0.0;
  velocity.val[1] = 0.0;
  velocity.val[2] = 0.0;
  velocity.val[3] = 0.0;
  
  for (int kk = 0; kk < 64; kk++)
    {
      double veckk = vec[kk];
      velocity.val[0] += veckk * tree[panel_index].moments[kk];
      velocity.val[1] += veckk * tree[panel_index].dxmoments[kk];
      velocity.val[2] += veckk * tree[panel_index].dymoments[kk];
      velocity.val[3] += veckk * tree[panel_index].dzmoments[kk];
    }
   
  return velocity;
}

//*******************************************************************//

vec_4d Call_Tricubic_B2(double px, double py, double pz,
			int panel_index,
			const Kernel& kernel)
{
  double bvec[64];
  init_vec(bvec);
  
  Compute_bvec_B2(bvec, px, py, pz, panel_index, kernel);
  
  return assemble_velocity(bvec, panel_index);
}

//*******************************************************************//

vec_4d Call_Tricubic_LM(double px, double py, double pz,
				      int panel_index,
                        const Kernel& kernel)
{
    double bvec[64];
    init_vec(bvec);
      
    Compute_bvec_LM(bvec, px, py, pz, panel_index, kernel);

    return assemble_velocity(bvec, panel_index);
}

//*******************************************************************//

vec_4d Call_Tricubic_B6(double px, double py, double pz,
                                      int panel_index,
                        const Kernel& kernel)
{
    double bvec[64];
    init_vec(bvec);

    Compute_bvec_B6(bvec, px, py, pz, panel_index, kernel);

    return assemble_velocity(bvec, panel_index);
}

//*************************************************************************//

// performs explicit multiplication by Transpose(Binv) matrix

void BinvMultiply_Transpose_B2(const double b[], double a[])
{
  double bfac = 0.037037037037037037037037;
  
  a[0]  = bfac*(
              - 93184.0* b[42]
              + 49856.0*(b[26] + b[38] + b[41])
              + 48128.0*(b[43] + b[46] + b[58])
              - 26353.0*(b[22] + b[25] + b[37])
              - 25840.0*(b[27] + b[30] + b[39] + b[45] + b[54] + b[57])
              - 24832.0*(b[47] + b[59] + b[62])
              + 13718.0*(b[21] + b[23] + b[29] + b[53])
              + 13376.0*(b[31] + b[55] + b[61])
              + 12800.0* b[63]
              -  4800.0*(b[10] + b[34] + b[40])
              +  2544.0*(b[11] + b[14] + b[35] + b[44] + b[50] + b[56])
              +  2337.0*(b[6]  + b[9]  + b[18] + b[24] + b[33] + b[36])
              -  1344.0*(b[15] + b[51] + b[60])
              -  1254.0*(b[7]  + b[13] + b[19] + b[28] + b[49] + b[52])
              -  1083.0*(b[5]  + b[17] + b[20])
              -    81.0*(b[2]  + b[8]  + b[32])
              +    54.0*(b[3]  + b[12] + b[48])
              +    27.0 *b[0] );
      a[1]  = bfac*(
                51200.0* b[42]
              - 48128.0* b[43]
              - 27664.0*(b[26] + b[38])
              - 26368.0*(b[46] + b[58])
              + 25840.0*(b[27] + b[39])
              + 24832.0*(b[47] + b[59])
              + 14801.0* b[22]
              + 14288.0*(b[30] + b[54])
              - 13718.0* b[23]
              + 13568.0* b[62]
              - 13376.0*(b[31] + b[55])
              - 12800.0* b[63]
              -  7872.0* b[41]
              +  4161.0*(b[25] + b[37])
              +  4080.0*(b[45] + b[57])
              +  2832.0*(b[10] + b[34])
              -  2544.0*(b[11] + b[35])
              -  2166.0*(b[21] + b[29] + b[53])
              -  2112.0* b[61]
              -  1488.0*(b[14] + b[50])
              -  1425.0*(b[6]  + b[18])
              +  1344.0*(b[15] + b[51])
              +  1254.0*(b[7]  + b[19])
              -   369.0*(b[9]  + b[33])
              +   198.0*(b[13] + b[49])
              +   171.0*(b[5]  + b[17])
              +    81.0* b[2]
              -    54.0* b[3] );
      a[2]  = bfac*(
                51200.0* b[42]
              - 48128.0* b[46]
              - 27664.0*(b[26] + b[41])
              - 26368.0*(b[43] + b[58])
              + 25840.0*(b[30] + b[45])
              + 24832.0*(b[47] + b[62])
              + 14801.0* b[25]
              + 14288.0*(b[27] + b[57])
              - 13718.0* b[29]
              + 13568.0* b[59]
              - 13376.0*(b[31] + b[61])
              - 12800.0* b[63]
              -  7872.0* b[38]
              +  4161.0*(b[22] + b[37])
              +  4080.0*(b[39] + b[54])
              +  2832.0*(b[10] + b[40])
              -  2544.0*(b[14] + b[44])
              -  2166.0*(b[21] + b[23] + b[53])
              -  2112.0* b[55]
              -  1488.0*(b[11] + b[56])
              -  1425.0*(b[9]  + b[24])
              +  1344.0*(b[15] + b[60])
              +  1254.0*(b[13] + b[28])
              -   369.0*(b[6]  + b[36])
              +   198.0*(b[7]  + b[52])
              +   171.0*(b[5]  + b[20])
              +    81.0* b[8]
              -    54.0* b[12] ); 
      a[3]  = bfac*(
              - 27904.0* b[42]
              + 26368.0*(b[43] + b[46])
              - 24832.0* b[47]
              + 15200.0* b[26]
              + 14336.0* b[58]
              - 14288.0*(b[27] + b[30])
              - 13568.0*(b[59] + b[62])
              + 13376.0* b[31]
              + 12800.0* b[63]
              +  4368.0*(b[38] + b[41])
              -  4080.0*(b[39] + b[45])
              -  2337.0*(b[22] + b[25])
              -  2256.0*(b[54] + b[57])
              +  2166.0*(b[23] + b[29])
              +  2112.0*(b[55] + b[61])
              -  1632.0* b[10]
              +  1488.0*(b[11] + b[14])
              -  1344.0* b[15]
              -   657.0* b[37]
              +   342.0*(b[21] + b[53])
              +   225.0*(b[6]  + b[9])
              -   198.0*(b[7]  + b[13])
              -    27.0* b[5] );
      a[4]  = bfac*(
                51200.0* b[42]
              - 48128.0* b[58]
              - 27664.0*(b[38] + b[41])
              - 26368.0*(b[43] + b[46])
              + 25840.0*(b[54] + b[57])
              + 24832.0*(b[59] + b[62])
              + 14801.0* b[37]
              + 14288.0*(b[39] + b[45])
              - 13718.0* b[53]
              + 13568.0* b[47]
              - 13376.0*(b[55] + b[61])
              - 12800.0* b[63]
              -  7872.0* b[26]
              +  4161.0*(b[22] + b[25])
              +  4080.0*(b[27] + b[30])
              +  2832.0*(b[34] + b[40])
              -  2544.0*(b[50] + b[56])
              -  2166.0*(b[21] + b[23] + b[29])
              -  2112.0* b[31]
              -  1488.0*(b[35] + b[44])
              -  1425.0*(b[33] + b[36])
              +  1344.0*(b[51] + b[60])
              +  1254.0*(b[49] + b[52])
              -   369.0*(b[18] + b[24])
              +   198.0*(b[19] + b[28])
              +   171.0*(b[17] + b[20])
              +    81.0* b[32]
              -    54.0* b[48] );
 
      a[5]  = bfac*(
              - 27904.0* b[42]
              + 26368.0*(b[43] + b[58])
              - 24832.0* b[59]
              + 15200.0* b[38]
              + 14336.0* b[46]
              - 14288.0*(b[39] + b[54])
              - 13568.0*(b[47] + b[62])
              + 13376.0* b[55]
              + 12800.0* b[63]
              +  4368.0*(b[26] + b[41])
              -  4080.0*(b[27] + b[57])
              -  2337.0*(b[22] + b[37])
              -  2256.0*(b[30] + b[45])
              +  2166.0*(b[23] + b[53])
              +  2112.0*(b[31] + b[61])
              -  1632.0* b[34]
              +  1488.0*(b[35] + b[50])
              -  1344.0* b[51]
              -   657.0* b[25]
              +   342.0*(b[21] + b[29])
              +   225.0*(b[18] + b[33])
              -   198.0*(b[19] + b[49])
              -    27.0* b[17] ); 
      a[6]  = bfac*(
              - 27904.0* b[42]
              + 26368.0*(b[46] + b[58])
              - 24832.0* b[62]
              + 15200.0* b[41]
              + 14336.0* b[43]
              - 14288.0*(b[45] + b[57])
              - 13568.0*(b[47] + b[59])
              + 13376.0* b[61]
              + 12800.0* b[63]
              +  4368.0*(b[26] + b[38])
              -  4080.0*(b[30] + b[54])
              -  2337.0*(b[25] + b[37])
              -  2256.0*(b[27] + b[39])
              +  2166.0*(b[29] + b[53])
              +  2112.0*(b[31] + b[55])
              -  1632.0* b[40]
              +  1488.0*(b[44] + b[56])
              -  1344.0* b[60]
              -   657.0* b[22]
              +   342.0*(b[21] + b[23])
              +   225.0*(b[24] + b[36])
              -   198.0*(b[28] + b[52])
              -    27.0* b[20] ); 
      a[7]  = bfac*(
                15104.0* b[42]
              - 14336.0*(b[43] + b[46] + b[58])
              + 13568.0*(b[47] + b[59] + b[62])
              - 12800.0* b[63]
              -  2400.0*(b[26] + b[38] + b[41])
              +  2256.0*(b[27] + b[30] + b[39] + b[45] + b[54] + b[57])
              -  2112.0*(b[31] + b[55] + b[61])
              +   369.0*(b[22] + b[25] + b[37])
              -   342.0*(b[23] + b[29] + b[53])
              -    54.0* b[21] );
      
      a[8]  = bfac*(
              -326144.0* b[42]
              +207872.0*(b[43] + b[46] + b[58])
              -131072.0*(b[47] + b[59] + b[62])
              +127680.0*(b[26] + b[38] + b[41])
              - 82176.0*(b[27] + b[30] + b[39] + b[45] + b[54] + b[57])
              + 81920.0* b[63]
              + 52224.0*(b[31] + b[55] + b[61])
              - 49536.0*(b[22] + b[25] + b[37])
              + 32256.0*(b[23] + b[29] + b[53])
              + 19008.0* b[21] );
      a[9]  = bfac*(
               297472.0* b[42]
              -207872.0* b[43]
              -185344.0*(b[46] + b[58])
              +131072.0*(b[47] + b[59])
              -118848.0*(b[26] + b[38])
              +114688.0* b[62]
              - 99008.0* b[41]
              + 82176.0*(b[27] + b[39])
              - 81920.0* b[63]
              + 74496.0*(b[30] + b[54])
              + 59648.0*(b[45] + b[57])
              - 52224.0*(b[31] + b[55])
              + 47232.0* b[22]
              + 40704.0*(b[25] + b[37])
              - 35840.0* b[61]
              - 32256.0* b[23]
              - 24576.0*(b[29] + b[53])
              - 16704.0* b[21]
              +  9408.0* b[40]
              -  5376.0*(b[44] + b[56])
              -  4032.0*(b[24] + b[36])
              +  3072.0* b[60]
              +  2304.0*(b[28] + b[52])
              +  1728.0* b[20] );
      a[10] = bfac*(
               297472.0* b[42]
              -207872.0* b[46]
              -185344.0*(b[43] + b[58])
              +131072.0*(b[47] + b[62])
              -118848.0*(b[26] + b[41])
              +114688.0* b[59]
              - 99008.0* b[38]
              + 82176.0*(b[30] + b[45])
              - 81920.0* b[63]
              + 74496.0*(b[27] + b[57])
              + 59648.0*(b[39] + b[54])
              - 52224.0*(b[31] + b[61])
              + 47232.0* b[25]
              + 40704.0*(b[22] + b[37])
              - 35840.0* b[55]
              - 32256.0* b[29]
              - 24576.0*(b[23] + b[53])
              - 16704.0* b[21]
              +  9408.0* b[34]
              -  5376.0*(b[35] + b[50])
              -  4032.0*(b[18] + b[33])
              +  3072.0* b[51]
              +  2304.0*(b[19] + b[49])
              +  1728.0* b[17] ); 
      a[11] = bfac*(
              -258560.0* b[42]
              +185344.0*(b[43] + b[46])
              +158720.0* b[58]
              -131072.0* b[47]
              -114688.0*(b[59] + b[62])
              +104640.0* b[26]
              + 81920.0* b[63]
              + 79936.0*(b[38] + b[41])
              - 74496.0*(b[27] + b[30])
              - 59648.0*(b[39] + b[45])
              + 52224.0* b[31]
              - 47872.0*(b[54] + b[57])
              + 35840.0*(b[55] + b[61])
              - 33024.0*(b[22] + b[25])
              + 24576.0*(b[23] + b[29])
              - 21632.0* b[37]
              + 12800.0* b[53]
              +  9024.0* b[21]
              -  6720.0*(b[34] + b[40])
              +  5376.0*(b[35] + b[44])
              +  3840.0*(b[50] + b[56])
              -  3072.0*(b[51] + b[60])
              +  2880.0*(b[18] + b[24])
              -  2304.0*(b[19] + b[28])
              +  1344.0*(b[33] + b[36])
              -   768.0*(b[49] + b[52])
              -   576.0*(b[17] + b[20]) ); 
      a[12] = bfac*(
               297472.0* b[42]
              -207872.0* b[58]
              -185344.0*(b[43] + b[46])
              +131072.0*(b[59] + b[62])
              -118848.0*(b[38] + b[41])
              +114688.0* b[47]
              - 99008.0* b[26]
              + 82176.0*(b[54] + b[57])
              - 81920.0* b[63]
              + 74496.0*(b[39] + b[45])
              + 59648.0*(b[27] + b[30])
              - 52224.0*(b[55] + b[61])
              + 47232.0* b[37]
              + 40704.0*(b[22] + b[25])
              - 35840.0* b[31]
              - 32256.0* b[53]
              - 24576.0*(b[23] + b[29])
              - 16704.0* b[21]
              +  9408.0* b[10]
              -  5376.0*(b[11] + b[14])
              -  4032.0*(b[6]  + b[9] )
              +  3072.0* b[15]
              +  2304.0*(b[7]  + b[13])
              +  1728.0* b[5] ); 
      a[13] = bfac*(
              -258560.0* b[42]
              +185344.0*(b[43] + b[58])
              +158720.0* b[46]
              -131072.0* b[59]
              -114688.0*(b[47] + b[62])
              +104640.0* b[38]
              + 81920.0* b[63]
              + 79936.0*(b[26] + b[41])
              - 74496.0*(b[39] + b[54])
              - 59648.0*(b[27] + b[57])
              + 52224.0* b[55]
              - 47872.0*(b[30] + b[45])
              + 35840.0*(b[31] + b[61])
              - 33024.0*(b[22] + b[37])
              + 24576.0*(b[23] + b[53])
              - 21632.0* b[25]
              + 12800.0* b[29]
              +  9024.0* b[21]
              -  6720.0*(b[10] + b[40])
              +  5376.0*(b[11] + b[56])
              +  3840.0*(b[14] + b[44])
              -  3072.0*(b[15] + b[60])
              +  2880.0*(b[6]  + b[36])
              -  2304.0*(b[7]  + b[52])
              +  1344.0*(b[9]  + b[24])
              -   768.0*(b[13] + b[28])
              -   576.0*(b[5]  + b[20]) ); 
      a[14] = bfac*(
              -258560.0* b[42]
              +185344.0*(b[46] + b[58])
              +158720.0* b[43]
              -131072.0* b[62]
              -114688.0*(b[47] + b[59])
              +104640.0* b[41]
              + 81920.0* b[63]
              + 79936.0*(b[26] + b[38])
              - 74496.0*(b[45] + b[57])
              - 59648.0*(b[30] + b[54])
              + 52224.0* b[61]
              - 47872.0*(b[27] + b[39])
              + 35840.0*(b[31] + b[55])
              - 33024.0*(b[25] + b[37])
              + 24576.0*(b[29] + b[53])
              - 21632.0* b[22]
              + 12800.0* b[23]
              +  9024.0* b[21]
              -  6720.0*(b[10] + b[34])
              +  5376.0*(b[14] + b[50])
              +  3840.0*(b[11] + b[35])
              -  3072.0*(b[15] + b[51])
              +  2880.0*(b[9]  + b[33])
              -  2304.0*(b[13] + b[49])
              +  1344.0*(b[6]  + b[18])
              -   768.0*(b[7]  + b[19])
              -   576.0*(b[5]  + b[17]) ); 
      a[15] = bfac*(
               217600.0* b[42]
              -158720.0*(b[43] + b[46] + b[58])
              +114688.0*(b[47] + b[59] + b[62])
              - 81920.0* b[63]
              - 63680.0*(b[26] + b[38] + b[41])
              + 47872.0*(b[27] + b[30] + b[39] + b[45] + b[54] + b[57])
              - 35840.0*(b[31] + b[55] + b[61])
              + 16768.0*(b[22] + b[25] + b[37])
              - 12800.0*(b[23] + b[29] + b[53])
              +  4800.0*(b[10] + b[34] + b[40])
              -  4160.0* b[21] 
              -  3840.0*(b[11] + b[14] + b[35] + b[44] + b[50] + b[56])
              +  3072.0*(b[15] + b[51] + b[60])
              -   960.0*(b[6]  + b[9]  + b[18] + b[24] + b[33] + b[36])
              +   768.0*(b[7]  + b[13] + b[19] + b[28] + b[49] + b[52])
              +   192.0*(b[5]  + b[17] + b[20]) );  
      ////////////////////

      a[16] = bfac*(
              -  6144.0* b[42]
              +  3648.0*(b[26] + b[38])
              +  3072.0*(b[41] + b[43] + b[46] + b[58])
              -  2166.0* b[22]
              -  1824.0*(b[25] + b[27] + b[30] + b[37] + b[39] + b[54])
              -  1536.0*(b[45] + b[47] + b[57] + b[59] + b[62])
              +  1083.0*(b[21] + b[23])
              +   912.0*(b[29] + b[31] + b[53] + b[55])
              +   768.0*(b[61] + b[63])
              -   576.0*(b[10] + b[34])
              +   342.0*(b[6]  + b[18])
              +   288.0*(b[9]  + b[11] + b[14] + b[33] + b[35] + b[50])
              -   171.0*(b[5]  + b[7]  + b[17] + b[19])
              -   144.0*(b[13] + b[15] + b[49] + b[51])
              -    54.0* b[2]
              +    27.0*(b[1]  + b[3]) );
      a[17] = bfac*(
              -  3072.0*(b[42] - b[43])
              +  1824.0*(b[26] - b[27] + b[38] - b[39])
              +  1536.0*(b[46] - b[47] + b[58] - b[59])
              -  1083.0*(b[22] - b[23])
              -   912.0*(b[30] - b[31] + b[54] - b[55])
              -   768.0*(b[62] - b[63])
              -   288.0*(b[10] - b[11] + b[34] - b[35])
              +   171.0*(b[6]  - b[7]  + b[18] - b[19])
              +   144.0*(b[14] - b[15] + b[50] - b[51])
              -    27.0*(b[2]  - b[3]) );
      a[18] = bfac*(
                 3072.0*(b[42] - b[46])
              -  1824.0*(b[26] - b[30])
              -  1536.0*(b[41] + b[43] - b[45] - b[47] + b[58] - b[62])
              +   912.0*(b[25] + b[27] - b[29] - b[31])
              +   768.0*(b[57] + b[59] - b[61] - b[63])
              -   576.0* b[38]
              +   342.0* b[22]
              +   288.0*(b[10] - b[14] + b[37] + b[39] + b[54])
              -   171.0*(b[21] + b[23])
              -   144.0*(b[9]  + b[11] - b[13] - b[15] + b[53] + b[55])
              -    54.0* b[6]
              +    27.0*(b[5]  + b[7]) );
      a[19] = bfac*(
                 1536.0*(b[42] - b[43] - b[46] + b[47])
              -   912.0*(b[26] - b[27] - b[30] + b[31])
              -   768.0*(b[58] - b[59] - b[62] + b[63])
              -   288.0*(b[38] - b[39])
              +   171.0*(b[22] - b[23])
              +   144.0*(b[10] - b[11] - b[14] + b[15] + b[54] - b[55])
              -    27.0*(b[6]  - b[7]) );
      a[20] = bfac*(
                 3072.0*(b[42] - b[58])
              -  1824.0*(b[38] - b[54])
              -  1536.0*(b[41] + b[43] + b[46] - b[57] - b[59] - b[62])
              +   912.0*(b[37] + b[39] - b[53] - b[55])
              +   768.0*(b[45] + b[47] - b[61] - b[63])
              -   576.0* b[26]
              +   342.0* b[22]
              +   288.0*(b[25] + b[27] + b[30] + b[34] - b[50])
              -   171.0*(b[21] + b[23])
              -   144.0*(b[29] + b[31] + b[33] + b[35] - b[49] - b[51])
              -    54.0* b[18]
              +    27.0*(b[17] + b[19]) ); 
      a[21] = bfac*(
                 1536.0*(b[42] - b[43] - b[58] + b[59])
              -   912.0*(b[38] - b[39] - b[54] + b[55])
              -   768.0*(b[46] - b[47] - b[62] + b[63])
              -   288.0*(b[26] - b[27])
              +   171.0*(b[22] - b[23])
              +   144.0*(b[30] - b[31] + b[34] - b[35] - b[50] + b[51])
              -    27.0*(b[18] - b[19]) ); 
      a[22] = bfac*(
              -  1536.0*(b[42] - b[46] - b[58] + b[62])
              +   768.0*(b[41] + b[43] - b[45] - b[47] - b[57] - b[59] + b[61] + b[63])
              +   288.0*(b[26] - b[30] + b[38] - b[54])
              -   144.0*(b[25] + b[27] - b[29] - b[31] + b[37] + b[39] - b[53] - b[55])
              -    54.0* b[22]
              +    27.0*(b[21] + b[23]) ); 
      a[23] = bfac*(
              -   768.0*(b[42] - b[43] - b[46] + b[47] - b[58] + b[59] + b[62] - b[63])
              +   144.0*(b[26] - b[27] - b[30] + b[31] + b[38] - b[39] - b[54] + b[55])
              -    27.0*(b[22] - b[23]) ); 
      a[24] = bfac*(
              - 65856.0* b[42]
              + 37632.0*(b[43] + b[46] + b[58])
              + 35280.0* b[41]
              + 28224.0*(b[26] + b[38])
              - 21504.0*(b[47] + b[59] + b[62])
              - 20160.0*(b[45] + b[57])
              - 16128.0*(b[27] + b[30] + b[39] + b[54])
              - 15120.0*(b[25] + b[37])
              + 12288.0* b[63]
              - 12096.0* b[22]
              + 11520.0* b[61]
              +  9216.0*(b[31] + b[55])
              +  8640.0*(b[29] + b[53])
              +  6912.0* b[23]
              +  6480.0* b[21]
              -  5292.0* b[40]
              +  3024.0*(b[44] + b[56])
              +  2268.0*(b[24] + b[36])
              -  1728.0* b[60]
              -  1296.0*(b[28] + b[52])
              -   972.0* b[20] );
      a[25] = bfac*(
              - 47040.0* b[42]
              + 37632.0* b[43]
              + 26880.0*(b[46] + b[58])
              - 21504.0*(b[47] + b[59])
              + 20160.0*(b[26] + b[38])
              + 16464.0* b[41]
              - 16128.0*(b[27] + b[39])
              - 15360.0* b[62]
              + 12288.0* b[63]
              - 11520.0*(b[30] + b[54])
              -  9408.0*(b[45] + b[57])
              +  9216.0*(b[31] + b[55])
              -  8640.0* b[22]
              -  7056.0*(b[25] + b[37])
              +  6912.0* b[23]
              +  5376.0* b[61]
              +  4032.0*(b[29] + b[53])
              +  3024.0* b[21]
              -  1764.0* b[40]
              +  1008.0*(b[44] + b[56])
              +   756.0*(b[24] + b[36])
              -   576.0* b[60]
              -   432.0*(b[28] + b[52])
              -   324.0* b[20] );
      a[26] = bfac*(
                47040.0* b[42]
              - 37632.0* b[46]
              - 26880.0*(b[43] + b[58])
              - 25200.0* b[41]
              + 21504.0*(b[47] + b[62])
              - 20160.0*(b[26] - b[45])
              + 16128.0* b[30]
              + 15360.0* b[59]
              + 14400.0* b[57]
              - 12288.0* b[63]
              + 11520.0*(b[27] - b[61])
              + 10800.0* b[25]
              -  9408.0* b[38]
              -  9216.0* b[31]
              -  8640.0* b[29]
              +  5376.0*(b[39] + b[54])
              +  5040.0* b[37]
              +  4032.0* b[22]
              +  3780.0* b[40]
              -  3072.0* b[55]
              -  3024.0* b[44]
              -  2880.0* b[53]
              -  2304.0* b[23]
              -  2160.0*(b[21] + b[56])
              +  1728.0* b[60]
              -  1620.0* b[24]
              +  1296.0* b[28]
              -   756.0* b[36]
              +   432.0* b[52]
              +   324.0* b[20] ); 
      a[27] = bfac*(
                33600.0* b[42]
              - 26880.0*(b[43] + b[46])
              + 21504.0* b[47]
              - 19200.0* b[58]
              + 15360.0*(b[59] + b[62])
              - 14400.0* b[26]
              - 12288.0* b[63]
              - 11760.0* b[41]
              + 11520.0*(b[27] + b[30])
              +  9408.0* b[45]
              -  9216.0* b[31]
              -  6720.0*(b[38] - b[57])
              +  5376.0*(b[39] - b[61])
              +  5040.0* b[25]
              -  4032.0* b[29]
              +  3840.0* b[54]
              -  3072.0* b[55]
              +  2880.0* b[22]
              +  2352.0* b[37]
              -  2304.0* b[23]
              -  1344.0* b[53]
              +  1260.0* b[40]
              -  1008.0*(b[21] + b[44])
              -   720.0* b[56]
              +   576.0* b[60]
              -   540.0* b[24]
              +   432.0* b[28]
              -   252.0* b[36]
              +   144.0* b[52]
              +   108.0* b[20] ); 
      a[28] = bfac*(
                47040.0* b[42]
              - 37632.0* b[58]
              - 26880.0*(b[43] + b[46])
              - 25200.0* b[41]
              + 21504.0*(b[59] + b[62])
              - 20160.0*(b[38] - b[57])
              + 16128.0* b[54]
              + 15360.0* b[47]
              + 14400.0* b[45]
              - 12288.0* b[63]
              + 11520.0*(b[39] - b[61])
              + 10800.0* b[37]
              -  9408.0* b[26]
              -  9216.0* b[55]
              -  8640.0* b[53]
              +  5376.0*(b[27] + b[30])
              +  5040.0* b[25]
              +  4032.0* b[22]
              +  3780.0* b[40]
              -  3072.0* b[31]
              -  3024.0* b[56]
              -  2880.0* b[29]
              -  2304.0* b[23]
              -  2160.0*(b[21] + b[44])
              +  1728.0* b[60]
              -  1620.0* b[36]
              +  1296.0* b[52]
              -   756.0* b[24]
              +   432.0* b[28]
              +   324.0* b[20] ); 
      a[29] = bfac*(
                33600.0* b[42]
              - 26880.0*(b[43] + b[58])
              + 21504.0* b[59]
              - 19200.0* b[46]
              + 15360.0*(b[47] + b[62])
              - 14400.0* b[38]
              - 12288.0* b[63]
              - 11760.0* b[41]
              + 11520.0*(b[39] + b[54])
              +  9408.0* b[57]
              -  9216.0* b[55]
              -  6720.0*(b[26] - b[45])
              +  5376.0*(b[27] - b[61])
              +  5040.0* b[37]
              -  4032.0* b[53]
              +  3840.0* b[30]
              -  3072.0* b[31]
              +  2880.0* b[22]
              +  2352.0* b[25]
              -  2304.0* b[23]
              -  1344.0* b[29]
              +  1260.0* b[40]
              -  1008.0*(b[21] + b[56])
              -   720.0* b[44]
              +   576.0* b[60]
              -   540.0* b[36]
              +   432.0* b[52]
              -   252.0* b[24]
              +   144.0* b[28]
              +   108.0* b[20] ); 
      a[30] = bfac*(
              - 33600.0* b[42]
              + 26880.0*(b[46] + b[58])
              - 21504.0* b[62]
              + 19200.0* b[43]
              + 18000.0* b[41]
              - 15360.0*(b[47] + b[59])
              - 14400.0*(b[45] + b[57])
              + 12288.0* b[63]
              + 11520.0* b[61]
              +  6720.0*(b[26] + b[38])
              -  5376.0*(b[30] + b[54])
              -  3840.0*(b[27] + b[39])
              -  3600.0*(b[25] + b[37])
              +  3072.0*(b[31] + b[55])
              +  2880.0*(b[29] + b[53])
              -  2700.0* b[40]
              +  2160.0*(b[44] + b[56])
              -  1728.0* b[60]
              -  1344.0* b[22]
              +   768.0* b[23]
              +   720.0* b[21]
              +   540.0*(b[24] + b[36])
              -   432.0*(b[28] + b[52])
              -   108.0* b[20] ); 
      a[31] = bfac*(
              - 24000.0* b[42]
              + 19200.0*(b[43] + b[46] + b[58])
              - 15360.0*(b[47] + b[59] + b[62])
              + 12288.0* b[63]
              +  8400.0* b[41]
              -  6720.0*(b[45] + b[57])
              +  5376.0* b[61]
              +  4800.0*(b[26] + b[38])
              -  3840.0*(b[27] + b[30] + b[39] + b[54])
              +  3072.0*(b[31] + b[55])
              -  1680.0*(b[25] + b[37])
              +  1344.0*(b[29] + b[53])
              -   960.0* b[22]
              -   900.0* b[40]
              +   768.0* b[23]
              +   720.0*(b[44] + b[56])
              -   576.0* b[60]
              +   336.0* b[21]
              +   180.0*(b[24] + b[36])
              -   144.0*(b[28] + b[52])
              -    36.0* b[20] ); 
      ////////////////////

      a[32] = bfac*(
              -  6144.0* b[42]
              +  3648.0*(b[26] + b[41])
              +  3072.0*(b[38] + b[43] + b[46] + b[58])
              -  2166.0* b[25]
              -  1824.0*(b[22] + b[27] + b[30] + b[37] + b[45] + b[57])
              -  1536.0*(b[39] + b[47] + b[54] + b[59] + b[62])
              +  1083.0*(b[21] + b[29])
              +   912.0*(b[23] + b[31] + b[53] + b[61])
              +   768.0*(b[55] + b[63])
              -   576.0*(b[10] + b[40])
              +   342.0*(b[9]  + b[24])
              +   288.0*(b[6]  + b[11] + b[14] + b[36] + b[44] + b[56])
              -   171.0*(b[5]  + b[13] + b[20] + b[28])
              -   144.0*(b[7]  + b[15] + b[52] + b[60])
              -    54.0* b[8]
              +    27.0*(b[4]  + b[12]) ); 
      a[33] = bfac*(
              +  3072.0*(b[42] - b[43])
              -  1824.0*(b[26] - b[27]) 
              -  1536.0*(b[38] - b[39] + b[46] - b[47] + b[58] - b[59])
              +   912.0*(b[22] - b[23] + b[30] - b[31])
              +   768.0*(b[54] - b[55] + b[62] - b[63])
              -   576.0* b[41]
              +   342.0* b[25]  
              +   288.0*(b[10] - b[11] + b[37] + b[45] + b[57])
              -   171.0*(b[21] + b[29])
              -   144.0*(b[6]  - b[7]  + b[14] - b[15] + b[53] + b[61])
              -    54.0* b[9]
              +    27.0*(b[5]  + b[13]) ); 
      a[34] = bfac*(
              -  3072.0*(b[42] - b[46])
              +  1824.0*(b[26] - b[30] + b[41] - b[45])
              +  1536.0*(b[43] - b[47] + b[58] - b[62])
              -  1083.0*(b[25] - b[29])
              -   912.0*(b[27] - b[31] + b[57] - b[61])
              -   768.0*(b[59] - b[63])
              -   288.0*(b[10] - b[14] + b[40] - b[44])
              +   171.0*(b[9]  - b[13] + b[24] - b[28])
              +   144.0*(b[11] - b[15] + b[56] - b[60])
              -    27.0*(b[8]  - b[12]) ); 
      a[35] = bfac*(
                 1536.0*(b[42] - b[43] - b[46] + b[47])
              -   912.0*(b[26] - b[27] - b[30] + b[31])
              -   768.0*(b[58] - b[59] - b[62] + b[63])
              -   288.0*(b[41] - b[45])
              +   171.0*(b[25] - b[29])
              +   144.0*(b[10] - b[11] - b[14] + b[15] + b[57] - b[61])
              -    27.0*(b[9]  - b[13]) ); 
      a[36] = bfac*(
              +  3072.0*(b[42] - b[58])
              -  1824.0*(b[41] - b[57])
              -  1536.0*(b[38] + b[43] + b[46] - b[54] - b[59] - b[62])
              +   912.0*(b[37] + b[45] - b[53] - b[61])
              +   768.0*(b[39] + b[47] - b[55] - b[63])
              -   576.0* b[26]
              +   342.0* b[25]
              +   288.0*(b[22] + b[27] + b[30] + b[40] - b[56])
              -   171.0*(b[21] + b[29])
              -   144.0*(b[23] + b[31] + b[36] + b[44] - b[52] - b[60])
              -    54.0* b[24]
              +    27.0*(b[20] + b[28]) ); 
      a[37] = bfac*(
              -  1536.0*(b[42] - b[43] - b[58] + b[59])
              +   768.0*(b[38] - b[39] + b[46] - b[47] - b[54] + b[55] - b[62] + b[63])
              +   288.0*(b[26] - b[27] + b[41] - b[57])
              -   144.0*(b[22] - b[23] + b[30] - b[31] + b[37] + b[45] - b[53] - b[61])
              -    54.0* b[25]
              +    27.0*(b[21] + b[29]) ); 
      a[38] = bfac*(
                 1536.0*(b[42] - b[46] - b[58] + b[62])
              -   912.0*(b[41] - b[45] - b[57] + b[61])
              -   768.0*(b[43] - b[47] - b[59] + b[63])
              -   288.0*(b[26] - b[30])
              +   171.0*(b[25] - b[29])
              +   144.0*(b[27] - b[31]  + b[40] - b[44] - b[56] + b[60])
              -    27.0*(b[24] - b[28]) ); 
      a[39] = bfac*(
              -  768.0*(b[42] - b[43] - b[46] + b[47] - b[58] + b[59] + b[62] - b[63])
              +  144.0*(b[26] - b[27] - b[30] + b[31] + b[41] - b[45] - b[57] + b[61])
              -   27.0*(b[25] - b[29]) ); 
      a[40] = bfac*(
              - 65856.0* b[42]
              + 37632.0*(b[43] + b[46] + b[58])
              + 35280.0* b[38]
              + 28224.0*(b[26] + b[41])
              - 21504.0*(b[47] + b[59] + b[62])
              - 20160.0*(b[39] + b[54])
              - 16128.0*(b[27] + b[30] + b[45] + b[57])
              - 15120.0*(b[22] + b[37])
              + 12288.0* b[63]
              - 12096.0* b[25]
              + 11520.0* b[55]
              +  9216.0*(b[31] + b[61])
              +  8640.0*(b[23] + b[53])
              +  6912.0* b[29]
              +  6480.0* b[21]
              -  5292.0* b[34]
              +  3024.0*(b[35] + b[50])
              +  2268.0*(b[18] + b[33])
              -  1728.0* b[51]
              -  1296.0*(b[19] + b[49])
              -   972.0* b[17] ); 
      a[41] = bfac*(
                47040.0* b[42]
              - 37632.0* b[43]
              - 26880.0*(b[46] + b[58])
              - 25200.0* b[38]
              + 21504.0*(b[47] + b[59])
              - 20160.0*(b[26] - b[39])
              + 16128.0* b[27]
              + 15360.0* b[62]
              + 14400.0* b[54]
              - 12288.0* b[63]
              + 11520.0*(b[30] - b[55])
              + 10800.0* b[22]
              -  9408.0* b[41]
              -  9216.0* b[31]
              -  8640.0* b[23]
              +  5376.0*(b[45] + b[57])
              +  5040.0* b[37]
              +  4032.0* b[25]
              +  3780.0* b[34]
              -  3072.0* b[61]
              -  3024.0* b[35]
              -  2880.0* b[53]
              -  2304.0* b[29]
              -  2160.0*(b[21] + b[50])
              +  1728.0* b[51]
              -  1620.0* b[18]
              +  1296.0* b[19]
              -   756.0* b[33]
              +   432.0* b[49]
              +   324.0* b[17] ); 
      a[42] = bfac*(
              - 47040.0* b[42]
              + 37632.0* b[46]
              + 26880.0*(b[43] + b[58])
              - 21504.0*(b[47] + b[62])
              + 20160.0*(b[26] + b[41])
              + 16464.0* b[38]
              - 16128.0*(b[30] + b[45])
              - 15360.0* b[59]
              + 12288.0* b[63]
              - 11520.0*(b[27] + b[57])
              -  9408.0*(b[39] + b[54])
              +  9216.0*(b[31] + b[61])
              -  8640.0* b[25]
              -  7056.0*(b[22] + b[37])
              +  6912.0* b[29]
              +  5376.0* b[55]
              +  4032.0*(b[23] + b[53])
              +  3024.0* b[21]
              -  1764.0* b[34]
              +  1008.0*(b[35] + b[50])
              +   756.0*(b[18] + b[33])
              -   576.0* b[51]
              -   432.0*(b[19] + b[49])
              -   324.0* b[17] ); 
      a[43] = bfac*(
                33600.0* b[42]
              - 26880.0*(b[43] + b[46])
              + 21504.0* b[47]
              - 19200.0* b[58]
              + 15360.0*(b[59] + b[62])
              - 14400.0* b[26]
              - 12288.0* b[63]
              - 11760.0* b[38]
              + 11520.0*(b[27] + b[30])
              +  9408.0* b[39]
              -  9216.0* b[31]
              -  6720.0*(b[41] - b[54])
              +  5376.0*(b[45] - b[55])
              +  5040.0* b[22]
              -  4032.0* b[23]
              +  3840.0* b[57]
              -  3072.0* b[61]
              +  2880.0* b[25]
              +  2352.0* b[37]
              -  2304.0* b[29]
              -  1344.0* b[53]
              +  1260.0* b[34]
              -  1008.0*(b[21] + b[35])
              -   720.0* b[50]
              +   576.0* b[51]
              -   540.0* b[18]
              +   432.0* b[19]
              -   252.0* b[33]
              +   144.0* b[49]
              +   108.0* b[17] ); 
      a[44] = bfac*(
                47040.0* b[42]
              - 37632.0* b[58]
              - 26880.0*(b[43] + b[46])
              - 25200.0* b[38]
              + 21504.0*(b[59] + b[62])
              - 20160.0*(b[41] - b[54])
              + 16128.0* b[57]
              + 15360.0* b[47]
              + 14400.0* b[39]
              - 12288.0* b[63]
              + 11520.0*(b[45] - b[55])
              + 10800.0* b[37]
              -  9408.0* b[26]
              -  9216.0* b[61]
              -  8640.0* b[53]
              +  5376.0*(b[27] + b[30])
              +  5040.0* b[22]
              +  4032.0* b[25]
              +  3780.0* b[34]
              -  3072.0* b[31]
              -  3024.0* b[50]
              -  2880.0* b[23]
              -  2304.0* b[29]
              -  2160.0*(b[21] + b[35])
              +  1728.0* b[51]
              -  1620.0* b[33]
              +  1296.0* b[49]
              -   756.0* b[18]
              +   432.0* b[19]
              +   324.0* b[17] ); 
      a[45] = bfac*(
              - 33600.0* b[42]
              + 26880.0*(b[43] + b[58])
              - 21504.0* b[59]
              + 19200.0* b[46]
              + 18000.0* b[38]
              - 15360.0*(b[47] + b[62])
              - 14400.0*(b[39] + b[54])
              + 12288.0* b[63]
              + 11520.0* b[55]
              +  6720.0*(b[26] + b[41])
              -  5376.0*(b[27] + b[57])
              -  3840.0*(b[30] + b[45])
              -  3600.0*(b[22] + b[37])
              +  3072.0*(b[31] + b[61])
              +  2880.0*(b[23] + b[53])
              -  2700.0* b[34]
              +  2160.0*(b[35] + b[50])
              -  1728.0* b[51]
              -  1344.0* b[25]
              +   768.0* b[29]
              +   720.0* b[21]
              +   540.0*(b[18] + b[33])
              -   432.0*(b[19] + b[49])
              -   108.0* b[17] ); 
      a[46] = bfac*(
                33600.0* b[42]
              - 26880.0*(b[46] + b[58])
              + 21504.0* b[62]
              - 19200.0* b[43]
              + 15360.0*(b[47] + b[59])
              - 14400.0* b[41]
              - 12288.0* b[63]
              - 11760.0* b[38]
              + 11520.0*(b[45] + b[57])
              +  9408.0* b[54]
              -  9216.0* b[61]
              -  6720.0*(b[26] - b[39])
              +  5376.0*(b[30] - b[55])
              +  5040.0* b[37]
              -  4032.0* b[53]
              +  3840.0* b[27]
              -  3072.0* b[31]
              +  2880.0* b[25]
              +  2352.0* b[22]
              -  2304.0* b[29]
              -  1344.0* b[23]
              +  1260.0* b[34]
              -  1008.0*(b[21] + b[50])
              -   720.0* b[35]
              +   576.0* b[51]
              -   540.0* b[33]
              +   432.0* b[49]
              -   252.0* b[18]
              +   144.0* b[19]
              +   108.0* b[17] ); 
      a[47] = bfac*(
              - 24000.0* b[42]
              + 19200.0*(b[43] + b[46] + b[58])
              - 15360.0*(b[47] + b[59] + b[62])
              + 12288.0* b[63]
              +  8400.0* b[38]
              -  6720.0*(b[39] + b[54])
              +  5376.0* b[55]
              +  4800.0*(b[26] + b[41])
              -  3840.0*(b[27] + b[30] + b[45] + b[57])
              +  3072.0*(b[31] + b[61])
              -  1680.0*(b[22] + b[37])
              +  1344.0*(b[23] + b[53])
              -   960.0* b[25]
              -   900.0* b[34]
              +   768.0* b[29]
              +   720.0*(b[35] + b[50])
              -   576.0* b[51]
              +   336.0* b[21]
              +   180.0*(b[18] + b[33])
              -   144.0*(b[19] + b[49])
              -    36.0* b[17] ); 
      ////////////////////

      a[48] = bfac*(
              -  6144.0* b[42]
              +  3648.0*(b[38] + b[41])
              +  3072.0*(b[26] + b[43] + b[46] + b[58])
              -  2166.0* b[37]
              -  1824.0*(b[22] + b[25] + b[39] + b[45] + b[54] + b[57])
              -  1536.0*(b[27] + b[30] + b[47] + b[59] + b[62])
              +  1083.0*(b[21] + b[53])
              +   912.0*(b[23] + b[29] + b[55] + b[61])
              +   768.0*(b[31] + b[63])
              -   576.0*(b[34] + b[40])
              +   342.0*(b[33] + b[36])
              +   288.0*(b[18] + b[24] + b[35] + b[44] + b[50] + b[56])
              -   171.0*(b[17] + b[20] + b[49] + b[52])
              -   144.0*(b[19] + b[28] + b[51] + b[60])
              -    54.0* b[32]
              +    27.0*(b[16] + b[48]) ); 
      a[49] = bfac*(
              +  3072.0*(b[42] - b[43])
              -  1824.0*(b[38] - b[39])
              -  1536.0*(b[26] - b[27] + b[46] - b[47] + b[58] - b[59])
              +   912.0*(b[22] - b[23] + b[54] - b[55])
              +   768.0*(b[30] - b[31] + b[62] - b[63])
              -   576.0* b[41]
              +   342.0* b[37]
              +  288.0*(b[25]  + b[34] - b[35] + b[45] + b[57])
              -   171.0*(b[21] + b[53])
              -   144.0*(b[18] - b[19] + b[29] + b[50] - b[51] + b[61])
              -    54.0* b[33]
              +    27.0*(b[17] + b[49]) ); 
      a[50] = bfac*(
              +  3072.0*(b[42] - b[46])
              -  1824.0*(b[41] - b[45])
              -  1536.0*(b[26] - b[30] + b[43] - b[47] + b[58] - b[62])
              +   912.0*(b[25] - b[29] + b[57] - b[61])
              +   768.0*(b[27] - b[31] + b[59] - b[63])
              -   576.0* b[38]
              +   342.0* b[37]
              +   288.0*(b[22]  + b[39] + b[40] - b[44] + b[54])
              -   171.0*(b[21]  + b[53])
              -   144.0*(b[23]  + b[24] - b[28] + b[55] + b[56] - b[60])
              -    54.0* b[36]
              +    27.0*(b[20] + b[52]) ); 
      a[51] = bfac*(
              -  1536.0*(b[42] - b[43] - b[46] + b[47])
              +   768.0*(b[26] - b[27] - b[30] + b[31] + b[58] - b[59] - b[62] + b[63])
              +   288.0*(b[38] - b[39] + b[41] - b[45])
              -   144.0*(b[22] - b[23] + b[25] - b[29] + b[54] - b[55] + b[57] - b[61])
              -    54.0* b[37]
              +    27.0*(b[21] + b[53]) ); 
      a[52] = bfac*(
              -  3072.0*(b[42] - b[58])
              +  1824.0*(b[38] + b[41] - b[54] - b[57])
              +  1536.0*(b[43] + b[46] - b[59] - b[62])
              -  1083.0*(b[37] - b[53])
              -   912.0*(b[39] + b[45] - b[55] - b[61])
              -   768.0*(b[47] - b[63])
              -   288.0*(b[34] + b[40] - b[50] - b[56])
              +   171.0*(b[33] + b[36] - b[49] - b[52])
              +   144.0*(b[35] + b[44] - b[51] - b[60])
              -    27.0*(b[32] - b[48]) ); 
      a[53] = bfac*(
                 1536.0*(b[42] - b[43] - b[58] + b[59])
              -   912.0*(b[38] - b[39] - b[54] + b[55])
              -   768.0*(b[46] - b[47] - b[62] + b[63])
              -   288.0*(b[41] - b[57])
              +   171.0*(b[37] - b[53])
              +   144.0*(b[34] - b[35] + b[45] - b[50] + b[51] - b[61])
              -    27.0*(b[33] - b[49]) ); 
      a[54] = bfac*(
                 1536.0*(b[42] - b[46] - b[58] + b[62])
              -   912.0*(b[41] - b[45] - b[57] + b[61])
              -   768.0*(b[43] - b[47] - b[59] + b[63])
              -   288.0*(b[38] - b[54])
              +   171.0*(b[37] - b[53])
              +   144.0*(b[39] + b[40] - b[44] - b[55] - b[56] + b[60])
              -    27.0*(b[36] - b[52]) ); 
      a[55] = bfac*(
              -   768.0*(b[42] - b[43] - b[46] + b[47] - b[58] + b[59] + b[62] - b[63])
              +   144.0*(b[38] - b[39] + b[41] - b[45] - b[54] + b[55] - b[57] + b[61])
              -    27.0*(b[37] - b[53]) ); 
      a[56] = bfac*(
              - 65856.0* b[42]
              + 37632.0*(b[43] + b[46] + b[58])
              + 35280.0* b[26]
              + 28224.0*(b[38] + b[41])
              - 21504.0*(b[47] + b[59] + b[62])
              - 20160.0*(b[27] + b[30])
              - 16128.0*(b[39] + b[45] + b[54] + b[57])
              - 15120.0*(b[22] + b[25])
              + 12288.0* b[63]
              - 12096.0* b[37]
              + 11520.0* b[31]
              +  9216.0*(b[55] + b[61])
              +  8640.0*(b[23] + b[29])
              +  6912.0* b[53]
              +  6480.0* b[21]
              -  5292.0* b[10]
              +  3024.0*(b[11] + b[14])
              +  2268.0*(b[6]  + b[9])
              -  1728.0* b[15]
              -  1296.0*(b[7]  + b[13])
              -   972.0* b[5] ); 
      a[57] =  bfac*(
                47040.0* b[42]
              - 37632.0* b[43]
              - 26880.0*(b[46] + b[58])
              - 25200.0* b[26]
              + 21504.0*(b[47] + b[59])
              + 20160.0*(b[27] - b[38])
              + 16128.0* b[39]
              + 15360.0* b[62]
              + 14400.0* b[30]
              - 12288.0* b[63]
              - 11520.0*(b[31] - b[54])
              + 10800.0* b[22]
              -  9408.0* b[41]
              -  9216.0* b[55]
              -  8640.0* b[23]
              +  5376.0*(b[45] + b[57])
              +  5040.0* b[25]
              +  4032.0* b[37]
              +  3780.0* b[10]
              -  3072.0* b[61]
              -  3024.0* b[11]
              -  2880.0* b[29]
              -  2304.0* b[53]
              -  2160.0*(b[14] + b[21])
              +  1728.0* b[15]
              -  1620.0* b[6]
              +  1296.0* b[7]
              -   756.0* b[9]
              +   432.0* b[13]
              +   324.0* b[5] ); 
      a[58] =  bfac*(
                47040.0* b[42]
              - 37632.0* b[46]
              - 26880.0*(b[43] + b[58])
              - 25200.0* b[26]
              + 21504.0*(b[47] + b[62])
              + 20160.0*(b[30] - b[41])
              + 16128.0* b[45]
              + 15360.0* b[59]
              + 14400.0* b[27]
              - 12288.0* b[63]
              - 11520.0*(b[31] - b[57])
              + 10800.0* b[25]
              -  9408.0* b[38]
              -  9216.0* b[61]
              -  8640.0* b[29]
              +  5376.0*(b[39] + b[54])
              +  5040.0* b[22]
              +  4032.0* b[37]
              +  3780.0* b[10]
              -  3072.0* b[55]
              -  3024.0* b[14]
              -  2880.0* b[23]
              -  2304.0* b[53]
              -  2160.0*(b[11] + b[21])
              +  1728.0* b[15]
              -  1620.0* b[9]
              +  1296.0* b[13]
              -   756.0* b[6]
              +   432.0* b[7]
              +   324.0* b[5] ); 
      a[59] = bfac*(
              - 33600.0* b[42]
              + 26880.0*(b[43] + b[46])
              - 21504.0* b[47]
              + 19200.0* b[58]
              + 18000.0* b[26]
              - 15360.0*(b[59] + b[62])
              - 14400.0*(b[27] + b[30])
              + 12288.0* b[63]
              + 11520.0* b[31]
              +  6720.0*(b[38] + b[41])
              -  5376.0*(b[39] + b[45])
              -  3840.0*(b[54] + b[57])
              -  3600.0*(b[22] + b[25])
              +  3072.0*(b[55] + b[61])
              +  2880.0*(b[23] + b[29])
              -  2700.0* b[10]
              +  2160.0*(b[11] + b[14])
              -  1728.0* b[15]
              -  1344.0* b[37]
              +   768.0* b[53]
              +   720.0* b[21]
              +   540.0*(b[6]  + b[9])
              -   432.0*(b[7]  + b[13])
              -   108.0* b[5] ); 
      a[60] = bfac*(
              - 47040.0* b[42]
              + 37632.0* b[58]
              + 26880.0*(b[43] + b[46])
              - 21504.0*(b[59] + b[62])
              + 20160.0*(b[38] + b[41])
              + 16464.0* b[26]
              - 16128.0*(b[54] + b[57])
              - 15360.0* b[47]
              + 12288.0* b[63]
              - 11520.0*(b[39] + b[45])
              -  9408.0*(b[27] + b[30])
              +  9216.0*(b[55] + b[61])
              -  8640.0* b[37]
              -  7056.0*(b[22] + b[25])
              +  6912.0* b[53]
              +  5376.0* b[31]
              +  4032.0*(b[23] + b[29])
              +  3024.0* b[21]
              -  1764.0* b[10]
              +  1008.0*(b[11] + b[14])
              +   756.0*(b[6]  + b[9])
              -   576.0* b[15]
              -   432.0*(b[7] + b[13])
              -   324.0* b[5] ); 
      a[61] = bfac*(
                33600.0* b[42]
              - 26880.0*(b[43] + b[58])
              + 21504.0* b[59]
              - 19200.0* b[46]
              + 15360.0*(b[47] + b[62])
              - 14400.0* b[38]
              - 12288.0* b[63]
              - 11760.0* b[26]
              + 11520.0*(b[39] + b[54])
              +  9408.0* b[27]
              -  9216.0* b[55]
              +  6720.0*(b[30] - b[41])
              -  5376.0*(b[31] - b[57])
              +  5040.0* b[22]
              -  4032.0* b[23]
              +  3840.0* b[45]
              -  3072.0* b[61]
              +  2880.0* b[37]
              +  2352.0* b[25]
              -  2304.0* b[53]
              -  1344.0* b[29]
              +  1260.0* b[10]
              -  1008.0*(b[11] + b[21])
              -   720.0* b[14]
              +   576.0* b[15]
              -   540.0* b[6]
              +   432.0* b[7]
              -   252.0* b[9]
              +   144.0* b[13]
              +   108.0* b[5] );
      a[62] = bfac*(
                33600.0* b[42]
              - 26880.0*(b[46] + b[58])
              + 21504.0* b[62]
              - 19200.0* b[43]
              + 15360.0*(b[47] + b[59])
              - 14400.0* b[41]
              - 12288.0* b[63]
              - 11760.0* b[26]
              + 11520.0*(b[45] + b[57])
              +  9408.0* b[30]
              -  9216.0* b[61]
              +  6720.0*(b[27] - b[38])
              -  5376.0*(b[31] - b[54])
              +  5040.0* b[25]
              -  4032.0* b[29]
              +  3840.0* b[39]
              -  3072.0* b[55]
              +  2880.0* b[37]
              +  2352.0* b[22]
              -  2304.0* b[53]
              -  1344.0* b[23]
              +  1260.0* b[10]
              -  1008.0*(b[14] + b[21])
              -   720.0* b[11]
              +   576.0* b[15]
              -   540.0* b[9]
              +   432.0* b[13]
              -   252.0* b[6]
              +   144.0* b[7]
              +   108.0* b[5] ); 
      a[63] = bfac*(
              - 24000.0* b[42]
              + 19200.0*(b[43] + b[46] + b[58])
              - 15360.0*(b[47] + b[59] + b[62])
              + 12288.0* b[63]
              +  8400.0* b[26]
              -  6720.0*(b[27] + b[30])
              +  5376.0* b[31]
              +  4800.0*(b[38] + b[41])
              -  3840.0*(b[39] + b[45] + b[54] + b[57])
              +  3072.0*(b[55] + b[61])
              -  1680.0*(b[22] + b[25])
              +  1344.0*(b[23] + b[29])
              -   960.0* b[37]
              -   900.0* b[10]
              +   768.0* b[53]
              +   720.0*(b[11] + b[14])
              -   576.0* b[15]
              +   336.0* b[21]
              +   180.0*(b[6]  + b[9])
              -   144.0*(b[7]  + b[13])
              -    36.0* b[5] ); 
}

//*************************************************************************//

// performs explicit multiplication by Binv matrix
// B2 

void BinvMultiply_B2(const double b[], double a[])
{

      a[0] =  27.0  * b[0];

      a[1] =  27.0  * b[16];
      
      a[2] =  81.0  * (b[1]-b[0])       
            - 54.0  * b[16] - 27.0  * b[17];

      a[3] =  54.0  * (b[0]-b[1])
            + 27.0  * (b[16]+b[17]);

      a[4] = 27.0  * b[32];

      a[5] = - 1083.0  * b[0]
             + 171.0   * (b[1]+b[2]-b[16]-b[32])
             + 27.0    * (b[18]+b[33]-b[3]) 
             + 1728.0  * b[12]
             - 576.0   * (b[13]+b[14])
             + 192.0   * b[15]
             - 972.0   * b[56]
             + 324.0   * (b[57]+b[58]-b[60])
             + 108.0   * (b[61]+b[62]-b[59])
             - 36.0    * b[63];

      a[6] =   2337.0  * b[0]  - 1425.0  * b[1]
             - 369.0   * b[2]  + 225.0   * b[3]
             - 4032.0  * b[12] + 2880.0  * b[13]
             + 1344.0  * b[14] - 960.0   * b[15]
             + 342.0   * b[16] + 171.0   * b[17]
             - 54.0    * b[18] - 27.0    * b[19]
             + 288.0   * b[32] - 144.0   * b[33]
             + 2268.0  * b[56] - 1620.0  * b[57]
             + 756.0   * (b[60]-b[58]) 
             + 540.0   * (b[59]-b[61])
             - 252.0   * b[62] + 180.0  * b[63];
      
      a[7] =   1254.0  * (b[1]-b[0])
             + 198.0   * (b[2]-b[3])
             + 2304.0  * (b[12]-b[13])
             + 768.0   * (b[15]-b[14])
             - 171.0   * (b[16]+b[17])
             + 27.0    * (b[18]+b[19])
             + 144.0   * (b[33]-b[32])
             + 1296.0  * (b[57]-b[56])
             + 432.0   * (b[58]-b[59]-b[60]+b[61])
             + 144.0   * (b[62]-b[63]);

      a[8] =   81.0  * (b[2] - b[0])
             - 54.0  * b[32] - 27.0  * b[34];

      a[9] =  2337.0  * b[0]  - 369.0  * b[1]
             - 1425.0  * b[2]  + 225.0  * b[3]
             - 4032.0  * b[12] + 1344.0  * b[13]
             + 2880.0  * b[14] - 960.0   * b[15]
             + 288.0   * b[16] - 144.0   * b[18]
             + 342.0   * b[32] - 54.0    * b[33]
             + 171.0   * b[34] - 27.0    * b[35]
             + 2268.0  * b[56] 
             + 756.0   * (b[60]-b[57])
             - 1620.0  * b[58] 
             + 540.0   * (b[59]-b[62])
             - 252.0   * b[61]  + 180.0 * b[63];

      a[10] =   4800.0  * (b[15]-b[0])  
                          + 2832.0  * (b[1]+b[2])
              - 1632.0  * b[3] + 9408.0  * b[12]
              - 6720.0  * (b[13]+b[14]) 
              - 576.0   * (b[16]+b[32]) 
              + 288.0   * (b[18]+b[33]-b[17]-b[34])
              + 144.0   * (b[19]+b[35]) 
              - 5292.0  * b[56]
              + 3780.0  * (b[57]+b[58]) 
              - 2700.0  * b[59]
              - 1764.0  * b[60]  
              + 1260.0  * (b[61]+b[62])
              - 900.0   * b[63];
     
      a[11] =   2544.0  * (b[0]-b[1]) 
              + 1488.0  * (b[3]-b[2])
              + 5376.0  * (b[13]-b[12]) 
              + 3840.0  * (b[14]-b[15])
              + 288.0   * (b[16]+b[17]+b[32]-b[33]) 
              + 144.0   * (b[34]-b[18]-b[19]-b[35])
              + 3024.0  * (b[56]-b[57]) 
              + 2160.0  * (b[59]-b[58])
              + 1008.0  * (b[60]-b[61]) 
              + 720.0  * (b[63]-b[62]);

      a[12] = 54.0  * (b[0]-b[2]) 
              + 27.0  * (b[32]+b[34]);

      a[13] =   1254.0  * (b[2]-b[0]) 
              + 198.0  * (b[1]-b[3])
              + 2304.0  * (b[12]-b[14]) 
              + 768.0  * (b[15]-b[13])
              + 144.0   * (b[18]-b[16] 
              + b[61]-b[63])
              - 171.0   * (b[32]+b[34]) 
              + 27.0  * (b[33]+b[35])
              + 1296.0  * (b[58]-b[56]) 
              + 432.0   * (b[57]-b[59]-b[60]+b[62]);

      a[14] =   2544.0  * (b[0]-b[2]) 
              + 1488.0  * (b[3]-b[1])
              + 5376.0  * (b[14]-b[12]) 
              + 3840.0  * (b[13]-b[15])
              + 288.0   * (b[16]-b[18]+b[32]+b[34])
              + 144.0   * (b[17]-b[19]-b[33]-b[35])
              + 3024.0  * (b[56]-b[58]) 
              + 2160.0  * (b[59]-b[57])
              + 1008.0  * (b[60]-b[62]) 
              + 720.0   * (b[63]-b[61]);

      a[15] =   1344.0  * (b[1]-b[0]+b[2]-b[3])
              + 3072.0  * (b[12]-b[13]-b[14]+b[15])
              + 144.0   * (b[18]-b[16]-b[17]+b[19]
                          -  b[32]+b[33]-b[34]+b[35])
              + 1728.0  * (b[57]-b[56]+b[58]-b[59])
              + 576.0   * (b[61]-b[60]+b[62]-b[63]);

      a[16] =   27.0  * b[48];

      a[17] =  -1083.0  * b[0] 
              + 171.0  * (b[1]+b[4]-b[16]-b[48])
              + 27.0  * (b[20]-b[5]+b[49]) 
              + 1728.0  * b[10]
              - 576.0  * (b[11]+b[14]) 
              + 192.0  * b[15] 
              - 972.0  * b[40] 
              + 324.0  * (b[41]-b[42]+b[44])
              + 108.0  * (b[43]-b[45]+b[46]) 
              - 36.0  * b[47];
     
      a[18] =   2337.0  * b[0] - 1425.0  * b[1] 
              - 369.0  * b[4]
              + 225.0   * b[5] - 4032.0  * b[10] 
              + 2880.  * b[11]
              + 1344.0  * b[14] - 960.0  * b[15] 
              + 342.0  * b[16]
              + 171.0   * b[17] - 54.0   * b[20] 
              - 27.0  * b[21]
              + 2268.0  * b[40] - 1620.0  * b[41] 
              + 756.0  * (b[42]-b[44])
              + 540.0   * (b[45]-b[43]) 
              - 252.0  * b[46] 
              + 180.0  * b[47]
              + 288.0   * b[48] - 144.0  * b[49];

      a[19] =   1254.0  * (b[1]-b[0]) 
              + 198.0  * (b[4]-b[5])
              + 2304.0  * (b[10]-b[11]) 
              + 768.0  * (b[15]-b[14])
              - 171.0   * (b[16]+b[17]) 
              + 27.0  * (b[20]+b[21])
              + 1296.0  * (b[41]-b[40]) 
              + 432.0  * (b[43]-b[42]+b[44]-b[45]) 
              + 144.0  * (b[46]-b[47]-b[48]+b[49]);

      a[20] = - 1083.0  * b[0] 
              + 171.0  * (b[2]+b[4]-b[32]-b[48])
              + 27.0    * (b[36]-b[6]+b[50]) 
              + 1728.0  * b[9]
              - 576.0   * (b[11]+b[13]) 
              + 192.0  * b[15] - 972.0  * b[24]
              + 324.0   * (b[26]-b[25]+b[28]) 
              + 108.0  * (b[27]+b[29]-b[30])
              - 36.0    * b[31];

      a[21] =   13718.0  * b[0] 
              - 2166.0  * (b[1]+b[2]+b[4])
              + 342.0    * (b[3]+b[5]+b[6]) 
              - 54.0  * b[7]
              + 19008.0  * b[8] 
              - 16704.0  * (b[9]+b[10]+b[12])
              + 9024.0   * (b[11]+b[13]+b[14]) 
              - 4160.0 *b[15]
              + 1083.0   * (b[16]+b[32]+b[48]) 
              - 171.0    * (b[18]+b[20]+b[33]
                           +b[36]+b[49]+b[50])
              + 27.0     * (b[22]+b[37]+b[51]) 
              + 6480.0  * (b[24]+b[40]+b[56])
              + 3024.0   * (b[25]+b[42]+b[60]) 
              - 2160.0   * (b[26]+b[28]+b[41]
                           +b[44]+b[57]+b[58])
              - 1008.0   * (b[27]+b[29]+b[43]
                           +b[46]+b[61]+b[62])
              + 720.0    * (b[30]+b[45]+b[59]) 
              + 336.0  * (b[31]+b[47]+b[63]);

      a[22] = - 26353.0  * b[0] + 14801.0  * b[1] 
              + 4161.0  * (b[2]+b[4])
              - 2337.0   * (b[3]+b[5]) 
              - 657.0  * b[6] + 369.0  * b[7]
              - 49536.0  * b[8] + 47232.0  * b[9] 
              + 40704.0  * (b[10]+b[12])
              - 33024.0  * (b[11]+b[13]) 
              - 21632.0  * b[14] + 16768.0  * b[15]
              - 2166.0  * b[16] - 1083.0  * b[17] 
              + 342.0  * (b[18]+b[20])
              + 171.0   * (b[19]+b[21]) - 54.0  * b[22] 
              - 27.0  * b[23]
              - 12096.0  * b[24] - 8640.0  * b[25] 
              + 4032.0  * (b[26]+b[28])
              + 2880.0  * (b[27]+b[29]) 
              - 1344.0  * b[30] - 960.0  * b[31]
              - 1824.0  * (b[32]+b[48]) 
              + 912.0  * (b[33]+b[49]) 
              + 288.0  * (b[36]+b[50]) 
              - 144.0   * (b[37]+b[51]) 
              - 15120.0  * (b[40]+b[56]) 
              + 10800.0  * (b[41]+b[57]) 
              - 7056.0  * (b[42]+b[60])  
              + 5040.0  * (b[43]+b[44]+b[58]+b[61]) 
              - 3600.0  * (b[45]+b[59])
              + 2352.0  * (b[46]+b[62]) 
              - 1680.0  * (b[47]+b[63]); 

      a[23] =   13718.0  * (b[0]-b[1]) 
              + 2166.0  * (b[3]-b[2]-b[4]+b[5])
              + 342.0  * (b[6]-b[7]) 
              + 32256.0  * (b[8]-b[9])
              + 24576.0  * (b[11]-b[10]-b[12]+b[13]) 
              + 12800.0  * (b[14]-b[15])
              + 1083.0  * (b[16]+b[17]) 
              - 171.0  * (b[18]+b[19]+b[20]+b[21])
              + 27.0  * (b[22]+b[23]) 
              + 6912.0  * (b[24]+b[25])
              - 2304.0  * (b[26]+b[27]+b[28]+b[29]) 
              + 768.0  * (b[30]+b[31])
              + 912.0  * (b[32]-b[33]+b[48]-b[49]) 
              + 144.0  * (b[37]-b[36]-b[50]+b[51])
              + 8640.0  * (b[40]-b[41]+b[56]-b[57]) 
              + 4032.0  * (b[42]-b[43]+b[60]-b[61]) 
              + 2880.0  * (b[45]-b[44]-b[58]+b[59]) 
              + 1344.0  * (b[47]-b[46]-b[62]+b[63]);

      a[24] = - 4032.0  * b[9] + 2880.0  * b[11] 
              + 2337.0  * b[0] 
              + 2268.0  * b[24] - 1620.0  * b[26] 
              - 1425.0  * b[2]
              + 1344.0  * b[13] - 960.0   * b[15] 
              + 756.0   * (b[25]-b[28])
              + 540.0   * (b[30]-b[27]) 
              - 369.0  * b[4] + 342.0  * b[32]
              + 288.0   * b[48] 
              - 252.0  * b[29] + 225.0  * b[6] 
              + 180.0   * b[31] 
              + 171.0  * b[34] - 144.0  * b[50]
              - 54.0  * b[36] - 27.0  * b[38];

      a[25] = - 49536.0  * b[8] + 47232.0  * b[10] 
              + 40704.0  * (b[9]+b[12])
              - 33024.0  * (b[11]+b[14]) 
              - 26353.0  * b[0] - 21632.0  * b[13]
              + 16768.0  * b[15] 
              - 15120.0  * (b[24]+b[56]) 
              + 14801.0  * b[2] - 12096.0  * b[40]
              + 10800.0  * (b[26]+b[58]) 
              - 8640.0  * b[42] 
              - 7056.0   * (b[25]+b[60]) 
              + 5040.0  * (b[27]+b[28]+b[57]+b[62])
              + 4161.0   * (b[1]+b[4]) 
              + 4032.0  * (b[41]+b[44])
              - 3600.0   * (b[30]+b[59]) 
              + 2880.0  * (b[43]+b[46])
              + 2352.0   * (b[29]+b[61]) 
              - 2337.0  * (b[3]+b[6])
              - 2166.0   * b[32] 
              - 1824.0  * (b[16]+b[48])
              - 1680.0   * (b[31]+b[63]) 
              - 1344.0  * b[45]
              - 1083.0   * b[34] - 960.0  * b[47] 
              + 912.0  * (b[18]+b[50])
              - 657.0  * b[5] + 369.0  * b[7] 
              + 342.0  * (b[33]+b[36])
              + 288.0  * (b[20]+b[49]) 
              + 171.0  * (b[35]+b[38])
              - 144.0  * (b[22]+b[51]) 
              - 54.0  * b[37] - 27.0  * b[39];

      a[26] = 127680.0  * b[8] 
              - 118848.0  * (b[9]+b[10]) 
              + 104640.0  * b[11]
              - 99008.0  * b[12] 
              + 79936.0  * (b[13]+b[14]) 
              - 63680.0  * b[15]
              + 49856.0  * b[0] + 35280.0  * b[56] 
              + 28224.0  * (b[24]+b[40])
              - 27664.0  * (b[1]+b[2]) 
              - 25200.0  * (b[57]+b[58]) 
              + 20160.0  * (b[25]-b[26]-b[41]+b[42]) 
              + 18000.0  * b[59]
              + 16464.0  * b[60] + 15200.0  * b[3] 
              - 14400.0  * (b[27]+b[43])
              - 11760.0  * (b[61]+b[62]) 
              - 9408.0  * (b[28]+b[44]) 
              + 8400.0  * b[63] - 7872.0  * b[4] 
              + 6720.0  * (b[30]-b[29]+b[45]-b[46]) 
              + 4800.0  * (b[31]+b[47])
              + 4368.0  * (b[5]+b[6]) 
              + 3648.0  * (b[16]+b[32]) 
              + 3072.0  * b[48] - 2400.0  * b[7] 
              + 1824.0  * (b[17]-b[18]-b[33]+b[34]) 
              - 1536.0  * (b[49]+b[50])
              - 912.0   * (b[19]+b[35]) 
              + 768.0  * b[51] 
              - 576.0  * (b[20]+b[36]) 
              + 288.0  * (b[22]-b[21]+b[37]-b[38])
              + 144.0  * (b[23]+b[39]);

      a[27] = 82176.0  * (b[9]-b[8]) 
              + 74496.0  * (b[10]-b[11])
              + 59648.0  * (b[12]-b[13]) 
              + 47872.0  * (b[15]-b[14])
              + 25840.0  * (b[1]-b[0]) 
              + 20160.0  * (b[57]-b[56])
              + 16128.0  * (b[41]-b[24]-b[25]-b[40]) 
              + 14400.0  * (b[58]-b[59])
              + 14288.0  * (b[2]-b[3]) 
              + 11520.0  * (b[26]+b[27]-b[42]+b[43])
              + 9408.0   * (b[61]-b[60]) 
              + 6720.0  * (b[62]-b[63])
              + 5376.0   * (b[28]+b[29]+b[44]-b[45]) 
              + 4080.0  * (b[4]-b[5])
              + 3840.0   * (b[46]-b[30]-b[31]-b[47]) 
              + 2256.0  * (b[7]-b[6])
              + 1824.0   * (b[33]-b[16]-b[17]-b[32]) 
              + 1536.0  * (b[49]-b[48])
              + 912.0    * (b[18]+b[19]-b[34]+b[35]) 
              + 768.0  * (b[50]-b[51])
              + 288.0    * (b[20]+b[21]+b[36]-b[37])
              + 144.0    * (b[38]-b[22]-b[23]-b[39]);

      a[28] = 2304.0  * (b[9]-b[11]) 
              + 1296.0  * (b[26]-b[24])
              + 1254.0  * (b[2]-b[0]) 
              + 768.0  * (b[15]-b[13])
              + 432.0   * (b[27]-b[25]+b[28]-b[30]) 
              + 198.0  * (b[4]-b[6])
              - 171.0   * (b[32]+b[34]) 
              + 144.0  * (b[29]-b[31]-b[48]+b[50])
              + 27.0    * (b[36]+b[38]);

      a[29] = 32256.0  * (b[8]-b[10]) 
              + 24576.0  * (b[11]-b[9]-b[12]+b[14])
              + 13718.0  * (b[0]-b[2]) 
              + 12800.0  * (b[13]-b[15])
              + 8640.0  * (b[24]-b[26]+b[56]-b[58]) 
              + 6912.0  * (b[40]+b[42])
              + 4032.0  * (b[25]-b[27]+b[60]-b[62])
              + 2880.0  * (b[30]-b[28]-b[57]+b[59])
              - 2304.0  * (b[41]+b[43]+b[44]+b[46])
              + 2166.0  * (b[3]-b[1]-b[4]+b[6])
              + 1344.0  * (b[31]-b[29]-b[61]+b[63]) 
              + 1083.0  * (b[32]+b[34])
              + 912.0   * (b[16]-b[18]+b[48]-b[50]) 
              + 768.0  * (b[45]+b[47])
              + 342.0   * (b[5]-b[7]) 
              - 171.0  * (b[33]+b[35]+b[36]+b[38])
              + 144.0   * (b[22]-b[20]-b[49]+b[51]) 
              + 27.0  * (b[37]+b[39]);

      a[30] = 82176.0  * (b[10]-b[8]) 
              + 74496.0  * (b[9]-b[11])
              + 59648.0  * (b[12]-b[14]) 
              + 47872.0  * (b[15]-b[13])
              + 25840.0  * (b[2]-b[0]) 
              + 20160.0  * (b[58]-b[56])
              + 16128.0  * (b[26]-b[24]-b[40]-b[42]) 
              + 14400.0  * (b[57]-b[59])
              + 14288.0  * (b[1]-b[3]) 
              + 11520.0  * (b[27]-b[25]+b[41]+b[43])
              + 9408.0   * (b[62]-b[60]) 
              + 6720.0  * (b[61]-b[63])
              + 5376.0   * (b[28]-b[30]+b[44]+b[46]) 
              + 4080.0  * (b[4]-b[6])
              + 3840.0   * (b[29]-b[31]-b[45]-b[47]) 
              + 2256.0  * (b[7]-b[5])
              + 1824.0   * (b[18]-b[16]-b[32]-b[34]) 
              + 1536.0  * (b[50]-b[48])
              + 912.0    * (b[19]-b[17]+b[33]+b[35]) 
              + 768.0  * (b[49]-b[51])
              + 288.0  * (b[20]-b[22]+b[36]+b[38]) 
              + 144.0  * (b[21]-b[23]-b[37]-b[39]);

      a[31] = 52224.0  * (b[8]-b[9]-b[10]+b[11]) 
              + 35840.0  * (b[13]-b[12]+b[14]-b[15])
              + 13376.0  * (b[0]-b[1]-b[2]+b[3])
              + 11520.0  * (b[56]-b[57]-b[58]+b[59])
              + 9216.0   * (b[24]+b[25]-b[26]-b[27]
                           +  b[40]-b[41]+b[42]-b[43])
              + 5376.0   * (b[60]-b[61]-b[62]+b[63])
              + 3072.0   * (b[30]-b[28]-b[29]+b[31]
                           - b[44]+b[45]-b[46]+b[47])
              + 2112.0   * (b[5]-b[4]+b[6]-b[7])
              + 912.0    * (b[16]+b[17]-b[18]-b[19]
                           + b[32]-b[33]+b[34]-b[35])
              + 768.0    * (b[48]-b[49]-b[50]+b[51])
              + 144.0    * (b[22]-b[20]-b[21]+b[23]
                           -  b[36]+b[37]-b[38]+b[39]);

      a[32] = 81.0  * (b[4]-b[0]) 
             - 54.0  * b[48] - 27.0  * b[52];

      a[33] = - 4032.0  * b[10] 
              + 2880.0  * b[14] + 2337.0  * b[0]
              + 2268.0  * b[40] 
              - 1620.0  * b[44] - 1425.0  * b[4]
              + 1344.0  * b[11] 
              - 960.0  * b[15] + 756.0  * (b[42]-b[41])
              + 540.0   * (b[45]-b[46]) 
              - 369.0  * b[1] + 342.0  * b[48]
              + 288.0  * b[16] - 252.0  * b[43] 
              + 225.0  * b[5]
              + 180.0  * b[47] + 171.0  * b[52] 
              - 144.0  * b[20]
              - 54.0   * b[49] - 27.0  * b[53];

      a[34] = 9408.0  * b[10] - 6720.0  * (b[11]+b[14]) 
              - 5292.0  * b[40]
              + 4800.0  * (b[15]-b[0]) 
              + 3780.0  * (b[41]+b[44])
              + 2832.0  * (b[1]+b[4]) - 2700.0  * b[45] 
              - 1764.0  * b[42]
              - 1632.0  * b[5] + 1260.0  * (b[43]+b[46]) 
              - 900.0  * b[47]
              - 576.0   * (b[16]+b[48]) 
              + 288.0  * (b[20]-b[17]+b[49]-b[52])
              + 144.0   * (b[21]+b[53]);

      a[35] = 5376.0  * (b[11]-b[10]) 
              + 3840.0  * (b[14]-b[15])
              + 3024.0  * (b[40]-b[41]) 
              + 2544.0  * (b[0]-b[1])
              + 2160.0  * (b[45]-b[44]) 
              + 1488.0  * (b[5]-b[4])
              + 1008.0  * (b[42]-b[43]) 
              + 720.0  * (b[47]-b[46])
              + 288.0   * (b[16]+b[17]+b[48]-b[49])
              + 144.0   * (b[52]-b[20]-b[21]-b[53]);

      a[36] = - 4032.0  * b[9] + 2880.0  * b[13] 
              + 2337.0  * b[0]
              + 2268.0  * b[24] - 1620.0  * b[28] 
              - 1425.0  * b[4]
              + 1344.0  * b[11] - 960.0  * b[15] 
              + 756.0  * (b[25]-b[26])
              + 540.0   * (b[30]-b[29]) - 369.0  * b[2] 
              + 342.0  * b[48]
              + 288.0   * b[32] - 252.0  * b[27] 
              + 225.0  * b[6]
              + 180.0   * b[31] + 171.0  * b[52] 
              - 144.0  * b[36]
              - 54.0    * b[50] - 27.0  * b[54];

      a[37] = - 49536.0  * b[8] + 47232.0  * b[12] 
              + 40704.0  * (b[9]+b[10])
              - 33024.0  * (b[13]+b[14]) 
              - 26353.0  * b[0] - 21632.0  * b[11]
              + 16768.0  * b[15] - 15120.0  * (b[24]+b[40]) 
              + 14801.0  * b[4]
              - 12096.0  * b[56] + 10800.0  * (b[28]+b[44]) 
              - 8640.0  * b[60]
              - 7056.0   * (b[25]+b[42]) 
              + 5040.0  * (b[26]+b[29]+b[41]+b[46])
              + 4161.0   * (b[1]+b[2]) 
              + 4032.0  * (b[57]+b[58])
              - 3600.0   * (b[30]+b[45]) 
              + 2880.0  * (b[61]+b[62])
              + 2352.0   * (b[27]+b[43]) 
              - 2337.0  * (b[5]+b[6]) 
              - 2166.0   * b[48] - 1824.0  * (b[16]+b[32]) 
              - 1680.0   * (b[31]+b[47]) 
              - 1344.0  * b[59] - 1083.0  * b[52]
              - 960.0    * b[63] + 912.0  * (b[20]+b[36]) 
              - 657.0  * b[3]
              + 369.0    * b[7] + 342.0  * (b[49]+b[50]) 
              + 288.0  * (b[18]+b[33]) 
              + 171.0  * (b[53]+b[54]) 
              - 144.0  * (b[22]+b[37]) 
              - 54.0  * b[51] - 27.0  * b[55];

      a[38] = 127680.0  * b[8] - 118848.0  * (b[9]+b[12]) 
              + 104640.0  * b[13]
              - 99008.0  * b[10] + 79936.0  * (b[11]+b[14]) 
              - 63680.0  * b[15]
              + 49856.0  * b[0] + 35280.0  * b[40] 
              + 28224.0  * (b[24]+b[56])
              - 27664.0  * (b[1]+b[4]) 
              - 25200.0  * (b[41]+b[44])
              + 20160.0  * (b[25]-b[28]-b[57]+b[60]) 
              + 18000.0  * b[45]
              + 16464.0  * b[42] + 15200.0  * b[5] 
              - 14400.0  * (b[29]+b[61])
              - 11760.0  * (b[43]+b[46]) 
              - 9408.0  * (b[26]+b[58])
              + 8400.0   * b[47] - 7872.0  * b[2] 
              + 6720.0  * (b[30]-b[27]+b[59]-b[62]) 
              + 4800.0  * (b[31]+b[63])
              + 4368.0  * (b[3]+b[6]) 
              + 3648.0  * (b[16]+b[48]) 
              + 3072.0  * b[32] - 2400.0  * b[7] 
              + 1824.0  * (b[17]-b[20]-b[49]+b[52]) 
              - 1536.0  * (b[33]+b[36])
              - 912.0   * (b[21]+b[53]) + 768.0  * b[37] 
              - 576.0  * (b[18]+b[50]) 
               + 288.0  * (b[22]-b[19]+b[51]-b[54])
              + 144.0  * (b[23]+b[55]);

      a[39] = 82176.0  * (b[9]-b[8]) 
              + 74496.0  * (b[12]-b[13])
              + 59648.0  * (b[10]-b[11]) 
              + 47872.0  * (b[15]-b[14])
              + 25840.0  * (b[1]-b[0]) 
              + 20160.0  * (b[41]-b[40])
              + 16128.0  * (b[57]-b[24]-b[25]-b[56]) 
              + 14400.0  * (b[44]-b[45])
              + 14288.0  * (b[4]-b[5]) 
              + 11520.0  * (b[28]+b[29]-b[60]+b[61]) 
              + 9408.0   * (b[43]-b[42]) 
              + 6720.0  * (b[46]-b[47])
              + 5376.0   * (b[26]+b[27]+b[58]-b[59]) 
              + 4080.0  * (b[2]-b[3])
              + 3840.0   * (b[62]-b[30]-b[31]-b[63]) 
              + 2256.0  * (b[7]-b[6])
              + 1824.0   * (b[49]-b[16]-b[17]-b[48]) 
              + 1536.0  * (b[33]-b[32])
              + 912.0    * (b[20]+b[21]-b[52]+b[53]) 
              + 768.0  * (b[36]-b[37])
              + 288.0    * (b[18]+b[19]+b[50]-b[51]) 
              + 144.0    * (b[54]-b[22]-b[23]-b[55]);

      a[40] = 9408.0  * b[9] - 6720.0  * (b[11]+b[13]) 
              - 5292.0  * b[24]
              + 4800.0  * (b[15]-b[0]) 
              + 3780.0  * (b[26]+b[28])
              + 2832.0  * (b[2]+b[4]) - 2700.0  * b[30] 
              - 1764.0  * b[25]
              - 1632.0  * b[6] + 1260.0  * (b[27]+b[29]) 
              - 900.0  * b[31]
              - 576.0   * (b[32]+b[48]) 
              + 288.0  * (b[36]-b[34]+b[50]-b[52])
              + 144.0   * (b[38]+b[54]);

      a[41] = 127680.0  * b[8] - 118848.0  * (b[10]+b[12]) 
              + 104640.0  * b[14] - 99008.0  * b[9]
              + 79936.0  * (b[11]+b[13]) - 63680.0  * b[15] 
              + 49856.0  * b[0]
              + 35280.0  * b[24] + 28224.0  * (b[40]+b[56]) 
              - 27664.0  * (b[2]+b[4]) 
              - 25200.0  * (b[26]+b[28]) 
              + 20160.0  * (b[42]-b[44]-b[58]+b[60]) 
              + 18000.0  * b[30]
              + 16464.0  * b[25] + 15200.0  * b[6] 
              - 14400.0  * (b[46]+b[62])
              - 11760.0  * (b[27]+b[29]) 
              - 9408.0  * (b[41]+b[57]) 
              + 8400.0   * b[31] - 7872.0  * b[1] 
              + 6720.0   * (b[45]-b[43]+b[59]-b[61]) 
              + 4800.0  * (b[47]+b[63])
              + 4368.0   * (b[3]+b[5]) 
              + 3648.0  * (b[32]+b[48]) 
              + 3072.0  * b[16] - 2400.0  * b[7] 
              + 1824.0  * (b[34]-b[36]-b[50]+b[52]) 
              - 1536.0  * (b[18]+b[20])
              - 912.0   * (b[38]+b[54]) + 768.0  * b[22] 
              - 576.0   * (b[33]+b[49]) 
              + 288.0  * (b[37]-b[35]+b[51]-b[53])
              + 144.0   * (b[39]+b[55]);

      a[42] = - 326144.0  * b[8] 
              + 297472.0  * (b[9]+b[10]+b[12])
              - 258560.0  * (b[11]+b[13]+b[14]) 
              + 217600.0  * b[15]
              - 93184.0   * b[0] 
              - 65856.0  * (b[24]+b[40]+b[56])
              + 51200.0   * (b[1]+b[2]+b[4])
              + 47040.0   * (b[26]-b[25]+b[28]+b[41]-b[42]
                            +  b[44]+b[57]+b[58]-b[60])
              + 33600.0   * (b[27]+b[29]-b[30]+b[43]-b[45]
                            +  b[46]-b[59]+b[61]+b[62])
              - 27904.0   * (b[3]+b[5]+b[6]) 
              - 24000.0  * (b[31]+b[47]+b[63])
              + 15104.0   * b[7] 
              - 6144.0  * (b[16]+b[32]+b[48])
              + 3072.0    * (b[18]-b[17]+b[20]+b[33]-b[34]
                            +  b[36]+b[49]+b[50]-b[52])
              + 1536.0    * (b[19]+b[21]-b[22]+b[35]-b[37]
                            + b[38]-b[51]+b[53]+b[54])
              - 768.0     * (b[23]+b[39]+b[55]);

      a[43] = 207872.0  * (b[8]-b[9]) 
              + 185344.0  * (b[11]-b[10]-b[12]+b[13])
              + 158720.0  * (b[14]-b[15]) 
              + 48128.0  * (b[0]-b[1])
              + 37632.0   * (b[24]+b[25]
                            + b[40]-b[41]+b[56]-b[57])
              + 26880.0   * (b[42]-b[26]-b[27]-b[28]-b[29]
                            -  b[43]-b[44]+b[45]-b[58]
                            +  b[59]+b[60]-b[61]) 
              + 26368.0   * (b[3]-b[2]-b[4]+b[5]) 
              + 19200.0   * (b[30]+b[31]-b[46]+b[47]-b[62]+b[63])
              + 14336.0   * (b[6]-b[7]) 
              + 3072.0    * (b[16]+b[17]+b[32]-b[33]+b[48]-b[49])  
              + 1536.0    * (b[34]-b[18]-b[19]-b[20]-b[21]
                            -  b[35]-b[36]+b[37]-b[50]
                            +  b[51]+b[52]-b[53])
              + 768.0     * (b[22]+b[23]-b[38]+b[39]-b[54]+b[55]);

      a[44] = 5376.0  * (b[11]-b[9]) 
              + 3840.0  * (b[13]-b[15])
              + 3024.0  * (b[24]-b[26]) 
              + 2544.0  * (b[0]-b[2]) 
              + 2160.0  * (b[30]-b[28]) 
              + 1488.0  * (b[6]-b[4])
              + 1008.0  * (b[25]-b[27]) 
              + 720.0  * (b[31]-b[29])
              + 288.0   * (b[32]+b[34]+b[48]-b[50]) 
              + 144.0  * (b[52]-b[36]-b[38]-b[54]); 

      a[45] = 82176.0  * (b[10]-b[8]) 
              + 74496.0  * (b[12]-b[14])
              + 59648.0  * (b[9]-b[11]) 
              + 47872.0  * (b[15]-b[13])
              + 25840.0  * (b[2]-b[0]) 
              + 20160.0  * (b[26]-b[24])
              + 16128.0  * (b[58]-b[40]-b[42]-b[56]) 
              + 14400.0  * (b[28]-b[30])
              + 14288.0  * (b[4]-b[6]) 
              + 11520.0  * (b[44]+b[46]-b[60]+b[62])
              + 9408.0   * (b[27]-b[25]) 
              + 6720.0  * (b[29]-b[31])
              + 5376.0   * (b[41]+b[43]+b[57]-b[59]) 
              + 4080.0  * (b[1]-b[3])
              + 3840.0   * (b[61]-b[45]-b[47]-b[63]) 
              + 2256.0  * (b[7]-b[5])
              + 1824.0   * (b[50]-b[32]-b[34]-b[48]) 
              + 1536.0  * (b[18]-b[16])
              + 912.0    * (b[36]+b[38]-b[52]+b[54]) 
              + 768.0  * (b[20]-b[22])
              + 288.0    * (b[33]+b[35]+b[49]-b[51])
              + 144.0    * (b[53]-b[37]-b[39]-b[55]);

      a[46] = 207872.0  * (b[8]-b[10]) 
              + 185344.0  * (b[11]-b[9]-b[12]+b[14])
              + 158720.0  * (b[13]-b[15]) 
              + 48128.0  * (b[0]-b[2])
              + 37632.0   * (b[24]-b[26]+b[40]+b[42]+b[56]-b[58])
              + 26880.0   * (b[25]-b[27]-b[28]+b[30]
                            -  b[41]-b[43]-b[44]-b[46]-b[57]
                            +  b[59]+b[60]-b[62]) 
              + 26368.0   * (b[3]-b[1]-b[4]+b[6]) 
              + 19200.0   * (b[31]-b[29]+b[45]+b[47]-b[61]+b[63])
              + 14336.0   * (b[5]-b[7]) 
              + 3072.0    * (b[16]-b[18]+b[32]+b[34]+b[48]-b[50])
              + 1536.0    * (b[17]-b[19]-b[20]+b[22]-b[33]-b[35]
                            -  b[36]-b[38]-b[49]
                            +  b[51]+b[52]-b[54])
              + 768.0     * (b[23]-b[21]+b[37]+b[39]-b[53]+b[55]);
  
      a[47] = 131072.0  * (b[9]-b[8]+b[10]-b[11]) 
              + 114688.0  * (b[12]-b[13]-b[14]+b[15])
              + 24832.0   * (b[1]-b[0]+b[2]-b[3]) 
              + 21504.0   * (b[26]-b[24]-b[25]+b[27]-b[40]
                            + b[41]-b[42]+b[43]-b[56]
                            +  b[57]+b[58]-b[59]) 
              + 15360.0   * (b[28]+b[29]-b[30]-b[31]+b[44]
                            -  b[45]+b[46]-b[47]-b[60]
                            +  b[61]+b[62]-b[63])
              + 13568.0   * (b[4]-b[5]-b[6]+b[7])
              + 1536.0    * (b[18]-b[16]-b[17]+b[19]-b[32]
                            + b[33]-b[34]+b[35]-b[48]
                            +  b[49]+b[50]-b[51])
              + 768.0     * (b[20]+b[21]-b[22]-b[23]+b[36]
                            -  b[37]+b[38]-b[39]-b[52]
                            +  b[53]+b[54]-b[55]);

      a[48] = 54.0  * (b[0]-b[4]) + 27.0  * (b[48]+b[52]);

      a[49] = 2304.0  * (b[10]-b[14]) + 1296.0  * (b[44]-b[40]) 
              + 1254.0  * (b[4]-b[0]) + 768.0  * (b[15]-b[11])
              + 432.0   * (b[41]-b[42]-b[45]+b[46]) 
              + 198.0  * (b[1]-b[5])
              - 171.0   * (b[48]+b[52]) 
              + 144.0  * (b[20]-b[16]+b[43]-b[47])
              + 27.0    * (b[49]+b[53]);

      a[50] = 5376.0  * (b[14]-b[10]) + 3840.0  * (b[11]-b[15])
              + 3024.0  * (b[40]-b[44]) + 2544.0  * (b[0]-b[4])
              + 2160.0  * (b[45]-b[41]) + 1488.0  * (b[5]-b[1])
              + 1008.0  * (b[42]-b[46]) + 720.0  * (b[47]-b[43])
              + 288.0   * (b[16]-b[20]+b[48]+b[52])
              + 144.0   * (b[17]-b[21]-b[49]-b[53]);

      a[51] = 3072.0  * (b[10]-b[11]-b[14]+b[15])
              + 1728.0  * (b[41]-b[40]+b[44]-b[45])
              + 1344.0  * (b[1]-b[0]+b[4]-b[5]) 
              + 576.0   * (b[43]-b[42]+b[46]-b[47])
              + 144.0   * (b[20]-b[16]-b[17]+b[21]
                          -  b[48]+b[49]-b[52]+b[53]);

      a[52] = 2304.0  * (b[9]-b[13]) + 1296.0  * (b[28]-b[24])
              + 1254.0  * (b[4]-b[0]) + 768.0  * (b[15]-b[11])
              + 432.0   * (b[26]-b[25]+b[29]-b[30]) 
              + 198.0  * (b[2]-b[6])
              - 171.0   * (b[48]+b[52]) 
              + 144.0   * (b[27]-b[31]-b[32]+b[36]) 
              + 27.0  * (b[50]+b[54]);

      a[53] = 32256.0  * (b[8]-b[12]) 
              + 24576.0  * (b[13]-b[9]-b[10]+b[14])
              + 13718.0  * (b[0]-b[4]) + 12800.0  * (b[11]-b[15])
              + 8640.0   * (b[24]-b[28]+b[40]-b[44]) 
              + 6912.0  * (b[56]+b[60])
              + 4032.0   * (b[25]-b[29]+b[42]-b[46])
              + 2880.0   * (b[30]-b[26]-b[41]+b[45])
              - 2304.0   * (b[57]+b[58]+b[61]+b[62])
              + 2166.0   * (b[5]-b[1]-b[2]+b[6])
              + 1344.0   * (b[31]-b[27]-b[43]+b[47]) 
              + 1083.0  * (b[48]+b[52])
              + 912.0    * (b[16]-b[20]+b[32]-b[36]) 
              + 768.0  * (b[59]+b[63])
              + 342.0    * (b[3]-b[7]) 
              - 171.0  * (b[49]+b[50]+b[53]+b[54])
              + 144.0    * (b[22]-b[18]-b[33]+b[37]) 
              + 27.0  * (b[51]+b[55]);

      a[54] = 82176.0  * (b[12]-b[8]) + 74496.0  * (b[9]-b[13])
              + 59648.0  * (b[10]-b[14]) 
              + 47872.0  * (b[15]-b[11])
              + 25840.0  * (b[4]-b[0]) + 20160.0  * (b[44]-b[40])
              + 16128.0  * (b[28]-b[24]-b[56]-b[60]) 
              + 14400.0  * (b[41]-b[45])
              + 14288.0  * (b[1]-b[5]) 
              + 11520.0  * (b[29]-b[25]+b[57]+b[61])
              + 9408.0   * (b[46]-b[42]) + 6720.0  * (b[43]-b[47])
              + 5376.0   * (b[26]-b[30]+b[58]+b[62]) 
              + 4080.0  * (b[2]-b[6])
              + 3840.0   * (b[27]-b[31]-b[59]-b[63]) 
              + 2256.0  * (b[7]-b[3])
              + 1824.0   * (b[20]-b[16]-b[48]-b[52]) 
              + 1536.0  * (b[36]-b[32])
              + 912.0    * (b[21]-b[17]+b[49]+b[53]) 
              + 768.0  * (b[33]-b[37])
              + 288.0    * (b[18]-b[22]+b[50]+b[54])
              + 144.0    * (b[19]-b[23]-b[51]-b[55]);

      a[55] = 52224.0  * (b[8]-b[9]-b[12]+b[13])
              + 35840.0  * (b[11]-b[10]+b[14]-b[15])
              + 13376.0  * (b[0]-b[1]-b[4]+b[5])
              + 11520.0  * (b[40]-b[41]-b[44]+b[45])
              + 9216.0   * (b[24]+b[25]-b[28]-b[29]
                           +  b[56]-b[57]+b[60]-b[61])
              + 5376.0   * (b[42]-b[43]-b[46]+b[47])
              + 3072.0   * (b[30]-b[26]-b[27]+b[31]
                           -  b[58]+b[59]-b[62]+b[63])
              + 2112.0   * (b[3]-b[2]+b[6]-b[7])
              + 912.0    * (b[16]+b[17]-b[20]-b[21]
                           +  b[48]-b[49]+b[52]-b[53])
              + 768.0    * (b[32]-b[33]-b[36]+b[37])
              + 144.0    * (b[22]-b[18]-b[19]+b[23]
                           -  b[50]+b[51]-b[54]+b[55]);

      a[56] = 5376.0  * (b[13]-b[9]) + 3840.0  * (b[11]-b[15])
              + 3024.0  * (b[24]-b[28]) + 2544.0  * (b[0]-b[4])
              + 2160.0  * (b[30]-b[26]) + 1488.0  * (b[6]-b[2])
              + 1008.0  * (b[25]-b[29]) + 720.0  * (b[31]-b[27])
              + 288.0   * (b[32]-b[36]+b[48]+b[52])
              + 144.0   * (b[34]-b[38]-b[50]-b[54]);

      a[57] = 82176.0  * (b[12]-b[8]) + 74496.0  * (b[10]-b[14])
              + 59648.0  * (b[9]-b[13]) + 47872.0  * (b[15]-b[11])
              + 25840.0  * (b[4]-b[0]) + 20160.0  * (b[28]-b[24])
              + 16128.0  * (b[44]-b[40]-b[56]-b[60]) 
              + 14400.0  * (b[26]-b[30])
              + 14288.0  * (b[2]-b[6]) 
              + 11520.0  * (b[46]-b[42]+b[58]+b[62])
              + 9408.0   * (b[29]-b[25]) + 6720.0  * (b[27]-b[31])
              + 5376.0   * (b[41]-b[45]+b[57]+b[61]) 
              + 4080.0  * (b[1]-b[5])
              + 3840.0   * (b[43]-b[47]-b[59]-b[63]) 
              + 2256.0  * (b[7]-b[3])
              + 1824.0   * (b[36]-b[32]-b[48]-b[52]) 
              + 1536.0  * (b[20]-b[16])
              + 912.0    * (b[38]-b[34]+b[50]+b[54]) 
              + 768.0  * (b[18]-b[22])
              + 288.0    * (b[33]-b[37]+b[49]+b[53])
              + 144.0    * (b[35]-b[39]-b[51]-b[55]);

      a[58] = 207872.0  * (b[8]-b[12]) 
              + 185344.0  * (b[13]-b[9]-b[10]+b[14])
              + 158720.0  * (b[11]-b[15]) 
              + 48128.0  * (b[0]-b[4])
              + 37632.0   * (b[24]-b[28]+b[40]
                            -b[44]+b[56]+b[60])
              + 26880.0   * (b[25]-b[26]-b[29]+b[30]-b[41]
                            +b[42]+b[45]-b[46]-b[57]
                            -  b[58]-b[61]-b[62]) 
              + 26368.0   * (b[5]-b[1]-b[2]+b[6]) 
              + 19200.0   * (b[31]-b[27]-b[43]+b[47]
                            +b[59]+b[63])
              + 14336.0   * (b[3]-b[7]) 
              + 3072.0    * (b[16]-b[20]+b[32]-b[36]
                            +b[48]+b[52])  
              + 1536.0    * (b[17]-b[18]-b[21]+b[22]-b[33]
                            +b[34]+b[37]-b[38]-b[49]
                            -  b[50]-b[53]-b[54])
              + 768.0     * (b[23]-b[19]-b[35]+b[39]+b[51]+b[55]);
 
      a[59] = 131072.0  * (b[9]-b[8]+b[12]-b[13]) 
              + 114688.0  * (b[10]-b[11]-b[14]+b[15])
              + 24832.0   * (b[1]-b[0]+b[4]-b[5]) 
              + 21504.0   * (b[28]-b[24]-b[25]+b[29]-b[40]
                            +b[41]+b[44]-b[45]-b[56]
                            +  b[57]-b[60]+b[61]) 
              + 15360.0   * (b[26]+b[27]-b[30]-b[31]-b[42]
                            +b[43]+b[46]-b[47]+b[58]
                            -  b[59]+b[62]-b[63])
              + 13568.0   * (b[2]-b[3]-b[6]+b[7])
              + 1536.0    * (b[20]-b[16]-b[17]+b[21]-b[32]
                            +b[33]+b[36]-b[37]-b[48]
                            +  b[49]-b[52]+b[53])
              + 768.0     * (b[18]+b[19]-b[22]-b[23]-b[34]
                            +b[35]+b[38]-b[39]+b[50]
                            -  b[51]+b[54]-b[55]);

      a[60] = 3072.0  * (b[9]-b[11]-b[13]+b[15])
              + 1728.0  * (b[26]-b[24]+b[28]-b[30])
              + 1344.0  * (b[2]-b[0]+b[4]-b[6]) 
              + 576.0  * (b[27]-b[25]+b[29]-b[31])      
              + 144.0  * (b[36]-b[32]-b[34]+b[38]-b[48]
                         +b[50]-b[52]+b[54]);

      a[61] = 52224.0  * (b[8]-b[10]-b[12]+b[14])
              + 35840.0  * (b[11]-b[9]+b[13]-b[15])
              + 13376.0  * (b[0]-b[2]-b[4]+b[6])
              + 11520.0  * (b[24]-b[26]-b[28]+b[30])
              + 9216.0   * (b[40]+b[42]-b[44]-b[46]+b[56]
              -b[58]+b[60]-b[62])
              + 5376.0   * (b[25]-b[27]-b[29]+b[31])
              + 3072.0   * (b[45]-b[41]-b[43]+b[47]
              -b[57]+b[59]-b[61]+b[63])
              + 2112.0   * (b[3]-b[1]+b[5]-b[7])
              + 912.0    * (b[32]+b[34]-b[36]-b[38]+b[48]
              -b[50]+b[52]-b[54])
              + 768.0    * (b[16]-b[18]-b[20]+b[22])
              + 144.0    * (b[37]-b[33]-b[35]+b[39]-b[49]
              +b[51]-b[53]+b[55]);

      a[62] = 131072.0  * (b[10]-b[8]+b[12]-b[14])
              + 114688.0  * (b[9]-b[11]-b[13]+b[15])
              + 24832.0   * (b[2]-b[0]+b[4]-b[6])
              + 21504.0   * (b[26]-b[24]+b[28]-b[30]-b[40]
                            -b[42]+b[44]+b[46]-b[56]
                            +  b[58]-b[60]+b[62])
              + 15360.0   * (b[27]-b[25]+b[29]-b[31]+b[41]
                            +b[43]-b[45]-b[47]+b[57]
                            -  b[59]+b[61]-b[63])
              + 13568.0   * (b[1]-b[3]-b[5]+b[7])
              + 1536.0    * (b[18]-b[16]+b[20]-b[22]-b[32]
                            -b[34]+b[36]+b[38]-b[48]
                            +  b[50]-b[52]+b[54])
              + 768.0     * (b[19]-b[17]+b[21]-b[23]+b[33]
                            +b[35]-b[37]-b[39]+b[49]
                            -  b[51]+b[53]-b[55]);

      a[63] =   81920.0  * (b[8]-b[9]-b[10]+b[11]-b[12]
                           +b[13]+b[14]-b[15])
              + 12800.0  * (b[0]-b[1]-b[2]+b[3]-b[4]+b[5]
                           +b[6]-b[7])
              + 12288.0  * (b[24]+b[25]-b[26]-b[27]-b[28]
                           -b[29]+b[30]+b[31]
                           +  b[40]-b[41]+b[42]-b[43]-b[44]
                           +b[45]-b[46]+b[47]
                           +  b[56]-b[57]-b[58]+b[59]+b[60]
                           -b[61]-b[62]+b[63])
              + 768.0    * (b[16]+b[17]-b[18]-b[19]-b[20]
                           -b[21]+b[22]+b[23]
                           +  b[32]-b[33]+b[34]-b[35]-b[36]
                           +b[37]-b[38]+b[39]
                           +  b[48]-b[49]-b[50]+b[51]+b[52]
                           -b[53]-b[54]+b[55]);
}

//*************************************************************************//

// performs explicit multiplication by Transpose(Binv) matrix
// Lekien&Marsden

void BinvMultiply_Transpose_LM(const double b[], double a[])
{
      a[0]  = - 27.0*b[42] 
              + 18.0*(b[43] + b[46] + b[58])
              - 12.0*(b[47] + b[59] + b[62]) 
              + 9.0 *(b[10] + b[34] + b[40]) 
              + 8.0 *b[63]
              - 6.0 *(b[11] + b[14] + b[35] + b[44] + b[50] + b[56])
              + 4.0 *(b[15] + b[51] + b[60])
              - 3.0 *(b[2]  + b[8]  + b[32])
              + 2.0 *(b[3]  + b[12] + b[48])
              + b[0];

      a[1]  =   27.0*b[42]
              - 18.0*(b[43] + b[46] + b[58])
              + 12.0*(b[47] + b[59] + b[62])
              - 9.0 *(b[10] + b[34])
              - 8.0 *b[63]
              + 6.0 *(b[11] + b[14] + b[35] + b[50])
              - 4.0 *(b[15] + b[51])
              + 3.0 *b[2]
              - 2.0 *b[3];

      a[2]  = 27.0*b[42]
              - 18.0*(b[43] + b[46] + b[58])
              + 12.0*(b[47] + b[59] + b[62])
              - 9.0 *(b[10] + b[40])
              - 8.0 *b[63]
              + 6.0 *(b[11] + b[14] + b[44] + b[56])
              - 4.0 *(b[15] + b[60])
              + 3.0 *b[8]
              - 2.0 *b[12]; 

      a[3]  = - 27.0*b[42]
              + 18.0*(b[43] + b[46] + b[58])
              - 12.0*(b[47] + b[59] + b[62]) 
              + 9.0 *b[10] 
              + 8.0 *b[63]
              - 6.0 *(b[11] + b[14])
              + 4.0 *b[15];

      a[4]  = 27.0*b[42]
              - 18.0*(b[43] + b[46] + b[58])
              + 12.0*(b[47] + b[59] + b[62])
              - 9.0 *(b[34] + b[40])
              - 8.0 *b[63]
              + 6.0 *(b[35] + b[44] + b[50] + b[56])
              - 4.0 *(b[51] + b[60])
              + 3.0 *b[32]
              - 2.0 *b[48];

      a[5]  = - 27.0*b[42]
              + 18.0*(b[43] + b[46] + b[58])
              - 12.0*(b[47] + b[59] + b[62])
              + 9.0 *b[34]
              + 8.0 *b[63]
              - 6.0 *(b[35] + b[50])
              + 4.0 *b[51];

      a[6]  = - 27.0*b[42]
              + 18.0*(b[43] + b[46] + b[58])
              - 12.0*(b[47] + b[59] + b[62])
              + 9.0 *b[40]
              + 8.0 *b[63]
              - 6.0 *(b[44] + b[56])
              + 4.0 *b[60];

      a[7]  =   27.0*b[42]
              - 18.0*(b[43] + b[46] + b[58])
              + 12.0*(b[47] + b[59] + b[62])
              - 8.0 *b[63];

      a[8]  = - 18.0*b[42]
              + 12.0*(b[46] + b[58])
              + 9.0 *(b[41] + b[43])
              - 8.0 *b[62]
              + 6.0 *(b[10] + b[34] - b[45] - b[47] - b[57] - b[59])
              - 4.0 *(b[14] + b[50] - b[61] - b[63])
              - 3.0 *(b[9]  + b[11] + b[33] + b[35])
              - 2.0 *(b[2]  - b[13] - b[15] - b[49] - b[51])
              + b[1] + b[3];

      a[9]  = - 9.0 *(b[42] - b[43])
              + 6.0 *(b[46] - b[47] + b[58] - b[59])
              - 4.0 *(b[62] - b[63])
              + 3.0 *(b[10] - b[11] + b[34] - b[35])
              - 2.0 *(b[14] - b[15] + b[50] - b[51])
              - b[2] + b[3];

      a[10] =   18.0*b[42]
              - 12.0*(b[46] + b[58])
              - 9.0 *(b[41] + b[43])
              + 8.0 *b[62]
              - 6.0 *(b[10] - b[45] - b[47] - b[57] - b[59])
              + 4.0 *(b[14] - b[61] - b[63])
              + 3.0 *(b[9]  + b[11])
              - 2.0 *(b[13] + b[15]);

      a[11] =   9.0 *(b[42] - b[43])
              - 6.0 *(b[46] - b[47] + b[58] - b[59])
              + 4.0 *(b[62] - b[63])
              - 3.0 *(b[10] - b[11])
              + 2.0 *(b[14] - b[15]); 

      a[12] =   18.0*b[42]
              - 12.0*(b[46] + b[58])
              - 9.0 *(b[41] + b[43])
              + 8.0 *b[62]
              - 6.0 *(b[34] - b[45] - b[47] - b[57] - b[59])
              + 4.0 *(b[50] - b[61] - b[63])
              + 3.0 *(b[33] + b[35])
              - 2.0 *(b[49] + b[51]);

      a[13] =   9.0 *(b[42] - b[43])
              - 6.0 *(b[46] - b[47] + b[58] - b[59])
              + 4.0 *(b[62] - b[63])
              - 3.0 *(b[34] - b[35])
              + 2.0 *(b[50] - b[51]);

      a[14] = - 18.0*b[42]
              + 12.0*(b[46] + b[58])
              + 9.0 *(b[41] + b[43])
              - 8.0 *b[62]
              - 6.0 *(b[45] + b[47] + b[57] + b[59])
              + 4.0 *(b[61] + b[63]);

      a[15] = - 9.0 *(b[42] - b[43])
              + 6.0 *(b[46] - b[47] + b[58] - b[59])
              - 4.0 *(b[62] - b[63]); 
      ////////////////////

      a[16] = - 18.0*b[42]
              + 12.0*(b[43] + b[58])
              + 9.0 *(b[38] + b[46])
              - 8.0 *b[59]
              + 6.0 *(b[10] - b[39] + b[40] - b[47] - b[54] - b[62])
              - 4.0 *(b[11] - b[55] + b[56] - b[63])
              - 3.0 *(b[6]  + b[14] + b[36] + b[44])
              + 2.0 *(b[7]  - b[8]  + b[15] + b[52] + b[60])
              + b[4] + b[12];
 
      a[17] =   18.0*b[42]
              - 12.0*(b[43] + b[58])
              - 9.0 *(b[38] + b[46])
              + 8.0 *b[59]
              - 6.0 *(b[10] - b[39] - b[47] - b[54] - b[62])
              + 4.0 *(b[11] - b[55] - b[63])
              + 3.0 *(b[6]  + b[14])
              - 2.0 *(b[7]  + b[15]);

      a[18] = - 9.0 *(b[42] - b[46])
              + 6.0 *(b[43] - b[47] + b[58] - b[62])
              - 4.0 *(b[59] - b[63])
              + 3.0 *(b[10] - b[14] + b[40] - b[44])
              - 2.0 *(b[11] - b[15] + b[56] - b[60])
              -     b[8]    + b[12];

      a[19] =   9.0 *(b[42] - b[46])
              - 6.0 *(b[43] - b[47] + b[58] - b[62])
              + 4.0 *(b[59] - b[63])
              - 3.0 *(b[10] - b[14])
              + 2.0 *(b[11] - b[15]);

      a[20] =   18.0*b[42]
              - 12.0*(b[43] + b[58])
              - 9.0 *(b[38] + b[46])
              + 8.0 *b[59]
              + 6.0 *(b[39] - b[40] + b[47] + b[54] + b[62])
              - 4.0 *(b[55] - b[56] + b[63])
              + 3.0 *(b[36] + b[44])
              - 2.0 *(b[52] + b[60]);

      a[21] = - 18.0*b[42]
              + 12.0*(b[43] + b[58])
              + 9.0 *(b[38] + b[46])
              - 8.0 *b[59]
              - 6.0 *(b[39] + b[47] + b[54] + b[62])
              + 4.0 *(b[55] + b[63]);

      a[22] =   9.0 *(b[42] - b[46])
              - 6.0 *(b[43] - b[47] + b[58] - b[62])
              + 4.0 *(b[59] - b[63])
              - 3.0 *(b[40] - b[44])
              + 2.0 *(b[56] - b[60]);

      a[23] = - 9.0 *(b[42] - b[46])
              + 6.0 *(b[43] - b[47] + b[58] - b[62])
              - 4.0 *(b[59] - b[63]);

      a[24] = - 18.0*b[42]
              + 12.0*(b[43] + b[46])
              + 9.0 *(b[26] + b[58])
              - 8.0 *b[47]
              - 6.0 *(b[27] + b[30] - b[34] - b[40] + b[59] + b[62])
              + 4.0 *(b[31] - b[35] - b[44] + b[63])
              - 3.0 *(b[18] + b[24] + b[50] + b[56])
              + 2.0 *(b[19] + b[28] - b[32] + b[51] + b[60])
              + b[16] + b[48];

      a[25] =   18.0*b[42]
              - 12.0*(b[43] + b[46])
              - 9.0 *(b[26] + b[58])
              + 8.0 *b[47]
              + 6.0 *(b[27] + b[30] - b[34] + b[59] + b[62])
              - 4.0 *(b[31] - b[35] + b[63])
              + 3.0 *(b[18] + b[50])
              - 2.0 *(b[19] + b[51]);
 
      a[26] =   18.0*b[42]
              - 12.0*(b[43] + b[46])
              - 9.0 *(b[26] + b[58])
              + 8.0 *b[47]
              + 6.0 *(b[27] + b[30] - b[40] + b[59] + b[62])
              - 4.0 *(b[31] - b[44] + b[63])
              + 3.0 *(b[24] + b[56])
              - 2.0 *(b[28] + b[60]); 

      a[27] = - 18.0*b[42]
              + 12.0*(b[43] + b[46])
              + 9.0 *(b[26] + b[58])
              - 8.0 *b[47]
              - 6.0 *(b[27] + b[30] + b[59] + b[62])
              + 4.0 *(b[31] + b[63]);

      a[28] = - 9.0 *(b[42] - b[58])
              + 6.0 *(b[43] + b[46] - b[59] - b[62])
              - 4.0 *(b[47] - b[63])
              + 3.0 *(b[34] + b[40] - b[50] - b[56])
              - 2.0 *(b[35] + b[44] - b[51] - b[60])
              -  b[32] + b[48];

      a[29] =   9.0 *(b[42] - b[58])
              - 6.0 *(b[43] + b[46] - b[59] - b[62])
              + 4.0 *(b[47] - b[63])
              - 3.0 *(b[34] - b[50])
              + 2.0 *(b[35] - b[51]);
      
      a[30] =   9.0 *(b[42] - b[58])
              - 6.0 *(b[43] + b[46] - b[59] - b[62])
              + 4.0 *(b[47] - b[63])
              - 3.0 *(b[40] - b[56])
              + 2.0 *(b[44] - b[60]); 

      a[31] = - 9.0 *(b[42] - b[58])
              + 6.0 *(b[43] + b[46] - b[59] - b[62])
              - 4.0 *(b[47] - b[63]); 
      ////////////////////

      a[32] = - 12.0*b[42]
              + 8.0 *b[58]
              + 6.0 *(b[38] + b[41] + b[43] + b[46])
              + 4.0 *(b[10] - b[54] - b[57] - b[59] - b[62])
              - 3.0 *(b[37] + b[39] + b[45] + b[47])
              - 2.0 *(b[6]  + b[9]  + b[11] + b[14] - b[53] - b[55] - b[61] - b[63])
              +       b[5]  + b[7]  + b[13] + b[15];
 
      a[33] = - 6.0 *(b[42] - b[43])
              + 4.0 *(b[58] - b[59])
              + 3.0 *(b[38] - b[39] + b[46] - b[47])
              + 2.0 *(b[10] - b[11] - b[54] + b[55] - b[62] + b[63])
              -       b[6]  + b[7]  - b[14] + b[15];

      a[34] = - 6.0 *(b[42] - b[46])
              + 4.0 *(b[58] - b[62])
              + 3.0 *(b[41] + b[43] - b[45] - b[47])
              + 2.0 *(b[10] - b[14] - b[57] - b[59] + b[61] + b[63])
              -       b[9]  - b[11] + b[13] + b[15];

      a[35] = - 3.0 *(b[42] - b[43] - b[46] + b[47])
              + 2.0 *(b[58] - b[59] - b[62] + b[63])
              +       b[10] - b[11] - b[14] + b[15]; 

      a[36] =   12.0*b[42]
              - 8.0 *b[58]
              - 6.0 *(b[38] + b[41] + b[43] + b[46])
              + 4.0 *(b[54] + b[57] + b[59] + b[62])
              + 3.0 *(b[37] + b[39] + b[45] + b[47])
              - 2.0 *(b[53] + b[55] + b[61] + b[63]);

      a[37] =   6.0 *(b[42] - b[43])
              - 4.0 *(b[58] - b[59])
              - 3.0 *(b[38] - b[39] + b[46] - b[47])
              + 2.0 *(b[54] - b[55] + b[62] - b[63]);

      a[38] =   6.0 *(b[42] - b[46])
              - 4.0 *(b[58] - b[62])
              - 3.0 *(b[41] + b[43] - b[45] - b[47])
              + 2.0 *(b[57] + b[59] - b[61] - b[63]);

      a[39] =   3.0 *(b[42] - b[43] - b[46] + b[47])
              - 2.0 *(b[58] - b[59] - b[62] + b[63]); 

      a[40] = - 12.0*b[42]
              + 8.0 *b[46]
              + 6.0 *(b[26] + b[41] + b[43] + b[58])
              - 4.0 *(b[30] - b[34] + b[45] + b[47] + b[62])
              - 3.0 *(b[25] + b[27] + b[57] + b[59])
              - 2.0 *(b[18] - b[29] - b[31] + b[33] + b[35] + b[50] - b[61] - b[63])
              +       b[17] + b[19] + b[49] + b[51]; 

      a[41] =  - 6.0 *(b[42] - b[43])
              + 4.0 *(b[46] - b[47])
              + 3.0 *(b[26] - b[27] + b[58] - b[59])
              - 2.0 *(b[30] - b[31] - b[34] + b[35] + b[62] - b[63])
              -       b[18] + b[19]  - b[50] + b[51];

      a[42] =   12.0*b[42]
              - 8.0 *b[46]
              - 6.0 *(b[26] + b[41] + b[43] + b[58])
              + 4.0 *(b[30] + b[45] + b[47] + b[62])
              + 3.0 *(b[25] + b[27] + b[57] + b[59])
              - 2.0 *(b[29] + b[31] + b[61] + b[63]);

      a[43] =   6.0 *(b[42] - b[43])
              - 4.0 *(b[46] - b[47])
              - 3.0 *(b[26] - b[27] + b[58] - b[59])
              + 2.0 *(b[30] - b[31] + b[62] - b[63]);


      a[44] =  - 6.0 *(b[42] - b[58])
              + 4.0 *(b[46] - b[62])
              + 3.0 *(b[41] + b[43] - b[57] - b[59])
              + 2.0 *(b[34] - b[45] - b[47] - b[50] + b[61] + b[63])
              -       b[33] - b[35] + b[49] + b[51];
 

      a[45] = - 3.0 *(b[42] - b[43] - b[58] + b[59])
              + 2.0 *(b[46] - b[47] - b[62] + b[63])
              +       b[34] - b[35] - b[50] + b[51];


      a[46] =   6.0 *(b[42] - b[58])
              - 4.0 *(b[46] - b[62])
              - 3.0 *(b[41] + b[43] - b[57] - b[59])
              + 2.0 *(b[45] + b[47] - b[61] - b[63]);

      a[47] =   3.0 *(b[42] - b[43] - b[58] + b[59])
              - 2.0 *(b[46] - b[47] - b[62] + b[63]);  

      ////////////////////

      a[48] = - 12.0*b[42]
              + 8.0 *b[43]
              + 6.0 *(b[26] + b[38] + b[46] + b[58])
              - 4.0 *(b[27] + b[39] - b[40] + b[47] + b[59])
              - 3.0 *(b[22] + b[30] + b[54] + b[62])
              + 2.0 *(b[23] - b[24] + b[31] - b[36] - b[44] + b[55] - b[56] + b[63])
              +       b[20] + b[28] + b[52] + b[60];
 
      a[49] =   12.0*b[42]
              - 8.0 *b[43]
              - 6.0 *(b[26] + b[38] + b[46] + b[58])
              + 4.0 *(b[27] + b[39] + b[47] + b[59])
              + 3.0 *(b[22] + b[30] + b[54] + b[62])
              - 2.0 *(b[23] + b[31] + b[55] + b[63]);

      a[50] = - 6.0 *(b[42] - b[46])
              + 4.0 *(b[43] - b[47])
              + 3.0 *(b[26] - b[30] + b[58] - b[62])
              - 2.0 *(b[27] - b[31] - b[40] + b[44] + b[59] - b[63])
              -       b[24] + b[28] - b[56] + b[60]; 

      a[51] =   6.0 *(b[42] - b[46])
              - 4.0 *(b[43] - b[47])
              - 3.0 *(b[26] - b[30] + b[58] - b[62])
              + 2.0 *(b[27] - b[31] + b[59] - b[63]);
	
      a[52] =  - 6.0 *(b[42] - b[58])
              + 4.0 *(b[43] - b[59])
              + 3.0 *(b[38] + b[46] - b[54] - b[62])
              - 2.0 *(b[39] - b[40] + b[47] - b[55] + b[56] - b[63])
              -       b[36] - b[44] + b[52] + b[60];

      a[53] =   6.0 *(b[42] - b[58])
              - 4.0 *(b[43] - b[59])
              - 3.0 *(b[38] + b[46] - b[54] - b[62])
              + 2.0 *(b[39] + b[47] - b[55] - b[63]);
	
      a[54] = - 3.0 *(b[42] - b[46] - b[58] + b[62])
              + 2.0 *(b[43] - b[47] - b[59] + b[63])
              +       b[40] - b[44] - b[56] + b[60]; 

      a[55] =   3.0 *(b[42] - b[46] - b[58] + b[62])
              - 2.0 *(b[43] - b[47] - b[59] + b[63]); 
 
      a[56] = - 8.0 *b[42]
              + 4.0 *(b[26] + b[38] + b[41] + b[43] + b[46] + b[58])
              - 2.0 *(b[22] + b[25] + b[27] + b[30] + b[37] + b[39] + b[45] + b[47] + b[54] + b[57] + b[59] + b[62])
              +       b[21] + b[23] + b[29] + b[31] + b[53] + b[55] + b[61] + b[63]; 

      a[57] = - 4.0 *(b[42] - b[43])
              + 2.0 *(b[26] - b[27] + b[38] - b[39] + b[46] - b[47] + b[58] - b[59])
              -       b[22] + b[23] - b[30] + b[31] - b[54] + b[55] - b[62] + b[63];

      a[58] = - 4.0 *(b[42] - b[46])
              + 2.0 *(b[26] - b[30] + b[41] + b[43] - b[45] - b[47] + b[58] - b[62])
              -       b[25] - b[27] + b[29] + b[31] - b[57] - b[59] + b[61] + b[63];
 
      a[59] = - 2.0 *(b[42] - b[43] - b[46] + b[47])
              +       b[26] - b[27] - b[30] + b[31] + b[58] - b[59] - b[62] + b[63]; 

      a[60] = - 4.0 *(b[42] - b[58])
              + 2.0 *(b[38] + b[41] + b[43] + b[46] - b[54] - b[57] - b[59] - b[62])
              -       b[37] - b[39] - b[45] - b[47] + b[53] + b[55] + b[61] + b[63];

      a[61] = - 2.0 *(b[42] - b[43] - b[58] + b[59])
              +       b[38] - b[39] + b[46] - b[47] - b[54] + b[55] - b[62] + b[63];
 
      a[62] = - 2.0 *(b[42] - b[46] - b[58] + b[62])
              +       b[41] + b[43] - b[45] - b[47] - b[57] - b[59] + b[61] + b[63];
 
      a[63] = -       b[42] + b[43] + b[46] - b[47] + b[58] - b[59] - b[62] + b[63]; 
}

//************************************************************************//

// performs explicit multiplication by Binv matrix

void BinvMultiply_LM(const double b[], double a[])

{

      a[0]  = b[0];
      a[1]  = b[8];
      a[2]  = -3.0*(b[0] - b[1]) - 2.0*b[8] - b[9];
      a[3]  =  2.0*(b[0] - b[1]) + b[8] + b[9];
      a[4]  = b[16];
      a[5]  = b[32];
      a[6]  = -3.0*(b[16] - b[17]) - 2.0*b[32] - b[33];
      a[7]  =  2.0*(b[16] - b[17]) + b[32] + b[33];
      
      a[8]  = -3.0*(b[0] - b[2])  - 2.0*b[16] - b[18];
      a[9]  = -3.0*(b[8] - b[10]) - 2.0*b[32] - b[34];
      a[10] = 9.0*(b[0] - b[1] - b[2] + b[3])
               + 6.0*(b[8] - b[10] + b[16] - b[17]) 
               + 3.0*(b[9] - b[11] + b[18] - b[19])
               + 4.0*b[32] + 2.0*(b[33] + b[34]) + b[35];
      
      a[11] = - 6.0*(b[0] - b[1] - b[2] + b[3])
              - 3.0*(b[8] + b[9] - b[10] - b[11])
              - 4.0*(b[16] - b[17]) - 2.0*(b[18] - b[19] + b[32] + b[33])
              - b[34] - b[35];

      a[12]  = 2.0*(b[0] - b[2])  + b[16] + b[18];
      a[13]  = 2.0*(b[8] - b[10]) + b[32] + b[34];

      a[14] = - 6.0*(b[0] - b[1] - b[2] + b[3])
              - 4.0*(b[8] - b[10]) - 2.0*(b[9] - b[11] + b[32] + b[34])
              - 3.0*(b[16] - b[17] + b[18] - b[19])
              - b[33] - b[35];

      a[15] = 4.0*(b[0] - b[1] - b[2] + b[3])
               + 2.0*(b[8] + b[9] - b[10] - b[11] + b[16] - b[17] + b[18] - b[19])
               + b[32] + b[33] + b[34] + b[35];
      ////////////////////

      a[16] = b[24];
      a[17] = b[40];
      a[18] = -3.0*(b[24] - b[25]) - 2.0*b[40] - b[41];
      a[19] =  2.0*(b[24] - b[25]) +     b[40] + b[41];
      a[20] = b[48];
      a[21] = b[56];
      a[22] = -3.0*(b[48] - b[49]) - 2.0*b[56] - b[57];
      a[23] =  2.0*(b[48] - b[49]) +   b[56] + b[57];

      a[24] = -3.0*(b[24] - b[26]) - 2.0*b[48] - b[50];
      a[25] = -3.0*(b[40] - b[42]) - 2.0*b[56] - b[58];
      a[26] = 9.0*(b[24] - b[25] - b[26] + b[27]) + 4.0*b[56]
              + 6.0*(b[40] - b[42] + b[48] - b[49]) 
              + 3.0*(b[41] - b[43] + b[50] - b[51])
              + 2.0*(b[57] + b[58]) + b[59];

      a[27] = - 6.0*(b[24] - b[25] - b[26] + b[27])
              - 4.0*(b[48] - b[49]) - 3.0*(b[40] + b[41] - b[42] - b[43])
              - 2.0*(b[50] - b[51] + b[56] + b[57])- b[58] - b[59];

      a[28] = 2.0*(b[24] - b[26]) + b[48] + b[50];
      a[29] = 2.0*(b[40] - b[42]) + b[56] + b[58];
      
      a[30] = - 6.0*(b[24] - b[25] - b[26] + b[27])
                 - 4.0*(b[40] - b[42]) 
                 - 3.0*(b[48] - b[49] + b[50] - b[51])
                 - 2.0*(b[41] - b[43] + b[56] + b[58]) - b[57] - b[59];

      a[31] = 4.0*(b[24] - b[25] - b[26] + b[27])
               + 2.0*(b[40] + b[41] - b[42] - b[43] + b[48] - b[49] + b[50] - b[51])
               + b[56] + b[57] + b[58] + b[59];
      ////////////////////

      a[32] = - 3.0*(b[0] - b[4])  - 2.0*b[24] - b[28];
      a[33] = - 3.0*(b[8] - b[12]) - 2.0*b[40] - b[44];

      a[34] = 9.0*(b[0]  - b[1]  - b[4]  + b[5])
               + 6.0*(b[8] - b[12] + b[24] - b[25])  + 4.0*b[40] 
               + 3.0*(b[9] - b[13] + b[28] - b[29])
               + 2.0*(b[41] + b[44]) + b[45];
      
      a[35] = - 6.0*(b[0]  - b[1]  - b[4]  + b[5])
                 - 3.0*(b[8]  + b[9]  - b[12] - b[13])
                 - 4.0*(b[24] - b[25]) 
                 - 2.0*(b[28] - b[29] + b[40] + b[41]) - b[44] - b[45];

      a[36] = - 3.0*(b[16] - b[20]) - 2.0*b[48] - b[52];
      a[37] = - 3.0*(b[32] - b[36]) - 2.0*b[56] - b[60];

      a[38] = 9.0*(b[16] - b[17] - b[20] + b[21])
               + 6.0*(b[32] - b[36] + b[48] - b[49]) + 4.0*b[56] 
               + 3.0*(b[33] - b[37] + b[52] - b[53])
               + 2.0*(b[57] + b[60]) + b[61];

      a[39] = - 6.0*(b[16] - b[17] - b[20] + b[21])
                 - 3.0*(b[32] + b[33] - b[36] - b[37])
                 - 4.0*(b[48] - b[49]) 
                 - 2.0*(b[52] - b[53] + b[56] + b[57]) - b[60] - b[61];

      a[40] = 9.0*(b[0]  - b[2]  - b[4]  + b[6])
               + 6.0*(b[16] - b[20] + b[24] - b[26]) + 4.0*b[48] 
               + 3.0*(b[18] - b[22] + b[28] - b[30])
               + 2.0*(b[50] + b[52]) + b[54];

      a[41] = 9.0*(b[8]  - b[10] - b[12] + b[14])
               + 6.0*(b[32] - b[36] + b[40] - b[42]) + 4.0*b[56] 
               + 3.0*(b[34] - b[38] + b[44] - b[46])
               + 2.0*(b[58] + b[60]) + b[62];

      a[42] = -27.0*(b[0] - b[1] - b[2] + b[3] - b[4] + b[5] + b[6] - b[7])
              -18.0*(b[8] - b[10] - b[12] + b[14] + b[16] - b[17] - b[20]
                   + b[21] + b[24] - b[25] - b[26] + b[27])
              -12.0*(b[32] - b[36] + b[40] - b[42] + b[48] - b[49])
              -9.0 *(b[9] - b[11] - b[13] + b[15] + b[18] - b[19] - b[22]
                   + b[23] + b[28] - b[29] - b[30] + b[31])
              -8.0 * b[56]
              -6.0 *(b[33] + b[34] - b[37] - b[38] + b[41] - b[43] + b[44] 
                   - b[46] + b[50] - b[51] + b[52] - b[53])
              -4.0 *(b[57] + b[58] + b[60])
              -3.0 *(b[35] - b[39] + b[45] - b[47] + b[54] - b[55])
              -2.0 *(b[59] + b[61] + b[62]) - b[63];

      a[43] = 18.0*(b[0] - b[1] - b[2] + b[3] - b[4] + b[5] + b[6] - b[7])
               + 9.0*(b[8] + b[9] - b[10] - b[11] - b[12] - b[13] + b[14] + b[15])
               + 12.0*(b[16] - b[17] - b[20] + b[21] + b[24] - b[25] - b[26] + b[27])
               + 6.0 *(b[18] - b[19] - b[22] + b[23] + b[28] - b[29] - b[30] + b[31] 
               + b[32] + b[33]  - b[36] - b[37] + b[40] + b[41] - b[42] - b[43]) 
               + 3.0*(b[34] + b[35] - b[38] - b[39] + b[44] + b[45] - b[46] - b[47])
               + 8.0*(b[48] - b[49]) + 4.0*(b[50] - b[51] + b[52] - b[53] + b[56] + b[57]) 
               + 2.0*(b[54] - b[55] + b[58] + b[59] + b[60] + b[61]) + b[62] + b[63];

      a[44] = - 6.0*(b[0]  - b[2]  - b[4]  + b[6])
                 - 3.0*(b[16] + b[18] - b[20] - b[22])
                 - 4.0*(b[24] - b[26]) - 2.0*(b[28] - b[30] + b[48] + b[50]) - b[52] - b[54];

      a[45] = - 6.0*(b[8]  - b[10] - b[12] + b[14])
                 - 3.0*(b[32] + b[34] - b[36] - b[38])
                 - 4.0*(b[40] - b[42]) - 2.0*(b[44] - b[46] + b[56] + b[58]) - b[60] - b[62];

      a[46] = 18.0*(b[0] - b[1] - b[2] + b[3] - b[4] + b[5] + b[6] - b[7])
               + 12.0*(b[8] - b[10] - b[12] + b[14] + b[24] - b[25] - b[26] + b[27]) 
               + 6.0*(b[9]  - b[11] - b[13] + b[15] + b[28] - b[29] - b[30] + b[31] 
               + b[32] + b[34] - b[36] - b[38] + b[48] - b[49] + b[50] - b[51])
               + 9.0*(b[16] - b[17] + b[18] - b[19] - b[20] + b[21] - b[22] + b[23])
               + 3.0*(b[33] + b[35] - b[37] - b[39] + b[52] - b[53] + b[54] - b[55])
               + 8.0*(b[40] - b[42]) + 4.0*(b[41] - b[43] + b[44] - b[46] + b[56] + b[58]) 
               + 2.0*(b[45] - b[47] + b[57]  + b[59] + b[60] + b[62]) + b[61] + b[63];

      a[47] = - 12.0*(b[0] - b[1] - b[2] + b[3] - b[4] + b[5] + b[6] - b[7])
                 - 6.0*(b[8] + b[9] - b[10] - b[11] - b[12] - b[13] + b[14] + b[15]
                 + b[16] - b[17] + b[18] - b[19] - b[20] + b[21] - b[22] + b[23])
                 - 8.0*(b[24] - b[25] - b[26] + b[27]) - 4.0*(b[28] - b[29] - b[30] + b[31]
                 + b[40] + b[41] - b[42] - b[43] + b[48] - b[49] + b[50] - b[51])
                 - 3.0*(b[32] + b[33] + b[34] + b[35] - b[36] - b[37] - b[38] - b[39])
                 - 2.0*(b[44] + b[45] - b[46] - b[47] + b[52] - b[53] + b[54] - b[55]
                 +      b[56] + b[57] + b[58] + b[59]) - b[60] - b[61] - b[62] - b[63];

      ////////////////////

      a[48] = 2.0*(b[0] - b[4])  + b[24] + b[28];
      a[49] = 2.0*(b[8] - b[12]) + b[40] + b[44];

      a[50] = - 6.0*(b[0] - b[1] - b[4] + b[5])
                 - 4.0*(b[8] - b[12]) - 3.0*(b[24] - b[25] + b[28] - b[29])
                 - 2.0*(b[9]  - b[13] + b[40] + b[44]) - b[41] - b[45];

      a[51] = 4.0*(b[0] - b[1] - b[4] + b[5])
               + 2.0*(b[8] + b[9] - b[12] - b[13] + b[24] - b[25] + b[28] - b[29])
               + b[40] + b[41] + b[44] + b[45];
        
      a[52] = 2.0*(b[16] - b[20]) + b[48] + b[52];
      a[53] = 2.0*(b[32] - b[36]) + b[56] + b[60];
        
      a[54] = - 6.0*(b[16] - b[17] - b[20] + b[21])
                 - 4.0*(b[32] - b[36]) - 3.0*(b[48] - b[49] + b[52] - b[53])
                 - 2.0*(b[33] - b[37] + b[56] + b[60]) - b[57] - b[61];

      a[55] = 4.0*(b[16] - b[17] - b[20] + b[21])
               + 2.0*(b[32] + b[33] - b[36] - b[37] + b[48] - b[49] + b[52] - b[53])
               + b[56] + b[57] + b[60] + b[61];

      a[56] = - 6.0*(b[0]  - b[2]  - b[4]  + b[6])
                 - 4.0*(b[16] - b[20]) - 3.0*(b[24] - b[26] + b[28] - b[30])
                 - 2.0*(b[18] - b[22] + b[48] + b[52]) - b[50] - b[54];

      a[57] = - 6.0*(b[8]  - b[10] - b[12] + b[14])
                 - 4.0*(b[32] - b[36]) - 3.0*(b[40] - b[42] + b[44] - b[46])
                 - 2.0*(b[34] - b[38] + b[56] + b[60]) - b[58]  - b[62];

      a[58] = 18.0*(b[0] - b[1] - b[2] + b[3] - b[4] + b[5] + b[6] - b[7])
               + 12.0*(b[8] - b[10] - b[12] + b[14] + b[16] - b[17] - b[20] + b[21])
               + 6.0 *(b[9] - b[11] - b[13] + b[15] + b[18] - b[19] - b[22] + b[23]
                     + b[40] - b[42] + b[44] - b[46] + b[48] - b[49] + b[52]  - b[53])
               + 9.0* (b[24]  - b[25]  - b[26]  + b[27] + b[28]  - b[29]  - b[30]  + b[31])
               + 8.0* (b[32]  - b[36]) + 4.0*(b[33]  + b[34]  - b[37]  - b[38] + b[56] + b[60])  
               + 3.0* (b[41] - b[43] + b[45]  - b[47] + b[50]  - b[51] + b[54]  - b[55])
               + 2.0* (b[35] - b[39] + b[57] + b[58] + b[61]  + b[62])  + b[59] + b[63];

      a[59] = - 12.0*(b[0] - b[1] - b[2] + b[3] - b[4] + b[5] + b[6] - b[7])
                 - 6.0*(b[8] + b[9] - b[10] - b[11] - b[12] - b[13] + b[14] + b[15]
                 +      b[24] - b[25] - b[26] + b[27] + b[28] - b[29] - b[30] + b[31])
                 - 8.0*(b[16] - b[17] - b[20] + b[21]) - 4.0*(b[18] - b[19] - b[22] + b[23] 
                 +      b[32] + b[33] - b[36] - b[37]  + b[48] - b[49] + b[52] - b[53])
                 - 3.0*(b[40] + b[41] - b[42] - b[43] +  b[44] + b[45] - b[46] - b[47])
                 - 2.0*(b[34] + b[35] - b[38] - b[39] +  b[50] - b[51] + b[54] - b[55]
                 +      b[56] + b[57] + b[60] + b[61]) - b[58] - b[59] - b[62] - b[63];

      a[60] = 4.0*(b[0] - b[2] - b[4] + b[6]) + 2.0*(b[16] + b[18] - b[20] - b[22]
               +   b[24] - b[26] + b[28] - b[30]) + b[48] + b[50] + b[52] + b[54];

      a[61] = 4.0*(b[8] - b[10] - b[12] + b[14]) + 2.0*(b[32] + b[34] - b[36] - b[38]
               + b[40] - b[42] + b[44] - b[46]) + b[56] + b[58] + b[60] + b[62];

      a[62] = - 12.0*(b[0] - b[1] - b[2] + b[3] - b[4] + b[5] + b[6] - b[7])
                 - 8.0*(b[8] - b[10] - b[12] + b[14])  - 4.0*(b[9] - b[11] - b[13] + b[15] 
                 +      b[32] + b[34] - b[36] - b[38] + b[40] - b[42] + b[44] - b[46]) 
                 - 6.0*(b[16] - b[17] + b[18] - b[19] - b[20] + b[21] - b[22] + b[23]
                 +      b[24] - b[25] - b[26] + b[27] + b[28] - b[29] - b[30] + b[31])
                 - 3.0*(b[48] - b[49] + b[50] - b[51] + b[52] - b[53] + b[54] - b[55])
                 - 2.0*(b[33] + b[35] - b[37] - b[39] + b[41] - b[43] + b[45] - b[47]
                 +      b[56] + b[58] + b[60] + b[62]) - b[57] - b[59] - b[61] - b[63];

      a[63] = 8.0*(b[0] - b[1] - b[2] + b[3] - b[4] + b[5] + b[6] - b[7])
               + 4.0*(b[8] + b[9] - b[10] - b[11] - b[12] - b[13] + b[14] + b[15]
               + b[16] - b[17] + b[18] - b[19] - b[20] + b[21] - b[22] + b[23]
               + b[24] - b[25] - b[26] + b[27] + b[28] - b[29] - b[30] + b[31])
               + 2.0*(b[32] + b[33] + b[34] + b[35] - b[36] - b[37] - b[38] - b[39]
               + b[40] + b[41] - b[42] - b[43] + b[44] + b[45] - b[46] - b[47]
               + b[48] - b[49] + b[50] - b[51] + b[52] - b[53] + b[54] - b[55])
               + b[56] + b[57] + b[58] + b[59] + b[60] + b[61] + b[62] + b[63];
}

//*****************************************************************************//

// performs explicit multiplication by Transpose(Binv6) matrix

void BinvMultiply_Transpose_B6(const double b[], double a[])
{
 
 a[0] = - 708.0*(b[43] + b[46] + b[58]) - 605.0*(b[26] + b[38] + b[41]) - 208.0*b[63]
        - 180.0*(b[31] + b[55] + b[61]) - 150.0*(b[23] + b[29] + b[53]) - 125.0*b[21]
        - 6.0  *(b[11] + b[14] + b[35] + b[44] + b[50] + b[56]) - 3.0*(b[2] + b[8] + b[32])
        + b[0] + 2.0*(b[3] + b[12] + b[48]) + 4.0*(b[15] + b[51] + b[60]) 
        + 9.0*(b[10] + b[34] + b[40]) + 275.0*(b[22] + b[25] + b[37]) 
        + 330.0*(b[27] + b[30] + b[39] + b[45] + b[54] + b[57]) + 384.0*(b[47] + b[59] + b[62])
        + 1304.0*b[42];

 a[1] = - 820.0*b[42] - 384.0*(b[47] + b[59]) - 330.0*(b[27] + b[39]) - 240.0*b[62] 
        - 210.0*(b[30] + b[54]) - 175.0*b[22] - 66.0*(b[45] + b[57]) - 55.0*(b[25] + b[37])
        -   9.0*(b[10] + b[34]) - 4.0*(b[15] + b[51]) - 2.0*b[3] + 3.0*b[2] + 6.0*(b[11] + b[14] + b[35] + b[50])
        +  25.0*b[21] + 30.0*(b[29] + b[53]) + 36.0*b[61] + 121.0*b[41] + 150.0*b[23]
        + 180.0*(b[31] + b[55]) + 208.0*b[63] + 385.0*(b[26] + b[38]) + 444.0*(b[46] + b[58])
        + 708.0*b[43]; 

 a[2] = - 820.0*b[42] - 384.0*(b[47] + b[62]) - 330.0*(b[30] + b[45]) - 240.0*b[59]   
        - 210.0*(b[27] + b[57]) - 175.0*b[25] - 66.0*(b[39] + b[54]) - 55.0*(b[22] + b[37])
        -   9.0*(b[10] + b[40]) - 4.0*(b[15] + b[60]) - 2.0*b[12] + 3.0*b[8] + 6.0*(b[11] + b[14] + b[44] + b[56])
        +  25.0*b[21] + 30.0*(b[23] + b[53]) + 36.0*b[55] + 121.0*b[38] + 150.0*b[29]
        + 180.0*(b[31] + b[61]) + 208.0*b[63] + 385.0*(b[26] + b[41]) + 444.0*(b[43] + b[58])
        + 708.0*b[46];

 a[3] = - 444.0*(b[43] + b[46]) - 276.0*b[58] - 245.0*b[26] - 208.0*b[63] - 180.0*b[31] 
        -  77.0*(b[38] + b[41]) - 36.0*(b[55] + b[61]) - 30.0*(b[23] + b[29]) - 6.0*(b[11] + b[14] + b[53])
        -   5.0*b[21] + 4.0*b[15] + 9.0*b[10] + 11.0*b[37] + 35.0*(b[22] + b[25]) + 42.0*(b[54] + b[57])
        +  66.0*(b[39] + b[45]) + 210.0*(b[27] + b[30]) + 240.0*(b[59]+b[62]) + 384.0*b[47] + 512.0*b[42];

 a[4] = - 820.0*b[42] - 384.0*(b[59] + b[62]) - 330.0*(b[54] + b[57]) - 240.0*b[47]   
        - 210.0*(b[39] + b[45]) - 175.0*b[37] - 66.0*(b[27] + b[30]) - 55.0*(b[22] + b[25])
        -   9.0*(b[34] + b[40]) - 4.0*(b[51] + b[60]) - 2.0*b[48] + 3.0*b[32] + 6.0*(b[35] + b[44] + b[50] + b[56])
        +  25.0*b[21] + 30.0*(b[23] + b[29]) + 36.0*b[31] + 121.0*b[26] + 150.0*b[53]
        + 180.0*(b[55] + b[61]) + 208.0*b[63] + 385.0*(b[38] + b[41]) + 444.0*(b[43] + b[46])
        + 708.0*b[58];

 a[5] = - 444.0*(b[43] + b[58]) - 276.0*b[46] - 245.0*b[38] - 208.0*b[63] - 180.0*b[55] 
        -  77.0*(b[26] + b[41]) - 36.0*(b[31] + b[61]) - 30.0*(b[23] + b[53]) - 6.0*(b[29] + b[35] + b[50])
        -   5.0*b[21] + 4.0*b[51] + 9.0*b[34] + 11.0*b[25] + 35.0*(b[22] + b[37]) + 42.0*(b[30] + b[45])
        +  66.0*(b[27] + b[57]) + 210.0*(b[39] + b[54]) + 240.0*(b[47]+b[62]) + 384.0*b[59] + 512.0*b[42];

 a[6] = - 444.0*(b[46] + b[58]) - 276.0*b[43] - 245.0*b[41] - 208.0*b[63] - 180.0*b[61] 
        -  77.0*(b[26] + b[38]) - 36.0*(b[31] + b[55]) - 30.0*(b[29] + b[53]) - 6.0*(b[23] + b[44] + b[56])
        -   5.0*b[21] + 4.0*b[60] + 9.0*b[40] + 11.0*b[22] + 35.0*(b[25] + b[37]) + 42.0*(b[27] + b[39])
        +  66.0*(b[30] + b[54]) + 210.0*(b[45] + b[57]) + 240.0*(b[47]+b[59]) + 384.0*b[62] + 512.0*b[42];

 a[7] = - 316.0*b[42] - 240.0*(b[47] + b[59] + b[62]) - 42.0*(b[27] + b[30] + b[39] + b[45] + b[54] + b[57])
        -   7.0*(b[22] + b[25] + b[37]) + b[21] + 6.0*(b[23] + b[29] + b[53]) + 36.0*(b[31] + b[55] + b[61])
        +  49.0*(b[26] + b[38] + b[41]) + 208.0*b[63] + 276.0*(b[43] + b[46] + b[58]);

 a[8] = - 64.0*(b[22] + b[25] + b[37] + b[42] - b[21] - b[26] - b[38] - b[41]);

 a[9] = - 120.0*(b[46] + b[58]) - 112.0*(b[41] + b[43]) - 110.0*(b[26] + b[38]) - 32.0*(b[61] + b[63])
        -  30.0*(b[29] + b[31] + b[53] + b[55]) - 25.0*(b[21] + b[23]) - 4.0*(b[14] + b[50])
        -   3.0*(b[9] + b[11] + b[33] + b[35]) - 2.0*(b[2] - b[13] - b[15] - b[49] - b[51]) + b[1] + b[3] 
        +   6.0*(b[10] + b[34]) + 50.0*b[22] + 55.0*(b[25] + b[27] + b[37] + b[39]) 
        +  60.0*(b[30] + b[45] + b[47] + b[54] + b[57] + b[59]) + 64.0*b[62] + 224.0*b[42];

 a[10] = - 112.0*(b[43] - b[42]) - 60.0*(b[46] + b[58] - b[47] - b[59]) - 55.0*(b[26] + b[38] - b[27] - b[39])
         -  32.0*(b[63] - b[62]) - 30.0*(b[31] + b[55] - b[30] - b[54]) - 25.0*(b[23] - b[22])
         -   3.0*(b[11] + b[35] - b[10] - b[34]) - 2.0*(b[14] + b[50] - b[15] - b[51]) - b[2] + b[3];

 a[11] = - 136.0*b[42] - 64.0*b[62] - 60.0*(b[30] + b[45] + b[47]) - 36.0*(b[57] + b[59]) - 35.0*(b[25] + b[27])
         -  12.0*b[54] - 11.0*(b[37] + b[39]) - 10.0*b[22] - 6.0*b[10] - 2.0*(b[13] + b[15]) + 3.0*(b[9] + b[11])
         +   4.0*b[14] + 5.0*(b[21] + b[23]) + 6.0*(b[53] + b[55]) + 22.0*b[38] + 30.0*(b[29] + b[31])
         +  32.0*(b[61] + b[63]) + 68.0*(b[41] + b[43]) + 70.0*b[26] + 72.0*b[58] + 120.0*b[46];

 a[12] = - 68.0*(b[42] - b[43]) - 60.0*(b[47] - b[46]) - 36.0*(b[59] - b[58]) - 35.0*(b[27] - b[26])
         - 32.0*(b[62] - b[63]) - 30.0*(b[30] - b[31]) - 11.0*(b[39] - b[38]) - 6.0*(b[54] - b[55])
         -  5.0*(b[22] - b[23]) -  3.0*(b[10] - b[11]) - 2.0 *(b[15] - b[14]);

 a[13] = -  2.0*(b[49] + b[51]) + 3.0*(b[33] + b[35]) + 4.0*b[50] + 5.0*(b[21] + b[23]) + 6.0*(b[29] + b[31] - b[34])
         - 10.0*b[22] - 11.0*(b[25] + b[27]) - 12.0*b[30] + 22.0*b[26] + 30.0*(b[53] + b[55]) + 32.0*(b[61] + b[63])
         - 35.0*(b[37] + b[39]) - 36.0*(b[45] + b[47]) - 60.0*(b[54] + b[57] + b[59]) - 64.0*b[62] + 68.0*(b[41] + b[43])
         + 70.0*b[38] + 72.0*b[46] + 120.0*b[58] - 136.0*b[42];
 
 a[14] = - 68.0*(b[42] - b[43]) - 60.0*(b[59] - b[58]) - 36.0*(b[47] - b[46]) - 35.0*(b[39] - b[38])
         - 32.0*(b[62] - b[63]) - 30.0*(b[54] - b[55]) - 11.0*(b[27] - b[26]) - 6.0*(b[30] - b[31])
         -  5.0*(b[22] - b[23]) -  3.0*(b[34] - b[35]) - 2.0 *(b[51] - b[50]);

 a[15] = - 72.0*(b[46] + b[58]) - 40.0*(b[41] + b[43]) - 32.0*(b[61] + b[63]) - 14.0*(b[26] + b[38])
         -  6.0*(b[29] + b[31] + b[53] + b[55]) - b[21] - b[23] + 2.0*b[22] + 7.0*(b[25] + b[27] + b[37] + b[39])
         + 12.0*(b[30] + b[54]) + 36.0*(b[45] + b[47] + b[57] + b[59]) + 64.0*b[62] + 80.0*b[42];

 a[16] = - 40.0*(b[43] - b[42]) - 36.0*(b[46] + b[58] -b[47] - b[59]) - 32.0*(b[63] - b[62])
         -  7.0*(b[26] + b[38] - b[27] - b[39]) - 6.0*(b[31] + b[55] - b[30] - b[54]) - b[23] + b[22];

 a[17] = - 96.0*(b[26] + b[38] - b[22] - b[42]) - 64.0*(b[23] + b[43] - b[27] - b[39]) 
         - 32.0*(b[21] + b[41] - b[25] - b[37]);

 a[18] = - 120.0*(b[43] + b[58]) - 112.0*(b[38] + b[46]) - 110.0*(b[26] + b[41]) - 32.0*(b[55] + b[63])
         -  30.0*(b[23] + b[31] + b[53] + b[61]) - 25.0*(b[21] + b[29]) - 4.0*(b[11] + b[56])
         -   3.0*(b[6] + b[14] + b[36] + b[44]) - 2.0*(b[8] - b[7] - b[15] - b[52] - b[60]) + b[4] + b[12] 
         +   6.0*(b[10] + b[40]) + 50.0*b[25] + 55.0*(b[22] + b[30] + b[37] + b[45]) 
         +  60.0*(b[27] + b[39] + b[47] + b[54] + b[57] + b[62]) + 64.0*b[59] + 224.0*b[42];

 a[19] = - 136.0*b[42] - 64.0*b[59] - 60.0*(b[27] + b[39] + b[47]) - 36.0*(b[54] + b[62]) - 35.0*(b[22] + b[30])
         -  12.0*b[57] - 11.0*(b[37] + b[45]) - 10.0*b[25] - 6.0*b[10] - 2.0*(b[7] + b[15]) + 3.0*(b[6] + b[14])
         +   4.0*b[11] + 5.0*(b[21] + b[29]) + 6.0*(b[53] + b[61]) + 22.0*b[41] + 30.0*(b[23] + b[31])
         +  32.0*(b[55] + b[63]) + 68.0*(b[38] + b[46]) + 70.0*b[26] + 72.0*b[58] + 120.0*b[43];

 a[20] = - 112.0*(b[46] - b[42]) - 60.0*(b[43] + b[58] - b[47] - b[62]) - 55.0*(b[26] + b[41] - b[30] - b[45])
         -  32.0*(b[63] - b[59]) - 30.0*(b[31] + b[61] - b[27] - b[57]) - 25.0*(b[29] - b[25])
         -   3.0*(b[14] + b[44] - b[10] - b[40]) - 2.0*(b[11] + b[56] - b[15] - b[60]) - b[8] + b[12];

 a[21] = - 68.0*(b[42] - b[46]) - 60.0*(b[47] - b[43]) - 36.0*(b[62] - b[58]) - 35.0*(b[30] - b[26])
         - 32.0*(b[59] - b[63]) - 30.0*(b[27] - b[31]) - 11.0*(b[45] - b[41]) - 6.0*(b[57] - b[61])
         -  5.0*(b[25] - b[29]) -  3.0*(b[10] - b[14]) - 2.0 *(b[15] - b[11]);

 a[22] = - 136.0*b[42] - 64.0*b[59] - 60.0*(b[54] + b[57] + b[62]) - 36.0*(b[39] + b[47]) - 35.0*(b[37] + b[45])
         -  12.0*b[27] - 11.0*(b[22] + b[30]) - 10.0*b[25] - 6.0*b[40] - 2.0*(b[52] + b[60]) + 3.0*(b[36] + b[44])
         +   4.0*b[56] + 5.0*(b[21] + b[29]) + 6.0*(b[23] + b[31]) + 22.0*b[26] + 30.0*(b[53] + b[61])
         +  32.0*(b[55] + b[63]) + 68.0*(b[38] + b[46]) + 70.0*b[41] + 72.0*b[43] + 120.0*b[58];

 a[23] = - 72.0*(b[43] + b[58]) - 40.0*(b[38] + b[46]) - 32.0*(b[55] + b[63]) - 14.0*(b[26] + b[41])
         -  6.0*(b[23] + b[31] + b[53] + b[61]) - b[21] - b[29] + 2.0*b[25] + 7.0*(b[22] + b[30] + b[37] + b[45])
         + 12.0*(b[27] + b[57]) + 36.0*(b[39] + b[47] + b[54] + b[62]) + 64.0*b[59] + 80.0*b[42];

 a[24] = - 68.0*(b[42] - b[46]) - 60.0*(b[62] - b[58]) - 36.0*(b[47] - b[43]) - 35.0*(b[45] - b[41])
         - 32.0*(b[59] - b[63]) - 30.0*(b[57] - b[61]) - 11.0*(b[30] - b[26]) - 6.0*(b[27] - b[31])
         -  5.0*(b[25] - b[29]) -  3.0*(b[40] - b[44]) - 2.0 *(b[60] - b[56]);

 a[25] =         b[25] - b[29] + 6.0*(b[27] + b[57] - b[31] - b[61]) + 7.0*(b[30] + b[45] - b[26] - b[41])
         + 32.0*(b[59] - b[63]) + 36.0*(b[47] + b[62] - b[43]  - b[58]) + 40.0*(b[42] - b[46]);

 a[26] =   32.0*(b[22] + b[37] - b[21] - b[38]) + 64.0*(b[30] + b[45] - b[29] - b[46]) 
         + 96.0*(b[25] + b[42] - b[26] - b[41]);

 a[27] =         b[16] + b[48] + 2.0*(b[19] + b[28] + b[51] + b[60] - b[32]) - 3.0*(b[18] + b[24] + b[50] + b[56])
         -  4.0*(b[35] + b[44]) + 6.0*(b[34] + b[40]) - 25.0*(b[21] + b[53]) - 30.0*(b[23] + b[29] + b[55] + b[61])
         - 32.0*(b[31] + b[63]) + 50.0*b[37] + 55.0*(b[22] + b[25] + b[54] + b[57]) 
         + 60.0*(b[27] + b[30] + b[39] + b[45] + b[59] + b[62]) + 64.0*b[47] - 110.0*(b[38] + b[41])
         - 112.0*(b[26] + b[58]) - 120.0*(b[43] + b[46]) + 224.0*b[42];

 a[28] = -   2.0*(b[19] + b[51]) + 3.0*(b[18] + b[50]) + 4.0*b[35] + 5.0*(b[21] + b[53]) + 6.0*(b[29] + b[61] - b[34])
         -  10.0*b[37] - 11.0*(b[25] + b[57]) - 12.0*b[45] + 22.0*b[41] + 30.0*(b[23] + b[55]) + 32.0*(b[31] + b[63])
         -  35.0*(b[22] + b[54]) - 36.0*(b[30] + b[62]) - 60.0*(b[27] + b[39] + b[59]) - 64.0*b[47] + 68.0*(b[26] + b[58])
         +  70.0*b[38] + 72.0*b[46] + 120.0*b[43] - 136.0*b[42];

 a[29] = -   2.0*(b[28] + b[60]) + 3.0*(b[24] + b[56]) + 4.0*b[44] + 5.0*(b[21] + b[53]) + 6.0*(b[23] + b[55] - b[40])
         -  10.0*b[37] - 11.0*(b[22] + b[54]) - 12.0*b[39] + 22.0*b[38] + 30.0*(b[29] + b[61]) + 32.0*(b[31] + b[63])
         -  35.0*(b[25] + b[57]) - 36.0*(b[27] + b[59]) - 60.0*(b[30] + b[45] + b[62]) - 64.0*b[47] + 68.0*(b[26] + b[58])
         +  70.0*b[41] + 72.0*b[43] + 120.0*b[46] - 136.0*b[42];

 a[30] = -  b[21] - b[53] + 2.0*b[37] - 6.0*(b[23] + b[29] + b[55] + b[61]) + 7.0*(b[22] + b[25] + b[54] + b[57])
         +  12.0*(b[39] + b[45]) - 14.0*(b[38] + b[41]) - 32.0*(b[31] + b[63]) + 36.0*(b[27] + b[30] + b[59] + b[62])
         -  40.0*(b[26] + b[58]) + 64.0*b[47] - 72.0*(b[43] + b[46]) + 80.0*b[42];

 a[31] =    b[48] - b[32] + 2.0*(b[51] + b[60] - b[35] - b[44]) + 3.0*(b[34] + b[40] - b[50] - b[56]) + 25.0*(b[37] - b[53])
         +  30.0*(b[39] + b[45] - b[55] - b[61]) + 32.0*(b[47] - b[63]) + 55.0*(b[54] + b[57] - b[38] - b[41])
         +  60.0*(b[59] + b[62] - b[43] - b[46]) + 112.0*(b[42] - b[58]);

 a[32] =     2.0*(b[35] - b[51]) + 3.0*(b[50] - b[34]) + 5.0*(b[53] - b[37]) + 6.0*(b[61] - b[45]) + 11.0*(b[41] - b[57])
         +  30.0*(b[55] - b[39]) + 32.0*(b[63] - b[47]) + 35.0*(b[38] - b[54]) + 36.0*(b[46] - b[62]) + 60.0*(b[43] - b[59])
         +  68.0*(b[58] - b[42]);

 a[33] =     2.0*(b[44] - b[60]) + 3.0*(b[56] - b[40]) + 5.0*(b[53] - b[37]) + 6.0*(b[55] - b[39]) + 11.0*(b[38] - b[54])
         +  30.0*(b[61] - b[45]) + 32.0*(b[63] - b[47]) + 35.0*(b[41] - b[57]) + 36.0*(b[43] - b[59]) + 60.0*(b[46] - b[62])
         +  68.0*(b[58] - b[42]);

 a[34] =    b[37] - b[53] + 6.0*(b[39] + b[45] - b[55] - b[61]) + 7.0*(b[54] + b[57] - b[38] - b[41]) + 32.0*(b[47] - b[63])
         +  36.0*(b[59] + b[62] - b[43] - b[46]) + 40.0*(b[42] - b[58]);

 a[35] =    32.0*(b[22] + b[25] - b[21] - b[26]) + 64.0*(b[54] + b[57] - b[53] - b[58]) + 96.0*(b[37] + b[42] - b[38] - b[41]);

 a[36] =    b[5] + b[7] + b[13] + b[15] - 2.0*(b[6] + b[9] + b[11] + b[14]) + 4.0*(b[10] - b[53] - b[55] - b[61] - b[63])
         -   5.0*(b[21] + b[23] + b[29] + b[31]) + 8.0*(b[37] + b[39] + b[45] + b[47] + b[54] + b[57] + b[59] + b[62])
         +  10.0*(b[22] + b[25] + b[27] + b[30]) - 16.0*(b[38] + b[41] + b[43] + b[46] + b[58]) - 20.0*b[26] + 32.0*b[42];

 a[37] =    b[7] + b[15] - b[6] - b[14] + 2.0*(b[10] - b[11]) + 4.0*(b[54] + b[62] - b[55] - b[63]) 
         +   5.0*(b[22] + b[30] - b[23] - b[31]) + 8.0*(b[39] + b[47] + b[59] - b[38] - b[46] - b[58]) + 10.0*(b[27] - b[26])
         +  16.0*(b[42] - b[43]);

 a[38] =    b[13] + b[15] - b[9] - b[11] + 2.0*(b[10] - b[14]) + 4.0*(b[57] + b[59] - b[61] - b[63]) 
         +   5.0*(b[25] + b[27] - b[29] - b[31]) + 8.0*(b[45] + b[47] + b[62] - b[41] - b[43] - b[58]) + 10.0*(b[30] - b[26])
         +  16.0*(b[42] - b[46]);

 a[39] =    b[10] + b[15] - b[11] - b[14] + 4.0*(b[59] + b[62] - b[58] - b[63])
         +   5.0*(b[27] + b[30] - b[26] - b[31]) + 8.0*(b[42] + b[47] - b[43] - b[46]);

 a[40] =    b[21] + b[23] + b[29] + b[31] - 2.0*(b[22] + b[25] + b[27] + b[30])
         +    4.0*(b[26] + b[53] + b[55] + b[61] + b[63] - b[37] - b[39] - b[45] - b[47])
         +    8.0*(b[38] + b[41] + b[43] + b[46] - b[54] - b[57] - b[59] - b[62]) + 16.0*(b[58] - b[42]);

 a[41] =    b[23] + b[31] - b[22] - b[30] + 2.0*(b[26] - b[27]) 
         +    4.0*(b[38] + b[46] + b[55] + b[63] - b[39] - b[47] - b[54] - b[62]) + 8.0*(b[43] + b[58] - b[42] - b[59]);

 a[42] =   b[29] + b[31] - b[25] - b[27] + 2.0*(b[26] - b[30]) 
         +    4.0*(b[41] + b[43] + b[61] + b[63] - b[45] - b[47] - b[57] - b[59]) + 8.0*(b[46] + b[58] - b[42] - b[62]);

 a[43] =   b[26] + b[31] - b[27] - b[30] + 4.0*(b[43] + b[46] + b[58] + b[63] - b[42] - b[47] - b[59] - b[62]);

 a[44] =     16.0*(b[21] - b[37]) + 32.0*(b[23] + b[29] - b[39] - b[45]) + 48.0*(b[38] + b[41] - b[22] - b[25])
         +   64.0*(b[31] - b[47]) + 96.0*(b[43] + b[46] - b[27] - b[30]) + 144.0*(b[26] - b[42]);

 a[45] =   b[17] + b[19] + b[49] + b[51] - 2.0*(b[18] + b[33] + b[35] + b[50]) + 4.0*(b[34] - b[29] - b[31] - b[61] - b[63])
         -   5.0*(b[21] + b[23] + b[53] + b[55]) + 8.0*(b[25] + b[27] + b[30] + b[45] + b[47] + b[57] + b[59] + b[62])
         +  10.0*(b[22] + b[37] + b[39] + b[54]) - 16.0*(b[26] + b[41] + b[43] + b[46] + b[58]) - 20.0*b[38] + 32.0*b[42];

 a[46] =    b[19] + b[51] - b[18] - b[50] + 2.0*(b[34] - b[35]) + 4.0*(b[30] + b[62] - b[31] - b[63]) 
         +   5.0*(b[22] + b[54] - b[23] - b[55]) + 8.0*(b[27] + b[47] + b[59] - b[26] - b[46] - b[58]) + 10.0*(b[39] - b[38])
         +  16.0*(b[42] - b[43]);

 a[47] =    b[21] + b[23] + b[53] + b[55] - 2.0*(b[22] + b[37] + b[39] + b[54])                    
         +    4.0*(b[29] + b[31] + b[38] + b[61] + b[63] - b[25] - b[27] - b[57] - b[59])
         +    8.0*(b[26] + b[41] + b[43] + b[58] - b[30] - b[45] - b[47] - b[62]) + 16.0*(b[46] - b[42]);

 a[48] =   b[23] + b[55] - b[22] - b[54] + 2.0*(b[38] - b[39])
         +    4.0*(b[26] + b[31] + b[58] + b[63] - b[27] - b[30] - b[59] - b[62]) + 8.0*(b[43] + b[46] - b[42] - b[47]);

 a[49] =    b[49] + b[51] - b[33] - b[35] + 2.0*(b[34] - b[50]) + 4.0*(b[45] + b[47] - b[61] - b[63])
         +   5.0*(b[37] + b[39] - b[53] - b[55]) + 8.0*(b[57] + b[59] + b[62] - b[41] - b[43] - b[46]) + 10.0*(b[54] - b[38])
         +  16.0*(b[42] - b[58]);

 a[50] =   b[34] + b[51] - b[35] - b[50] + 4.0*(b[47] + b[62] - b[46] - b[63])
         +    5.0*(b[39] + b[54] - b[38] - b[55]) + 8.0*(b[42] + b[59] - b[43] - b[58]);

 a[51] =   b[53] + b[55] - b[37] - b[39] + 2.0*(b[38] - b[54])
         +    4.0*(b[41] + b[43] + b[61] + b[63] - b[45] - b[47] - b[57] - b[59]) + 8.0*(b[46] + b[58] - b[42] - b[62]);

 a[52] =  b[38] + b[55] - b[39] - b[54] + 4.0*(b[43] + b[46] + b[58] + b[63] - b[42] - b[47] - b[59] - b[62]);

 a[53] =     16.0*(b[21] - b[25]) + 32.0*(b[23] + b[53] - b[27] - b[57]) + 48.0*(b[26] + b[41] - b[22] - b[37])
         +   64.0*(b[55] - b[59]) + 96.0*(b[43] + b[58] - b[39] - b[54]) + 144.0*(b[38] - b[42]);
 
 a[54] =   b[20] + b[28] + b[52] + b[60] - 2.0*(b[24] + b[36] + b[44] + b[56]) + 4.0*(b[40] - b[23] - b[31] - b[55] - b[63])
         -   5.0*(b[21] + b[29] + b[53] + b[61]) + 8.0*(b[22] + b[27] + b[30] + b[39] + b[47] + b[54] + b[59] + b[62])
         +  10.0*(b[25] + b[37] + b[45] + b[57]) - 16.0*(b[26] + b[38] + b[43] + b[46] + b[58]) - 20.0*b[41] + 32.0*b[42];

 a[55] =    b[21] + b[29] + b[53] + b[61] - 2.0*(b[25] + b[37] + b[45] + b[57])
         +    4.0*(b[23] + b[31] + b[41] + b[55] + b[63] - b[22] - b[30] - b[54] - b[62])
         +    8.0*(b[26] + b[38] + b[46] + b[58] - b[27] - b[39] - b[47] - b[59]) + 16.0*(b[43] - b[42]);

 a[56] =    b[28] + b[60] - b[24] - b[56] + 2.0*(b[40] - b[44]) + 4.0*(b[27] + b[59] - b[31] - b[63])
         +   5.0*(b[25] + b[57] - b[29] - b[61]) + 8.0*(b[30] + b[47] + b[62] - b[26] - b[43] - b[58]) + 10.0*(b[45] - b[41])
         +  16.0*(b[42] - b[46]);
 
 a[57] =   b[29] + b[61] - b[25] - b[57] + 2.0*(b[41] - b[45])
         +    4.0*(b[26] + b[31] + b[58] + b[63] - b[27] - b[30] - b[59] - b[62]) + 8.0*(b[43] + b[46] - b[42] - b[47]);

 a[58] =    b[52] + b[60] - b[36] - b[44] + 2.0*(b[40] - b[56]) + 4.0*(b[39] + b[47] - b[55] - b[63])
         +   5.0*(b[37] + b[45] - b[53] - b[61]) + 8.0*(b[54] + b[59] + b[62] - b[38] - b[43] - b[46]) + 10.0*(b[57] - b[41])
         +  16.0*(b[42] - b[58]);

 a[59] =   b[53] + b[61] - b[37] - b[45] + 2.0*(b[41] - b[57])
         +    4.0*(b[38] + b[46] + b[55] + b[63] - b[39] - b[47] - b[54] - b[62]) + 8.0*(b[43] + b[58] - b[42] - b[59]);

 a[60] =   b[40] + b[60] - b[44] - b[56] + 4.0*(b[47] + b[59] - b[43] - b[63])
         +    5.0*(b[45] + b[57] - b[41] - b[61]) + 8.0*(b[42] + b[62] - b[46] - b[58]);

 a[61] =  b[41] + b[61] - b[45] - b[57] + 4.0*(b[43] + b[46] + b[58] + b[63] - b[42] - b[47] - b[59] - b[62]);

 a[62] =     16.0*(b[21] - b[22]) + 32.0*(b[29] + b[53] - b[30] - b[54]) + 48.0*(b[26] + b[38] - b[25] - b[37])
         +   64.0*(b[61] - b[62]) + 96.0*(b[46] + b[58] - b[45] - b[57]) + 144.0*(b[41] - b[42]);

 a[63] = -   8.0*b[21] - 16.0*(b[23] + b[29] + b[53]) + 24.0*(b[22] + b[25] + b[37]) - 32.0*(b[31] + b[55] + b[61])
         +  48.0*(b[27] + b[30] + b[39] + b[45] + b[54] + b[57]) - 64.0*b[63] - 72.0*(b[26] + b[38] + b[41])
         +  96.0*(b[47] + b[59] + b[62]) -144.0*(b[43] + b[46] + b[58]) + 216.0*b[42];

}

//*************************************************************************//

// performs explicit multiplication by Binv matrix
// B6 

void BinvMultiply_B6(const double b[], double a[])
{

 a[0] =    b[0];
 a[1] =    b[9];
 a[2] = -  b[10] - 2.0*b[9] + 3.0*(b[1] - b[0]);
 a[3] =    b[9] + b[10] + 2.0*(b[0] - b[1]);
 a[4] =    b[18];
 a[5] =    b[36];
 a[6] = -  b[37] - 2.0*b[36] + 3.0*(b[19] - b[18]);
 a[7] =    b[36] + b[37] + 2.0*(b[18] - b[19]);
 a[8] = -  b[20] - 2.0*b[18] + 3.0*(b[2] - b[0]);
 a[9] = -  b[38] - 2.0*b[36] + 3.0*(b[11] - b[9]);
 
 a[10] =    b[39] + 2.0*(b[37] + b[38]) + 3.0*(b[10] + b[20] - b[12] - b[21]) + 4.0*b[36]
        +   6.0*(b[9] + b[18] - b[11] - b[19]) + 9.0*(b[0] + b[3] - b[1] - b[2]);

 a[11] = -  b[38] - b[39] + 2.0*(b[21] - b[20] - b[36] - b[37]) + 3.0*(b[11] + b[12] - b[9] - b[10])
         +  4.0*(b[19] - b[18]) + 6.0*(b[1] + b[2] - b[0] - b[3]);

 a[12] =    b[18] + b[20] + 2.0*(b[0] - b[2]);
 
 a[13] =    b[36] + b[38] + 2.0*(b[9] - b[11]);

 a[14] = -  b[37] - b[39] + 2.0*(b[12] - b[10] - b[36] - b[38]) + 3.0*(b[19] + b[21] - b[18] - b[20])
         +  4.0*(b[11] - b[9]) + 6.0*(b[1] + b[2] - b[0] - b[3]);

 a[15] =    b[36] + b[37] + b[38] + b[39] + 2.0*(b[9] + b[10] + b[18] + b[20] - b[11] - b[12] - b[19] - b[21])
         +  4.0*(b[0] + b[3] - b[1] - b[2]);

 a[16] =    b[27];

 a[17] =    b[45];

 a[18] = -  b[46] - 2.0*b[45] + 3.0*(b[28] - b[27]);

 a[19] =    b[45] + b[46] + 2.0*(b[27] - b[28]);

 a[20] =    b[54];

 a[21] =    b[7] + b[40] + b[47] + b[55] - b[15] - b[23] - b[30]
         +  5.0*(b[11] + b[13] + b[19] + b[22] + b[28] + b[29] - b[3] - b[5] - b[6] - b[36] - b[45] - b[54])
         -  8.0*b[63] + 16.0*(b[44] + b[53] + b[62]) + 25.0*(b[1] + b[2] + b[4] - b[9] - b[18] - b[27])
         -  32.0*(b[17] + b[26] + b[35]) + 64.0*b[8] - 125.0*b[0];

 a[22] =    b[16] - b[41] - b[48] + 2.0*(b[15] - b[40] - b[47]) - 4.0*b[55] + 5.0*(b[37] + b[46] - b[12] - b[14])
         +  7.0*(b[23] + b[30] - b[7]) + 8.0*b[54] + 10.0*(b[36] + b[45] - b[11] - b[13]) + 11.0*(b[6] - b[22] - b[29])
         -  16.0*b[62] + 24.0*b[63] + 25.0*b[10] + 32.0*(b[26] + b[35]) + 35.0*(b[3] + b[5] - b[19] - b[28])
         -  48.0*(b[44] + b[53]) + 50.0*b[9] + 55.0*(b[18] + b[27] - b[2] - b[4]) - 64.0*b[8] + 96.0*b[17] - 175.0*b[1]
         +  275.0*b[0];

 a[23] =    b[40] + b[41] + b[47] + b[48] - b[15] - b[16] + 4.0*(b[55] - b[54])
         +  5.0*(b[11] + b[12] + b[13] + b[14] - b[36] - b[37] - b[45] - b[46]) 
         +  6.0*(b[7] + b[22] + b[29] - b[6] - b[23] - b[30]) - 16.0*b[63] - 25.0*(b[9] + b[10])
         +  30.0*(b[2] + b[4] + b[19] + b[28] - b[3] - b[5] - b[18] - b[27]) + 32.0*(b[44] + b[53]) - 64.0*b[17]
         +  150.0*(b[1] - b[0]);

 a[24] = -  b[56] - 2.0*b[54] + 3.0*(b[29] - b[27]);

 a[25] =   b[25] - b[42] - b[57] + 2.0*(b[23] - b[40] - b[55]) - 4.0*b[47] + 5.0*(b[38] + b[56] - b[21] - b[24])
         +  7.0*(b[15] + b[30] - b[7]) + 8.0*b[45] + 10.0*(b[36] + b[54] - b[19] - b[22]) + 11.0*(b[5] - b[13] - b[28])
         -  16.0*b[53] + 24.0*b[63] + 25.0*b[20] + 32.0*(b[17] + b[35]) + 35.0*(b[3] + b[6] - b[11] - b[29])
         -  48.0*(b[44] + b[62]) + 50.0*b[18] + 55.0*(b[9] + b[27] - b[1] - b[4]) - 64.0*b[8] + 96.0*b[26] - 175.0*b[2]
         +  275.0*b[0];

 a[26] =   b[43] + 2.0*(b[41] + b[42]) + 4.0*(b[40] + b[48] + b[57]) - 5.0*b[39] - 7.0*(b[16] + b[25])
         + 8.0*(b[47] + b[55] - b[46] - b[56]) - 10.0*(b[37] + b[38]) + 11.0*(b[14] + b[24]) - 14.0*(b[15] + b[23])
         - 16.0*(b[45] + b[54]) - 20.0*b[36] + 22.0*(b[13] + b[22]) - 32.0*b[35] + 35.0*(b[12] + b[21]) - 40.0*b[30]
         + 48.0*(b[53] + b[62]) + 49.0*b[7] - 55.0*(b[10] + b[20]) + 64.0*b[8] + 68.0*(b[28] + b[29]) + 70.0*(b[11] + b[19])
         - 72.0*b[63] - 77.0*(b[5] + b[6]) - 96.0*(b[17] + b[26]) - 110.0*(b[9] + b[18]) - 112.0*b[27] + 121.0*b[4]
         + 144.0*b[44] - 245.0*b[3] + 385.0*(b[1] + b[2]) - 605.0*b[0];

 a[27] = - b[42] - b[43] - 2.0*(b[40] + b[41]) + 4.0*(b[56] - b[47] - b[48] - b[57]) + 5.0*(b[38] + b[39])
         + 6.0*(b[25] - b[24]) + 7.0*(b[15] + b[16]) + 8.0*(b[45] + b[46] + b[54] - b[55]) + 10.0*(b[36] + b[37])
         - 11.0*(b[13] + b[14]) + 12.0*(b[23] - b[22]) + 30.0*(b[20] - b[21]) - 32.0*b[53] - 35.0*(b[11] + b[12])
         + 36.0*(b[30] - b[29]) + 42.0*(b[6] - b[7]) + 48.0*b[63] + 55.0*(b[9] + b[10]) + 60.0*(b[18] + b[27] - b[19] - b[28])
         + 64.0*b[17] + 66.0*(b[5] - b[4]) - 96.0*b[44] + 210.0*(b[3] - b[2]) + 330.0*(b[0] - b[1]);

 a[28] =   b[54] + b[56] + 2.0*(b[27] - b[29]);

 a[29] =    b[40] + b[42] + b[55] + b[57] - b[23] - b[25] + 4.0*(b[47] - b[45])
         +  5.0*(b[19] + b[21] + b[22] + b[24] - b[36] - b[38] - b[54] - b[56]) 
         +  6.0*(b[7] + b[13] + b[28] - b[5] - b[15] - b[30]) - 16.0*b[63] - 25.0*(b[18] + b[20])
         +  30.0*(b[1] + b[4] + b[11] + b[29] - b[3] - b[6] - b[9] - b[27]) + 32.0*(b[44] + b[62]) - 64.0*b[26]
         +  150.0*(b[2] - b[0]);

 a[30] = - b[41] - b[43] - 2.0*(b[40] + b[42]) + 4.0*(b[46] - b[48] - b[55] - b[57]) + 5.0*(b[37] + b[39])
         + 6.0*(b[16] - b[14]) + 7.0*(b[23] + b[25]) + 8.0*(b[45] + b[54] + b[56] - b[47]) + 10.0*(b[36] + b[38])
         - 11.0*(b[22] + b[24]) + 12.0*(b[15] - b[13]) + 30.0*(b[10] - b[12]) - 32.0*b[62] - 35.0*(b[19] + b[21])
         + 36.0*(b[30] - b[28]) + 42.0*(b[5] - b[7]) + 48.0*b[63] + 55.0*(b[18] + b[20]) + 60.0*(b[9] + b[27] - b[11] - b[29])
         + 64.0*b[26] + 66.0*(b[6] - b[4]) - 96.0*b[44] + 210.0*(b[3] - b[1]) + 330.0*(b[0] - b[2]);

 a[31] =   b[40] + b[41] + b[42] + b[43] + 4.0*(b[47] + b[48] + b[55] + b[57] - b[45] - b[46] - b[54] - b[56])
         - 5.0*(b[36] + b[37] + b[38] + b[39]) + 6.0*(b[13] + b[14] + b[22] + b[24] - b[15] - b[16] - b[23] - b[25])
         + 30.0*(b[11] + b[12] + b[19] + b[21] - b[9] - b[10] - b[18] - b[20]) + 32.0*(b[28] + b[29] - b[27] - b[30] - b[63])
         + 36.0*(b[4] + b[7] - b[5] - b[6]) + 64.0*b[44] + 180.0*(b[1] + b[2] - b[0] - b[3]);

 a[32] = - b[31] - 2.0*b[27] + 3.0*(b[4] - b[0]);

 a[33] = - b[49] - 2.0*b[45] + 3.0*(b[13] - b[9]);

 a[34] =   b[50] + 2.0*(b[46] + b[49]) + 3.0*(b[10] + b[31] - b[14] - b[32]) + 4.0*b[45] + 6.0*(b[9] + b[27] - b[13] - b[28])
         + 9.0*(b[0] + b[5] - b[1] - b[4]);

 a[35] = - b[49] - b[50] + 2.0*(b[32] - b[31] - b[45] - b[46]) + 3.0*(b[13] + b[14] - b[9] - b[10]) + 4.0*(b[28] - b[27])
         + 6.0*(b[1] + b[4] - b[0] - b[5]);

 a[36] = - b[58] - 2.0*b[54] + 3.0*(b[22] - b[18]);

 a[37] =   b[34] - b[51] - b[59] + 2.0*(b[30] - b[47] - b[55]) - 4.0*b[40] + 5.0*(b[49] + b[58] - b[32] - b[33])
         +  7.0*(b[15] + b[23] - b[7]) + 8.0*b[36] + 10.0*(b[45] + b[54] - b[28] - b[29]) + 11.0*(b[3] - b[11] - b[19])
         -  16.0*b[44] + 24.0*b[63] + 25.0*b[31] + 32.0*(b[17] + b[26]) + 35.0*(b[5] + b[6] - b[13] - b[22])
         -  48.0*(b[53] + b[62]) + 50.0*b[27] + 55.0*(b[9] + b[18] - b[1] - b[2]) - 64.0*b[8] + 96.0*b[35] - 175.0*b[4]
         +  275.0*b[0];

 a[38] =   b[52] + 2.0*(b[48] + b[51]) + 4.0*(b[41] + b[47] + b[59]) - 5.0*b[50] - 7.0*(b[16] + b[34])
         + 8.0*(b[40] + b[55] - b[37] - b[58]) - 10.0*(b[46] + b[49]) + 11.0*(b[12] + b[33]) - 14.0*(b[15] + b[30])
         - 16.0*(b[36] + b[54]) - 20.0*b[45] + 22.0*(b[11] + b[29]) - 32.0*b[26] + 35.0*(b[14] + b[32]) - 40.0*b[23]
         + 48.0*(b[44] + b[62]) + 49.0*b[7] - 55.0*(b[10] + b[31]) + 64.0*b[8] + 68.0*(b[19] + b[22]) + 70.0*(b[13] + b[28])
         - 72.0*b[63] - 77.0*(b[3] + b[6]) - 96.0*(b[17] + b[35]) - 110.0*(b[9] + b[27]) - 112.0*b[18] + 121.0*b[2]
         + 144.0*b[53] - 245.0*b[5] + 385.0*(b[1] + b[4]) - 605.0*b[0];

 a[39] = - b[51] - b[52] - 2.0*(b[47] + b[48]) + 4.0*(b[58] - b[40] - b[41] - b[59]) + 5.0*(b[49] + b[50])
         + 6.0*(b[34] - b[33]) + 7.0*(b[15] + b[16]) + 8.0*(b[36] + b[37] + b[54] - b[55]) + 10.0*(b[45] + b[46])
         - 11.0*(b[11] + b[12]) + 12.0*(b[30] - b[29]) + 30.0*(b[31] - b[32]) - 32.0*b[44] - 35.0*(b[13] + b[14])
         + 36.0*(b[23] - b[22]) + 42.0*(b[6] - b[7]) + 48.0*b[63] + 55.0*(b[9] + b[10]) + 60.0*(b[18] + b[27] - b[19] - b[28])
         + 64.0*b[17] + 66.0*(b[3] - b[2]) - 96.0*b[53] + 210.0*(b[5] - b[4]) + 330.0*(b[0] - b[1]);

 a[40] =   b[60] + 2.0*(b[56] + b[58]) + 3.0*(b[20] + b[31] - b[24] - b[33]) + 4.0*b[54] + 6.0*(b[18] + b[27] - b[22] - b[29])
         + 9.0*(b[0] + b[6] - b[2] - b[4]);

 a[41] =   b[61] + 2.0*(b[57] + b[59]) + 4.0*(b[42] + b[51] + b[55]) - 5.0*b[60] - 7.0*(b[25] + b[34])
         + 8.0*(b[40] + b[47] - b[38] - b[49]) - 10.0*(b[56] + b[58]) + 11.0*(b[21] + b[32]) - 14.0*(b[23] + b[30])
         - 16.0*(b[36] + b[45]) - 20.0*b[54] + 22.0*(b[19] + b[28]) - 32.0*b[17] + 35.0*(b[24] + b[33]) - 40.0*b[15]
         + 48.0*(b[44] + b[53]) + 49.0*b[7] - 55.0*(b[20] + b[31]) + 64.0*b[8] + 68.0*(b[11] + b[13]) + 70.0*(b[22] + b[29])
         - 72.0*b[63] - 77.0*(b[3] + b[5]) - 96.0*(b[26] + b[35]) - 110.0*(b[18] + b[27]) - 112.0*b[9] + 121.0*b[1]
         + 144.0*b[62] - 245.0*b[6] + 385.0*(b[2] + b[4]) - 605.0*b[0];

 a[42] = -   4.0*(b[43] + b[52] + b[61]) + 8.0*(b[39] + b[50] + b[60] - b[41] - b[42] - b[48] - b[51] - b[57] - b[59])
         +  16.0*(b[37] + b[38] + b[46] + b[49] + b[56] + b[58] - b[40] - b[47] - b[55]) + 32.0*(b[36] + b[45] + b[54])
         +  40.0*(b[16] + b[25] + b[34]) - 64.0*b[8] - 68.0*(b[12] + b[14] + b[21] + b[24] + b[32] + b[33]) + 80.0*(b[15] + b[23] + b[30])
         +  96.0*(b[17] + b[26] + b[35]) + 112.0*(b[10] + b[20] + b[31]) - 136.0*(b[11] + b[13] + b[19] + b[22] + b[28] + b[29])
         - 144.0*(b[44] + b[53] + b[62]) + 216.0*b[63] + 224.0*(b[9] + b[18] + b[27]) - 316.0*b[7] + 512.0*(b[3] + b[5] + b[6])
         - 820.0*(b[1] + b[2] + b[4]) + 1304.0*b[0];

 a[43] =     4.0*(b[42] + b[43] + b[51] + b[52] + b[61] - b[60])
         +   8.0*(b[40] + b[41] + b[47] + b[48] + b[57] + b[59] - b[38] - b[39] - b[49] - b[50] - b[56] - b[58])
         +  16.0*(b[55] - b[36] - b[37] - b[45] - b[46] - b[54])
         +  36.0*(b[24] + b[33] - b[25] - b[34])
         -  40.0*(b[15] + b[16])
         +  60.0*(b[21] + b[32] - b[20] - b[31])
         -  64.0*b[17]
         +  68.0*(b[11] + b[12] + b[13] + b[14])
         +  72.0*(b[22] + b[29] - b[23] - b[30])
         +  96.0*(b[44] + b[53])
         - 112.0*(b[9]  + b[10])
         + 120.0*(b[19] + b[28] - b[18] - b[27])
         - 144.0* b[63]
         + 276.0*(b[7] - b[6])
         + 444.0*(b[2] + b[4] - b[3] - b[5])
         + 708.0*(b[1] - b[0]);

 a[44] = -  b[58] - b[60] + 2.0*(b[33] - b[31] - b[54] - b[56]) + 3.0*(b[22] + b[24] - b[18] - b[20]) + 4.0*(b[29] - b[27])
         +  6.0*(b[2] + b[4] - b[0] - b[6]);

 a[45] = - b[59] - b[61] - 2.0*(b[55] + b[57]) + 4.0*(b[49] - b[40] - b[42] - b[51]) + 5.0*(b[58] + b[60])
         + 6.0*(b[34] - b[32]) + 7.0*(b[23] + b[25]) + 8.0*(b[36] + b[38] + b[45] - b[47]) + 10.0*(b[54] + b[56])
         - 11.0*(b[19] + b[21]) + 12.0*(b[30] - b[28]) + 30.0*(b[31] - b[33]) - 32.0*b[44] - 35.0*(b[22] + b[24])
         + 36.0*(b[15] - b[13]) + 42.0*(b[5] - b[7]) + 48.0*b[63] + 55.0*(b[18] + b[20]) + 60.0*(b[9] + b[27] - b[11] - b[29])
         + 64.0*b[26] + 66.0*(b[3] - b[1]) - 96.0*b[62] + 210.0*(b[6] - b[4]) + 330.0*(b[0] - b[2]);

 a[46] =     4.0*(b[41] + b[43] + b[52] + b[59] + b[61] - b[50])
         +   8.0*(b[40] + b[42] + b[48] + b[51] + b[55] + b[57] - b[37] - b[39] - b[46] - b[49] - b[58] - b[60])
         +  16.0*(b[47] - b[36] - b[38] - b[45] - b[54] - b[56])
         +  36.0*(b[14] + b[32] - b[16] - b[34])
         -  40.0*(b[23] + b[25])
         +  60.0*(b[12] + b[33] - b[10] - b[31])
         -  64.0*b[26]
         +  68.0*(b[19] + b[21] + b[22] + b[24])
         +  72.0*(b[13] + b[28] - b[15] - b[30])
         +  96.0*(b[44] + b[62])
         - 112.0*(b[18] + b[20])
         + 120.0*(b[11] + b[29] - b[9]  - b[27])
         - 144.0* b[63]
         + 276.0*(b[7] - b[5])
         + 444.0*(b[1] + b[4] - b[3] - b[6])
         + 708.0*(b[2] - b[0]);

 a[47] =     4.0*(b[49] + b[50] + b[58] + b[60] - b[40] - b[41] - b[42] - b[43] - b[51] - b[52] - b[59] - b[61])
         +   8.0*(b[36] + b[37] + b[38] + b[39] + b[45] + b[46] + b[54] + b[56] - b[47] - b[48] - b[55] - b[57])
         +  32.0*(b[31] + b[34] - b[32] - b[33])
         +  36.0*(b[15] + b[16] + b[23] + b[25] - b[13] - b[14] - b[22] - b[24])
         +  60.0*(b[9]  + b[10] + b[18] + b[20] - b[11] - b[12] - b[19] - b[21])
         +  64.0*(b[27] + b[30] - b[28] - b[29] - b[44])
         +  96.0* b[63]
         + 240.0*(b[5] + b[6] - b[4] - b[7])
         + 384.0*(b[0] + b[3] - b[1] - b[2]);

 a[48] =     b[27] + b[31] + 2.0*(b[0] - b[4]);

 a[49] =     b[45] + b[49] + 2.0*(b[9] - b[13]);

 a[50] = -  b[46] - b[50] + 2.0*(b[14] - b[10] - b[45] - b[49]) + 3.0*(b[28] + b[32] - b[27] - b[31]) + 4.0*(b[13] - b[9])
         +  6.0*(b[1] + b[4] - b[0] - b[5]);

 a[51] =    b[45] + b[46] + b[49] + b[50] + 2.0*(b[9] + b[10] + b[27] + b[31] - b[13] - b[14] - b[28] - b[32])
         + 4.0*(b[0] + b[5] - b[1] - b[4]);

 a[52] =    b[54] + b[58] + 2.0*(b[18] - b[22]);

 a[53] =    b[47] + b[51] + b[55] + b[59] - b[30] - b[34] + 4.0*(b[40] - b[36])
         +  5.0*(b[28] + b[29] + b[32] + b[33] - b[45] - b[49] - b[54] - b[58])
         +  6.0*(b[7] + b[11] + b[19] - b[3] - b[15] - b[23]) - 16.0*b[63] - 25.0*(b[27] + b[31])
         +  30.0*(b[1] + b[2] + b[13] + b[22] - b[5] - b[6] - b[9] - b[18]) + 32.0*(b[53] + b[62]) - 64.0*b[35]
         +  150.0*(b[4] - b[0]);

 a[54] = - b[48] - b[52] - 2.0*(b[47] + b[51]) + 4.0*(b[37] - b[41] - b[55] - b[59]) + 5.0*(b[46] + b[50])
         + 6.0*(b[16] - b[12]) + 7.0*(b[30] + b[34]) + 8.0*(b[36] + b[54] + b[58] - b[40]) + 10.0*(b[45] + b[49])
         - 11.0*(b[29] + b[33]) + 12.0*(b[15] - b[11]) + 30.0*(b[10] - b[14]) - 32.0*b[62] - 35.0*(b[28] + b[32])
         + 36.0*(b[23] - b[19]) + 42.0*(b[3] - b[7]) + 48.0*b[63] + 55.0*(b[27] + b[31]) + 60.0*(b[9] + b[18] - b[13] - b[22])
         + 64.0*b[35] + 66.0*(b[6] - b[2]) - 96.0*b[53] + 210.0*(b[5] - b[1]) + 330.0*(b[0] - b[4]);

 a[55] =   b[47] + b[48] + b[51] + b[52] + 4.0*(b[40] + b[41] + b[55] + b[59] - b[36] - b[37] - b[54] - b[58])
         - 5.0*(b[45] + b[46] + b[49] + b[50]) + 6.0*(b[11] + b[12] + b[29] + b[33] - b[15] - b[16] - b[30] - b[34])
         + 30.0*(b[13] + b[14] + b[28] + b[32] - b[9] - b[10] - b[27] - b[31]) + 32.0*(b[19] + b[22] - b[18] - b[23] - b[63])
         + 36.0*(b[2] + b[7] - b[3] - b[6]) + 64.0*b[53] + 180.0*(b[1] + b[4] - b[0] - b[5]);

 a[56] = -  b[56] - b[60] + 2.0*(b[24] - b[20] - b[54] - b[58]) + 3.0*(b[29] + b[33] - b[27] - b[31]) + 4.0*(b[22] - b[18])
         +  6.0*(b[2] + b[4] - b[0] - b[6]);

 a[57] = - b[57] - b[61] - 2.0*(b[55] + b[59]) + 4.0*(b[38] - b[42] - b[47] - b[51]) + 5.0*(b[56] + b[60])
         + 6.0*(b[25] - b[21]) + 7.0*(b[30] + b[34]) + 8.0*(b[36] + b[45] + b[49] - b[40]) + 10.0*(b[54] + b[58])
         - 11.0*(b[28] + b[32]) + 12.0*(b[23] - b[19]) + 30.0*(b[20] - b[24]) - 32.0*b[53] - 35.0*(b[29] + b[33])
         + 36.0*(b[15] - b[11]) + 42.0*(b[3] - b[7]) + 48.0*b[63] + 55.0*(b[27] + b[31]) + 60.0*(b[9] + b[18] - b[13] - b[22])
         + 64.0*b[35] + 66.0*(b[5] - b[1]) - 96.0*b[62] + 210.0*(b[6] - b[2]) + 330.0*(b[0] - b[4]);
 
 a[58] =     4.0*(b[43] + b[48] + b[52] + b[57] + b[61] - b[39])
         +   8.0*(b[41] + b[42] + b[47] + b[51] + b[55] + b[59] - b[37] - b[38] - b[46] - b[50] - b[56] - b[60])
         +  16.0*(b[40] - b[36] - b[45] - b[49] - b[54] - b[58])
         +  36.0*(b[12] + b[21] - b[16] - b[25])
         -  40.0*(b[30] + b[34])
         +  60.0*(b[14] + b[24] - b[10] - b[20])
         -  64.0*b[35]
         +  68.0*(b[28] + b[29] + b[32] + b[33])
         +  72.0*(b[11] + b[19] - b[15] - b[23])
         +  96.0*(b[53] + b[62])
         - 112.0*(b[27] + b[31])
         + 120.0*(b[13] + b[22] - b[9]  - b[18])
         - 144.0* b[63]
         + 276.0*(b[7] - b[3])
         + 444.0*(b[1] + b[2] - b[5] - b[6])
         + 708.0*(b[4] - b[0]);
 
 a[59] =     4.0*(b[38] + b[39] + b[56] + b[60] - b[42] - b[43] - b[47] - b[48] - b[51] - b[52] - b[57] - b[61])
         +   8.0*(b[36] + b[37] + b[45] + b[46] + b[49] + b[50] + b[54] + b[58] - b[40] - b[41] - b[55] - b[59])
         +  32.0*(b[20] + b[25] - b[21] - b[24])
         +  36.0*(b[15] + b[16] + b[30] + b[34] - b[11] - b[12] - b[29] - b[33])
         +  60.0*(b[9]  + b[10] + b[27] + b[31] - b[13] - b[14] - b[28] - b[32])
         +  64.0*(b[18] + b[23] - b[19] - b[22] - b[53])
         +  96.0* b[63]
         + 240.0*(b[3] + b[6] - b[2] - b[7])
         + 384.0*(b[0] + b[5] - b[1] - b[4]);
 
 a[60] =   b[54] + b[56] + b[58] + b[60] + 2.0*(b[18] + b[20] + b[27] + b[31] - b[22] - b[24] - b[29] - b[33])
         + 4.0*(b[0] + b[6] - b[2] - b[4]);

 a[61] =   b[55] + b[57] + b[59] + b[61] + 4.0*(b[40] + b[42] + b[47] + b[51] - b[36] - b[38] - b[45] - b[49])
         - 5.0*(b[54] + b[56] + b[58] + b[60]) + 6.0*(b[19] + b[21] + b[28] + b[32] - b[23] - b[25] - b[30] - b[34])
         + 30.0*(b[22] + b[24] + b[29] + b[33] - b[18] - b[20] - b[27] - b[31]) + 32.0*(b[11] + b[13] - b[9] - b[15] - b[63])
         + 36.0*(b[1] + b[7] - b[3] - b[5]) + 64.0*b[62] + 180.0*(b[2] + b[4] - b[0] - b[6]);

 a[62] =     4.0*(b[37] + b[39] + b[46] + b[50] - b[41] - b[43] - b[48] - b[52] - b[55] - b[57] - b[59] - b[61])
         +   8.0*(b[36] + b[38] + b[45] + b[49] + b[54] + b[56] + b[58] + b[60] - b[40] - b[42] - b[47] - b[51])
         +  32.0*(b[10] + b[16] - b[12] - b[14])
         +  36.0*(b[23] + b[25] + b[30] + b[34] - b[19] - b[21] - b[28] - b[32])
         +  60.0*(b[18] + b[20] + b[27] + b[31] - b[22] - b[24] - b[29] - b[33])
         +  64.0*(b[9]  + b[15] - b[11] - b[13] - b[62])
         +  96.0* b[63]
         + 240.0*(b[3] + b[5] - b[1] - b[7])
         + 384.0*(b[0] + b[6] - b[2] - b[4]);

 a[63] =      4.0*(b[40] + b[41] + b[42] + b[43] + b[47] + b[48] + b[51] + b[52] + b[55] + b[57] + b[59] + b[61]
                 - b[36] - b[37] - b[38] - b[39] - b[45] - b[46] - b[49] - b[50] - b[54] - b[56] - b[58] - b[60])
          +  32.0*(b[11] + b[12] + b[13] + b[14] + b[19] + b[21] + b[22] + b[24] + b[28] + b[29] + b[32] + b[33]
                 - b[9]  - b[10] - b[15] - b[16] - b[18] - b[20] - b[23] - b[25] - b[27] - b[30] - b[31] - b[34])
          -  64.0* b[63]
          + 208.0*(b[1] + b[2] + b[4] + b[7] - b[0] - b[3] - b[5] - b[6]);

}

