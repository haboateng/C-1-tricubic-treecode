// Tricubic B2 version with or without finite differences
// Moments multiplied by B^(-1)^T

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "tricubic_utils.h"

using namespace std;

//*****************************************************************************//

void Compute_CP2(size_t panel_index,
                           struct xyz &particles)
{
 if (tree[panel_index].moment_flag == 1)
 {

   double tm[Pflat];
   double dm[Pflat];

// Multiply coefficient by B inverse

   BinvMultiply_LM(tree[panel_index].moments, tm);



//============================================================
/*
   for (int i=0; i<64; i++)
     {
       tm[i] = 0.0;
       for (int j=0; j<64; j++)
         {
           tm[i] += Binv[i][j]*tree[panel_index].moments[j];
         }
     }
*/
//============================================================

  double xmin = tree[panel_index].xinterval[0];
  double ymin = tree[panel_index].yinterval[0];
  double zmin = tree[panel_index].zinterval[0];

  double rdx = tree[panel_index].rxl; 
  double rdy = tree[panel_index].ryl; 
  double rdz = tree[panel_index].rzl;

  double tp0 = tree[panel_index].members[0];
  double tp1 = tree[panel_index].members[1];

  for (size_t tp_j = tp0; tp_j <= tp1; tp_j++)
    {
      size_t old_j = particles.old_index[tp_j];

      double x = (particles.x[tp_j] - xmin) * rdx;
      double y = (particles.y[tp_j] - ymin) * rdy;
      double z = (particles.z[tp_j] - zmin) * rdz;

      double peng = 0.0;
      double tvx  = 0.0; double tvy = 0.0; double tvz = 0.0;
      double si   = 1.0;
      for (int i = 0; i < P + 1; i++)
       {
          double sij = si;
          for (int j = 0; j < P + 1; j++)
            {
              int i4j = i + 4*j;

              double s = sij;
              for (int k = 0; k < P + 1; k++)
                {
                  int ii = i4j + 16*k;
                  peng  += tm[ii] * s;
                  dm[ii] = s;

                  s *= z;
                }
                sij *= y;
            }
            si *= x;
        }
          cpvelo[old_j] += peng;

      for (int i = 1; i < P + 1; i++)
        { 

          double di = static_cast<double>(i);

          int fourI = 4*i; int sixtI = 16*i;
          int fII = 4*(i-1); int sII = 16*(i-1);

          for (int j = 0; j < P + 1; j++)
            { 

              int fourJ = 4*j;
              int i4j   = i + fourJ;
              int j4i   = j + fourI;

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

                  tvx   += tm[ii] * di * dm[ii-1];
                  tvy   += tm[jj] * di * dm[njj];
                  tvz   += tm[kk] * di * dm[nkk];

                }
            }
        }

         cpdvx[old_j]  -= rdx*tvx;
         cpdvy[old_j]  -= rdy*tvy;
         cpdvz[old_j]  -= rdz*tvz; 
    }
 }

// Loop over children
   size_t length = tree[panel_index].children.size();
   for (size_t i = 0; i < length; i++)
     {
       size_t index = tree[panel_index].children[i];
       Compute_CP2(index, particles);
     }
}

//*****************************************************************************//

void Compute_CP1(struct xyz &particles,
                 double p_x, double p_y, double p_z,
                 size_t panel_index,
                 const Kernel& kernel)
{
    size_t limit_1 = tree[panel_index].members[0];
    size_t limit_2 = tree[panel_index].members[1];   
        
    double tpx = p_x - tree[panel_index].xc; 
    double tpy = p_y - tree[panel_index].yc; 
    double tpz = p_z - tree[panel_index].zc;
        
    double R_sq = tpx * tpx + tpy * tpy + tpz * tpz;
    
    if (tree[panel_index].MAC < R_sq)
      {
        tree[panel_index].moment_flag = 1;
        Compute_Coeff_LM(p_x, p_y, p_z, panel_index, kernel);

        maxCradius = (tree[panel_index].radius > maxCradius) ? tree[panel_index].radius : maxCradius;
      }
    
    else
      {
        if (limit_2 - limit_1 < N0) //  otherwise, if cluster is a leaf, use direct sum
          {
            Call_Ds_CP(limit_1, limit_2, p_x,  p_y, p_z, particles,kernel);
          }
        else // othervise, if cluster is not a leaf, look at children
          {
            size_t length = tree[panel_index].children.size();
            for (size_t i = 0; i < length; i++)
              {
                size_t index = tree[panel_index].children[i];
                Compute_CP1(particles, p_x,  p_y, p_z, index, kernel);
              }            
          }
      }
}

//*****************************************************************************//

vec_4d Compute_Velocity(double *lambda,
			struct xyz &particles,
			double p_x, double p_y, double p_z,
			size_t panel_index,
			const Kernel& kernel)
{
    size_t limit_1 = tree[panel_index].members[0];
    size_t limit_2 = tree[panel_index].members[1];   
	
    double tpx = p_x - tree[panel_index].xc; 
    double tpy = p_y - tree[panel_index].yc; 
    double tpz = p_z - tree[panel_index].zc;
	
    double R_sq = tpx * tpx + tpy * tpy + tpz * tpz;
    
    vec_4d velocity;
    velocity.val[0] = 0.0; velocity.val[1] = 0.0;
    velocity.val[2] = 0.0; velocity.val[3] = 0.0;
	
    if (tree[panel_index].MAC < R_sq)
      {
        vec_4d tree_result = Call_Tricubic_LM(p_x, p_y, p_z,
					      panel_index,
					      kernel);
        velocity.val[0] += tree_result.val[0];
        velocity.val[1] += tree_result.val[1]*tree[panel_index].rxl;
        velocity.val[2] += tree_result.val[2]*tree[panel_index].ryl;
        velocity.val[3] += tree_result.val[3]*tree[panel_index].rzl;

        maxCradius = (tree[panel_index].radius > maxCradius) ? tree[panel_index].radius : maxCradius;
      }
    
    else
      {
	if (limit_2 - limit_1 < N0) //if cluster is a leaf, use direct sum
	  {
            vec_4d DS_result = Call_Ds_PC(limit_1, limit_2,
				       p_x,  p_y, p_z,
				       particles,
				       lambda,
				       kernel);
	    
            velocity.val[0] += DS_result.val[0];
            velocity.val[1] += DS_result.val[1];
            velocity.val[2] += DS_result.val[2];
            velocity.val[3] += DS_result.val[3];
            
	  }
	else //if cluster is not a leaf, look at children
	  {
	    velocity.val[0] = 0.0; velocity.val[1] = 0.0;
            velocity.val[2] = 0.0; velocity.val[3] = 0.0;
 
	    size_t length = tree[panel_index].children.size();
	    for (size_t i = 0; i < length; i++)
	      {
		size_t index = tree[panel_index].children[i];
                vec_4d temp_result = Compute_Velocity(lambda,
						      particles,
						      p_x,  p_y, p_z,
						      index,
						      kernel);
                velocity.val[0] += temp_result.val[0];
                velocity.val[1] += temp_result.val[1];
                velocity.val[2] += temp_result.val[2];
                velocity.val[3] += temp_result.val[3];

	      }            
	  }
      }
    return velocity;
}

//*****************************************************************************//
void lookup_table(double dr, const Kernel& kernel)
{
 fval[0]=0.0;
 
 for (int i = 1; i < mgrid; i++)
  {
    fval[i] = kernel.formula(dr * i, 0, 0);
  }
}
//*****************************************************************************//

int main()
{
    tree.reserve(5000);
    leaf.reserve(5000);

    double theta;
    string treeMethod;
    string kernelName;
    read_tree_params(treeMethod, kernelName, N_cube, N0, theta);

    cout<<" "<<endl;
    cout<<"Using "<<treeMethod<<" treecode"<<endl;
    cout<<" "<<endl;

    Kernel* kernel;
    if (fdiff==1) {
      if (kernelName == "Coulomb") {
	kernel = new FDiffCoulombKernel();
      } else if (kernelName == "ScreenedCoulomb") {
	kernel = new FDiffScreenedCoulombKernel();
      } else if (kernelName == "RSEwaldSum") {
	kernel = new FDiffRSEwaldSumKernel();
      } else {
	cerr << "Invalid kernel name: " << kernelName << endl;
	exit(-1);
      }
    }
    else {
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
    }

    mtheta = theta;
    sq_theta = theta*theta; // theta^2

    struct xyz particles(N_cube);
    double *lambda = new double[N_cube];

    cout << "===== Tricubic LM version =====" << endl;
    if (fdiff==1) cout << "===== with finite differences =====" << endl;
    cout << "Kernel = " << kernelName << endl;
    cout << "P = " << P << endl;
    cout << "N = " << N_cube << endl;
    cout << "theta = " << theta << endl;
    cout << "N0 = " << N0 << endl;

    // read particle coordinates and strengths from files
    read_particle_data(N_cube, particles, lambda);

    if (fdiff==1) {
      twoh = 2.0*hgrid; rtwoh = 1.0/twoh; hsq = hgrid*hgrid; 
      rfhsq = 0.25/hsq; r8hcb = rtwoh*rfhsq;
   
      // grid size
      
      double xlmax = particles.x[0];
      double xlmin = particles.x[0];
      double ylmax = particles.y[0];
      double ylmin = particles.y[0];
      double zlmax = particles.z[0];
      double zlmin = particles.z[0];
      
      for (int i = 1; i < N_cube; i++)
	{
	  double xi, yi, zi;
	  xi=particles.x[i]; yi=particles.y[i]; zi=particles.z[i];
	  
	  if (xi > xlmax) xlmax = xi;
	  if (xi < xlmin) xlmin = xi;
	  if (yi > ylmax) ylmax = yi;
	  if (yi < ylmin) ylmin = yi;
	  if (zi > zlmax) zlmax = zi;
	  if (zi < zlmin) zlmin = zi;
	}
      
      double maxlen = max(xlmax-xlmin,ylmax-ylmin);
      maxlen = sqrt(3.0)*max(maxlen,zlmax-zlmin);
      double dr = maxlen/static_cast<double>(mgrid);
      rdr = 1.0/dr;
      
      // Call look-up table to compute function values at grid points
      fval = new double[mgrid]; //Allocate mgrid doubles and save ptr in fval
      lookup_table(dr, *kernel);
    }
    
 /*
    if (treeMethod == "Cluster-Particle"){

    // Read the inverse matrix from a file
       Binv = new double*[64];
       for (int i = 0; i < 64; i++)
         Binv[i] = new double[64];

       ifstream BinvFile("BinvMatrix_LM.txt");
       for (int i=0; i<64; i++)
         {
           for (int j=0; j<64; j++)
             {
               double value;
               BinvFile >> value;
               Binv[i][j] = value;
             }
         }
       BinvFile.close();
    }
 */

    cout << "Starting treecode" << endl;
	
    //***************** Set up tree *******************************
    long Start_total, Start_btree;
    long End_total, End_btree;
    
    Start_total = getTickCount(); // Get currenct CPU time
    Start_btree = getTickCount();

    int newMAC = MACflag;
    build_tree_init(newMAC, particles);
    build_tree_3D_Recursive(newMAC, 0, lambda, particles, 0);

    if(treeMethod=="Particle-Cluster"){    

    //***************** Compute moment for each panel **************
       size_t size = tree.size();
    
       for (size_t k = 1; k < size; k++) // skip root
         {
	   double tm[Pflat], tmx[Pflat], tmy[Pflat], tmz[Pflat];

	   Panel_Moment_Tricubic(k, lambda, particles, tm, tmx, tmy, tmz);
           BinvMultiply_Transpose_LM(tm, tree[k].moments);
           BinvMultiply_Transpose_LM(tmx, tree[k].dxmoments);
           BinvMultiply_Transpose_LM(tmy, tree[k].dymoments);
           BinvMultiply_Transpose_LM(tmz, tree[k].dzmoments);
         }
    
    }
    End_btree = getTickCount();
    
    cout << "build tree time (incl. moments) = "
	 << End_btree - Start_btree << endl;

    //***************** Run Treecode *****************

       double *velo = new double[N_cube];
       double *dvx  = new double[N_cube];
       double *dvy  = new double[N_cube];
       double *dvz  = new double[N_cube];

    cpvelo = new double[N_cube];
    cpdvx  = new double[N_cube];
    cpdvy  = new double[N_cube];
    cpdvz  = new double[N_cube];

    if(treeMethod == "Particle-Cluster"){ 
     
      //********Particle-Cluster Treecode**************   
/*
       double *velo = new double[N_cube];
       double *dvx  = new double[N_cube];
       double *dvy  = new double[N_cube];
       double *dvz  = new double[N_cube];
*/
	
       for (int i = 0; i < N_cube; i++)
         {
	   double p_x = particles.x[i];
	   double p_y = particles.y[i];
	   double p_z = particles.z[i];
	
	   int old_i = particles.old_index[i];
	   double temp_lambda = lambda[i];
	   lambda[i] = 0.0;
	
	   double temp_x = particles.x[i];
	   particles.x[i] += 1000.0;
	
	   vec_4d vdv = Compute_Velocity(lambda,
	 			         particles,
				         p_x, p_y, p_z,
				         0,
	  			         *kernel);

           velo[old_i] = vdv.val[0];
           dvx[old_i]  = vdv.val[1];
           dvy[old_i]  = vdv.val[2];
           dvz[old_i]  = vdv.val[3];
 	
	   lambda[i] = temp_lambda;
   	   particles.x[i] = temp_x;          
	
         }
    }else if(treeMethod == "Cluster-Particle"){

      //********Cluster-Particle Treecode**************
/*
    cpvelo = new double[N_cube];
    cpdvx  = new double[N_cube];
    cpdvy  = new double[N_cube];
    cpdvz  = new double[N_cube];
*/
    
    for (int i = 0; i < N_cube; i++)
      {
        cpvelo[i] = 0.0;
      }

    double bfac = 0.037037037037037037037037;
    for (int i = 0; i < N_cube; i++)
      {
        double p_x = particles.x[i];
        double p_y = particles.y[i];
        double p_z = particles.z[i];

        int old_i = particles.old_index[i];
        sweight   = lambda[i];
        sfweight  = bfac * sweight;
        sID       = i;

        Compute_CP1(particles, p_x, p_y, p_z, 0, *kernel);

      }

      Compute_CP2(0, particles);

     }else{
      // No tree method specified
      cout << "Please specify a tree method" <<endl;
      exit (EXIT_FAILURE);
    }

    End_total = getTickCount(); // Time for all treecode computing
    long treecode_cpu_time;
    treecode_cpu_time = End_total - Start_total;
    
    cout << "treecode_cpu_time = " << treecode_cpu_time << endl;
    
    //***************** End Run Treecode *****************
    
    // ********* read exact data from a file *********

    double *v_true = new double[N_cube];
    double *fx_true = new double[N_cube];
    double *fy_true = new double[N_cube];
    double *fz_true = new double[N_cube];

    long ds_cpu_time;
    read_direct_sum(kernelName,
		    N_cube,
		    v_true,
                    fx_true,fy_true,fz_true,
		    ds_cpu_time);
    cout << "ds time = " << ds_cpu_time << endl;

    // ********* write treecode data to a file *********
    string method;
    method = "LM";

   if(treeMethod == "Particle-Cluster"){    
       write_treecode_sum(kernelName,
		       method,
		       N_cube,
		       ds_cpu_time,
		       treecode_cpu_time,
		       v_true,
                       fx_true,fy_true,fz_true,
		       velo,dvx,dvy,dvz);
    }

   if(treeMethod == "Cluster-Particle"){
       write_treecode_sum(kernelName,
                       method,
                       N_cube,
                       ds_cpu_time,
                       treecode_cpu_time,
                       v_true,
                       fx_true,fy_true,fz_true,
                       cpvelo,cpdvx,cpdvy,cpdvz);
    }

    // compute error
    double e_n, e_d, e_d_ex, E, fd,fe, FE, eL1, fL1;

    if(treeMethod == "Particle-Cluster"){
       compute_error(v_true, fx_true, fy_true, fz_true, 
                  velo, dvx, dvy, dvz, e_n, e_d, e_d_ex, E, fd, fe, FE, eL1, fL1);
    }

    if(treeMethod == "Cluster-Particle"){
       compute_error(v_true, fx_true, fy_true, fz_true,
        cpvelo, cpdvx, cpdvy, cpdvz, e_n, e_d, e_d_ex, E, fd, fe, FE, eL1, fL1);
    }

// Sum up forces

    double *lambda_old = new double[N_cube];

  // ********* read particle weights from a file *********
  
  FILE * fp;
  char lambda_Str_data_file[64] = {0};
  sprintf(lambda_Str_data_file, "./lambda_%d.txt", N_cube);
  fp = fopen(lambda_Str_data_file, "r");

  double x1;
  int count = -1;  
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
      
      lambda_old[count] = x1;
    }


    double dfxs, dfys, dfzs, tfxs, tfys, tfzs;
    dfxs=0.0; dfys=0.0; dfzs=0.0;
    tfxs=0.0; tfys=0.0; tfzs=0.0;

    if(treeMethod == "Particle-Cluster"){

       for (int i = 0; i < N_cube; i++)
           {

            double lm = lambda_old[i];
            dfxs += lm*fx_true[i]; dfys += lm*fy_true[i]; dfzs += lm*fz_true[i];
            tfxs += lm*dvx[i];     tfys += lm*dvy[i];     tfzs += lm*dvz[i];
           }
    }

    if(treeMethod == "Cluster-Particle"){

       for (int i = 0; i < N_cube; i++)
           {
            double lm = lambda_old[i];
            dfxs += lm*fx_true[i]; dfys += lm*fy_true[i]; dfzs += lm*fz_true[i];
            tfxs += lm*cpdvx[i];   tfys += lm*cpdvy[i];   tfzs += lm*cpdvz[i];
           }
    }

    cout << "L2 velocity Error = " << E << endl;
    cout << "Numerator of L2 error = " << e_n << endl;
    cout << "Denominator of L2 error (uses exact value) = " << e_d_ex << endl;
    cout << " " << endl;

    cout << "L2 force error = " << FE << endl;   
    cout << "Numerator of L2 force error = " << fd << endl;
    cout << "Denominator of L2 force error  = " << fe << endl;
    cout << " " << endl;
    cout << "Sum of direct sum forces "<< endl;
    cout << dfxs <<" "<< dfys <<" "<< dfzs << endl;
    cout << " " << endl;
    cout << "Sum of tree forces "<< endl;
    cout << tfxs <<" "<< tfys <<" "<< tfzs << endl;

    cout << " " << endl;
// Output to file
   
    ofstream output_file;
    output_file.open("output.txt", ios::out | ios::ate | ios::app);
    output_file << N_cube <<"  "<< N0 << "  " << theta <<"  "<< setprecision(16) 
    << E <<"  "<< FE << "  "<< eL1 <<"  "<< fL1 <<"  "<<ds_cpu_time <<"  "<<treecode_cpu_time << endl;

    output_file.close();
   //**************************************************
    
    delete [] lambda;
    delete [] v_true;
    delete [] fx_true;
    delete [] fy_true;
    delete [] fz_true;
   
    delete [] velo;
    delete [] dvx;
    delete [] dvy;
    delete [] dvz;

    delete [] cpvelo;
    delete [] cpdvx;
    delete [] cpdvy;
    delete [] cpdvz;

   // if (treeMethod == "Cluster-Particle") delete *Binv;
    if (fdiff==1) delete [] fval;  

    delete kernel;
    
    return 0;
}
