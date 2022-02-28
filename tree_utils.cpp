#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "tree_utils.h"

using namespace std;

size_t node_count = 0;

int N0; // leaf size
double mtheta;
double sq_theta; // theta^2
int max_level = 0;
int MACflag;
double maxCradius = 0.0;
vector<panel> tree;
vector<size_t> leaf;

//*********************************************************************//

void read_tree_params(string &treeMethod, string &kernelName, int &N_cube, int &N0, double &theta)
{
  ifstream ifile("input_params.txt");
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

      getline(ifile,line);
      getline(ifile,line);
      stream << line << endl;
      stream >> N0;
      getline(ifile,line);

      getline(ifile,line);
      getline(ifile,line);
      stream << line << endl;
      stream >> MACflag;
      getline(ifile,line);

      getline(ifile,line);
      getline(ifile,line);
      stream << line << endl;
      stream >> theta;
      getline(ifile,line);

      getline(ifile,line);
      getline(ifile,line);
      stream << line << endl;
      stream >> fdiff;
      getline(ifile,line);

      getline(ifile,line);
      getline(ifile,line);
      stream << line << endl;
      stream >> mgrid;
      getline(ifile,line);

      getline(ifile,line);
      getline(ifile,line);
      stream << line << endl;
      stream >> hgrid;
    }
  ifile.close();
}

//************************************************************************//

// read direct sum (exact data) from a file

void read_direct_sum(const string& kernelName,
		     int N_cube,
		     double v_true[],
                     double fx_true[],
                     double fy_true[],
                     double fz_true[],
		     long &ds_cpu_time)
{  
  int idp;
  char file_name[256];
  sprintf(file_name, "exact_sum_%s_N%d", kernelName.c_str(), N_cube);
  ifstream rfile(file_name);
  string line;
  stringstream stream;
  
  for (int i = 0; i < N_cube; i++)
    {
      getline(rfile,line);
      stream << line << endl;      
      stream >> idp >> v_true[i] >> fx_true[i] >> fy_true[i] >> fz_true[i];
    }
  getline(rfile,line);
  stream << line << endl;
  stream >>  ds_cpu_time;
  
  rfile.close();
}

//*********************************************************************//

void build_tree_init(int newMAC, struct xyz &particles)
{
  panel temp_panel;
  
  // indices of particles belonging to panel
  temp_panel.members[0] = 0;
  temp_panel.members[1] = N_cube - 1;

  double min_x = particles.x[0];
  double max_x = particles.x[0];
  double min_y = particles.y[0];
  double max_y = particles.y[0];
  double min_z = particles.z[0];
  double max_z = particles.z[0];
  
  for (size_t index = 1; index < N_cube-1; index++)
    {
      min_x = min(min_x, particles.x[index]);
      max_x = max(max_x, particles.x[index]);
      min_y = min(min_y, particles.y[index]);
      max_y = max(max_y, particles.y[index]);
      min_z = min(min_z, particles.z[index]);
      max_z = max(max_z, particles.z[index]);
    }

  temp_panel.xinterval[0] = min_x; // interval defining the panel
  temp_panel.xinterval[1] = max_x;
  temp_panel.yinterval[0] = min_y;
  temp_panel.yinterval[1] = max_y;
  temp_panel.zinterval[0] = min_z;
  temp_panel.zinterval[1] = max_z;
  temp_panel.xc = 0.5 * (temp_panel.xinterval[0] + temp_panel.xinterval[1]);
  temp_panel.yc = 0.5 * (temp_panel.yinterval[0] + temp_panel.yinterval[1]);
  temp_panel.zc = 0.5 * (temp_panel.zinterval[0] + temp_panel.zinterval[1]);
  temp_panel.rxl = 1.0/(max_x-min_x);
  temp_panel.ryl = 1.0/(max_y-min_y);
  temp_panel.rzl = 1.0/(max_z-min_z);
  if (newMAC == 0)
    {
//      temp_panel.MAC = (3 * 10 * 10/ 4) / sq_theta; // MAC = r^2 / theta^2
       temp_panel.MAC = ((max_x-temp_panel.xc)*(max_x-temp_panel.xc) +
                         (max_y-temp_panel.yc)*(max_y-temp_panel.yc) +
                         (max_z-temp_panel.zc)*(max_z-temp_panel.zc)) / sq_theta;
    }
  else
    {
       temp_panel.radius = sqrt((max_x-temp_panel.xc)*(max_x-temp_panel.xc) +
                         (max_y-temp_panel.yc)*(max_y-temp_panel.yc) +
                         (max_z-temp_panel.zc)*(max_z-temp_panel.zc));
       temp_panel.MAC = (temp_panel.radius + mtheta)*(temp_panel.radius + mtheta);
    }
  tree.push_back(temp_panel);
  node_count = 1;
}

//**********************************************************************//

void Swap(size_t i, size_t j, double *lambda, struct xyz &s)
{
  if (i == j)
    return;
  
  double x = s.x[i];
  double y = s.y[i];
  double z = s.z[i];
  size_t index = s.index[i];
  size_t old_index = s.old_index[i];
  double lam = lambda[i];
  
  s.x[i] = s.x[j];
  s.y[i] = s.y[j];
  s.z[i] = s.z[j];
  s.index[i] = s.index[j];
  s.old_index[i] = s.old_index[j];
  lambda[i] = lambda[j];
  
  s.x[j] = x;
  s.y[j] = y;
  s.z[j] = z;
  s.index[j] = index;
  s.old_index[j] = old_index;
  lambda[j] = lam;
}

//**********************************************************************//

void split_tree_node(int newMAC,
		     size_t panel_index,
		     double *lambda,
		     struct xyz &particles)
{
  panel child[8];
  
  double tp_x0 = tree[panel_index].xinterval[0];
  double tp_x1 = tree[panel_index].xinterval[1];
  double tp_y0 = tree[panel_index].yinterval[0];
  double tp_y1 = tree[panel_index].yinterval[1];
  double tp_z0 = tree[panel_index].zinterval[0];
  double tp_z1 = tree[panel_index].zinterval[1];

  double midpointx = (tp_x0 + tp_x1) / 2.0;
  double midpointy = (tp_y0 + tp_y1) / 2.0;
  double midpointz = (tp_z0 + tp_z1) / 2.0;

  double xc0 = (tp_x0 + midpointx) / 2.0;
  double xc1 = (tp_x1 + midpointx) / 2.0;
  double yc0 = (tp_y0 + midpointy) / 2.0;
  double yc1 = (tp_y1 + midpointy) / 2.0;
  double zc0 = (tp_z0 + midpointz) / 2.0;
  double zc1 = (tp_z1 + midpointz) / 2.0;

  child[0].xinterval[0] = tp_x0;
  child[0].xinterval[1] = midpointx;
  child[0].yinterval[0] = tp_y0;
  child[0].yinterval[1] = midpointy;
  child[0].zinterval[0] = tp_z0;
  child[0].zinterval[1] = midpointz;
  child[0].xc = xc0;
  child[0].yc = yc0;
  child[0].zc = zc0;
  child[0].rxl = 1.0/(midpointx-tp_x0);
  child[0].ryl = 1.0/(midpointy-tp_y0);
  child[0].rzl = 1.0/(midpointz-tp_z0);
  if (newMAC == 0)
    {
      child[0].MAC = ((midpointx - xc0) * (midpointx - xc0) + (midpointy - yc0) * (midpointy - yc0) + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;
    }
  else
    {
      child[0].radius = sqrt((midpointx - xc0) * (midpointx - xc0) + (midpointy - yc0) * (midpointy - yc0) + (midpointz - zc0) * (midpointz - zc0));
      child[0].MAC = (child[0].radius + mtheta)*(child[0].radius + mtheta);
    }

  child[1].xinterval[0] = midpointx;
  child[1].xinterval[1] = tp_x1;
  child[1].yinterval[0] = tp_y0;
  child[1].yinterval[1] = midpointy;
  child[1].zinterval[0] = tp_z0;
  child[1].zinterval[1] = midpointz;
  child[1].xc = xc1;
  child[1].yc = yc0;
  child[1].zc = zc0;
  child[1].rxl = 1.0/(tp_x1-midpointx);
  child[1].ryl = 1.0/(midpointy-tp_y0);
  child[1].rzl = 1.0/(midpointz-tp_z0);
  if (newMAC == 0)
    {
      child[1].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1) + (midpointy - yc0) * (midpointy - yc0) + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;
    }
  else
    {
      child[1].radius = sqrt((tp_x1 - xc1) * (tp_x1 - xc1) + (midpointy - yc0) * (midpointy - yc0) + (midpointz - zc0) * (midpointz - zc0));
      child[1].MAC = (child[1].radius + mtheta)*(child[1].radius + mtheta);
    }
  
  child[2].xinterval[0] = tp_x0;
  child[2].xinterval[1] = midpointx;
  child[2].yinterval[0] = midpointy;
  child[2].yinterval[1] = tp_y1;
  child[2].zinterval[0] = tp_z0;
  child[2].zinterval[1] = midpointz;
  child[2].xc = xc0;
  child[2].yc = yc1;
  child[2].zc = zc0;
  child[2].rxl = 1.0/(midpointx-tp_x0);
  child[2].ryl = 1.0/(tp_y1-midpointy);
  child[2].rzl = 1.0/(midpointz-tp_z0);
  if (newMAC == 0)
    {
      child[2].MAC = ((midpointx - xc0) * (midpointx - xc0) + (tp_y1 - yc1) * (tp_y1 - yc1) + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;
    }
  else
    {
      child[2].radius = ((midpointx - xc0) * (midpointx - xc0) + (tp_y1 - yc1) * (tp_y1 - yc1) + (midpointz - zc0) * (midpointz - zc0));
      child[2].MAC = (child[2].radius + mtheta)*(child[2].radius + mtheta);
    }

  child[3].xinterval[0] = midpointx;
  child[3].xinterval[1] = tp_x1;
  child[3].yinterval[0] = midpointy;
  child[3].yinterval[1] = tp_y1;
  child[3].zinterval[0] = tp_z0;
  child[3].zinterval[1] = midpointz;
  child[3].xc = xc1;
  child[3].yc = yc1;
  child[3].zc = zc0;
  child[3].rxl = 1.0/(tp_x1-midpointx);
  child[3].ryl = 1.0/(tp_y1-midpointy);
  child[3].rzl = 1.0/(midpointz-tp_z0);
  if (newMAC == 0)
    {
      child[3].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1) + (tp_y1 - yc1) * (tp_y1 - yc1) + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;
    }
  else
    {
      child[3].radius = sqrt((tp_x1 - xc1) * (tp_x1 - xc1) + (tp_y1 - yc1) * (tp_y1 - yc1) + (midpointz - zc0) * (midpointz - zc0));
      child[3].MAC = (child[3].radius + mtheta)*(child[3].radius + mtheta);
    }
  
  child[4].xinterval[0] = tp_x0;
  child[4].xinterval[1] = midpointx;
  child[4].yinterval[0] = tp_y0;
  child[4].yinterval[1] = midpointy;
  child[4].zinterval[0] = midpointz;
  child[4].zinterval[1] = tp_z1;
  child[4].xc = xc0;
  child[4].yc = yc0;
  child[4].zc = zc1;
  child[4].rxl = 1.0/(midpointx-tp_x0);
  child[4].ryl = 1.0/(midpointy-tp_y0);
  child[4].rzl = 1.0/(tp_z1-midpointz);
  if (newMAC == 0)
    {
      child[4].MAC = ((midpointx - xc0) * (midpointx - xc0) + (midpointy - yc0) * (midpointy - yc0) + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;
    }
  else
    {
      child[4].radius = sqrt((midpointx - xc0) * (midpointx - xc0) + (midpointy - yc0) * (midpointy - yc0) + (tp_z1 - zc1) * (tp_z1 - zc1));
      child[4].MAC = (child[4].radius + mtheta)*(child[4].radius + mtheta);
    }
  
  child[5].xinterval[0] = midpointx;
  child[5].xinterval[1] = tp_x1;
  child[5].yinterval[0] = tp_y0;
  child[5].yinterval[1] = midpointy;
  child[5].zinterval[0] = midpointz;
  child[5].zinterval[1] = tp_z1;
  child[5].xc = xc1;
  child[5].yc = yc0;
  child[5].zc = zc1;
  child[5].rxl = 1.0/(tp_x1-midpointx);
  child[5].ryl = 1.0/(midpointy-tp_y0);
  child[5].rzl = 1.0/(tp_z1-midpointz);
  if (newMAC == 0)
    {
      child[5].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1) + (midpointy - yc0) * (midpointy - yc0) + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;
    }
  else
    {
      child[5].radius = sqrt((tp_x1 - xc1) * (tp_x1 - xc1) + (midpointy - yc0) * (midpointy - yc0) + (tp_z1 - zc1) * (tp_z1 - zc1));
      child[5].MAC = (child[5].radius + mtheta)*(child[5].radius + mtheta);
    }
  
  child[6].xinterval[0] = tp_x0;
  child[6].xinterval[1] = midpointx;
  child[6].yinterval[0] = midpointy;
  child[6].yinterval[1] = tp_y1;
  child[6].zinterval[0] = midpointz;
  child[6].zinterval[1] = tp_z1;
  child[6].xc = xc0;
  child[6].yc = yc1;
  child[6].zc = zc1;
  child[6].rxl = 1.0/(midpointx-tp_x0);
  child[6].ryl = 1.0/(tp_y1-midpointy);
  child[6].rzl = 1.0/(tp_z1-midpointz);
  if (newMAC == 0)
    {
      child[6].MAC = ((midpointx - xc0) * (midpointx - xc0) + (tp_y1 - yc1) * (tp_y1 - yc1) + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;
    }
  else
    {
      child[6].radius = sqrt((midpointx - xc0) * (midpointx - xc0) + (tp_y1 - yc1) * (tp_y1 - yc1) + (tp_z1 - zc1) * (tp_z1 - zc1));
      child[6].MAC = (child[6].radius + mtheta)*(child[6].radius + mtheta);
    }
  
  child[7].xinterval[0] = midpointx;
  child[7].xinterval[1] = tp_x1;
  child[7].yinterval[0] = midpointy;
  child[7].yinterval[1] = tp_y1;
  child[7].zinterval[0] = midpointz;
  child[7].zinterval[1] = tp_z1;
  child[7].xc = xc1;
  child[7].yc = yc1;
  child[7].zc = zc1;
  child[7].rxl = 1.0/(tp_x1-midpointx);
  child[7].ryl = 1.0/(tp_y1-midpointy);
  child[7].rzl = 1.0/(tp_z1-midpointz);
  if (newMAC == 0)
    {
      child[7].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1) + (tp_y1 - yc1) * (tp_y1 - yc1) + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;
    }
  else
    {
      child[7].radius = sqrt((tp_x1 - xc1) * (tp_x1 - xc1) + (tp_y1 - yc1) * (tp_y1 - yc1) + (tp_z1 - zc1) * (tp_z1 - zc1));
      child[7].MAC = (child[7].radius + mtheta)*(child[7].radius + mtheta);
    }
  
  vector<size_t> v[8];
  size_t start = tree[panel_index].members[0];
  size_t end = tree[panel_index].members[1];
  size_t* addr_table = new size_t[end - start + 1];

  size_t index;
  for (index = start; index <= end; index++)
    {
      particles.index[index] = index;
      addr_table[index - start] = index;

      if (particles.x[index] <= midpointx && particles.y[index] <= midpointy &&
	  particles.z[index] <= midpointz)
	v[0].push_back(index);
      else if (particles.x[index] > midpointx && particles.y[index] <= midpointy &&
	       particles.z[index] <= midpointz )
	v[1].push_back(index);
      else if (particles.x[index] <= midpointx && particles.y[index] > midpointy &&
	       particles.z[index]<= midpointz)
	v[2].push_back(index);
      else if (particles.x[index] > midpointx && particles.y[index] > midpointy &&
	       particles.z[index] <= midpointz)
	v[3].push_back(index);
      else if(particles.x[index] <= midpointx && particles.y[index] <= midpointy &&
	      particles.z[index] > midpointz )
	v[4].push_back(index);
      else if (particles.x[index] > midpointx && particles.y[index] <= midpointy &&
	       particles.z[index] > midpointz)
	v[5].push_back(index);
      else if (particles.x[index] <= midpointx && particles.y[index] > midpointy &&
	       particles.z[index] > midpointz)
	v[6].push_back(index);
      else if (particles.x[index] > midpointx && particles.y[index] > midpointy &&
	       particles.z[index] > midpointz)
	v[7].push_back(index);
    }

  size_t seq = start;
  for (size_t j = 0; j < 8; j++)
    {
      size_t size = v[j].size();

      if (size >= 1)
	{
	  for (size_t k = 0; k < size; k++)
	    {
	      if (k == 0)
		child[j].members[0] = seq;
	      if (k == size - 1)
		child[j].members[1] = seq;

	      index = v[j][k];
	      // This uses an address table
	      size_t pos = addr_table[index - start];
	      size_t out = particles.index[seq];
	      Swap(pos, seq, lambda, particles);
	      addr_table[index - start] = seq;
	      addr_table[out - start] = pos;

	      seq++;
	    }

	  node_count++;
	  tree[panel_index].children.push_back(node_count - 1);
	  tree.push_back(child[j]);
	  v[j].clear();
	}
    }

  delete[] addr_table;
}

//***********************************************************************//

void build_tree_3D_Recursive(int newMAC,
			     size_t panel_index,
			     double *lambda,
			     struct xyz &particles,
			     int level)
{
  if (level > max_level)
    max_level = level;
	
  size_t n = tree[panel_index].members[1] - tree[panel_index].members[0] + 1;

  if (n >= (size_t)N0)
    {
      split_tree_node(newMAC, panel_index, lambda, particles);

      for (size_t i = 0; i < tree[panel_index].children.size(); i++)
	{
	  size_t panel_index_new = tree[panel_index].children[i];
	  build_tree_3D_Recursive(newMAC,
				  panel_index_new,
				  lambda,
				  particles,
				  level + 1);
	}
    }
  else
    leaf.push_back(panel_index);
}

//************************************************************************//

// direct sum used in particle-cluster treecode

vec_4d Call_Ds_PC(int limit_1, int limit_2, 
	       double p_x, double p_y, double p_z,
	       struct xyz &particles,
	       double *lambda, const Kernel& kernel)
{    
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
        
    vec_4d velocity;
    velocity.val[0]=0.0; velocity.val[1]=0.0;
    velocity.val[2]=0.0; velocity.val[3]=0.0;

    for (size_t j = limit_1; j <= limit_2; j++)
    {
        x = p_x - particles.x[j];
        y = p_y - particles.y[j];
        z = p_z - particles.z[j];

        double wj = lambda[j];
        double kerF = wj * kernel.Fformula(x, y, z);

	velocity.val[0] += wj * kernel.formula(x, y, z);
	velocity.val[1] -= x * kerF;
        velocity.val[2] -= y * kerF;
        velocity.val[3] -= z * kerF;

    }
    return velocity;
}

//************************************************************************//

// direct sum used in cluster-particle treecode

void Call_Ds_CP(int limit_1, int limit_2, 
               double p_x, double p_y, double p_z,
               struct xyz &particles,
               const Kernel& kernel)

{    
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
      
   if (sID >= limit_1 && sID <= limit_2)
   { 
    for (size_t jj = limit_1; jj < sID; jj++)
    {
        x = particles.x[jj] - p_x;
        y = particles.y[jj] - p_y;
        z = particles.z[jj] - p_z;
        
        int old_j = particles.old_index[jj];
 
        double kerF    = sweight * kernel.Fformula(x, y, z);

        cpvelo[old_j] += sweight * kernel.formula(x, y, z);         
        cpdvx[old_j]  -= x * kerF;
        cpdvy[old_j]  -= y * kerF;
        cpdvz[old_j]  -= z * kerF;
        
    }

    for (size_t jj = sID+1; jj <= limit_2; jj++)
    {   
        x = particles.x[jj] - p_x;
        y = particles.y[jj] - p_y;
        z = particles.z[jj] - p_z;
        
        int old_j = particles.old_index[jj];
        
        double kerF    = sweight * kernel.Fformula(x, y, z);

        cpvelo[old_j] += sweight * kernel.formula(x, y, z);
        cpdvx[old_j]  -= x * kerF;
        cpdvy[old_j]  -= y * kerF;
        cpdvz[old_j]  -= z * kerF;
       
    }
   
  }
  else {

    for (size_t jj = limit_1; jj <= limit_2; jj++)
    {
        x = particles.x[jj] - p_x;
        y = particles.y[jj] - p_y;
        z = particles.z[jj] - p_z;
    
        int old_j = particles.old_index[jj];

        double kerF    = sweight * kernel.Fformula(x, y, z);

        cpvelo[old_j] += sweight * kernel.formula(x, y, z);
        cpdvx[old_j]  -= x * kerF;
        cpdvy[old_j]  -= y * kerF;
        cpdvz[old_j]  -= z * kerF;

    }

  }
}

//********************************************************************//

// save treecode computed data to a file

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
                        const double dvz[])
{
    char ofile_name[256];
    sprintf(ofile_name,
	    "tricubic_sum_%s_%s_N%d",
	    kernelName.c_str(),
	    type.c_str(),
	    N_cube);
    ofstream output_file(ofile_name);

    output_file << ds_cpu_time << "   " << treecode_cpu_time << endl;
    output_file << "============================== " << endl;
    for (size_t i = 0; i < N_cube; i++)
      {
        output_file << i << "   " << setprecision(16)
		    << v_true[i] << "   " << velo[i] << endl;
        output_file << setprecision(16) << fx_true[i] <<" "<< fy_true[i] <<" "<< fz_true[i] <<endl;
        output_file << setprecision(16) << dvx[i] <<" "<< dvy[i] <<" "<< dvz[i] <<endl;
        output_file << " " << endl;     
      }    
    output_file.close();
}

//********************************************************************//

// compute error

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
                   double &fL1)
{
  e_n = 0.0;
  e_d = 0.0;
  e_d_ex = 0.0;
  
  double fd2 = 0.0;
  double fe2 = 0.0;

  double en1 = 0.0;
  double ed1 = 0.0;
  double fd1 = 0.0;
  double fe1 = 0.0;

  for (size_t i = 0; i < N_cube; i++)
    {
      e_n += (v_true[i] - velo[i]) * (v_true[i] - velo[i]) ;
      e_d += velo[i] * velo[i];
      e_d_ex += v_true[i] * v_true[i];

      en1 += abs(v_true[i] - velo[i]);
      ed1 += abs(v_true[i]);

      double fdx = fx_true[i] - dvx[i];
      double fdy = fy_true[i] - dvy[i];
      double fdz = fz_true[i] - dvz[i];

      fd2 += fdx*fdx + fdy*fdy + fdz*fdz;
      fe2 += fx_true[i]*fx_true[i] + fy_true[i]*fy_true[i] + fz_true[i]*fz_true[i];

      fd1 += abs(fdx)+abs(fdy)+abs(fdz);
      fe1 += abs(fx_true[i])+abs(fy_true[i])+abs(fz_true[i]);
    }
  e_n = sqrt(e_n);
  e_d = sqrt(e_d);
  e_d_ex = sqrt(e_d_ex);
  
  E = e_n/e_d_ex;

  fd=sqrt(fd2); fe=sqrt(fe2);
  FE = fd/fe;

  eL1 = en1/ed1;
  fL1 = fd1/fe1;
}
