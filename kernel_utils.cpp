
#include "kernel_utils.h"

//************************************************************************//

// Lekien&Marsden: Coulomb only
void kerEval_LM(double x, double y, double z, double b[])
{
  
  double rr2 = 1./(x*x + y*y + z*z);
  double D1 = sqrt(rr2);
  double D3 = D1*rr2; 
  double D5 = D3*rr2; 

  b[0] = D1;

  b[1] = x*D3;
  b[2] = y*D3;
  b[3] = z*D3;
  
  b[4] = 3.0*x*y*D5;
  b[5] = 3.0*x*z*D5;
  b[6] = 3.0*y*z*D5;
  
  b[7] = 5.0*x*b[6]*rr2;

}
