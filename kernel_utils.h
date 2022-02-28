
#include "utilities.h"

class Kernel {
 public:
  virtual ~Kernel() = default;
  virtual void eval_bvec(double x, double y, double z, double b[]) const = 0;
  virtual double formula(double x, double y, double z) const = 0;
  virtual double Fformula(double x, double y, double z) const = 0;
};

//***********************************************************************//

class FDiffKernel : public Kernel {
 public:
  ~FDiffKernel() override = default;
  void eval_bvec(double x, double y, double z, double b[]) const override {
    // second order finite differences
    double r2 = x*x + y*y + z*z; 
    double r2h = r2 + hsq;
    double r4h = r2 + 4.0*hsq;
    
    b[0] = fEval(sqrt(r2));
    
    b[1] = f1derivOrd2(r2h, x);
    b[2] = f1derivOrd2(r2h, y);
    b[3] = f1derivOrd2(r2h, z);
  }
};

//***********************************************************************//

class CoulombKernel : public Kernel {
 public:
  ~CoulombKernel() override = default;
  void eval_bvec(double x, double y, double z, double b[]) const override {
    double rr2 = 1./(x*x + y*y + z*z);
    double D1  = sqrt(rr2);
    double D3  = D1*rr2;
    double D5  = D3*rr2;
    
    b[0] = D1;
    
    b[1] = x*D3;
    b[2] = y*D3;
    b[3] = z*D3;
    
    b[4] = 3.0*x*y*D5;
    b[5] = 3.0*x*z*D5;
    b[6] = 3.0*y*z*D5;

    b[7] = 5.0*x*b[6]*rr2;
  }
  double formula(double x, double y, double z) const override {
    return 1.0 / sqrt(x*x + y*y + z*z); //Coulomb kernel
  }
  double Fformula(double x, double y, double z) const override {
    double r2 = x*x + y*y + z*z;
    return -1.0 / (r2*sqrt(r2)); //Coulomb kernel
  }
};

class FDiffCoulombKernel : public FDiffKernel {
 public:
  ~FDiffCoulombKernel() override = default;
  double formula(double x, double y, double z) const override {
    return 1.0 / sqrt(x*x + y*y + z*z); //Coulomb kernel
  }
  double Fformula(double x, double y, double z) const override {
    double r2 = x*x + y*y + z*z;
    return -1.0 / (r2*sqrt(r2)); //Coulomb kernel
  }
};

//***********************************************************************//

class ScreenedCoulombKernel : public Kernel {
 public:
  ~ScreenedCoulombKernel() override = default;
  void eval_bvec(double x, double y, double z, double b[]) const override {
    double r2 = x*x + y*y + z*z;
    double r = sqrt(r2);
    double F = exp(-kappa*r)/r;
    double rr = 1.0/r;
    double rr2 = rr*rr;
    double k2  = kappa*kappa;
    double k1r = kappa + rr;
    double dF = F*rr*k1r;
    double dF2 = F*rr2*(k2 + 3.0*rr*k1r);
    b[0] = F;
    
    b[1] = x*dF;
    b[2] = y*dF;
    b[3] = z*dF;

    b[4] = x*y*dF2;
    b[5] = x*z*dF2;
    b[6] = y*z*dF2;

    b[7] = x*y*z*F*rr2*rr*(k2*kappa + 3.0*rr*(2*k2+5.0*rr*k1r));
  }
  double formula(double x, double y, double z) const override {
    double R = sqrt(x*x + y*y + z*z);	
    return exp(-kappa*R)/R;   //Screened Coulomb kernel
  }
  double Fformula(double x, double y, double z) const override {
    double R2 = x*x + y*y + z*z;
    double R  = sqrt(R2);             
    return -exp(-kappa*R)*(kappa*R + 1)/(R2*R); //Screened Coulomb kernel
  }
};

class FDiffScreenedCoulombKernel : public FDiffKernel {
 public:
  ~FDiffScreenedCoulombKernel() override = default;
  double formula(double x, double y, double z) const override {
    double R  = sqrt(x*x + y*y + z*z);	
    return exp(-kappa*R)/R; //Screened Coulomb kernel
  }
  double Fformula(double x, double y, double z) const override {
    double R2 = x*x + y*y + z*z;
    double R  = sqrt(R2);             
    return -exp(-kappa*R)*(kappa*R + 1)/(R2*R); //Screened Coulomb kernel
  }
};

//***********************************************************************//

class RSEwaldSumKernel : public Kernel {
 public:
  ~RSEwaldSumKernel() override = default;
  void eval_bvec(double x, double y, double z, double b[]) const override {
    double r2 = x*x + y*y + z*z;
    double r = sqrt(r2);
    double rr = 1.0/r;
    double rr2 = rr*rr;
    double F = erfc(kappa*r)*rr;
    double k2 = kappa*kappa;
    double fk3 = kappa*k2;
    double  ekr2pi = exp(-k2*r2)/rootPI; 
    double G = rr2*(2.0*kappa*ekr2pi + F);
    double H = rr2*(3.0*G + fk3*ekr2pi);

    b[0] = F;
    
    b[1] = x*G;
    b[2] = y*G;
    b[3] = z*G;

    b[4] = x*y*H;
    b[5] = x*z*H;
    b[6] = y*z*H;

    b[7] = x*y*z*rr2*(15.0*rr2*G + fk3*ekr2pi*(2.0*k2+5.0*rr2)); 
  }
  double formula(double x, double y, double z) const override {
    double R = sqrt(x*x + y*y + z*z);
    return erfc(kappa*R)/R; //Real Space Ewald Sum kernel
  }
  double Fformula(double x, double y, double z) const override {
    double R2 = x*x + y*y + z*z;
    double R  = sqrt(x*x + y*y + z*z);
    return -(2.0*kappa/rootPI*R*exp(-kappa*kappa*R2) + erfc(kappa*R))/(R2*R); //Real Space Ewald Sum kernel
  }
};

class FDiffRSEwaldSumKernel : public FDiffKernel {
 public:
  ~FDiffRSEwaldSumKernel() override = default;
  double formula(double x, double y, double z) const override {
    double R = sqrt(x*x + y*y + z*z);	
    return erfc(R)/R; //Real Space Ewald Sum kernel
  }
  double Fformula(double x, double y, double z) const override {
    double R2 = x*x + y*y + z*z;
    double R  = sqrt(x*x + y*y + z*z);
    return -(2.0*kappa/rootPI*R*exp(-kappa*kappa*R2) + erfc(kappa*R))/(R2*R); //Real Space Ewald Sum kernel
  }
};

//*******************************************************************//

void kerEval_LM(double x, double y, double z, double b[]);
