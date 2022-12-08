#include "headcpp.h"
#include "array.h"
#include "point.h"
#include "dim_array.h"
#include "polar.h"
#include "system_of_eqs.h"
#include "name_tools.h"
#include "matrice.h"
#include <math.h>
#include "scalar.h"
extern "C"  {
 extern int readfromC(int , int ,double *, double *, double *,
                     double *, double *, double *, double *, double *,
                     double *,double *,double *, double *, double *,
                    double *, double *, double *, double *);
}
int readfromC(int first, int nsize, double *posr, double *post,
              double *g11, double *g22,double *g33,
              double *alpha, double *beta,
              double* dg11dr,double* dg11dt,
              double* dg22dr,double* dg22dt,
              double* dg33dr,
              double* dg33dt , double* dalphadr, double* dalphadt,
              double* dbetadr, double*  dbetadt){
	// Read the input data file
	  // Constants :
      double M = 1. ;
      double kerrparam = 0.5;
      
      // Go on...
      double aa = kerrparam*M ;
      double hh = sqrt(M*M - aa*aa) ; // Eq. 101

      // The new config :
      int dim = 2 ;
      // Number of points
      int type_coloc = CHEB_TYPE ;
      int nr = 13 ;
      int np = 9 ;
      Dim_array res (dim) ;
      res.set(0) = nr ; res.set(1) = np ;

      Point center (2) ;	
      for (int i=1 ; i<=dim ; i++)
	center.set(i) = 0 ;

      // Domains
      int ndom = 8 ; 
      Array<double> bounds(ndom-1) ;
      bounds.set(0) = hh/2. ; // Radius of the horizon Eq. 10
      for (int d=1 ; d<ndom-1 ; d++)
	  bounds.set(d) = bounds(d-1)*2  ;
      Space_polar space(type_coloc, center, res, bounds) ;
  
      // Auxiliary quantities
      Scalar one (space) ;
      one = 1 ;
      one.std_base() ;
      Scalar cos (space) ;
      for (int d=0 ; d<ndom ; d++) {
	cos.set_domain(d) = one(d).mult_cos_theta() ;
      }
      
      Scalar rad(space) ;
      rad.set_domain(0) = 1 ;
      for (int d=1 ; d<ndom ; d++) {
	rad.set_domain(d) = space.get_domain(d)->get_radius() ;
      }
      
      Scalar bigR (space) ;
      bigR = rad + (M*M-aa*aa)/4./rad + M ; // Eq. (111)
      
      Scalar sigma (space) ;
      sigma = bigR*bigR + aa*aa*cos*cos ; // Eq. (93)

      // 3+1 quantities 
      Scalar Asquare (space) ;
      Asquare = 1.0 + 2.0*M/rad + 
          (3*M*M+aa*aa *(2*cos*cos-1))/2./rad/rad + 
          hh*hh*M/2/rad/rad/rad + hh*hh*hh*hh/16/rad/rad/rad/rad ; // Eq. 121
      Asquare.set_val_inf(1) ;
      Scalar bigA (space) ;
      bigA = sqrt(Asquare) ;
      bigA.std_base() ;
      bigA.set_domain(0) = 0 ;
    
      Scalar Nsquare (space) ;
      Nsquare = 1 - 2*M*bigR / sigma + (4*aa*aa*M*M*bigR*bigR*(1-cos*cos))/(sigma*sigma*(bigR*bigR+aa*aa) + 2*aa*aa*sigma*M*bigR*(1-cos*cos)) ; // Eq. 124
      Nsquare.set_val_inf(1) ;
      Scalar lapse (space) ;
      lapse = sqrt(Nsquare) ;
      // Regularity on the horizon is enforced
      Index pos (space.get_domain(1)->get_nbr_points()) ;
      do {
	if(pos(0)==0)
	  lapse.set_domain(1).set(pos) = 0 ;
      }
      while(pos.inc()) ;
      
      lapse.std_base() ;
      lapse.set_domain(0) = 0 ;
      
     
      Scalar Bsquare(space) ;
      Bsquare = (bigR*bigR+aa*aa+2*aa*aa*M*bigR*(1-cos*cos)/sigma)/rad/rad ; // Eq. 125
      Bsquare.set_val_inf(1) ;
      Scalar bigB (space) ;
      bigB = sqrt(Bsquare) ;
      bigB.std_base() ;
      bigB.set_domain(0) = 0 ;
        
      Scalar Np (space) ;
      Np = (2*aa*M*bigR)/(sigma*(bigR*bigR+aa*aa)+2*aa*aa*M*bigR*(1-cos*cos)) ; // Eq. 126
      Np.set_val_inf(0) ;
      Np.std_base() ;
      Np.set_domain(0) = 0 ;

	//cout << "Domain 0 : r= a x avec a = " << bounds(0) << endl ;
	//for (int d=1 ; d<ndom-1 ; d++) 
	 // cout << "Domain " << d << " : r = a x + b avec a = " << (bounds(d)-bounds(d-1))/2. << " et b = " << (bounds(d)+bounds(d-1))/2. << endl ;
	//cout << "Domain " << ndom-1 << " : 1/r = a (x-1) avec a = " << -0.5/bounds(ndom-2) << endl ;
	
        Scalar lapse_der_r(lapse.der_r());
        Scalar lapse_der_t(lapse.der_var(2));
        Scalar beta_der_r(Np.der_r());
        Scalar beta_der_t(Np.der_var(2));
        Scalar bigA_der_r(bigA.der_r());
        Scalar bigA_der_t(bigA.der_var(2));
        Scalar bigB_der_r(bigB.der_r());
        Scalar bigB_der_t(bigB.der_var(2));

	// Reconstruction mapping :
        //cout<<"lolo"<<nsize<<endl;
        for(int i=0;i<nsize;i++){
	double xx = posr[i];
	double zz = post[i];
	Point MM(2) ;
	MM.set(1) = xx ; MM.set(2) = zz;
	// Reconstruction mapping :
        g11[i]      = pow(bigA.val_point(MM),2.0);
        g22[i]      = pow(bigA.val_point(MM),2.0);
        g33[i]      = pow(bigB.val_point(MM),2.0);
        alpha[i]    = lapse.val_point(MM);
        beta[i]     = -Np.val_point(MM);
        dalphadr[i] = lapse_der_r.val_point(MM);
        dalphadt[i] = lapse_der_t.val_point(MM);
        dbetadr[i]  = -beta_der_r.val_point(MM);
        dbetadt[i]  = -beta_der_t.val_point(MM);
        dg11dr[i]   = 2.0*bigA.val_point(MM)*bigA_der_r.val_point(MM);
        dg11dt[i]   = 2.0*bigA.val_point(MM)*bigA_der_t.val_point(MM);
        dg22dr[i]   = 2.0*bigA.val_point(MM)*bigA_der_r.val_point(MM);
        dg22dt[i]   = 2.0*bigA.val_point(MM)*bigA_der_t.val_point(MM);
        dg33dr[i]   = 2.0*bigB.val_point(MM)*bigB_der_r.val_point(MM);
        dg33dt[i]   = 2.0*bigB.val_point(MM)*bigB_der_t.val_point(MM);
        }
}

