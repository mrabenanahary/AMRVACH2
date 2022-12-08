/*#include "headcpp.h"
#include "array.h"
#include "point.h"
#include "dim_array.h"
#include "polar.h"
#include "system_of_eqs.h"
#include "name_tools.h"
#include "matrice.h"
#include <math.h>
#include "scalar.h"*/

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
#include <iostream>
#include <fstream>
extern "C"  {
 extern int readfromC(int , int ,double *, double *, double *,
                     double *, double *, double *, double *,double *, 
                     double *,double *,double *, double *, double *,
                    double *, double *,double *, double *);
}
int readfromC(int first, int nsize, double *posr, double *post,
              double *g11, double *g22,double *g33, 
              double *alpha, double *beta, 
              double* dg11dr,double* dg11dt,
              double* dg22dr,double* dg22dt,
              double* dg33dr,
              double* dg33dt , double* dalphadr, double* dalphadt, 
              double* dbetadr, double*  dbetadt) {
	// Read the input data file
	int kk ;
	double omega ;
  
	FILE* fin = fopen ("bos_1_0.800000.dat", "r") ;
	//FILE* fin = fopen (argv[1], "r") ;
        //cout<<"well true"<<endl; 
	Space_polar space(fin) ;
	//cout << space << endl ;
	fread_be (&kk, sizeof(int), 1, fin) ;
	fread_be (&omega, sizeof(double), 1, fin) ;
	Scalar nu (space,fin) ;
	Scalar incA (space, fin) ;
	Scalar incB (space, fin) ;
	Scalar incNphi (space, fin) ;
	Scalar phi (space, fin) ;
	fclose(fin) ;
	phi.affect_parameters() ;
        phi.set_parameters()->set_m_quant() = kk ;
      //cout<<"well"<<endl;
	
	int ndom = space.get_nbr_domains() ;
	int nr = space.get_domain(0)->get_nbr_coefs()(0) ;
	int nt = space.get_domain(0)->get_nbr_coefs()(1) ;
	
	// Reconstruction mapping :
	Array<double> bounds (ndom-1) ;
	for (int d=0 ; d<ndom-1 ; d++) {
       //   cout<<d<<endl;
	  Index pos (space.get_domain(d)->get_nbr_points()) ;
	  pos.set(0) = space.get_domain(d)->get_nbr_points()(0)-1 ; // Last point in r 
	  bounds.set(d) = space.get_domain(d)->get_radius()(pos) ;
	}

	//cout << bounds << endl ;
	//cout << "Domain 0 : r= a x avec a = " << bounds(0) << endl ;
	//for (int d=1 ; d<ndom-1 ; d++) 
	 // cout << "Domain " << d << " : r = a x + b avec a = " << (bounds(d)-bounds(d-1))/2. << " et b = " << (bounds(d)+bounds(d-1))/2. << endl ;
	//cout << "Domain " << ndom-1 << " : 1/r = a (x-1) avec a = " << -0.5/bounds(ndom-2) << endl ;

	// Compute the "true" fields
	Scalar lapse (exp(nu)) ;
        lapse.std_base() ;
	
        Scalar bigA (exp(incA-nu)) ;
        bigA.std_base() ;

        Scalar bigB ((incB.div_rsint()+1)/lapse) ;
        bigB.std_base() ;
	
      	Scalar Np(incNphi.div_rsint()) ;	

        
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
        //cout<<beta[i]<<"  "<<posr[i]<<"  "<<rad.val_point(MM)<<endl;
        }
}

