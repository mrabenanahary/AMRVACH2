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
      double keRRparam = 0.5;0.5 ;
      
      // Go on...
      double aa = keRRparam*M ;
      double hh = pow(M,2.0)+sqrt(M*M - aa*aa) ; // Eq. 101

      // The new config :
      int dim = 2 ;
      // Number of points
      int type_coloc = CHEB_TYPE ;
      int nr = 129 ;
      int np = 33 ;
      Dim_array res (dim) ;
      res.set(0) = nr ; res.set(1) = np ;

      Point center (2) ;	
      for (int i=1 ; i<=dim ; i++)
	center.set(i) = 0 ;

      // Domains
      int ndom = 8 ; 
      Array<double> bounds(ndom-1) ;
      bounds.set(0) = hh *(1.1); // Radius of the horizon Eq. 10
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
      Scalar sin (space) ;
      for (int d=0 ; d<ndom ; d++) {
        sin.set_domain(d) = one(d).mult_sin_theta() ;
      }
      
      Scalar rad(space) ;
      rad.set_domain(0) = 1 ;
      for (int d=1 ; d<ndom ; d++) {
	rad.set_domain(d) = space.get_domain(d)->get_radius() ;
      }
      rad.std_base() ;      

      Scalar bigR2 (space) ;
      bigR2 = 1.0+ pow(aa/rad*cos,2) ; // Eq. (3)

      Scalar delta (space) ; 
      delta=1.0-2.0*M/rad+pow(aa/rad,2);

      Scalar sigma (space) ;
      sigma = pow(1.0+ pow(aa/rad,2),2) - pow(aa/rad*sin,2)*delta ; // Eq. (4)

      Scalar gRR (space) ;
      gRR = bigR2/delta;
      gRR.set_domain(0) = 1 ;
 
      Scalar gtt (space) ;
      gtt = bigR2;
      gtt.set_domain(0) = 1 ;

      Scalar gff (space) ;
      gff = sigma/bigR2;
      gff.set_domain(0) = 1 ;
      gff.coef();
      Scalar dgrr (space) ;
      dgrr = -((2.0*M/rad-2.0*pow(aa/rad,2))+pow(aa/rad,2)
        *(2.0-2.0*M/rad)*cos*cos)/(rad*delta*delta);
      dgrr.set_domain(0) = 0 ;
      // Check if one can compute the coef
      dgrr.coef() ; 

      Scalar dgff (space) ;
      dgff = -(pow(aa/rad,2))/rad*(2.0+3.0*2.0*M/rad);
      dgff.set_domain(0) = 0 ;
      // Check if one can compute the coef
             dgff.coef() ;
      Scalar Nf (space) ;
      Nf =-  2.0*M*aa/(pow(rad,3)*sigma);//*sin;
      Nf.std_base() ;
      Nf.set_domain(0) = 0 ;
      Nf.coef() ;
      Scalar lapse (space) ;
      lapse = sqrt(delta*bigR2/sigma);
      lapse.std_base() ;
      lapse.set_domain(0) = 0 ;
      
	
      Scalar dlapse (space) ;
      /*dlapse =0.5/lapse*2.0*M/(rad*pow(sigma,2))*
               (pow(pow((aa/rad),2)+1.0,2)-2.0*
                pow(aa,2)*2.0*M/rad);*/ 
       dlapse =0.5/sqrt(delta*bigR2/sigma)*
               2.0*M/(pow(rad,2)*pow(sigma,2))*(pow(pow((aa/rad),2)+1.0,2)-2.0*
                pow(aa/rad,2)*2.0*M/rad);
        dlapse.set_domain(0) = 0 ;

      Scalar dbeta (space) ;
       dbeta =Nf/(rad*sigma)*(2.0-pow(aa/rad,2)*2.0*M/rad+pow(aa/rad,2)*(2.0*M/rad-2.0*pow(aa/rad,2))*cos);
        dbeta.set_domain(0) = 0 ;
        Scalar lapse_der_r(lapse.der_r());
        Scalar lapse_der_t(lapse.der_var(2));
        Scalar beta_der_r(Nf.der_r());
        Scalar beta_der_t(Nf.der_var(2));
        Scalar gRR_der_r(gRR.der_r());
        Scalar gRR_der_t(gRR.der_var(2));
        Scalar gtt_der_r(gtt.der_r());
        Scalar gtt_der_t(gtt.der_var(2));
        Scalar gff_der_r(gff.der_r());
        Scalar gff_der_t(gff.der_var(2));

	// Reconstruction mapping :
        //cout<<"lolo"<<nsize<<endl;
        for(int i=0;i<nsize;i++){
	double xx = posr[i];
	double zz = post[i];
	Point MM(2) ;
	MM.set(1) = xx ; MM.set(2) = zz;
	// Reconstruction mapping :
        g11[i]      = gRR.val_point(MM);
        g22[i]      = gtt.val_point(MM);
        g33[i]      = gff.val_point(MM);
        alpha[i]    = lapse.val_point(MM);
        beta[i]     = Nf.val_point(MM);
        dalphadr[i] = lapse_der_r.val_point(MM);
        dalphadt[i] = lapse_der_t.val_point(MM);
        dbetadr[i]  = beta_der_r.val_point(MM);
        dbetadt[i]  = beta_der_t.val_point(MM);
        dg11dr[i]   = gRR_der_r.val_point(MM);
        dg11dt[i]   = gRR_der_t.val_point(MM);
        dg22dr[i]   = gtt_der_r.val_point(MM);
        dg22dt[i]   = gtt_der_t.val_point(MM);
        dg33dr[i]   = gff_der_r.val_point(MM);
        dg33dt[i]   = gff_der_t.val_point(MM);
       /*cout<<"all test "<<dg11dr[i]-dgrr.val_point(MM)
        <<endl;*/
        }
/*for(int i=0;i<24-1;i++){
        double xx = posr[i];
        double zz = post[i];
        Point MM(2) ;
        MM.set(1) = xx ; MM.set(2) = zz;
        cout<<" beta "<<beta[i]<<" "<<sigma.val_point(MM)<<" "<<sin.val_point(MM)<<" "<<rad.val_point(MM)<<" "<<posr[i]<< "  "<<post[i]<<endl;
//	cout<<"best test "<<g11[i]<<"gifR "<<bigR2.val_point(MM)<<"delta"<<delta.val_point(MM)<<" check"<<bigR2.val_point(MM)/delta.val_point(MM)<<" r "<<rad.val_point(MM)<<" "<<posr[i]<< "  "<<post[i]<< endl;

}*/
//  abort();
}

