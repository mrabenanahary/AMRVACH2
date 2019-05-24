!============================================================================
!##############################################################################
! SUB module METRIC - KERR -
!=============================================================================
! Project : KERR -STATIC METRIC : Conformal Decomposition 
! Aim     :
! Ref     :
! pre-compile: setamrvac -p=gr
! par file : typemric='kerr'
! update : 20/12/2012,  Zakaria

!============================================================================
subroutine setmetric_initialize_kerrbl_spherical(pm)
  include 'amrvacdef.f'
  type(themetric), intent(inout)    :: pm
  ! .. local variables
  integer         :: i,j,k
 !----------------------------------------------------
! element of the metric

! in the kerr-spherical coordinate case
 Loop_i : do i = 1,setgr%ndir
  pm%elemsp(i,i)%elpt%elm_on=.true. ! diagonal
  pm%Qsp(i,i)%elpt%elm_on=.true.
  if(pm%cell_center) then
   Loop_k : do k=1,setgr%ndir
      if(k==r_{^ZIN .or.k==theta_})then
        pm%elemsp(i,i)%elpt%drv(k)%elm_on=.true.
      end if
   end do Loop_k
  end if
 end do Loop_i
! the lapse function
 pm%alfa%elm_on              =.true.
 if(pm%cell_center)then
  pm%alfa%drv(r_)%elm_on      =.true.
  {^ZIN 
  pm%alfa%drv(theta_)%elm_on  =.true.}
 end if
! the shift vector
 pm%bt_cont(phi_)%sftpt%elm_on            =.true.
 if(pm%cell_center)then
  pm%bt_cont(phi_)%sftpt%drv(r_)%elm_on     =.true.
  {^ZIN 
  pm%bt_cont(phi_)%sftpt%drv(theta_)%elm_on =.true.}
 end if
end subroutine setmetric_initialize_kerrbl_spherical
!=============================================================================
subroutine setmetric_spacetime_kerrbl_spherical(ixI^L,ixO^L,x,pm)
 use mod_mat
 include 'amrvacdef.f'

 integer, intent(in)               :: ixI^L,ixO^L
 double precision, intent(in)      :: x(ixI^S,1:ndim)
 type(themetric), intent(inout)    :: pm
 ! .. local variables ..
 double precision, dimension(ixG^T):: ator,rstor,grrho2,grdelta,grsigma,&
                                      sint,sin2t,cost
 integer                           :: ixd^L
 !-----------------------------------------------------------------------------
 setgr%diag         = .false.
 setgr%diagspace    = .true.
 ! .. setting the local variables ..
 ator(ixO^S)    = setgr%a/x(ixO^S,r_)
 rstor(ixO^S)   = setgr%rs/x(ixO^S,r_)
 grdelta(ixO^S) = one-rstor(ixO^S)+ator(ixO^S)**2.0d0

{^ZIN
 sint(ixO^S)=dsin(x(ixO^S,theta_))
 sin2t(ixO^S)=dsin(2.0d0*x(ixO^S,theta_))
 cost(ixO^S)=dcos(x(ixO^S,theta_))
 grrho2(ixO^S)  = one+(ator(ixO^S)*cost(ixO^S))**2.0d0
 grsigma(ixO^S) = (one+ator(ixO^S)**2.0d0)**2.0d0&
                  -ator(ixO^S)**2.0d0*grdelta(ixO^S)*sint(ixO^S)**2.0d0
 }
 {^NOZIN
 sint(ixO^S)=one;
 sin2t(ixO^S)=zero;
 cost(ixO^S) =zero;
 grrho2(ixO^S)  = one
 grsigma(ixO^S) = (one+ator(ixO^S)**2.0d0)**2.0d0&
                  -ator(ixO^S)**2.0d0*grdelta(ixO^S)
 }
 
 !..................................

 ixd^L=ixI^L^LSUB2;

 call setmetric_lapse_kerrbl_spherical(pm%alfa,pm%cell_center)
 call setmetric_shiftvector_contr_kerrbl_spherical(pm%bt_cont,pm%cell_center)
 call setmetric_space_kerrbl_spherical(pm%alfa,pm%elemsp,pm%cell_center)
 contains
!===============================================================================
!-----------------------------------------------------------------------------
!============================================================================
 subroutine setmetric_shiftvector_contr_kerrbl_spherical(pbeta_cont,drv_needed)
 type(pshift_metric), intent(inout) :: pbeta_cont(:)
 logical, intent(in)                :: drv_needed

 double precision  :: dw(ixI^S) 
!-----------------------------------------------------------
 pbeta_cont(phi_)%sftpt%wg(ixO^S)=-ator(ixO^S)*rstor(ixO^S)*sint(ixO^S)/grsigma(ixO^S);
 if(drv_needed)then
   ! set derivative of the shift vector
   pbeta_cont(phi_)%sftpt%drv(r_)%dwg(ixO^S)=pbeta_cont(phi_)%sftpt%wg(ixO^S)&
                     /(grsigma(ixO^S)*x(ixO^S,r_))&
                     *(-2.0D0+ator(ixO^S)**2.0d0*rstor(ixO^S)&
                     -ator(ixO^S)**2.0d0*(rstor(ixO^S)&
                      -2.0D0*ator(ixO^S)**2.0d0)*cost(ixO^S)**2.0D0 )
{^ZIN
   pbeta_cont(phi_)%sftpt%drv(theta_)%dwg(ixO^S)=pbeta_cont(phi_)%sftpt%wg(ixO^S)&
 *ator(ixO^S)**2.0D0*grdelta(ixO^S)*sin2t(ixO^S)/(x(ixO^S,r_)*grsigma(ixO^S))&
 -ator(ixO^S)*rstor(ixO^S)*cost(ixO^S)/(x(ixO^S,r_)*grsigma(ixO^S))
}
{^NOZIN
!  pbeta_cont(phi_)%sftpt%drv(theta_)%dwg(ixO^S)=zero;
}
 end if
end subroutine setmetric_shiftvector_contr_kerrbl_spherical
!============================================================================
 subroutine setmetric_lapse_kerrbl_spherical(palpha,drv_needed)
  type(lapse_metric), intent(inout) :: palpha
  logical, intent(in)               :: drv_needed
!------------------------------------------------------------
 ! alpha=Sqrt(rho^2 *delta/Sima^2)
  palpha%wg(ixO^S)        = dsqrt(max(grrho2(ixO^S)*grdelta(ixO^S)&
                                  /grsigma(ixO^S),smalldouble))
  if(drv_needed)then
    ! set derivative of the laps fonction
    palpha%drv(r_)%dwg(ixO^S)=half/(palpha%wg(ixO^S))&
        *(rstor(ixO^S)/(x(ixO^S,r_)*grsigma(ixO^S)**2.0D0))*&
        ((1.0+ator(ixO^S)**2.0D0)**2.0D0-2.0d0*ator(ixO^S)**2.0D0*rstor(ixO^S)&
        -ator(ixO^S)**2.0D0*cost(ixO^S)**2.0D0*&
        ((1.0+ator(ixO^S)**2.0D0)**2.0D0-2.0D0*rstor(ixO^S)))
  {^ZIN
    palpha%drv(theta_)%dwg(ixO^S) = -half*rstor(ixO^S)&
            *(1.0d0+ator(ixO^S)**2.0d0)*ator(ixO^S)**2.0d0*grdelta(ixO^S)&
            *sin2t(ixO^S)/(x(ixO^S,r_)*palpha%wg(ixO^S)*grsigma(ixO^S)**2.0D0)
   }
  {^NOZIN
!  palpha%drv(theta_)%dwg(ixO^S) =0.0D0;
  }
  end if
 end subroutine setmetric_lapse_kerrbl_spherical


!============================================================================

 subroutine setmetric_space_kerrbl_spherical(palpha,pmc,drv_needed)
  type(lapse_metric), intent(in)        :: palpha
  type(pelements_metric), intent(inout) :: pmc(:,:)
  logical, intent(in)                   :: drv_needed
! .. local variable ...

!-----------------------------------------------------------------------------

 ! fill diagonal elements
 pmc(r_,r_)%elpt%wg(ixO^S)         = max(grrho2(ixO^S)/grdelta(ixO^S),smalldouble)

 pmc(theta_,theta_)%elpt%wg(ixO^S) = max(grrho2(ixO^S),smalldouble)
 pmc(phi_,phi_)%elpt%wg(ixO^S)     = max(grsigma(ixO^S)/grrho2(ixO^S),smalldouble)


 if(drv_needed)then
 ! fill derivative of space metric
 ! fill diaogonal first
 pmc(r_,r_)%elpt%drv(r_)%dwg(ixO^S)= -((rstor(ixO^S)-2.0*ator(ixO^S)**2.0d0)&
       +ator(ixO^S)**2.0d0*(2.0d0-rstor(ixO^S))*(cost(ixO^S))**2.0d0)&
       / (x(ixO^S,r_)*grdelta(ixO^S)**2.0d0)
 {^IFTHETA

 pmc(theta_,theta_)%elpt%drv(r_)%dwg(ixO^S)= -2.0d0*(ator(ixO^S)&
            *cost(ixO^S))**2.0d0/x(ixO^S,r_)
 }
 {^IFPHI
 pmc(phi_,phi_)%elpt%drv(r_)%dwg(ixO^S)= (-2.0D0*ator(ixO^S)**2.d0&
 +rstor(ixO^S)-(rstor(ixO^S)*((1.0D0+3.0D0*ator(ixO^S)**2.d0)&
   -ator(ixO^S)**2.d0*(1.0D0-ator(ixO^S)**2.d0)*cost(ixO^S)**2.0D0))&
   / grrho2(ixO^S)**2.0D0)/x(ixO^S,r_)
 }
 ! derive to theta

{^ZIN

 pmc(r_,r_)%elpt%drv(theta_)%dwg(ixO^S)=-ator(ixO^S)**2.0d0*sin2t(ixO^S)&
               / (x(ixO^S,r_)*grdelta(ixO^S))

 pmc(theta_,theta_)%elpt%drv(theta_)%dwg(ixO^S)=-ator(ixO^S)**2.0d0*sin2t(ixO^S)&
               /x(ixO^S,r_)

 pmc(phi_,phi_)%elpt%drv(theta_)%dwg(ixO^S)= rstor(ixO^S)*ator(ixO^S)**2.0d0&
                       * sin2t(ixO^S)*(one+ator(ixO^S)**2.0d0)&
                       / (x(ixO^S,r_)*grrho2(ixO^S)**2.0d0)

}
{^NOZIN
!pmc(r_,r_)%elpt%drv(theta_)%dwg(ixO^S)=zero;
!pmc(theta_,theta_)%elpt%drv(theta_)%dwg(ixO^S)=zero;
!pmc(phi_,phi_)%elpt%drv(theta_)%dwg(ixO^S)=zero;
}
 end if
 end subroutine setmetric_space_kerrbl_spherical
! ============================================================================
end subroutine setmetric_spacetime_kerrbl_spherical


