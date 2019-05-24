!============================================================================
!##############################################################################
! SUB module METRIC - kerr - KERRSCHILD
!=============================================================================
! Project : KERR -STATIC kerr KERRSCHILD : Conformal Decomposition 
! Aim     :
! Ref     :
! pre-compile: setamrvac -p=gr
! par file : typemric='kerrKS'
! update : 21/03/2015,  Zakaria

!============================================================================
subroutine setmetric_initialize_kerrKS_spherical(pm)
 ! Set the metric configuration
  include 'amrvacdef.f'
  type(themetric), intent(inout)    :: pm
  ! .. local variables
  integer         :: i,j,k
 !----------------------------------------------------
! element of the metric

! in the kerr KERRSCHILD -spherical coordinate case
 
 Loop_i : do i = 1,setgr%ndir
  pm%elemsp(i,i)%elpt%elm_on=.true. ! diagonal
  pm%Qsp(i,i)%elpt%elm_on=.true.
  cond_CC_1 : if(pm%cell_center) then
   Loop_k : do k=1,setgr%ndir
      if(k==r_.or.k==theta_)then
        pm%elemsp(i,i)%elpt%drv(k)%elm_on=.true.
      else
        pm%elemsp(i,i)%elpt%drv(k)%elm_on=.false.
      end if
   end do Loop_k
  end if cond_CC_1
 end do Loop_i

 pm%elemsp(r_,phi_)%elpt%elm_on=.true.
   if(pm%cell_center) then
   Loop_k2 : do k=1,setgr%ndir
      if(k==r_.or.k==theta_)then
        pm%elemsp(r_,phi_)%elpt%drv(k)%elm_on=.true.
      else
        pm%elemsp(r_,phi_)%elpt%drv(k)%elm_on=.false.
      end if
   end do Loop_k2
  end if
 pm%elemsp(r_,r_)%elpt%elm_isone          = .false.
 pm%elemsp(theta_,theta_)%elpt%elm_isone  = .false.
 pm%elemsp(phi_,phi_)%elpt%elm_isone      = .false.
 pm%elemsp(r_,phi_)%elpt%elm_isone        = .false.

! the lapse function
 pm%alfa%elm_on              =.true.
 pm%alfa%elm_isone           =.false.
 if(pm%cell_center)then
  if(setgr%nostatic)pm%alfa%drv(time_)%elm_on   =.false.
  pm%alfa%drv(r_)%elm_on      =.true.
  pm%alfa%drv(theta_)%elm_on  =.true.
  pm%alfa%drv(phi_)%elm_on    =.false.

  pm%alfa%drv(r_)%elm_isone     =.false.
  pm%alfa%drv(theta_)%elm_isone =.false.
 end if
! the shift vector

 
 pm%bt_cont(r_)%sftpt%elm_on              =.true.
 pm%bt_cont(theta_)%sftpt%elm_on          =.false.
 pm%bt_cont(phi_)%sftpt%elm_on            =.false.
 if(pm%cell_center)then
  pm%bt_cont(r_)%sftpt%drv(r_)%elm_on         =.true.
  pm%bt_cont(r_)%sftpt%drv(theta_)%elm_on     =.true.
  pm%bt_cont(r_)%sftpt%drv(phi_)%elm_on       =.false.
  pm%bt_cont(phi_)%sftpt%drv(r_)%elm_on     =.false.
  pm%bt_cont(phi_)%sftpt%drv(theta_)%elm_on =.false.
  pm%bt_cont(phi_)%sftpt%drv(phi_)%elm_on   =.false.
 end if

end subroutine setmetric_initialize_kerrKS_spherical
!=============================================================================
subroutine setmetric_spacetime_kerrKS_spherical(ixI^L,ixO^L,x,pm)
 use mod_mat
 include 'amrvacdef.f'
 integer, intent(in)               :: ixI^L,ixO^L
 double precision, intent(in)      :: x(ixI^S,1:ndim)
 type(themetric), intent(inout)    :: pm
 ! .. local variables ..
 double precision, dimension(ixG^T):: rstor,grdelta,grrho2,ator,sint,sin2t,cost
 !-----------------------------------------------------------------------------
 setgr%diag         = .false.
 setgr%diagspace    = .false.
 ! .. setting the local variables ..
 ator(ixO^S)    = setgr%a/x(ixO^S,r_)
 rstor(ixO^S)   = setgr%rs/x(ixO^S,r_)
 {^ZIN
 sint(ixO^S)    = dsin(x(ixO^S,theta_))
 sin2t(ixO^S)   = dsin(2.0d0*x(ixO^S,theta_))
 cost(ixO^S)    = dcos(x(ixO^S,theta_))
 grrho2(ixO^S)  = one+(ator(ixO^S)*cost(ixO^S))**2.0d0
 grdelta(ixO^S) = one+rstor(ixO^S)/grrho2(ixO^S)
 }
 {^NOZIN
 sint(ixO^S)=one;
 sin2t(ixO^S)=zero;
 cost(ixO^S) =zero;

 grrho2(ixO^S)=1.0d0
 grdelta(ixO^S) = one+rstor(ixO^S)
 }
 
 !..................................


 call setmetric_lapse_kerrKS_spherical(pm%alfa,pm%cell_center)
 call setmetric_shiftvector_contr_kerrKS_spherical(pm%bt_cont,pm%cell_center)
 call setmetric_space_kerrKS_spherical(pm%alfa,pm%elemsp,pm%cell_center)

 contains
!============================================================================

!-----------------------------------------------------------------------------
!============================================================================
 subroutine setmetric_lapse_kerrKS_spherical(palpha,drv_needed)
  type(lapse_metric), intent(inout) :: palpha
  logical, intent(in)               :: drv_needed
  double precision:: fprim(ixI^S)
!------------------------------------------------------------
  palpha%wg(ixO^S)        = one/dsqrt(grdelta(ixO^S))
  if(drv_needed)then
   if(palpha%drv(r_)%elm_needed)then
    ! set derivative of the laps fonction
    palpha%drv(r_)%dwg(ixO^S)=half*rstor(ixO^S)/&
                             (x(ixO^S,r_)*grrho2(ixO^S)**2.0D0)&
                             *(1-(ator(ixO^S)*cost(ixO^S))**2.0d0)&
                             * palpha%wg(ixO^S)**3.0d0

    palpha%drv(theta_)%dwg(ixO^S)=-half*&
    rstor(ixO^S)/x(ixO^S,r_)*ator(ixO^S)**2.0d0*sin2t(ixO^S)&
    * palpha%wg(ixO^S)**3.0d0/grrho2(ixO^S)**2.0D0
   end if
  end if
 end subroutine setmetric_lapse_kerrKS_spherical


!============================================================================
subroutine setmetric_shiftvector_contr_kerrKS_spherical(pbeta_cont,&
                                                                 drv_needed)
 type(pshift_metric), intent(inout) :: pbeta_cont(:)
 logical, intent(in)                :: drv_needed
!-----------------------------------------------------------
 pbeta_cont(r_)%sftpt%wg(ixO^S)=rstor(ixO^S)/(grrho2(ixO^S)+rstor(ixO^S))
 if(drv_needed)then
  if(pbeta_cont(r_)%sftpt%drv(r_)%elm_needed)then
   ! set derivative of the shift vector
   pbeta_cont(r_)%sftpt%drv(r_)%dwg(ixO^S)=-rstor(ixO^S)/x(ixO^S,r_)&
               *(one-(ator(ixO^S)*cost(ixO^S))**2.0D0)&
               / (grrho2(ixO^S)+rstor(ixO^S))**2.0d0
   pbeta_cont(r_)%sftpt%drv(theta_)%dwg(ixO^S)=rstor(ixO^S)/x(ixO^S,r_)*&
               ator(ixO^S)**2.0d0*sin2t(ixO^S)&
               /  (grrho2(ixO^S)+rstor(ixO^S))**2.0d0
  end if
 end if
end subroutine setmetric_shiftvector_contr_kerrKS_spherical
!============================================================================


 subroutine setmetric_space_kerrKS_spherical(palpha,pmc,drv_needed)
  type(lapse_metric), intent(in)        :: palpha
  type(pelements_metric), intent(inout) :: pmc(:,:)
  logical, intent(in)                   :: drv_needed
! .. local variable ...

!-----------------------------------------------------------------------------
 ! fill diagonal elements
 pmc(r_,r_)%elpt%wg(ixO^S)         = grdelta(ixO^S) 
 pmc(theta_,theta_)%elpt%wg(ixO^S) = grrho2(ixO^S)
 pmc(phi_,phi_)%elpt%wg(ixO^S)     = grrho2(ixO^S)&
                        +ator(ixO^S)**2.0D0*grdelta(ixO^S)*sint(ixO^S)**2.0D0
 pmc(r_,phi_)%elpt%wg(ixO^S)       = -grdelta(ixO^S)*ator(ixO^S)*sint(ixO^S);
 if(drv_needed)then 
 ! fill derivative of space metric
  pmc(r_,r_)%elpt%drv(r_)%dwg(ixO^S)      = -rstor(ixO^S)/x(ixO^S,r_)&
   / grrho2(ixO^S)**2.0d0*(1.0-(ator(ixO^S)*cost(ixO^S))**2.0d0)
  pmc(r_,r_)%elpt%drv(theta_)%dwg(ixO^S)  =  rstor(ixO^S)/x(ixO^S,r_)&
   / grrho2(ixO^S)**2.0d0*ator(ixO^S)**2.0D0 * sin2t(ixO^S)


  pmc(theta_,theta_)%elpt%drv(r_)%dwg(ixO^S)  = -2.0d0*ator(ixO^S)**2.0D0&
   /x(ixO^S,r_)*cost(ixO^S)**2.0d0           
  pmc(theta_,theta_)%elpt%drv(theta_)%dwg(ixO^S)=-&
   ator(ixO^S)**2.0D0/x(ixO^S,r_)*sin2t(ixO^S)

  pmc(phi_,phi_)%elpt%drv(r_)%dwg(ixO^S)      = -ator(ixO^S)**2.0D0&
     / x(ixO^S,r_)*(2+rstor(ixO^S)*(sint(ixO^S)/grrho2(ixO^S))**2.0D0&
       *(3.0D0+(ator(ixO^S)*cost(ixO^S))**2.0D0))

 pmc(phi_,phi_)%elpt%drv(theta_)%dwg(ixO^S) = ator(ixO^S)**2.0D0&
    /x(ixO^S,r_)*sin2t(ixO^S)&
  *rstor(ixO^S)*(1.0d0+ator(ixO^S)**2.0D0)/grrho2(ixO^S)**2.0D0

 pmc(r_,phi_)%elpt%drv(r_)%dwg(ixO^S)      = ator(ixO^S)&
    /x(ixO^S,r_)*sint(ixO^S)*&
   (one+2.0D0*rstor(ixO^S)/grrho2(ixO^S)**2.0D0)
 pmc(r_,phi_)%elpt%drv(theta_)%dwg(ixO^S) = -ator(ixO^S)&
    /x(ixO^S,r_)*cost(ixO^S)&
  *(1.0d0+rstor(ixO^S)/grrho2(ixO^S)**2.0D0*(1.0D0+ator(ixO^S)**2.0D0)) 
 end if                                      
 end subroutine setmetric_space_kerrKS_spherical
! ============================================================================
end subroutine setmetric_spacetime_kerrKS_spherical

!============================================================================

