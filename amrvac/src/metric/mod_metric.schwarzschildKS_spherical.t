!============================================================================
!##############################################################################
! SUB module METRIC - schwarzschild - KERRSCHILD
!=============================================================================
! Project : KERR -STATIC schwarzschild KERRSCHILD : Conformal Decomposition 
! Aim     :
! Ref     :
! pre-compile: setamrvac -p=gr
! par file : typemric='schwarzschildKS'
! update : 20/07/2014,  Zakaria

!============================================================================
subroutine setmetric_initialize_schwarzschildKS_spherical(pm)
  include 'amrvacdef.f'
  type(themetric), intent(inout)    :: pm
  ! .. local variables
  integer         :: i,j,k
 !----------------------------------------------------
! element of the metric

! in the schwarzschild KERRSCHILD -spherical coordinate case
 
 Loop_i : do i = 1,setgr%ndir
  pm%elemsp(i,i)%elpt%elm_on=.true. ! diagonal
  pm%Qsp(i,i)%elpt%elm_on=.true.
  if(pm%cell_center) then
   Loop_k : do k=1,setgr%ndir
      if(k==r_.and.i==r_)then
        pm%elemsp(i,i)%elpt%drv(k)%elm_on=.true.
      else
        pm%elemsp(i,i)%elpt%drv(k)%elm_on=.false.
      end if
   end do Loop_k
  end if
  if(i<ndim)then
   Loop_j1 : do j = i+1,setgr%ndir
    pm%elemsp(i,j)%elpt%elm_on= .false.
    pm%Qsp(i,j)%elpt%elm_on=.false.
    if(pm%cell_center)pm%elemsp(i,j)%elpt%drv(:)%elm_on=.false.
   end do Loop_j1
  end if
 end do Loop_i

 pm%elemsp(theta_,theta_)%elpt%elm_isone=.true.
 pm%elemsp(phi_,phi_)%elpt%elm_isone=.true.
! the lapse function
 pm%alfa%elm_on              =.true.
 pm%alfa%elm_isone           =.false.
 if(pm%cell_center)then
  if(setgr%nostatic)pm%alfa%drv(time_)%elm_on   =.false.
  pm%alfa%drv(r_)%elm_on      =.true.
  pm%alfa%drv(theta_)%elm_on  =.false.
  pm%alfa%drv(phi_)%elm_on    =.false.

  pm%alfa%drv(r_)%elm_isone    =.false.
 end if
! the shift vector

 
 pm%bt_cont(r_)%sftpt%elm_on              =.true.
 pm%bt_cont(theta_)%sftpt%elm_on          =.false.
 pm%bt_cont(phi_)%sftpt%elm_on            =.false.
 if(pm%cell_center)then
  pm%bt_cont(r_)%sftpt%drv(r_)%elm_on     =.true.
  pm%bt_cont(r_)%sftpt%drv(theta_)%elm_on =.false.
  pm%bt_cont(r_)%sftpt%drv(phi_)%elm_on   =.false.
 end if

end subroutine setmetric_initialize_schwarzschildKS_spherical
!=============================================================================
subroutine setmetric_spacetime_schwarzschildKS_spherical(ixI^L,ixO^L,x,pm)
 include 'amrvacdef.f'
 integer, intent(in)               :: ixI^L,ixO^L
 double precision, intent(in)      :: x(ixI^S,1:ndim)
 type(themetric), intent(inout)    :: pm
 ! .. local variables ..
 double precision, dimension(ixG^T):: rstor,grdelta
 !-----------------------------------------------------------------------------
 setgr%diag         = .true.
 setgr%diagspace    = .true.
 ! .. setting the local variables ..

 rstor(ixO^S)   = setgr%rs/x(ixO^S,r_)
 grdelta(ixO^S) = one+rstor(ixO^S)
 
 
 !..................................


 call setmetric_lapse_schwarzschildKS_spherical(pm%alfa,pm%cell_center)
 call setmetric_shiftvector_contr_schwarzschildKS_spherical(pm%bt_cont,pm%cell_center)
 call setmetric_space_schwarzschildKS_spherical(pm%alfa,pm%elemsp,pm%cell_center)

 contains
!===============================================================================

!-----------------------------------------------------------------------------
!============================================================================
 subroutine setmetric_lapse_schwarzschildKS_spherical(palpha,drv_needed)
  type(lapse_metric), intent(inout) :: palpha
  logical, intent(in)               :: drv_needed
!------------------------------------------------------------
  palpha%wg(ixO^S)        = one/dsqrt(grdelta(ixO^S))
  if(drv_needed)then
   if(palpha%drv(r_)%elm_needed)then
    ! set derivative of the laps fonction
    palpha%drv(r_)%dwg(ixO^S)=half*rstor(ixO^S)/x(ixO^S,r_)&
                               / (grdelta(ixO^S))**(3.0/2.0)
   end if
  end if
 end subroutine setmetric_lapse_schwarzschildKS_spherical


!============================================================================
subroutine setmetric_shiftvector_contr_schwarzschildKS_spherical(pbeta_cont,drv_needed)
 type(pshift_metric), intent(inout) :: pbeta_cont(:)
 logical, intent(in)                :: drv_needed
!-----------------------------------------------------------
 pbeta_cont(r_)%sftpt%wg(ixO^S)=rstor(ixO^S)/grdelta(ixO^S)
 if(drv_needed)then
  if(pbeta_cont(r_)%sftpt%drv(r_)%elm_needed)then
   ! set derivative of the shift vector
   pbeta_cont(r_)%sftpt%drv(r_)%dwg(ixO^S)=-rstor(ixO^S)/x(ixO^S,r_)&
               /(grdelta(ixO^S))**2.0d0
  end if
 end if
end subroutine setmetric_shiftvector_contr_schwarzschildKS_spherical
!============================================================================


 subroutine setmetric_space_schwarzschildKS_spherical(palpha,pmc,drv_needed)
  type(lapse_metric), intent(in)        :: palpha
  type(pelements_metric), intent(inout) :: pmc(:,:)
  logical, intent(in)                   :: drv_needed
! .. local variable ...

!-----------------------------------------------------------------------------
 ! fill diagonal elements
 pmc(r_,r_)%elpt%wg(ixO^S)         = one+rstor(ixO^S)

 if(drv_needed)then 
 ! fill derivative of space metric
  pmc(r_,r_)%elpt%drv(r_)%dwg(ixO^S)= -rstor(ixO^S)/x(ixO^S,r_)
 end if                                      
 end subroutine setmetric_space_schwarzschildKS_spherical
! ============================================================================
end subroutine setmetric_spacetime_schwarzschildKS_spherical


