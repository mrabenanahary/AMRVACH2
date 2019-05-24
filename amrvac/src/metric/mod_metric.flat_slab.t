!============================================================================
!##############################################################################
! SUB module METRIC - FLAT -
!=============================================================================
! Project : GR -STATIC FLAT
! coordinate : cartesian
! Aim     :
! Ref     :
! pre-compile: setamrvac -mtr=flat -gm=spherical
! par file : typemetric='flat'
! update : 1/02/2013,  Zakaria

!============================================================================
subroutine setmetric_initialize_flat_slab(pm)
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

  pm%elemsp(i,i)%elpt%elm_isone=.true.
  pm%Qsp(i,i)%elpt%elm_isone=.true.
  if(i<ndim)then
   Loop_j1 : do j = i+1,setgr%ndir
    pm%elemsp(i,j)%elpt%elm_on= .false.
    pm%Qsp(i,j)%elpt%elm_on=.false.
    if(pm%cell_center)pm%elemsp(i,j)%elpt%drv(:)%elm_on=.false.
   end do Loop_j1
  end if
 end do Loop_i
 
! the lapse function
 pm%alfa%elm_on              =.true.
 pm%alfa%elm_isone           =.true.
 if(pm%cell_center)then
  if(setgr%nostatic)pm%alfa%drv(time_)%elm_on   =.false.
  pm%alfa%drv(r_)%elm_on      =.false.
  pm%alfa%drv(theta_)%elm_on  =.false.
  pm%alfa%drv(phi_)%elm_on    =.false.
 end if
! the shift vector

 
 pm%bt_cont(r_)%sftpt%elm_on              =.false. 
 pm%bt_cont(theta_)%sftpt%elm_on          =.false.
 pm%bt_cont(phi_)%sftpt%elm_on            =.false.
end subroutine setmetric_initialize_flat_slab
!=============================================================================
subroutine setmetric_spacetime_flat_slab(ixI^L,ixO^L,x,pm)
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

 
 
 !..................................



!============================================================================
! ============================================================================
end subroutine setmetric_spacetime_flat_slab


