!===========================================================================
!##############################################################################
! module METRIC - GENERAL RELATIVISTY -STATIC METRIC
!=============================================================================
! Project : GENERAL RELATIVISTY -STATIC METRIC : Conformal Decomposition 
! Aim     :
! Ref     :
! pre-compile: setamrvac -p=gr
! update : 05/12/2012,  Zakaria
! based on  geometry.t 
! configuration
! typemetric
! typeaxial
! parameters :
! variable :
! main variables
! convention for the metric : (-,+,+,+)
!============================================================================
Module mod_metric
{^IFGR 
!============================================================================
! global variable in the module
 use mod_interface_metric
 use interface_kadath
 use mod_forest
 implicit none
!-------------------------------------

integer, target, save    :: LC3_EPS(1:3,1:3,1:3),Eps3(1:3,1:3)
integer,save             :: i1_gr(1:3,1:3),i2_gr(1:3,1:3)
!------------------------------------- 

type logical_needed
 logical    :: el
 logical    :: drvel
end type logical_needed
type setspacetime
 logical                 :: normalized
 logical                 :: diag,diagspace
 logical                 :: check_on,check_isone
 logical                 :: nostatic
 double precision        :: rs,a
 integer                 :: ndir
 type(logical_needed)    :: elem_needed,elem_inv_needed
 type(logical_needed)    :: elemsp_needed,elemsp_inv_needed,elemsp_updwn_needed
 type(logical_needed)    :: Qsp_needed,Qsp_inv_needed
 type(logical_needed)    :: bt_needed,bt_cont_needed
 type(logical_needed)    :: alfa_needed
 type(logical_needed)    :: detr_g_needed,detr_Gama_needed
end type setspacetime
type(setspacetime),save  :: setgr
!=============================================================================
}
{^IFGR
contains
}
{^IFGR
INCLUDE: Metric/mod_metric.^MTR^COORD_^GM.t
}
{^IFGR
!=============================================================================
subroutine setmetric_initphys_spacetime
 include 'amrvacdef.f'


! .. local variable ...
 
!-----------------------------------------------------------------------------
  setgr%a                  = eqpar(gr_a_)
  setgr%rs                 = eqpar(gr_rs_)
end subroutine setmetric_initphys_spacetime
!=============================================================================
subroutine setmetric_init_spacetime
 include 'amrvacdef.f'
!-----------------------------------------------------------------------------
  setgr%normalized             = mtrnormalized
  setgr%diag                   = .false.
  setgr%diagspace              = .false.
  setgr%check_on               = .true.
  setgr%check_isone            = .true.
  setgr%nostatic               = .false.

  setgr%ndir                   = ^NC

  setgr%elem_needed%el         = .false.
  setgr%elem_inv_needed%el     = .false.


  setgr%elemsp_needed%el       = .true.
  setgr%elemsp_inv_needed%el   = .true.
  setgr%elemsp_updwn_needed%el = .true.

  setgr%Qsp_needed%el          = .true.
  setgr%Qsp_inv_needed%el      = .true.

  setgr%bt_needed%el           = .false.
  setgr%bt_cont_needed%el      = .true.

  setgr%alfa_needed%el         = .true.

  setgr%detr_g_needed%el       = .true.
  setgr%detr_Gama_needed%el    = .true.

  call Levi_Civita3  
  call Kronecker_symbol3
  call index_sym_metric
  call setmetric_initalloc_spacetime
  
end subroutine setmetric_init_spacetime
!=============================================================================
subroutine setmetric_initalloc_spacetime
  include 'amrvacdef.f'
  allocate(pm_cell(1:ngridshi),pm_face^D(1:ngridshi),&
           pmCoarse(1:ngridshi),pmCoarse_face^D(1:ngridshi))
end subroutine setmetric_initalloc_spacetime
!=============================================================================
subroutine setmetric_clean_spacetime
 deallocate(pm_cell,pm_face^D,pmCoarse,pmCoarse_face^D)
end subroutine setmetric_clean_spacetime
!=============================================================================
subroutine setmetric_grid(igrid)
 include 'amrvacdef.f'
 integer, intent(in) :: igrid
 ! ... local variables ...
 integer             :: ixCoG^L 
!-----------------------------------------------------------------------------
 saveigrid=igrid
 pm_cell(igrid)%names='pm_cell'
 {^D&write(pm_face^D(igrid)%names,'(a,I)')'pm_face',^D ;}
 pmCoarse(igrid)%names='pmCoarse'
 call setmetric_themetric(ixG^LL,ixG^LL,px(igrid)%x,pm_cell(igrid))
 call setmetric_face(igrid,px(igrid)%x)
 ixCoGmin^D=1; ixCoGmax^D=ixGhi^D/2+dixB;
 call setmetric_themetric(ixCoG^L,ixCoG^L,pxCoarse(igrid)%x,pmCoarse(igrid))
end subroutine setmetric_grid
! ============================================================================
subroutine setmetric_face(igrid,x)
 include 'amrvacdef.f'

 integer, intent(in)          :: igrid
 double precision, intent(in) :: x(ixG^T,1:ndim)

 double precision :: xC(ixG^T,1:ndim),xCoC(ixG^T,1:ndim)
 double precision :: dx^D,xmin^D,xshift^D,xCoshift^D
 integer          :: idims, ixC^L, ixCoG^L,ixCoC^L, ix, idims2,ishift
!-----------------------------------------------------------------------------
 dx^D=rnode(rpdx^D_,igrid);
 xmin^D=rnode(rpxmin^D_,igrid);
 ishift=max(1,dixB/2)
 do idims=1,ndim
   ixCmin^D=ixMlo^D-1-ishift*kr(^D,idims); ixCmax^D=ixMhi^D+1+ishift*kr(^D,idims);
   xshift^D=half*(one-kr(^D,idims));
   do idims2=1,ndim
      select case(idims2)
      {case(^D)
        do ix = ixC^LIM^D
          xC(ix^D%ixC^S,^D)=xmin^D+(dble(ix-dixB)-xshift^D)*dx^D
        end do\}
      end select
   end do

   select case(idims)
   {case (^D)
      call setmetric_themetric(ixG^LL,ixC^L,xC,pm_face^D(igrid)) 
   \}
   end select
  
   ixCoGmin^D=ixGlo^D; ixCoGmax^D=ixGhi^D/2+dixB;
   ixCoCmin^D=ixCoGmin^DD+dixB/2-kr(^DD,^D); ixCoCmax^DD=ixCoGmax^DD-dixB/2;
   xCoshift^D=(one-kr(^D,idims));
   do idims2=1,ndim
      select case(idims2)
      {case(^D)
        do ix = ixCoC^LIM^D
          xCoC(ix^D%ixCoC^S,^D)=xmin^D+(dble(ix-dixB/2)-xCoshift^D)*dx^D
        end do\}
      end select
   end do

   select case(idims)
   {case (^D)
      call setmetric_themetric(ixG^LL,ixCoC^L,xCoC,pmCoarse_face^D(igrid))
    \}
   end select


 end do

end subroutine setmetric_face
 
!=============================================================================
subroutine setmetric_themetric(ixI^L,ixO^L,x,pm)
 use thecom
 include 'amrvacdef.f'
 integer, intent(in)               :: ixI^L,ixO^L
 double precision, intent(in)      :: x(ixI^S,1:^ND)
 type(themetric), intent(inout)    :: pm
!----------------------------------------
 call setmetric_default_init(pm)
 call setmetric_initialize_^MTR^COORD_^GM(pm)
 call setmetric_complete_metric(pm)
 call setmetric_allocstructure_grid(ixI^L,pm)
 call setmetric_spacetime_^MTR^COORD_^GM(ixI^L,ixO^L,x,pm)
 if(setgr%check_on)call setmetric_check_on(ixI^L,ixO^L,pm)
 if(setgr%check_isone)call setmetric_check_isone(ixI^L,ixO^L,pm)
 call setmetric_sqrt_space_determinant_nondiag(ixI^L,ixO^L,&
                                 pm%elemsp,pm%sqrtdetr%Gama)
 if(setgr%detr_g_needed%el)call setmetric_sqrt_4D_determinant_nondiag(ixI^L,&
                 ixO^L,pm%alfa,pm%sqrtdetr%Gama,pm%sqrtdetr%g)

 call setmetric_space_inv(ixI^L,ixO^L,pm%elemsp,pm%sqrtdetr%Gama,&
                          pm%elemsp_inv)
 call setmetric_check_needed_element(pm%elemsp_inv,1,&
                                  setgr%elemsp_inv_needed&
                                    ,.false.)
 if(setgr%elem_inv_needed%el)then
      call setmetric_spacetime_inv(ixI^L,ixO^L,pm%elem,&
                              pm%elemsp_inv,pm%alfa%wg,pm%bt,&
                              pm%elem_inv)
      call setmetric_check_needed_element(pm%elem_inv,1,&
                                  setgr%elem_inv_needed&
                                    ,.false.)
  end if
 if(setgr%diagspace)then
   call setmeric_transfer_matrixsp_diag(ixI^L,ixO^L,pm%elemsp,pm%Qsp)
   call setmeric_transfer_matrixspinv_diag(ixI^L,ixO^L,pm%Qsp,pm%Qsp_inv)
 else
   !call mpistop("--transfer_matrix is not implimented for no diagonal metric--")
 end if

 if(all(dabs(pm%sqrtdetr%Gama(ixO^S)-one)<smalldouble))then
  pm%sqrtdetr%Gama_isone=.true.
 else
  pm%sqrtdetr%Gama_isone=.false.
 end if
 if(all(dabs(pm%sqrtdetr%g(ixO^S)-one)<smalldouble))then
  pm%sqrtdetr%g_isone=.true.
 else
  pm%sqrtdetr%g_isone=.false.
 end if
end subroutine setmetric_themetric
!=============================================================================
subroutine setmetric_default_init(pm)
 type(themetric) , intent(inout) :: pm
! .. local variables ..
 integer                         :: i,j,k
!-------------------------------------
 pm%alfa%elm_on=.false.
 pm%alfa%elm_isone=.false.
 pm%alfa%elm_needed=.true.
Loop_i : do i=1,setgr%ndir
  pm%bt_cont(i)%sftpt%elm_on=.false.
  pm%bt_cont(i)%sftpt%elm_isone=.false.
  pm%bt_cont(i)%sftpt%elm_needed=.true.
  if(pm%cell_center) then
   pm%alfa%drv(i)%elm_on=.false.
   pm%alfa%drv(i)%elm_isone=.false.
   pm%alfa%drv(i)%elm_needed=.true.
  end if
  Loop_j : do j=i,setgr%ndir
   pm%elemsp(i,j)%elpt%elm_on=.false.
   pm%elemsp(i,j)%elpt%elm_isone=.false.
   pm%elemsp(i,j)%elpt%elm_needed=.true.
   pm%Qsp(i,j)%elpt%elm_on=.false.
   pm%Qsp(i,j)%elpt%elm_isone=.false.
   if(pm%cell_center) then
    Loop_k1 : do k=1,setgr%ndir
     pm%elemsp(i,j)%elpt%drv(k)%elm_on=.false.
     pm%elemsp(i,j)%elpt%drv(k)%elm_isone=.false.
     pm%elemsp(i,j)%elpt%drv(k)%elm_needed=.true.
    end do Loop_k1
   end if
  end do Loop_j
  if(pm%cell_center) then
   Loop_k : do k=1,setgr%ndir
    pm%bt_cont(i)%sftpt%drv(k)%elm_on=.false.
    pm%bt_cont(i)%sftpt%drv(k)%elm_isone=.false.
    pm%bt_cont(i)%sftpt%drv(k)%elm_needed=.true.
   end do Loop_k
  end if
end do Loop_i
end subroutine setmetric_default_init
!=============================================================================
subroutine setmetric_complete_metric(pm)
 type(themetric) , intent(inout) :: pm
! .. local variables
 integer                         :: i,j
!------------------------------------- 

 if((pm%alfa%elm_on .and. .not.pm%alfa%elm_isone).and.setgr%alfa_needed%el) then
  pm%alfa%elm_needed=.true.
  if(pm%cell_center)then
    Loop_i : do i=1,setgr%ndir
     if(pm%alfa%elm_needed)then
      if(pm%alfa%drv(i)%elm_on .and. .not.pm%alfa%drv(i)%elm_isone)then
        pm%alfa%drv(i)%elm_needed=.true.
       else
       pm%alfa%drv(i)%elm_needed=.false.
      end if
     end if
    end do Loop_i
   end if
 else
  pm%alfa%elm_needed=.false.
 end if
 Loop_i2 : do i=1,setgr%ndir
  if((pm%bt_cont(i)%sftpt%elm_on .and. .not.pm%bt_cont(i)%sftpt%elm_isone)&
     .and.setgr%bt_cont_needed%el) then
   pm%bt_cont(i)%sftpt%elm_needed=.true.
   if(pm%cell_center)then
    Loop_j : do j=1,setgr%ndir
     if(pm%bt_cont(i)%sftpt%drv(j)%elm_on .and. &
       .not.pm%bt_cont(i)%sftpt%drv(j)%elm_isone)then
      pm%bt_cont(i)%sftpt%drv(j)%elm_needed=.true.
     else
      pm%bt_cont(i)%sftpt%drv(j)%elm_needed=.false.
     end if
    end do Loop_j
   end if
  else
   pm%bt_cont(i)%sftpt%elm_needed=.false.
  end if
 end do Loop_i2


 call setmetric_check_needed_element(pm%elemsp,1,setgr%elemsp_needed&
                                    ,pm%cell_center)
 call setmetric_setdefault_element(pm%elemsp_inv,1,setgr%elemsp_inv_needed&
                                  ,.false.)
 call setmetric_setdefault_Kronecker(pm%elemsp_updwn,1&
                                    ,setgr%elemsp_updwn_needed,.false.)
 call setmetric_check_needed_element(pm%Qsp,1,setgr%Qsp_needed,.false.)
 call setmetric_setdefault_element(pm%Qsp_inv,1,setgr%Qsp_inv_needed,.false.)

end subroutine setmetric_complete_metric
!============================================================================
subroutine setmetric_clean_noneed_element(sub_elem,istart,kernel_needed)
  type(pelements_metric), pointer, intent(inout)     :: sub_elem(:,:)
 integer, intent(in)                                :: istart
 type(logical_needed), intent(in)                   :: kernel_needed
 ! .. local variables ..
 integer                                            :: i,j,k
!------------------------------------------------
  Loop_i : do i=istart,setgr%ndir
  Loop_j : do j=i,setgr%ndir

if(associated(sub_elem(i,j)%elpt%drv))then
 if(all(.not.sub_elem(i,j)%elpt%drv(:)%elm_needed))deallocate(sub_elem(i,j)%elpt%drv)
end if
  end do Loop_j
 end do Loop_i

end subroutine setmetric_clean_noneed_element
!============================================================================

subroutine setmetric_setdefault_Kronecker(sub_elem,istart,kernel_needed,drv_needed)
 type(pelements_metric), pointer, intent(inout)     :: sub_elem(:,:)
 integer, intent(in)                                :: istart
 type(logical_needed), intent(in)                   :: kernel_needed
 logical, intent(in)                                :: drv_needed

 ! .. local variables ..
 integer                                            :: i,j,k
!------------------------------------------------

 Loop_i : do i=istart,setgr%ndir
  Loop_j : do j=i,setgr%ndir
   sub_elem(i,j)%elpt%elm_needed=.false.
   if(i==j)then
    sub_elem(i,j)%elpt%elm_isone=.true.
   else
    sub_elem(i,j)%elpt%elm_on=.false.
   end if
   if(drv_needed)then
    Loop_k: do k=istart,setgr%ndir
     sub_elem(i,j)%elpt%drv(k)%elm_needed=.false.
    end do Loop_k
   end if
  end do  Loop_j
 end do Loop_i
end subroutine setmetric_setdefault_Kronecker
!============================================================================
subroutine setmetric_setdefault_element(sub_elem,istart,kernel_needed,drv_needed)
 type(pelements_metric), pointer, intent(inout)     :: sub_elem(:,:)
 integer, intent(in)                                :: istart
 type(logical_needed), intent(in)                   :: kernel_needed
 logical, intent(in)                                :: drv_needed

 ! .. local variables ..
 integer                                            :: i,j,k
!------------------------------------------------
if(kernel_needed%el) then
 Loop_i : do i=istart,setgr%ndir
  Loop_j : do j=i,setgr%ndir
   sub_elem(i,j)%elpt%elm_needed=.true.
   sub_elem(i,j)%elpt%elm_on=.true.
   if(drv_needed)then
    Loop_k : do k=istart,setgr%ndir   
     sub_elem(i,j)%elpt%drv(k)%elm_needed=.true. 
     sub_elem(i,j)%elpt%drv(k)%elm_on=.true.
    end do Loop_k
   end if
  end do Loop_j
 end do Loop_i
end if
end subroutine setmetric_setdefault_element
!============================================================================
subroutine setmetric_check_needed_element(sub_elem,istart,kernel_needed,drv_needed)
 type(pelements_metric), pointer, intent(inout)     :: sub_elem(:,:)
 integer, intent(in)                                :: istart
 type(logical_needed), intent(in)                   :: kernel_needed 
 logical, intent(in)                                :: drv_needed

 ! .. local variables ..
 integer                                            :: i,j,k
!------------------------------------------------
if(kernel_needed%el) then
 Loop_i : do i=istart,setgr%ndir
  Loop_j : do j=i,setgr%ndir
   if(sub_elem(i,j)%elpt%elm_on .and. .not.sub_elem(i,j)%elpt%elm_isone) then
    sub_elem(i,j)%elpt%elm_needed=.true.
    if(drv_needed)then
     Loop_k : do k=istart,setgr%ndir
      if(sub_elem(i,j)%elpt%drv(k)%elm_on .and. .not.sub_elem(i,j)%elpt%drv(k)%elm_isone) then
        sub_elem(i,j)%elpt%drv(k)%elm_needed=.true.
      else
        sub_elem(i,j)%elpt%drv(k)%elm_needed=.false.
      end if
     end do Loop_k
    end if
   else
    sub_elem(i,j)%elpt%elm_needed=.false.
   end if
  end do Loop_j
 end do Loop_i
end if
end subroutine setmetric_check_needed_element
!============================================================================
subroutine setmetric_check_on(ixI^L,ixO^L,pm)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixI^L,ixO^L
 type(themetric), intent(inout)  :: pm
 ! local variables
 integer                         ::i,j,k                   
!------------------------------------------------
if(pm%alfa%elm_needed)then
 if(maxval(dabs(pm%alfa%wg(ixO^S)))<smalldouble)then
  pm%alfa%elm_on=.false.
  pm%alfa%elm_needed=.false.
  deallocate(pm%alfa%wg)
 end if
 Loop_i :do i=1,setgr%ndir
 
  if(pm%cell_center)then
   if(pm%alfa%drv(i)%elm_on)then
    if(maxval(dabs(pm%alfa%drv(i)%dwg(ixO^S)))<smalldouble)then
      pm%alfa%drv(i)%elm_on=.false.
      pm%alfa%drv(i)%elm_needed=.false.
     deallocate(pm%alfa%drv(i)%dwg)
    end if
   end if
  end if
  if(pm%bt_cont(i)%sftpt%elm_needed) then
   Loop_jbt: do j=1,setgr%ndir
    if(pm%cell_center)then
     if(pm%bt_cont(i)%sftpt%drv(j)%elm_on) then
      if(maxval(dabs(pm%bt_cont(i)%sftpt%drv(j)%dwg(ixO^S)))<smalldouble)then
       pm%bt_cont(i)%sftpt%drv(j)%elm_on=.false.
       pm%bt_cont(i)%sftpt%drv(j)%elm_needed=.false.
       deallocate(pm%bt_cont(i)%sftpt%drv(j)%dwg)
      end if
     end if
    end if
   end do Loop_jbt
   if(maxval(dabs(pm%bt_cont(i)%sftpt%wg(ixO^S)))<smalldouble)then
    pm%bt_cont(i)%sftpt%elm_on=.false.
    pm%bt_cont(i)%sftpt%elm_needed=.false.
    deallocate(pm%bt_cont(i)%sftpt%wg)
   end if
  end if
  
  Loop_j : do j=i,setgr%ndir
   if(pm%elemsp(i,j)%elpt%elm_needed) then
    if(pm%cell_center)then
     Loop_k : do k=1,setgr%ndir
      if(pm%elemsp(i,j)%elpt%drv(k)%elm_needed) then
       if(maxval(dabs(pm%elemsp(i,j)%elpt%drv(k)%dwg(ixO^S)))<smalldouble)then
        pm%elemsp(i,j)%elpt%drv(k)%elm_on=.false.
        pm%elemsp(i,j)%elpt%drv(k)%elm_needed=.false.
        deallocate(pm%elemsp(i,j)%elpt%drv(k)%dwg)
       end if
      end if
     end do Loop_k
    end if
    if(maxval(dabs(pm%elemsp(i,j)%elpt%wg(ixO^S)))<smalldouble)then
     pm%elemsp(i,j)%elpt%elm_on=.false.
     pm%elemsp(i,j)%elpt%elm_needed=.false.
     deallocate(pm%elemsp(i,j)%elpt%wg)
    end if
   end if
  end do Loop_j
 end do Loop_i
end if
end subroutine setmetric_check_on
!============================================================================
subroutine setmetric_check_isone(ixI^L,ixO^L,pm)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixI^L,ixO^L
 type(themetric), intent(inout)  :: pm
 ! .. local variables ..
 integer                         :: i,j,k               
!------------------------------------------------
if(pm%alfa%elm_needed.and.pm%alfa%elm_on)then
 if(maxval(dabs(pm%alfa%wg(ixO^S)-1.0d0))<smalldouble)then
  pm%alfa%elm_isone=.true.
  pm%alfa%elm_needed=.false.
  deallocate(pm%alfa%wg)
 end if
 Loop_i : do i=1,setgr%ndir
  if(pm%bt_cont(i)%sftpt%elm_needed.and.pm%bt_cont(i)%sftpt%elm_on) then
   if(maxval(dabs(pm%bt_cont(i)%sftpt%wg(ixO^S)-1.0d0))<smalldouble)then
    pm%bt_cont(i)%sftpt%elm_isone=.true.
    pm%bt_cont(i)%sftpt%elm_needed=.false.
    deallocate(pm%bt_cont(i)%sftpt%wg)
   end if
  end if
  Loop_j : do j=i,setgr%ndir
   if(pm%cell_center)then
    Loop_k : do k=1,setgr%ndir
     if(pm%elemsp(i,j)%elpt%drv(k)%elm_needed.and.pm%elemsp(i,j)%elpt%drv(k)%elm_on) then
      if(maxval(dabs(pm%elemsp(i,j)%elpt%drv(k)%dwg(ixO^S)-1.0d0))&
         <smalldouble)then
       pm%elemsp(i,j)%elpt%drv(k)%elm_isone=.true.
       pm%elemsp(i,j)%elpt%drv(k)%elm_needed=.false.
       deallocate(pm%elemsp(i,j)%elpt%drv(k)%dwg)
      end if
     end if
    end do Loop_k
   end if
   if(pm%elemsp(i,j)%elpt%elm_needed.and.pm%elemsp(i,j)%elpt%elm_on) then
    if(maxval(dabs(pm%elemsp(i,j)%elpt%wg(ixO^S)-1.0d0))<smalldouble)then
     pm%elemsp(i,j)%elpt%elm_isone=.true.
     pm%elemsp(i,j)%elpt%elm_needed=.false.
     deallocate(pm%elemsp(i,j)%elpt%wg)
    end if
   end if
  end do Loop_j
 end do Loop_i
end if
    
end subroutine setmetric_check_isone
!============================================================================
!============================================================================
subroutine setmetric_multipelement(ixI^L,ixO^L,elem,ldetr)
 integer, intent(in)                :: ixI^L,ixO^L
 type(elements_metric), intent(in)  :: elem(:)
 double precision, intent(inout)    :: ldetr(:^D&)
 ! .. local variables ..
 integer                            :: i, j
!-------------------------------------------------------------------
  if(all(elem(:)%elm_isone))then
    ldetr(ixO^S)=1.0d0
  else
   Loop_i : do i=1,setgr%ndir
    if(.not.elem(i)%elm_isone)then
     ldetr(ixO^S)=elem(i)%wg(ixO^S)
     do j=i+1,setgr%ndir
       if(.not.elem(j)%elm_isone)then
        ldetr(ixO^S)=ldetr(ixO^S)*elem(i)%wg(ixO^S)
       end if
     end do
     exit Loop_i
    end if
   end do Loop_i
  end if
end subroutine setmetric_multipelement

!============================================================================
subroutine setmetric_sqrt_4D_determinant_nondiag(ixI^L,ixO^L,alpha,sqrtdet_gamma,sqrtdet_g)
 include 'amrvacdef.f'
 integer, intent(in)                :: ixI^L,ixO^L
 type(lapse_metric), intent(in)     :: alpha
 double precision, intent(in)       :: sqrtdet_gamma(ixI^S)
 double precision, intent(inout)    :: sqrtdet_g(ixI^S)
!-----------------------------------------------------------------------------
if(.not.alpha%elm_isone)then
  sqrtdet_g(ixO^S)=-dsqrt(alpha%wg(ixO^S))*sqrtdet_gamma(ixO^S)
else
  sqrtdet_g(ixO^S)=-sqrtdet_gamma(ixO^S)
end if
end subroutine setmetric_sqrt_4D_determinant_nondiag
}
{^IFGR
!============================================================================
subroutine setmetric_sqrt_space_determinant_nondiag(ixI^L,ixO^L,metsp,sqrtdetr)
 include 'amrvacdef.f'
 integer, intent(in)                 :: ixI^L,ixO^L
 type(pelements_metric), intent(in)  :: metsp(:,:)
 double precision, intent(out)       :: sqrtdetr(ixI^S)
! .. local variable ...
 double precision                    :: ldetr(ixG^T),detr(ixG^T)
 logical                             :: first
 integer                             :: ii,jj,kk,idir,mdir(setgr%ndir)
!-----------------------------------------------------------------------------
first=.true.
do ii=1,setgr%ndir;mdir(1)=ii
 do jj=1,setgr%ndir-1; mdir(2)=jkdir(mdir(1),jj);
   if(setgr%ndir>2)then
     kk=merge(jj+1,jj+2-ndir,jj+2<=ndir);
     mdir(3)=jkdir(mdir(1),kk)
    end if
  if(all((/(metsp(i1_gr(idir,mdir(idir)),i2_gr(idir,mdir(idir)))%elpt%elm_on,idir=1,setgr%ndir)/))&
       .and.LC3_eps(mdir(1),mdir(2),mdir(3))/=0)then
  call setmetric_multipelement_detr(ixI^L,ixO^L,setgr%ndir,&
     (/ (metsp(i1_gr(idir,mdir(idir)),i2_gr(idir,mdir(idir)))%elpt,idir=1,setgr%ndir)/),ldetr)
 if(LC3_eps(mdir(1),mdir(2),mdir(3))==1)then
     if(first)then
          detr(ixO^S)=ldetr(ixO^S)
     else
          detr(ixO^S)=detr(ixO^S)+ldetr(ixO^S)
     end if
  elseif(LC3_eps(mdir(1),mdir(2),mdir(3))==-1)then
        if(first)then
          detr(ixO^S)=-ldetr(ixO^S)
        else
         detr(ixO^S)=detr(ixO^S)-ldetr(ixO^S)
        end if
  end if
      first=.false.
    end if
 end do
end do
sqrtdetr(ixO^S)=dsqrt(detr(ixO^S))
end subroutine setmetric_sqrt_space_determinant_nondiag
!============================================================================
 subroutine setmetric_multipelement_detr(ixI^L,ixO^L,nelem,elem,ldetr)
 include 'amrvacdef.f'
 integer, intent(in)                  :: ixI^L ,ixO^L,nelem
 type(elements_metric), intent(in)   :: elem(:)
 double precision, intent(inout)      :: ldetr(ixG^T)
 ! .. local variables ..
 integer                              :: i,j
!-------------------------------------------------------------------
  if(all(elem(:)%elm_isone))then
    ldetr=one
  else
   Loop_i : do i=1,nelem
    if(.not.elem(i)%elm_isone)then
     ldetr(ixO^S)=elem(i)%wg(ixO^S)
     Loop_j : do j=i+1,nelem
       if(.not.elem(j)%elm_isone)then
        ldetr(ixO^S)=ldetr(ixO^S)*elem(j)%wg(ixO^S)
       end if
     end do Loop_j
     exit Loop_i
    end if
   end do Loop_i
  end if
end subroutine setmetric_multipelement_detr

!============================================================================
subroutine setmetric_space_inv(ixI^L,ixO^L,metsp,sqrtGama,metsp_inv)
 include 'amrvacdef.f'
 integer, intent(in)                  :: ixI^L,ixO^L
 type(pelements_metric), intent(in)   :: metsp(:,:)
 double precision, intent(in)         :: sqrtGama(ixI^S)
 type(pelements_metric), intent(inout):: metsp_inv(:,:)
 ! .. local variables ..
 integer                              :: i,j,k
 integer                              :: ii(1:setgr%ndir-1),jj(1:setgr%ndir-1)
 double precision,dimension(ixG^T)    :: ldetr,ldetr2
!---------------------------------------
Loop_i : do i=1,setgr%ndir; ii=(/(jkdir(i,k),k=1,setgr%ndir-1)/)
 Loop_j : do j=i,setgr%ndir; jj=(/(jkdir(j,k),k=1,setgr%ndir-1)/)

  if((.not.metsp(i1_gr(ii(1),jj(1)),i2_gr(ii(1),jj(1)))%elpt%elm_on&
        .or..not.metsp(i1_gr(ii(2),jj(2)),i2_gr(ii(2),jj(2)))%elpt%elm_on)&
      .and.(.not.metsp(i1_gr(ii(1),jj(2)),i2_gr(ii(1),jj(2)))%elpt%elm_on&
       .or..not.metsp(i1_gr(ii(2),jj(1)),i2_gr(ii(2),jj(1)))%elpt%elm_on))then
       metsp_inv(i,j)%elpt%elm_on=.false.
       cycle
  end if
  if(metsp(i1_gr(ii(1),jj(1)),i2_gr(ii(1),jj(1)))%elpt%elm_on&
   .and.metsp(i1_gr(ii(2),jj(2)),i2_gr(ii(2),jj(2)))%elpt%elm_on.and.&
     metsp(i1_gr(ii(2),jj(1)),i2_gr(ii(2),jj(1)))%elpt%elm_on&
     .and.metsp(i1_gr(ii(1),jj(2)),i2_gr(ii(1),jj(2)))%elpt%elm_on)then
   call setmetric_multipelement_detr(ixI^L,ixO^L,2,&
      (/metsp(i1_gr(ii(1),jj(1)),i2_gr(ii(1),jj(1)))%elpt,&
       metsp(i1_gr(ii(2),jj(2)),i2_gr(ii(2),jj(2)))%elpt/),ldetr)
   call setmetric_multipelement_detr(ixI^L,ixO^L,2,&
                    (/metsp(i1_gr(ii(2),jj(1)),i2_gr(ii(2),jj(1)))%elpt,&
                    metsp(i1_gr(ii(1),jj(2)),i2_gr(ii(1),jj(2)))%elpt/),ldetr2)
   metsp_inv(i,j)%elpt%wg(ixO^S)=(ldetr(ixO^S)-ldetr2(ixO^S))/sqrtGama(ixO^S)**2.0D0
  else if(.not.metsp(i1_gr(ii(1),jj(1)),i2_gr(ii(1),jj(1)))%elpt%elm_on.or.&
     .not.metsp(i1_gr(ii(2),jj(2)),i2_gr(ii(2),jj(2)))%elpt%elm_on)then
    call setmetric_multipelement_detr(ixI^L,ixO^L,2,&
                 (/metsp(i1_gr(ii(2),jj(1)),i2_gr(ii(2),jj(1)))%elpt,&
                 metsp(i1_gr(ii(1),jj(2)),i2_gr(ii(1),jj(2)))%elpt/),ldetr)
    metsp_inv(i,j)%elpt%wg(ixO^S)=-ldetr(ixO^S)/sqrtGama(ixO^S)**2.0d0
  else if(.not.metsp(i1_gr(ii(1),jj(2)),i2_gr(ii(1),jj(2)))%elpt%elm_on&
    .or..not.metsp(i1_gr(ii(2),jj(1)),i2_gr(ii(2),jj(1)))%elpt%elm_on)then
    call setmetric_multipelement_detr(ixI^L,ixO^L,2,&
              (/metsp(i1_gr(ii(1),jj(1)),i2_gr(ii(1),jj(1)))%elpt,&
                metsp(i1_gr(ii(2),jj(2)),i2_gr(ii(2),jj(2)))%elpt/),ldetr)
    metsp_inv(i,j)%elpt%wg(ixO^S)=ldetr(ixO^S)/sqrtGama(ixO^S)**2.d0
  end if
  if(any(dabs(metsp_inv(i,j)%elpt%wg(ixO^S))>smalldouble))then
   metsp_inv(i,j)%elpt%elm_on=.true.
  else
   metsp_inv(i,j)%elpt%elm_on=.false.
   deallocate(metsp_inv(i,j)%elpt%wg)
   if(allocated(metsp_inv(i,j)%elpt%drv))then
    Loop_k : do k=1,setgr%ndir
     if(allocated(metsp_inv(i,j)%elpt%drv(k)%dwg))&
              deallocate(metsp_inv(i,j)%elpt%drv(k)%dwg)
    end do Loop_k
    deallocate(metsp_inv(i,j)%elpt%drv)
   end if
  end if
 end do Loop_j
end do Loop_i
end subroutine setmetric_space_inv
!============================================================================
subroutine setmetric_spacetime_inv(ixI^L,ixO^L,met4D,metsp_inv,&
                                   alphaw,beta_inv,met4D_inv)
 
 include 'amrvacdef.f'
 integer, intent(in)                   :: ixI^L,ixO^L
 type(pelements_metric), intent(in)    :: met4D(0:,0:)
 type(pelements_metric), intent(in)    :: metsp_inv(:,:)
 type(pshift_metric), intent(in)       :: beta_inv(:)
 double precision, intent(in)          :: alphaw(ixI^S)
 type(pelements_metric), intent(inout) :: met4D_inv(0:,0:)

 ! .. local variables ..
 integer                               :: i, j
!---------------------------------------
met4D_inv(0,0)%elpt%wg(ixO^S)=-one/alphaw(ixO^S)
{^D&if(met4D(0,^D)%elpt%elm_on)then
   if(beta_inv(^D)%sftpt%elm_on)then
     met4D_inv(0,^D)%elpt%elm_on=.true.
     allocate(met4D_inv(0,^D)%elpt%wg(ixG^T))
     met4D_inv(0,^D)%elpt%wg(ixO^S)=beta_inv(^D)%sftpt%wg(ixO^S)/alphaw(ixO^S);
   else
     met4D_inv(0,^D)%elpt%elm_on=.false.
   end if
 end if\}
 Loop_i : do i=1,setgr%ndir
  Loop_j : do j=i,setgr%ndir
   if(.not.metsp_inv(i,j)%elpt%elm_on) then
    if(.not.beta_inv(i)%sftpt%elm_on.or..not.beta_inv(j)%sftpt%elm_on)then
     met4D_inv(i,j)%elpt%elm_on=.false.
     cycle Loop_j
    elseif(beta_inv(i)%sftpt%elm_on.and.beta_inv(j)%sftpt%elm_on)then
     met4D_inv(i,j)%elpt%elm_on=.true.
     allocate(met4D_inv(i,j)%elpt%wg(ixG^T))
      met4d_inv(i,j)%elpt%wg(ixO^S)=-beta_inv(i)%sftpt%wg(ixO^S)&
                                    *beta_inv(j)%sftpt%wg(ixO^S)&
                                     / alphaw(ixO^S);
    end if
   else
     met4D_inv(i,j)%elpt%elm_on=.true.
     allocate(met4D_inv(i,j)%elpt%wg(ixG^T))
     met4d_inv(i,j)%elpt%wg(ixO^S)=metsp_inv(i,j)%elpt%wg(ixO^S)
     if(beta_inv(i)%sftpt%elm_on.and.beta_inv(j)%sftpt%elm_on)then
     met4d_inv(i,j)%elpt%wg(ixO^S)=     met4d_inv(i,j)%elpt%wg(ixO^S)&
       -beta_inv(i)%sftpt%wg(ixO^S)*beta_inv(j)%sftpt%wg(ixO^S)/alphaw(ixO^S);
     end if
   end if
  end do Loop_j
 end do Loop_i
end subroutine setmetric_spacetime_inv
!============================================================================
subroutine setmetric_space_updown(ixI^L,ixO^L,metsp,met_inv,metsp_updwn)
include 'amrvacdef.f'
 integer, intent(in)                    :: ixI^L,ixO^L
 type(pelements_metric), intent(in)     :: metsp(:,:)
 type(pelements_metric), intent(in)     :: met_inv(:,:)
 type(pelements_metric), intent(inout)  :: metsp_updwn(:,:)

 ! .. local variables ..
 integer                                :: iup,jdown,k,addk
 integer :: iupkdir,jdownkdir,kiupdir,kjdowndir,idirup,jdirdown
 integer :: iupaddkdir,jdownaddkdir,addkiupdir,addkjdowndir
!---------------------------------------
 Loop_iup : do iup=1,setgr%ndir
  Loop_jdown : do jdown=1,setgr%ndir
   Loop_k : do k=1,setgr%ndir
    iupkdir=i1_gr(iup,k);jdownkdir=i2_gr(k,jdown);
    kiupdir=i2_gr(iup,k);kjdowndir=i1_gr(k,jdown)
    idirup=i1_gr(iup,jdown);jdirdown=i2_gr(iup,jdown)
    if(met_inv(iupkdir,kiupdir)%elpt%elm_on.and.metsp(kjdowndir,jdownkdir)%elpt%elm_on)then
     if(met_inv(iupkdir,kiupdir)%elpt%elm_isone.and.metsp(kjdowndir,jdownkdir)%elpt%elm_isone)then
      metsp_updwn(idirup,jdirdown)%elpt%wg(ixO^S)=one
     else if(met_inv(iupkdir,kiupdir)%elpt%elm_isone.and..not.metsp(kjdowndir,jdownkdir)%elpt%elm_isone)then
      metsp_updwn(idirup,jdirdown)%elpt%wg(ixO^S)=metsp(kjdowndir,jdownkdir)%elpt%wg(ixO^S)
     else if(.not.met_inv(iupkdir,kiupdir)%elpt%elm_isone.and.metsp(kjdowndir,jdownkdir)%elpt%elm_isone)then
      metsp_updwn(idirup,jdirdown)%elpt%wg(ixO^S)=met_inv(iupkdir,kiupdir)%elpt%wg(ixO^S)
     else if(.not.met_inv(iupkdir,kiupdir)%elpt%elm_isone.and..not.metsp(kjdowndir,jdownkdir)%elpt%elm_isone)then
      metsp_updwn(idirup,jdirdown)%elpt%wg(ixO^S)=met_inv(iupkdir,kiupdir)%elpt%wg(ixO^S)*metsp(k,jdown)%elpt%wg(ixO^S)
     end if
     if(k==ndir)exit Loop_k
     Loop_addk : do addk=k+1,setgr%ndir
      iupaddkdir=i1_gr(iup,k);jdownaddkdir=i2_gr(k,jdown);
      addkiupdir=i2_gr(iup,k);addkjdowndir=i1_gr(k,jdown)
      idirup=i1_gr(iup,jdown);jdirdown=i2_gr(iup,jdown)
      if(met_inv(iupaddkdir,addkiupdir)%elpt%elm_on&
         .and.metsp(addkjdowndir,jdownaddkdir)%elpt%elm_on)then
         if(met_inv(iupaddkdir,addkiupdir)%elpt%elm_isone&
          .and.metsp(addkjdowndir,jdownaddkdir)%elpt%elm_isone)then
      metsp_updwn(idirup,jdirdown)%elpt%wg(ixO^S)=metsp_updwn(idirup,jdirdown)%elpt%wg(ixO^S)+one
     else if(met_inv(iupaddkdir,addkiupdir)%elpt%elm_isone&
          .and..not.metsp(addkjdowndir,jdownaddkdir)%elpt%elm_isone)then
      metsp_updwn(idirup,jdirdown)%elpt%wg(ixO^S)&
                =metsp_updwn(idirup,jdirdown)%elpt%wg(ixO^S)+&
                 metsp(addkjdowndir,jdownaddkdir)%elpt%wg(ixO^S)
     else if(.not.met_inv(iupaddkdir,addkiupdir)%elpt%elm_isone&
             .and.metsp(addkjdowndir,jdownaddkdir)%elpt%elm_isone)then
      metsp_updwn(iup,jdown)%elpt%wg(ixO^S)=&
                 metsp_updwn(idirup,jdirdown)%elpt%wg(ixO^S)+&
                 met_inv(iupaddkdir,addkiupdir)%elpt%wg(ixO^S)
     else if(.not.met_inv(iupaddkdir,addkiupdir)%elpt%elm_isone&
          .and..not.metsp(addkjdowndir,jdownaddkdir)%elpt%elm_isone)then
       metsp_updwn(idirup,jdirdown)%elpt%wg(ixO^S)= &
                 metsp_updwn(idirup,jdirdown)%elpt%wg(ixO^S)+&
                 met_inv(iupaddkdir,addkiupdir)%elpt%wg(ixO^S)&
                 *metsp(addkjdowndir,jdownaddkdir)%elpt%wg(ixO^S)
     end if

      end if
     end do Loop_addk
     exit Loop_k
    end if
   end do Loop_k

  end do Loop_jdown
 end do Loop_iup

end subroutine setmetric_space_updown
!============================================================================
subroutine setmeric_transfer_matrixsp_diag(ixI^L,ixO^L,metsp,Qsp)
 
 integer, intent(in)                    :: ixI^L,ixO^L
 type(pelements_metric), intent(in)     :: metsp(:,:)
 type(pelements_metric), intent(out)    :: Qsp(:,:)
 ! .. local ..
 integer   :: iup,jdwn
 !--------------------------------------------
 Loop_iup : do iup=1,setgr%ndir
  Loop_jdwn : do jdwn=iup,setgr%ndir
   if(iup==jdwn)then
    if(metsp(iup,jdwn)%elpt%elm_on)then
     Qsp(iup,jdwn)%elpt%elm_on     =.true.
     if(metsp(iup,jdwn)%elpt%elm_isone)then
      Qsp(iup,jdwn)%elpt%elm_isone =.true.
     else
       Qsp(iup,jdwn)%elpt%elm_isone =.false.   
       Qsp(iup,jdwn)%elpt%wg(ixO^S) =1.0d0/dsqrt(metsp(iup,jdwn)%elpt%wg(ixO^S))
     end if
    else
     Qsp(iup,jdwn)%elpt%elm_on      = .false.
    end if
   else
    Qsp(iup,jdwn)%elpt%elm_on       = .false.
   end if
  end do Loop_jdwn
 end do Loop_iup
end subroutine setmeric_transfer_matrixsp_diag
!============================================================================
subroutine setmeric_transfer_matrixspinv_diag(ixI^L,ixO^L,Qsp,Qsp_inv)

 integer, intent(in)                    :: ixI^L,ixO^L
 type(pelements_metric), intent(in)     :: Qsp(:,:)
 type(pelements_metric), intent(out)    :: Qsp_inv(:,:)
 ! .. local ..
 integer   :: iup,jdwn
 !--------------------------------------------
 Loop_iup : do iup=1,setgr%ndir
  Loop_jdwn : do jdwn=iup,setgr%ndir
   if(iup==jdwn)then
    if(Qsp(iup,jdwn)%elpt%elm_on)then
     Qsp_inv(iup,jdwn)%elpt%elm_on     =.true.
     if(Qsp(iup,jdwn)%elpt%elm_isone)then
      Qsp_inv(iup,jdwn)%elpt%elm_isone =.true.
     else
       Qsp_inv(iup,jdwn)%elpt%elm_isone =.false.
       Qsp_inv(iup,jdwn)%elpt%wg(ixO^S) =1.0/(Qsp(iup,jdwn)%elpt%wg(ixO^S))
     end if
    else
     Qsp_inv(iup,jdwn)%elpt%elm_on      = .false.
    end if
   else
    Qsp_inv(iup,jdwn)%elpt%elm_on       = .false.
   end if
  end do Loop_jdwn
 end do Loop_iup
end subroutine setmeric_transfer_matrixspinv_diag

!============================================================================
subroutine setmetric_alloc_grid(igrid)

 integer, intent(in)                     :: igrid
 type(themetric)                         :: pm
 !----------------------------------
 pm_cell(igrid)%cell_center=.true.
 call setmetric_alloc_pm(pm_cell(igrid))
 {^D&
 pm_face^D(igrid)%cell_center=.false.;
 call setmetric_alloc_pm(pm_face^D(igrid));\}
 pmCoarse(igrid)%cell_center=.true.
 call setmetric_alloc_pm(pmCoarse(igrid))
 {^D&
 pmCoarse_face^D(igrid)%cell_center=.false.;
 call setmetric_alloc_pm(pmCoarse_face^D(igrid));\}
end subroutine setmetric_alloc_grid
!============================================================================
subroutine setmetric_alloc_pm(pm)
 type(themetric), intent(inout)          :: pm
 ! .. local variable ..
 integer                                 :: i,j
 !----------------------------------
  
 call setmetric_alloc_element(pm%elem,0,setgr%elem_needed,pm%cell_center)
 call setmetric_alloc_element(pm%elem_inv,0,setgr%elem_inv_needed,.false.)
 call setmetric_alloc_element(pm%elemsp,1,setgr%elemsp_needed,pm%cell_center)
 call setmetric_alloc_element(pm%elemsp_inv,1,setgr%elemsp_inv_needed,.false.)
 call setmetric_alloc_element(pm%elemsp_updwn,1,setgr%elemsp_updwn_needed,.false.)

 call setmetric_alloc_element(pm%Qsp,1,setgr%Qsp_needed,.false.)
 call setmetric_alloc_element(pm%Qsp_inv,1,setgr%Qsp_inv_needed,.false.)
 ! the lapse function
 if(setgr%alfa_needed%el)then 
   nullify(pm%alfa)
   allocate(pm%alfa)
   if(pm%cell_center)then
    if(setgr%nostatic)then 
       allocate(pm%alfa%drv(0:setgr%ndir))
    else
       allocate(pm%alfa%drv(1:setgr%ndir))
    end if
    
   end if
 end if
 ! the shift vector
 if(setgr%bt_needed%el)allocate(pm%bt(1:setgr%ndir))
 if(setgr%bt_cont_needed%el)allocate(pm%bt_cont(1:setgr%ndir))
 Loop_i : do i=1,setgr%ndir
  if(setgr%bt_needed%el)then
   nullify(pm%bt(i)%sftpt)
   allocate(pm%bt(i)%sftpt)
   if(pm%cell_center)then
     allocate(pm%bt(i)%sftpt%drv(1:setgr%ndir))
   end if
  end if
  if(setgr%bt_cont_needed%el)then
   nullify(pm%bt_cont(i)%sftpt)
   allocate(pm%bt_cont(i)%sftpt)
   if(pm%cell_center)then 
     allocate(pm%bt_cont(i)%sftpt%drv(1:setgr%ndir))
   end if
  end if
 
 end do Loop_i
 nullify(pm%sqrtdetr)
 allocate(pm%sqrtdetr)
end subroutine setmetric_alloc_pm
!============================================================================
subroutine setmetric_alloc_element(sub_elm,istart,kernel_needed,drv_needed)
 type(pelements_metric), intent(inout),pointer    :: sub_elm(:,:)
 integer, intent(in)                              :: istart
 type(logical_needed), intent(in)                 :: kernel_needed
 logical, intent(in)                              :: drv_needed
 ! .. local ..
 integer                                          :: i,j
 !-----------------------------------------
 if(kernel_needed%el)then
   allocate(sub_elm(istart:setgr%ndir,istart:setgr%ndir))
   Loop_i : do i=istart,setgr%ndir
    Loop_j : do j=i,setgr%ndir
      nullify(sub_elm(i,j)%elpt)
      allocate(sub_elm(i,j)%elpt)
      if(j>=i.and.drv_needed)then
       if(setgr%nostatic)then
        allocate(sub_elm(i,j)%elpt%drv(0:setgr%ndir))
       else
        allocate(sub_elm(i,j)%elpt%drv(1:setgr%ndir))
       end if
      end if
    end do Loop_j
   end do Loop_i
 end if
end subroutine setmetric_alloc_element
!============================================================================
subroutine setmetric_allocstructure_grid(ixI^L,pm)
 integer, intent(in)                     :: ixI^L
 type(themetric), intent(inout)          :: pm
 ! .. local variables ..
 integer         :: i,j,k
!----------------------------------------------------
! allocate memory for the grid in Kerr case
 call setmetric_allocstructure_element(ixI^L,pm%elemsp,1,&
                                    setgr%elemsp_needed,pm%cell_center)
 call setmetric_allocstructure_element(ixI^L,pm%elemsp_inv,1,&
                                    setgr%elemsp_inv_needed,.false.)
 call setmetric_allocstructure_element(ixI^L,pm%elemsp_updwn,1,&
                                    setgr%elemsp_updwn_needed,.False.)
 call setmetric_allocstructure_element(ixI^L,pm%Qsp,1,&
                                    setgr%Qsp_needed,.False.)
 call setmetric_allocstructure_element(ixI^L,pm%Qsp_inv,1,&
                                    setgr%Qsp_inv_needed,.False.)
 call setmetric_allocstructure_shift_vector(ixI^L,pm%bt,setgr%bt_needed,pm%cell_center)
 call setmetric_allocstructure_shift_vector(ixI^L,pm%bt_cont,setgr%bt_cont_needed,pm%cell_center)


 allocate(pm%alfa%wg(ixI^S))
 if(pm%cell_center)then
  do i=1,setgr%ndir
   if(pm%alfa%drv(i)%elm_needed)allocate(pm%alfa%drv(i)%dwg(ixI^S))
  end do
 end if
 allocate(pm%sqrtdetr%Gama(ixI^S),pm%sqrtdetr%g(ixI^S))
end subroutine setmetric_allocstructure_grid
!============================================================================
subroutine setmetric_allocstructure_element(ixI^L,sub_elem,istart,kernel_needed,drv_needed)
 integer, intent(in)                              :: ixI^L
 type(pelements_metric), intent(inout),pointer    :: sub_elem(:,:)
 integer, intent(in)                              :: istart
 type(logical_needed), intent(in)                 :: kernel_needed
 logical, intent(in)                              :: drv_needed
 ! .. local ..
 integer                                          :: i,j,k
 !-----------------------------------------
 if(kernel_needed%el)then
  Loop_i : do i=istart,setgr%ndir
   Loop_j : do j=i,setgr%ndir
   if(sub_elem(i,j)%elpt%elm_needed)then
    allocate(sub_elem(i,j)%elpt%wg(ixI^S))
    if (drv_needed) then
     Loop_k : do k=istart,setgr%ndir
       if(sub_elem(i,j)%elpt%drv(k)%elm_needed)then
        allocate(sub_elem(i,j)%elpt%drv(k)%dwg(ixI^S))
       end if
      end do Loop_k
     end if
    end if
   end do Loop_j
  end do Loop_i
 end if
end subroutine setmetric_allocstructure_element
!============================================================================
subroutine setmetric_allocstructure_shift_vector(ixI^L,sub_beta,kernel_needed,drv_needed)
 integer, intent(in)                              :: ixI^L
 type(pshift_metric), intent(inout)               :: sub_beta(:)
 type(logical_needed), intent(in)                 :: kernel_needed
 logical, intent(in)                              :: drv_needed
 ! .. local ..
 integer                                          :: i,k
 !-----------------------------------------
 if(kernel_needed%el)then
  Loop_i : do i=1,setgr%ndir
    if(sub_beta(i)%sftpt%elm_needed)then
     allocate(sub_beta(i)%sftpt%wg(ixI^S))
     if(drv_needed)then
      Loop_k : do k=1,setgr%ndir
       if(sub_beta(i)%sftpt%drv(k)%elm_needed)then
        allocate(sub_beta(i)%sftpt%drv(k)%dwg(ixI^S))
       end if
      end do Loop_k
     end if
    end if
   end do Loop_i
 end if
end subroutine setmetric_allocstructure_shift_vector
!============================================================================
subroutine setmetric_dealloc_grid(igrid)
 integer, intent(in) :: igrid
 !----------------------------------
 call setmetric_deallocstructure_grid(pm_cell(igrid))
 {^D&call setmetric_deallocstructure_grid(pm_face^D(igrid));}
 call setmetric_deallocstructure_grid(pmCoarse(igrid))
  {^D&call setmetric_deallocstructure_grid(pmCoarse_face^D(igrid));} 

end subroutine setmetric_dealloc_grid
!============================================================================
subroutine setmetric_deallocstructure_grid(pm)
 type(themetric)       :: pm 
 integer               :: i,j,k
!----------------------------------------------------
! clean memory for the grid in Kerr case
 !call setmetric_deassociate_element(pm%elem,0,setgr%elem_needed)
 !call setmetric_deassociate_element(pm%elem_inv,0,setgr%elem_inv_needed)
 !call setmetric_deassociate_element(pm%elemsp,1,setgr%elemsp_needed)
 !call setmetric_deassociate_element(pm%elemsp_inv,1,setgr%elemsp_inv_needed)
 !call setmetric_deassociate_element(pm%elemsp_updwn,1,setgr%elemsp_updwn_needed)

 call setmetric_dealloc_element(pm%elem,0,setgr%elem_needed,pm%cell_center)
 call setmetric_dealloc_element(pm%elem_inv,0,setgr%elem_inv_needed,.false.)
 call setmetric_dealloc_element(pm%elemsp,1,setgr%elemsp_needed,pm%cell_center)
 call setmetric_dealloc_element(pm%elemsp_inv,1,setgr%elemsp_inv_needed,.false.)
 call setmetric_dealloc_element(pm%elemsp_updwn,1,setgr%elemsp_updwn_needed,.false.)

 call setmetric_dealloc_element(pm%Qsp,1,setgr%Qsp_needed,.false.)
 call setmetric_dealloc_element(pm%Qsp_inv,1,setgr%Qsp_inv_needed,.false.)
 ! the shift vector

 call setmetric_clean_shift_vector(pm%bt,setgr%bt_needed,pm%cell_center)
 call setmetric_clean_shift_vector(pm%bt_cont,setgr%bt_cont_needed,pm%cell_center)

 ! the lapse function
 if(allocated(pm%alfa%wg))deallocate(pm%alfa%wg)
 if(pm%cell_center)then
  do i=1,setgr%ndir
   if(allocated(pm%alfa%drv(i)%dwg))deallocate(pm%alfa%drv(i)%dwg)
  end do
 ! the lapse function
  deallocate(pm%alfa%drv)
 end if
 nullify(pm%alfa)
 ! determinant
 deallocate(pm%sqrtdetr%Gama,pm%sqrtdetr%g)
 nullify(pm%sqrtdetr)



end subroutine setmetric_deallocstructure_grid
!============================================================================
subroutine setmetric_deassociate_element(sub_elem,istart,kernel_needed)
 type(pelements_metric), intent(inout),pointer    :: sub_elem(:,:)
 integer, intent(in)                              :: istart
 type(logical_needed), intent(in)                 :: kernel_needed
 ! .. local ..
 integer                                          :: i,j
 !-----------------------------------------
 if(kernel_needed%el)then
  Loop_i1 : do i=istart+1,setgr%ndir
   Loop_j1 : do j=istart,i-1
    nullify(sub_elem(i,j)%elpt) 
    !if(associated(sub_elem(i,j)%elpt))deallocate(sub_elem(i,j)%elpt)
   end do Loop_j1
  end do Loop_i1
 end if 
end subroutine setmetric_deassociate_element
!============================================================================
subroutine setmetric_dealloc_element(sub_elm,istart,kernel_needed,drv_used)
 type(pelements_metric), intent(inout),pointer    :: sub_elm(:,:)
 integer, intent(in)                              :: istart
 type(logical_needed), intent(in)                 :: kernel_needed
 logical, intent(in)                              :: drv_used
 ! .. local ..
 integer                                          :: i,j,k
 !-----------------------------------------
 if(kernel_needed%el)then
   Loop_i : do i=istart,setgr%ndir
    Loop_j : do j=i,setgr%ndir
     if(associated(sub_elm(i,j)%elpt))then
      if(.not.sub_elm(i,j)%elpt%elm_on)cycle Loop_j
      if(drv_used)then
       if(associated(sub_elm(i,j)%elpt%drv))then
        Loop_k : do k=1,setgr%ndir
          if(.not.sub_elm(i,j)%elpt%drv(k)%elm_on)cycle Loop_k
          if(allocated(sub_elm(i,j)%elpt%drv(k)%dwg))&
                   deallocate(sub_elm(i,j)%elpt%drv(k)%dwg)
        end do Loop_k
        deallocate(sub_elm(i,j)%elpt%drv)
      end if
     end if 
       if(sub_elm(i,j)%elpt%elm_needed)deallocate(sub_elm(i,j)%elpt%wg)
       nullify(sub_elm(i,j)%elpt)
     end if
    end do Loop_j
   end do Loop_i
 end if
 if(associated(sub_elm))deallocate(sub_elm)
end subroutine setmetric_dealloc_element
!============================================================================
subroutine setmetric_clean_shift_vector(sub_beta,kernel_needed,drv_used)
 type(pshift_metric),pointer,intent(inout)        :: sub_beta(:)
 type(logical_needed), intent(in)                 :: kernel_needed
 logical, intent(in)                              :: drv_used
 ! .. local variables ..
 integer                                          :: k,i
!--------------------------------------
 if(kernel_needed%el)then
  Loop_i: do i=1,setgr%ndir
   if(allocated(sub_beta(i)%sftpt%wg))then
    deallocate(sub_beta(i)%sftpt%wg)
    if(drv_used) then
     Loop_k : do k=1,setgr%ndir
      if(allocated(sub_beta(i)%sftpt%drv(k)%dwg))deallocate(sub_beta(i)%sftpt%drv(k)%dwg)
     end do Loop_k
     deallocate(sub_beta(i)%sftpt%drv)
    end if
    nullify(sub_beta(i)%sftpt)
   end if
  end do Loop_i
 end if
 if(associated(sub_beta))deallocate(sub_beta)
end subroutine setmetric_clean_shift_vector
!============================================================================
subroutine index_sym_metric

! .. local variable
 integer  :: i,j
!-----------------------------------------------------------------------------

do i=1,3
 do j=i,3
    i1_gr(i,j)=i;i2_gr(i,j)=j
 end do
 if(i>1)then
  do j=1,i-1
    i1_gr(i,j)=j;i2_gr(i,j)=i
  end do
 end if
end do
end subroutine index_sym_metric
!============================================================================
subroutine Levi_Civita3

! .. local variable
 integer  :: i,j,k
!-----------------------------------------------------------------------------
do i=1,3
 do j=1,3
  do k=1,3
    LC3_eps(i,j,k)=((i-j)*(k-i)*(j-k))/2
  end do
 end do
end do
end subroutine Levi_Civita3
!============================================================================
subroutine Kronecker_symbol3
! .. local variable
 integer  :: i,j
!---------------------------------------------------------------
Loop_i : do i=1,setgr%ndir
 Loop_j :do j=1,setgr%ndir
   if(i==j)then
      Eps3(i,j)=1
   else
      Eps3(i,j)=0
   end if
 end do Loop_j
end do Loop_i
end subroutine Kronecker_symbol3

}
!============================================================================
End module mod_metric
!============================================================================

