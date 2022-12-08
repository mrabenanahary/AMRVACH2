subroutine setmetric_init_MblTOks
! .. local ..
 integer :: idir
!-----------------------------------
 allocate(MblTOks(0:^NC,0:^NC))
 MblTOks(:,:)%elm_on    =.false.
 MblTOks(0,0)%elm_on    = .true.
 MblTOks(0,0)%elm_isone = .true.

 do idir=1,^NC
  MblTOks(idir,idir)%elm_on    = .true.
  MblTOks(idir,idir)%elm_isone = .true.
 end do
 MblTOks(0,r_)%elm_on    = .true.
 MblTOks(0,r_)%elm_isone = .false.
 MblTOks(phi_,r_)%elm_on    = .true.
 MblTOks(phi_,r_)%elm_isone = .false.

end subroutine setmetric_init_MblTOks
!=======================================================================
subroutine setmetric_allocMblTOks
include 'amrvacdef.f'
integer :: i,j
!-----------------------------------------------
do i=0,ndir
 do j=0,ndir
   if(MblTOks(i,j)%elm_on)then
    if(.not.MblTOks(i,j)%elm_isone)then
     allocate(MblTOks(i,j)%wg(ixG^T))
    end if
   end if
 end do
end do
end subroutine setmetric_allocMblTOks
!===========================================================================
subroutine setmetric_setMblTOks(ixI^L,ixO^L,patchw,x)
 include 'amrvacdef.f'

 integer, intent(in)             :: ixI^L, ixO^L
 logical, intent(in)             :: patchw(ixI^S)
 double precision,intent(in)     :: x(ixI^S,ndir)
 ! .. local ..
 double precision                :: delta(ixG^T) 
!--------------------------------
 delta(ixO^S) =1.0-eqpar(gr_rs_)/x(ixO^S,r_)+&
               eqpar(gr_a_)**2.0D0/x(ixO^S,r_)**2.0D0
 MblTOks(0,r_)%wg(ixO^S) = eqpar(gr_rs_)/x(ixO^S,r_)/delta(ixO^S)
 MblTOks(phi_,r_)%wg(ixO^S) =eqpar(gr_a_)/delta(ixO^S)

end subroutine setmetric_setMblTOks
!==========================================================================
subroutine setmetric_cleanMblTOks
integer :: i,j
do i=0,ndir
 do j=0,ndir
  if(allocated(MblTOks(i,j)%wg))deallocate(MblTOks(i,j)%wg)
 end do
end do
end subroutine setmetric_cleanMblTOks
!==========================================================================
subroutine setmetric_cleanMblTOks_structure
 deallocate(MblTOks)
end subroutine setmetric_cleanMblTOks_structure
!==========================================================================
{^IFCOORDKS
INCLUDE: Metric/mod_metric.^MTRbl_^GM.t
}
{^IFCOORDKS
subroutine setmetric_BL(ixI^L,ixO^L,x,pm)
 use thecom
 include 'amrvacdef.f'
 integer, intent(in)               :: ixI^L,ixO^L
 double precision, intent(in)      :: x(ixI^S,1:^ND)
 type(themetric), intent(inout)    :: pm
!----------------------------------------
 call setmetric_alloc_pm(pm)
 call setmetric_default_init(pm)
 call setmetric_initialize_^MTRbl_^GM(pm)
 call setmetric_complete_metric(pm)
 call setmetric_allocstructure_grid(ixI^L,pm)
 call setmetric_spacetime_^MTRbl_^GM(ixI^L,ixO^L,x,pm)
 call setmetric_sqrt_space_determinant_nondiag(ixI^L,ixO^L,&
                                 pm%elemsp,pm%sqrtdetr%Gama)
 if(setgr%detr_g_needed%el)call setmetric_sqrt_4D_determinant_nondiag(ixI^L,&
                 ixO^L,pm%alfa,pm%sqrtdetr%Gama,pm%sqrtdetr%g)
 call setmetric_space_inv(ixI^L,ixO^L,pm%elemsp,pm%sqrtdetr%Gama,&
                          pm%elemsp_inv)
end subroutine setmetric_BL
}
!==========================================================================
{^IFCOORDBL
INCLUDE: Metric/mod_metric.^MTRks_^GM.t
}
{^IFCOORDBL
subroutine setmetric_KS(ixI^L,ixO^L,x,pm)
 use thecom
 include 'amrvacdef.f'
 integer, intent(in)               :: ixI^L,ixO^L
 double precision, intent(in)      :: x(ixI^S,1:^ND)
 type(themetric), intent(inout)    :: pm
!----------------------------------------
 call setmetric_alloc_pm(pm)
 call setmetric_default_init(pm)
 call setmetric_initialize_^MTRks_^GM(pm)
 call setmetric_complete_metric(pm)
 call setmetric_allocstructure_grid(ixI^L,pm)
 call setmetric_spacetime_^MTRks_^GM(ixI^L,ixO^L,x,pm)
 call setmetric_sqrt_space_determinant_nondiag(ixI^L,ixO^L,&
                                 pm%elemsp,pm%sqrtdetr%Gama)
 if(setgr%detr_g_needed%el)call setmetric_sqrt_4D_determinant_nondiag(ixI^L,&
                 ixO^L,pm%alfa,pm%sqrtdetr%Gama,pm%sqrtdetr%g)

 call setmetric_space_inv(ixI^L,ixO^L,pm%elemsp,pm%sqrtdetr%Gama,&
                          pm%elemsp_inv)
end subroutine setmetric_KS
}
!==========================================================================
subroutine usemetric_speed4_BLtoKS(ixI^L,ixO^L,patchw,x,subname,uBL,uKS,lfacKS)
 include 'amrvacdef.f'

 integer, intent(in)             :: ixI^L, ixO^L
 logical, intent(in)             :: patchw(ixI^S)
 character(len=*), intent(in)    :: subname
 double precision,intent(in)     :: x(ixI^S,ndir)
 double precision,intent(in)     :: uBL(ixI^S,0:ndir)
 double precision,intent(inout)  :: uKS(ixI^S,0:ndir)
 double precision, intent(out)   :: lfacKS(ixG^T)  
 ! .. local variables ..
 integer                            :: i,j
 double precision, dimension(ixG^T) :: delta
 !----------------------------------------------
 call setmetric_allocMblTOks
 call setmetric_setMblTOks(ixI^L,ixO^L,patchw,x)
 uKS(ixO^S,:) =0.0d0
 do i=0,ndir
  do j=0,ndir
   if(MblTOks(i,j)%elm_on) then
    if(.not.MblTOks(i,j)%elm_isone)then
     uKS(ixO^S,i) = uKS(ixO^S,i)+&
                         MblTOks(i,j)%wg(ixO^S)*uBL(ixO^S,j)
    else
     uKS(ixO^S,i) = uKS(ixO^S,i) + uBL(ixO^S,j)
    end if 
   end if
  end do
 end do
 call usemetric_multiby_alpha(ixO^L,uKS(ixO^S,0),lfacKS)
 call setmetric_cleanMblTOks

end subroutine usemetric_speed4_BLtoKS
