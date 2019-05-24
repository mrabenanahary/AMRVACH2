

subroutine set_B0_grid(igrid,qt)
  use mod_global_parameters

  integer, intent(in)             :: igrid
  double precision, intent(in)    :: qt

  call set_B0_cell(pw(igrid)%B0(:,:,:,0),pw(igrid)%x,ixGlo1,ixGlo2,ixGhi1,&
     ixGhi2,ixGlo1,ixGlo2,ixGhi1,ixGhi2,qt)
  call set_J0_cell(igrid,pw(igrid)%J0,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1-1,&
     ixMlo2-1,ixMhi1+1,ixMhi2+1,qt)
  call set_B0_face(igrid,pw(igrid)%x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,&
     ixMhi1,ixMhi2,qt)

end subroutine set_B0_grid

subroutine set_B0_cell(wB0,x,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,&
   ixmax1,ixmax2,qt)
  use mod_constants
  use mod_usr_methods, only: usr_set_B0
  use mod_global_parameters
  use mod_geometry
  implicit none
  integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,&
     ixmin2,ixmax1,ixmax2
  double precision, intent(inout) :: wB0(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndir)
  double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
  double precision, intent(in)    :: qt
  integer                         :: idims
  real(dp)                        :: mu0_sphe(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim)
  real(dp)                        :: mu0_cart(1:ndim)
  real(dp)                        :: dmutheta,dmuphi
  wB0(ixmin1:ixmax1,ixmin2:ixmax2,1:ndir)=zero

  ! approximate cell-averaged B0 as cell-centered B0
  select case (typeaxial)
  case ("spherical")
     
     if (dabs(Bdip)>smalldouble) then
        wB0(ixmin1:ixmax1,ixmin2:ixmax2,1)=2.0d0*Bdip*dcos(x(ixmin1:ixmax1,&
           ixmin2:ixmax2,2))/x(ixmin1:ixmax1,ixmin2:ixmax2,1)**3
        wB0(ixmin1:ixmax1,ixmin2:ixmax2,2)=Bdip*dsin(x(ixmin1:ixmax1,&
           ixmin2:ixmax2,2))/x(ixmin1:ixmax1,ixmin2:ixmax2,1)**3
     end if
  
     if (abs(Bquad)>smalldouble) then
        wB0(ixmin1:ixmax1,ixmin2:ixmax2,1)=wB0(ixmin1:ixmax1,ixmin2:ixmax2,&
           1) +Bquad*0.5d0*(1.0d0+3.0d0*dcos(2.0d0*x(ixmin1:ixmax1,&
           ixmin2:ixmax2,2)))/x(ixmin1:ixmax1,ixmin2:ixmax2,1)**4
        wB0(ixmin1:ixmax1,ixmin2:ixmax2,2)=wB0(ixmin1:ixmax1,ixmin2:ixmax2,&
           2)+Bquad*dsin(2.0d0*x(ixmin1:ixmax1,ixmin2:ixmax2,&
           2))/x(ixmin1:ixmax1,ixmin2:ixmax2,1)**4
     end if
     if (abs(Boct)>smalldouble) then
        wB0(ixmin1:ixmax1,ixmin2:ixmax2,1)=wB0(ixmin1:ixmax1,ixmin2:ixmax2,&
           1) +Boct*(10.0d0*dcos(2.0d0*x(ixmin1:ixmax1,ixmin2:ixmax2,&
           2))-2.0d0) *dcos(x(ixmin1:ixmax1,ixmin2:ixmax2,2))/x(ixmin1:ixmax1,&
           ixmin2:ixmax2,1)**5
        wB0(ixmin1:ixmax1,ixmin2:ixmax2,2)=wB0(ixmin1:ixmax1,ixmin2:ixmax2,&
           2) +Boct*1.5d0*(3.0d0+5.0d0*dcos(2.0d0*x(ixmin1:ixmax1,&
           ixmin2:ixmax2,2))) *dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,&
           2))/x(ixmin1:ixmax1,ixmin2:ixmax2,1)**5
     end if
     if(dabs(mu0dip)>smalldouble)then
     
     end if
    
  end select
  
  if (associated(usr_set_B0)) call usr_set_B0(ixImin1,ixImin2,ixImax1,ixImax2,&
     ixmin1,ixmin2,ixmax1,ixmax2,x,wB0,qt)

end subroutine set_B0_cell

subroutine set_J0_cell(igrid,wJ0,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,&
   ixmax1,ixmax2,qt)
  use mod_usr_methods, only: usr_set_J0
  use mod_global_parameters
  use mod_geometry

  integer, intent(in)             :: igrid,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixmin1,ixmin2,ixmax1,ixmax2
  double precision, intent(inout) :: wJ0(ixImin1:ixImax1,ixImin2:ixImax2,&
     7-2*ndir:3)
  double precision, intent(in)    :: qt
  integer :: idirmin0, idirmin

  if(associated(usr_set_J0)) then
    call usr_set_J0(ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
       ixmax2,pw(igrid)%x,wJ0,qt)
  else
    idirmin0 = 7-2*ndir
    call curlvector(pw(igrid)%B0(:,:,:,0),ixImin1,ixImin2,ixImax1,ixImax2,&
       ixmin1,ixmin2,ixmax1,ixmax2,wJ0,idirmin,idirmin0,ndir)
  end if

end subroutine set_J0_cell

subroutine set_B0_face(igrid,x,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,&
   ixmax1,ixmax2,qt)
  use mod_global_parameters

  integer, intent(in) :: igrid, ixImin1,ixImin2,ixImax1,ixImax2, ixmin1,ixmin2,&
     ixmax1,ixmax2
  double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
  double precision, intent(in)    :: qt

  double precision :: delx(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
  double precision :: xC(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),xshift1,&
     xshift2
  integer :: idims, ixCmin1,ixCmin2,ixCmax1,ixCmax2, ix, idims2

  if(slab)then
   delx(ixImin1:ixImax1,ixImin2:ixImax2,1)=rnode(rpdx1_,igrid)
   delx(ixImin1:ixImax1,ixImin2:ixImax2,2)=rnode(rpdx2_,igrid)
  else
   ! for all non-cartesian and stretched coordinate(s)
   delx(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)=pw(igrid)%dx(ixImin1:ixImax1,&
      ixImin2:ixImax2,1:ndim)
  endif

  do idims=1,ndim
     ixCmin1=ixmin1-kr(1,idims);ixCmin2=ixmin2-kr(2,idims); ixCmax1=ixmax1
     ixCmax2=ixmax2;
     ! always xshift=0 or 1/2
     xshift1=half*(one-kr(1,idims));xshift2=half*(one-kr(2,idims));
     do idims2=1,ndim
       select case(idims2)
       case(1)
         do ix = ixCmin1,ixCmax1
           ! xshift=half: this is the cell center coordinate
           ! xshift=0: this is the cell edge i+1/2 coordinate
           xC(ix,ixCmin2:ixCmax2,1)=x(ix,ixCmin2:ixCmax2,&
              1)+(half-xshift1)*delx(ix,ixCmin2:ixCmax2,1)
         end do
       case(2)
         do ix = ixCmin2,ixCmax2
           ! xshift=half: this is the cell center coordinate
           ! xshift=0: this is the cell edge i+1/2 coordinate
           xC(ixCmin1:ixCmax1,ix,2)=x(ixCmin1:ixCmax1,ix,&
              2)+(half-xshift2)*delx(ixCmin1:ixCmax1,ix,2)
         end do
       end select
     end do
     call set_B0_cell(pw(igrid)%B0(:,:,:,idims),xC,ixImin1,ixImin2,ixImax1,&
        ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,qt)
  end do

end subroutine set_B0_face

subroutine alloc_B0_grid(igrid)
  use mod_global_parameters

  integer, intent(in) :: igrid

  if(.not. allocated(pw(igrid)%B0)) then
    allocate(pw(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndir,0:ndim))
    allocate(pw(igrid)%J0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,7-2*ndir:3))
  end if

end subroutine alloc_B0_grid

subroutine dealloc_B0_grid(igrid)
  use mod_global_parameters

  integer, intent(in) :: igrid

  deallocate(pw(igrid)%B0)
  deallocate(pw(igrid)%J0)

end subroutine dealloc_B0_grid
