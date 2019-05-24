

subroutine set_B0_grid(igrid,qt)
  use mod_global_parameters

  integer, intent(in)             :: igrid
  double precision, intent(in)    :: qt

  call set_B0_cell(pw(igrid)%B0(:,:,0),pw(igrid)%x,ixGlo1,ixGhi1,ixGlo1,ixGhi1,&
     qt)
  call set_J0_cell(igrid,pw(igrid)%J0,ixGlo1,ixGhi1,ixMlo1-1,ixMhi1+1,qt)
  call set_B0_face(igrid,pw(igrid)%x,ixGlo1,ixGhi1,ixMlo1,ixMhi1,qt)

end subroutine set_B0_grid

subroutine set_B0_cell(wB0,x,ixImin1,ixImax1,ixmin1,ixmax1,qt)
  use mod_constants
  use mod_usr_methods, only: usr_set_B0
  use mod_global_parameters
  use mod_geometry
  implicit none
  integer, intent(in)             :: ixImin1,ixImax1,ixmin1,ixmax1
  double precision, intent(inout) :: wB0(ixImin1:ixImax1,1:ndir)
  double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
  double precision, intent(in)    :: qt
  integer                         :: idims
  real(dp)                        :: mu0_sphe(ixImin1:ixImax1,1:ndim)
  real(dp)                        :: mu0_cart(1:ndim)
  real(dp)                        :: dmutheta,dmuphi
  wB0(ixmin1:ixmax1,1:ndir)=zero

  ! approximate cell-averaged B0 as cell-centered B0
  select case (typeaxial)
  case ("spherical")
     
  end select
  
  if (associated(usr_set_B0)) call usr_set_B0(ixImin1,ixImax1,ixmin1,ixmax1,x,&
     wB0,qt)

end subroutine set_B0_cell

subroutine set_J0_cell(igrid,wJ0,ixImin1,ixImax1,ixmin1,ixmax1,qt)
  use mod_usr_methods, only: usr_set_J0
  use mod_global_parameters
  use mod_geometry

  integer, intent(in)             :: igrid,ixImin1,ixImax1,ixmin1,ixmax1
  double precision, intent(inout) :: wJ0(ixImin1:ixImax1,7-2*ndir:3)
  double precision, intent(in)    :: qt
  integer :: idirmin0, idirmin

  if(associated(usr_set_J0)) then
    call usr_set_J0(ixImin1,ixImax1,ixmin1,ixmax1,pw(igrid)%x,wJ0,qt)
  else
    idirmin0 = 7-2*ndir
    call curlvector(pw(igrid)%B0(:,:,0),ixImin1,ixImax1,ixmin1,ixmax1,wJ0,&
       idirmin,idirmin0,ndir)
  end if

end subroutine set_J0_cell

subroutine set_B0_face(igrid,x,ixImin1,ixImax1,ixmin1,ixmax1,qt)
  use mod_global_parameters

  integer, intent(in) :: igrid, ixImin1,ixImax1, ixmin1,ixmax1
  double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
  double precision, intent(in)    :: qt

  double precision :: delx(ixImin1:ixImax1,1:ndim)
  double precision :: xC(ixImin1:ixImax1,1:ndim),xshift1
  integer :: idims, ixCmin1,ixCmax1, ix, idims2

  if(slab)then
   delx(ixImin1:ixImax1,1)=rnode(rpdx1_,igrid)
  else
   ! for all non-cartesian and stretched coordinate(s)
   delx(ixImin1:ixImax1,1:ndim)=pw(igrid)%dx(ixImin1:ixImax1,1:ndim)
  endif

  do idims=1,ndim
     ixCmin1=ixmin1-kr(1,idims); ixCmax1=ixmax1;
     ! always xshift=0 or 1/2
     xshift1=half*(one-kr(1,idims));
     do idims2=1,ndim
       select case(idims2)
       case(1)
         do ix = ixCmin1,ixCmax1
           ! xshift=half: this is the cell center coordinate
           ! xshift=0: this is the cell edge i+1/2 coordinate
           xC(ix,1)=x(ix,1)+(half-xshift1)*delx(ix,1)
         end do
       end select
     end do
     call set_B0_cell(pw(igrid)%B0(:,:,idims),xC,ixImin1,ixImax1,ixCmin1,&
        ixCmax1,qt)
  end do

end subroutine set_B0_face

subroutine alloc_B0_grid(igrid)
  use mod_global_parameters

  integer, intent(in) :: igrid

  if(.not. allocated(pw(igrid)%B0)) then
    allocate(pw(igrid)%B0(ixGlo1:ixGhi1,1:ndir,0:ndim))
    allocate(pw(igrid)%J0(ixGlo1:ixGhi1,7-2*ndir:3))
  end if

end subroutine alloc_B0_grid

subroutine dealloc_B0_grid(igrid)
  use mod_global_parameters

  integer, intent(in) :: igrid

  deallocate(pw(igrid)%B0)
  deallocate(pw(igrid)%J0)

end subroutine dealloc_B0_grid
