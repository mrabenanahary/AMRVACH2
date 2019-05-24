

subroutine set_B0_grid(igrid,qt)
  use mod_global_parameters

  integer, intent(in)             :: igrid
  double precision, intent(in)    :: qt

  call set_B0_cell(pw(igrid)%B0(:^D&,:,0),pw(igrid)%x,ixG^LL,ixG^LL,qt)
  call set_J0_cell(igrid,pw(igrid)%J0,ixG^LL,ixM^LL^LADD1,qt)
  call set_B0_face(igrid,pw(igrid)%x,ixG^LL,ixM^LL,qt)

end subroutine set_B0_grid

subroutine set_B0_cell(wB0,x,ixI^L,ix^L,qt)
  use mod_constants
  use mod_usr_methods, only: usr_set_B0
  use mod_global_parameters
  use mod_geometry
  implicit none
  integer, intent(in)             :: ixI^L,ix^L
  double precision, intent(inout) :: wB0(ixI^S,1:ndir)
  double precision, intent(in)    :: x(ixI^S,1:ndim)
  double precision, intent(in)    :: qt
  integer                         :: idims
  real(dp)                        :: mu0_sphe(ixI^S,1:ndim)
  real(dp)                        :: mu0_cart(1:ndim)
  real(dp)                        :: dmutheta,dmuphi
  wB0(ix^S,1:ndir)=zero

  ! approximate cell-averaged B0 as cell-centered B0
  select case (typeaxial)
  case ("spherical")
     {^NOONED
     if (dabs(Bdip)>smalldouble) then
        wB0(ix^S,1)=2.0d0*Bdip*dcos(x(ix^S,2))/x(ix^S,1)**3
        wB0(ix^S,2)=Bdip*dsin(x(ix^S,2))/x(ix^S,1)**3
     end if
  
     if (abs(Bquad)>smalldouble) then
        wB0(ix^S,1)=wB0(ix^S,1) &
             +Bquad*0.5d0*(1.0d0+3.0d0*dcos(2.0d0*x(ix^S,2)))/x(ix^S,1)**4
        wB0(ix^S,2)=wB0(ix^S,2)+Bquad*dsin(2.0d0*x(ix^S,2))/x(ix^S,1)**4
     end if
     if (abs(Boct)>smalldouble) then
        wB0(ix^S,1)=wB0(ix^S,1) &
                     +Boct*(10.0d0*dcos(2.0d0*x(ix^S,2))-2.0d0) &
                          *dcos(x(ix^S,2))/x(ix^S,1)**5
        wB0(ix^S,2)=wB0(ix^S,2) &
                     +Boct*1.5d0*(3.0d0+5.0d0*dcos(2.0d0*x(ix^S,2))) &
                          *dsin(x(ix^S,2))/x(ix^S,1)**5
     end if
     if(dabs(mu0dip)>smalldouble)then
     {^IFTHREED
       if(dabs(mutheta_period)<smalldouble)then
        dmutheta = 0.0_dp
       else
        dmutheta = 2.0*dpi*qt/mutheta_period
       end if
       if(dabs(muphi_period)<smalldouble)then
        dmuphi = 0.0_dp
       else
        dmuphi = 2.0*dpi*qt/muphi_period
       end if

       mu0_cart(z_)=mu0dip*dcos(mu0theta+dmutheta)
       mu0_cart(x_)=mu0dip*dsin(mu0theta+dmutheta)&
                *dcos(mu0phi+dmuphi)
       mu0_cart(y_)=mu0dip*dsin(mu0theta+dmutheta)&
                *dsin(mu0phi+dmuphi)
       call transV_cartesian_to_spherical(ixI^L,ix^L,x,mu0_cart,mu0_sphe)
       wB0(ix^S,r_) = 2.0_dp*mu0_sphe(ix^S,r_)/x(ix^S,r_)**3.0_dp
       Loop_idims_mu : do idims=2,ndim
        wB0(ix^S,idims) = -mu0_sphe(ix^S,idims)/x(ix^S,r_)**3.0_dp
       end do Loop_idims_mu
      }
     end if
     }
  end select
  
  if (associated(usr_set_B0)) call usr_set_B0(ixI^L,ix^L,x,wB0,qt)

end subroutine set_B0_cell

subroutine set_J0_cell(igrid,wJ0,ixI^L,ix^L,qt)
  use mod_usr_methods, only: usr_set_J0
  use mod_global_parameters
  use mod_geometry

  integer, intent(in)             :: igrid,ixI^L,ix^L
  double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:3)
  double precision, intent(in)    :: qt
  integer :: idirmin0, idirmin

  if(associated(usr_set_J0)) then
    call usr_set_J0(ixI^L,ix^L,pw(igrid)%x,wJ0,qt)
  else
    idirmin0 = 7-2*ndir
    call curlvector(pw(igrid)%B0(:^D&,:,0),ixI^L,ix^L,wJ0,idirmin,idirmin0,ndir)
  end if

end subroutine set_J0_cell

subroutine set_B0_face(igrid,x,ixI^L,ix^L,qt)
  use mod_global_parameters

  integer, intent(in) :: igrid, ixI^L, ix^L
  double precision, intent(in) :: x(ixI^S,1:ndim)
  double precision, intent(in)    :: qt

  double precision :: delx(ixI^S,1:ndim)
  double precision :: xC(ixI^S,1:ndim),xshift^D
  integer :: idims, ixC^L, ix, idims2

  if(slab)then
   ^D&delx(ixI^S,^D)=rnode(rpdx^D_,igrid)\
  else
   ! for all non-cartesian and stretched coordinate(s)
   delx(ixI^S,1:ndim)=pw(igrid)%dx(ixI^S,1:ndim)
  endif

  do idims=1,ndim
     ixCmin^D=ixmin^D-kr(^D,idims); ixCmax^D=ixmax^D;
     ! always xshift=0 or 1/2
     xshift^D=half*(one-kr(^D,idims));
     do idims2=1,ndim
       select case(idims2)
       {case(^D)
         do ix = ixC^LIM^D
           ! xshift=half: this is the cell center coordinate
           ! xshift=0: this is the cell edge i+1/2 coordinate
           xC(ix^D%ixC^S,^D)=x(ix^D%ixC^S,^D)+(half-xshift^D)*delx(ix^D%ixC^S,^D)
         end do\}
       end select
     end do
     call set_B0_cell(pw(igrid)%B0(:^D&,:,idims),xC,ixI^L,ixC^L,qt)
  end do

end subroutine set_B0_face

subroutine alloc_B0_grid(igrid)
  use mod_global_parameters

  integer, intent(in) :: igrid

  if(.not. allocated(pw(igrid)%B0)) then
    allocate(pw(igrid)%B0(ixG^T,1:ndir,0:ndim))
    allocate(pw(igrid)%J0(ixG^T,7-2*ndir:3))
  end if

end subroutine alloc_B0_grid

subroutine dealloc_B0_grid(igrid)
  use mod_global_parameters

  integer, intent(in) :: igrid

  deallocate(pw(igrid)%B0)
  deallocate(pw(igrid)%J0)

end subroutine dealloc_B0_grid
