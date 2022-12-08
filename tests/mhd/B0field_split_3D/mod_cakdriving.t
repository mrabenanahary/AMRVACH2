!#############################################################################
! module amrvacusr - cakdriving
!=============================================================================
! Project : line driven winds
! pre-compile: include this routine in the amrvacusr.t.****** file
! update : 15/07/2009, Allard Jan
! state  : testing
! parameters :
! main variables
!   rho          : density
!   m^C=rho v^C  : momentum
!   e            : total energy
!============================================================================
subroutine cak_pointgrav(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

!
!  Calculate velocity change due to point source gravity
!  and CAK line driving force.
!  As it stands now, the points source has to be in the same
!  'dimensionality' as the grid. E.g. It is not possible to
!  have an xy grid, with the point source localized in the
!  z-direction.
!  On the other hand, it IS possible to have the point source
!  in a location that is not part of the grid. This will in
!  fact often be the case in spherical symmetries with the
!  gravitational source at r=0
!
!  the scaled value for the central mass: eqpar(Mstar_) has to be set in
!  the parfile, together with the coordinates x1ptms,x2ptms,x3ptms
!


include 'amrvacdef.f'


integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(in)    :: wCT(ixI^S,1:nw)
double precision, intent(inout) :: w(ixI^S,1:nw)

double precision                :: gforce(ixG^T)


!-----------------------------------------------------------------------------
if (star_mass == zero) return

call cak_accel(ixI^L,ixO^L,wCT,x,gforce)
!
!  update momentum and energy
!

w(ixO^S,m1_) =  w(ixO^S,m1_) + qdt*gforce(ixO^S)*wCT(ixO^S,rho_)
w(ixO^S,e_)  =  w(ixO^S,e_)  + qdt*gforce(ixO^S)*wCT(ixO^S,m1_)
call phys_check_smallvalues(w,ixI^L,ixO^L,"pointgravcak")

end subroutine cak_pointgrav
!==============================================================================
subroutine cak_pointgrav_getdt(w,ixG^L,ixO^L,dtnew,dx^D,x)

!
! Limits timestep for gravitational pointsource cak force
!
include 'amrvacdef.f'

integer, intent(in)             :: ixG^L, ixO^L
double precision, intent(in)    :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
double precision                :: gforce(ixG^T)

! .. local ..
double precision                :: r2i(ixG^T,1:ndim)
double precision                :: dtgravcak
integer                         :: idims
!-----------------------------------------------------------------------------

if( star_Mass == zero ) return
!
! call cak_accel with w instead of wCT, because you need values for the next time step
!
call cak_accel(ixG^LL,ixO^L,w,x,gforce)
dtnew = min(bigdouble,{^D&minval(sqrt(dabs(block%dx(ixO^S,^D)/max(gforce(ixO^S),smalldouble))))|,})

!dtnew = dtgravcak

end subroutine cak_pointgrav_getdt
!===========================================================================
ubroutine cak_accel(ixI^L,ixO^L,wCT,x,gforce)
include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(in)    :: wCT(ixI^S,1:nw)

double precision, intent(out)   :: gforce(ixG^T)

! ... local ...
double precision  :: opa, oma
double precision,dimension(ixG^T)  :: v1, dv1dr
double precision,dimension(ixG^T)  :: Tlocal,tt, gthin, sigma, &
                                      fdisk, gcak, kqbareff, ptherm
double precision  :: r2i(ixG^T,1:ndim)
integer           :: idims
!-----------------------------------------------------------------------------
call getgvector(ixI^L,ixO^L,x,r2i)
idims=1

call getv(wCT,ixI^L,ixI^L,idims,v1)

select case(typegrad)
case("central")
     call gradient(v1,ixI^L,ixO^L,idims,dv1dr)
case("limited")
     call gradientS(v1,ixI^L,ixO^L,idims,dv1dr)
end select

dv1dr(ixO^S) = max(dv1dr(ixO^S), smalldouble)


call getpthermal(wCT,ixI^L,ixO^L,ptherm)
cak_opa=1.0_dp + cak_alpha
cak_oma=1.0_dp - cak_alpha)

Tlocal(ixO^S)   = ptherm(ixO^S)/wCT(ixO^S,rho_)
tt(ixO^S)= Tlocal(ixO^S)/star_temperature
! exp(-100)==0 here for numerical stability, compiler dependent, fairly arbitrary
! cutoff.
where(tt(ixO^S)<100.0)
 kqbareff(ixO^S)=merge(dexp(-tt(ixO^S)+1.0d0),1.0d0,tt(ixO^S)>1.0_dp)
else where
 kqbareff(ixO^S)=zero
end where
! Eq. 12 Asif
gthin(ixO^S) = (cak_Qbar*cak_kappae)**cak_oma &
              *start_luminosity*r2i(ixO^S,1) &
              /(4.0d0*dpi*cak_oma*const_c**opa)

sigma(ixO^S) = (1.0_dp-v1(ixO^S)*dsqrt(r2i(ixO^S,1))/dv1dr(ixO^S)) &
               * (star_radius**2.0d0) *r2i(ixO^S,1)

!
!  Finite disk correction factor
!


where(sigma(ixO^S) >= 1.0_dp)
      fdisk(ixO^S)= 1.0_dp/cak_opa
else where (sigma(ixO^S) < -1.0d10)
      fdisk(ixO^S) = ((-sigma(ixO^S))**cak_alpha/cak_opa)
else where (dabs(sigma(ixO^S)) > 1.0d-3)
      fdisk(ixO^S) = (1.0_dp-(1.0_dp-sigma(ixO^S))**cak_opa)&
                     /(sigma(ixO^S)*cak_opa)
else where
      fdisk(ixO^S) = 1.0_dp-0.5_dp*cak_alpha*sigma(ixO^S)   &
            * (1.0_dp+(1.0_dp-cak_alpha))*sigma(ixO^S)/3.0d0)
end where

!
!  Calculates cak and gravity at the same time.
!

gforce(ixO^S) = gthin(ixO^S)*fdisk(ixO^S)&
                *(dv1dr(ixO^S)/wCT(ixO^S,rho_))**cak_alpha

gforce(ixO^S) = gforce(ixO^S)*kqbareff(ixO^S) &
                + (eqpar(Newt_)*star_mass*r2i(ixO^S,1))*(star_Eddington-1.0_dp)

return
end subroutine cak_accel

!============================================================================

