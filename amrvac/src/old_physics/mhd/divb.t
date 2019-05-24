{#IFDEF GLM
!=============================================================================
subroutine addsource_glm1(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add divB related sources to w within ixO
! corresponding to Dedner JCP 2002, 175, 645 _equation 24_
! giving the EGLM-MHD scheme
use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(in) :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)
integer :: iw
double precision:: divb(ixI^S)
integer          :: idims
double precision :: gradPsi(ixI^S)
!-----------------------------------------------------------------------------
! We calculate now div B
call getdivb(wCT,ixI^L,ixO^L,divb)

{#IFDEF FCT
! When using FCT, we treat this as a source term
w(ixO^S,psi_) = wCT(ixO^S,psi_) -  cmax_global**2 * qdt * divb(ixO^S) 
}


! Psi = Psi - qdt Ch^2/Cp^2 Psi
if (eqpar(Cr_) < zero) then
  w(ixO^S,psi_) = abs(eqpar(Cr_))*wCT(ixO^S,psi_)
else 
  ! implicit update of psi variable
  w(ixO^S,psi_) = dexp(-qdt*(cmax_global/eqpar(Cr_)))*wCT(ixO^S,psi_)
end if

! gradient of Psi
do idims=1,ndim
   select case(typegrad)
   case("central")
      call gradient(wCT(ixI^S,psi_),ixI^L,ixO^L,idims,gradPsi)
   case("limited")
      call gradientS(wCT(ixI^S,psi_),ixI^L,ixO^L,idims,gradPsi)
   end select
{#IFDEF ENERGY
! e  = e  -qdt (b . grad(Psi))
   w(ixO^S,e_) = w(ixO^S,e_)&
        -qdt*wCT(ixO^S,b0_+idims)*gradPsi(ixO^S)
}
end do

! m = m - qdt b div b
if (B0field) then
{^C&
   w(ixO^S,m^C_)=w(ixO^S,m^C_)-qdt*(wCT(ixO^S,b^C_)+myB0_cell%w(ixO^S,^C)) &
        *divb(ixO^S)\}
else
{^C&
   w(ixO^S,m^C_)=w(ixO^S,m^C_)-qdt*wCT(ixO^S,b^C_)*divb(ixO^S)\}
end if
{#IFDEF ENERGY
! since this option changes energy: smallvalues call
if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"addsource_glm")
}
end subroutine addsource_glm1
!=============================================================================
subroutine addsource_glm2(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add divB related sources to w within ixO
! corresponding to Dedner JCP 2002, 175, 645 _equation 38_
! giving the non conservative EGLM-MHD scheme.
use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(in) :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)
integer :: iw
double precision:: divb(ixI^S)
integer          :: idims
double precision :: gradPsi(ixI^S)
!-----------------------------------------------------------------------------
! We calculate now div B
call getdivb(wCT,ixI^L,ixO^L,divb)

{#IFDEF FCT
! When using FCT, we treat this as a source term
! We calculate now div B consistent with FCT
w(ixO^S,psi_) = wCT(ixO^S,psi_) -  cmax_global**2 * qdt * divb(ixO^S) 
}


! Psi = Psi - qdt Ch^2/Cp^2 Psi
if (eqpar(Cr_) < zero) then
  w(ixO^S,psi_) = abs(eqpar(Cr_))*wCT(ixO^S,psi_)
else 
  ! implicit update of psi variable
  w(ixO^S,psi_) = dexp(-qdt*(cmax_global/eqpar(Cr_)))*wCT(ixO^S,psi_)
end if

! gradient of Psi
do idims=1,ndim
   select case(typegrad)
   case("central")
      call gradient(wCT(ixI^S,psi_),ixI^L,ixO^L,idims,gradPsi)
   case("limited")
      call gradientS(wCT(ixI^S,psi_),ixI^L,ixO^L,idims,gradPsi)
   end select


! Psi=Psi - qdt (v . grad(Psi))
   w(ixO^S,psi_) = w(ixO^S,psi_)&
                   -qdt*wCT(ixO^S,m0_+idims)/wCT(ixO^S,rho_)*gradPsi(ixO^S)

{#IFDEF ENERGY
! e  = e  - qdt (b . grad(Psi))
   w(ixO^S,e_) = w(ixO^S,e_)&
        -qdt*wCT(ixO^S,b0_+idims)*gradPsi(ixO^S)
}
end do

{#IFDEF ENERGY
! e = e - qdt (v . b) * div b
   w(ixO^S,e_)=w(ixO^S,e_)-&
      qdt*(^C&wCT(ixO^S,m^C_)*wCT(ixO^S,b^C_)+)/wCT(ixO^S,rho_)*divb(ixO^S)
}
! b = b - qdt v * div b
{^C&w(ixO^S,b^C_)=w(ixO^S,b^C_)-qdt*wCT(ixO^S,m^C_)/wCT(ixO^S,rho_)*divb(ixO^S)\}


! m = m - qdt b * div b
if (B0field) then
{^C&
   w(ixO^S,m^C_)=w(ixO^S,m^C_)-qdt*(wCT(ixO^S,b^C_)+myB0_cell%w(ixO^S,^C)) &
        *divb(ixO^S)\}
else
{^C&
   w(ixO^S,m^C_)=w(ixO^S,m^C_)-qdt*wCT(ixO^S,b^C_)*divb(ixO^S)\}
end if

{#IFDEF ENERGY
! since this option changes energy: smallvalues call
if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"addsource_glm2")
}
end subroutine addsource_glm2
!=============================================================================
subroutine addsource_glm3(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add divB related sources to w within ixO
! corresponding to Dedner JCP 2002, 175, 645 _equation (1a), (1b), (4), (1d), 19
! conservative hyperbolic mixed GLM-MHD with no additional source terms.
use mod_global_parameters

integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(in)    :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)
{#IFDEF FCT
double precision                :: divb(ixI^S)
}
!-----------------------------------------------------------------------------

{#IFDEF FCT
! When using FCT, we treat this as a source term
! We calculate now div B consistent with FCT
call getdivb(wCT,ixI^L,ixO^L,divb)
w(ixO^S,psi_) = wCT(ixO^S,psi_) -  cmax_global**2 * qdt * divb(ixO^S) 
}


! Psi = Psi - qdt Ch^2/Cp^2 Psi
if (eqpar(Cr_) < zero) then
  w(ixO^S,psi_) = abs(eqpar(Cr_))*w(ixO^S,psi_)
else 
  ! implicit update of psi variable
  w(ixO^S,psi_) = dexp(-qdt*(cmax_global/eqpar(Cr_)))*w(ixO^S,psi_)
end if

end subroutine addsource_glm3
}
!=============================================================================
subroutine addsource_powel(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add divB related sources to w within ixO
! corresponding to Powel
use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(in) :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)
integer :: iw
double precision:: divb(ixI^S)
!-----------------------------------------------------------------------------
! We calculate now div B
call getdivb(wCT,ixI^L,ixO^L,divb)
{#IFDEF ENERGY
! e = e - qdt (v . b) * div b
w(ixO^S,e_)=w(ixO^S,e_)-&
     qdt*(^C&wCT(ixO^S,m^C_)*wCT(ixO^S,b^C_)+)/wCT(ixO^S,rho_)*divb(ixO^S)
}

! b = b - qdt v * div b
{^C&w(ixO^S,b^C_)=w(ixO^S,b^C_)-qdt*wCT(ixO^S,m^C_)/wCT(ixO^S,rho_)*divb(ixO^S)\}

! m = m - qdt b * div b
if (B0field) then
{^C&
   w(ixO^S,m^C_)=w(ixO^S,m^C_)-qdt*(wCT(ixO^S,b^C_)+myB0_cell%w(ixO^S,^C)) &
        *divb(ixO^S)\}
else
{^C&
   w(ixO^S,m^C_)=w(ixO^S,m^C_)-qdt*wCT(ixO^S,b^C_)*divb(ixO^S)\}
end if

! since this option changes energy: smallvalues call
if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"addsource_powel")

end subroutine addsource_powel
!=============================================================================
subroutine addsource_janhunen(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add divB related sources to w within ixO
! corresponding to Janhunen, just the term in the induction equation.
use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(in) :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)
!.. local ..
double precision:: divb(ixI^S)
!-----------------------------------------------------------------------------
! We calculate now div B
call getdivb(wCT,ixI^L,ixO^L,divb)

! b = b - qdt v * div b
{^C&w(ixO^S,b^C_)=w(ixO^S,b^C_)-qdt*wCT(ixO^S,m^C_)/wCT(ixO^S,rho_)*divb(ixO^S)\}
end subroutine addsource_janhunen
!=============================================================================
subroutine addsource_linde(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add Linde's divB related sources to wnew within ixO
use mod_global_parameters

integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(in)    :: dx^D
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
integer :: iw, idims, ix^L, ixp^L, i^D, iside
double precision :: divb(ixI^S),graddivb(ixI^S)
!-----------------------------------------------------------------------------

! Calculate div B
ix^L=ixO^L^LADD1;
call getdivb(wCT,ixI^L,ix^L,divb)
! for AMR stability, retreat one cell layer from the boarders of level jump
ixp^L=ixO^L;
do idims=1,ndim
  select case(idims)
   {case(^D)
      do iside=1,2
        i^DD=kr(^DD,^D)*(2*iside-3);
        if(leveljump(i^DD)) then
          if(iside==1) then
            ixpmin^D=ixOmin^D-i^D
          else
            ixpmax^D=ixOmax^D-i^D
          end if 
        end if
      end do
   \}
  end select
end do

! Add Linde's diffusive terms
do idims=1,ndim
   ! Calculate grad_idim(divb)
   select case(typegrad)
   case("central")
     call gradient(divb,ixI^L,ixp^L,idims,graddivb)
   case("limited")
     call gradientS(divb,ixI^L,ixp^L,idims,graddivb)
   end select
   !ixmin^D=ixpmin^D;
   !ixmax^D=merge(ixmin^D,ixmax^D,kr(idims,^D)==1);
   !call gradient(divb,ixI^L,ix^L,idims,graddivb)
   !ixmin^D=merge(ixmax^D,ixmin^D,kr(idims,^D)==1);
   !ixmax^D=ixpmax^D;
   !call gradient(divb,ixI^L,ix^L,idims,graddivb)
   !ixmin^D=ixpmin^D+kr(idims,^D);
   !ixmax^D=ixpmax^D-kr(idims,^D);
   !call gradientS(divb,ixI^L,ix^L,idims,graddivb)

   ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
   if (slab) then
      graddivb(ixp^S)=graddivb(ixp^S)*divbdiff/(^D&1.0d0/dxlevel(^D)**2+)
   else
      graddivb(ixp^S)=graddivb(ixp^S)*divbdiff &
                      /(^D&1.0d0/mygeo%dx(ixp^S,^D)**2+)
   end if
   do iw= iw^LIM
      if (iw==b0_+idims) then
         ! B_idim += eta*grad_idim(divb)
         w(ixp^S,iw)=w(ixp^S,iw)+graddivb(ixp^S)
{#IFDEF ENERGY
      else if (iw==e_ .and. typedivbdiff=='all') then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixp^S,iw)=w(ixp^S,iw)+wCT(ixp^S,b0_+idims)*graddivb(ixp^S)
}
      end if
   end do
end do

if(fixsmall) call smallvalues(w,x,ixI^L,ixp^L,"addsource_linde")

end subroutine addsource_linde
!=============================================================================
subroutine getdivb(w,ixI^L,ixO^L,divb)

! Calculate div B within ixO

use mod_global_parameters

integer, intent(in)                :: ixI^L, ixO^L
double precision, intent(in)       :: w(ixI^S,1:nw)
double precision                   :: divb(ixI^S)

double precision                   :: bvec(ixI^S,1:ndir)

{#IFDEF FCT
integer                            :: ixC^L, idir, idim, ixJp^L, ic^D, ix^L
integer                            :: ixKp^L, ixJpKp^L, ixJm^L, ixJmKm^L
double precision                   :: divb_corner(ixI^S), sign
}
!-----------------------------------------------------------------------------

bvec(ixI^S,1:ndir)=w(ixI^S,b0_+1:b0_+ndir)

{#IFNDEF FCT
select case(typediv)
case("central")
  call divvector(bvec,ixI^L,ixO^L,divb)
case("limited")
  call divvectorS(bvec,ixI^L,ixO^L,divb)
end select
}{#IFDEF FCT

! For fct, we calculate the divB on the corners according to Toth (2000), 
! eq. (27) and average to the cell centers for output.

{ixCmax^D=ixOmax^D;}
{ixCmin^D=ixOmin^D-1;} ! Extend range by one


! Get the corner-centered divb:

divb_corner(ixC^S) = zero
do idir = 1, ndim ! idir is the component of the field to consider (j)
{do ic^DB=0,1\}
   {ix^L=ixC^L+ic^D;}

   select case(idir)
      {^D&   
           case(^D)
      sign = dble(ic^D*2 - 1)
      \}
   end select

   divb_corner(ixC^S) = divb_corner(ixC^S) &
        + sign * bvec(ix^S,idir)/dxlevel(idir)
   
   {end do\}
end do
divb_corner(ixC^S) = divb_corner(ixC^S) / 2.0d0**(ndim-1)



! Now average back to the cell centers:

divb(ixO^S) = zero
{do ic^DB=-1,0\}
   {ixC^L=ixO^L+ic^D;}

   divb(ixO^S) = divb(ixO^S) &
        + divb_corner(ixC^S)

{end do\}
divb(ixO^S) = divb(ixO^S) / 2.0d0**(ndim)

}
end subroutine getdivb
!=============================================================================
