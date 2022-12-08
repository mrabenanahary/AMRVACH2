!##############################################################################
! module amrvacphys/- srhd  -- version april 2009
! added positions to fluxes etc. -- march 2012
!=============================================================================

INCLUDE:amrvacnul/addsource.t
INCLUDE:amrvacnul/getdt.t
INCLUDE:amrvacphys/rhdroe.t
INCLUDE:amrvacphys/srhdhllc.t
INCLUDE:amrvacphys/correctauxsrhd.t
!=============================================================================
subroutine checkglobaldata
use mod_global_parameters
!-----------------------------------------------------------------------------
minp    = max(zero,smallp)
minrho  = max(zero,smallrho)
smalltau= minp/(eqpar(gamma_)-one)
smallxi = minrho + minp*eqpar(gamma_)/(eqpar(gamma_)-one)

end subroutine checkglobaldata
!=============================================================================
subroutine initglobaldata

! set default values 

use mod_global_parameters

!-----------------------------------------------------------------------------
eqpar(gamma_)=4./3.

end subroutine initglobaldata
!=============================================================================
subroutine checkw(checkprimitive,ixI^L,ixO^L,w,flag)

use mod_global_parameters
  
logical, intent(in)          :: checkprimitive
integer, intent(in)          :: ixI^L, ixO^L
double precision, intent(in) :: w(ixI^S,nw)
logical, intent(out)         :: flag(ixG^T)
!-----------------------------------------------------------------------------

flag(ixG^T)=.true.

if (checkprimitive) then
     if(useprimitiveRel)then
     ! check   rho>=0, p>=smallp
     flag(ixO^S) = (w(ixO^S,rho_) >= minrho).and. &
                   (w(ixO^S,pp_)  >= minp)
     else
     ! check  v^2 < 1, rho>=0, p>=smallp
     flag(ixO^S) = ({^C&w(ixO^S,v^C_)**2.0d0+} < one).and. &
                   (w(ixO^S,rho_) >= minrho).and. &
                   (w(ixO^S,pp_)  >= minp)
     endif
else
     ! Check D>=0 and lower limit for tau
     flag(ixO^S) = (w(ixO^S,d_)    >= minrho).and. &
                   (w(ixO^S,tau_)  >= smalltau)
endif

end subroutine checkw
!=============================================================================
subroutine conserve(ixI^L,ixO^L,w,x,patchw)

! Transform primitive variables into conservative ones
! (rho,v,p) ---> (D,S,tau)
! call to smallvalues
! --> latter only used for correcting procedure in correctaux 
! --> input array patchw for spatially selective transformation

use mod_global_parameters

integer, intent(in)               :: ixI^L, ixO^L
double precision, intent(inout)   :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
logical, intent(in)               :: patchw(ixG^T)

double precision,dimension(ixG^T) :: xi
!-----------------------------------------------------------------------------

if(useprimitiveRel)then
  ! assumes four velocity computed in primitive (rho u p) with u=lfac*v
  xi(ixO^S)=one+({^C&w(ixO^S,u^C_)**2.0d0+})
  ! fill the auxiliary variable lfac_ (Lorentz Factor) and p_ (pressure)
  w(ixO^S,lfac_)=dsqrt(xi(ixO^S))
  w(ixO^S,p_)=w(ixO^S,pp_)
else
  ! assumes velocity in primitive (rho v p) 
  xi(ixO^S)=one-({^C&w(ixO^S,v^C_)**2.0d0+})
  ! fill the auxiliary variable lfac_ (Lorentz Factor) and p_ (pressure)
  w(ixO^S,lfac_)=one/dsqrt(xi(ixO^S))
  w(ixO^S,p_)=w(ixO^S,pp_)
endif

if(useprimitiveRel)then
  ! compute xi=Lfac w  (enthalphy w)
  xi(ixO^S)=w(ixO^S,lfac_)* &
        (w(ixO^S,rho_)+w(ixO^S,p_)*eqpar(gamma_)/(eqpar(gamma_)-one))
else
  ! compute xi=Lfac^2 w  (enthalphy w)
  xi(ixO^S)=w(ixO^S,lfac_)*w(ixO^S,lfac_)* &
        (w(ixO^S,rho_)+w(ixO^S,p_)*eqpar(gamma_)/(eqpar(gamma_)-one))
endif

if(useprimitiveRel)then
  w(ixO^S,d_)=w(ixO^S,rho_)*w(ixO^S,lfac_) 
  ^C&w(ixO^S,s^C_)=xi(ixO^S)*w(ixO^S,u^C_);
  w(ixO^S,tau_)=xi(ixO^S)*w(ixO^S,lfac_) - w(ixO^S,p_) - w(ixO^S,d_)
else
  w(ixO^S,d_)=w(ixO^S,rho_)*w(ixO^S,lfac_) 
  ^C&w(ixO^S,s^C_)=xi(ixO^S)*w(ixO^S,v^C_);
  w(ixO^S,tau_)=xi(ixO^S) - w(ixO^S,p_) - w(ixO^S,d_)
endif

{#IFDEF TRACER
! We got D, now we can get the conserved tracers:
{^FL&w(ixO^S,tr^FL_) = w(ixO^S,d_)*w(ixO^S,tr^FL_)\}
}

if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"conserve")

end subroutine conserve
!=============================================================================
subroutine conserven(ixI^L,ixO^L,w,patchw)

! Transform primitive variables into conservative ones
! (rho,v,p) ---> (D,S,tau)
! no call to smallvalues
! --> latter only used for correcting procedure in correctaux 
! --> input array patchw for spatially selective transformation

use mod_global_parameters

integer, intent(in)               :: ixI^L, ixO^L
double precision, intent(inout)   :: w(ixI^S,nw)
logical, intent(in)               :: patchw(ixG^T)

double precision,dimension(ixG^T) :: xi
!-----------------------------------------------------------------------------

if(useprimitiveRel)then
  ! assumes four velocity computed in primitive (rho u p) with u=lfac*v
  where(.not.patchw(ixO^S))
     xi(ixO^S)=one+({^C&w(ixO^S,u^C_)**2.0d0+})
     ! fill the auxiliary variable lfac_ (Lorentz Factor) and p_ (pressure)
     w(ixO^S,lfac_)=dsqrt(xi(ixO^S))
     w(ixO^S,p_)=w(ixO^S,pp_)
  endwhere
else
  ! assumes velocity in primitive (rho v p) 
  where(.not.patchw(ixO^S))
     xi(ixO^S)=one-({^C&w(ixO^S,v^C_)**2.0d0+})
     ! fill the auxiliary variable lfac_ (Lorentz Factor) and p_ (pressure)
     w(ixO^S,lfac_)=one/dsqrt(xi(ixO^S))
     w(ixO^S,p_)=w(ixO^S,pp_)
  endwhere
endif

if(useprimitiveRel)then
   ! compute xi=Lfac w  (enthalphy w)
   where(.not.patchw(ixO^S))
      xi(ixO^S)=w(ixO^S,lfac_)* &
            (w(ixO^S,rho_)+w(ixO^S,p_)*eqpar(gamma_)/(eqpar(gamma_)-one))
   endwhere
else
   ! compute xi=Lfac^2 w  (enthalphy w)
   where(.not.patchw(ixO^S))
      xi(ixO^S)=w(ixO^S,lfac_)*w(ixO^S,lfac_)* &
            (w(ixO^S,rho_)+w(ixO^S,p_)*eqpar(gamma_)/(eqpar(gamma_)-one))
   endwhere
endif

if(useprimitiveRel)then
  where(.not.patchw(ixO^S))
    w(ixO^S,d_)=w(ixO^S,rho_)*w(ixO^S,lfac_) 
    ^C&w(ixO^S,s^C_)=xi(ixO^S)*w(ixO^S,u^C_);
    w(ixO^S,tau_)=xi(ixO^S)*w(ixO^S,lfac_) - w(ixO^S,p_) - w(ixO^S,d_)
  endwhere
else
  where(.not.patchw(ixO^S))
    w(ixO^S,d_)=w(ixO^S,rho_)*w(ixO^S,lfac_) 
    ^C&w(ixO^S,s^C_)=xi(ixO^S)*w(ixO^S,v^C_);
    w(ixO^S,tau_)=xi(ixO^S) - w(ixO^S,p_) - w(ixO^S,d_)
  endwhere
endif

{#IFDEF TRACER
! We got D, now we can get the conserved tracers:
where(.not.patchw(ixO^S))
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,d_)*w(ixO^S,tr^FL_)\}
end where
}

end subroutine conserven
!=============================================================================
subroutine primitive(ixI^L,ixO^L,w,x)

! Transform conservative variables into primitive ones
! (D,S,tau) ---> (rho,v,p)

use mod_global_parameters

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(inout) :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

! Calculate pressure and Lorentz factor from conservative variables only
! these are put in lfac_ and p_ auxiliaries
call getaux(.true.,w,x,ixI^L,ixO^L,'primitive')
! note: on exit from getaux: gauranteed 
!    xi=(d+tau+p)>smallp*gamma/(gamma-1), lfac>=1, p>smallp

! replace conservative with primitive variables
! compute velocity
if(useprimitiveRel)then
^C&w(ixO^S,u^C_)=w(ixO^S,lfac_)*w(ixO^S,s^C_)/(w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_));
else
^C&w(ixO^S,v^C_)=w(ixO^S,s^C_)/(w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_));
endif
! compute density
w(ixO^S,rho_)=w(ixO^S,d_)/w(ixO^S,lfac_)
! fill pressure
w(ixO^S,pp_)=w(ixO^S,p_)

{#IFDEF TRACER
! We got lor, rho, Dtr, now we can get the tracers:
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,Dtr^FL_)/w(ixO^S,lfac_)/w(ixO^S,rho_)\}
}

end subroutine primitive
!=============================================================================
subroutine primitiven(ixI^L,ixO^L,w,patchw)

! Transform conservative variables into primitive ones
! (D,S,tau) ---> (rho,v,p)
! similar to subroutine primitive, but no call to getaux here!!
! --> only used for correcting procedure in correctaux 
! --> input array patchw for spatially selective transformation

use mod_global_parameters

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(inout) :: w(ixI^S,nw)
logical, intent(in),dimension(ixG^T)   :: patchw
!-----------------------------------------------------------------------------

if(useprimitiveRel)then
 where(.not.patchw(ixO^S))
    ^C&w(ixO^S,u^C_)=w(ixO^S,lfac_)*w(ixO^S,s^C_)/(w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_));
 end where
else
 where(.not.patchw(ixO^S))
    ^C&w(ixO^S,v^C_)=w(ixO^S,s^C_)/(w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_));
 end where
endif

where(.not.patchw(ixO^S))
  ! compute density
  w(ixO^S,rho_)=w(ixO^S,d_)/w(ixO^S,lfac_)
  ! fill pressure
  w(ixO^S,pp_)=w(ixO^S,p_)
end where

{#IFDEF TRACER
! We got lor, rho, Dtr, now we can get the tracers:
where(.not.patchw(ixO^S))
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,Dtr^FL_)/w(ixO^S,lfac_)/w(ixO^S,rho_)\}
endwhere
}
end subroutine primitiven
!=============================================================================
subroutine e_to_rhos(ixI^L,ixO^L,w,x)

! replace tau by DS

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

! Calculate pressure and Lorentz factor from conservative variables only
! these are put in lfac_ and p_ auxiliaries
call getaux(.true.,w,x,ixI^L,ixO^L,'e_to_rhos')

w(ixO^S,rhos_)=w(ixO^S,d_)*w(ixO^S,p_)*(w(ixO^S,lfac_)/w(ixO^S,d_))**eqpar(gamma_)

end subroutine e_to_rhos
!=============================================================================
subroutine rhos_to_e(ixI^L,ixO^L,w,x)

! replace DS by tau

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision :: w(ixI^S,nw), xi(ixG^T)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

! Calculate pressure and Lorentz factor from conservative variables BUT
! with tau replaced by DS, these are put in lfac_ and p_ auxiliaries
call getaux2(.true.,w,x,ixI^L,ixO^L,'rhos_to_e')

xi(ixO^S)=w(ixO^S,lfac_)*w(ixO^S,d_)+w(ixO^S,lfac_)**2*w(ixO^S,p_)*eqpar(gamma_)/(eqpar(gamma_)-one)
w(ixO^S,tau_)=xi(ixO^S)-w(ixO^S,p_)-w(ixO^S,d_)

end subroutine rhos_to_e
!=============================================================================
subroutine ppmflatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dp)

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L,ixL^L,ixR^L
double precision, intent(in)  :: w(ixI^S,nw),d2w(ixG^T,1:nwflux)

double precision, intent(out) :: drho(ixG^T),dp(ixG^T)
!-----------------------------------------------------------------------------

if(useprimitive)then
 drho(ixO^S) = eqpar(gamma_)*dabs(d2w(ixO^S,rho_))&
               /min(w(ixL^S,rho_),w(ixR^S,rho_))
 dp(ixO^S)   = dabs(d2w(ixO^S,pp_))/min(w(ixL^S,pp_),w(ixR^S,pp_))
end if

end subroutine ppmflatcd
!=============================================================================
subroutine ppmflatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp,dv)

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L
integer, intent(in)           :: idims
double precision, intent(in)  :: w(ixI^S,nw)

double precision, intent(out) :: drho(ixG^T),dp(ixG^T),dv(ixG^T)

double precision :: v(ixG^T)
!-----------------------------------------------------------------------------

if(useprimitive)then
   ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
   where (dabs(w(ixRR^S,pp_)-w(ixLL^S,pp_))>smalldouble)
      drho(ixO^S) = dabs((w(ixR^S,pp_)-w(ixL^S,pp_))&
                   /(w(ixRR^S,pp_)-w(ixLL^S,pp_)))
   elsewhere
      drho(ixO^S) = zero
   end where

   !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26 
   !  use "dp" to save sound speed squared
   dp(ixO^S)=(eqpar(gamma_)*w(ixO^S,pp_) &
         /(w(ixO^S,rho_)+eqpar(gamma_)*w(ixO^S,pp_)/(eqpar(gamma_)-one)))

   dp(ixO^S) = dabs(w(ixR^S,pp_)-w(ixL^S,pp_))&
                /(w(ixO^S,rho_)*dp(ixO^S))
   if (useprimitiveRel) then
     v(ixI^S)  = w(ixI^S,u0_+idims)/w(ixI^S,lfac_)
   else
     v(ixI^S)  = w(ixI^S,v0_+idims)
   end if
   call gradient(v,ixI^L,ixO^L,idims,dv)
end if

end subroutine ppmflatsh
!=============================================================================
subroutine getv(w,x,ixI^L,ixO^L,idim,v)

! Calculate v_idim=m_idim/rho within ixO^L

use mod_global_parameters

integer, intent(in)           :: ixI^L, ixO^L, idim
double precision, intent(in)  :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(out) :: v(ixG^T)
!-----------------------------------------------------------------------------
! assuming getv FOLLOWS a getaux call for updated p_ entry and
! well-behaved xi=d+tau+p

! first case only happens when d=0 and p/tau are at enforced minimal
! values, namely p=smallp and tau=smallp/(gamma-1) (no velocities)
! This will typically NEVER be the case, but still...
where(w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_)==smallxi)
   v(ixO^S)=zero
elsewhere
   v(ixO^S)=w(ixO^S,s0_+idim)/ &
        (w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_))
endwhere

end subroutine getv
!=============================================================================
subroutine getcmax(new_cmax,w,x,ixI^L,ixO^L,idim,cmax,cmin,needcmin)

! Calculate cmax_idim=csound+abs(v_idim) within ixO^L

use mod_global_parameters

logical, intent(in)                               :: new_cmax,needcmin
integer, intent(in)                               :: ixI^L, ixO^L, idim
double precision, dimension(ixI^S,nw), intent(in) :: w
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, dimension(ixG^T), intent(out)   :: cmax,cmin

double precision, dimension(ixG^T)                :: csound2,rhoh,vidim2,v2,vidim
!-----------------------------------------------------------------------------
rhoh(ixO^S) = w(ixO^S,d_)/w(ixO^S,lfac_) + &
         eqpar(gamma_)*w(ixO^S,p_)/(eqpar(gamma_)-one)
csound2(ixO^S)=eqpar(gamma_)*w(ixO^S,p_)/rhoh(ixO^S)
if(.not.needcmin)then
   v2(ixO^S)=({^C&w(ixO^S,s^C_)**2.0d0+})/ &
             ((rhoh(ixO^S)*w(ixO^S,lfac_)*w(ixO^S,lfac_))**2.0d0)
else
   vidim(ixO^S)=(w(ixO^S,s0_+idim)/ &
                  (rhoh(ixO^S)*w(ixO^S,lfac_)*w(ixO^S,lfac_)))
endif

if(.not.needcmin)then
  if (ndir==1) then
     cmax(ixO^S)= (dsqrt(v2(ixO^S))+dsqrt(csound2(ixO^S)))/ & 
                  (one+dsqrt(csound2(ixO^S)*v2(ixO^S)))
  else
     vidim2(ixO^S)=(w(ixO^S,s0_+idim)/ &
                  (rhoh(ixO^S)*w(ixO^S,lfac_)*w(ixO^S,lfac_)))**2.0d0

     cmax(ixO^S)=( dsqrt(vidim2(ixO^S))*(one-csound2(ixO^S)) + &
                dsqrt(csound2(ixO^S)*(one-v2(ixO^S))*( &
        one-v2(ixO^S)*csound2(ixO^S)-vidim2(ixO^S)*(one-csound2(ixO^S))) &
               ) ) / (one-v2(ixO^S)*csound2(ixO^S))
  endif
else
  if (ndir==1) then
     cmax(ixO^S)= min(one,max(zero,(vidim(ixO^S)+dsqrt(csound2(ixO^S)))/ & 
                  (one+dsqrt(csound2(ixO^S))*vidim(ixO^S))))
     cmin(ixO^S)= max(-one,min(zero,(vidim(ixO^S)-dsqrt(csound2(ixO^S)))/ & 
                  (one-dsqrt(csound2(ixO^S))*vidim(ixO^S))))
  else
     v2(ixO^S)=({^C&w(ixO^S,s^C_)**2.0d0+})/ &
             (rhoh(ixO^S)*w(ixO^S,lfac_)*w(ixO^S,lfac_))**2.0d0

     cmax(ixO^S)=min(one,max(zero,(vidim(ixO^S)*(one-csound2(ixO^S)) + &
      dsqrt(csound2(ixO^S)*(one-v2(ixO^S))*( &
      one-v2(ixO^S)*csound2(ixO^S)-vidim(ixO^S)**2.0d0*(one-csound2(ixO^S))) &
       ) ) / (one-v2(ixO^S)*csound2(ixO^S))))

     cmin(ixO^S)=max(-one,min(zero,(vidim(ixO^S)*(one-csound2(ixO^S)) - &
       dsqrt(csound2(ixO^S)*(one-v2(ixO^S))*( &
       one-v2(ixO^S)*csound2(ixO^S)-vidim(ixO^S)**2.0d0*(one-csound2(ixO^S))) &
                 ) ) / (one-v2(ixO^S)*csound2(ixO^S))))
  endif
endif

end subroutine getcmax
!=============================================================================
subroutine getflux(w,x,ixI^L,ixO^L,iw,idim,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.

use mod_global_parameters

integer, intent(in)           :: ixI^L, ixO^L, iw, idim
double precision, intent(in)  :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(out) ::  f(ixG^T)
logical, intent(out)          :: transport
!-----------------------------------------------------------------------------
! assuming getflux FOLLOWS a getaux call for updated p_ entry and
! well-behaved xi=d+tau+p, write:

if(iw==s0_+idim)then
     ! f_i[s_i]=v_i*s_i + p
     f(ixO^S)=w(ixO^S,p_) 
{#IFDEF TRACER
{else if (iw==tr^FL_) then 
      f(ixO^S)=zero\}
}
else if(iw==tau_)then
     ! f_i[tau]=v_i*tau + v_i*p
     ! first case only happens when d=0 and p/tau are at enforced minimal
     ! values, namely p=smallp and tau=smallp/(gamma-1) (no velocities)
     ! This will typically NEVER be the case, but still...
     where(w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_)==smallxi)
        f(ixO^S)=zero
     elsewhere
        f(ixO^S)=w(ixO^S,s0_+idim)*w(ixO^S,p_)/ &
             (w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_))
     endwhere
else
     f(ixO^S)=zero
endif

transport=.true.

end subroutine getflux
!=============================================================================
subroutine getfluxforhllc(w,x,ixI^L,ixO^L,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L,iw,idims
double precision, intent(in)  :: w(ixI^S,1:nw)
double precision, intent(out) :: f(ixG^T,1:nwflux)
double precision, intent(in)    :: x(ixI^S,1:ndim)
logical, intent(out)          :: transport
!----------------------------------------------

! assuming getflux FOLLOWS a getaux call for updated p_ entry and
! well-behaved xi=d+tau+p, write:

if(iw==s0_+idims)then
     ! f_i[s_i]=v_i*s_i + p
     f(ixO^S,iw)=w(ixO^S,p_)
{#IFDEF TRACER
{else if (iw==tr^FL_) then 
      f(ixO^S,iw)=zero\}
}
else if(iw==tau_)then
     ! f_i[tau]=v_i*tau + v_i*p
     ! first case only happens when d=0 and p/tau are at enforced minimal
     ! values, namely p=smallp and tau=smallp/(gamma-1) (no velocities)
     ! This will typically NEVER be the case, but still...
     where(w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_)==smallxi)
        f(ixO^S,iw)=zero
     elsewhere
        f(ixO^S,iw)=w(ixO^S,s0_+idims)*w(ixO^S,p_)/ &
             (w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_))
     endwhere
else
     f(ixO^S,iw)=zero
endif

transport=.true.

return
end subroutine getfluxforhllc
!=============================================================================
subroutine con2prim(pressure,lfac,d,s^C,tau,ierror)
!use ieee_arithmetic
use mod_global_parameters

double precision:: pressure,lfac
double precision:: d,s^C,tau
integer:: ierror
      
integer:: ni,niiter
double precision:: govergminone,pcurrent,pnew
double precision:: er,er1,ff,df,dp,v^C
double precision:: pmin,lfac2inv,pLabs,pRabs,pprev
double precision:: s2overcubeG2rh,sqrs
double precision:: xicurrent
double precision:: oldff1,oldff2
double precision:: Nff
double precision:: pleft,pright,pnewi
integer::nit,n2it,ni2,ni3
!-----------------------------------------------------------------------------

ierror=0
! ierror=0 : ok
!
! ierror<>0
!
! ierror=1 : error on entry: must have D>=minrho, tau>=smalltau
! ierror=2 : maxitnr reached without convergence
! ierror=3 : final pressure value < smallp or xi<smallxi during iteration
! ierror=4 : final v^2=1 hence problem as lfac=1/0
! ierror=5 : nonmonotonic function f?
! ierror=7 : stop due to strictnr violation

if(d<minrho .or. tau<smalltau) then
  ierror=1
  return
endif

sqrs={s^C**2.0d0+}

if(sqrs==zero)then
pressure = (eqpar(gamma_)-one)*tau
lfac =one
return
endif


! left and right brackets for p-range
pmin=dsqrt(sqrs)/(one-dmaxvel)-tau-d
pLabs=max(minp,pmin)


if(.false..and.mype==0) then
 print *,'pmin,sqrs,dmaxvel,tau,d,minp'
 print *,pmin,sqrs,dmaxvel,tau,d,minp
endif

pRabs=1.0d99
! start value from input
pcurrent=pLabs
govergminone = eqpar(gamma_)/(eqpar(gamma_)-one)

er1=one
pprev=pcurrent

! Fudge Parameters 
oldff1=1.0d7  ! High number
oldff2=1.0d9  ! High number bigger then oldff1
n2it = 0
nit  = 0

LoopNR:  do ni=1,maxitnr
     nit = nit + 1
     !============= Controle ~1~=============!
     if(nit>maxitnr/4)then
        ! mix pressure value for convergence
        pcurrent=half*(pcurrent+pprev)
        ! relax accuracy requirement
        er1=10.*er1
        nit = nit - maxitnr/10
     endif
     !=======================================!

     niiter=ni  
     xicurrent=tau+d+pcurrent

     if(xicurrent<smallxi) then
      if(strictgetaux)then       
       print*,'!--- amrvacphys/t.srhd-- con2prim ---!'
       print *,'stop: too small xi iterate:',xicurrent
       print *,'for pressure iterate p',pcurrent
       print *,'pressure bracket pLabs pRabs',pLabs,pRabs
       print *,'iteration number:',ni
       print *,'values for d,s,tau,s^2:',d,s^C,tau,sqrs
      end if
      ierror = 3
      return
     endif

     {v^C=s^C/xicurrent\}
     lfac2inv=one - ({v^C**2.0d0+})
     if(lfac2inv>zero) then
       lfac=one/dsqrt(lfac2inv)
     else
      if(strictgetaux)then
       print*,'!--- amrvacphys/t.srhd-- con2prim ---!'
       print *,'stop: negative or zero factor 1-v^2:',lfac2inv
       print *,'for pressure iterate p',pcurrent
       print *,'absolute pressure bracket pLabs pRabs',pLabs,pRabs
       print *,'iteration number:',ni
       print *,'values for d,s,tau,s^2:',d,s^C,tau,sqrs
       print *,'values for v,xi:',v^C,xicurrent
      end if
      ierror=4
      return
     endif
       

     s2overcubeG2rh=sqrs/(xicurrent**3.0d0)
     ff=(d*(lfac-one)-pcurrent-tau)*lfac2inv + pcurrent*govergminone 
     df=(d*lfac - 2.0d0*xicurrent)*s2overcubeG2rh  &
          + govergminone - lfac2inv

     if (ff*df==zero) then
        if (ff==zero) then
            exit ! zero found
        else
            if(strictgetaux)print *,'stop: df becomes zero, non-monotonic f(p)!!'
            ierror=5
            return
        endif
     else 
        pnew=pcurrent-ff/df
        if (ff*df>zero) then
            ! pressure iterate has decreased
            ! restrict to left 
            pnew=max(pnew,pLabs)
        else  ! ff*df<0
            ! pressure iterate has increased
            ! restrict to right 
            pnew=min(pnew,pRabs)
        endif
     endif
        

     ! handle special case where NR incorrectly believes in convergence
     if(pnew == pLabs .and. pcurrent==pnew .and. &
        abs(ff)> absaccnr .and. sqrs > zero)then
        print *,'found pnew=pcurrent=pLabs=',pnew,pcurrent,pLabs
        print *,' while abs(ff)=',abs(ff)
        print *,' while absaccnr=',absaccnr
        print *,' while sqrs=',sqrs
        pnewi=pnew
        ! try 2 higher pressure values to locate a sign change for f(p)
LoopCor:  do ni2=1,2
	   !=====================!
	   pcurrent=pnewi*500.0d0
  	   xicurrent=tau+d+pcurrent
     	   {v^C=s^C/xicurrent\}
     	   lfac2inv=one - ({v^C**2.0d0+})
 	   !=====================!
        
	   !=====================!
    	   if(lfac2inv>zero)then
	      lfac=one/dsqrt(lfac2inv)
	   else
            print *,'factor 1-v^2:',lfac2inv
            print *,'within special case handling',pcurrent
            ierror=4
            return
	   endif
 	   !=====================!

	   s2overcubeG2rh=-sqrs/(xicurrent**3.0d0)

	   !==== Calculate enthalpy and derivative ====!
	   Nff=(d*(lfac-one)-pcurrent-tau)*lfac2inv + pcurrent*govergminone 

	   !== Save old value of pressure ==!
	   pnewi=pcurrent
	   !================================!

	   !== find the interval where is the root ==!
	   if(Nff * ff <=zero)then
	      pnew=pcurrent
	      exit LoopCor
	   endif
	   !=========================================!
        enddo LoopCor

        !== No possible solution, correct all including the conservatives ==!
        if( Nff*ff>zero)then
     
           ! following is in accord with trick done in smallvalues
           d   = 2.0d0*(one + 10.0d0 * minrho) * minrho
           tau = 2.0d0*(one + 10.0d0 * smalltau) * smalltau
           {^C&s^C =zero;}
           pressure     = (eqpar(gamma_)-one)*tau
           lfac = one
 
           if(strictnr)ierror=7 
           ! leave the do loop here
           return
        endif
     endif
     !===============================================!
     dp=pcurrent-pnew
     er=2.0d0*dabs(dp)/(pnew+pcurrent)
     if(((er<tolernr*er1).or.(dabs(dp)<absaccnr))) exit LoopNR
     !===============================================!

     ! For very small values of pressure, NR algorithm is not efficient to
     ! find root, use Euler algorithm to find precise value of pressure
     if((dabs(oldff2-ff) < 1.0d-8 .or. niiter >= maxitnr-maxitnr/20).and.&
	   ff * oldff1 < zero    .and.	dabs(ff)>absaccnr)then

       n2it=n2it+1
       if(n2it <= 3) pcurrent=half*(pnew+pcurrent)
       if(n2it >3)then
         pright =pcurrent
         pleft=pprev
         pcurrent=half*(pleft+pright)
 Dicho:  do ni3=1,maxitnr
	   !===================!
     	   xicurrent=tau+d+pcurrent
     	   {v^C=s^C/xicurrent\}
     	   lfac2inv=one - ({v^C**2.0d0+})
     	   if(lfac2inv>zero)then
	     lfac=one/dsqrt(lfac2inv)
	   else
             print *,'factor 1-v^2:',lfac2inv
             print *,'within Euler case handling',pcurrent
	     ierror=4
             return
	   endif
	   !===================!

	   !== ZM calculation done using the EOS ==!
    	   Nff=(d*(lfac-one)-pcurrent-tau)*lfac2inv + pcurrent*govergminone
 	   !=======================================!
	   !==== Iterate ====!
	   if(ff * Nff < zero)then
     		pleft=pcurrent
  	   else
		pright=pcurrent
    	   endif

	   pcurrent=half*(pleft+pright)
	   !==================!

	   !=== The iteration converge ===!
   	   if(2.0d0*dabs(pleft-pright)/(pleft+pright) < absaccnr &
	      .or. dabs(ff)<absaccnr)then
              pnew=pcurrent
	      exit LoopNR
            endif
	    !==============================!

	    !=== conserve the last value of Nff ===!
	    ff=Nff
	    !======================================!
         enddo    Dicho
       endif

     else
       !====== There is no problems, continue the NR iteration ======!
       pprev=pcurrent
       pcurrent=pnew
       !=============================================================!
     endif 

 
     !=== keep the values of the 2 last ff ===!
     oldff2=oldff1
     oldff1=ff
     !========================================!
enddo LoopNR

if(niiter==maxitnr)then
   !print*,' ff = ',ff,' df = ',df
   !print*,'reachs maxitnr = ', niiter
   ierror=2
   return
endif

if(pcurrent<minp) then 
   ierror=3
   return
endif

!--end result for pressure and lorentz factor------!
pressure=pcurrent
xicurrent=tau+d+pressure
{v^C = s^C/xicurrent\}
lfac2inv=one - ({v^C**2.0d0+})
if(lfac2inv>zero) then
    lfac=one/dsqrt(lfac2inv)
else 
    print *,'factor 1-v^2:',lfac2inv
    print *,'at final getaux check for',pcurrent
    ierror=4
    return
endif
!------------------------------!

end subroutine con2prim
!=============================================================================
subroutine getaux2(clipping,w,x,ixI^L,ixO^L,subname)

! Calculate auxilary variables ixO^L from non-auxiliary entries in w
! however, this time tau replaced by rhos (=DS)

use mod_global_parameters

logical, intent(in)             :: clipping
integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(inout) :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
character(len=*), intent(in)    :: subname

integer          :: err,ix^D
double precision :: dold,rhosold,{^C&sold^C_}
!-----------------------------------------------------------------------------

! we compute auxiliaries p, lfac from D,S,rhos (=DS)
! put the p and lfac in the auxiliary fields lfac_ and p_
! on entry: p_ field may contain previous value for iteration
! however, for filling ghost cells, this does not exist, hence fix here
if(subname=="bc")then
  w(ixO^S,p_)=minp
  w(ixO^S,lfac_)=one
endif

{do ix^D= ixO^LIM^D\}
     dold=w(ix^D,d_)
     rhosold=w(ix^D,rhos_)
    { ^C&sold^C_=w(ix^D,s^C_);}

    call con2prim2(w(ix^D,p_),w(ix^D,lfac_), &
            w(ix^D,d_),{^C&w(ix^D,s^C_)},w(ix^D,rhos_),err)
    if (err/=0) then
       print*,'Getaux2 error:',err,'ix^D=',ix^D,it
       print*,'start value for p=',w(ix^D,p_)
       print*,'end value for lfac=',w(ix^D,lfac_)
       print*,'input d=',w(ix^D,d_),'s=',{^C&w(ix^D,s^C_)},'rhos=',w(ix^D,rhos_)
       print*,'input d=',dold,'s=',{^C&sold^C_},'rhos=',rhosold
       print*,'Called from: ',subname
       call mpistop("problem in getaux2")
   endif
{enddo^D&\}

end subroutine getaux2
!=============================================================================
subroutine con2prim2(pressure,lfac,d,s^C,rhos,ierror)

use mod_global_parameters

double precision:: pressure,lfac
double precision:: d,s^C,rhos
integer:: ierror

double precision:: sqrs,sqrsond2,gmax,gLabs,gRabs,gnew,dg,er
double precision:: gcurrent,govergminone,gprev,er1,xksi,ff,df
integer:: nit, ni, niiter
!-----------------------------------------------------------------------------

ierror=0
! ierror=0 : ok
!
! ierror<>0
!
! ierror=1 : error on entry: must have D>=minrho, rhos>0
! ierror=2 : maxitnr reached without convergence
! ierror=3 : final lfac value < 1 or > maximal value
! ierror=5 : nonmonotonic function f?

if(d<minrho.or.rhos>zero) then
  ierror=1
  return
endif

sqrs={s^C**2.0d0+}

if(sqrs==zero)then
pressure = rhos*d**(eqpar(gamma_)-one)
lfac =one
return
endif

sqrsond2=sqrs/d**2

! left and right bracket for lorentz factor-range
gmax=dsqrt(sqrs/d**2+one)
gLabs=one
gRabs=gmax
! start value 
gcurrent=one
govergminone = eqpar(gamma_)/(eqpar(gamma_)-one)

er1=one
gprev=gcurrent

nit  = 0

LoopNR:  do ni=1,maxitnr
     nit = nit + 1
     !============= Controle ~1~=============!
     if(nit>maxitnr/4)then
        ! mix value for convergence
        gcurrent=half*(gcurrent+gprev)
        ! relax accuracy requirement
        er1=10.*er1
        nit = nit - maxitnr/10
     endif
     !=======================================!

     niiter=ni  

     xksi=govergminone*rhos*d**(eqpar(gamma_)-two)*gcurrent**(one-eqpar(gamma_))
     ff=(one+xksi)**2*(one-gcurrent**2)+sqrsond2
     df=two*gcurrent*(one+xksi)*(xksi*((one-eqpar(gamma_))*(one-gcurrent**2)/gcurrent**2-one)-one)

     if (ff*df==zero) then
        if (ff==zero) then
            exit ! zero found
        else
            print *,'stop: df becomes zero, non-monotonic f(p)!!'
            ierror=5
            return
        endif
     else 
        gnew=gcurrent-ff/df
        if (ff*df>zero) then
            ! iterate has decreased
            ! restrict to left 
            gnew=max(gnew,gLabs)
        else  ! ff*df<0
            ! iterate has increased
            ! restrict to right 
            gnew=min(gnew,gRabs)
        endif
     endif
        

     !===============================================!
     dg=gcurrent-gnew
     er=2.0d0*dabs(dg)/(gnew+gcurrent)
     if(((er<tolernr*er1).or.(dabs(dg)<absaccnr))) exit LoopNR
     !===============================================!

     !======  continue the NR iteration ======!
     gprev=gcurrent
     gcurrent=gnew
     !=============================================================!
 
enddo LoopNR

if(niiter==maxitnr)then
   ierror=2
   return
endif

if(gcurrent<one.or.gcurrent>gmax) then 
   ierror=3
   return
endif

!--end result for pressure and lorentz factor------!
pressure = rhos*d**(eqpar(gamma_)-one)*gcurrent**eqpar(gamma_)
lfac=gcurrent
      
end subroutine con2prim2
!=============================================================================
subroutine addgeometry(qdt,ixI^L,ixO^L,wCT,w,x)

! Add geometrical source terms to w

use mod_global_parameters

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(in)    :: qdt
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(in)    :: wCT(ixI^S,1:nw)
double precision, intent(inout) ::  w(ixI^S,1:nw)

integer          :: iw,idir, h1x^L{^NOONED, h2x^L}
double precision :: tmp(ixG^T)
logical          :: angmomfix=.false.
!-----------------------------------------------------------------------------
!!! this getaux call now ensured in tvdmusclf/hancock/hll variants
!!!if(typeaxial /= 'slab') call getaux(.true.,wCT,ixI^L,ixO^L,'addgeometry')

select case (typeaxial)
case ('slab')
   ! No source terms in slab symmetry
case ('spherical')
     h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
     do iw=1,nwflux
         select case(iw)
         ! s[s1]=((Stheta**2+Sphi**2)/xi+2*p)/r
         case(s1_)
            tmp(ixO^S)=wCT(ixO^S,p_)*x(ixO^S,1) &
                    *(mygeo%surfaceC1(ixO^S)-mygeo%surfaceC1(h1x^S)) &
                    /mygeo%dvolume(ixO^S){&^CE&
                  +wCT(ixO^S,s^CE_)**2.0d0&
                 /(wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_)) }
{^NOONEC
         ! s[s2]=-(Sr*Stheta/xi)/r
         !       + cot(theta)*(Sphi**2/xi+p)/r
         case(s2_)
}
{^NOONED
            tmp(ixO^S) = +wCT(ixO^S,p_)
            w(ixO^S,iw)=w(ixO^S,iw) &
                  +qdt*tmp(ixO^S)*(mygeo%surfaceC2(ixO^S)-mygeo%surfaceC2(h2x^S)) &
                      /mygeo%dvolume(ixO^S)
}
{^NOONEC
            tmp(ixO^S)=-(wCT(ixO^S,s1_)*wCT(ixO^S,s2_)&
                      /(wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_)))
}
{^IFTHREEC
{^NOONED
            tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,s3_)**2.0&
                      /(wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_))) &
                           *dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
}
         ! s[s3]=-(sphi*sr/xi)/r
         !       -cot(theta)*(stheta*sphi/xi)/r
         case(s3_)
            if(.not.angmomfix) &
            tmp(ixO^S)=-(wCT(ixO^S,s3_)*wCT(ixO^S,s1_)&
                      /(wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_))){^NOONED &
                         -(wCT(ixO^S,s2_)*wCT(ixO^S,s3_)&
                      /(wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_))) &
                      *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
}
         end select
         ! Divide by radius and add to w
         if(iw==s1_{^NOONEC.or.iw==s2_}{^IFTHREEC  &
               .or.(iw==s3_.and..not.angmomfix)}) &
            w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
      end do
case ('cylindrical')
   do iw=1,nwflux
      select case (iw)
      ! source[sr]=sphi*vphi/radius + p/radius
      case (sr_)
         w(ixO^S,sr_)=w(ixO^S,sr_)+qdt*wCT(ixO^S,p_)/x(ixO^S,1)
         tmp(ixO^S)=zero
{^IFPHI
         tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,sphi_)**2/ &
                (wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_))
      ! source[sphi]=(-sphi*vr)/radius
      case (sphi_)
         tmp(ixO^S)=-wCT(ixO^S,sphi_)*wCT(ixO^S,sr_)/ &
                (wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_))
}
      end select
      ! Divide by radius and add to w
      if (iw==sr_.or.iw==sphi_) then
           w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
      end if
   end do
end select

end subroutine addgeometry
!=============================================================================
! end module amrvacphys/- srhd
!##############################################################################
