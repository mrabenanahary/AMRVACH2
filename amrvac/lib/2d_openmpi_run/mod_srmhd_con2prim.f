module mod_srmhd_con2prim


 use mod_global_parameters
 use mod_srmhd_parameters
 use mod_srmhd_eos
 implicit none

contains

subroutine srmhd_con2prim(d,tau,flag_use_old,Bsqr,Ssqr,sdotb,sb2,E2,lfac,xi,&
   ierror)
  ! made by Z. MELIANI 14/02/2018
! This subroutine srmhd_con2prim is re-writed by Z. Meliani

! (D,S,tau,B) --> compute auxiliaries lfac and xi
use mod_srmhd_parameters
implicit none

real(kind=dp), intent(inout):: lfac,xi
real(kind=dp), intent(inout):: d,tau
logical      , intent(in)   :: flag_use_old
real(kind=dp), intent(in)   :: Ssqr,Bsqr,sdotb,sb2,E2
integer, intent(out)        :: ierror

! .. local ..
logical         :: flag_use_oldin
real(kind=dp)   :: fin,dfin

real(kind=dp)   :: keepxi1,xi1,xi2,xii,xih,xil,dxi,f,df,xi1bsqr
real(kind=dp)   :: temp,vsqr,fl,fh
real(kind=dp)   :: lfacl,lfach,lfaci,lfacin,lfacmax,lfacmin
real(kind=dp)   :: er,er1,xiprev,dplus
real(kind=dp)   :: fv,dfv,v(1:ndir),pressure
real(kind=dp)   :: lastxiok,lastlfacok
real(kind=dp)   :: lastf,oldF
real(kind=dp)   :: dlfac,p,dpdxi
real(kind=dp)   :: taud
integer         :: ni,nit,niiter
logical         :: finished,tests,testsl,testsh
!-----------------------------------------------------------------------------

ierror=0

! ierror=1 : error on entry: must have D>=0, tau>=smallp/(gamma-1)
! ierror=2 : srmhd_maxiterationNR reached without convergence
! ierror=3 : no range for solution was find
! ierror=4 : final xi <small_xi
! ierror=5 : lower bound non-negative f(xi) value
! ierror=6 : v^2>1
! ierror=7 :  solution out from the range?
! ierror=8 : nonmonotonic function f?

if(d<small_density .or. tau<small_e)then
  ierror=1
  return
endif
taud=tau+d


! bsqr=sum(b**2.0d0)
! ssqr=sum(s**2.0d0)

! handle the case of no flow first
cond_noflow : if(ssqr<=small_vec2)then
  call srmhd_get_h_noflux(d,taud-0.5*bsqr,xi)
  lfac=1.0_dp
  return
endif cond_noflow





! Hydro case: handle exactly as before (using p to iterate on)
!{#IFDEF ENERGY
! opedit: I don't know how this would be d1.0_dp if we don't have energy so removing for now.
! cond_srhd: if(bsqr<=small_vec2)then
!   pressure=xi-taud
!   call srmhd_con2prim_hd(d,tau,taud,flag_use_old,Ssqr,&
!                           lfac,pressure,ierror)
!   xi=pressure+taud
!   return
! end if cond_srhd
!}




! temp=sum(s*b)
! if(temp<0.0_dp)then
!   sdotb=-min(dsqrt(ssqr*bsqr),-temp)
! else
!   sdotb=min(dsqrt(ssqr*bsqr),temp)
! endif
!sb2=sdotb**2.0




cond_useold_0 : if(flag_use_old)then
 if (xi>smalldouble) then
    flag_use_oldin=.true.

    call srmhd_funcd(xi,fin,dfin,lfacin,d,ssqr,tau,taud,bsqr,sdotb,sb2,E2,&
       ierror)
    if(dabs(fin)<srmhd_absaccnr.and.ierror/=0)then
      return ! no change
    end if
 else
   flag_use_oldin=.false.
 end if
else
  flag_use_oldin =.false.
end if cond_useold_0


fin=1.0d99

! the maximal allowed Lorentz factor
lfacmax=1.0_dp/dsqrt(1.0_dp-srmhd_maxspeed2)
! starting point for NR on xi
! xi1 is meant to be the lower bound for the bracket to be used in NR
! estimate inyial value for xi
! solve cubic equation

dplus=d+srmhd_gamma*small_pressure/srmhd_gamma_1
xi1=dplus
lfacmin=1.0d0
lfacl=lfacmax



call srmhd_initialise_con2prim(d,tau,taud,small_pressure,bsqr,ssqr,sdotb,E2,&
   lfacmin,lfacmax,dplus,xi1,lfacl,ierror)

 xi1bsqr=xi1+bsqr
 vsqr = ((ssqr + ((xi1+xi1bsqr)*sdotb**2.0_dp)/xi1**2.0_dp)/ &
    (xi1bsqr**2.0_dp))
if(vsqr<srmhd_maxspeed2)then
    call srmhd_funcd(xi1,F,dF,lfac,d,ssqr,tau,taud,bsqr,sdotb,sb2,E2,ierror)
    if(f>0.0_dp)xi1=dplus
end if

! compute v^2(xi1)
keepxi1=-1.0_dp;oldF=1.0d99
cond_sdotB : if(dabs(sdotb)<srmhd_limitvalue)then
 xi1=max(xi1,dsqrt(ssqr)/srmhd_maxspeed-bsqr)
else cond_sdotB
 xi1bsqr=xi1+bsqr
 vsqr = ((ssqr + ((xi1+xi1bsqr)*sdotb**2.0d0)/xi1**2.0d0)/ (xi1bsqr**2.0d0))
 !=== find new xi1, in the case v2(xi1) > maxvel^2 = (1-dmaxvel)^2
 ! locate xi corresponding to maximal velocity allowed, namely 1-dmaxvel
 niiter=-1
 cond_speed1 : if(vsqr > srmhd_maxspeed2 )then
   er1=1.0_dp
   xiprev=xi1
  LoopVmax:  do ni = 1,srmhd_maxiterationNR

      if(ni>srmhd_maxiterationNR/2)then
        xi1=half*(xi1+xiprev)
        er1=10.0d0*er1
      endif
      xi1bsqr=xi1+bsqr
      ! v^2(xi1) - maxvel^2
      fv= (ssqr + (xi1+xi1bsqr)*sb2/xi1**2.0d0)/xi1bsqr**2.0d0-srmhd_maxspeed2
      ! d(v^2(xi)-maxvel^2)/dxi
      dfv= -two * (sb2*(3.0d0*xi1*xi1bsqr+bsqr**2.0d0)+ssqr*xi1**3)/ &
         ((xi1*xi1bsqr)**3.0d0)

      cond_isone : if(fv*dfv==0.0_dp) then
         if(fv==0.0_dp)then
            exit LoopVmax
         else
            ierror=8
            return
         endif
      else cond_isone
        xiprev=xi1
        xi1   =xi1 -fv/dfv
        cond_Rfv : if(fv*dfv>0.0_dp)then
          ! xi-iterate decreased
          ! restrict to left
          xi1=max(xi1,dplus)
        else  cond_Rfv ! fv*dfv <0
          ! xi-iterate increased
          ! restrict to right
          xi1=min(xi1,taud)
        endif cond_Rfv
      endif cond_isone
     ! TESSTSSS=========
     vsqr = ((ssqr + ((xi1+xi1bsqr)*sb2)/xi1**2.0d0)/(xi1bsqr**2.0d0))
     if(vsqr<srmhd_maxspeed2)then

      call srmhd_funcd(xi1,F,dF,lfac,d,ssqr,tau,taud,bsqr,sdotb,sb2,E2,ierror)
      if(dabs(F)<oldF.and.f<0.0_dp)then
       keepxi1=xi1
       oldF=dabs(F)
      end if
     end if
     !================
      er=dabs(fv/dfv)/xi1
      if((er<srmhd_tolernr*er1).or.(dabs(fv/dfv)<srmhd_absaccnr))exit LoopVmax
      niiter=ni
   enddo LoopVmax
 endif cond_speed1
 !=============================================!


 cond_noconverge0 :if(niiter==srmhd_maxiterationNR)then

  ! could not find consistent value of lower bound for xi compliant with maxvel
  !print *,' could not find value of lower bound for xi compliant with maxvel'
  ! print *,'THE FIRST ','xi1=',xi1,'dplus=',dplus,' tau+d=',tau+d,oldF,xi1
  ierror=2
  return
 endif cond_noconverge0

 ! we now compute f(xi1) and lfac(xi1)
 if(keepxi1>0.0_dp)xi1=keepxi1
end if cond_sdotB


call srmhd_funcd(xi1,fl,df,lfacl,d,ssqr,tau,taud,bsqr,sdotb,sb2,E2,ierror)
if(ierror /=0) then
 !Print*,'START IERROR', IERROR,niiter,er,srmhd_tolernr*er1,er1,dabs(fv/dfv),srmhd_absaccnr
 return
end if
testsl=(xi1>=dplus*lfacl.and.lfacl<=lfacmax)

if(fl>0.0_dp)then
  ierror=5
  return
endif

!--------------------------------------------------------!
! xi2 is meant to be the maximal bound for the bracket to be used in NR
! for this we take the value of E=tau+d, increased by smallp

xi2= max(taud+small_pressure - half *bsqr,10.0d0*xi1)

niiter=-1

LoopxiMax : do ni=1,srmhd_maxiterationNR
    ! we now compute f(xi2) and lfac(xi2)

    call srmhd_funcd(xi2,fh,df,lfach,d,ssqr,tau,taud,bsqr,sdotb,sb2,E2,ierror)

    if(ierror /=0)then
      return
    end if
    testsh=(xi2>=dplus*lfach.and.lfach<=lfacmax)

    ! maximal bound found when f(xi1) opposite sign of f(xi2)
    !, enforce consistency tests on xi2
    if (testsh.and.fh *fl <=0.0_dp) exit LoopxiMax
    !==== Zak 17/05 fast convergence ====!
    xi1=xi2
    fl=fh
    lfacl=lfach
    testsl=testsh
    !====================================!
    xi2=two*xi2
    niiter=ni

end do   LoopxiMax

if(niiter>1)call srmhd_funcd(xi1,fl,df,lfacl,d,ssqr,tau,taud,bsqr,sdotb,sb2,E2,&
   ierror)
if(flag_use_oldin)then
 if(fin*fl<0.0_dp)then
  if(xi2>xi)then
    xi2=xi;lfach=lfacin;fh=fin
  end if
 else
  if(xi1<xi)then
    xi1=xi;lfacl=lfacin;fl=fin
  end if
 end if
end if

!--------------------------------------------------------!

if(niiter == srmhd_maxiterationNR .or. (fh*fl>0.0_dp))then
  ierror = 3
  return
end if

finished=.false.
if(fl==0.0_dp)then
  xii=xi1
  lfaci=lfacl
  finished=(xi1>=dplus*lfacl.and.lfacl<=lfacmax)
endif

if(fh==0.0_dp)then
  xii=xi2
  lfaci=lfach
  finished=(xi2>=dplus*lfach.and.lfach<=lfacmax)
endif

if(finished)then
  xi=xii
  lfac=lfaci
  return
end if


xil=xi1
xih=xi2
xii=(dabs(fl)*xih+dabs(fh)*xil)/(dabs(fl)+dabs(fh)) !Initialize the guess for rootfinder
call srmhd_funcd(xii,f,df,lfaci,d,ssqr,tau,taud,bsqr,sdotb,sb2,E2,ierror)
if(ierror /=0)then

  return
end if

call  Newton_Ralphson(srmhd_funcd,d,ssqr,tau,taud,bsqr,sdotb,sb2,E2,lfacmax,&
   dplus, xil,xih,xii,fl,fh,xi,lfac,ierror)

if(ierror/=0)then
 return
end if





!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






end subroutine srmhd_con2prim

!================================================================================
subroutine Newton_Ralphson(fmhd_l,d,ssqr,tau,taud,bsqr,sdotb,sb2,E2,lfacmax,&
   dplus, xil,xih,xii,fl,fh,xi,lfac,ierror)

use mod_global_parameters
use mod_srmhd_parameters
implicit none
procedure(gene_funcd)           :: fmhd_l
real(kind=dp)   , intent(in)    :: d,ssqr,tau,taud,bsqr,sdotb,sb2,E2,lfacmax,&
   dplus
real(kind=dp)   , intent(inout) :: xil,xih,fl,fh
real(kind=dp)   , intent(inout) :: xii
real(kind=dp)   , intent(out)   :: xi,lfac
integer, intent(inout)          :: ierror

integer          :: nit,niiter,ni
real(kind=dp)    :: er1,xiprev,lastxiok,lastlfacok,lastf,lfaci,er
real(kind=dp)    :: f,df
logical          :: tests
! ------------------------------------------------------
er1 = 1.0_dp
nit = 0
niiter=-1
xiprev=xii
lastxiok=-1.0_dp
lastlfacok=-1.0_dp
lastf=-1.0_dp
!--- Start iteration ---!

LoopNRRMHD :  do ni=1,srmhd_maxiterationNR
      nit = nit + 1

      if(nit>srmhd_maxiterationNR/2)then
        ! mix the last  value for convergence
        xii=(dabs(fl)*xih+dabs(fh)*xil)/(dabs(fl)+dabs(fh))
        ! relax accuracy requirement
        er1=10.0d0*er1
        ! following avoids decrease of accuracy requirement
        ! *every* iteration step beyond srmhd_maxiterationNR/2
        nit = nit - srmhd_maxiterationNR/10
      endif

      call fmhd_l(xii,f,df,lfaci,d,ssqr,tau,taud,bsqr,sdotb,sb2,E2,ierror)

      if(ierror /=0)then
              return
      end if
      tests=(xii>=dplus*lfaci.and.lfaci<=lfacmax)
      if(tests)then
        lastxiok=xii
        lastlfacok=lfaci
        lastf=f
      endif


      cond_iszero : if(f*df==0.0_dp) then
         cond_fiszero : if(f==0.0_dp) then
            tests=(xii>=dplus*lfaci.and.lfaci<=lfacmax)
            cond_zerophys0 : if(tests)then
               exit ! 0.0_dp found and consistency checks fullfilled
            else cond_zerophys0
              if(lfaci<=lfacmax)then
               ! use df as enthalpy
               df=xi/lfaci
               xiprev=d/lfaci
               call srmhd_get_pressure_primitive_scalar(xiprev,xi/lfaci,f)
               !call FuncPressure(xi,lfaci,d,ssqr,tau,1.0d0,f,df)
               if(f>0)exit ! 0.0_dp found and consistency checks fullfilled
              end if
              ierror=7
              return
            endif cond_zerophys0
         else cond_fiszero
            ierror=8
            return
         endif cond_fiszero
      else cond_iszero
        xiprev=xii
        xii   =xii -f/df
        cond_limitf1: if(f*df>0.0_dp)then
          ! xi-iterate decreased
          ! restrict to left
          xii=max(xii,xil)
        else cond_limitf1
         ! xi-iterate increased
         ! restrict to right
          xii=min(xii,xih)
        endif cond_limitf1
        cond_signf : if(fl*f>0.0)then
          if(xii>xil)then
            xil=xiprev
            fl=f
          end if
        else cond_signf
          if(xii<xih)then
            xih=xiprev
            fh=f
          end if
        end if cond_signf
      endif cond_iszero

      er=dabs(f/df)/xii
      cond_converge : if((er<srmhd_tolernr*er1).or.&
         (dabs(f/df)<srmhd_absaccnr))then
            call fmhd_l(xii,f,df,lfaci,d,ssqr,tau,taud,bsqr,sdotb,sb2,E2,&
               ierror)
            tests=(xii>=dplus*lfaci.and.lfaci<=lfacmax)
            if(tests)then
               exit LoopNRRMHD ! converged solution with tests ensured
            else
               if(lfaci<=lfacmax)then
               df=xi/lfaci
               xiprev=d/lfaci
               call srmhd_get_pressure_primitive_scalar(xiprev,xi/lfaci,f)
               if(f>0)exit  LoopNRRMHD ! converged solution with tests ensured
              end if
              ierror=7
              return
            endif
      endif cond_converge
      niiter=ni

enddo LoopNRRMHD
srmhd_maxiter_nr=max(srmhd_maxiter_nr,ni)
if(niiter==srmhd_maxiterationNR) then
     ierror=2
     return
endif ! niiter==srmhd_maxiterationNR

 !===============================!
 ! final values for auxiliary variables are now passed to w-array
  xi=xii
  lfac=lfaci
 !===============================!
 ! we now perform some additional consistency checks

  if(xi<small_xi)then
    ierror=4
    return
  endif
end subroutine Newton_Ralphson

!----------------------------------------------------------------------------
subroutine srmhd_funcd(xi,F,dF,lfac,d,ssqr,tau,taud,bsqr,sdotb,sb2,e2,ierror)

real(kind=dp)   , intent(in)  :: xi,d,ssqr,tau,taud,bsqr,sdotb,sb2,e2
real(kind=dp)   , intent(out) :: F,dF,lfac
integer, intent(inout)        :: ierror
! .. local ..
real(kind=dp)     :: dlfacdxi,rhoh,dhdxi,rho,drhodxi,xibsqr
real(kind=dp)     :: vsqr,p,dpdxi,invlfac
!-----------------------------------------------------------------------------
xibsqr=xi+bsqr

vsqr = (ssqr + (2.0_dp*(xibsqr)*sb2)/(xi**2.0_dp))/(xibsqr**2.0_dp)
if (vsqr<1.0_dp) then
   lfac    = 1.0_dp/dsqrt(1.0_dp-vsqr)
   invlfac = 1.0_dp/lfac
   rho     = d/lfac
   rhoh    = xi/lfac**2.0_dp


   dlfacdxi = -lfac**3.0_dp*(sb2*(3.0d0*xi*xibsqr+bsqr**2.0_dp)+&
      ssqr*xi**3.0_dp)/((xi*xibsqr)**3.0_dp)

   drhodxi  = -rho*dlfacdxi*invlfac
   dhdxi    = invlfac**2.0_dp*(1.0_dp-2.0_dp*xi*invlfac*dlfacdxi)

   call srmhd_get_val_p_dpdxi(rho,rhoh,drhodxi,dhdxi,p,dpdxi)



   F  = xibsqr-(taud+p+half*bsqr)+half*E2/xibsqr**2.0_dp
   dF = 1.0_dp + E2/xibsqr**3.0_dp  -dpdxi


!   F  = xibsqr-tau-d+half*bsqr*vsqr-half*sb2/xi**2.0d0-p
!   dF = 1.0_dp + sb2/xi**3.0d0  -dpdxi+dlfacdxi*bsqr/lfac**3.0d0
else
  ! print *,'Warning: err1.0_dpous input to funcd since vsrq=',vsqr,' >=1'
   print *,'input values d, ssqr, tau, bsqr, sdotb:',d,ssqr,tau,bsqr,sdotb
   print*,'ierror ==6 ',it,mype,global_time,saveigrid
   ierror =6
   call mpistop('at srmhd_funcd')
   return
end if

end subroutine srmhd_funcd
!=============================================================================

!=============================================================================
subroutine srmhd_initialise_con2prim(d,tau,taud,pressure,bsqr,ssqr,sdotb,E2,&
   lfacmin,lfacmax,dplus,xi,lfac,ierror)

use mod_global_parameters
implicit none


real(kind=dp)   , intent(in)    :: d,tau, pressure,taud
real(kind=dp)   , intent(in)    :: ssqr,bsqr,sdotb,E2,lfacmin,lfacmax,dplus
real(kind=dp)   , intent(inout) :: xi,lfac
integer, intent(inout)          :: ierror
! .. local ..
real(kind=dp)                   :: xibsqr,taudp,F,dF,lfacn,dxidlfac,er,er1,rho,&
   rhoh,a0,b0,x0,dx0
integer                         :: ni,n_h
real(kind=dp)   , parameter     :: smallratio=1.0d-6
!-----------------------------------------------------------------------------
taudp=taud+pressure
er1=1.0
n_h=0
lfacn=lfac


a0=E2/(2.0_dp*(taudp+bsqr/2.0_dp))
b0=a0/((taudp+bsqr/2.0_dp)**2.0_dp)
if(a0<smallratio)then
  xi=max(dsqrt(a0)-bsqr,dplus)
  lfacn=max(1.0_dp,min(xi/d,lfacmax))
end if

if(1.0_dp/a0<smallratio)then
 xi=max(taudp-half*bsqr,dplus)
 lfacn=max(1.0_dp,min(xi/d,lfacmax))
end if

Loop_lfac_max : do ni=1,srmhd_maxiterationNR
 if(ni>srmhd_maxiterationNR/2)then
    er1=10.0d0*er1
 endif
 rho  = d/lfacn

 call srmhd_get_enthalpy_scalar(rho,pressure,rhoh)

 xi=lfacn**2.0*rhoh
 dxidlfac=2.0*lfacn*rhoh-d
 xibsqr=xi+bsqr
 x0=xibsqr/(taud+half*bsqr)
 dx0=dxidlfac/(taud+half*bsqr)
 F = x0-1.0_dp+b0/x0**2.0d0
 dF = dx0* (1.0_dp- two*b0/x0**3.0d0)


 cond_fdfzero : if(f*df==0.0_dp) then
   if(f==0.0_dp)then
     exit Loop_lfac_max
   else
     !print *,'stop: dfv becomes 0.0_dp, non-monotonic function of xi!!'
     ierror=10
     return
   endif
 else cond_fdfzero
    lfacn   =lfacn -f/df
    if(f*df>0.0_dp)then
     ! lfac-iterate decreased
     ! restrict to left
     if(lfacn<=lfacmin)then
       lfacn=lfacmin
       n_h=n_h+1
     else
       n_h=0
     end if
     lfacn=max(lfacn,lfacmin)
    else ! fv*dfv <0
     ! lfac-iterate increased
     ! restrict to right
     if(lfacn>=lfacmax)then
       lfacn=lfacmax
       n_h=n_h+1
     else
       n_h=0
     end if
    endif
 end if cond_fdfzero
  if(n_h>=3)exit Loop_lfac_max
  er=dabs(f/df)/lfacn
  if((er<srmhd_tolernr*er1).or.(dabs(f/df)<srmhd_absaccnr))exit Loop_lfac_max

end do Loop_lfac_max
!maxiter_ini=max(ni,maxiter_ini)
end subroutine srmhd_initialise_con2prim
!=============================================================================

!=============================================================================
subroutine srmhd_funcd_hd(xi,F,dF,lfac,d,ssqr,tau,taud,ierror)

real(kind=dp)   , intent(in)  :: xi,d,ssqr,tau,taud
real(kind=dp)   , intent(out) :: F,dF,lfac
integer, intent(inout)        :: ierror
! .. local ..
real(kind=dp)     :: dlfacdxi,rhoh,dhdxi,rho,drhodxi
real(kind=dp)     :: vsqr,p,dpdxi,invlfac
!-----------------------------------------------------------------------------

vsqr = ssqr /xi**2.0_dp
if (vsqr<1.0_dp) then

   lfac    = 1.0_dp/dsqrt(1.0_dp-vsqr)
   invlfac = 1.0/lfac
   rho     = d/lfac
   rhoh    = xi/lfac**2.0_dp


   dlfacdxi = -(lfac/xi)**3.0_dp*ssqr*xi**3.0_dp
   drhodxi  = -rho*dlfacdxi*invlfac
   dhdxi    = invlfac**2.0_dp*(1.0_dp-2.0_dp*xi*invlfac*dlfacdxi)

   call srmhd_get_val_p_dpdxi(rho,rhoh,drhodxi,dhdxi,p,dpdxi)



   F  = xi-(taud+p)
   dF = 1.0_dp-dpdxi

else
  ! print *,'Warning: err1.0_dpous input to funcd since vsrq=',vsqr,' >=1'
   print *,'input values d, ssqr, tau=:',d,ssqr,tau
   print*,'ierror ==6 ',it,mype,global_time,saveigrid
   ierror =6
   !call mpistop('at srmhd_funcd')
   return
end if

end subroutine srmhd_funcd_hd
!=============================================================================


subroutine srmhd_con2prim_hd(d,tau,taud,flag_use_old,Ssqr,lfac,pressure,&
   ierror)
!use ieee_arithmetic
use mod_global_parameters
implicit none
real(kind=dp)   , intent(inout) :: pressure,lfac
real(kind=dp)   , intent(in)    :: d,tau,taud,Ssqr
logical, intent(in)             :: flag_use_old
integer, intent(inout)          :: ierror

! .. local ..
integer:: ni,niiter
real(kind=dp)   :: pcurrent,pnew,pL,pR
real(kind=dp)   :: er,er1,ff,df,dpressure
real(kind=dp)   :: pmin,lfac2inv,pLabs,pRabs,pprev
real(kind=dp)   :: xicurrent,xi
real(kind=dp)   :: oldff1,oldff2
real(kind=dp)   :: Nff,fin
real(kind=dp)   :: pleft,pright,pnewi

integer::nit,n2it,ni2,ni3

logical                     :: flag_use_oldin
!-----------------------------------------------------------------------------

ierror=0
! ierror=0 : ok
!
! ierror<>0
!
! ierror=1 : error on entry: must have D>=small_density, tau>=small_e
! ierror=2 : srmhd_maxiterationNR reached without convergence
! ierror=3 : final pressure value < smallp or xi<smallxi during iteration
! ierror=4 : final v^2=1 hence problem as lfac=1/0
! ierror=5 : nonmonotonic function f?
! ierror=7 : stop due to strictnr violation

 if(d<small_density .or. tau<small_e) then
  ierror=1
  return
 endif

if(ssqr<=small_vec2)then
  call srmhd_get_h_noflux(d,taud,xi)
  lfac=1.0_dp
  return
endif



 fin=1.0d99

! check old auxilairies
 flag_use_oldin=.false.
 if(flag_use_old)then
  if (pressure>small_pressure) then
     xicurrent=taud+pressure
     flag_use_oldin=.true.
     call  srmhd_funcd_hd(xicurrent,Fin,dF,lfac,d,ssqr,tau,taud,ierror)
     if(dabs(fin)<srmhd_absaccnr)return ! no change
  end if
 end if

 ! left and right brackets for p-range
 pmin=dsqrt(Ssqr/srmhd_maxspeed2)-taud
 pLabs=max(small_pressure,pmin)
 pRabs=1.0d99
 ! start value from input
 pcurrent=pLabs
 xicurrent=taud+pLabs

 if(pcurrent<small_pressure) then
       ierror=3
       return
 endif


 call  srmhd_funcd_hd(xicurrent,fin,dF,lfac,d,ssqr,tau,taud,ierror)

 if(dabs(fin)<dabs(ff))then
   pcurrent=pressure
 else
   pcurrent=pLabs
 end if
 er1=1.0_dp
 pprev=pcurrent

 ! Fudge Parameters
 oldff1=1.0d7  ! High number
 oldff2=1.0d9  ! High number bigger then oldff1
 n2it  = 0
 nit   = 0

 LoopNR:  do ni=1,srmhd_maxiterationNR
     nit = nit + 1
     !============= Controle ~1~=============!
     if(nit>srmhd_maxiterationNR/4)then
        ! mix pressure value for convergence
        pcurrent=half*(pcurrent+pprev)
        ! relax accuracy requirement
        er1=10.*er1
        nit = nit - srmhd_maxiterationNR/10
     endif
     !=======================================!

     niiter=ni
     if(pcurrent<small_pressure) then
      ierror = 3
      return
     endif
     xicurrent=taud+pcurrent
     call  srmhd_funcd_hd(xicurrent,ff,dF,lfac,d,ssqr,tau,taud,ierror)
     if (ff*df==0.0_dp) then
        if (ff==0.0_dp) then
            exit ! 0.0_dp found
        else
            print *,'stop: df becomes 0.0_dp, non-monotonic f(p)!!'
            ierror=5
            return
        endif
     else
        pnew=pcurrent-ff/df
        if (ff*df>0.0_dp) then
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
  cond_noconverge1: if(pnew == pLabs .and. pcurrent==pnew .and. abs(ff)> &
     srmhd_absaccnr .and. Ssqr > 0.0_dp)then
        pnewi=pnew
        ! try 2 higher pressure values to locate a sign change for f(p)
    LoopCor:  do ni2=1,2
     !=====================!
     pcurrent=pnewi*500.0d0
     xicurrent=taud+pcurrent


     !==== Calculate enthalpy and derivative ====!
     call  srmhd_funcd_hd(xicurrent,Nff,dF,lfac,d,ssqr,tau,taud,ierror)

     !== Save old value of pressure ==!
     pnewi=pcurrent
     !================================!

     !== find the interval where is the root ==!
     if(Nff * ff <=0.0_dp)then
        pnew=pcurrent
        exit LoopCor
     endif
     !=========================================!
    enddo LoopCor

        !== No possible solution, correct all including the conservatives ==!
        if( Nff*ff>0.0_dp)then
           ! following is in accord with trick d1.0_dp in smallvalues
           pressure = (srmhd_gamma-1.0_dp)*2.0d0*(1.0_dp + 10.0d0 * small_e) * &
              small_e
           lfac     = 1.0_dp
           ierror=7
           ! leave the do loop here
           return
        endif
  endif cond_noconverge1
  !===============================================!
  dpressure=pcurrent-pnew
  er=2.0d0*dabs(dpressure)/(pnew+pcurrent)
  if(((er<srmhd_tolernr*er1).or.(dabs(dpressure)<srmhd_absaccnr))) exit LoopNR
  !===============================================!

  ! For very small values of pressure, NR algorithm is not efficient to
  ! find root, use Euler algorithm to find precise value of pressure
  if((dabs(oldff2-ff) < 1.0d-8 .or. niiter >= &
     srmhd_maxiterationNR-srmhd_maxiterationNR/20).and.ff * oldff1 < 0.0_dp    &
     .and.  dabs(ff)>srmhd_absaccnr)then

   n2it=n2it+1
   if(n2it <= 3) pcurrent=half*(pnew+pcurrent)
   if(n2it >3)then
     pright =pcurrent
     pleft=pprev
     pcurrent=half*(pleft+pright)
     Dicho:  do ni3=1,srmhd_maxiterationNR
      !===================!
      xicurrent=taud+pcurrent
      call  srmhd_funcd_hd(xicurrent,Nff,dF,lfac,d,ssqr,tau,taud,ierror)
       !=======================================!
      !==== Iterate ====!
      if(ff * Nff < 0.0_dp)then
          pleft=pcurrent
      else
         pright=pcurrent
      endif

      pcurrent=half*(pleft+pright)
     !==================!

     !=== The iteration converge ===!
       if(2.0d0*dabs(pleft-pright)/(pleft+pright) < srmhd_absaccnr .or. &
          dabs(ff)<srmhd_absaccnr)then
              pnew=pcurrent
        exit LoopNR
       endif
      !==============================!

      !=== conserve the last value of Nff ===!
      ff=Nff
      !======================================!
    enddo Dicho
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
!maxiter_nr=max(maxiter_nr,ni)
if(niiter==srmhd_maxiterationNR)then
   !print*,' ff = ',ff,' df = ',df
   !print*,'reachs srmhd_maxiterationNR = ', niiter
   ierror=2
   return
endif

if(pcurrent<small_pressure) then
   ierror=3
   return
endif

!--end result for pressure and lorentz factor------!
pressure=pcurrent
xicurrent=taud+pressure

lfac2inv=1.0_dp - Ssqr/xicurrent
if(lfac2inv>0.0_dp) then
    lfac=1.0_dp/dsqrt(lfac2inv)
else
    ierror=4
    return
endif
!------------------------------!

end subroutine srmhd_con2prim_hd
!
end module mod_srmhd_con2prim
