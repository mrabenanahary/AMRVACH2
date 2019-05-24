!===========================================================================
!##############################################################################
! module USE METRIC - GENERAL RELATIVISTY -STATIC METRIC
!=============================================================================
! Project : GENERAL RELATIVISTY -STATIC METRIC : Conformal Decomposition 
! Aim     :
! Ref     :
! pre-compile: setamrvac -p=gr
! update : 10/02/2015,  Zakaria
! based on  geometry.t 
! configuration
! typemetric
! typeaxial
! parameters :
! variable :
! main variables

!============================================================================
Module mod_usemetric
!============================================================================
! global variable in the module
!Implicit none
{^IFGR
!... use modules ...
 use mod_interface_metric
 use mod_metric
!...................
include 'amrvacdef.f'
!-------------------------------------
type vector_phys
 double precision :: w(ixG^T,1:^NC)
end type vector_phys
type vecnorm
 double precision         :: norm2(ixG^T)
end type vecnorm
type sptmetric
 double precision         :: wg(ixG^T)
end type sptmetric
type tensor_vf
 double precision         :: w(ixG^T)
 logical                  :: elem_needed
end type tensor_vf

!------------------------------------- 
!=============================================================================
contains
}
{^IFGR
INCLUDE: mod_metric.BL_to_KS.t
}
{^IFGR
!=============================================================================
!=============================================================================
subroutine usemetric_conserve(ixI^L,ixO^L,patchw,subname,w,nwv,iwv,it_cov,keep_withoutG)
! w should be contravariant

 integer, intent(in)             :: ixI^L,ixO^L
 logical, intent(in)             :: patchw(ixG^T)
 character(len=*), intent(in)    :: subname
 double precision, intent(inout) :: w(ixI^S,1:nw)
 integer, intent(in)             :: nwv,iwv(nwv)
 logical, optional               :: it_cov
 logical, optional               :: keep_withoutG
 ! .. local variables ..
 integer              :: idir,jdir,kdir,iidir,jjdir,kkdir,iv
 type(vector_phys)    :: vcont(nwv)
 logical              :: it_cov_sub,nozero
!-----------------------------------------------------------------------------
cond_cov1: if(present(it_cov))then
 it_cov_sub=it_cov
else
 it_cov_sub=.false.
end if cond_cov1
cond_cov: if(.not.it_cov_sub)then
 Loop_iniv : do iv=1,nwv
  where(patchw(ixO^S))
   ^C&vcont(iv)%w(ixO^S,^C)=w(ixO^S,iwv(iv)+^C);
  end where
 end do Loop_iniv
 Loop_idir : do iidir=1,setgr%ndir
   nozero=.false.
   if(all((/(.not.mypm%elemsp(i1_gr(iidir,jdir),i2_gr(iidir,jdir))%elpt%elm_on&
            ,jdir=1,setgr%ndir)/)))then
    Loop_setv0 : do iv=1,nwv
     where(patchw(ixO^S))
      w(ixO^S,iwv(iv)+iidir)=zero;
     end where
    end do Loop_setv0
    cycle Loop_idir
   end if

  Loop_jdir : do jjdir=1,setgr%ndir;idir=i1_gr(iidir,jjdir);
   jdir=i2_gr(iidir,jjdir)
   cond_on : if(mypm%elemsp(idir,jdir)%elpt%elm_on)then
     cond_isone : if(.not.mypm%elemsp(idir,jdir)%elpt%elm_isone)then
      cond_nozero1 : if(nozero)then
       Loop_sumv1add : do iv=1,nwv
        where(patchw(ixO^S))
         w(ixO^S,iwv(iv)+iidir)=w(ixO^S,iwv(iv)+iidir)+&
                          vcont(iv)%w(ixO^S,jjdir)&
                         *mypm%elemsp(idir,jdir)%elpt%wg(ixO^S);
        end where
        end do Loop_sumv1add
      else cond_nozero1
       Loop_sumv1 : do iv=1,nwv
        where(patchw(ixO^S))
         w(ixO^S,iwv(iv)+iidir)=vcont(iv)%w(ixO^S,jjdir)&
                          *mypm%elemsp(idir,jdir)%elpt%wg(ixO^S);
        end where
       end do Loop_sumv1
      end if  cond_nozero1
     else cond_isone
      cond_nozero2 : if(nozero)then
        Loop_sumv1badd : do iv=1,nwv
         where(patchw(ixO^S))
          w(ixO^S,iwv(iv)+iidir)=w(ixO^S,iwv(iv)+iidir)+&
                                  vcont(iv)%w(ixO^S,jjdir);
         end where
        end do Loop_sumv1badd
      else cond_nozero2
       Loop_sumv1b : do iv=1,nwv
        where(patchw(ixO^S))
         w(ixO^S,iwv(iv)+iidir)=vcont(iv)%w(ixO^S,jjdir);
        end where
       end do Loop_sumv1b
     end if cond_nozero2
    end if cond_isone
    nozero=.true.
   end if cond_on
  end do Loop_jdir
 end do Loop_idir
 
end if cond_cov
if(present(keep_withoutG))then
 if(.not.keep_withoutG)then
  where(patchw(ixO^S))
   ^FL&w(ixO^S,^FL)=w(ixO^S,^FL)*mypm%sqrtdetr%Gama(ixO^S);
  end where
 end if
else
  where(patchw(ixO^S))
   ^FL&w(ixO^S,^FL)=w(ixO^S,^FL)*mypm%sqrtdetr%Gama(ixO^S);
  end where
end if
end subroutine usemetric_conserve

!=============================================================================
!=============================================================================
subroutine  usemetric_getaux_ready(ixI^L,ixO^L,patchw,subname,w,&
                                   nwv_cont,iwv_cont,nwv_cov,iwv_cov,&
                                   pscal_cov,pscal_cont,&
                                   wphys,withoutG)
 ! w is covariant
 include 'amrvacdef.f'
 integer, intent(in)             :: ixI^L,ixO^L
 logical, intent(in)             :: patchw(ixG^T)
 character(len=*), intent(in)    :: subname
 double precision, intent(inout) :: w(ixI^S,1:nw)
 integer, intent(in)             :: nwv_cont,iwv_cont(:)
 integer, intent(in)             :: nwv_cov,iwv_cov(:)
 type(vecnorm), optional         :: pscal_cov(nwv_cov)
 type(vecnorm), optional         :: pscal_cont(nwv_cont)
 double precision, optional      :: wphys(ixG^T,1:nw)
 logical,optional                :: withoutG
! .. local variables ..
 integer                          :: idir,jdir,iw,iidir,jjdir,iv
 type(vector_phys),allocatable    :: vcov(:),vcont(:)
 logical                          :: withouG_sub,nozero
!-----------------------------------------------------------------------------
if(present(withoutG))then
 withouG_sub=withoutG
else
 withouG_sub=.false.
end if

cond_wphys: if(present(wphys).and. .not.withouG_sub) then
 where(patchw(ixO^S))
  ^FL&wphys(ixO^S,^FL)=w(ixO^S,^FL)/mypm%sqrtdetr%Gama(ixO^S);
 end where
 cond_cov1 : if(nwv_cov/=0)then
  allocate(vcont(nwv_cov))
  Loop_inivcov : do iv=1,nwv_cov
   where(patchw(ixO^S))
    ^C&vcov(iv)%w(ixO^S,^C)=wphys(ixO^S,iwv_cov(iv)+^C);
   end where
  end do Loop_inivcov
 end if cond_cov1
 cond_cont1 : if(nwv_cont/=0)then
  allocate(vcont(nwv_cont))
  Loop_inivcont : do iv=1,nwv_cont
   where(patchw(ixO^S))
    ^C&vcont(iv)%w(ixO^S,^C)=wphys(ixO^S,iwv_cont(iv)+^C);
   end where
  end do Loop_inivcont
 end if cond_cont1
 cond_cov2: if(nwv_cov/=0) then
 Loop_idircov : do iidir=1,setgr%ndir
   nozero=.false.
   if(all((/(.not.mypm%elemsp_inv(i1_gr(iidir,jdir),&
                 i2_gr(iidir,jdir))%elpt%elm_on&
            ,jdir=1,setgr%ndir)/)))then
    Loop_setvcov0 : do iv=1,nwv_cov
     where(patchw(ixO^S))
      wphys(ixO^S,iwv_cov(iv)+iidir)=zero;
     end where
    end do Loop_setvcov0
    cycle Loop_idircov
   end if
 Loop_jdircov : do jjdir=1,setgr%ndir;idir=i1_gr(iidir,jjdir);
   jdir=i2_gr(iidir,jjdir)
   cond_on : if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on)then
     cond_isone : if(.not.mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
      cond_nozero1 : if(nozero)then
       Loop_sumv1add : do iv=1,nwv_cov
        where(patchw(ixO^S))
         wphys(ixO^S,iwv_cov(iv)+iidir)=wphys(ixO^S,iwv_cov(iv)+iidir)+&
                         vcov(iv)%w(ixO^S,jjdir)&
                         *mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S);
        end where
       end do Loop_sumv1add
      else cond_nozero1
       Loop_sumvcov1 : do iv=1,nwv_cov
        where(patchw(ixO^S))
         wphys(ixO^S,iwv_cov(iv)+iidir)=vcov(iv)%w(ixO^S,jjdir)&
                         *mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S);
        end where
       end do Loop_sumvcov1
      end if cond_nozero1
     else cond_isone
      cond_nozero2 : if(nozero)then
       Loop_sumv1bAadd : do iv=1,nwv_cov
        where(patchw(ixO^S))
         wphys(ixO^S,iwv_cov(iv)+iidir)=wphys(ixO^S,iwv_cov(iv)+iidir)&
                                    + vcov(iv)%w(ixO^S,jjdir);
        end where
       end do Loop_sumv1bAadd
      else cond_nozero2
       Loop_sumvcov1bA : do iv=1,nwv_cov
        where(patchw(ixO^S))
         wphys(ixO^S,iwv_cov(iv)+iidir)=vcov(iv)%w(ixO^S,jjdir);
        end where
       end do Loop_sumvcov1bA
     end if cond_nozero2
     nozero=.true.
    end if cond_isone
   end if cond_on
 end do Loop_jdircov
 end do Loop_idircov
 if(present(pscal_cov)) then
 Loop_getvcov2 : do iv=1,nwv_cov
  where(patchw(ixO^S))
   pscal_cov(iv)%norm2(ixO^S)=&
          {^C&wphys(ixO^S,iwv_cov(iv)+^C)*vcov(iv)%w(ixO^S,^C)+};
  end where
 end do Loop_getvcov2
 end if
end if cond_cov2
else cond_wphys

if(.not.withouG_sub)then
 where(patchw(ixO^S))
  ^FL&w(ixO^S,^FL)=w(ixO^S,^FL)/mypm%sqrtdetr%Gama(ixO^S);
 end where
end if
  Loop_inivcovB : do iv=1,nwv_cov
   where(patchw(ixO^S))
    ^C&vcov(iv)%w(ixO^S,^C)=w(ixO^S,iwv_cov(iv)+^C);
   end where
  end do Loop_inivcovB
  Loop_inivcontB : do iv=1,nwv_cont
   where(patchw(ixO^S))
    ^C&vcont(iv)%w(ixO^S,^C)=w(ixO^S,iwv_cont(iv)+^C);
   end where
  end do Loop_inivcontB
 cond_cov2B : if(nwv_cov/=0) then
 nozero=.false.
 Loop_idircovB : do iidir=1,setgr%ndir
  nozero=.false.
   if(all((/(.not.mypm%elemsp_inv(i1_gr(iidir,jdir),&
                 i2_gr(iidir,jdir))%elpt%elm_on&
            ,jdir=1,setgr%ndir)/)))then
    Loop_setvcov0B : do iv=1,nwv_cov
     where(patchw(ixO^S))
      w(ixO^S,iwv_cov(iv)+iidir)=zero;
     end where
    end do Loop_setvcov0B
    cycle Loop_idircovB
   end if
 Loop_jdircovB : do jjdir=1,setgr%ndir;idir=i1_gr(iidir,jjdir);
   jdir=i2_gr(iidir,jjdir)
   cond_oncov1B : if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on)then
     cond_isonecov1B : if(.not.mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
      cond_nozerocov1B : if(nozero)then
       Loop_sumvcov1B : do iv=1,nwv_cov
        where(patchw(ixO^S))
          w(ixO^S,iwv_cov(iv)+iidir)=w(ixO^S,iwv_cov(iv)+iidir)+&
                         vcov(iv)%w(ixO^S,jjdir)&
                         *mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S);
        end where
       end do Loop_sumvcov1B
      else cond_nozerocov1B
       Loop_sumvcov2B : do iv=1,nwv_cov
        where(patchw(ixO^S))
          w(ixO^S,iwv_cov(iv)+iidir)=vcov(iv)%w(ixO^S,jjdir)&
                         *mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S);
        end where
       end do Loop_sumvcov2B
      end if cond_nozerocov1B 
     else cond_isonecov1B
      cond_nozerocov1bB : if(nozero)then
       Loop_sumvcov1bB : do iv=1,nwv_cov
        where(patchw(ixO^S))
         w(ixO^S,iwv_cov(iv)+iidir)=w(ixO^S,iwv_cov(iv)+iidir)+&
                                    vcov(iv)%w(ixO^S,jjdir);
        end where
       end do Loop_sumvcov1bB
      else cond_nozerocov1bB
       Loop_sumvcov2bB : do iv=1,nwv_cov
        where(patchw(ixO^S))
         w(ixO^S,iw_vector(iwv_cov(iv))+iidir)=vcov(iv)%w(ixO^S,jjdir);
        end where
       end do Loop_sumvcov2bB
      end if cond_nozerocov1bB
     end if cond_isonecov1B
     nozero=.true.
    end if cond_oncov1B
 end do Loop_jdircovB
end do Loop_idircovB
 if(present(pscal_cov))then
 Loop_getvcov2B : do iv=1,nwv_cov
  where(patchw(ixO^S))
   pscal_cov(iv)%norm2(ixO^S)={^C&w(ixO^S,iwv_cov(iv)+^C)*vcov(iv)%w(ixO^S,^C)+};
  end where
 end do Loop_getvcov2B
 end if
 end if cond_cov2B
end if cond_wphys

cond_cont2 : if(nwv_cont/=0.and.present(pscal_cont)) then
 nozero=.false.
 Loop_idircont : do iidir=1,setgr%ndir
   if(all((/(.not.mypm%elemsp(i1_gr(iidir,jdir),&
                 i2_gr(iidir,jdir))%elpt%elm_on&
            ,jdir=1,setgr%ndir)/)))cycle
 Loop_jdircont : do jjdir=1,setgr%ndir;idir=i1_gr(iidir,jjdir);
   jdir=i2_gr(iidir,jjdir)
   cond_on_cont : if(mypm%elemsp(idir,jdir)%elpt%elm_on)then
     cond_isone_cont :if(.not.mypm%elemsp(idir,jdir)%elpt%elm_isone)then
      cond_nozeroB1 : if(nozero)then
       Loop_nwv1a : do iv=1,nwv_cont
        where(patchw(ixO^S))
         pscal_cont(iv)%norm2(ixO^S)=pscal_cont(iv)%norm2(ixO^S)+&
                           vcont(iv)%w(ixO^S,jjdir)&
                           *vcont(iv)%w(ixO^S,iidir)*&
                           mypm%elemsp(idir,jdir)%elpt%wg(ixO^S)
        end where
       end do Loop_nwv1a
      else cond_nozeroB1
       Loop_nwv1b : do iv=1,nwv_cont
        where(patchw(ixO^S))
          pscal_cont(iv)%norm2(ixO^S)=vcont(iv)%w(ixO^S,jjdir)&
                  *vcont(iv)%w(ixO^S,iidir)*&
                  mypm%elemsp(idir,jdir)%elpt%wg(ixO^S)
        end where
       end do Loop_nwv1b
      end if cond_nozeroB1

     else cond_isone_cont
      cond_nozeroB2 : if(nozero)then
       Loop_nwv2a : do iv=1,nwv_cont
        where(patchw(ixO^S))
         pscal_cont(iv)%norm2(ixO^S)=pscal_cont(iv)%norm2(ixO^S)+&
                           vcont(iv)%w(ixO^S,jjdir)&
                           *vcont(iv)%w(ixO^S,iidir)
        end where
       end do Loop_nwv2a
      else cond_nozeroB2
     Loop_nwv2b : do iv=1,nwv_cont
        where(patchw(ixO^S))
          pscal_cont(iv)%norm2(ixO^S)=vcont(iv)%w(ixO^S,jjdir)&
                  *vcont(iv)%w(ixO^S,iidir)
        end where
       end do Loop_nwv2b
      end if cond_nozeroB2
     end if cond_isone_cont
    end if cond_on_cont
  end do Loop_jdircont
 end do Loop_idircont
end if cond_cont2

end subroutine  usemetric_getaux_ready
!=============================================================================
subroutine usemetric_primitive(ixI^L,ixO^L,patchw,subname,w,nwv,iwv,it_cov)
! w should be covariant            
   
 integer, intent(in)             :: ixI^L,ixO^L
 logical, intent(in)             :: patchw(ixG^T)
 character(len=*), intent(in)    :: subname
 double precision, intent(inout) :: w(ixI^S,1:nw)
 integer, intent(in)             :: nwv,iwv(nwv)
 logical, optional               :: it_cov
 ! .. local variables ..
 integer              :: idir,jdir,kdir,iidir,jjdir,kkdir,iv
 type(vector_phys)    :: vcov(nwv)
 logical              :: nozero
!-----------------------------------------------------------------------------
if(.not.present(it_cov))then
 Loop_iniv : do iv=1,nwv
  where(patchw(ixO^S))
   ^C&vcov(iv)%w(ixO^S,^C)=w(ixO^S,iwv(iv)+^C);
  end where
 end do Loop_iniv
 Loop_idir : do iidir=1,setgr%ndir
   nozero=.false.
   if(all((/(.not.mypm%elemsp_inv(i1_gr(iidir,jdir),&
                 i2_gr(iidir,jdir))%elpt%elm_on&
            ,jdir=1,setgr%ndir)/)))then
    Loop_setv0 : do iv=1,nwv
     where(patchw(ixO^S))
      w(ixO^S,iwv(iv)+iidir)=zero;
     end where
    end do Loop_setv0
    cycle Loop_idir
   end if
 Loop_jdir : do jjdir=1,setgr%ndir;idir=i1_gr(iidir,jjdir);
   jdir=i2_gr(iidir,jjdir)
   if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on)then
     if(.not.mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
      if(nozero)then
       Loop_sumv1add : do iv=1,nwv
        where(patchw(ixO^S))
         w(ixO^S,iwv(iv)+iidir)= w(ixO^S,iwv(iv)+iidir)+vcov(iv)%w(ixO^S,jjdir)&
                          *mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S);
        end where
       end do Loop_sumv1add
      else
       Loop_sumv1 : do iv=1,nwv
        where(patchw(ixO^S))
         w(ixO^S,iwv(iv)+iidir)=vcov(iv)%w(ixO^S,jjdir)&
                          *mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S);
        end where
       end do Loop_sumv1
      end if 
     else
      if(nozero)then
       Loop_sumv1badd : do iv=1,nwv
        where(patchw(ixO^S))
         w(ixO^S,iwv(iv)+iidir)=w(ixO^S,iwv(iv)+iidir)+vcov(iv)%w(ixO^S,jjdir);
        end where
       end do Loop_sumv1badd
      else                                   
       Loop_sumv1b : do iv=1,nwv
        where(patchw(ixO^S))
         w(ixO^S,iwv(iv)+iidir)=vcov(iv)%w(ixO^S,jjdir);
        end where
       end do Loop_sumv1b
      end if
     end if
     nozero=.true.
    end if
 end do Loop_jdir
 end do Loop_idir
end if

end subroutine usemetric_primitive
!=============================================================================
subroutine  usemetric_w_withoutspacedetr(ixI^L,ixO^L,patchw,w,wphys,iwneeded)
 include 'amrvacdef.f'
 integer, intent(in)                    :: ixI^L,ixO^L
 logical, intent(in)                    :: patchw(ixG^T)
 double precision, intent(inout)        :: w(ixI^S,1:nw)
 double precision, optional, intent(out):: wphys(ixG^T,1:nw)
 integer, optional                      :: iwneeded(:)
 ! .. local variables ..
 integer                                 :: iw,iiw,nwneeded,nwneeded_flux
!-----------------------------------------------------------------------------
 if(present(wphys)) then
 cond_iwnd1 : if(present(iwneeded))then
  nwneeded=size(iwneeded)
  nwneeded_flux=nwneeded-count(iwneeded>nwflux)
  Loop_iwf1 : do iiw=1,nwneeded_flux;iw=iwneeded(iiw)
   where(patchw(ixO^S))
    wphys(ixO^S,iw)=w(ixO^S,iw)/ mypm%sqrtdetr%Gama(ixO^S);
   end where
  end do Loop_iwf1

  cond_aux1 : if(nwaux>0.and.any(iwneeded>nwflux))then
   Loop_iwaux : do iiw=nwneeded_flux+1,nwneeded;iw=iwneeded(iiw)
    where(patchw(ixO^S))
     wphys(ixO^S,iw)=w(ixO^S,iw)
    end where
   end do Loop_iwaux
  end if cond_aux1
 else cond_iwnd1
   where(patchw(ixO^S))
    {^FL&wphys(ixO^S,^FL)=w(ixO^S,^FL)/mypm%sqrtdetr%Gama(ixO^S);\}
   end where
  if(nwaux>0)then
   {^FLAUX&wphys(ixO^S,^FLAUX)=w(ixO^S,^FLAUX);}
  end if
 end if cond_iwnd1
else
 cond_iwnd2 : if(present(iwneeded))then
  nwneeded=size(iwneeded)
  nwneeded_flux=nwneeded-count(iwneeded>nwflux)
  Loop_iwf : do iiw=1,nwneeded_flux;iw=iwneeded(iw)
   where(patchw(ixO^S))
    w(ixO^S,iw)=w(ixO^S,iw)/mypm%sqrtdetr%Gama(ixO^S);
   end where
  end do Loop_iwf

 else cond_iwnd2
   where(patchw(ixO^S))
    {^FL&w(ixO^S,^FL)=w(ixO^S,^FL)/mypm%sqrtdetr%Gama(ixO^S);\}
   end where
 end if cond_iwnd2

end if
end subroutine  usemetric_w_withoutspacedetr
!===========================================================================
subroutine  usemetric_vec_withoutspacedetr(ixI^L,ixO^L,patchw,vec,vecphys)
 include 'amrvacdef.f'
 integer, intent(in)                    :: ixI^L,ixO^L
 logical, intent(in)                    :: patchw(ixG^T)
 double precision, intent(inout)        :: vec(ixI^S,1:ndir)
 double precision, optional, intent(out):: vecphys(ixG^T,1:ndir)
!-----------------------------------------------------------------------------
if(present(vecphys)) then
   where(patchw(ixO^S))
    {^C&vecphys(ixO^S,^C)=vec(ixO^S,^C)/mypm%sqrtdetr%Gama(ixO^S);}
   end where
else
   where(patchw(ixO^S))
    {^C&vec(ixO^S,^C)=vec(ixO^S,^C)/mypm%sqrtdetr%Gama(ixO^S);}
   end where
end if
end subroutine  usemetric_vec_withoutspacedetr
!=============================================================================
!=============================================================================
subroutine  usemetric_covtobasis(ixI^L,ixO^L,wvcov,wvbasis)
 include 'amrvacdef.f'
 integer, intent(in)              :: ixI^L,ixO^L
 double precision, intent(in)     :: wvcov(ixI^S,1:ndir)
 double precision, intent(out)    :: wvbasis(ixG^T,1:ndir)
 ! .. local variables ..
 integer                          :: idir,idir2,jdir,iidir,iidir2,jjdir
!-----------------------------------------------------------------------------
 Loop_jdir : do jjdir=1,setgr%ndir
  Loop_idir : do iidir=1,setgr%ndir;idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
   if(mypm%Qsp(idir,jdir)%elpt%elm_on)then
     if(mypm%Qsp(idir,jdir)%elpt%elm_isone)then
      wvbasis(ixO^S,jjdir)=wvcov(ixO^S,iidir)
     else
      wvbasis(ixO^S,jjdir)=mypm%Qsp(idir,jdir)%elpt%wg(ixO^S)*wvcov(ixO^S,iidir)
     end if
     exit Loop_idir
   end if
  end do Loop_idir

  if(idir<setgr%ndir)then
   Loop_idir2 : do iidir2=idir+1,setgr%ndir;idir2=i1_gr(iidir2,jjdir);
                                            jdir=i2_gr(iidir2,jjdir)
    if(mypm%Qsp(idir2,jdir)%elpt%elm_on)then
      if(mypm%Qsp(idir2,jdir)%elpt%elm_isone)then
       wvbasis(ixO^S,jjdir)=wvbasis(ixO^S,jjdir)+wvcov(ixO^S,iidir2)
      else
       wvbasis(ixO^S,jjdir)=wvbasis(ixO^S,jjdir)+&
                          mypm%Qsp(idir2,jdir)%elpt%wg(ixO^S)*wvcov(ixO^S,iidir2)
      end if
    end if
   end do Loop_idir2
  end if
 end do Loop_jdir
end subroutine  usemetric_covtobasis
!=============================================================================
subroutine  usemetric_conttobasis(ixI^L,ixO^L,wvcont,wvbasis)
 include 'amrvacdef.f'
 integer, intent(in)              :: ixI^L,ixO^L
 double precision, intent(in)     :: wvcont(ixI^S,1:ndir)
 double precision, intent(out)    :: wvbasis(ixG^T,1:ndir)
 ! .. local variables ..
 integer                          :: idir,idir2,jdir,iidir,iidir2,jjdir
!-----------------------------------------------------------------------------
 Loop_jdir : do jjdir=1,setgr%ndir
  Loop_idir : do iidir=1,setgr%ndir;idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
   if(mypm%Qsp_inv(idir,jdir)%elpt%elm_on)then
     if(mypm%Qsp_inv(idir,jdir)%elpt%elm_isone)then
      wvbasis(ixO^S,jjdir)=wvcont(ixO^S,iidir)
     else
      wvbasis(ixO^S,jjdir)=mypm%Qsp_inv(idir,jdir)%elpt%wg(ixO^S)&
                          *wvcont(ixO^S,iidir)
     end if
     exit Loop_idir
   end if
  end do Loop_idir

  if(idir<setgr%ndir)then
   Loop_idir2 : do iidir2=idir+1,setgr%ndir;idir2=i1_gr(iidir2,jjdir);
                                            jdir=i2_gr(iidir2,jjdir)
    if(mypm%Qsp_inv(idir2,jdir)%elpt%elm_on)then
      if(mypm%Qsp_inv(idir2,jdir)%elpt%elm_isone)then
       wvbasis(ixO^S,jjdir)=wvbasis(ixO^S,jjdir)+wvcont(ixO^S,iidir2)
      else
       wvbasis(ixO^S,jjdir)=wvbasis(ixO^S,jjdir)+&
                   mypm%Qsp_inv(idir2,jdir)%elpt%wg(ixO^S)*wvcont(ixO^S,iidir2)
      end if
    end if
   end do Loop_idir2
  end if
 end do Loop_jdir
end subroutine  usemetric_conttobasis
!=============================================================================
subroutine  usemetric_basistocont(ixI^L,ixO^L,wvbasis,wvcont)
 include 'amrvacdef.f'
 integer, intent(in)              :: ixI^L,ixO^L
 double precision, intent(in)     :: wvbasis(ixI^S,1:ndir)
 double precision, intent(out)    :: wvcont(ixG^T,1:ndir)
 ! .. local variables ..
 integer                          :: idir,idir2,jdir,iidir,iidir2,jjdir
!-----------------------------------------------------------------------------
 Loop_jdir : do jjdir=1,setgr%ndir
  Loop_idir : do iidir=1,setgr%ndir;idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
   if(mypm%Qsp(idir,jdir)%elpt%elm_on)then
     if(mypm%Qsp(idir,jdir)%elpt%elm_isone)then
      wvcont(ixO^S,jjdir)=wvbasis(ixO^S,iidir)
     else
      wvcont(ixO^S,jjdir)=mypm%Qsp(idir,jdir)%elpt%wg(ixO^S)&
                         *wvbasis(ixO^S,iidir)
     end if
     exit Loop_idir
   end if
  end do Loop_idir

  if(idir<setgr%ndir)then
   Loop_idir2 : do iidir2=idir+1,setgr%ndir;idir2=i1_gr(iidir2,jjdir);&
                                            jdir=i2_gr(iidir2,jjdir)
    if(mypm%Qsp(idir2,jdir)%elpt%elm_on)then
      if(mypm%Qsp(idir2,jdir)%elpt%elm_isone)then
       wvcont(ixO^S,jjdir)=wvcont(ixO^S,jjdir)+wvbasis(ixO^S,iidir2)
      else
       wvcont(ixO^S,jjdir)=wvcont(ixO^S,jjdir)+&
            mypm%Qsp(idir2,jdir)%elpt%wg(ixO^S)*wvbasis(ixO^S,iidir2)
      end if
    end if
   end do Loop_idir2
  end if
 end do Loop_jdir
end subroutine  usemetric_basistocont
!=============================================================================
subroutine  usemetric_basistocov(ixI^L,ixO^L,wvbasis,wvcov)
 include 'amrvacdef.f'
 integer, intent(in)              :: ixI^L,ixO^L
 double precision, intent(in)     :: wvbasis(ixI^S,1:ndir)
 double precision, intent(out)    :: wvcov(ixG^T,1:ndir)
 ! .. local variables ..
 integer                          :: idir,idir2,jdir,iidir,iidir2,jjdir
!-----------------------------------------------------------------------------
 Loop_jdir : do jjdir=1,setgr%ndir
  Loop_idir : do iidir=1,setgr%ndir;idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
   if(mypm%Qsp_inv(idir,jdir)%elpt%elm_on)then
     if(mypm%Qsp_inv(idir,jdir)%elpt%elm_isone)then
      wvcov(ixO^S,jjdir)=wvbasis(ixO^S,iidir)
     else
      wvcov(ixO^S,jjdir)=mypm%Qsp_inv(idir,jjdir)%elpt%wg(ixO^S)&
                        *wvbasis(ixO^S,iidir)
     end if
     exit Loop_idir
   end if
  end do Loop_idir

  if(idir<setgr%ndir)then
   Loop_idir2 : do iidir2=idir+1,setgr%ndir;idir2=i1_gr(iidir2,jjdir);
                                            jdir=i2_gr(iidir2,jjdir)
    if(mypm%Qsp_inv(idir2,jdir)%elpt%elm_on)then
      if(mypm%Qsp_inv(idir2,jdir)%elpt%elm_isone)then
       wvcov(ixO^S,jjdir)=wvcov(ixO^S,jjdir)+wvbasis(ixO^S,iidir2)
      else
       wvcov(ixO^S,jjdir)=wvcov(ixO^S,jjdir)+&
            mypm%Qsp_inv(idir2,jdir)%elpt%wg(ixO^S)*wvbasis(ixO^S,iidir2)
      end if
    end if
   end do Loop_idir2
  end if
 end do Loop_jdir
end subroutine  usemetric_basistocov
!=============================================================================
{subroutine  usemetric_bcvcot_switchG^D(ixI^L,ixO^L,ixB^L,lshift,niw,iws,wvct)
 include 'amrvacdef.f'
 integer, intent(in)                     :: ixI^L,ixO^L,ixB^L,niw,iws(nw)
 logical, intent(in)                     :: lshift(1:nw)
 double precision, intent(inout)         :: wvct(ixI^S,1:nw)
 !.. local ..
 double precision, dimension(ixG^T,1:nw) :: wvbasis_symm, wvbasis_cont
 integer                                 :: idir,idir2,jdir,iidir,iidir2,jjdir
 integer                                 :: iiw,iw,iv0_,niv
 logical                                 :: ivect(1:nw),ivs(nvector)
!----------------------------------------------------------
 niv=0
 ivect=.false.
 Loop_iw_0 : do iw=1,niw
  iv0_=iw-idir_iw(iw)
  if(ivector_iw(iw)>0.and..not.ivect(iv0_))then
    niv=niv+1
    ivs(niv)=iv0_
    ivect(iv0_)=.true.
  end if
 end do Loop_iw_0
 Loop_jdir : do jjdir=1,setgr%ndir
   Loop_idir : do iidir=1,setgr%ndir;idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
    if(mypm%Qsp_inv(idir,jdir)%elpt%elm_on)then
     if(mypm%Qsp_inv(idir,jdir)%elpt%elm_isone)then
      Loop_vect1_1 :do iw=1,niv;iv0_=ivs(iw)
       if(any(lshift(iv0_+1:iv0_+^NC)))then
        wvbasis_cont(ixO^S,iv0_+jjdir)=wvct(ixO^S,iv0_+iidir)
       end if
       if(any(.not.lshift(iv0_+1:iv0_+^NC)))then
        wvbasis_symm(ixO^S,iv0_+jjdir)=wvct(ixO^S,iv0_+iidir)
       end if

      end do Loop_vect1_1
     else
      Loop_vect1_2 : do iw=1,niv;iv0_=ivs(iw)
       if(any(lshift(iv0_+1:iv0_+^NC)))then
        wvbasis_cont(ixO^S,iv0_+jjdir)=mypm%Qsp_inv(idir,jdir)%elpt%wg(ixB^S)&
                                *wvct(ixO^S,iv0_+iidir)
       end if
       if(any(.not.lshift(iv0_+1:iv0_+^NC)))then
        wvbasis_symm(ixO^S,iv0_+jjdir)=wvct(ixO^S,iv0_+iidir)*&
                   mypm%Qsp_inv(idir,jdir)%elpt%wg(ixBmax^D:ixBmin^D:-1^D%ixB^S)
       end if
      end do Loop_vect1_2
     end if
     exit Loop_idir
    end if
   end do Loop_idir
  if(idir<setgr%ndir)then
    Loop_idir2 : do iidir2=iidir+1,setgr%ndir;idir2=i1_gr(iidir2,jjdir);
                                              jdir=i2_gr(iidir2,jjdir)
     if(mypm%Qsp_inv(idir2,jdir)%elpt%elm_on)then
      if(mypm%Qsp_inv(idir2,jdir)%elpt%elm_isone)then
       Loop_vect2_1 : do iw=1,niv;iv0_=ivs(iw)
        if(any(lshift(iv0_+1:iv0_+^NC)))then
         wvbasis_cont(ixO^S,iv0_+jjdir)=wvbasis_cont(ixO^S,jjdir)+&
                                       wvct(ixO^S,iv0_+iidir2)
        end if
       if(any(.not.lshift(iv0_+1:iv0_+^NC)))then
         wvbasis_symm(ixO^S,iv0_+jjdir)=wvbasis_symm(ixO^S,jjdir)+&
                                       wvct(ixO^S,iv0_+iidir2)
       end if
       end do Loop_vect2_1
      else
       Loop_vect2_2 : do iw=1,niv;iv0_=ivs(iw)
        if(any(lshift(iv0_+1:iv0_+^NC)))then
         wvbasis_cont(ixO^S,iv0_+jjdir)=wvbasis_cont(ixO^S,iv0_+jjdir)+&
                mypm%Qsp_inv(idir2,jdir)%elpt%wg(ixB^S)*wvct(ixO^S,iv0_+iidir2)
        end if
        if(any(.not.lshift(iv0_+1:iv0_+^NC)))then
         wvbasis_symm(ixO^S,iv0_+jjdir)=wvbasis_symm(ixO^S,iv0_+jjdir)+&
               mypm%Qsp_inv(idir2,jdir)%elpt%wg(ixBmax^D:ixBmin^D:-1^D%ixB^S)*&
                    wvct(ixO^S,iv0_+iidir2)
        end if
       end do Loop_vect2_2
      end if
     end if
    end do Loop_idir2
   end if
 end do Loop_jdir
 Loop_iw : do iiw=1,niw;iw=iws(iiw);jjdir=idir_iw(iw);iv0_=iw-jjdir
   Loop_idir2_1 : do iidir=1,setgr%ndir;idir=i1_gr(iidir,jjdir);
                                        jdir=i2_gr(iidir,jjdir)
    if(mypm%Qsp(idir,jdir)%elpt%elm_on)then
     if(mypm%Qsp(idir,jdir)%elpt%elm_isone)then
      if(lshift(iw))then
       wvct(ixO^S,iw)=wvbasis_cont(ixO^S,iv0_+iidir)
      else
       wvct(ixO^S,iw)=wvbasis_symm(ixO^S,iv0_+iidir)
      end if
     else
      if(lshift(iw))then
       wvct(ixO^S,iw)=mypm%Qsp(idir,jdir)%elpt%wg(ixO^S)*&
                       wvbasis_cont(ixO^S,iv0_+iidir)
      else
       wvct(ixO^S,iw)=mypm%Qsp(idir,jdir)%elpt%wg(ixO^S)*&
                       wvbasis_symm(ixO^S,iv0_+iidir)
      end if
     end if
     exit Loop_idir2_1
    end if
   end do Loop_idir2_1

   if(idir<setgr%ndir)then
    Loop_idir2_2 : do iidir2=idir+1,setgr%ndir;idir2=i1_gr(iidir2,jjdir);
                                               jdir=i2_gr(iidir2,jjdir)
     if(mypm%Qsp(idir2,jdir)%elpt%elm_on)then
      if(mypm%Qsp(idir2,jdir)%elpt%elm_isone)then
       if(lshift(iw))then
        wvct(ixO^S,iw)=wvct(ixO^S,iw)+wvbasis_cont(ixO^S,iv0_+iidir2)
       else
        wvct(ixO^S,iw)=wvct(ixO^S,iw)+wvbasis_symm(ixO^S,iv0_+iidir2)
       end if

      else
       if(lshift(iw))then
        wvct(ixO^S,iw)=wvct(ixO^S,iw)+&
           mypm%Qsp(idir2,jdir)%elpt%wg(ixO^S)*wvbasis_cont(ixO^S,iv0_+iidir2)
       else
        wvct(ixO^S,iw)=wvct(ixO^S,iw)+&
          mypm%Qsp(idir2,jdir)%elpt%wg(ixO^S)*wvbasis_symm(ixO^S,iv0_+iidir2)
       end if
      end if
     end if
    end do Loop_idir2_2
   end if
 end do Loop_iw
end subroutine  usemetric_bcvcot_switchG^D
\}
!=============================================================================
{subroutine  usemetric_bcvcov_switchG^D(ixI^L,ixO^L,ixB^L,lshift,niw,iws,wvcov)
 include 'amrvacdef.f'
 integer, intent(in)                     :: ixI^L,ixO^L,ixB^L,niw,iws(nw)
 logical, intent(in)                     :: lshift(1:nw)
 double precision, intent(inout)         :: wvcov(ixI^S,1:nw)
 !.. local ..
 double precision, dimension(ixG^T,1:nw) :: wvbasis_symm, wvbasis_cont
 integer                                 :: idir,idir2,jdir,iidir,iidir2,jjdir
 integer                                 ::iiw,iw,iv0_,niv
 logical                                 :: ivect(1:nw),ivs(nvector)
!----------------------------------------------------------
 niv=0
 ivect=.false.
 do iw=1,niw
  iv0_=iw-idir_iw(iw)
  if(ivector_iw(iw)>0.and..not.ivect(iv0_))then
    niv=niv+1
    ivs(niv)=iv0_
    ivect(iv0_)=.true.
  end if
 end do
 Loop_jdir : do jjdir=1,setgr%ndir
   Loop_idir : do iidir=1,setgr%ndir;idir=i1_gr(iidir,jjdir)
                                     jdir=i2_gr(iidir,jjdir)
    if(mypm%Qsp(idir,jdir)%elpt%elm_on)then
     if(mypm%Qsp(idir,jdir)%elpt%elm_isone)then
      Loop_vect1_1 :do iw=1,niv;iv0_=ivs(iw)
       if(any(lshift(iv0_+1:iv0_+^NC)))then
        wvbasis_cont(ixO^S,iv0_+jjdir)=wvcov(ixO^S,iv0_+iidir)
       end if
       if(any(.not.lshift(iv0_+1:iv0_+^NC)))then
        wvbasis_symm(ixO^S,iv0_+jjdir)=wvcov(ixO^S,iv0_+iidir)
       end if

      end do Loop_vect1_1
     else
      Loop_vect1_2 : do iw=1,niv;iv0_=ivs(iw)
       if(any(lshift(iv0_+1:iv0_+^NC)))then
        wvbasis_cont(ixO^S,iv0_+jjdir)=mypm%Qsp(idir,jdir)%elpt%wg(ixB^S)&
                                *wvcov(ixO^S,iv0_+iidir)
       end if
       if(any(.not.lshift(iv0_+1:iv0_+^NC)))then
        wvbasis_symm(ixO^S,iv0_+jjdir)=wvcov(ixO^S,iv0_+iidir)*&
                   mypm%Qsp(idir,jdir)%elpt%wg(ixBmax^D:ixBmin^D:-1^D%ixB^S)
       end if
      end do Loop_vect1_2
     end if
     exit Loop_idir
    end if
   end do Loop_idir
  if(iidir<setgr%ndir)then
    Loop_idir2 : do iidir2=iidir+1,setgr%ndir;idir2=i1_gr(iidir2,jjdir);
                                              jdir=i2_gr(iidir2,jjdir)
     if(mypm%Qsp(idir2,jdir)%elpt%elm_on)then
      if(mypm%Qsp(idir2,jdir)%elpt%elm_isone)then
       Loop_vect2_1 : do iw=1,niv;iv0_=ivs(iw)
        if(any(lshift(iv0_+1:iv0_+^NC)))then
         wvbasis_cont(ixO^S,iv0_+jjdir)=wvbasis_cont(ixO^S,jjdir)+&
                                       wvcov(ixO^S,iv0_+iidir2)
        end if
       if(any(.not.lshift(iv0_+1:iv0_+^NC)))then
         wvbasis_symm(ixO^S,iv0_+jjdir)=wvbasis_symm(ixO^S,jdir)+&
                                       wvcov(ixO^S,iv0_+iidir2)
       end if
       end do Loop_vect2_1
      else
       Loop_vect2_2 : do iw=1,niv;iv0_=ivs(iw)
        if(any(lshift(iv0_+1:iv0_+^NC)))then
         wvbasis_cont(ixO^S,iv0_+jjdir)=wvbasis_cont(ixO^S,iv0_+jjdir)+&
                     mypm%Qsp(idir2,jdir)%elpt%wg(ixB^S)*wvcov(ixO^S,iv0_+iidir2)
        end if
        if(any(.not.lshift(iv0_+1:iv0_+^NC)))then
         wvbasis_symm(ixO^S,iv0_+jjdir)=wvbasis_symm(ixO^S,iv0_+jjdir)+&
                    mypm%Qsp(idir2,jdir)%elpt%wg(ixBmax^D:ixBmin^D:-1^D%ixB^S)*&
                    wvcov(ixO^S,iv0_+iidir2)
        end if
       end do Loop_vect2_2
      end if
     end if
    end do Loop_idir2
   end if
 end do Loop_jdir
 Loop_iw : do iiw=1,niw;iw=iws(iiw);jjdir=idir_iw(iw);iv0_=iw-jjdir
   Loop_idir2_1 : do iidir=1,setgr%ndir;idir=i1_gr(iidir,jjdir)
                                        jdir=i2_gr(iidir,jjdir)
    if(mypm%Qsp_inv(idir,jdir)%elpt%elm_on)then
     if(mypm%Qsp_inv(idir,jdir)%elpt%elm_isone)then
      if(lshift(iw))then
       wvcov(ixO^S,iw)=wvbasis_cont(ixO^S,iv0_+iidir)
      else
       wvcov(ixO^S,iw)=wvbasis_symm(ixO^S,iv0_+iidir)
      end if
     else
      if(lshift(iw))then
       wvcov(ixO^S,iw)=mypm%Qsp_inv(idir,jdir)%elpt%wg(ixO^S)*&
                       wvbasis_cont(ixO^S,iv0_+iidir)
      else
       wvcov(ixO^S,iw)=mypm%Qsp_inv(idir,jdir)%elpt%wg(ixO^S)*&
                       wvbasis_symm(ixO^S,iv0_+iidir)
      end if
     end if
     exit Loop_idir2_1
    end if
   end do Loop_idir2_1

   if(idir<setgr%ndir)then
    Loop_idir2_2 : do iidir2=idir+1,setgr%ndir;idir2=i1_gr(iidir2,jdir);
                                               jdir=i2_gr(iidir2,jdir)
     if(mypm%Qsp_inv(idir2,jdir)%elpt%elm_on)then
      if(mypm%Qsp_inv(idir2,jdir)%elpt%elm_isone)then
       if(lshift(iw))then
        wvcov(ixO^S,iw)=wvcov(ixO^S,iw)+wvbasis_cont(ixO^S,iv0_+iidir2)
       else
        wvcov(ixO^S,iw)=wvcov(ixO^S,iw)+wvbasis_symm(ixO^S,iv0_+iidir2)
       end if

      else
       if(lshift(iw))then
        wvcov(ixO^S,iw)=wvcov(ixO^S,iw)+&
           mypm%Qsp_inv(idir2,jdir)%elpt%wg(ixO^S)*wvbasis_cont(ixO^S,iv0_+iidir2)
       else
        wvcov(ixO^S,iw)=wvcov(ixO^S,iw)+&
           mypm%Qsp_inv(idir2,jdir)%elpt%wg(ixO^S)*wvbasis_symm(ixO^S,iv0_+iidir2)
       end if
      end if
     end if
    end do Loop_idir2_2
   end if
 end do Loop_iw
end subroutine  usemetric_bcvcov_switchG^D
\}

!=============================================================================
{subroutine  usemetric_bc_switchG^D(ixI^L,ixO^L,ixB^L,lshift,niw,iws,w)
 include 'amrvacdef.f'
 integer, intent(in)              :: ixI^L,ixO^L,ixB^L,niw,iws(nw)
 logical, intent(in)              :: lshift(nw)
 double precision, intent(inout)  :: w(ixI^S,1:nw)
 ! .. local ..
 integer                          :: iiw,iw
!----------------------------------------------------------
 do iiw=1,niw;iw=iws(iiw)
  if(lshift(iw))then
   w(ixO^S,iw)=w(ixO^S,iw)*(mypm%sqrtdetr%Gama(ixO^S)/mypm%sqrtdetr%Gama(ixB^S))
  else
   w(ixO^S,iw)=w(ixO^S,iw)*(mypm%sqrtdetr%Gama(ixO^S)&
                 /mypm%sqrtdetr%Gama(ixBmax^D:ixBmin^D:-1^D%ixB^S))
  end if
 end do
end subroutine  usemetric_bc_switchG^D
\}

!=============================================================================
subroutine  usemetric_multiby_sqrtspacedetr(ixO^L,patchw,w,wphys,powerG)
 include 'amrvacdef.f'
 integer, intent(in)                        :: ixO^L
 logical, intent(in)                        :: patchw(ixG^T)
 double precision, intent(inout)            :: w(ixO^S)
 double precision, optional, intent(out)    :: wphys(ixG^T)
 double precision, optional                 :: powerG
!-----------------------------------------------------------------------------
 if(present(powerG)) then
  if(present(wphys))then
   where(patchw(ixO^S))
    wphys(ixO^S)=w(ixO^S)*mypm%sqrtdetr%Gama(ixO^S)**powerG
   end where
  else
   where(patchw(ixO^S))
    w(ixO^S)=w(ixO^S)*mypm%sqrtdetr%Gama(ixO^S)**powerG
   end where
  end if
 else
  if(present(wphys))then
   where(patchw(ixO^S))
    wphys(ixO^S)=w(ixO^S)*mypm%sqrtdetr%Gama(ixO^S)
   end where
  else
   where(patchw(ixO^S))
    w(ixO^S)=w(ixO^S)*mypm%sqrtdetr%Gama(ixO^S)
   end where
  end if
 end if
end subroutine  usemetric_multiby_sqrtspacedetr
!=============================================================================
subroutine  usemetric_devideby_sqrtspacedetr(ixO^L,patchw,w,wphys,powerG)
 include 'amrvacdef.f'
 integer, intent(in)                       :: ixO^L
 logical, intent(in)                       :: patchw(ixG^T)
 double precision, intent(inout)           :: w(ixO^S)
 double precision, optional,intent(out)    :: wphys(ixG^T)
 double precision, optional                :: powerG
!-----------------------------------------------------------------------------
if(present(wphys))then
 if(present(powerG)) then
  where(patchw(ixO^S))
   wphys(ixO^S)=w(ixO^S)/mypm%sqrtdetr%Gama(ixO^S)**powerG
 end where
 else
  where(patchw(ixO^S))
    wphys(ixO^S)=w(ixO^S)/mypm%sqrtdetr%Gama(ixO^S)
  end where
 end if
else
 if(present(powerG)) then
  where(patchw(ixO^S))
   w(ixO^S)=w(ixO^S)/mypm%sqrtdetr%Gama(ixO^S)**powerG
 end where
 else
  where(patchw(ixO^S))
    w(ixO^S)=w(ixO^S)/mypm%sqrtdetr%Gama(ixO^S)
  end where
 end if

end if
end subroutine  usemetric_devideby_sqrtspacedetr

!=============================================================================
subroutine  usemetric_getvec2_fromcov(ixI^L,ixO^L,nwa,wa,patchw,nwv,iwv,v,&
                                      iwe,ip0,p0,tkout_G,vcovphys)
 include 'amrvacdef.f'
 integer, intent(in)                  :: ixI^L,ixO^L
 integer, intent(in)                  :: nwa
 double precision, intent(in)         :: wa(ixI^S,1:nwa)
 logical, intent(in)                  :: patchw(ixG^T)
 integer, intent(in)                  :: nwv,iwv(nwv)
 type(vecnorm), intent(inout)         :: v(:)
 integer, optional                    :: iwe(2)
 integer, optional                    :: ip0(:)
 type(vector_phys), optional          :: p0(:)
 logical, optional                    :: tkout_G
 type(vector_phys), pointer, optional :: vcovphys(:)
! .. local variables ..
 type(vector_phys),target, allocatable   :: wcov(:)
 logical                         :: nozero
 integer                         :: idir,jdir,iidir,jjdir,i,nwv_sub
!-----------------------------------------------------------------------------
 nwv_sub=nwv
 if(present(p0)) then
   nwv_sub=nwv_sub+size(p0)
 end if
 allocate(wcov(nwv_sub))
 cond_ptkG : if(present(tkout_G)) then
  cond_tkG : if(tkout_G)then
   Loop_nwv0awG : do i=1,nwv
    where(patchw(ixO^S))
     ^C&wcov(i)%w(ixO^S,^C)=wa(ixO^S,iwv(i)+^C)/mypm%sqrtdetr%Gama(ixO^S);
    end where
   end do Loop_nwv0awG
   cond_p0wG : if(present(p0))then
    Loop_nwv0bwG: do i=nwv+1,nwv_sub
     where(patchw(ixO^S))
      {^C&wcov(i)%w(ixO^S,^C)=(wa(ixO^S,iwv(ip0(i-nwv))+^C)&
                            +p0(ip0(i-nwv))%w(ixO^S,^C))&
                             / mypm%sqrtdetr%Gama(ixO^S);}
     end where
    end do Loop_nwv0bwG
   end if cond_p0wG
   cond_vphys : if(present(vcovphys)) then
    vcovphys=>wcov
   end if cond_vphys
  else cond_tkG
   Loop_nwv0a : do i=1,nwv
    where(patchw(ixO^S))
     ^C&wcov(i)%w(ixO^S,^C)=wa(ixO^S,iwv(i)+^C);
    end where
   end do Loop_nwv0a
   if(present(p0))then
    Loop_nwv0b: do i=nwv+1,nwv_sub
     where(patchw(ixO^S))
      {^C&wcov(i)%w(ixO^S,^C)=wa(ixO^S,iwv(ip0(i-nwv))+^C)&
                            +p0(ip0(i-nwv))%w(ixO^S,^C);}
     end where
     end do Loop_nwv0b
    end if
  end if cond_tkG
 else cond_ptkG
   Loop_nwv0ae : do i=1,nwv
    where(patchw(ixO^S))
     ^C&wcov(i)%w(ixO^S,^C)=wa(ixO^S,iwv(i)+^C);
    end where
   end do Loop_nwv0ae
   if(present(p0))then
    Loop_nwv0be: do i=nwv+1,nwv_sub
     where(patchw(ixO^S))
      {^C&wcov(i)%w(ixO^S,^C)=wa(ixO^S,iwv(ip0(i-nwv))+^C)&
                            +p0(ip0(i-nwv))%w(ixO^S,^C);}
     end where
     end do Loop_nwv0be
    end if
 end if cond_ptkG
 nozero=.false.
 Loop_idir : do iidir=1,setgr%ndir
   if(all((/(.not.mypm%elemsp_inv(i1_gr(iidir,jjdir),&
              i2_gr(iidir,jjdir))%elpt%elm_on&
             ,jjdir=1,setgr%ndir)/))) cycle
   Loop_jdir : do jjdir=1,setgr%ndir;idir=i1_gr(iidir,jjdir);&
                                     jdir=i2_gr(iidir,jjdir)
    if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on)then
     cond_isone : if(.not.mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
      cond_nozeroA : if(nozero)then
       Loop_nwv1a : do i=1,nwv_sub
        where(patchw(ixO^S))
         v(i)%norm2(ixO^S)=v(i)%norm2(ixO^S)+ &
                   wcov(i)%w(ixO^S,jjdir)*wcov(i)%w(ixO^S,iidir)*&
                    mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S)
        end where
       end do Loop_nwv1a
       cond_iwe1 : if(present(iwe))then
        where(patchw(ixO^S))
         v(nwv_sub+1)%norm2(ixO^S)=v(nwv_sub+1)%norm2(ixO^S)+ &
           wcov(iwe(1))%w(ixO^S,jjdir)*wcov(iwe(1))%w(ixO^S,iidir)*&
           mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S)
        end where
       end if cond_iwe1
      else cond_nozeroA
       Loop_nwv1b : do i=1,nwv_sub
        where(patchw(ixO^S))
         v(i)%norm2(ixO^S)=wcov(i)%w(ixO^S,jjdir)**2.0d0&
                   *mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S)
        end where
       end do Loop_nwv1b
       cond_iwe2 : if(present(iwe))then
        where(patchw(ixO^S))
         v(nwv_sub+1)%norm2(ixO^S)=wcov(iwe(1))%w(ixO^S,iidir)&
                             *wcov(iwe(2))%w(ixO^S,jjdir)&
                    *mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S)
        end where
       end if cond_iwe2
      end if cond_nozeroA
     else cond_isone
      cond_nozeroB : if(nozero)then
       Loop_nwv1c : do i=1,nwv_sub
        where(patchw(ixO^S))
         v(i)%norm2(ixO^S)=v(i)%norm2(ixO^S)+wcov(i)%w(ixO^S,jjdir)*&
                                     wcov(i)%w(ixO^S,iidir)
        end where
       end do Loop_nwv1c
       cond_iweB1 : if(present(iwe))then
        where(patchw(ixO^S))
         v(nwv_sub+1)%norm2(ixO^S)=v(nwv_sub+1)%norm2(ixO^S)+&
                                   wcov(iwe(1))%w(ixO^S,jjdir)&
				   *wcov(iwe(2))%w(ixO^S,iidir)
        end where
       end if cond_iweB1
      else cond_nozeroB
       Loop_nwv1d : do i=1,nwv_sub
        where(patchw(ixO^S))
         v(i)%norm2(ixO^S)=wcov(i)%w(ixO^S,jjdir)*wcov(i)%w(ixO^S,iidir)
        end where
       end do Loop_nwv1d
       cond_iweB2 : if(present(iwe))then
        where(patchw(ixO^S))
         v(i+1)%norm2(ixO^S)=wcov(iwe(1))%w(ixO^S,iidir)*&
                             wcov(iwe(2))%w(ixO^S,jjdir)
        end where
       end if cond_iweB2
      end if cond_nozeroB
     end if cond_isone
     nozero=.true.
    end if
   end do Loop_jdir
 end do Loop_idir
 deallocate(wcov)
end subroutine  usemetric_getvec2_fromcov
!=============================================================================
subroutine  usemetric_getvdotb_fromcont(ixI^L,ixO^L,nwa,w,patchw,v,iwe,&
                                       tkout_G)
 include 'amrvacdef.f'
 integer, intent(in)                  :: ixI^L,ixO^L
 integer, intent(in)                  :: nwa
 double precision, intent(in)         :: w(ixI^S,1:nwa)
 logical, intent(in)                  :: patchw(ixG^T)
 integer, intent(in)                  :: iwe(2)
 type(vecnorm), intent(inout)         :: v
 logical, optional                    :: tkout_G
! .. local variables ..
 type(vector_phys),target, allocatable :: wcont(:)
 logical                               :: nozero
 integer                               :: idir,jdir,&
                                          iidir,jjdir,i,nwv_sub
!-----------------------------------------------------------------------------
 nwv_sub=2

 allocate(wcont(nwv_sub))
 cond_ptkG : if(present(tkout_G)) then
  cond_tkG : if(tkout_G)then
  Loop_nwv0aGT  : do i=1,nwv_sub
    where(patchw(ixO^S))
     ^C&wcont(i)%w(ixO^S,^C)=w(ixO^S,iwe(i)+^C)/mypm%sqrtdetr%Gama(ixO^S);
    end where
  end do Loop_nwv0aGT
  else cond_tkG
  Loop_nwv0aGf  : do i=1,nwv_sub
    where(patchw(ixO^S))
     ^C&wcont(i)%w(ixO^S,^C)=w(ixO^S,iwe(i)+^C);
    end where
  end do Loop_nwv0aGf
  end if cond_tkG
 else cond_ptkG
   Loop_nwv0a : do i=1,nwv_sub
    where(patchw(ixO^S))
     ^C&wcont(i)%w(ixO^S,^C)=w(ixO^S,iwe(i)+^C);
    end where
   end do Loop_nwv0a
end if cond_ptkG
 nozero=.false.
Loop_idir : do iidir=1,setgr%ndir
   if(all((/(.not.mypm%elemsp(i1_gr(iidir,jjdir),&
             i2_gr(iidir,jjdir))%elpt%elm_on,jjdir=1,setgr%ndir)/) ) ) cycle
   Loop_jdir : do jjdir=1,setgr%ndir;idir=i1_gr(iidir,jjdir)&
                                    ;jdir=i2_gr(iidir,jjdir)
    if(mypm%elemsp(idir,jdir)%elpt%elm_on)then
     cond_noone1 : if(.not.mypm%elemsp(idir,jdir)%elpt%elm_isone)then
      cond_nozero1 : if(nozero)then
        where(patchw(ixO^S))
         v%norm2(ixO^S)=v%norm2(ixO^S)+wcont(1)%w(ixO^S,jjdir)&
                           *wcont(2)%w(ixO^S,iidir)&
                           *mypm%elemsp(idir,jdir)%elpt%wg(ixO^S)
        end where
      else cond_nozero1
        where(patchw(ixO^S))
          v%norm2(ixO^S)=wcont(1)%w(ixO^S,jjdir)*wcont(2)%w(ixO^S,iidir)*&
                  mypm%elemsp(idir,jdir)%elpt%wg(ixO^S)
        end where
      end if cond_nozero1
     else cond_noone1
      cond_nozero2 : if(nozero)then
        where(patchw(ixO^S))
         v%norm2(ixO^S)=v%norm2(ixO^S)+wcont(1)%w(ixO^S,jjdir)*&
                        wcont(2)%w(ixO^S,iidir)
        end where
      else cond_nozero2
        where(patchw(ixO^S))
         v%norm2(ixO^S)=wcont(1)%w(ixO^S,jjdir)*wcont(2)%w(ixO^S,iidir)
        end where
      end if cond_nozero2
     end if cond_noone1
     nozero=.true.
    end if
   end do Loop_jdir
 end do Loop_idir


end subroutine  usemetric_getvdotb_fromcont
!=============================================================================
subroutine  usemetric_getvec2_fromcont(ixI^L,ixO^L,nwa,w,patchw,nwv,iwv,v,iwe,&
                                       ip0,p0,tkout_G,vcontphys)
 include 'amrvacdef.f'
 integer, intent(in)                  :: ixI^L,ixO^L
 integer, intent(in)                  :: nwa
 double precision, intent(in)         :: w(ixI^S,1:nwa)
 logical, intent(in)                  :: patchw(ixG^T)
 integer, intent(in)                  :: nwv,iwv(:)
 type(vecnorm), intent(inout)         :: v(:)
 integer, optional                    :: iwe(:)
 integer, optional                    :: ip0(:)
 type(vector_phys), optional          :: p0(:)
 logical, optional                    :: tkout_G
 type(vector_phys),  optional         :: vcontphys(nwv)
! .. local variables ..
 type(vector_phys),target, allocatable :: wcont(:)
 logical                               :: nozero
 integer                               :: idir,jdir,&
                                          iidir,jjdir,i,nwv_sub
!-----------------------------------------------------------------------------
 nwv_sub=nwv
 if(present(p0)) then
   nwv_sub=nwv_sub+size(p0)
 end if

 allocate(wcont(nwv_sub))
 cond_ptkG : if(present(tkout_G)) then
  cond_tkG : if(tkout_G)then
  Loop_nwv0aGT  : do i=1,nwv
    where(patchw(ixO^S))
     ^C&wcont(i)%w(ixO^S,^C)=w(ixO^S,iwv(i)+^C)/mypm%sqrtdetr%Gama(ixO^S);
    end where
  end do Loop_nwv0aGT
  cond_p0GT : if(present(p0)) then
   Loop_nwv0bGT :do i=nwv+1,nwv_sub
    where(patchw(ixO^S))
     {^C&wcont(i)%w(ixO^S,^C)=(w(ixO^S,iwv(ip0(i-nwv))+^C)+p0(i)%w(ixO^S,^C))&
                              / mypm%sqrtdetr%Gama(ixO^S);}
    end where
   end do Loop_nwv0bGT
  end if  cond_p0GT
  else cond_tkG
  Loop_nwv0aGf  : do i=1,nwv
    where(patchw(ixO^S))
     ^C&wcont(i)%w(ixO^S,^C)=w(ixO^S,iwv(i)+^C);
    end where
  end do Loop_nwv0aGf
  cond_p0Gf : if(present(p0)) then
   Loop_nwv0bGf :do i=nwv+1,nwv_sub
    where(patchw(ixO^S))
     {^C&wcont(i)%w(ixO^S,^C)=w(ixO^S,iwv(ip0(i-nwv))+^C)+p0(i)%w(ixO^S,^C);}
    end where
   end do Loop_nwv0bGf
  end if  cond_p0Gf
  end if cond_tkG
 else cond_ptkG
   Loop_nwv0a : do i=1,nwv
    where(patchw(ixO^S))
     ^C&wcont(i)%w(ixO^S,^C)=w(ixO^S,iwv(i)+^C);
    end where
   end do Loop_nwv0a
  cond_p0 : if(present(p0)) then
   Loop_nwv0b :do i=nwv+1,nwv_sub
    where(patchw(ixO^S))
     {^C&wcont(i)%w(ixO^S,^C)=w(ixO^S,iwv(ip0(i-nwv))+^C)+p0(i)%w(ixO^S,^C);}
    end where
   end do Loop_nwv0b
  end if  cond_p0
end if cond_ptkG
cond_vphys : if (present(vcontphys))then
  Loop_nwvphys  : do i=1,nwv

   where(patchw(ixO^S))
    ^C&vcontphys(i)%w(ixO^S,^C)=wcont(i)%w(ixO^S,^C);
   end where
 end do Loop_nwvphys
end if cond_vphys

nozero=.false.
Loop_idir : do iidir=1,setgr%ndir
   if(all((/(.not.mypm%elemsp(i1_gr(iidir,jjdir),&
             i2_gr(iidir,jjdir))%elpt%elm_on,jjdir=1,setgr%ndir)/) ) ) cycle
   Loop_jdir : do jjdir=1,setgr%ndir;idir=i1_gr(iidir,jjdir)&
                                    ;jdir=i2_gr(iidir,jjdir)
    if(mypm%elemsp(idir,jdir)%elpt%elm_on)then
     if(.not.mypm%elemsp(idir,jdir)%elpt%elm_isone)then
      cond_nozero1 : if(nozero)then
       Loop_nwv1a : do i=1,nwv_sub
        where(patchw(ixO^S))
         v(i)%norm2(ixO^S)=v(i)%norm2(ixO^S)+wcont(i)%w(ixO^S,jjdir)*&
                           wcont(i)%w(ixO^S,iidir)*&
                           mypm%elemsp(idir,jdir)%elpt%wg(ixO^S)
        end where
       end do Loop_nwv1a
       if(present(iwe))then
         v(nwv_sub+1)%norm2(ixO^S)=v(nwv_sub+1)%norm2(ixO^S)+&
                    wcont(iwe(1))%w(ixO^S,jjdir)*&
                    wcont(iwe(2))%w(ixO^S,iidir)*&
                    mypm%elemsp(idir,jdir)%elpt%wg(ixO^S)
       end if
      else cond_nozero1
       Loop_nwv1b : do i=1,nwv_sub
        where(patchw(ixO^S))
          v(i)%norm2(ixO^S)=wcont(i)%w(ixO^S,jjdir)*wcont(i)%w(ixO^S,iidir)*&
                  mypm%elemsp(idir,jdir)%elpt%wg(ixO^S)
        end where
       end do Loop_nwv1b
       if(present(iwe))then
         v(nwv_sub+1)%norm2(ixO^S)=wcont(iwe(1))%w(ixO^S,jjdir)*&
                    wcont(iwe(2))%w(ixO^S,iidir)*&
                    mypm%elemsp(idir,jdir)%elpt%wg(ixO^S)
       end if

      end if cond_nozero1
     else
      if(nozero)then
       Loop_nwv1c : do i=1,nwv_sub
        where(patchw(ixO^S))
         v(i)%norm2(ixO^S)=v(i)%norm2(ixO^S)+wcont(i)%w(ixO^S,jjdir)*&
                                            wcont(i)%w(ixO^S,iidir)
        end where
       end do Loop_nwv1c
       if(present(iwe))then
         v(nwv_sub+1)%norm2(ixO^S)=v(nwv_sub+1)%norm2(ixO^S)+&
                    wcont(iwe(1))%w(ixO^S,jjdir)*&
                    wcont(iwe(2))%w(ixO^S,iidir)
       end if

      else
       Loop_nwv1d : do i=1,nwv_sub
        where(patchw(ixO^S))
         v(i)%norm2(ixO^S)=wcont(i)%w(ixO^S,jjdir)*wcont(i)%w(ixO^S,iidir)
        end where
       end do Loop_nwv1d
       if(present(iwe))then
         v(nwv_sub+1)%norm2(ixO^S)=wcont(iwe(1))%w(ixO^S,jjdir)*&
                    wcont(iwe(2))%w(ixO^S,iidir)
       end if

      end if
     end if
     nozero=.true.
    end if
   end do Loop_jdir
 end do Loop_idir


end subroutine  usemetric_getvec2_fromcont
!=============================================================================
subroutine usemetric_multiby_sqrtelemsp(ixO^L,iidir,jjdir,patchw,&
                                        scalw_in,scalw_out,&
                                        nozero)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixO^L,iidir,jjdir
 logical, intent(in)             :: patchw(ixG^T)
 double precision, intent(in)    :: scalw_in(ixO^S)
 double precision, intent(out)   :: scalw_out(ixG^T)
 logical, intent(out)            :: nozero
 ! .. local ..
 integer                         :: idir,jdir
!-----------------------------------------------------------------------------
 idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
 if(mypm%elemsp(idir,jdir)%elpt%elm_on)then
  nozero=.true.
  if(mypm%elemsp(idir,jdir)%elpt%elm_isone)then
   where(patchw(ixO^S))
    scalw_out(ixO^S)=scalw_in(ixO^S)
   end where
  else 
   where(patchw(ixO^S))
    scalw_out(ixO^S)=scalw_in(ixO^S)&   
                     *dsqrt(mypm%elemsp(idir,jdir)%elpt%wg(ixO^S))
   end where
  end if
 else
  nozero=.false.
 end if

end subroutine usemetric_multiby_sqrtelemsp

!=============================================================================
subroutine usemetric_multiby_sqrtelemsp_inv(ixO^L,iidir,jjdir,patchw,&
                                            scalw_in,scalw_out,nozero)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixO^L,iidir,jjdir
 logical, intent(in)             :: patchw(ixG^T)
 double precision, intent(in)    :: scalw_in(ixO^S)
 double precision, intent(out)   :: scalw_out(ixG^T)
 logical, intent(out)            :: nozero
 ! .. local ..
 integer                         :: idir,jdir
!-----------------------------------------------------------------------------
 idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
 if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on)then
  nozero=.true.
  if(mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
   where(patchw(ixO^S))
    scalw_out(ixO^S)=scalw_in(ixO^S)
   end where
  else
   where(patchw(ixO^S))
    scalw_out(ixO^S)=scalw_in(ixO^S)*&
                     dsqrt(mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S))
   end where
  end if
 else
  nozero=.false.
 end if

end subroutine usemetric_multiby_sqrtelemsp_inv
!=============================================================================
subroutine usemetric_multiby_elemsp(ixO^L,iidir,jjdir,patchw,&
                                    scalw_in,scalw_out,nozero)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixO^L,iidir,jjdir
 logical, intent(in)             :: patchw(ixG^T)
 double precision, intent(in)    :: scalw_in(ixO^S)
 double precision, intent(out)   :: scalw_out(ixG^T)
 logical, intent(out)            :: nozero
 ! .. local ..
 integer                         :: idir,jdir
!-----------------------------------------------------------------------------
 idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
 if(mypm%elemsp(idir,jdir)%elpt%elm_on)then
  nozero=.true.
  if(mypm%elemsp(idir,jdir)%elpt%elm_isone)then
   where(patchw(ixO^S))
    scalw_out(ixO^S)=scalw_in(ixO^S)
   end where
  else
   where(patchw(ixO^S))
    scalw_out(ixO^S)=scalw_in(ixO^S)*mypm%elemsp(idir,jdir)%elpt%wg(ixO^S)
   end where
  end if
 else         
  nozero=.false.
 end if
    
end subroutine usemetric_multiby_elemsp
    
!=============================================================================
subroutine usemetric_multiby_elemsp_inv(ixO^L,iidir,jjdir,patchw,scalw_in,&
                                        scalw_out,nozero)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixO^L,iidir,jjdir
 logical, intent(in)             :: patchw(ixG^T)
 double precision, intent(in)    :: scalw_in(ixO^S)
 double precision, intent(out)   :: scalw_out(ixG^T)
 logical, intent(out)            :: nozero
 ! .. local ..
 integer                         :: idir,jdir
!-----------------------------------------------------------------------------
 idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
 if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on)then
  nozero=.true.
  if(mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
   where(patchw(ixO^S))
    scalw_out(ixO^S)=scalw_in(ixO^S)
   end where
  else
   where(patchw(ixO^S))
    scalw_out(ixO^S)=scalw_in(ixO^S)*mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S)
   end where
  end if
 else
  nozero=.false.
 end if

end subroutine usemetric_multiby_elemsp_inv
!=============================================================================
subroutine usemetric_add_elemsp_inv(ixO^L,iidir,jjdir,patchw,scalw_in,&
                                        scalw_out)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixO^L,iidir,jjdir
 logical, intent(in)             :: patchw(ixG^T)
 double precision, intent(in)    :: scalw_in(ixO^S)
 double precision, intent(out)   :: scalw_out(ixG^T)
 ! .. local ..
 integer                         :: idir,jdir
!-----------------------------------------------------------------------------
 idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
 if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on)then
  if(mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
   where(patchw(ixO^S))
    scalw_out(ixO^S)=scalw_in(ixO^S)+1
   end where
  else
   where(patchw(ixO^S))
    scalw_out(ixO^S)=scalw_in(ixO^S)+mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S)
   end where
  end if
 else
  scalw_out(ixO^S)=scalw_in(ixO^S)
 end if

end subroutine usemetric_add_elemsp_inv
!=============================================================================
subroutine usemetric_multiby_elemspt_inv(ixI^L,ixO^L,iidir,jjdir,patchw,scalw_in,&
                                        scalw_out,nozero)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixI^L,ixO^L,iidir,jjdir
 logical, intent(in)             :: patchw(ixG^T)
 double precision, intent(in)    :: scalw_in(ixI^S)
 double precision, intent(out)   :: scalw_out(ixG^T)
 logical, intent(out)            :: nozero
 ! .. local ..
 integer                         :: idir,jdir
!-----------------------------------------------------------------------------
 idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
 cond_1 : if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on& 
         .or.(mypm%bt_cont(iidir)%sftpt%elm_on.and.&
             mypm%bt_cont(jjdir)%sftpt%elm_on))then
  nozero=.true.
  condbi1 : if(mypm%bt_cont(iidir)%sftpt%elm_on.and.&
    mypm%bt_cont(jjdir)%sftpt%elm_on)then
    condbi1_isone :if(mypm%bt_cont(iidir)%sftpt%elm_isone)then
      if(mypm%bt_cont(jjdir)%sftpt%elm_isone)then
        if(mypm%alfa%elm_isone)then
           if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on) then
            if(mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
              nozero=.false.
            else
             scalw_out(ixO^S)=scalw_in(ixO^S)*(mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S)&
                               -1.0D0)
            end if
          else
           scalw_out(ixO^S)=-scalw_in(ixO^S)
          end if
        else
         if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on) then
          if(mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
            scalw_out(ixO^S)=scalw_in(ixO^S)*(1.0d0-1.0D0/mypm%alfa%wg(ixO^S))
          else
           scalw_out(ixO^S)=scalw_in(ixO^S)*(mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S)&
                               -1.0D0/mypm%alfa%wg(ixO^S))
          end if
         else
          scalw_out(ixO^S)=-scalw_in(ixO^S)/mypm%alfa%wg(ixO^S)
         end if
        end if
      else
        if(mypm%alfa%elm_isone)then
         if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on) then
          if(mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
            scalw_out(ixO^S)=scalw_in(ixO^S)*(1.0d0-mypm%bt_cont(jjdir)%sftpt%wg(ixO^S))
          else
           scalw_out(ixO^S)=scalw_in(ixO^S)*(mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S)&
                               -mypm%bt_cont(jjdir)%sftpt%wg(ixO^S))
          end if
         else
          scalw_out(ixO^S)=-scalw_in(ixO^S)*mypm%bt_cont(jjdir)%sftpt%wg(ixO^S)
         end if
        else
         if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on) then
          if(mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
            scalw_out(ixO^S)=scalw_in(ixO^S)*(1.0d0-mypm%bt_cont(jjdir)%sftpt%wg(ixO^S)&
                                         /mypm%alfa%wg(ixO^S))
          else
           scalw_out(ixO^S)=scalw_in(ixO^S)*(mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S)&
                               -mypm%bt_cont(jjdir)%sftpt%wg(ixO^S)/mypm%alfa%wg(ixO^S))
          end if
         else
          scalw_out(ixO^S)=-scalw_in(ixO^S)*mypm%bt_cont(jjdir)%sftpt%wg(ixO^S)&
                            /mypm%alfa%wg(ixO^S)
         end if
        end if
      end if
    else condbi1_isone
      if(mypm%bt_cont(jjdir)%sftpt%elm_isone)then
        if(mypm%alfa%elm_isone)then
         if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on) then
          if(mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
            scalw_out(ixO^S)=scalw_in(ixO^S)*(1.0d0-mypm%bt_cont(iidir)%sftpt%wg(ixO^S))
          else
           scalw_out(ixO^S)=scalw_in(ixO^S)*(mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S)&
                               -mypm%bt_cont(iidir)%sftpt%wg(ixO^S))
          end if
         else
          scalw_out(ixO^S)=-scalw_in(ixO^S)*mypm%bt_cont(iidir)%sftpt%wg(ixO^S)
         end if
        else
         if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on) then
          if(mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
            scalw_out(ixO^S)=scalw_in(ixO^S)*(1.0d0-mypm%bt_cont(iidir)%sftpt%wg(ixO^S)&
                                         /mypm%alfa%wg(ixO^S))
          else
           scalw_out(ixO^S)=scalw_in(ixO^S)*(mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S)&
                               -mypm%bt_cont(iidir)%sftpt%wg(ixO^S)/mypm%alfa%wg(ixO^S))
          end if
         else
          scalw_out(ixO^S)=-scalw_in(ixO^S)*mypm%bt_cont(iidir)%sftpt%wg(ixO^S)&
                            /mypm%alfa%wg(ixO^S)
         end if
        end if
      else
        if(mypm%alfa%elm_isone)then
         if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on) then
          if(mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
            scalw_out(ixO^S)=scalw_in(ixO^S)*(1.0d0-mypm%bt_cont(iidir)%sftpt%wg(ixO^S)&
                               *mypm%bt_cont(jjdir)%sftpt%wg(ixO^S))
          else
           scalw_out(ixO^S)=scalw_in(ixO^S)*(mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S)&
                               -mypm%bt_cont(iidir)%sftpt%wg(ixO^S)&
                                *mypm%bt_cont(jjdir)%sftpt%wg(ixO^S))
          end if
         else
          scalw_out(ixO^S)=-scalw_in(ixO^S)*mypm%bt_cont(iidir)%sftpt%wg(ixO^S)&
                            *mypm%bt_cont(jjdir)%sftpt%wg(ixO^S)
         end if
        else
         if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on) then
          if(mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
            scalw_out(ixO^S)=scalw_in(ixO^S)*(1.0d0-mypm%bt_cont(iidir)%sftpt%wg(ixO^S)&
                 *mypm%bt_cont(jjdir)%sftpt%wg(ixO^S) /mypm%alfa%wg(ixO^S))
          else
           scalw_out(ixO^S)=scalw_in(ixO^S)*(mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S)&
                      -mypm%bt_cont(iidir)%sftpt%wg(ixO^S)&
                      *mypm%bt_cont(jjdir)%sftpt%wg(ixO^S)/mypm%alfa%wg(ixO^S))
          end if
         else
          scalw_out(ixO^S)=-scalw_in(ixO^S)*mypm%bt_cont(iidir)%sftpt%wg(ixO^S)&
                            *mypm%bt_cont(jjdir)%sftpt%wg(ixO^S)/mypm%alfa%wg(ixO^S)
         end if
        end if
      end if
    end if condbi1_isone
  else condbi1
   scalw_out(ixO^S)=scalw_in(ixO^S)*mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S)
  end if condbi1
 else cond_1
  nozero=.false.
 end if cond_1
end subroutine usemetric_multiby_elemspt_inv
!=============================================================================
subroutine usemetric_multiby_elemsp_updown(ixO^L,iidir,jjdir,patchw,scalw_in,&
                                           nozero,scalw_out)
 include 'amrvacdef.f'
 integer, intent(in)                       :: ixO^L,iidir,jjdir
 logical, intent(in)                       :: patchw(ixG^T)
 double precision, intent(inout)           :: scalw_in(ixO^S)
 logical, intent(out)                      :: nozero
 double precision, optional, intent(out)   :: scalw_out(ixG^T)
 ! .. local ..
 integer                                   :: idir,jdir
!-----------------------------------------------------------------------------
 idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
 if(mypm%elemsp_updwn(idir,jdir)%elpt%elm_on)then
  nozero=.true.
  if(mypm%elemsp_updwn(idir,jdir)%elpt%elm_isone)then
   if(present(scalw_out))then
    where(patchw(ixO^S))
     scalw_out(ixO^S)=scalw_in(ixO^S)
    end where
   end if
  else
   if(present(scalw_out))then
    where(patchw(ixO^S))
     scalw_out(ixO^S)=scalw_in(ixO^S)*&
                       mypm%elemsp_updwn(idir,jdir)%elpt%wg(ixO^S)
    end where
   else
    where(patchw(ixO^S))
     scalw_in(ixO^S)=scalw_in(ixO^S)*&
         mypm%elemsp_updwn(idir,jdir)%elpt%wg(ixO^S)  
    end where
   end if
  end if
 else
  nozero=.false.
 end if

end subroutine usemetric_multiby_elemsp_updown
!=============================================================================
subroutine usemetric_multiby_sqrtelemsp_updown(ixO^L,iidir,jjdir,patchw,&
                                               scalw_in,nozero,scalw_out)
 include 'amrvacdef.f'
 integer, intent(in)                     :: ixO^L,iidir,jjdir
 logical, intent(in)                     :: patchw(ixG^T)
 double precision, intent(inout)         :: scalw_in(ixO^S)
 double precision, optional, intent(out) :: scalw_out(ixG^T)
 logical, intent(out)                    :: nozero
 ! .. local ..
 integer                                 :: idir,jdir
!-----------------------------------------------------------------------------
 idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
 if(mypm%elemsp_updwn(idir,jdir)%elpt%elm_on)then
  nozero=.true.
  if(mypm%elemsp_updwn(idir,jdir)%elpt%elm_isone)then
   if(present(scalw_out))then
    where(patchw(ixO^S))
     scalw_out(ixO^S)=scalw_in(ixO^S)
    end where
   end if
  else
   if(present(scalw_out))then
    where(patchw(ixO^S))
     scalw_out(ixO^S)=scalw_in(ixO^S)*&
                      dsqrt(mypm%elemsp_updwn(idir,jdir)%elpt%wg(ixO^S))
    end where
   else
    where(patchw(ixO^S))
     scalw_in(ixO^S)=scalw_in(ixO^S)*&
                      dsqrt(mypm%elemsp_updwn(idir,jdir)%elpt%wg(ixO^S))
    end where
   end if
  end if
 else
  nozero=.false.
 end if

end subroutine usemetric_multiby_sqrtelemsp_updown
!=============================================================================
subroutine usemetric_multiby_deriveelemsp(ixO^L,iidir,jjdir,drv_idir,patchw,&
                                          scalw_in,scalw_out,nozero)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixO^L,iidir,jjdir,drv_idir
 logical, intent(in)             :: patchw(ixG^T)
 double precision, intent(in)    :: scalw_in(ixO^S)
 double precision, intent(out)   :: scalw_out(ixG^T)
 logical, intent(out)            :: nozero
 ! .. local ..
 integer                         :: idir,jdir
!-----------------------------------------------------------------------------
 idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
 if(mypm%elemsp(idir,jdir)%elpt%drv(drv_idir)%elm_on)then
  nozero=.true.
  if(mypm%elemsp(idir,jdir)%elpt%drv(drv_idir)%elm_isone)then
   where(patchw(ixO^S))
    scalw_out(ixO^S)=scalw_in(ixO^S)
   end where
  else
   where(patchw(ixO^S))
    scalw_out(ixO^S)=scalw_in(ixO^S)&
                     *mypm%elemsp(idir,jdir)%elpt%drv(drv_idir)%dwg(ixO^S)
   end where
  end if
 else
  nozero=.false.
 end if

end subroutine usemetric_multiby_deriveelemsp
!=============================================================================
subroutine usemetric_multiby_derivealpha(ixO^L,drv_idir,scalw_in,scalw_out,nozero)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixO^L,drv_idir
 double precision, intent(in)    :: scalw_in(ixO^S)
 double precision, intent(out)   :: scalw_out(ixG^T)
 logical, intent(out)            :: nozero
!-----------------------------------------------------------------------------
 if(mypm%alfa%drv(drv_idir)%elm_on)then
  nozero=.true.
  if(mypm%alfa%drv(drv_idir)%elm_isone)then
   scalw_out(ixO^S)=scalw_in(ixO^S)
  else
   scalw_out(ixO^S)=scalw_in(ixO^S)*mypm%alfa%drv(drv_idir)%dwg(ixO^S)
  end if
 else
  nozero=.false.
 end if

end subroutine usemetric_multiby_derivealpha
!=============================================================================
!=============================================================================
subroutine usemetric_multiby_alpha(ixO^L,scalw_in,scalw_out)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixO^L
 double precision, intent(inout) :: scalw_in(ixO^S)
 double precision, optional      :: scalw_out(ixG^T)
!-----------------------------------------------------------------------------
  if(.not.mypm%alfa%elm_isone)then
   if(present(scalw_out))then
    scalw_out(ixO^S)=scalw_in(ixO^S)*mypm%alfa%wg(ixO^S)
   else
    scalw_in(ixO^S)=scalw_in(ixO^S)*mypm%alfa%wg(ixO^S)
   end if
  else   
  if(present(scalw_out))then
    scalw_out(ixO^S)=scalw_in(ixO^S)
   end if
  end if


end subroutine usemetric_multiby_alpha
!=============================================================================
subroutine usemetric_devideby_alpha(ixO^L,scalw_in,scalw_out)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixO^L
 double precision, intent(inout) :: scalw_in(ixO^S)
 double precision, optional      :: scalw_out(ixG^T)
!-----------------------------------------------------------------------------
  if(.not.mypm%alfa%elm_isone)then
   if(present(scalw_out))then
    scalw_out(ixO^S)=scalw_in(ixO^S)/mypm%alfa%wg(ixO^S)
   else
    scalw_in(ixO^S)=scalw_in(ixO^S)/mypm%alfa%wg(ixO^S)
   end if
  else
  if(present(scalw_out))then
    scalw_out(ixO^S)=scalw_in(ixO^S)
   end if
  end if


end subroutine usemetric_devideby_alpha
!=============================================================================

subroutine usemetric_devideby_elemsp(ixO^L,iidir,jjdir,scalw_in,scalw_out,nozero)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixO^L,iidir,jjdir
 double precision, intent(in)    :: scalw_in(ixO^S)
 double precision, intent(out)   :: scalw_out(ixG^T)
 logical, intent(out)            :: nozero
 ! .. local ..
 integer                         :: idir,jdir
!-----------------------------------------------------------------------------
 idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)

 if(mypm%elemsp(idir,jdir)%elpt%elm_on)then
  nozero=.true.
  if(mypm%elemsp(idir,jdir)%elpt%elm_isone)then
   scalw_out(ixO^S)=scalw_in(ixO^S)
  else
   scalw_out(ixO^S)=scalw_in(ixO^S)/mypm%elemsp(idir,jdir)%elpt%wg(ixO^S)
  end if
 else
  nozero=.false.
 end if

end subroutine usemetric_devideby_elemsp
!=============================================================================
subroutine usemetric_devideby_sqrtelemsp(ixO^L,iidir,jjdir,&
                                         scalw_in,scalw_out,nozero)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixO^L,iidir,jjdir
 double precision, intent(in)    :: scalw_in(ixO^S)
 double precision, intent(out)   :: scalw_out(ixG^T)
 logical, intent(out)            :: nozero
 ! .. local ..
 integer                         :: idir,jdir
!-----------------------------------------------------------------------------
 idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
 if(mypm%elemsp(idir,jdir)%elpt%elm_on)then
  nozero=.true.
  if(mypm%elemsp(idir,jdir)%elpt%elm_isone)then
   scalw_out(ixO^S)=scalw_in(ixO^S)
  else
   scalw_out(ixO^S)=scalw_in(ixO^S)/dsqrt(mypm%elemsp(idir,jdir)%elpt%wg(ixO^S))
  end if
 else
  nozero=.false.
 end if

end subroutine usemetric_devideby_sqrtelemsp
!=============================================================================
!=============================================================================
subroutine usemetric_multiby_beta(ixO^L,idir,scalw_in,scalw_out,nozero)
 include 'amrvacdef.f'
 integer, intent(in)                       :: ixO^L,idir
 double precision, intent(inout)           :: scalw_in(ixO^S)
 double precision, optional, intent(out)   :: scalw_out(ixG^T)
 logical, intent(out)                      :: nozero
 ! .. local ..
 logical                                         ::beta_on(setgr%ndir)
 double precision, dimension(ixO^S,1:setgr%ndir) :: beta_cov
!-----------------------------------------------------------------------------
 call usemetric_betacov(ixO^L,ixO^L,idir,patchtrue,beta_cov,beta_on)
 if(beta_on(idir))then
  nozero=.true.
   if(present(scalw_out))then
    scalw_out(ixO^S)=scalw_in(ixO^S)*beta_cov(ixO^S,idir)
   else
    scalw_in(ixO^S)=scalw_in(ixO^S)*beta_cov(ixO^S,idir)
   end if
 else
  nozero=.false.
 end if

end subroutine usemetric_multiby_beta
!=============================================================================
subroutine usemetric_multiby_beta_cont(ixO^L,idir,scalw_in,scalw_out,nozero)
 include 'amrvacdef.f'
 integer, intent(in)                      :: ixO^L,idir
 double precision, intent(inout)          :: scalw_in(ixO^S)
 double precision,optional, intent(out)   :: scalw_out(ixG^T)
 logical, optional,intent(out)            :: nozero
!-----------------------------------------------------------------------------
 if(mypm%bt_cont(idir)%sftpt%elm_on)then
  if(present(nozero))nozero=.true.
  if(mypm%bt_cont(idir)%sftpt%elm_isone)then
   if(present(scalw_out))scalw_out(ixO^S)=scalw_in(ixO^S)
  else
   if(present(scalw_out))then
    scalw_out(ixO^S)=scalw_in(ixO^S)*mypm%bt_cont(idir)%sftpt%wg(ixO^S)
   else
    scalw_in(ixO^S)=scalw_in(ixO^S)*mypm%bt_cont(idir)%sftpt%wg(ixO^S)
   end if
  end if
 else
  if(present(nozero))nozero=.false.
 end if

end subroutine usemetric_multiby_beta_cont
!=============================================================================
subroutine usemetric_multiby_derivebeta_cont(ixO^L,idir,drv_idir,&
                                             scalw_in,scalw_out,nozero)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixO^L,idir,drv_idir
 double precision, intent(in)    :: scalw_in(ixO^S)
 double precision, intent(out)   :: scalw_out(ixG^T)
 logical, intent(out)            :: nozero
!-----------------------------------------------------------------------------
 if(mypm%bt_cont(idir)%sftpt%drv(drv_idir)%elm_on)then
  nozero=.true.
  if(mypm%bt_cont(idir)%sftpt%drv(drv_idir)%elm_isone)then
   scalw_out(ixO^S)=scalw_in(ixO^S)
  else
   scalw_out(ixO^S)=scalw_in(ixO^S)&
                    *mypm%bt_cont(idir)%sftpt%drv(drv_idir)%dwg(ixO^S)
  end if
 else
  nozero=.false.
 end if

end subroutine usemetric_multiby_derivebeta_cont

!=============================================================================
subroutine usemetric_covtocontr(ixI^L,ixO^L,nwv,iwv,iidir_out,patchw,nwl,w,wvcont&
                                ,tkout_G,wcovphys,allvec)
 ! u_i to u^i
 include 'amrvacdef.f' 
 integer, intent(in)                   :: ixI^L,ixO^L,nwl,iidir_out
 logical, intent(in)                   :: patchw(ixG^T)
 double precision, intent(in)          :: w(ixI^S,1:nwl)
 integer, intent(in)                   :: nwv,iwv(:)
 type(vector_phys), intent(inout)      :: wvcont(nwv)
 logical, optional                     :: tkout_G
 type(vector_phys) ,optional           :: wcovphys(nwv)
 logical, optional                     :: allvec
! .. local variables ..

 integer                               :: iw0,i,jdir,kdir,jjdir,kkdir,idir
 integer                               :: idir_min,idir_max,iidir
 type(vector_phys), target             :: wcov(nwv) 
 logical                               :: nozero
!----------------------------------------------------------------------------
if(present(allvec))then
 if(allvec)then
  idir_min=1;idir_max=ndir
 else
  idir_min=iidir_out;idir_max=iidir_out
 end if
else
 idir_min=iidir_out;idir_max=iidir_out
end if
cond_ptkG : if(present(tkout_G)) then 
 cond_tkG : if(tkout_G)then
  Loop_ivect0tkG : do i=1,nwv;iw0=iwv(i)
   where(patchw(ixO^S))
    ^C&wcov(i)%w(ixO^S,^C)=w(ixO^S,iw0+^C)/mypm%sqrtdetr%Gama(ixO^S);
   end where
  end do Loop_ivect0tkG
  cond_wp : if(present(wcovphys))then
 Loop_ivect0tkGphys : do i=1,nwv
   where(patchw(ixO^S))
   ^C&wcovphys(i)%w(ixO^S,^C)=^C&wcov(i)%w(ixO^S,^C);
   end where
 end do Loop_ivect0tkGphys
!i   wcovphys=>wcov
  end if cond_wp

 else cond_tkG
  Loop_ivect0kpG : do i=1,nwv;iw0=iwv(i)
   where(patchw(ixO^S))
    ^C&wcov(i)%w(ixO^S,^C)=w(ixO^S,iw0+^C);
   end where
  end do Loop_ivect0kpG
 endif cond_tkG
else cond_ptkG
  Loop_ivect0wG : do i=1,nwv;iw0=iwv(i)
   where(patchw(ixO^S))
    ^C&wcov(i)%w(ixO^S,^C)=w(ixO^S,iw0+^C);
   end where
  end do Loop_ivect0wG
end if cond_ptkG

Loop_allvec : do iidir = idir_min,idir_max
 nozero=.false.
 cond_alof :if(all((/(.not.mypm%elemsp_inv(i1_gr(iidir,jjdir),&
            i2_gr(iidir,jjdir))%elpt%elm_on,jjdir=1,setgr%ndir)/)))then
  Loop_ivect1 : do i=1,nwv;iw0=iwv(i)
   where(patchw(ixO^S))
     wvcont(i)%w(ixO^S,iidir)=zero
   end where
  end do Loop_ivect1
  return;
 end if cond_alof
 Loop_jdir : do jjdir=1,setgr%ndir;idir=i1_gr(iidir,jjdir)
                                   jdir=i2_gr(iidir,jjdir)
  cond_ijon : if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on)then
   cond_ijone : if(.not.mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
    cond_nozero2: if(nozero)then
     Loop_ivect2add : do i=1,nwv;iw0=iwv(i)
      where(patchw(ixO^S))
       wvcont(i)%w(ixO^S,iidir)=wvcont(i)%w(ixO^S,iidir)+&
                                wcov(i)%w(ixO^S,jjdir)*&
                                mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S);
      end where
     end do Loop_ivect2add
    else cond_nozero2
     Loop_ivect2 : do i=1,nwv;iw0=iwv(i)
      where(patchw(ixO^S))
       wvcont(i)%w(ixO^S,iidir)=wcov(i)%w(ixO^S,jjdir)*&
                             mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S);
      end where
     end do Loop_ivect2
    end if cond_nozero2
   else cond_ijone
    cond_nozero3: if(nozero)then
     Loop_ivect3add : do i=1,nwv
      where(patchw(ixO^S))
       wvcont(i)%w(ixO^S,iidir)=wvcont(i)%w(ixO^S,iidir)+&
                                wcov(i)%w(ixO^S,jjdir);
      end where
     end do Loop_ivect3add
    else cond_nozero3
     Loop_ivect3 : do i=1,nwv
      where(patchw(ixO^S))
       wvcont(i)%w(ixO^S,iidir)=wcov(i)%w(ixO^S,jjdir);
      end where
     end do Loop_ivect3
    end if cond_nozero3
   end if cond_ijone
   nozero=.true.
   end if cond_ijon
  end do Loop_jdir
end do Loop_allvec
end subroutine usemetric_covtocontr
!=============================================================================
subroutine usemetric_contrtocov(ixI^L,ixO^L,nwv,iwv,iidir_out,patchw,nwl,w,wvcov,&
                                tkout_G,wcontphys,allvec)
 ! u_^ to u_i
 include 'amrvacdef.f'
 integer, intent(in)                   :: ixI^L,ixO^L,iidir_out,nwl
 logical, intent(in)                   :: patchw(ixG^T)
 double precision, intent(in)          :: w(ixI^S,1:nwl)
 integer, intent(in)                   :: nwv,iwv(:)
 type(vector_phys), intent(inout)      :: wvcov(nwv)
 logical, optional                     :: tkout_G
 type(vector_phys),          optional  :: wcontphys(nwv)
 logical, optional                     :: allvec
! .. local variables ..
 integer                               :: iw0,i,jdir,idir,jjdir
 integer                               :: idir_min,idir_max,iidir
 type(vector_phys),target              :: wvcont(nwv)
 logical                               :: nozero
!----------------------------------------------------------------------------
if(present(allvec))then
 if(allvec)then
  idir_min=1;idir_max=ndir
 else
  idir_min=iidir_out;idir_max=iidir_out
 end if
else
 idir_min=iidir_out;idir_max=iidir_out
end if
cond_ptkG : if(present(tkout_G)) then
 cond_tkG : if(tkout_G)then
  Loop_ivect0tkG : do i=1,nwv;iw0=iwv(i)
   where(patchw(ixO^S))
    ^C&wvcont(i)%w(ixO^S,^C)=w(ixO^S,iw0+^C)/mypm%sqrtdetr%Gama(ixO^S);
   end where
  end do Loop_ivect0tkG
  cond_wp : if(present(wcontphys)) then
  Loop_ivect0tkGphys : do i=1,nwv;iw0=iwv(i)
   where(patchw(ixO^S))
       ^C&wcontphys(i)%w(ixO^S,^C)=wvcont(i)%w(ixO^S,^C);
   end where
  end do Loop_ivect0tkGphys

  end if cond_wp
 else cond_tkG
  Loop_ivect0wG : do i=1,nwv;iw0=iwv(i)
  where(patchw(ixO^S))
   ^C&wvcont(i)%w(ixO^S,^C)=w(ixO^S,iw0+^C);
  end where
 end do Loop_ivect0wG
  end if cond_tkG
else cond_ptkG
  Loop_ivect0wGn : do i=1,nwv;iw0=iwv(i)
  where(patchw(ixO^S))
   ^C&wvcont(i)%w(ixO^S,^C)=w(ixO^S,iw0+^C);
  end where
 end do Loop_ivect0wGn
end if cond_ptkG
Loop_allvec : do iidir = idir_min,idir_max
 nozero=.false.
cond_alof : if(all((/(.not.mypm%elemsp(i1_gr(iidir,jjdir),i2_gr(iidir,jjdir))&
                      %elpt%elm_on,jjdir=1,setgr%ndir)/)))then
 Loop_ivect1 : do i=1,nwv;iw0=iwv(i)
   where(patchw(ixO^S))
    wvcov(i)%w(ixO^S,iidir)=zero
   end where
 end do Loop_ivect1
 return;
end if cond_alof

Loop_jdir : do jjdir=1,setgr%ndir;idir=i1_gr(iidir,jjdir);
                                  jdir=i2_gr(iidir,jjdir)
 cond_ison :if(mypm%elemsp(idir,jdir)%elpt%elm_on)then
     cond_isone : if(.not.mypm%elemsp(idir,jdir)%elpt%elm_isone)then
      if(nozero)then
      Loop_ivect2add : do i=1,nwv;iw0=iwv(i)
        where(patchw(ixO^S))
         wvcov(i)%w(ixO^S,iidir)=wvcov(i)%w(ixO^S,iidir)+&
                                wvcont(i)%w(ixO^S,jjdir)*&
                               mypm%elemsp(idir,jdir)%elpt%wg(ixO^S);
        end where
       end do Loop_ivect2add
      else
       Loop_ivect2 : do i=1,nwv;iw0=iwv(i)
        where(patchw(ixO^S))
         wvcov(i)%w(ixO^S,iidir)=wvcont(i)%w(ixO^S,jjdir)*&
                               mypm%elemsp(idir,jdir)%elpt%wg(ixO^S);
        end where
       end do Loop_ivect2
      end if
     else cond_isone
      if(nozero)then
       Loop_ivect3add : do i=1,nwv
        where(patchw(ixO^S))
         wvcov(i)%w(ixO^S,iidir)=wvcov(i)%w(ixO^S,iidir)+&
                                 wvcont(i)%w(ixO^S,jjdir);
        end where
       end do Loop_ivect3add
      else
       Loop_ivect3 : do i=1,nwv
        where(patchw(ixO^S))
         wvcov(i)%w(ixO^S,iidir)=wvcont(i)%w(ixO^S,jjdir);
        end where
       end do Loop_ivect3
      end if 
     end if cond_isone
     nozero=.true.
   end if cond_ison
 end do Loop_jdir
end do Loop_allvec
end subroutine usemetric_contrtocov
!=============================================================================
subroutine usemetric_shiftspeed_contr(ixI^L,ixO^L,idir,patchw,vshift,vshift_out)
 integer, intent(in)             :: ixI^L,ixO^L,idir
 logical, intent(in)             :: patchw(ixG^T)
 double precision, intent(inout) :: vshift(ixG^T)
 double precision, optional      :: vshift_out(ixG^T)
! .. local variables ..
!----------------------------------------------------------------------------
 if(mypm%bt_cont(idir)%sftpt%elm_on) then
  if(mypm%alfa%elm_isone) then
   if(present(vshift_out)) then
    where(patchw(ixO^S))
     vshift_out(ixO^S)=vshift(ixO^S)-mypm%bt_cont(idir)%sftpt%wg(ixO^S)
    end where
   else 
    where(patchw(ixO^S))
     vshift(ixO^S)=vshift(ixO^S)-mypm%bt_cont(idir)%sftpt%wg(ixO^S)
    end where
   end if
  else
   if(present(vshift_out)) then
    where(patchw(ixO^S))
     vshift_out(ixO^S)=vshift(ixO^S)-mypm%bt_cont(idir)%sftpt%wg(ixO^S)&
                                 /mypm%alfa%wg(ixO^S)
    end where
   else
    where(patchw(ixO^S))
     vshift(ixO^S)=vshift(ixO^S)-mypm%bt_cont(idir)%sftpt%wg(ixO^S)&
                                 /mypm%alfa%wg(ixO^S)
    end where
   end if
  end if
 end if
end subroutine usemetric_shiftspeed_contr
!=============================================================================
subroutine usemetric_addshiftspeed_contr(ixI^L,ixO^L,idir,patchw,vshift,vout)
 integer, intent(in)                        :: ixI^L,ixO^L,idir
 logical, intent(in)                        :: patchw(ixG^T)
 double precision, intent(inout)            :: vshift(ixI^S)
 double precision, optional, intent(out)    :: vout(ixG^T)
! .. local variables ..
!----------------------------------------------------------------------------
if(present(vout))then
 if(mypm%bt_cont(idir)%sftpt%elm_on) then
  if(mypm%alfa%elm_isone) then
    where(patchw(ixO^S))
     vout(ixO^S)=vshift(ixO^S)+mypm%bt_cont(idir)%sftpt%wg(ixO^S)
    end where
  else
    where(patchw(ixO^S))
     vout(ixO^S)=vshift(ixO^S)+mypm%bt_cont(idir)%sftpt%wg(ixO^S)&
                                 /mypm%alfa%wg(ixO^S)
    end where
  end if
 else
  vout(ixO^S)=vshift(ixO^S)
 end if
else
 if(mypm%bt_cont(idir)%sftpt%elm_on) then
  if(mypm%alfa%elm_isone) then
    where(patchw(ixO^S))
     vshift(ixO^S)=vshift(ixO^S)+mypm%bt_cont(idir)%sftpt%wg(ixO^S)
    end where
  else
    where(patchw(ixO^S))
     vshift(ixO^S)=vshift(ixO^S)+mypm%bt_cont(idir)%sftpt%wg(ixO^S)&
                                 /mypm%alfa%wg(ixO^S)
    end where
  end if
 end if
end if
end subroutine usemetric_addshiftspeed_contr
!=============================================================================
subroutine usemetric_getgdd00(ixI^L,ixO^L,gdd00)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixI^L,ixO^L
 double precision                :: gdd00(ixI^S)
 ! .. local ..
 integer                         :: idir,jdir,jjdir
 double precision                :: beta_cov(ixI^S,setgr%ndir)
 logical                         :: beta_on(setgr%ndir)
!----------------------------------------------------------------------------
 if(mypm%alfa%elm_isone)then
  gdd00(ixO^S)=-1.0d0
 else
  gdd00(ixO^S)=-mypm%alfa%wg(ixO^S)**2.0d0
 end if
 call usemetric_betacov(ixI^L,ixO^L,1,patchtrue,beta_cov,beta_on,allvec=.true.)
 Loop_idir: do idir=1,^NC
  if(mypm%bt_cont(idir)%sftpt%elm_on.and.beta_on(idir))then
    if(mypm%bt_cont(idir)%sftpt%elm_isone) then
     gdd00(ixO^S)= gdd00(ixO^S)+beta_cov(ixO^S,idir)
    else
     gdd00(ixO^S)= gdd00(ixO^S)+&
       beta_cov(ixO^S,idir)&
       *mypm%bt_cont(idir)%sftpt%wg(ixO^S)
    end if
   end if
 end do Loop_idir
 
end subroutine usemetric_getgdd00
!=============================================================================
subroutine usemetric_getgddij(ixI^L,ixO^L,iidir,jjdir,gddij,nozero)
 use mod_geometry
 include 'amrvacdef.f'
 
 integer, intent(in)             :: ixI^L,ixO^L
 integer, intent(in)             :: iidir, jjdir
 double precision                :: gddij(ixI^S)
 logical, intent(out)            :: nozero
 ! .. local ..
 integer                         :: idir,jdir
!----------------------------------------------------------------------------
 idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
 cond_gii_on : if(mypm%elemsp(idir,jdir)%elpt%elm_on)then
  nozero=.true.
  cond_gii_isone : if(mypm%elemsp(idir,jdir)%elpt%elm_isone)then
    gddij(ixO^S)=1.0d0
  else cond_gii_isone
    gddij(ixO^S)=mypm%elemsp(idir,jdir)%elpt%wg(ixO^S)
  end if cond_gii_isone
 else cond_gii_on
  nozero=.false.
 end if cond_gii_on
end subroutine usemetric_getgddij
!======================================================================
subroutine usemetric_getgdd0i(ixI^L,ixO^L,idir,gdd0i,nozero)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixI^L,ixO^L
 integer, intent(in)             :: idir
 double precision                :: gdd0i(ixI^S)
 logical, intent(out)            :: nozero
 ! .. local ..
 logical                         :: beta_on(setgr%ndir)
 double precision                :: beta_cov(ixI^S,1:setgr%ndir)
!----------------------------------------------------------------------------
 call usemetric_betacov(ixI^L,ixO^L,idir,patchtrue,beta_cov,beta_on)
 cond_beta_on : if(beta_on(idir))then
   nozero=.true.
   gdd0i(ixO^S) = beta_cov(ixO^S,idir)
 else
  nozero=.false.
 end if cond_beta_on
end subroutine usemetric_getgdd0i
!=============================================================================
subroutine usemetric_getguu00(ixI^L,ixO^L,guu00)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixI^L,ixO^L
 double precision                :: guu00(ixI^S)
 ! .. local ..
 integer                         :: idir,jdir,jjdir
!----------------------------------------------------------------------------
 if(mypm%alfa%elm_isone)then
  guu00(ixO^S)=-1.0d0
 else
  guu00(ixO^S)=-mypm%alfa%wg(ixO^S)**(-2.0d0)
 end if
end subroutine usemetric_getguu00
!=============================================================================
subroutine usemetric_getguuij(ixI^L,ixO^L,iidir,jjdir,guuij,nozero)
 use mod_geometry
 include 'amrvacdef.f'

 integer, intent(in)             :: ixI^L,ixO^L
 integer, intent(in)             :: iidir, jjdir
 double precision                :: guuij(ixI^S)
 logical, intent(out)            :: nozero
 ! .. local ..
 integer                         :: idir,jdir
 double precision                :: bibj(ixI^S)
!----------------------------------------------------------------------------
 idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
 nozero=.false.
 
 cond_gii_on : if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on)then
  nozero=.true.
  cond_gii_isone : if(mypm%elemsp_inv(idir,jdir)%elpt%elm_isone)then
    guuij(ixO^S)=1.0d0
  else cond_gii_isone
    guuij(ixO^S)=mypm%elemsp_inv(idir,jdir)%elpt%wg(ixO^S)
  end if cond_gii_isone
 end if cond_gii_on 
 if(mypm%bt_cont(iidir)%sftpt%elm_on.and.&
    mypm%bt_cont(jjdir)%sftpt%elm_on)then
    nozero=.true.
    if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on)then
     if(mypm%bt_cont(iidir)%sftpt%elm_isone) then
      if(mypm%bt_cont(jjdir)%sftpt%elm_isone) then
       bibj(ixO^S) = 1.0D0
      else
       bibj(ixO^S)=mypm%bt_cont(jjdir)%sftpt%wg(ixO^S)
      end if      
     else
      if(mypm%bt_cont(jjdir)%sftpt%elm_isone) then
       bibj(ixO^S)=mypm%bt_cont(iidir)%sftpt%wg(ixO^S)
      else
       bibj(ixO^S)=mypm%bt_cont(iidir)%sftpt%wg(ixO^S)&
                                *mypm%bt_cont(jjdir)%sftpt%wg(ixO^S)
      end if
     end if
    end if
  if(mypm%alfa%elm_isone) then
   if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on)then
    guuij(ixO^S)=guuij(ixO^S)-bibj(ixO^S)
   else
    guuij(ixO^S)=-bibj(ixO^S)
   end if
  else
   if(mypm%elemsp_inv(idir,jdir)%elpt%elm_on)then
    guuij(ixO^S)=guuij(ixO^S)-bibj(ixO^S)/mypm%alfa%wg(ixO^S)
   else
    guuij(ixO^S)=-bibj(ixO^S)/mypm%alfa%wg(ixO^S)
   end if
  end if
 end if

end subroutine usemetric_getguuij
!======================================================================
subroutine usemetric_getguu0i(ixI^L,ixO^L,idir,guu0i,nozero)
 include 'amrvacdef.f'
 integer, intent(in)             :: ixI^L,ixO^L
 integer, intent(in)             :: idir
 double precision                :: guu0i(ixI^S)
 ! .. local ..
 logical, intent(out)            :: nozero
!----------------------------------------------------------------------------
 cond_beta_on : if(mypm%bt_cont(idir)%sftpt%elm_on)then
  nozero=.true.
  cond_beta_isone : if(.not.mypm%bt_cont(idir)%sftpt%elm_isone)then
   cond_alpha_isone : if(mypm%alfa%elm_isone)then
     guu0i(ixO^S) = mypm%bt_cont(idir)%sftpt%wg(ixO^S)
   else cond_alpha_isone
     guu0i(ixO^S) = mypm%bt_cont(idir)%sftpt%wg(ixO^S)/mypm%alfa%wg(ixO^S)
   end if cond_alpha_isone
  else cond_beta_isone
   cond_alpha_isone2 : if(mypm%alfa%elm_isone)then
    guu0i(ixO^S)  = one
   else cond_alpha_isone2
    guu0i(ixO^S) = one/mypm%alfa%wg(ixO^S)
   end if cond_alpha_isone2
  end if cond_beta_isone
 else
  nozero=.false.
 end if cond_beta_on
end subroutine usemetric_getguu0i

!======================================================================

subroutine usemetric_contrtocov_tensor_1(ixI^L,ixO^L,patchw,Tup,iidir,&
                                         ikeep,wdwn)
 ! T_(ik,i) to u_ik^i
 include 'amrvacdef.f'
 integer, intent(in)             :: ixI^L,ixO^L,ikeep,iidir
 logical, intent(in)             :: patchw(ixG^T)
 type(tensor_vf), intent(in)     :: Tup(1:^NC)
 double precision, intent(inout) :: wdwn(ixG^T)
! .. local variables ..
 integer                         :: jdir,kdir,idir,jjdir,kkdir
!----------------------------------------------------------------------------
 

if(all((/(.not.mypm%elemsp(i1_gr(iidir,jjdir),i2_gr(iidir,jjdir))&
                  %elpt%elm_on,jjdir=1,setgr%ndir)/)))then
 where(patchw(ixO^S))
  wdwn(ixO^S)=zero;
 end where
 return;
end if
Loop_jdir: do jjdir=1,setgr%ndir;idir=i1_gr(iidir,jjdir);jdir=i2_gr(iidir,jjdir)
 if(mypm%elemsp(idir,jdir)%elpt%elm_on)then
  if(.not.mypm%elemsp(idir,jdir)%elpt%elm_isone)then
    where(patchw(ixO^S))
       wdwn(ixO^S)=Tup(jjdir)%w(ixO^S)*&
                              mypm%elemsp(idir,jdir)%elpt%wg(ixO^S);
    end where
  else
    where(patchw(ixO^S))
      wdwn(ixO^S)=Tup(jjdir)%w(ixO^S);
    end where
  end if
     if(jjdir+1>setgr%ndir)exit Loop_jdir
      Loop_kdir : do kkdir=jjdir+1,setgr%ndir;idir=i1_gr(iidir,kkdir);
                                            kdir=i2_gr(iidir,kkdir)
       if(mypm%elemsp(idir,kdir)%elpt%elm_on)then
        if(.not.mypm%elemsp_inv(idir,kdir)%elpt%elm_isone)then
          where(patchw(ixO^S))
           wdwn(ixO^S)=wdwn(ixO^S)+&
             Tup(kkdir)%w(ixO^S)*mypm%elemsp(idir,kdir)%elpt%wg(ixO^S);
          end where
        else
          where(patchw(ixO^S))  
           wdwn(ixO^S)=wdwn(ixO^S)+Tup(kkdir)%w(ixO^S);
          end where
        end if
      
       end if
      end do Loop_kdir
      exit Loop_jdir
     end if
end do Loop_jdir
 

end subroutine usemetric_contrtocov_tensor_1
!=============================================================================
subroutine usemetric_addsource_ready(wFuu,wFdu)
 type(tensor_vF), intent(inout)     :: wFuu(:,:),wFdu(:,:)   
 ! .. local variables ..
 integer                            :: i,j,k,idir,jdir,kdir
 !----------------------------------------------------
wFdu(:,:)%elem_needed=.false.
wFuu(:,:)%elem_needed=.false.

Loop_j : do j=1,^NC
 Loop_i : do i=1,^NC
  idir=i1_gr(i,j);jdir=i2_gr(i,j)
  if(mypm%bt_cont(i)%sftpt%elm_on)then
   if(.not.mypm%bt_cont(i)%sftpt%elm_isone)then
    if(mypm%bt_cont(i)%sftpt%drv(j)%elm_on)then
    wFdu(i,j)%elem_needed=.true.
    end if
   end if
  end if
  Loop_k : do k=1,ndir;idir=i1_gr(i,k);kdir=i2_gr(i,k)
   if(mypm%elemsp(idir,kdir)%elpt%elm_on)then
    if(.not.mypm%elemsp(idir,kdir)%elpt%elm_isone)then
     if(any((/(mypm%elemsp(idir,kdir)%elpt%drv(j)%elm_on,j=1,setgr%ndir)/)))then
      wFuu(idir,kdir)%elem_needed=.true.
     end if
    end if
   end if
  end do Loop_k
 end do Loop_i
end do Loop_j
end subroutine usemetric_addsource_ready
!============================================================================
}
!===============================================================================
{^IFGR

subroutine usemetric_allg(ixI^L,ixO^L,x,patchw,G,alpha,beta)
use mod_geometry
include 'amrvacdef.f'
integer, intent(in)          :: ixI^L,ixO^L
double precision, intent(in) :: x(ixI^S,ndir)
logical, intent(in)          :: patchw(ixG^T)
type(sptmetric), intent(out) :: G(0:^NC,0:^NC)
double precision, intent(out):: alpha(ixI^S),beta(ixI^S,3)
! .. local ..
integer                    :: i,j,idir,jjdir,jdir
double precision           :: h(ixI^S)
type(vector_phys)          :: betacov(1)
!---------------------------------------
Loop_ibeta : do i=1,3
 if(mypm%bt_cont(i)%sftpt%elm_on)then
  if(mypm%bt_cont(i)%sftpt%elm_isone)then
   beta(ixO^S,i) = one
  else
   beta(ixO^S,i) = mypm%bt_cont(i)%sftpt%wg(ixO^S)
  end if
 else
  beta(ixO^S,i) =zero
 end if
end do Loop_ibeta
Loop_ibeta_cov : do i=1,3
 call usemetric_contrtocov(ixI^L,ixO^L,1,(/0/),i,patchw,3,beta,betacov)
 call getcoordinatehi(ixI^L,ixO^L,i,x,h)
 where(patchw(ixO^S))
  G(0,i)%wg(ixO^S)=betacov(1)%w(ixO^S,i)*h(ixO^S)
  G(i,0)%wg(ixO^S)= G(0,i)%wg(ixO^S)
 end where
end do Loop_ibeta_cov
if(mypm%alfa%elm_isone)then
 where(patchw(ixO^S))
  G(0,0)%wg(ixO^S)=-1.0d0
  alpha(ixO^S) =1.0D0
 end where
else
 where(patchw(ixO^S))
  G(0,0)%wg(ixO^S)=-mypm%alfa%wg(ixO^S)**2.0d0
  alpha(ixO^S) =mypm%alfa%wg(ixO^S)
 end where
end if
 Loop_idir: do idir=1,3
  where(patchw(ixO^S))
   G(0,0)%wg(ixO^S)= G(0,0)%wg(ixO^S)+beta(ixO^S,idir)*betacov(1)%w(ixO^S,idir)
  end where
 end do Loop_idir
Loop_i :do i=1,3
 if(mypm%elemsp(i,i)%elpt%elm_on)then
  call getcoordinatehi(ixI^L,ixO^L,i,x,h)
  where(patchw(ixO^S))
   G(i,i)%wg(ixO^S)=h(ixO^S)**2.0d0
  end where
  if(.not.mypm%elemsp(i,i)%elpt%elm_isone)then
   where(patchw(ixO^S))
    G(i,i)%wg(ixO^S)=G(i,i)%wg(ixO^S)*mypm%elemsp(i,i)%elpt%wg(ixO^S)
   end where
  end if
 else
  where(patchw(ixO^S))
   G(i,i)%wg(ixO^S)=zero
  end where
 end if
 if(i<3)then
  Loop_j :do j=i+1,3
   if(mypm%elemsp(i,j)%elpt%elm_on)then
     call getcoordinatehi(ixI^L,ixO^L,i,x,h)
     where(patchw(ixO^S))
      G(i,j)%wg(ixO^S)=h(ixO^S)
     end where
     call getcoordinatehi(ixI^L,ixO^L,j,x,h)
     where(patchw(ixO^S))
      G(i,j)%wg(ixO^S)=G(i,j)%wg(ixO^S)*h(ixO^S)
     end where
    if(.not.mypm%elemsp(i,j)%elpt%elm_isone)then
      where(patchw(ixO^S))
       G(i,j)%wg(ixO^S)=G(i,j)%wg(ixO^S)*mypm%elemsp(i,j)%elpt%wg(ixO^S)
      end where
    end if
   else
    where(patchw(ixO^S))
     G(i,j)%wg(ixO^S)=zero
    end where
   end if
  end do Loop_j
 end if

end do Loop_i
do i=1,ndir
 call getcoordinatehi(ixI^L,ixO^L,i,x,h)
 beta(ixO^S,i) = beta(ixO^S,i)/h(ixO^S)
end do
end subroutine usemetric_allg
!=========================================================================
subroutine usemetric_allg_inv(ixI^L,ixO^L,x,patchw,G_inv)
use mod_geometry
include 'amrvacdef.f'
integer, intent(in)          :: ixI^L,ixO^L
double precision, intent(in) :: x(ixI^S,ndir)
logical, intent(in)          :: patchw(ixG^T)
type(sptmetric), intent(out) :: G_inv(0:^NC,0:^NC)

! .. local ..
integer                    :: i,j,idir,jjdir,jdir
double precision           :: beta(ixI^S,3),h(ixI^S)
type(vector_phys)          :: betacov(1)
!---------------------------------------
Loop_ibeta : do i=1,3
 if(mypm%bt_cont(i)%sftpt%elm_on)then
  if(mypm%bt_cont(i)%sftpt%elm_isone)then
   beta(ixO^S,i) = one
  else
   beta(ixO^S,i) = mypm%bt_cont(i)%sftpt%wg(ixO^S)
  end if
 else
  beta(ixO^S,i) =zero
 end if
end do Loop_ibeta
Loop_ibeta_cov : do i=1,3
 call usemetric_contrtocov(ixI^L,ixO^L,1,(/0/),i,patchw,3,beta,betacov)
 call getcoordinatehi(ixI^L,ixO^L,i,x,h)
 G_inv(0,i)%wg(ixO^S)=beta(ixO^S,i)/h(ixO^S)
 if(.not.mypm%alfa%elm_isone)then
  G_inv(0,i)%wg(ixO^S)=G_inv(0,i)%wg(ixO^S)/mypm%alfa%wg(ixO^S)**2.0d0
 end if
 G_inv(i,0)%wg(ixO^S)=G_inv(0,i)%wg(ixO^S)
end do Loop_ibeta_cov
if(mypm%alfa%elm_isone)then
 G_inv(0,0)%wg(ixO^S)=-1.0d0
else
 G_inv(0,0)%wg(ixO^S)=-mypm%alfa%wg(ixO^S)**(-2.0d0)
end if
Loop_i :do i=1,3
 if(mypm%elemsp_inv(i,i)%elpt%elm_on)then
  call getcoordinatehi(ixI^L,ixO^L,i,x,h)
  G_inv(i,i)%wg(ixO^S)=h(ixO^S)**(-2.0d0)
  if(.not.mypm%elemsp_inv(i,i)%elpt%elm_isone)then
   G_inv(i,i)%wg(ixO^S)=G_inv(i,i)%wg(ixO^S)*mypm%elemsp_inv(i,i)%elpt%wg(ixO^S)
  end if
 else
  G_inv(i,i)%wg(ixO^S)=zero
 end if
 if(mypm%alfa%elm_isone)then
  if(mypm%bt_cont(i)%sftpt%elm_on)then
   G_inv(i,i)%wg(ixO^S)=G_inv(i,i)%wg(ixO^S)-beta(ixO^S,i)*betacov(1)%w(ixO^S,i)
  end if
 else
  if(mypm%bt_cont(i)%sftpt%elm_on)then
    G_inv(i,i)%wg(ixO^S)=G_inv(i,i)%wg(ixO^S)-&
               beta(ixO^S,i)*betacov(1)%w(ixO^S,i)/mypm%alfa%wg(ixO^S)**2.0D0
  end if
 end if
 if(i<3)then
  Loop_j :do j=i+1,3
   if(mypm%elemsp_inv(i,j)%elpt%elm_on)then
     call getcoordinatehi(ixI^L,ixO^L,i,x,h)
     G_inv(i,j)%wg(ixO^S)=1.0d0/h(ixO^S)
     call getcoordinatehi(ixI^L,ixO^L,j,x,h)
     G_inv(i,j)%wg(ixO^S)=G_inv(i,j)%wg(ixO^S)/h(ixO^S)
    if(.not.mypm%elemsp_inv(i,j)%elpt%elm_isone)then
      G_inv(i,j)%wg(ixO^S)=G_inv(i,j)%wg(ixO^S)*mypm%elemsp_inv(i,j)%elpt%wg(ixO^S)
    end if
   else
    G_inv(i,j)%wg(ixO^S)=zero
   end if
   if(mypm%alfa%elm_isone)then
    if(mypm%bt_cont(i)%sftpt%elm_on)then
     G_inv(i,j)%wg(ixO^S)=G_inv(i,j)%wg(ixO^S)-beta(ixO^S,i)*betacov(1)%w(ixO^S,j)
    end if
   else
    if(mypm%bt_cont(i)%sftpt%elm_on)then
     G_inv(i,j)%wg(ixO^S)=G_inv(i,j)%wg(ixO^S)-&
               beta(ixO^S,i)*betacov(1)%w(ixO^S,j)/mypm%alfa%wg(ixO^S)**2.0D0
    end if
   end if
  end do Loop_j
 end if
end do Loop_i
end subroutine usemetric_allg_inv
!=============================================================================
subroutine usemetric_betacov(ixI^L,ixO^L,idir_out,patchw,beta_cov,betacov_on,allvec)
! w should be covariant

 integer, intent(in)             :: ixI^L,ixO^L,idir_out
 logical, intent(in)             :: patchw(ixG^T)
 double precision, intent(out)   :: beta_cov(ixI^S,setgr%ndir)
 logical, intent(out)            :: betacov_on(setgr%ndir)
 logical, optional               :: allvec
 ! .. local variables ..
 logical                         :: nozero
 integer                         :: idir,jdir,iidir,jjdir
 integer                         :: idir_min,idir_max
!-----------------------------------------------------------------------------
 if(present(allvec))then
  if(allvec)then
   idir_min = 1
   idir_max = setgr%ndir
  else
   idir_min =idir_out;idir_max = idir_out
  end if
 else
  idir_min =idir_out;idir_max = idir_out
 end if
 betacov_on=.false.
 Loop_idir : do iidir=idir_min,idir_max
   nozero=.false.
   if(all((/((.not.mypm%elemsp(i1_gr(iidir,jdir),i2_gr(iidir,jdir))%elpt%elm_on&
             .or..not.mypm%bt_cont(jdir)%sftpt%elm_on)&
             ,jdir=1,setgr%ndir)/)))then
    beta_cov(ixO^S,iidir)=zero;
    cycle Loop_idir
   else
    betacov_on(iidir)=.true.
   end if
 Loop_jdir : do jjdir=1,setgr%ndir;idir=i1_gr(iidir,jjdir);
                                    jdir=i2_gr(iidir,jjdir)
   if(mypm%elemsp(idir,jdir)%elpt%elm_on&
      .and.mypm%bt_cont(jjdir)%sftpt%elm_on)then
     if(.not.mypm%elemsp(idir,jdir)%elpt%elm_isone)then

       if(.not.mypm%bt_cont(jjdir)%sftpt%elm_isone)then
        if(nozero)then
         where(patchw(ixO^S))
          beta_cov(ixO^S,iidir)=beta_cov(ixO^S,iidir)+&
                            mypm%bt_cont(jjdir)%sftpt%wg(ixO^S)&
                            *mypm%elemsp(idir,jdir)%elpt%wg(ixO^S)
         end where
        else
         where(patchw(ixO^S))
          beta_cov(ixO^S,iidir)=mypm%bt_cont(jjdir)%sftpt%wg(ixO^S)&
                            *mypm%elemsp(idir,jdir)%elpt%wg(ixO^S)
         end where
        end if
       else
        if(nozero)then
         where(patchw(ixO^S))
          beta_cov(ixO^S,iidir)=beta_cov(ixO^S,iidir)+&
                                mypm%elemsp(idir,jdir)%elpt%wg(ixO^S)
         end where
        else
         where(patchw(ixO^S))
          beta_cov(ixO^S,iidir)=mypm%elemsp(idir,jdir)%elpt%wg(ixO^S)
         end where
        end if
       end if
     else
       if(.not.mypm%bt_cont(jjdir)%sftpt%elm_isone)then
        if(nozero)then
         where(patchw(ixO^S))
          beta_cov(ixO^S,iidir)=beta_cov(ixO^S,iidir)+&
                                mypm%bt_cont(jjdir)%sftpt%wg(ixO^S)
         end where
        else
         where(patchw(ixO^S))
          beta_cov(ixO^S,iidir)=mypm%bt_cont(jjdir)%sftpt%wg(ixO^S)
         end where
        end if
       else
        if(nozero)then
         where(patchw(ixO^S))
          beta_cov(ixO^S,iidir)=beta_cov(ixO^S,iidir)+1.0D0
         end where
        else
         where(patchw(ixO^S))
          beta_cov(ixO^S,iidir)=1.0D0
         end where
        end if
       end if

     end if
    else
     if(.not.nozero)then
      where(patchw(ixO^S))
       beta_cov(ixO^S,iidir)=zero
      end where
     end if
    end if
    nozero=.true.
    cycle
  end do Loop_jdir
 end do Loop_idir
end subroutine usemetric_betacov
!===========================================================================
subroutine usemetric_getGabCabJ(ixI^L,ixO^L,jdir,patchw,x,GabCabJ,Gnozero)
 use mod_geometry
 include 'amrvacdef.f'
 integer, intent(in)           :: ixI^L,ixO^L,jdir
 logical, intent(in)           :: patchw(ixG^T)
 double precision, intent(in)  :: x(ixI^S,1:ndim)
 double precision, intent(out) :: GabCabJ(ixI^S)
 logical, intent(out)          :: Gnozero
 ! .. local ..
 logical                       :: nozero
 integer                       :: idims, ixC^L,hxO^L
 double precision              :: fC(ixG^T),guuij(ixG^T),hC(ixG^T)
!----------------------------------------------------
Loop_idim :  do idims=1,ndim
 hxO^L=ixO^L-kr(idims,^D);
 ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
 call usemetric_getguuij(ixI^L,ixC^L,idims,jdir,guuij,nozero)
 if(nozero)then
  call usemetric_multiby_sqrtspacedetr(ixC^L,patchw,guuij,fC)
  call usemetric_multiby_alpha(ixC^L,fC)
  select case (idims)
         {case (^D)
          call getcoordinatehiC(ixI^L,ixC^L,jdir,x,hC)
          if(Gnozero) then
           GabCabJ(ixO^S) = GabCabJ(ixO^S) + &
                           (fC(ixO^S) * mygeo%surfaceC^D(ixO^S)/hC(ixO^S)&
                          - fC(hxO^S) * mygeo%surfaceC^D(hxO^S)/hC(hxO^S))&
                          / mygeo%dvolume(ixO^S)
          else
           GabCabJ(ixO^S) = (fC(ixO^S) * mygeo%surfaceC^D(ixO^S)/hC(ixO^S)&
                          - fC(hxO^S) * mygeo%surfaceC^D(hxO^S)/hC(hxO^S))&
                          / mygeo%dvolume(ixO^S)
          end if
         \}
  end select
  Gnozero = .true.
 end if
end do Loop_idim
if(Gnozero)call usemetric_multiby_alpha(ixO^L,GabCabJ)
end subroutine usemetric_getGabCabJ
!=========================================================================
subroutine usemetric_getGabCab0(ixI^L,ixO^L,patchw,x,GabCab0,Gnozero)
 use mod_geometry
 include 'amrvacdef.f'
 integer, intent(in)           :: ixI^L,ixO^L
 logical, intent(in)           :: patchw(ixG^T)
 double precision, intent(in)  :: x(ixI^S,1:ndim)
 double precision, intent(out) :: GabCab0(ixI^S)
 logical, intent(out)          :: Gnozero
 ! .. local ..
 logical                       :: nozero
 integer                       :: idims, ixC^L,hxO^L
 double precision              :: fC(ixG^T),guu0i(ixG^T),hC(ixG^T)
!----------------------------------------------------
Loop_idim :  do idims=1,ndim
 hxO^L=ixO^L-kr(idims,^D);
 ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
 call usemetric_getguu0i(ixI^L,ixC^L,idims,guu0i,nozero)
 if(nozero)then
  call usemetric_multiby_sqrtspacedetr(ixC^L,patchw,guu0i,fC)
  call usemetric_multiby_alpha(ixC^L,fC)
  select case (idims)
         {case (^D)
          if(Gnozero) then
           GabCab0(ixO^S) = GabCab0(ixO^S) + &
                           (fC(ixO^S) * mygeo%surfaceC^D(ixO^S)&
                          - fC(hxO^S) * mygeo%surfaceC^D(hxO^S))&
                          / mygeo%dvolume(ixO^S)
          else
           GabCab0(ixO^S) = (fC(ixO^S) * mygeo%surfaceC^D(ixO^S)&
                          - fC(hxO^S) * mygeo%surfaceC^D(hxO^S))&
                          / mygeo%dvolume(ixO^S)
          end if
         \}
  end select
  Gnozero = .true.
 end if

end do Loop_idim
if(Gnozero)call usemetric_multiby_alpha(ixO^L,GabCab0)
end subroutine usemetric_getGabCab0
}
!============================================================================
End module mod_usemetric 
!============================================================================
