!============================================================================
subroutine eigenvectors(idims,ixI^L,ixO^L,w,leftm,righm,lamda)
! calculate eigenvectors and eigenvalues for flux Jacobian matrix of x_idims
! (Jiang, Chaowei et al. 2011, ApJ, 727, 101)
! Chun Xia 
! 13 May 2016
use mod_global_parameters

integer, intent(in) :: idims,ixI^L,ixO^L
double precision, intent(in) :: w(ixI^S,1:nw)
double precision :: leftm(ixI^S,1:nwwave,1:nwprim)
double precision :: righm(ixI^S,1:nwprim,1:nwwave)
double precision :: lamda(ixI^S,1:nwwave), betad(ixI^S,1:ndir)
double precision, dimension(ixI^S) :: a2,betau,signbi,va,vf,vs,af,as,tmp
double precision :: sqrthalf
integer :: idirs, ii,is(ndir-1),si(ndir-1)

!----------------------------------------------------------------------------
sqrthalf=dsqrt(0.5d0)
{#IFDEF ENERGY
! sound speed squared
a2(ixO^S)=eqpar(gamma_)*w(ixO^S,p_)/w(ixO^S,rho_)
}
{#IFDEF ISO
a2(ixO^S)=eqpar(gamma_)*eqpar(adiab_)*w(ixO^S,rho_)**(eqpar(gamma_)-one)
}
betau(ixO^S)=dsqrt((^C&w(ixO^S,b^C_)**2+ )-w(ixO^S,b0_+idims)**2) 
ii=0
do idirs=1,ndir
  if(idirs==idims) cycle 
  ii=ii+1
  is(ii)=idirs
  where(betau(ixO^S)<smalldouble)
    betad(ixO^S,idirs)=sqrthalf
  elsewhere
    betad(ixO^S,idirs)=w(ixO^S,b0_+idirs)/betau(ixO^S)
  end where
end do
si(1)=is(2)
si(2)=is(1)
signbi(ixO^S)=sign(1.d0,w(ixO^S,b0_+idims))
! idims dimensional Alfven speed 
va(ixO^S)=dabs(w(ixO^S,b0_+idims))/dsqrt(w(ixO^S,rho_))
vf(ixO^S)=a2(ixO^S)+ (^C&w(ixO^S,b^C_)**2+ )/w(ixO^S,rho_)
tmp(ixO^S)=dsqrt(vf(ixO^S)**2-4.d0*a2(ixO^S)*va(ixO^S)**2)
vs(ixO^S)=0.5d0*(vf(ixO^S)-tmp(ixO^S))
vf(ixO^S)=0.5d0*(vf(ixO^S)+tmp(ixO^S))
tmp(ixO^S)=vf(ixO^S)-vs(ixO^S)
af(ixO^S)=a2(ixO^S)-vs(ixO^S)
where(af(ixO^S)<0.d0)
  af(ixO^S)=0.d0
end where
af(ixO^S)=dsqrt(af(ixO^S)/tmp(ixO^S))
as(ixO^S)=vf(ixO^S)-a2(ixO^S)
where(as(ixO^S)<0.d0)
  as(ixO^S)=0.d0
end where
as(ixO^S)=dsqrt(as(ixO^S)/tmp(ixO^S))
! fast magneto-acoustic wave speed 
vf(ixO^S)=dsqrt(vf(ixO^S))
! slow magneto-acoustic wave speed 
vs(ixO^S)=dsqrt(vs(ixO^S))
leftm(ixO^S,:,:)=0.d0
righm(ixO^S,:,:)=0.d0
{#IFDEF ENERGY
! entropy wave
lamda(ixO^S,entroW_)=w(ixO^S,v0_+idims)
righm(ixO^S,rho_,entroW_)=1.d0
leftm(ixO^S,entroW_,rho_)=1.d0
leftm(ixO^S,entroW_,  p_)=-1.d0/a2(ixO^S)
}
! magnetic flux wave
lamda(ixO^S,diverW_)=w(ixO^S,v0_+idims)
righm(ixO^S,b0_+idims,diverW_)=1.d0
leftm(ixO^S,diverW_,b0_+idims)=1.d0
! alfven wave down-stream
lamda(ixO^S,alfvRW_)=w(ixO^S,v0_+idims)+va(ixO^S)
righm(ixO^S,rho_+is(1),alfvRW_)= sqrthalf*betad(ixO^S,si(1))
righm(ixO^S,rho_+is(2),alfvRW_)=-sqrthalf*betad(ixO^S,si(2))
righm(ixO^S, b0_+is(1),alfvRW_)=-sqrthalf*betad(ixO^S,si(1))*dsqrt(w(ixO^S,rho_))*signbi(ixO^S)
righm(ixO^S, b0_+is(2),alfvRW_)= sqrthalf*betad(ixO^S,si(2))*dsqrt(w(ixO^S,rho_))*signbi(ixO^S)
leftm(ixO^S,alfvRW_,rho_+is(1))=righm(ixO^S,rho_+is(1),alfvRW_)
leftm(ixO^S,alfvRW_,rho_+is(2))=righm(ixO^S,rho_+is(2),alfvRW_)
leftm(ixO^S,alfvRW_, b0_+is(1))=righm(ixO^S, b0_+is(1),alfvRW_)/w(ixO^S,rho_)
leftm(ixO^S,alfvRW_, b0_+is(2))=righm(ixO^S, b0_+is(2),alfvRW_)/w(ixO^S,rho_)
! alfven wave up-stream
lamda(ixO^S,alfvLW_)=w(ixO^S,v0_+idims)-va(ixO^S)
righm(ixO^S,rho_+is(1),alfvLW_)=-righm(ixO^S,rho_+is(1),alfvRW_)
righm(ixO^S,rho_+is(2),alfvLW_)=-righm(ixO^S,rho_+is(2),alfvRW_)
righm(ixO^S, b0_+is(1),alfvLW_)= righm(ixO^S, b0_+is(1),alfvRW_)
righm(ixO^S, b0_+is(2),alfvLW_)= righm(ixO^S, b0_+is(2),alfvRW_)
leftm(ixO^S,alfvLW_,rho_+is(1))=-leftm(ixO^S,alfvRW_,rho_+is(1))
leftm(ixO^S,alfvLW_,rho_+is(2))=-leftm(ixO^S,alfvRW_,rho_+is(2))
leftm(ixO^S,alfvLW_, b0_+is(1))= leftm(ixO^S,alfvRW_, b0_+is(1))
leftm(ixO^S,alfvLW_, b0_+is(2))= leftm(ixO^S,alfvRW_, b0_+is(2))
! fast wave down-stream
lamda(ixO^S,fastRW_)=w(ixO^S,v0_+idims)+vf(ixO^S)
righm(ixO^S,rho_      ,fastRW_)=w(ixO^S,rho_)*af(ixO^S)
righm(ixO^S,rho_+idims,fastRW_)=af(ixO^S)*vf(ixO^S)
righm(ixO^S,rho_+is(1),fastRW_)=-as(ixO^S)*vs(ixO^S)*signbi(ixO^S)
righm(ixO^S,rho_+is(2),fastRW_)=righm(ixO^S,rho_+is(1),fastRW_)*betad(ixO^S,is(2))
righm(ixO^S,rho_+is(1),fastRW_)=righm(ixO^S,rho_+is(1),fastRW_)*betad(ixO^S,is(1))
{#IFDEF ENERGY
righm(ixO^S,  p_      ,fastRW_)=af(ixO^S)*eqpar(gamma_)*w(ixO^S,p_)
}
righm(ixO^S, b0_+is(1),fastRW_)=as(ixO^S)*dsqrt(w(ixO^S,rho_)*a2(ixO^S))
righm(ixO^S, b0_+is(2),fastRW_)=righm(ixO^S, b0_+is(1),fastRW_)*betad(ixO^S,is(2))
righm(ixO^S, b0_+is(1),fastRW_)=righm(ixO^S, b0_+is(1),fastRW_)*betad(ixO^S,is(1))
leftm(ixO^S,fastRW_,rho_+idims)=righm(ixO^S,rho_+idims,fastRW_)/(2.d0*a2(ixO^S))
leftm(ixO^S,fastRW_,rho_+is(1))=righm(ixO^S,rho_+is(1),fastRW_)/(2.d0*a2(ixO^S))
leftm(ixO^S,fastRW_,rho_+is(2))=righm(ixO^S,rho_+is(2),fastRW_)/(2.d0*a2(ixO^S))
{#IFDEF ENERGY
leftm(ixO^S,fastRW_,  p_      )=af(ixO^S)/(2.d0*a2(ixO^S)*w(ixO^S,rho_))
}
leftm(ixO^S,fastRW_, b0_+is(1))=righm(ixO^S, b0_+is(1),fastRW_)/(2.d0*a2(ixO^S)*w(ixO^S,rho_))
leftm(ixO^S,fastRW_, b0_+is(2))=righm(ixO^S, b0_+is(2),fastRW_)/(2.d0*a2(ixO^S)*w(ixO^S,rho_))
! fast wave up-stream
lamda(ixO^S,fastLW_)=w(ixO^S,v0_+idims)-vf(ixO^S)
righm(ixO^S,rho_      ,fastLW_)= righm(ixO^S,rho_      ,fastRW_)
righm(ixO^S,rho_+idims,fastLW_)=-righm(ixO^S,rho_+idims,fastRW_)
righm(ixO^S,rho_+is(1),fastLW_)=-righm(ixO^S,rho_+is(1),fastRW_)
righm(ixO^S,rho_+is(2),fastLW_)=-righm(ixO^S,rho_+is(2),fastRW_)
{#IFDEF ENERGY
righm(ixO^S,  p_      ,fastLW_)= righm(ixO^S,  p_      ,fastRW_)
}
righm(ixO^S, b0_+is(1),fastLW_)= righm(ixO^S, b0_+is(1),fastRW_)
righm(ixO^S, b0_+is(2),fastLW_)= righm(ixO^S, b0_+is(2),fastRW_)
leftm(ixO^S,fastLW_,rho_      )= leftm(ixO^S,fastRW_,rho_      )
leftm(ixO^S,fastLW_,rho_+idims)=-leftm(ixO^S,fastRW_,rho_+idims)
leftm(ixO^S,fastLW_,rho_+is(1))=-leftm(ixO^S,fastRW_,rho_+is(1))
leftm(ixO^S,fastLW_,rho_+is(2))=-leftm(ixO^S,fastRW_,rho_+is(2))
{#IFDEF ENERGY
leftm(ixO^S,fastLW_,  p_      )= leftm(ixO^S,fastRW_,  p_      )
}
leftm(ixO^S,fastLW_, b0_+is(1))= leftm(ixO^S,fastRW_, b0_+is(1))
leftm(ixO^S,fastLW_, b0_+is(2))= leftm(ixO^S,fastRW_, b0_+is(2))
! slow wave down-stream
lamda(ixO^S,slowRW_)=w(ixO^S,v0_+idims)+vs(ixO^S)
righm(ixO^S,rho_      ,slowRW_)=w(ixO^S,rho_)*as(ixO^S)
righm(ixO^S,rho_+idims,slowRW_)=as(ixO^S)*vs(ixO^S)
righm(ixO^S,rho_+is(1),slowRW_)=af(ixO^S)*vf(ixO^S)*signbi(ixO^S)
righm(ixO^S,rho_+is(2),slowRW_)=righm(ixO^S,rho_+is(1),slowRW_)*betad(ixO^S,is(2))
righm(ixO^S,rho_+is(1),slowRW_)=righm(ixO^S,rho_+is(1),slowRW_)*betad(ixO^S,is(1))
{#IFDEF ENERGY
righm(ixO^S,  p_      ,slowRW_)=as(ixO^S)*eqpar(gamma_)*w(ixO^S,p_)
}
righm(ixO^S, b0_+is(1),slowRW_)=-af(ixO^S)*dsqrt(w(ixO^S,rho_)*a2(ixO^S))
righm(ixO^S, b0_+is(2),slowRW_)=righm(ixO^S, b0_+is(1),slowRW_)*betad(ixO^S,is(2))
righm(ixO^S, b0_+is(1),slowRW_)=righm(ixO^S, b0_+is(1),slowRW_)*betad(ixO^S,is(1))
leftm(ixO^S,slowRW_,rho_+idims)=righm(ixO^S,rho_+idims,slowRW_)/(2.d0*a2(ixO^S))
leftm(ixO^S,slowRW_,rho_+is(1))=righm(ixO^S,rho_+is(1),slowRW_)/(2.d0*a2(ixO^S))
leftm(ixO^S,slowRW_,rho_+is(2))=righm(ixO^S,rho_+is(2),slowRW_)/(2.d0*a2(ixO^S))
{#IFDEF ENERGY
leftm(ixO^S,slowRW_,  p_      )=as(ixO^S)/(2.d0*a2(ixO^S)*w(ixO^S,rho_))
}
leftm(ixO^S,slowRW_, b0_+is(1))=righm(ixO^S, b0_+is(1),slowRW_)/(2.d0*a2(ixO^S)*w(ixO^S,rho_))
leftm(ixO^S,slowRW_, b0_+is(2))=righm(ixO^S, b0_+is(2),slowRW_)/(2.d0*a2(ixO^S)*w(ixO^S,rho_))
! slow wave up-stream
lamda(ixO^S,slowLW_)=w(ixO^S,v0_+idims)-vs(ixO^S)
righm(ixO^S,rho_      ,slowLW_)= righm(ixO^S,rho_      ,slowRW_)
righm(ixO^S,rho_+idims,slowLW_)=-righm(ixO^S,rho_+idims,slowRW_)
righm(ixO^S,rho_+is(1),slowLW_)=-righm(ixO^S,rho_+is(1),slowRW_)
righm(ixO^S,rho_+is(2),slowLW_)=-righm(ixO^S,rho_+is(2),slowRW_)
{#IFDEF ENERGY
righm(ixO^S,  p_      ,slowLW_)= righm(ixO^S,  p_      ,slowRW_)
}
righm(ixO^S, b0_+is(1),slowLW_)= righm(ixO^S, b0_+is(1),slowRW_)
righm(ixO^S, b0_+is(2),slowLW_)= righm(ixO^S, b0_+is(2),slowRW_)
leftm(ixO^S,slowLW_,rho_      )= leftm(ixO^S,slowRW_,rho_      )
leftm(ixO^S,slowLW_,rho_+idims)=-leftm(ixO^S,slowRW_,rho_+idims)
leftm(ixO^S,slowLW_,rho_+is(1))=-leftm(ixO^S,slowRW_,rho_+is(1))
leftm(ixO^S,slowLW_,rho_+is(2))=-leftm(ixO^S,slowRW_,rho_+is(2))
{#IFDEF ENERGY
leftm(ixO^S,slowLW_,  p_      )= leftm(ixO^S,slowRW_,  p_      )
}
leftm(ixO^S,slowLW_, b0_+is(1))= leftm(ixO^S,slowRW_, b0_+is(1))
leftm(ixO^S,slowLW_, b0_+is(2))= leftm(ixO^S,slowRW_, b0_+is(2))

end subroutine eigenvectors
!============================================================================
subroutine characteristic_project(idims,iside,ixI^L,ixO^L,w,x,dxndim,qdt)
! prepare time dependent boundary conditions using projected-characteristics
! method (Jiang, Chaowei et al. 2011, ApJ, 727, 101)
! Chun Xia
! 13 May 2016
use mod_global_parameters

integer, intent(in) :: idims, iside, ixI^L, ixO^L
double precision :: w(ixI^S,1:nw),x(ixI^S,1:ndim),dxndim(1:ndim),qdt
double precision :: jcobi(ixI^S,1:nwprim,1:nwprim)
double precision :: leftm(ixI^S,1:nwwave,1:nwprim)
double precision :: righm(ixI^S,1:nwprim,1:nwwave)
double precision :: ritmp(ixI^S,1:nwprim,1:nwwave)
double precision, dimension(ixI^S,1:nwprim) :: Sx, dw
double precision, dimension(ixI^S,1:nwwave) :: lamda, RHS
double precision :: invdx(ndir), dxall(ndir)
double precision :: divb(ixI^S)
integer :: idirs, ix^D, ixR^L,iwave,iprim,ii,is(ndir-1)
logical :: patchw(ixI^S)
!----------------------------------------------------------------------------
dxall(1:ndim)=dxndim(1:ndim)
if(ndir>ndim) dxall(ndim+1:ndir)=0.d0
invdx(1:ndim)=1.d0/dxall(1:ndim)
!print*,'cp ok 0'
call primitive(ixI^L,ixI^L,w,x)
!print*,'cp ok 1'
ii=0
do idirs=1,ndir
  if(idirs==idims) cycle
  ii=ii+1
  is(ii)=idirs
end do
Sx=0.d0
do ii=1,ndir-1
  idirs=is(ii)
  if(dxall(idirs)==0) cycle
  call eigenvectors(idirs,ixI^L,ixO^L,w,leftm,righm,lamda)
  ! negative part of flux Jacobian matrix
  ixR^L=ixO^L+kr(idirs,^D);
  dw(ixO^S,1:nwprim)=(w(ixR^S,1:nwprim)-w(ixO^S,1:nwprim))*invdx(idirs)
  lamda(ixO^S,1:nwwave)=0.5d0*(lamda(ixO^S,1:nwwave)-dabs(lamda(ixO^S,1:nwwave)))
  do iwave=1,nwwave
    do iprim=1,nwprim
      ritmp(ixO^S,iprim,iwave)=lamda(ixO^S,iwave)*righm(ixO^S,iprim,iwave)
    end do
  end do
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
    jcobi(ix^D,1:nwprim,1:nwprim)=matmul(ritmp(ix^D,1:nwprim,1:nwwave),leftm(ix^D,1:nwwave,1:nwprim))
  {end do\}
  do iprim=1,nwprim
    Sx(ixO^S,iprim)=Sx(ixO^S,iprim)-sum(jcobi(ixO^S,1:nwprim,iprim)*dw(ixO^S,1:nwprim),^ND+1)
  end do
  ! positive part of flux Jacobian matrix
  ixR^L=ixO^L-kr(idirs,^D);
  dw(ixO^S,1:nwprim)=(w(ixO^S,1:nwprim)-w(ixR^S,1:nwprim))*invdx(idirs)
  lamda(ixO^S,1:nwwave)=0.5d0*(lamda(ixO^S,1:nwwave)+dabs(lamda(ixO^S,1:nwwave)))
  do iwave=1,nwwave
    do iprim=1,nwprim
      ritmp(ixO^S,iprim,iwave)=lamda(ixO^S,iwave)*righm(ixO^S,iprim,iwave)
    end do
  end do
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
    jcobi(ix^D,1:nwprim,1:nwprim)=matmul(ritmp(ix^D,1:nwprim,1:nwwave),leftm(ix^D,1:nwwave,1:nwprim))
  {end do\}
  do iprim=1,nwprim
    Sx(ixO^S,iprim)=Sx(ixO^S,iprim)+sum(jcobi(ixO^S,1:nwprim,iprim)*dw(ixO^S,1:nwprim),^ND+1)
  end do
end do
!if(iside==1) then
!  ixR^L=ixO^L;
!  ixRmin^D=ixOmin^D+kr(idims,^D)\
!  call getdivb(w,ixI^L,ixR^L,divb)
!  select case (idims)
!  {case (^D)
!     divb(ixOmin^D^D%ixO^S) = divb(ixOmin^D+1^D%ixO^S) \}
!  end select
!else
!  ixR^L=ixO^L;
!  ixRmax^D=ixOmax^D-kr(idims,^D)\
!  call getdivb(w,ixI^L,ixR^L,divb)
!  select case (idims)
!  {case (^D)
!     divb(ixOmax^D^D%ixO^S) = divb(ixOmax^D-1^D%ixO^S)\}
!  end select
!end if
!do idirs=1,ndir
!   Sx(ixO^S,b0_+idirs)=Sx(ixO^S,b0_+idirs)-w(ixO^S,v0_+idirs)*divb(ixO^S)
!end do

call eigenvectors(idims,ixI^L,ixO^L,w,leftm,righm,lamda)
if(iside==1) then
  ixR^L=ixO^L+kr(idims,^D);
  dw(ixO^S,1:nwprim)=(w(ixR^S,1:nwprim)-w(ixO^S,1:nwprim))*invdx(idims)
else
  ixR^L=ixO^L-kr(idims,^D);
  dw(ixO^S,1:nwprim)=(w(ixO^S,1:nwprim)-w(ixR^S,1:nwprim))*invdx(idims)
end if
do iwave=1,nwwave
  RHS(ixO^S,iwave)=-lamda(ixO^S,iwave)*sum(leftm(ixO^S,iwave,1:nwprim)*&
    dw(ixO^S,1:nwprim),^ND+1)+sum(leftm(ixO^S,iwave,1:nwprim)*Sx(ixO^S,1:nwprim),^ND+1)
  where(lamda(ixO^S,iwave)>=0.d0 .and. iside==1 .or. lamda(ixO^S,iwave)<=0.d0 .and. iside==2)
    RHS(ixO^S,iwave)=0.d0
  end where
end do
do iprim=1,nwprim
  Sx(ixO^S,iprim)=sum(righm(ixO^S,iprim,1:nwwave)*RHS(ixO^S,1:nwwave),^ND+1)
end do
w(ixO^S,1:nwprim)=w(ixO^S,1:nwprim)+qdt*Sx(ixO^S,1:nwprim)
patchw=.false.
call conserve(ixI^L,ixI^L,w,x,patchw)

end subroutine characteristic_project
!============================================================================
