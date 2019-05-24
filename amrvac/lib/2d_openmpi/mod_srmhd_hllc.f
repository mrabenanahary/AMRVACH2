module mod_srmhd_hllc
  use mod_srmhd_parameters
  use mod_srmhd_phys 
  implicit none
  private

  public :: srmhd_hllc_init

contains

  subroutine srmhd_hllc_init()
    use mod_physics_hllc

    phys_diffuse_hllcd => srmhd_diffuse_hllcd
    phys_get_lCD       => srmhd_get_lCD
    phys_get_wCD       => srmhd_get_wCD

  end subroutine srmhd_hllc_init

  subroutine srmhd_diffuse_hllcd(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim,wLC,wRC,fLC,fRC,patchf)
  ! made by Z. MELIANI 14/02/2018
  ! when method is hllcd or hllcd1 then: 
  ! this subroutine is to enforce regions where we AVOID HLLC
  ! and use TVDLF instead: this is achieved by setting patchf to 4 in
  ! certain regions. An additional input parameter is nxdiffusehllc
  ! which sets the size of the fallback region.
    use mod_global_parameters
    
    integer, intent(in)                                     :: ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        intent(in)     :: wRC,wLC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
       intent(in)  :: fLC, fRC
    integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
        intent(inout) :: patchf
    
    integer :: ixOO1,ixOO2,TxOOmin1,TxOOmin2,TxOOmax1,TxOOmax2

    
    ! In a user-controlled region around any point with flux sign change between
    ! left and right, ensure fallback to TVDLF
    do ixOO1= ixOmin1,ixOmax1
    do ixOO2= ixOmin2,ixOmax2
      
      TxOOmin1= max(ixOO1 - nxdiffusehllc*kr(idim,1), ixOmin1);
      TxOOmax1= min(ixOO1 + nxdiffusehllc*kr(idim,1), ixOmax1);
      
      
      TxOOmin2= max(ixOO2 - nxdiffusehllc*kr(idim,2), ixOmin2);
      TxOOmax2= min(ixOO2 + nxdiffusehllc*kr(idim,2), ixOmax2);
      
      if(abs(patchf(ixOO1,ixOO2)) == 1 .or. abs(patchf(ixOO1,ixOO2)) == 4)Then
         if(any(fRC(ixOO1,ixOO2,1:nwflux)*fLC(ixOO1,ixOO2,&
            1:nwflux)<-smalldouble))Then
           where(Abs(patchf(TxOOmin1:TxOOmax1,TxOOmin2:TxOOmax2))==1)
             patchf(TxOOmin1:TxOOmax1,TxOOmin2:TxOOmax2) = 4
           endwhere
         endif
      endif
    enddo
    enddo
  
  end subroutine srmhd_diffuse_hllcd

  subroutine srmhd_get_lCD(wLC,wRC,fLC,fRC,cmin,cmax,idim,ixImin1,ixImin2,&
     ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2, whll,Fhll,lambdaCD,&
     patchf)
  ! made by Z. MELIANI 14/02/2018
  
  ! Calculate lambda at CD and set the patchf to know the orientation
  ! of the riemann fan and decide on the flux choice
  ! We also compute here the HLL flux and w value, for fallback strategy
  
    use mod_global_parameters
    implicit none    

    integer, intent(in)                                     :: ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        intent(in)     :: wLC,wRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
        intent(in) :: fLC,fRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
        intent(in)          :: cmax,cmin
    integer         , dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
        intent(inout)       :: patchf
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
        intent(out):: Fhll,whll
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
        intent(out)         :: lambdaCD
    

    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       mom(1):mom(ndir))     :: vCD
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2)     :: Aco,Bco,&
       Cco,Delta, BpdotfBp, Bp2,fBp2
    logical         , dimension(ixImin1:ixImax1,&
       ixImin2:ixImax2)     :: Cond_patchf, Cond_Bidimhll
    double precision                       :: Epsilon
    integer                                :: iw

    ! on entry, patch is preset to contain values from -2,1,2,4
    !      -2: take left flux, no computation here
    !      +2: take right flux, no computation here
    !      +4: take TVDLF flux, no computation here
    !       1: compute the characteristic speed for the CD
    
    Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=abs(patchf(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))==1
    
    do iw=1,nwflux
      where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
      !============= compute HLL flux ==============!
      Fhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)= (cmax(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         iw)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*fRC(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,iw) + cmin(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*cmax(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*(wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         iw)-wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)))/(cmax(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
      !======== compute intermediate HLL state =======!
      whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = (cmax(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         iw)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*wLC(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,iw)+fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         iw)-fRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw))/(cmax(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
      endwhere
    enddo


!Calculate the Characteristic speed at the contact
! part specific to SRMHD and HLLC
! Eq. 41, Mignone & Bodo, MNRAS 2006
!  write this eq as:  Aco lambda^2 - Bco lambda + Cco = 0
where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
  ! only store the HD part first, sufficient for case where normal B vanishes
  Aco(ixOmin1:ixOmax1,ixOmin2:ixOmax2)    = Fhll(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,e_)+Fhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,d_)
  Bco(ixOmin1:ixOmax1,ixOmin2:ixOmax2)    = Fhll(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,mom(idim))+whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     e_)+whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,d_)
  Cco(ixOmin1:ixOmax1,ixOmin2:ixOmax2)    = whll(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,mom(idim))
endwhere

! condition on the normal magnetic field
Cond_Bidimhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (dabs(whll(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,mag(idim)))<=smalldouble)

! Case With Normal Magnetic field
if(any(.not.Cond_Bidimhll(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2).and.Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))then
  where(.not.Cond_Bidimhll(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2).and.Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
   !---- Initialisation ----!
    BpdotfBp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = zero
    Bp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = zero
    fBp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = zero
  endwhere

  !The calculation of the Transverse components part
  do iw=mag(1),mag(ndir)
    if(iw /= mag(idim))then
      where(.not.Cond_Bidimhll(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2) .and. Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
        BpdotfBp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = BpdotfBp(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2) + whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           iw)*Fhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)
        Bp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = Bp2(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2) + whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)**2.0d0
        fBp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = fBp2(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2) + Fhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)**2.0d0
      endwhere
     endif
    enddo

    where(.not.Cond_Bidimhll(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) .and. Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
     Aco(ixOmin1:ixOmax1,ixOmin2:ixOmax2)    = Aco(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) - BpdotfBp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
     Bco(ixOmin1:ixOmax1,ixOmin2:ixOmax2)    = Bco(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)- Bp2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) - fBp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
     Cco(ixOmin1:ixOmax1,ixOmin2:ixOmax2)    = Cco(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) - BpdotfBp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    endwhere
   endif


   where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    Delta(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = Bco(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)**2.0d0- 4.0d0*Aco(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * Cco(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    where(Aco(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/=zero .and. &
       Delta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)>=zero)
     !Calculate the Characteristic speed at the contact
     ! only the minus sign is between [-1,1]
     lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (Bco(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) - dsqrt(Delta(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)))/(2.0d0*Aco(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    elsewhere(Aco(ixOmin1:ixOmax1,ixOmin2:ixOmax2)==zero .and.  &
       Bco(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/=zero)
     lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = Cco(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)/Bco(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    elsewhere(Delta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<zero)
     lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = zero
     ! we will fall back to HLL flux case in this degeneracy
     patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) =  3
    endwhere
   endwhere

   where(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)==3)
    Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=.false.
   end where

   where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    ! double check whether obtained speed is in between min and max speeds given
    ! and identify in which part of the Riemann fan the time-axis is
    where(cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<zero.and.&
       lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)>zero.and.&
       lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<cmax(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))
     patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = -1
    elsewhere(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)>zero.and.&
       lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<zero.and.&
       lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)>cmin(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))
     patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) =  1
    elsewhere(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)>=cmax(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2).or.lambdaCD(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) <= cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
     lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = zero
     ! we will fall back to HLL flux case in this degeneracy
     patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) =  3
    endwhere
   endwhere

   where(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)== 3)
    Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=.false.
   end where


   ! handle the specific case where the time axis is exactly on the CD
   if(any(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)==&
      zero.and.Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))then
    !determine which sector (forward or backward) of the Riemann fan is smallest
    ! and select left or right flux accordingly
    where(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)==&
       zero.and.Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
     where(-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)>=cmax(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2))
      patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) =  1
     elsewhere
      patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = -1
     endwhere
    endwhere
   endif



   ! eigenvalue lambda for contact is near zero: decrease noise by this trick
   if(flathllc)then
    Epsilon=1.0d-6
    where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2).and. &
       dabs(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2))/max(cmax(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2),Epsilon)< Epsilon  .and. &
       dabs(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2))/max(dabs(cmin(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2)),Epsilon)< Epsilon)
     lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2) =  zero
    end where
   end if
    

   if(any(dabs(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2))>one .and. &
      Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))then
    call mpistop("problems with lambdaCD>1")
   endif

   ! next should never happen
   if(any(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)==0))then
    call mpistop("patchf=0")
   endif


  end subroutine srmhd_get_lCD

  subroutine srmhd_get_wCD(x,wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,cmin,&
     cmax,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,&
     f)
  ! made by Z. MELIANI 14/02/2018
  ! compute the intermediate state U*
  ! only needed where patchf=-1/1
  
  ! reference Li S., JCP, 203, 2005, 344-357
  ! reference T. Miyoski, Kusano JCP, 2008, 2005.
    use mod_global_parameters
    
    integer, intent(in)                                      :: ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        intent(in)      :: wRC,wLC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
        intent(in)    :: x
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
        intent(in)  :: whll, Fhll
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
        intent(in)           :: lambdaCD
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
        intent(in)           :: cmax,cmin
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
        intent(in)  :: fRC,fLC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
       intent(out)  :: f
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)         :: wCD,wSub
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux)     :: fSub
    double precision, dimension(ixImin1:ixImax1,&
       ixImin2:ixImax2)              :: vSub,cspeed,pCD,VdotBCD
    integer         , dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
       intent(inout):: patchf

    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       mom(1):mom(ndir))      :: vCD
    double precision, dimension(ixImin1:ixImax1,&
       ixImin2:ixImax2)             :: Ratio_CD
    double precision, dimension(ixImin1:ixImax1,&
       ixImin2:ixImax2)             :: vRC,vLC
    logical         , dimension(ixImin1:ixImax1,&
       ixImin2:ixImax2)             :: Cond_Bidimhll
    integer                                        :: n, iw, idir,ix1,ix2
    double precision, dimension(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)             :: VdotB, B2

    !-------------- auxiliary Speed and array-------------!
    call srmhd_get_v_idim(wRC,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,idim,vRC)
    call srmhd_get_v_idim(wLC,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,idim,vLC)


    where(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) == 1)
     cspeed(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cmax(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)
     vSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2)   = vRC(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)
    elsewhere(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) == -1)
     cspeed(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cmin(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)
     vSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2)   = vLC(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)
    endwhere

    do iw=1,nwflux
     if(iw /= mag(idim))then
      where(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) == 1)
       wSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) =  wRC(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,iw)
       fSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) =  fRC(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,iw)
      elsewhere(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) == -1)
       wSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) =  wLC(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,iw)
       fSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) =  fLC(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,iw)
      endwhere
     endif
    enddo

    wSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idim)) = whll(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mag(idim))
    Cond_Bidimhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
       (dabs(wSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idim)))<=smalldouble)

    where(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)
     Ratio_CD(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (cspeed(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)-vSub(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2))/(cspeed(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)-lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
     wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,d_)   = wSub(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,d_)*Ratio_CD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end where
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
      if(abs(patchf(ix1,ix2))==1) then
       do n=1,srmhd_n_tracer
          iw = tracer(n)
          wCD(ix1,ix2,iw) = wSub(ix1,ix2,iw)*Ratio_CD(ix1,ix2)
        end do
      end if
    end do
    end do







     ! in case we have somewhere a normal component, need to distinguish
     bxnozero : if(any(.not.Cond_Bidimhll(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)))then

     !==== Magnetic field ====!
     do iw =mag(1),mag(ndir)
      if(iw /= mag(idim))then
       ! Transverse components
       where(.not.Cond_Bidimhll(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)  .and. abs(patchf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2))==1)
         ! case from eq 37
         wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = whll(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw)
       elsewhere(Cond_Bidimhll(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) .and. abs(patchf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2))==1)
         ! case from eq 53
         wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = wSub(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw) * Ratio_CD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       endwhere
      else  !  Normal component
       wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = wSub(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,iw)
      endif
     enddo

     !====== velocities ========!
     do idir = 1,ndir
      if(idir /=idim)then
       ! Transverse components
       where(.not.Cond_Bidimhll(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)  .and. abs(patchf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2))==1)
        ! case from eq 38
        vCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=(wCD(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mag(idir))*lambdaCD(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-Fhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(idir)))/wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idim))
       elsewhere(Cond_Bidimhll(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)  .and. abs(patchf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2))==1)
        ! unused case
        vCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=zero
       endwhere
      else ! Normal component
       where(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)
        vCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir))=lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       endwhere
      endif
     enddo
     ! enforce fallback strategy for case where characteristic velocity 
     !  at contact
     ! is unphysical
     where(.not. Cond_Bidimhll(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)  .and. sum(vCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        mom(:))**2.0d0,dim=ndim+1) > one .and. abs(patchf(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2))==1)
      patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 3
     endwhere

     where(.not.Cond_Bidimhll(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2).and. abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)
      wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         lfac_) = one/dsqrt(one-sum(vCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(:))**2.0d0,dim=ndim+1))
      VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = sum(vCD(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(:))*wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:)),&
         dim=ndim+1)
      !--- total Pressure from eq 40 ---!
      pCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  = -lambdaCD(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*(Fhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_)+Fhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,d_)-VdotB(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mag(idim)))+Fhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idim))+(wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mag(idim))/wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_))**2.0d0

     elsewhere(Cond_Bidimhll(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) .and. abs(patchf(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2))==1)
      !------ total Pressure from 48 ------!
      pCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  = -(Fhll(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,e_)+Fhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         d_))*lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2) + Fhll(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim))
     endwhere

     ! enforce fallback strategy for case where total pressure at contact
     ! is unphysical (should never happen?)
     where(pCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2) <= zero.and. &
        abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)
      patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 3
     endwhere
     !------- Momentum ------!
     do iw =mom(1),mom(ndir)
      if(iw /= mom(idim))then
       where(.not.Cond_Bidimhll(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)   .and.  abs(patchf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2))==1)
        ! eq. 44 45
        wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = (cspeed(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)*wSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           iw)-fSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)-wCD(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mag(idim))*(wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           iw-mom(1)+mag(1))/wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           lfac_)**2.0d0+VdotB(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2) * vCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           iw))) /(cspeed(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
       elsewhere(Cond_Bidimhll(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) .and. abs(patchf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2))==1)
        ! eq. 51
        wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = wSub(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,iw) * Ratio_CD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       endwhere
      endif
     enddo

     where(.not.Cond_Bidimhll(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)   .and. abs(patchf(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2))==1)
      !---- Tau Right combine 46 43----!
      wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = (cspeed(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2) * wSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_) + vSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2) * wSub(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,d_) - wSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idim))+lambdaCD(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*pCD(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)-VdotB(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mag(idim)))/(cspeed(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)-lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
      !--- Sidim Right from eq 33 ---!
      wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idim)) = (wCD(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,e_)+wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         d_)+pCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2))*lambdaCD(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)-VdotB(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idim))
     elsewhere(Cond_Bidimhll(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)  .and. abs(patchf(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2))==1)
      !---- Tau Right combine 50 52----!
      wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = (cspeed(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2) * wSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_) + vSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2) * wSub(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,d_) - wSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idim))+lambdaCD(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*pCD(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))/(cspeed(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)-lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
      !--- sidim Right from 49 ---!
      wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idim)) = (wCD(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,e_)+wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         d_)+pCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2))*lambdaCD(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)
      !-- Not real lfac, but fill it anyway --!
      wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_) = one
      !-------------------!
     endwhere
else ! bxnozero
  ! in case we have everywhere NO normal B component, no need to distinguish
  !---- Magnetic field ----!
  do iw =mag(1),mag(ndir)
    if(iw /= mag(idim))then
      !Transverse components
      where(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)
        wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = wSub(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,iw) * Ratio_CD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      endwhere
    else   !Normal components
      where(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)
        wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = wSub(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,iw)
      endwhere
    endif
  enddo

  !---- momentum ----!
  do iw =mom(1),mom(ndir)
    if(iw /= mom(idim))Then
      where(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)
        wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = wSub(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,iw) * Ratio_CD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      endwhere
    endif
  enddo

  where(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)
   !------- Pressure -------!
   pCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)      = Fhll(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,mom(idim))- (Fhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      e_)+Fhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,d_))* lambdaCD(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)
   !---- Tau Right ---!
   wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = (cspeed(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2) * wSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      e_) + vSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2) * wSub(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,d_) - wSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      mom(idim))+lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*pCD(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2))/(cspeed(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)-lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
   !---- Sx_ Right ---!
   wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idim)) =(wCD(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,e_)+wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      d_)+pCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2))*lambdaCD(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)
   !-- Not real lfac, but fill it anyway --!
   wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_) = one
  endwhere

 ! In case of B_idim = 0, we need only the normal velocity
  do iw = mom(1),mom(ndir)
    if(iw /= mom(idim))Then
      where(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)
        vCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = zero
      endwhere
    else
      where(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)
        vCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = lambdaCD(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)
      endwhere
    endif
  enddo

endif bxnozero

do iw=1,nwflux
 if(iw /= mag(idim).and. iw/=psi_)then
   where(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)
       ! f_i=fSub+lambda (wCD-wSub)
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=fSub(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,iw)+cspeed(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)*(wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          iw)-wsub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw))
   endwhere
 else
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=zero
 end if
end do

  end subroutine srmhd_get_wCD

end module mod_srmhd_hllc
