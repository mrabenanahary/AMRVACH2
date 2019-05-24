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

  subroutine srmhd_diffuse_hllcd(ixImin1,ixImax1,ixOmin1,ixOmax1,idim,wLC,wRC,&
     fLC,fRC,patchf)
  ! made by Z. MELIANI 14/02/2018
  ! when method is hllcd or hllcd1 then: 
  ! this subroutine is to enforce regions where we AVOID HLLC
  ! and use TVDLF instead: this is achieved by setting patchf to 4 in
  ! certain regions. An additional input parameter is nxdiffusehllc
  ! which sets the size of the fallback region.
    use mod_global_parameters
    
    integer, intent(in)                                     :: ixImin1,ixImax1,&
       ixOmin1,ixOmax1,idim
    double precision, dimension(ixImin1:ixImax1,1:nw), intent(in)     :: wRC,&
       wLC
    double precision, dimension(ixImin1:ixImax1,1:nwflux),intent(in)  :: fLC,&
        fRC
    integer, dimension(ixImin1:ixImax1), intent(inout) :: patchf
    
    integer :: ixOO1,TxOOmin1,TxOOmax1

    
    ! In a user-controlled region around any point with flux sign change between
    ! left and right, ensure fallback to TVDLF
    do ixOO1= ixOmin1,ixOmax1
      
      TxOOmin1= max(ixOO1 - nxdiffusehllc*kr(idim,1), ixOmin1);
      TxOOmax1= min(ixOO1 + nxdiffusehllc*kr(idim,1), ixOmax1);
      
      if(abs(patchf(ixOO1)) == 1 .or. abs(patchf(ixOO1)) == 4)Then
         if(any(fRC(ixOO1,1:nwflux)*fLC(ixOO1,1:nwflux)<-smalldouble))Then
           where(Abs(patchf(TxOOmin1:TxOOmax1))==1)
             patchf(TxOOmin1:TxOOmax1) = 4
           endwhere
         endif
      endif
    enddo
  
  end subroutine srmhd_diffuse_hllcd

  subroutine srmhd_get_lCD(wLC,wRC,fLC,fRC,cmin,cmax,idim,ixImin1,ixImax1,&
     ixOmin1,ixOmax1, whll,Fhll,lambdaCD,patchf)
  ! made by Z. MELIANI 14/02/2018
  
  ! Calculate lambda at CD and set the patchf to know the orientation
  ! of the riemann fan and decide on the flux choice
  ! We also compute here the HLL flux and w value, for fallback strategy
  
    use mod_global_parameters
    implicit none    

    integer, intent(in)                                     :: ixImin1,ixImax1,&
       ixOmin1,ixOmax1,idim
    double precision, dimension(ixImin1:ixImax1,1:nw), intent(in)     :: wLC,&
       wRC
    double precision, dimension(ixImin1:ixImax1,1:nwflux), intent(in) :: fLC,&
       fRC
    double precision, dimension(ixImin1:ixImax1), intent(in)          :: cmax,&
       cmin
    integer         , dimension(ixImin1:ixImax1),&
        intent(inout)       :: patchf
    double precision, dimension(ixImin1:ixImax1,1:nwflux), intent(out):: Fhll,&
       whll
    double precision, dimension(ixImin1:ixImax1),&
        intent(out)         :: lambdaCD
    

    double precision, dimension(ixGlo1:ixGhi1,mom(1):mom(ndir))     :: vCD
    double precision, dimension(ixGlo1:ixGhi1)     :: Aco,Bco,Cco,Delta,&
        BpdotfBp, Bp2,fBp2
    logical         , dimension(ixImin1:ixImax1)     :: Cond_patchf,&
        Cond_Bidimhll
    double precision                       :: Epsilon
    integer                                :: iw

    ! on entry, patch is preset to contain values from -2,1,2,4
    !      -2: take left flux, no computation here
    !      +2: take right flux, no computation here
    !      +4: take TVDLF flux, no computation here
    !       1: compute the characteristic speed for the CD
    
    Cond_patchf(ixOmin1:ixOmax1)=abs(patchf(ixOmin1:ixOmax1))==1
    
    do iw=1,nwflux
      where(Cond_patchf(ixOmin1:ixOmax1))
      !============= compute HLL flux ==============!
      Fhll(ixOmin1:ixOmax1,iw)= (cmax(ixOmin1:ixOmax1)*fLC(ixOmin1:ixOmax1,&
         iw)-cmin(ixOmin1:ixOmax1)*fRC(ixOmin1:ixOmax1,&
         iw) + cmin(ixOmin1:ixOmax1)*cmax(ixOmin1:ixOmax1)*(wRC(&
         ixOmin1:ixOmax1,iw)-wLC(ixOmin1:ixOmax1,&
         iw)))/(cmax(ixOmin1:ixOmax1)-cmin(ixOmin1:ixOmax1))
      !======== compute intermediate HLL state =======!
      whll(ixOmin1:ixOmax1,iw) = (cmax(ixOmin1:ixOmax1)*wRC(ixOmin1:ixOmax1,&
         iw)-cmin(ixOmin1:ixOmax1)*wLC(ixOmin1:ixOmax1,iw)+fLC(ixOmin1:ixOmax1,&
         iw)-fRC(ixOmin1:ixOmax1,iw))/(cmax(ixOmin1:ixOmax1)-&
         cmin(ixOmin1:ixOmax1))
      endwhere
    enddo


!Calculate the Characteristic speed at the contact
! part specific to SRMHD and HLLC
! Eq. 41, Mignone & Bodo, MNRAS 2006
!  write this eq as:  Aco lambda^2 - Bco lambda + Cco = 0
where(Cond_patchf(ixOmin1:ixOmax1))
  ! only store the HD part first, sufficient for case where normal B vanishes
  Aco(ixOmin1:ixOmax1)    = Fhll(ixOmin1:ixOmax1,e_)+Fhll(ixOmin1:ixOmax1,d_)
  Bco(ixOmin1:ixOmax1)    = Fhll(ixOmin1:ixOmax1,&
     mom(idim))+whll(ixOmin1:ixOmax1,e_)+whll(ixOmin1:ixOmax1,d_)
  Cco(ixOmin1:ixOmax1)    = whll(ixOmin1:ixOmax1,mom(idim))
endwhere

! condition on the normal magnetic field
Cond_Bidimhll(ixOmin1:ixOmax1) = (dabs(whll(ixOmin1:ixOmax1,&
   mag(idim)))<=smalldouble)

! Case With Normal Magnetic field
if(any(.not.Cond_Bidimhll(ixOmin1:ixOmax1).and.&
   Cond_patchf(ixOmin1:ixOmax1)))then
  where(.not.Cond_Bidimhll(ixOmin1:ixOmax1).and.Cond_patchf(ixOmin1:ixOmax1))
   !---- Initialisation ----!
    BpdotfBp(ixOmin1:ixOmax1) = zero
    Bp2(ixOmin1:ixOmax1) = zero
    fBp2(ixOmin1:ixOmax1) = zero
  endwhere

  !The calculation of the Transverse components part
  do iw=mag(1),mag(ndir)
    if(iw /= mag(idim))then
      where(.not.Cond_Bidimhll(ixOmin1:ixOmax1) .and. &
         Cond_patchf(ixOmin1:ixOmax1))
        BpdotfBp(ixOmin1:ixOmax1) = BpdotfBp(ixOmin1:ixOmax1) + &
           whll(ixOmin1:ixOmax1,iw)*Fhll(ixOmin1:ixOmax1,iw)
        Bp2(ixOmin1:ixOmax1) = Bp2(ixOmin1:ixOmax1) + whll(ixOmin1:ixOmax1,&
           iw)**2.0d0
        fBp2(ixOmin1:ixOmax1) = fBp2(ixOmin1:ixOmax1) + Fhll(ixOmin1:ixOmax1,&
           iw)**2.0d0
      endwhere
     endif
    enddo

    where(.not.Cond_Bidimhll(ixOmin1:ixOmax1) .and. &
       Cond_patchf(ixOmin1:ixOmax1))
     Aco(ixOmin1:ixOmax1)    = Aco(ixOmin1:ixOmax1) - &
        BpdotfBp(ixOmin1:ixOmax1)
     Bco(ixOmin1:ixOmax1)    = Bco(ixOmin1:ixOmax1)- Bp2(ixOmin1:ixOmax1) - &
        fBp2(ixOmin1:ixOmax1)
     Cco(ixOmin1:ixOmax1)    = Cco(ixOmin1:ixOmax1) - &
        BpdotfBp(ixOmin1:ixOmax1)
    endwhere
   endif


   where(Cond_patchf(ixOmin1:ixOmax1))
    Delta(ixOmin1:ixOmax1) = Bco(ixOmin1:ixOmax1)**2.0d0- &
       4.0d0*Aco(ixOmin1:ixOmax1) * Cco(ixOmin1:ixOmax1)
    where(Aco(ixOmin1:ixOmax1)/=zero .and. Delta(ixOmin1:ixOmax1)>=zero)
     !Calculate the Characteristic speed at the contact
     ! only the minus sign is between [-1,1]
     lambdaCD(ixOmin1:ixOmax1) = (Bco(ixOmin1:ixOmax1) - &
        dsqrt(Delta(ixOmin1:ixOmax1)))/(2.0d0*Aco(ixOmin1:ixOmax1))
    elsewhere(Aco(ixOmin1:ixOmax1)==zero .and.  Bco(ixOmin1:ixOmax1)/=zero)
     lambdaCD(ixOmin1:ixOmax1) = Cco(ixOmin1:ixOmax1)/Bco(ixOmin1:ixOmax1)
    elsewhere(Delta(ixOmin1:ixOmax1)<zero)
     lambdaCD(ixOmin1:ixOmax1) = zero
     ! we will fall back to HLL flux case in this degeneracy
     patchf(ixOmin1:ixOmax1) =  3
    endwhere
   endwhere

   where(patchf(ixOmin1:ixOmax1)==3)
    Cond_patchf(ixOmin1:ixOmax1)=.false.
   end where

   where(Cond_patchf(ixOmin1:ixOmax1))
    ! double check whether obtained speed is in between min and max speeds given
    ! and identify in which part of the Riemann fan the time-axis is
    where(cmin(ixOmin1:ixOmax1)<zero.and.lambdaCD(ixOmin1:ixOmax1)>zero.and.&
       lambdaCD(ixOmin1:ixOmax1)<cmax(ixOmin1:ixOmax1))
     patchf(ixOmin1:ixOmax1) = -1
    elsewhere(cmax(ixOmin1:ixOmax1)>zero.and.lambdaCD(ixOmin1:ixOmax1)<zero.and.&
       lambdaCD(ixOmin1:ixOmax1)>cmin(ixOmin1:ixOmax1))
     patchf(ixOmin1:ixOmax1) =  1
    elsewhere(lambdaCD(ixOmin1:ixOmax1)>=cmax(ixOmin1:ixOmax1).or.&
       lambdaCD(ixOmin1:ixOmax1) <= cmin(ixOmin1:ixOmax1))
     lambdaCD(ixOmin1:ixOmax1) = zero
     ! we will fall back to HLL flux case in this degeneracy
     patchf(ixOmin1:ixOmax1) =  3
    endwhere
   endwhere

   where(patchf(ixOmin1:ixOmax1)== 3)
    Cond_patchf(ixOmin1:ixOmax1)=.false.
   end where


   ! handle the specific case where the time axis is exactly on the CD
   if(any(lambdaCD(ixOmin1:ixOmax1)==zero.and.&
      Cond_patchf(ixOmin1:ixOmax1)))then
    !determine which sector (forward or backward) of the Riemann fan is smallest
    ! and select left or right flux accordingly
    where(lambdaCD(ixOmin1:ixOmax1)==zero.and.Cond_patchf(ixOmin1:ixOmax1))
     where(-cmin(ixOmin1:ixOmax1)>=cmax(ixOmin1:ixOmax1))
      patchf(ixOmin1:ixOmax1) =  1
     elsewhere
      patchf(ixOmin1:ixOmax1) = -1
     endwhere
    endwhere
   endif



   ! eigenvalue lambda for contact is near zero: decrease noise by this trick
   if(flathllc)then
    Epsilon=1.0d-6
    where(Cond_patchf(ixOmin1:ixOmax1).and. &
       dabs(lambdaCD(ixOmin1:ixOmax1))/max(cmax(ixOmin1:ixOmax1),&
       Epsilon)< Epsilon  .and. dabs(lambdaCD(ixOmin1:ixOmax1))/max(dabs(cmin(&
       ixOmin1:ixOmax1)),Epsilon)< Epsilon)
     lambdaCD(ixOmin1:ixOmax1) =  zero
    end where
   end if
    

   if(any(dabs(lambdaCD(ixOmin1:ixOmax1))>one .and. &
      Cond_patchf(ixOmin1:ixOmax1)))then
    call mpistop("problems with lambdaCD>1")
   endif

   ! next should never happen
   if(any(patchf(ixOmin1:ixOmax1)==0))then
    call mpistop("patchf=0")
   endif


  end subroutine srmhd_get_lCD

  subroutine srmhd_get_wCD(x,wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,cmin,&
     cmax,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,f)
  ! made by Z. MELIANI 14/02/2018
  ! compute the intermediate state U*
  ! only needed where patchf=-1/1
  
  ! reference Li S., JCP, 203, 2005, 344-357
  ! reference T. Miyoski, Kusano JCP, 2008, 2005.
    use mod_global_parameters
    
    integer, intent(in)                                      :: ixImin1,&
       ixImax1,ixOmin1,ixOmax1,idim
    double precision, dimension(ixImin1:ixImax1,1:nw), intent(in)      :: wRC,&
       wLC
    double precision, dimension(ixImin1:ixImax1,1:ndim), intent(in)    :: x
    double precision, dimension(ixImin1:ixImax1,1:nwflux), intent(in)  :: whll,&
        Fhll
    double precision, dimension(ixImin1:ixImax1),&
        intent(in)           :: lambdaCD
    double precision, dimension(ixImin1:ixImax1), intent(in)           :: cmax,&
       cmin
    double precision, dimension(ixImin1:ixImax1,1:nwflux), intent(in)  :: fRC,&
       fLC
    double precision, dimension(ixImin1:ixImax1,1:nwflux),intent(out)  :: f
    double precision, dimension(ixImin1:ixImax1,1:nw)         :: wCD,wSub
    double precision, dimension(ixImin1:ixImax1,1:nwflux)     :: fSub
    double precision, dimension(ixImin1:ixImax1)              :: vSub,cspeed,&
       pCD,VdotBCD
    integer         , dimension(ixImin1:ixImax1),intent(inout):: patchf

    double precision, dimension(ixGlo1:ixGhi1,mom(1):mom(ndir))      :: vCD
    double precision, dimension(ixImin1:ixImax1)             :: Ratio_CD
    double precision, dimension(ixImin1:ixImax1)             :: vRC,vLC
    logical         , dimension(ixImin1:ixImax1)             :: Cond_Bidimhll
    integer                                        :: n, iw, idir,ix1
    double precision, dimension(ixOmin1:ixOmax1)             :: VdotB, B2

    !-------------- auxiliary Speed and array-------------!
    call srmhd_get_v_idim(wRC,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,vRC)
    call srmhd_get_v_idim(wLC,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,vLC)


    where(patchf(ixOmin1:ixOmax1) == 1)
     cspeed(ixOmin1:ixOmax1) = cmax(ixOmin1:ixOmax1)
     vSub(ixOmin1:ixOmax1)   = vRC(ixOmin1:ixOmax1)
    elsewhere(patchf(ixOmin1:ixOmax1) == -1)
     cspeed(ixOmin1:ixOmax1) = cmin(ixOmin1:ixOmax1)
     vSub(ixOmin1:ixOmax1)   = vLC(ixOmin1:ixOmax1)
    endwhere

    do iw=1,nwflux
     if(iw /= mag(idim))then
      where(patchf(ixOmin1:ixOmax1) == 1)
       wSub(ixOmin1:ixOmax1,iw) =  wRC(ixOmin1:ixOmax1,iw)
       fSub(ixOmin1:ixOmax1,iw) =  fRC(ixOmin1:ixOmax1,iw)
      elsewhere(patchf(ixOmin1:ixOmax1) == -1)
       wSub(ixOmin1:ixOmax1,iw) =  wLC(ixOmin1:ixOmax1,iw)
       fSub(ixOmin1:ixOmax1,iw) =  fLC(ixOmin1:ixOmax1,iw)
      endwhere
     endif
    enddo

    wSub(ixOmin1:ixOmax1,mag(idim)) = whll(ixOmin1:ixOmax1,mag(idim))
    Cond_Bidimhll(ixOmin1:ixOmax1) = (dabs(wSub(ixOmin1:ixOmax1,&
       mag(idim)))<=smalldouble)

    where(abs(patchf(ixOmin1:ixOmax1))==1)
     Ratio_CD(ixOmin1:ixOmax1) = (cspeed(ixOmin1:ixOmax1)-&
        vSub(ixOmin1:ixOmax1))/(cspeed(ixOmin1:ixOmax1)-&
        lambdaCD(ixOmin1:ixOmax1))
     wCD(ixOmin1:ixOmax1,d_)   = wSub(ixOmin1:ixOmax1,&
        d_)*Ratio_CD(ixOmin1:ixOmax1)
    end where
    do ix1=ixOmin1,ixOmax1
      if(abs(patchf(ix1))==1) then
       do n=1,srmhd_n_tracer
          iw = tracer(n)
          wCD(ix1,iw) = wSub(ix1,iw)*Ratio_CD(ix1)
        end do
      end if
    end do







     ! in case we have somewhere a normal component, need to distinguish
     bxnozero : if(any(.not.Cond_Bidimhll(ixOmin1:ixOmax1)))then

     !==== Magnetic field ====!
     do iw =mag(1),mag(ndir)
      if(iw /= mag(idim))then
       ! Transverse components
       where(.not.Cond_Bidimhll(ixOmin1:ixOmax1)  .and. &
          abs(patchf(ixOmin1:ixOmax1))==1)
         ! case from eq 37
         wCD(ixOmin1:ixOmax1,iw) = whll(ixOmin1:ixOmax1,iw)
       elsewhere(Cond_Bidimhll(ixOmin1:ixOmax1) .and. &
          abs(patchf(ixOmin1:ixOmax1))==1)
         ! case from eq 53
         wCD(ixOmin1:ixOmax1,iw) = wSub(ixOmin1:ixOmax1,&
            iw) * Ratio_CD(ixOmin1:ixOmax1)
       endwhere
      else  !  Normal component
       wCD(ixOmin1:ixOmax1,iw) = wSub(ixOmin1:ixOmax1,iw)
      endif
     enddo

     !====== velocities ========!
     do idir = 1,ndir
      if(idir /=idim)then
       ! Transverse components
       where(.not.Cond_Bidimhll(ixOmin1:ixOmax1)  .and. &
          abs(patchf(ixOmin1:ixOmax1))==1)
        ! case from eq 38
        vCD(ixOmin1:ixOmax1,mom(idir))=(wCD(ixOmin1:ixOmax1,&
           mag(idir))*lambdaCD(ixOmin1:ixOmax1)-Fhll(ixOmin1:ixOmax1,&
           mag(idir)))/wCD(ixOmin1:ixOmax1,mag(idim))
       elsewhere(Cond_Bidimhll(ixOmin1:ixOmax1)  .and. &
          abs(patchf(ixOmin1:ixOmax1))==1)
        ! unused case
        vCD(ixOmin1:ixOmax1,mom(idir))=zero
       endwhere
      else ! Normal component
       where(abs(patchf(ixOmin1:ixOmax1))==1)
        vCD(ixOmin1:ixOmax1,mom(idir))=lambdaCD(ixOmin1:ixOmax1)
       endwhere
      endif
     enddo
     ! enforce fallback strategy for case where characteristic velocity 
     !  at contact
     ! is unphysical
     where(.not. Cond_Bidimhll(ixOmin1:ixOmax1)  .and. sum(vCD(ixOmin1:ixOmax1,&
        mom(:))**2.0d0,dim=ndim+1) > one .and. &
        abs(patchf(ixOmin1:ixOmax1))==1)
      patchf(ixOmin1:ixOmax1) = 3
     endwhere

     where(.not.Cond_Bidimhll(ixOmin1:ixOmax1).and. &
        abs(patchf(ixOmin1:ixOmax1))==1)
      wCD(ixOmin1:ixOmax1,lfac_) = one/dsqrt(one-sum(vCD(ixOmin1:ixOmax1,&
         mom(:))**2.0d0,dim=ndim+1))
      VdotB(ixOmin1:ixOmax1) = sum(vCD(ixOmin1:ixOmax1,&
         mom(:))*wCD(ixOmin1:ixOmax1,mag(:)),dim=ndim+1)
      !--- total Pressure from eq 40 ---!
      pCD(ixOmin1:ixOmax1)  = -lambdaCD(ixOmin1:ixOmax1)*(Fhll(ixOmin1:ixOmax1,&
         e_)+Fhll(ixOmin1:ixOmax1,d_)-VdotB(ixOmin1:ixOmax1)*wCD(&
         ixOmin1:ixOmax1,mag(idim)))+Fhll(ixOmin1:ixOmax1,&
         mom(idim))+(wCD(ixOmin1:ixOmax1,mag(idim))/wCD(ixOmin1:ixOmax1,&
         lfac_))**2.0d0

     elsewhere(Cond_Bidimhll(ixOmin1:ixOmax1) .and. &
        abs(patchf(ixOmin1:ixOmax1))==1)
      !------ total Pressure from 48 ------!
      pCD(ixOmin1:ixOmax1)  = -(Fhll(ixOmin1:ixOmax1,e_)+Fhll(ixOmin1:ixOmax1,&
         d_))*lambdaCD(ixOmin1:ixOmax1) + Fhll(ixOmin1:ixOmax1,mom(idim))
     endwhere

     ! enforce fallback strategy for case where total pressure at contact
     ! is unphysical (should never happen?)
     where(pCD(ixOmin1:ixOmax1) <= zero.and. abs(patchf(ixOmin1:ixOmax1))==1)
      patchf(ixOmin1:ixOmax1) = 3
     endwhere
     !------- Momentum ------!
     do iw =mom(1),mom(ndir)
      if(iw /= mom(idim))then
       where(.not.Cond_Bidimhll(ixOmin1:ixOmax1)   .and.  &
          abs(patchf(ixOmin1:ixOmax1))==1)
        ! eq. 44 45
        wCD(ixOmin1:ixOmax1,iw) = (cspeed(ixOmin1:ixOmax1)*wSub(&
           ixOmin1:ixOmax1,iw)-fSub(ixOmin1:ixOmax1,iw)-wCD(ixOmin1:ixOmax1,&
           mag(idim))*(wCD(ixOmin1:ixOmax1,&
           iw-mom(1)+mag(1))/wCD(ixOmin1:ixOmax1,&
           lfac_)**2.0d0+VdotB(ixOmin1:ixOmax1) * vCD(ixOmin1:ixOmax1,&
           iw))) /(cspeed(ixOmin1:ixOmax1)-lambdaCD(ixOmin1:ixOmax1))
       elsewhere(Cond_Bidimhll(ixOmin1:ixOmax1) .and. &
          abs(patchf(ixOmin1:ixOmax1))==1)
        ! eq. 51
        wCD(ixOmin1:ixOmax1,iw) = wSub(ixOmin1:ixOmax1,&
           iw) * Ratio_CD(ixOmin1:ixOmax1)
       endwhere
      endif
     enddo

     where(.not.Cond_Bidimhll(ixOmin1:ixOmax1)   .and. &
        abs(patchf(ixOmin1:ixOmax1))==1)
      !---- Tau Right combine 46 43----!
      wCD(ixOmin1:ixOmax1,e_) = (cspeed(ixOmin1:ixOmax1) * &
         wSub(ixOmin1:ixOmax1,e_) + vSub(ixOmin1:ixOmax1) * &
         wSub(ixOmin1:ixOmax1,d_) - wSub(ixOmin1:ixOmax1,&
         mom(idim))+lambdaCD(ixOmin1:ixOmax1)*pCD(ixOmin1:ixOmax1)-&
         VdotB(ixOmin1:ixOmax1)*wCD(ixOmin1:ixOmax1,&
         mag(idim)))/(cspeed(ixOmin1:ixOmax1)-lambdaCD(ixOmin1:ixOmax1))
      !--- Sidim Right from eq 33 ---!
      wCD(ixOmin1:ixOmax1,mom(idim)) = (wCD(ixOmin1:ixOmax1,&
         e_)+wCD(ixOmin1:ixOmax1,d_)+pCD(ixOmin1:ixOmax1))*lambdaCD(&
         ixOmin1:ixOmax1)-VdotB(ixOmin1:ixOmax1)*wCD(ixOmin1:ixOmax1,&
         mag(idim))
     elsewhere(Cond_Bidimhll(ixOmin1:ixOmax1)  .and. &
        abs(patchf(ixOmin1:ixOmax1))==1)
      !---- Tau Right combine 50 52----!
      wCD(ixOmin1:ixOmax1,e_) = (cspeed(ixOmin1:ixOmax1) * &
         wSub(ixOmin1:ixOmax1,e_) + vSub(ixOmin1:ixOmax1) * &
         wSub(ixOmin1:ixOmax1,d_) - wSub(ixOmin1:ixOmax1,&
         mom(idim))+lambdaCD(ixOmin1:ixOmax1)*pCD(ixOmin1:ixOmax1))/(cspeed(&
         ixOmin1:ixOmax1)-lambdaCD(ixOmin1:ixOmax1))
      !--- sidim Right from 49 ---!
      wCD(ixOmin1:ixOmax1,mom(idim)) = (wCD(ixOmin1:ixOmax1,&
         e_)+wCD(ixOmin1:ixOmax1,d_)+pCD(ixOmin1:ixOmax1))*lambdaCD(&
         ixOmin1:ixOmax1)
      !-- Not real lfac, but fill it anyway --!
      wCD(ixOmin1:ixOmax1,lfac_) = one
      !-------------------!
     endwhere
else ! bxnozero
  ! in case we have everywhere NO normal B component, no need to distinguish
  !---- Magnetic field ----!
  do iw =mag(1),mag(ndir)
    if(iw /= mag(idim))then
      !Transverse components
      where(abs(patchf(ixOmin1:ixOmax1))==1)
        wCD(ixOmin1:ixOmax1,iw) = wSub(ixOmin1:ixOmax1,&
           iw) * Ratio_CD(ixOmin1:ixOmax1)
      endwhere
    else   !Normal components
      where(abs(patchf(ixOmin1:ixOmax1))==1)
        wCD(ixOmin1:ixOmax1,iw) = wSub(ixOmin1:ixOmax1,iw)
      endwhere
    endif
  enddo

  !---- momentum ----!
  do iw =mom(1),mom(ndir)
    if(iw /= mom(idim))Then
      where(abs(patchf(ixOmin1:ixOmax1))==1)
        wCD(ixOmin1:ixOmax1,iw) = wSub(ixOmin1:ixOmax1,&
           iw) * Ratio_CD(ixOmin1:ixOmax1)
      endwhere
    endif
  enddo

  where(abs(patchf(ixOmin1:ixOmax1))==1)
   !------- Pressure -------!
   pCD(ixOmin1:ixOmax1)      = Fhll(ixOmin1:ixOmax1,&
      mom(idim))- (Fhll(ixOmin1:ixOmax1,e_)+Fhll(ixOmin1:ixOmax1,&
      d_))* lambdaCD(ixOmin1:ixOmax1)
   !---- Tau Right ---!
   wCD(ixOmin1:ixOmax1,e_) = (cspeed(ixOmin1:ixOmax1) * wSub(ixOmin1:ixOmax1,&
      e_) + vSub(ixOmin1:ixOmax1) * wSub(ixOmin1:ixOmax1,&
      d_) - wSub(ixOmin1:ixOmax1,mom(idim))+&
      lambdaCD(ixOmin1:ixOmax1)*pCD(ixOmin1:ixOmax1))/(cspeed(ixOmin1:ixOmax1)-&
      lambdaCD(ixOmin1:ixOmax1))
   !---- Sx_ Right ---!
   wCD(ixOmin1:ixOmax1,mom(idim)) =(wCD(ixOmin1:ixOmax1,&
      e_)+wCD(ixOmin1:ixOmax1,d_)+pCD(ixOmin1:ixOmax1))*lambdaCD(&
      ixOmin1:ixOmax1)
   !-- Not real lfac, but fill it anyway --!
   wCD(ixOmin1:ixOmax1,lfac_) = one
  endwhere

 ! In case of B_idim = 0, we need only the normal velocity
  do iw = mom(1),mom(ndir)
    if(iw /= mom(idim))Then
      where(abs(patchf(ixOmin1:ixOmax1))==1)
        vCD(ixOmin1:ixOmax1,iw) = zero
      endwhere
    else
      where(abs(patchf(ixOmin1:ixOmax1))==1)
        vCD(ixOmin1:ixOmax1,iw) = lambdaCD(ixOmin1:ixOmax1)
      endwhere
    endif
  enddo

endif bxnozero

do iw=1,nwflux
 if(iw /= mag(idim).and. iw/=psi_)then
   where(abs(patchf(ixOmin1:ixOmax1))==1)
       ! f_i=fSub+lambda (wCD-wSub)
       f(ixOmin1:ixOmax1,iw)=fSub(ixOmin1:ixOmax1,&
          iw)+cspeed(ixOmin1:ixOmax1)*(wCD(ixOmin1:ixOmax1,&
          iw)-wsub(ixOmin1:ixOmax1,iw))
   endwhere
 else
       f(ixOmin1:ixOmax1,iw)=zero
 end if
end do

  end subroutine srmhd_get_wCD

end module mod_srmhd_hllc
