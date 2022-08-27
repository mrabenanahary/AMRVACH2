!> Hydrodynamics HLLC module
module mod_hd_hllc
  use mod_hd_grackle_phys
  use mod_dust
  implicit none
  private

  public :: hd_hllc_init

contains

  subroutine hd_hllc_init()
    use mod_physics_hllc

    phys_diffuse_hllcd => hd_diffuse_hllcd
    phys_get_lCD => hd_get_lCD
    phys_get_wCD => hd_get_wCD
  end subroutine hd_hllc_init

  subroutine hd_diffuse_hllcd(ixI^L,ixO^L,idim,wLC,wRC,fLC,fRC,patchf)

    ! when method is hllcd or hllcd1 then:

    ! this subroutine is to enforce regions where we AVOID HLLC
    ! and use TVDLF instead: this is achieved by setting patchf to 4 in
    ! certain regions. An additional input parameter is nxdiffusehllc
    ! which sets the size of the fallback region.

    use mod_global_parameters

    integer, intent(in)                                      :: ixI^L,ixO^L,idim
    double precision, dimension(ixI^S,1:nw), intent(in)      :: wRC,wLC
    double precision, dimension(ixI^S,1:nwflux),intent(in) :: fLC, fRC
    integer         , dimension(ixI^S), intent(inout)        :: patchf

    integer                                           :: ixOO^D,TxOO^L

    ! In a user-controlled region around any point with flux sign change between
    ! left and right, ensure fallback to TVDLF
    {do ixOO^D= ixO^LIM^D\}
    {
    TxOOmin^D= max(ixOO^D - nxdiffusehllc*kr(idim,^D), ixOmin^D);
    TxOOmax^D= min(ixOO^D + nxdiffusehllc*kr(idim,^D), ixOmax^D);
    \}
    if(abs(patchf(ixOO^D)) == 1 .or. abs(patchf(ixOO^D)) == 4)Then
       if(any(fRC(ixOO^D,1:nwflux)*fLC(ixOO^D,1:nwflux)<-smalldouble))Then
          where(Abs(patchf(TxOO^S))==1)
             patchf(TxOO^S) = 4
          endwhere
       endif
    endif
    {enddo^D&\}

  end subroutine hd_diffuse_hllcd

  subroutine hd_get_lCD(wLC,wRC,fLC,fRC,cmin,cmax,idim,ixI^L,ixO^L, &
       whll,Fhll,lambdaCD,patchf)

    ! Calculate lambda at CD and set the patchf to know the orientation
    ! of the riemann fan and decide on the flux choice
    ! We also compute here the HLL flux and w value, for fallback strategy

    ! In this nul version, we simply compute nothing and ensure TVDLF fallback

    use mod_global_parameters

    integer, intent(in)                                      :: ixI^L,ixO^L,idim
    double precision, dimension(ixI^S,1:nw), intent(in)      :: wLC,wRC
    double precision, dimension(ixI^S,1:nwflux), intent(in):: fLC,fRC
    double precision, dimension(ixI^S), intent(in)           :: cmax,cmin
    integer         , dimension(ixI^S), intent(inout)        :: patchf
    double precision, dimension(ixI^S,1:nwflux), intent(out) :: Fhll,whll
    double precision, dimension(ixI^S), intent(out)          :: lambdaCD

    logical         , dimension(ixI^S)     :: Cond_patchf
    double precision                       :: Epsilon
    integer                                :: iw

    !--------------------------------------------
    ! on entry, patch is preset to contain values from -2,1,2,4
    !      -2: take left flux, no computation here
    !      +2: take right flux, no computation here
    !      +4: take TVDLF flux, no computation here
    !       1: compute the characteristic speed for the CD

    Cond_patchf(ixO^S)=(abs(patchf(ixO^S))==1)

    do iw=1,nwflux
       where(Cond_patchf(ixO^S))
          !============= compute HLL flux ==============!
          Fhll(ixO^S,iw)= (cmax(ixO^S)*fLC(ixO^S,iw)-cmin(ixO^S)*fRC(ixO^S,iw) &
               + cmin(ixO^S)*cmax(ixO^S)*(wRC(ixO^S,iw)-wLC(ixO^S,iw)))/&
               (cmax(ixO^S)-cmin(ixO^S))
          !======== compute intermediate HLL state =======!
          whll(ixO^S,iw) = (cmax(ixO^S)*wRC(ixO^S,iw)-cmin(ixO^S)*wLC(ixO^S,iw)&
               +fLC(ixO^S,iw)-fRC(ixO^S,iw))/(cmax(ixO^S)-cmin(ixO^S))
       endwhere
    enddo

    ! deduce the characteristic speed at the CD
    where(Cond_patchf(ixO^S))
       lambdaCD(ixO^S)=whll(ixO^S,mom(idim))/whll(ixO^S,rho_)
    elsewhere
       lambdaCD(ixO^S)=  0.0_dp
    end where

    where(Cond_patchf(ixO^S))
       ! double check whether obtained speed is in between min and max speeds given
       ! and identify in which part of the Riemann fan the time-axis is
       where(cmin(ixO^S) < zero .and. lambdaCD(ixO^S)>zero&
            .and.lambdaCD(ixO^S)<cmax(ixO^S))
          patchf(ixO^S) = -1
       elsewhere(cmax(ixO^S) > zero .and. lambdaCD(ixO^S) < zero&
            .and.lambdaCD(ixO^S)>cmin(ixO^S))
          patchf(ixO^S) =  1
       elsewhere(lambdaCD(ixO^S) >= cmax(ixO^S) .or. &
            lambdaCD(ixO^S) <= cmin(ixO^S))
          lambdaCD(ixO^S) = zero
          ! we will fall back to HLL flux case in this degeneracy
          patchf(ixO^S) =  3
       endwhere
    endwhere

    where(patchf(ixO^S)== 3)
       Cond_patchf(ixO^S)=.false.
    end where

    ! handle the specific case where the time axis is exactly on the CD
    if(any(lambdaCD(ixO^S)==zero.and.Cond_patchf(ixO^S)))then
       ! determine which sector (forward or backward) of the Riemann fan is smallest
       ! and select left or right flux accordingly
       where(lambdaCD(ixO^S)==zero.and.Cond_patchf(ixO^S))
          where(-cmin(ixO^S)>=cmax(ixO^S))
             patchf(ixO^S) =  1
          elsewhere
             patchf(ixO^S) = -1
          endwhere
       endwhere
    endif


    ! eigenvalue lambda for contact is near zero: decrease noise by this trick
    if(flathllc)then
       Epsilon=1.0d-6
       where(Cond_patchf(ixO^S).and. &
            dabs(lambdaCD(ixO^S))/max(cmax(ixO^S),Epsilon)< Epsilon  .and. &
            dabs(lambdaCD(ixO^S))/max(dabs(cmin(ixO^S)),Epsilon)< Epsilon)
          lambdaCD(ixO^S) =  zero
       end where
    end if

  end subroutine hd_get_lCD

  subroutine hd_get_wCD(x,wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,cmin,cmax,&
       ixI^L,ixO^L,idim,f)
    ! compute the intermediate state U*
    ! only needed where patchf=-1/1

    ! reference Li S., JCP, 203, 2005, 344-357
    ! reference T. Miyoski, Kusano JCP, 2008, 2005.

    use mod_global_parameters
    use mod_dust
    integer, intent(in)                                      :: ixI^L,ixO^L,idim
    double precision, dimension(ixI^S,1:nw), intent(in)      :: wRC,wLC
    double precision, dimension(ixI^S,1:ndim), intent(in)    :: x
    double precision, dimension(ixI^S,1:nwflux), intent(in)  :: whll, Fhll
    double precision, dimension(ixI^S), intent(in)           :: lambdaCD
    double precision, dimension(ixI^S), intent(in)           :: cmax,cmin
    double precision, dimension(ixI^S,1:nwflux), intent(in)  :: fRC,fLC
    double precision, dimension(ixI^S,1:nwflux),intent(out)  :: f
    double precision, dimension(ixI^S,1:nw)                  :: wCD,wSub
    double precision, dimension(ixI^S,1:nwflux)              :: fSub
    double precision, dimension(ixI^S)                       :: vSub,cspeed,pCD
    integer         , dimension(ixI^S), intent(inout)        :: patchf

    integer                                      :: n, iw, ix^D

    !-------------- auxiliary Speed and array-------------!
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
       if(patchf(ix^D)==1) then
         cspeed(ix^D)=cmax(ix^D)
         vSub(ix^D)=wRC(ix^D,mom(idim))/wRC(ix^D,rho_)
         wSub(ix^D,:)=wRC(ix^D,:)
         fSub(ix^D,:)=fRC(ix^D,:)
       else if(patchf(ix^D)==-1) then
         cspeed(ix^D)=cmin(ix^D)
         vSub(ix^D)=wLC(ix^D,mom(idim))/wLC(ix^D,rho_)
         wSub(ix^D,:)=wLC(ix^D,:)
         fSub(ix^D,:)=fLC(ix^D,:)
       end if
    {end do\}

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      if(abs(patchf(ix^D))==1) then
        wCD(ix^D,rho_) = wSub(ix^D,rho_)&
                         *(cspeed(ix^D)-vSub(ix^D))/(cspeed(ix^D)-lambdaCD(ix^D))
        do n=1,hd_n_tracer
          iw = tracer(n)
          wCD(ix^D,iw) = wSub(ix^D,iw)*(cspeed(ix^D)-vSub(ix^D))&
             /(cspeed(ix^D)-lambdaCD(ix^D))
        end do


        ! Mialy s adding:
        if(hd_config%use_grackle)then
          wCD(ix^D,rhoHI_) = wSub(ix^D,rhoHI_)*(cspeed(ix^D)-vSub(ix^D))&
          /(cspeed(ix^D)-lambdaCD(ix^D))
          wCD(ix^D,rhoHII_) = wSub(ix^D,rhoHII_)*(cspeed(ix^D)-vSub(ix^D))&
          /(cspeed(ix^D)-lambdaCD(ix^D))
          wCD(ix^D,rhoHM_) = wSub(ix^D,rhoHM_)*(cspeed(ix^D)-vSub(ix^D))&
          /(cspeed(ix^D)-lambdaCD(ix^D))
          wCD(ix^D,rhoHeI_) = wSub(ix^D,rhoHeI_)*(cspeed(ix^D)-vSub(ix^D))&
          /(cspeed(ix^D)-lambdaCD(ix^D))
          wCD(ix^D,rhoHeII_) = wSub(ix^D,rhoHeII_)*(cspeed(ix^D)-vSub(ix^D))&
          /(cspeed(ix^D)-lambdaCD(ix^D))
          wCD(ix^D,rhoHeIII_) = wSub(ix^D,rhoHeIII_)*(cspeed(ix^D)-vSub(ix^D))&
          /(cspeed(ix^D)-lambdaCD(ix^D))
          wCD(ix^D,rhoH2I_) = wSub(ix^D,rhoH2I_)*(cspeed(ix^D)-vSub(ix^D))&
          /(cspeed(ix^D)-lambdaCD(ix^D))
          wCD(ix^D,rhoH2II_) = wSub(ix^D,rhoH2II_)*(cspeed(ix^D)-vSub(ix^D))&
          /(cspeed(ix^D)-lambdaCD(ix^D))
          wCD(ix^D,rhoe_) = wSub(ix^D,rhoe_)*(cspeed(ix^D)-vSub(ix^D))&
          /(cspeed(ix^D)-lambdaCD(ix^D))
          wCD(ix^D,rhoDI_) = wSub(ix^D,rhoDI_)*(cspeed(ix^D)-vSub(ix^D))&
          /(cspeed(ix^D)-lambdaCD(ix^D))
          wCD(ix^D,rhoDII_) = wSub(ix^D,rhoDII_)*(cspeed(ix^D)-vSub(ix^D))&
          /(cspeed(ix^D)-lambdaCD(ix^D))
          wCD(ix^D,rhoHDI_) = wSub(ix^D,rhoHDI_)*(cspeed(ix^D)-vSub(ix^D))&
          /(cspeed(ix^D)-lambdaCD(ix^D))
          wCD(ix^D,rhometal_) = wSub(ix^D,rhometal_)*(cspeed(ix^D)-vSub(ix^D))&
          /(cspeed(ix^D)-lambdaCD(ix^D))
          wCD(ix^D,rhodust_) = wSub(ix^D,rhodust_)*(cspeed(ix^D)-vSub(ix^D))&
          /(cspeed(ix^D)-lambdaCD(ix^D))
        end if

        !------- Momentum ------!
        do iw=1, ndir
          if(iw /= idim)then
            ! eq. 21 22
            wCD(ix^D,mom(iw))=(cspeed(ix^D)*wSub(ix^D,mom(iw))-fSub(ix^D,mom(iw)))/&
                (cspeed(ix^D)-lambdaCD(ix^D))
          else
            ! eq. 20
            wCD(ix^D,mom(iw)) =  wCD(ix^D,rho_) * lambdaCD(ix^D)
          endif
        enddo
        if(hd_energy) then
          pCD(ix^D) = wSub(ix^D,rho_)*(cspeed(ix^D)-vSub(ix^D))&
                        *(lambdaCD(ix^D)-vSub(ix^D))&
                        +fSub(ix^D,mom(idim))-wSub(ix^D,mom(idim))*vSub(ix^D)
          ! Eq 31
          wCD(ix^D,e_) = (cspeed(ix^D)*wSub(ix^D,e_) &
                          -fSub(ix^D,e_)+lambdaCD(ix^D)*pCD(ix^D))&
                          /(cspeed(ix^D)-lambdaCD(ix^D))
        end if

        Loop_iw_flux_gas: do iw=1,hd_config%nwhllc
          ! f_i=fsub+lambda (wCD-wSub)
          f(ix^D,iw)=fsub(ix^D,iw)+cspeed(ix^D)*(wCD(ix^D,iw)-wSub(ix^D,iw))
        end do Loop_iw_flux_gas

      end if
    {end do\}
    if(hd_config%dust_on)then
      call dust_get_wCD(x,wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,cmin,cmax,&
       ixI^L,ixO^L,idim,f)
    end if
  end subroutine hd_get_wCD

end module mod_hd_hllc
