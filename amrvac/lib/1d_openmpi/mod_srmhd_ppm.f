module mod_srmhd_ppm
  use mod_srmhd_parameters
  use mod_srmhd_phys
  implicit none
  private

  public :: srmhd_ppm_init

contains

  subroutine srmhd_ppm_init()
    use mod_physics_ppm

    phys_ppm_flatcd => srmhd_ppm_flatcd
    phys_ppm_flatsh => srmhd_ppm_flatsh
  end subroutine srmhd_ppm_init

  subroutine srmhd_ppm_flatcd(ixImin1,ixImax1,ixOmin1,ixOmax1,ixLmin1,ixLmax1,&
     ixRmin1,ixRmax1,w,d2w,drho,dpressure)
  ! made by Z. MELIANI 14/02/2018
    use mod_global_parameters
    use mod_srmhd_phys

    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1,ixLmin1,&
       ixLmax1,ixRmin1,ixRmax1
    double precision, intent(in)    :: w(ixImin1:ixImax1,nw),&
       d2w(ixImin1:ixImax1,1:nwflux)
    double precision, intent(inout) :: drho(ixImin1:ixImax1),&
       dpressure(ixImin1:ixImax1)

    double precision                ::  pmag(ixImin1:ixImax1),&
       dpressuremag(ixImin1:ixImax1)

    if(srmhd_energy) then
      call srmhd_get_p_mag(w,ixImin1,ixImax1,ixImin1,ixImax1,pmag)
      dpressuremag(ixOmin1:ixOmax1) = pmag(ixLmin1:ixLmax1)-&
         pmag(ixRmin1:ixRmax1)

      drho(ixOmin1:ixOmax1) =dabs(d2w(ixOmin1:ixOmax1,&
         rho_))/min(w(ixLmin1:ixLmax1,rho_),w(ixRmin1:ixRmax1,rho_))

      dpressure(ixOmin1:ixOmax1) = dabs(d2w(ixOmin1:ixOmax1,&
         p_)+dpressuremag(ixOmin1:ixOmax1))/min(w(ixLmin1:ixLmax1,p_),&
         w(ixRmin1:ixRmax1,p_))
    else
      call mpistop("PPM with flatcd=.true. can not be used without energy &
         equation !")
    end if
  end subroutine srmhd_ppm_flatcd

  ! based on Mignone and Miller and Collela 2002
  ! PPM flattening at shocks: we use total pressure and not thermal pressure
  subroutine srmhd_ppm_flatsh(ixImin1,ixImax1,ixOmin1,ixOmax1,ixLLmin1,&
     ixLLmax1,ixLmin1,ixLmax1,ixRmin1,ixRmax1,ixRRmin1,ixRRmax1,idims,w,drho,&
     dpressure,dv)
  ! made by Z. MELIANI 14/02/2018
    use mod_global_parameters
    use mod_srmhd_phys
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1,&
       ixLLmin1,ixLLmax1,ixLmin1,ixLmax1,ixRmin1,ixRmax1,ixRRmin1,ixRRmax1
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,nw)
    double precision, intent(inout) :: drho(ixImin1:ixImax1),&
       dpressure(ixImin1:ixImax1),dv(ixImin1:ixImax1)
    double precision                :: ptot(ixImin1:ixImax1)

    if(srmhd_energy) then
      ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
      call srmhd_get_p_total(w,pw(saveigrid)%x,ixImin1,ixImax1,ixImin1,ixImax1,&
         ptot)
      where (dabs(ptot(ixRRmin1:ixRRmax1)-&
         ptot(ixLLmin1:ixLLmax1))>smalldouble)
         drho(ixOmin1:ixOmax1) = dabs((ptot(ixRmin1:ixRmax1)-&
            ptot(ixLmin1:ixLmax1))/(ptot(ixRRmin1:ixRRmax1)-&
            ptot(ixLLmin1:ixLLmax1)))
      elsewhere
         drho(ixOmin1:ixOmax1) = zero
      end where

      !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
      !  use "dpressure" to save squared sound speed, assume primitive in w
      call srmhd_get_csound_prim(w,pw(saveigrid)%x,ixImin1,ixImax1,ixOmin1,&
         ixOmax1,idims,dpressure)

      dpressure(ixOmin1:ixOmax1)  = dabs(ptot(ixRmin1:ixRmax1)-&
         ptot(ixLmin1:ixLmax1))/(w(ixOmin1:ixOmax1,&
         rho_)*dpressure(ixOmin1:ixOmax1))
      ! recycle ptot to store v
      ptot(ixImin1:ixImax1)= w(ixImin1:ixImax1,mom(idims))/w(ixImin1:ixImax1,&
         lfac_)
      call gradient(ptot,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,dv)
    else
      call mpistop("PPM with flatsh=.true. can not be used without energy &
         equation !")
    end if

  end subroutine srmhd_ppm_flatsh

end module mod_srmhd_ppm
