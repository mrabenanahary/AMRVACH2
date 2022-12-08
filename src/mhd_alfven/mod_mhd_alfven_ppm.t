module mod_mhd_alfven_ppm
  use mod_mhd_alfven_phys

  implicit none
  private

  public :: mhd_alfven_ppm_init

contains

  subroutine mhd_alfven_ppm_init()
    use mod_physics_ppm

    phys_ppm_flatcd => mhd_alfven_ppm_flatcd
    phys_ppm_flatsh => mhd_alfven_ppm_flatsh
  end subroutine mhd_alfven_ppm_init

  subroutine mhd_alfven_ppm_flatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dpres)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L,ixL^L,ixR^L
    double precision, intent(in)    :: w(ixI^S,nw),d2w(ixI^S,1:nwflux)
    double precision, intent(inout) :: drho(ixI^S),dpres(ixI^S)

    if(mhd_alfven_energy) then
      drho(ixO^S) =mhd_alfven_gamma*dabs(d2w(ixO^S,rho_))&
           /min(w(ixL^S,rho_),w(ixR^S,rho_))
      dpres(ixO^S) = dabs(d2w(ixO^S,p_))/min(w(ixL^S,p_),w(ixR^S,p_))
    else
      call mpistop("PPM with flatcd=.true. can not be used without energy equation!")
    end if
  end subroutine mhd_alfven_ppm_flatcd

  ! based on Mignone and Miller and Collela 2002
  ! PPM flattening at shocks: we use total pressure and not thermal pressure
  subroutine mhd_alfven_ppm_flatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dpres,dv)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(inout) :: drho(ixI^S),dpres(ixI^S),dv(ixI^S)
    double precision                :: ptot(ixI^S)

    if(mhd_alfven_energy) then
      ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
      ptot(ixO^S)=w(ixO^S,p_)+half*sum(w(ixO^S,mag(:))**2,dim=ndim+1)
      where (dabs(ptot(ixRR^S)-ptot(ixLL^S))>smalldouble)
         drho(ixO^S) = dabs((ptot(ixR^S)-ptot(ixL^S))&
              /(ptot(ixRR^S)-ptot(ixLL^S)))
      elsewhere
         drho(ixO^S) = zero
      end where

      !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
      !  use "dpres" to save squared sound speed, assume primitive in w
      dpres(ixO^S)=(mhd_alfven_gamma*w(ixO^S,p_)/w(ixO^S,rho_))

      dpres(ixO^S)  = dabs(ptot(ixR^S)-ptot(ixL^S))&
           /(w(ixO^S,rho_)*dpres(ixO^S))
      ! recycle ptot to store v
      ptot(ixI^S)= w(ixI^S,mom(idims))
      call gradient(ptot,ixI^L,ixO^L,idims,dv)
    else
      call mpistop("PPM with flatsh=.true. can not be used without energy equation!")
    end if

  end subroutine mhd_alfven_ppm_flatsh

end module mod_mhd_alfven_ppm
