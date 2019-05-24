module mod_srmhd_eos
  use mod_global_parameters
  use mod_srmhd_parameters
  implicit none
  !> gamma minus one and its inverse
  real(kind=dp)   , save :: srmhd_gamma_1, inv_srmhd_gamma_1,&
     gamma_to_srmhd_gamma_1

 contains

    !>compute the enthalpy
  subroutine srmhd_get_enthalpy(ixOmin1,ixOmin2,ixOmax1,ixOmax2,rho,p,rhoh)
  ! made by Z. MELIANI 14/02/2018
    implicit none
    integer, intent(in)                :: ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)       :: rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
       p(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    real(kind=dp)   , intent(inout)    :: rhoh(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    ! .. local ..
    real(kind=dp)   , dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2) :: E_th,E

    if(srmhd_eos) then
     E_th = p*inv_srmhd_gamma_1
     E    = E_th+dsqrt(E_Th**2.0_dp+rho**2.0_dp)
     rhoh = half*((srmhd_gamma+one) * E-srmhd_gamma_1* rho*(rho/E))
    else
     rhoh = (rho+gamma_to_srmhd_gamma_1*p)
    end if
  end subroutine srmhd_get_enthalpy
    !>compute the enthalpy from scalar
  subroutine srmhd_get_enthalpy_scalar(rho,p,rhoh)
  ! made by Z. MELIANI 14/02/2018
    implicit none

    real(kind=dp)   , intent(in)       :: rho,p
    real(kind=dp)   , intent(inout)    :: rhoh
    ! .. local ..
    real(kind=dp)                      :: E_th,E

    if(srmhd_eos) then
     E_th = p*inv_srmhd_gamma_1
     E    = E_th+dsqrt(E_Th**2.0_dp+rho**2.0_dp)
     rhoh = half*((srmhd_gamma+one) * E-srmhd_gamma_1* rho*(rho/E))
    else
     rhoh = (rho+gamma_to_srmhd_gamma_1*p)
    end if
  end subroutine srmhd_get_enthalpy_scalar

  !> Calculate thermal pressure for enthalpy and density
  subroutine srmhd_get_pressure_primitive_eos(ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,rho,rhoh,pth)
  ! made by Z. MELIANI 14/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)            :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)   :: rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
       rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    real(kind=dp)   , intent(out)  :: pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    real(kind=dp)                  :: E(ixOmin1:ixOmax1,ixOmin2:ixOmax2)


    cond_iseos : if(srmhd_eos) then
     E = (rhoh+dsqrt(rhoh**2.0d0+(srmhd_gamma**2.0d0-&
        1.0_dp)*rho**2.0d0))/(srmhd_gamma+1.0_dp)
     pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.5d0*srmhd_gamma_1* &
        (E-rho*(rho/E))
    else cond_iseos
     pth = (rhoh-rho)/gamma_to_srmhd_gamma_1
    end if cond_iseos
  end subroutine srmhd_get_pressure_primitive_eos

 !==========================================================================
 !> Calculate thermal pressure for enthalpy and density
 subroutine srmhd_get_pressure_primitive_scalar(rho,rhoh,pth)
 ! made by Z. MELIANI 14/02/2018
   use mod_global_parameters
   implicit none
   real(kind=dp)   , intent(in)   :: rho,rhoh
   real(kind=dp)   , intent(out)  :: pth
   ! .. local ..
   real(kind=dp)                  :: E
   !--------------------------------------------

   cond_iseos : if(srmhd_eos) then
    E = (rhoh+dsqrt(rhoh**2.0d0+(srmhd_gamma**2.0d0-&
       1.0_dp)*rho**2.0d0))/(srmhd_gamma+1.0_dp)
    pth = 0.5d0*srmhd_gamma_1* (E-rho*(rho/E))
   else cond_iseos
    pth = (rhoh-rho)/gamma_to_srmhd_gamma_1
   end if cond_iseos
 end subroutine srmhd_get_pressure_primitive_scalar
 !==========================================================================
  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine srmhd_get_pthermal_eos(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,x,rho,rhoh,e_in,pth)
  ! made by Z. MELIANI 14/02/2018
    use mod_global_parameters
    implicit none

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in) :: rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
       rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2),e_in(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    real(kind=dp)   , intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(out):: pth(ixImin1:ixImax1,ixImin2:ixImax2)

    real(kind=dp)                :: e(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    is_notdiab: if(srmhd_energy) then

      is_einternal : if(block%e_is_internal) then
        is_eos_einter : if(srmhd_eos)then
         pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
            0.5d0*srmhd_gamma_1*(e_in-rho*(rho/e_in))
        else is_eos_einter
         pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=srmhd_gamma_1*e_in
        end if is_eos_einter
      else is_einternal
       is_eos : if(srmhd_eos) then
         e  = (rhoh+dsqrt(rhoh**2.0d0+(srmhd_gamma**2.0d0-&
            one)*rho**2.0d0))/(srmhd_gamma+one)

         pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.5d0*srmhd_gamma_1* &
            (e-rho*(rho/e))
       else is_eos
         pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (rhoh - &
            rho)/gamma_to_srmhd_gamma_1
       end if is_eos

      end if is_einternal

    else is_notdiab

      pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=srmhd_adiab*rho**srmhd_gamma

    end if is_notdiab
  end subroutine srmhd_get_pthermal_eos


  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamm*p/rho
  subroutine srmhd_get_csound2_eos(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,x,rho,rhoh,csound2)
  ! made by Z. MELIANI 14/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)    :: rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
       rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    real(kind=dp)   , intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    real(kind=dp)   , intent(out)   :: csound2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    real(kind=dp)                   :: h(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
       E(ixOmin1:ixOmax1,ixOmin2:ixOmax2),p(ixOmin1:ixOmax1,ixOmin2:ixOmax2)


    if(srmhd_energy) then
      if(srmhd_eos) then
       E(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(rhoh(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)+dsqrt(rhoh(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)**2.0_dp+(srmhd_gamma**2.0_dp-&
          one)*rho**2.0_dp))/(srmhd_gamma+one)

       p=0.5_dp*(E-rho*(rho/E))

       csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(p*((srmhd_gamma+one)+&
          srmhd_gamma_1*(rho/E)**2.0d0))/(2.0_dp*rhoh)

      else
    !       csound2(ixO^S)= (eqpar(gamma_)-one) &
    !          *(w(ixO^S,xi_)-w(ixO^S,d_)*w(ixO^S,lfac_))/w(ixO^S,xi_)
    !   p =   (rhoh - rho)/gamma_to_srmhd_gamma_1
       csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(1.0_dp - &
          rho/rhoh)*srmhd_gamma_1


      end if
    else
      csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=srmhd_gamma*srmhd_adiab*rho**&
         srmhd_gamma_1
    end if
  end subroutine srmhd_get_csound2_eos


  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamm*p/rho
  subroutine srmhd_get_csound2_prim_eos(ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,x,rho,rhoh,p,csound2)
  ! made by Z. MELIANI 14/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)    :: rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
       rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2),p(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    real(kind=dp)   , intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    real(kind=dp)   , intent(out)   :: csound2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    real(kind=dp)                   :: h(ixImin1:ixImax1,ixImin2:ixImax2),&
       E(ixImin1:ixImax1,ixImin2:ixImax2)
    if(srmhd_energy) then
      if(srmhd_eos) then
       E(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(rhoh(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)+dsqrt(rhoh(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)**2.0d0+(srmhd_gamma**2.0d0-&
          one)*rho**2.0d0))/(srmhd_gamma+one)
       csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(p*((srmhd_gamma+one)+&
          srmhd_gamma_1*(rho/E)**2.0d0))/(2.0*rhoh)
      else
       csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=srmhd_gamma*p/rhoh
      end if
    else
      csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=srmhd_gamma*srmhd_adiab*rho**&
         srmhd_gamma_1
    end if
  end subroutine srmhd_get_csound2_prim_eos

  !> Calculate the Enthalpy and dhdp from pressure
  subroutine srmhd_get_val_h_dhdp(rho,p,drhodp,h,dhdp)
  ! made by Z. MELIANI 14/02/2018
    use mod_global_parameters
    implicit none

    real(kind=dp)   ,           intent(in) :: rho,p,drhodp
    real(kind=dp)   ,           intent(out):: h
    real(kind=dp)   , optional, intent(out):: dhdp

    real(kind=dp)                :: Eth,E,dEthdp,dEdp,sqrtEt2rho2,rhotoE

    is_energy : if(srmhd_energy) then
      is_eos : if(srmhd_eos) then

       Eth = p*inv_srmhd_gamma_1
       sqrtEt2rho2 = dsqrt(Eth**2.0d0+rho**2.0d0)
       E = (Eth + sqrtEt2rho2)
       rhotoE= rho/E
       h = 0.5 *((srmhd_gamma+one) * E-srmhd_gamma_1 * rho*rhotoE)

       if(present(dhdp))then
        dEthdp = inv_srmhd_gamma_1
        dEdp   = dEthdp +(dEthdp*Eth+rho*drhodp)/sqrtEt2rho2
        dhdp   = 0.5 *(((srmhd_gamma+one) + rhotoE**2.0)* dEdp- &
           2.0*srmhd_gamma_1*rhotoE* drhodp)
       end if
      else is_eos
       h=rho+gamma_to_srmhd_gamma_1*p
       if(present(dhdp))dhdp=drhodp+gamma_to_srmhd_gamma_1
      end if is_eos
    else is_energy
    end if is_energy
  end  subroutine srmhd_get_val_h_dhdp

   !> Calculate the pressure and dpdxi from xi
  subroutine srmhd_get_val_p_dpdxi(rho,h,drhodxi,dhdxi,p,dpdxi)
  ! made by Z. MELIANI 14/02/2018
    use mod_global_parameters
    implicit none

    real(kind=dp)   ,           intent(in) :: rho,h,drhodxi,dhdxi
    real(kind=dp)   ,           intent(out):: p
    real(kind=dp)   ,           intent(out):: dpdxi

    real(kind=dp)                          :: Eth,E,dEdxi,rhotoE

    is_energy : if(srmhd_energy) then
      is_eos : if(srmhd_eos) then
       E =(h+dsqrt(h**2.0d0+(srmhd_gamma**2.0d0-one)*rho**2.0d0))/(srmhd_gamma+&
          one)
       rhotoE = rho/E
       ! output pressure
       p=srmhd_gamma_1/2.0_dp* (E-rho*rhotoE)

       dEdxi=(dhdxi+(h*dhdxi+(srmhd_gamma**2.0d0-&
          one)*rho*drhodxi)/dsqrt(h**2.0d0+&
          (srmhd_gamma**2.0d0-one)*rho**2.0d0))/ (srmhd_gamma+one)
       dpdxi=srmhd_gamma_1/2.0d0*((1.0+rhotoE**2.0)*dEdxi-2.0*rhotoE*drhodxi)
      else is_eos
       p = (h - rho)/gamma_to_srmhd_gamma_1
       dpdxi = (dhdxi-drhodxi)/gamma_to_srmhd_gamma_1
      end if is_eos
    else is_energy

!
    end if is_energy
  end  subroutine srmhd_get_val_p_dpdxi


  subroutine srmhd_get_h_noflux(rho,h_p,h)
  ! made by Z. MELIANI 14/02/2018
    use mod_global_parameters
    implicit none

    real(kind=dp)   ,           intent(in) :: rho,h_p

    real(kind=dp)   ,           intent(out):: h

    is_energy : if(srmhd_energy) then
      is_eos : if(srmhd_eos) then
       h=0.5*((srmhd_gamma+1.0)*h_p-srmhd_gamma_1*rho*(rho/h_p))
      else is_eos
       h=rho+(h_p-rho)*srmhd_gamma
      end if is_eos
    else is_energy

!
    end if is_energy

  end subroutine srmhd_get_h_noflux
end module mod_srmhd_eos
