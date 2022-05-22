module mod_grackle_chemistry_solver
  use mod_constants
  use mod_global_parameters
  use mod_physics
  use mod_obj_dust, only : dust
  use mod_obj_global_parameters
  use mod_obj_usr_unit
  use mod_obj_mat
  use grackle_header
  use mod_grackle_chemistry




contains

subroutine set_mean_mup(self,densitymethod,iobject)
  use mod_global_parameters
  implicit none
  class(gr_objects),TARGET                         :: self
  CHARACTER(LEN=30), intent(in)                    :: densitymethod
  integer                    :: iobject
  ! local -------
  integer :: n_species
  integer :: n_coeff_tot
  integer :: n_coeff_nH_tot

  !----------------------------------------------------------

  n_species = 0
  n_coeff_tot = 0
  select case(densitymethod)
    case('chemical_coefficient')
      if(self%myconfig%gr_primordial_chemistry >= 1)then
        n_species = 6
        n_coeff_tot = gr_HI_density(gr_patches_indices_global(iobject))
        n_coeff_tot = n_coeff_tot + gr_HII_density(gr_patches_indices_global(iobject))
        n_coeff_tot = n_coeff_tot +gr_e_density(gr_patches_indices_global(i_all))
        ! Densities as fraction of n_H=n(H+H^+)+2n(H_2)
        gr_HeI_density(gr_patches_indices_global(i_all)) = 0.1
        gr_HeII_density(gr_patches_indices_global(i_all)) = 0.
        gr_HeIII_density(gr_patches_indices_global(i_all)) = 0.

        n_coeff_tot = n_coeff_tot + gr_HM_density(gr_patches_indices_global(i_all)) = 0.
        gr_H2I_density(gr_patches_indices_global(i_all)) = 2.0
        gr_H2II_density(gr_patches_indices_global(i_all)) = 0.
        gr_DI_density(gr_patches_indices_global(i_all)) = 0.
        gr_DII_density(gr_patches_indices_global(i_all)) = 0.
        gr_HDI_density(gr_patches_indices_global(i_all)) = 0.

        ! metal density : treat separately ==> density = chem species + metals + dust
        ! Densities in the fraction of ISM/JET in fraction of rho_tot = mean_mass * n_H
        gr_metal_density(gr_patches_indices_global(i_all)) = 0.
        ! dust density : treat separately ==> density = chem species + metals + dust
        ! Densities in the fraction of ISM/JET in fraction of rho_tot = mean_mass * n_H
        gr_dust_density(gr_patches_indices_global(i_all)) = 0.
        if(self%myconfig%gr_primordial_chemistry >= 2)then
          n_species = n_species + 3
          if(self%myconfig%gr_primordial_chemistry >= 3)then
            n_species = n_species + 3
          end if
        end if
      end if
    case default
      write(*,*) 'Method ',densitymethod, ' unknown!'
      call mpistop('The code stops now!')
  end select

  if(He_abundance_sub>0.0_dp)then
    mean_nall_to_nH = (1.0_dp+2.0_dp*He_abundance_sub)
    mean_mass       = (1.0_dp+4.0_dp*He_abundance_sub)
    select case(trim(hd_chemical_gas_type))
      case('fullymolecular')
        mean_mup = mean_mass/(0.5_dp+He_abundance_sub)
        mean_ne_to_nH =0.0_dp  ! is the value used in implimentation of low temperarure cooling table DM2
      case('fullyatomic')
        mean_mup = mean_mass/(1.0_dp+He_abundance_sub)
        mean_ne_to_nH =0.0_dp ! is the value used in implimentation of low temperarure cooling table DM2
      case('ionised')
        mean_mup = mean_mass/(2.0_dp+2.0_dp*He_abundance_sub)
        mean_ne_to_nH = (1.0_dp+He_abundance_sub)
      case('fullyionised')
        mean_mup = mean_mass/(2.0_dp+3.0_dp*He_abundance_sub)
        mean_ne_to_nH = (1.0_dp+2.0_dp*He_abundance_sub)
     case default
       write(*,*) 'The chemical gas type : ', trim(hd_chemical_gas_type)
        write (*,*) "Undefined gas chemical type entered in mod_hd_phys.t "
        call mpistop('The stops at hd_fill_chemical_ionisation in src/hd/mod_hd_phys.t')
    end select
    !hd_config%mean_mup          = (2.0_dp+3.0_dp*He_abundance_sub)

  else if(dabs(He_abundance_sub)<=smalldouble)then
    mean_mass         = 1.0_dp
    mean_mup          = 1.0_dp
    mean_ne_to_nH     = 1.0_dp


end subroutine set_mean_mup



end module mod_grackle_chemistry_solver
