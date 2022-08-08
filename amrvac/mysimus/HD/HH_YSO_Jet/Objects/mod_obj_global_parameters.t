module mod_obj_global_parameters
 use mod_constants
 use mod_physics
 use mod_srmhd
 implicit none
 integer,public  :: itr,phys_n_tracer
   type unit_expand
   real(dp)  :: volum
   real(dp)  :: mass
   real(dp)  :: mass_flux
   real(dp)  :: energy_flux
   real(dp)  :: luminosity
   real(dp)  :: energy
   end type  unit_expand

   type(unit_expand),public :: unit_user

   logical,public   :: use1Dfile
   ! physical constantes:
   type const_expand
   real(dp)  :: G       = 6.6743d-8 !OLD : 6.67259D-8, cm^3 g^-1 s^-2
   real(dp)  :: clight  = 2.99792458d10   ! cm s^-1
   real(dp)  :: z_solar = 0.012_dp
   end type const_expand
   type(const_expand),public :: constusr

   ! solar constantes
   real(dp),public :: solar_mass       = 1.9892d+33 !> solar mass cgs
   real(dp),public :: solar_luminosity = 3.826d+33  !> solar luminosity (erg/s)
   real(dp),public :: solar_radius     = 6.95987d+10!> solar radius (cm)
   ! geometry
   character(len=30),public:: coordinate_system
contains

  !--------------------------------------------------------------------
    subroutine usr_physical_unit

     unit_user%volum          = unit_length**3.0_dp
     unit_user%mass           = unit_density*unit_user%volum
     unit_user%energy         = unit_velocity**2.0_dp*unit_density
     unit_user%luminosity     = unit_length**2.0_dp*unit_user%mass&
                                *unit_time**(-3.0_dp)!*2.81166770085091
     unit_user%mass_flux      = unit_user%mass/ unit_time
     unit_user%energy_flux    = unit_user%mass*unit_time**(-3.0_dp)



    end subroutine usr_physical_unit
  !------------------------------------------------------------------------
end module mod_obj_global_parameters
