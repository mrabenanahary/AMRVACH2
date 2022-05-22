module mod_grackle_parameters
  use mod_constants
  use mod_global_parameters
  use mod_obj_global_parameters
  use grackle_header

!=================================================================================

! Define constants

real*8,parameter    :: mh=mp_cgs, kboltz=kB_cgs
real(kind=gr_rpknd) :: tiny_number= 1.0e-20_gr_RKIND,huge_number= 1.0e+20_gr_RKIND


! Grackle parameters for this solver time step
type gr_params
  real*8               :: initial_redshift
  real*8               :: temperature_units, pressure_units, dtchem
  real*8               :: grid_dx
  integer              :: number_of_objects
  real(kind=dp)   :: He_abundance
  real(kind=dp)   :: mean_nall_to_nH
  real(kind=dp)   :: mean_mass
  real(kind=dp)   :: mean_mup
  real(kind=dp)   :: mean_ne_to_nH

  !     Grid size and dimension
  !     grid_start and grid_end are used to ignore ghost zones.
  !     grid_dx is used in H2 self-shielding approximation only
  integer                          :: field_size(ndim)
  !real(kind=gr_rpknd) :: density(1:2)
  !real(kind=gr_rpknd) :: energy(1:2) !=internal energy
  !INTEGER                          :: grid_rank(1:2)
  !INTEGER                          :: grid_dimension(1:2,3)
  !INTEGER                          :: grid_start(1:2,3)
  !INTEGER                          :: grid_end(1:2,3)
  CHARACTER(LEN=257)   :: data_file
end type gr_params

! Grackle fields
type gr_fields
  character(len=80):: filename
  integer          :: field_size(ndim)
  real(kind=gr_rpknd), allocatable :: gr_density(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_energy(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_x_velocity(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_y_velocity(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_z_velocity(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_HI_density(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_HII_density(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_HM_density(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_HeI_density(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_HeII_density(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_HeIII_density(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_H2I_density(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_H2II_density(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_DI_density(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_DII_density(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_HDI_density(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_e_density(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_metal_density(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_dust_density(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_volumetric_heating_rate(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_specific_heating_rate(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_RT_HI_ionization_rate(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_RT_HeI_ionization_rate(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_RT_HeII_ionization_rate(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_RT_H2_dissociation_rate(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_RT_heating_rate(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_cooling_time(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_gamma(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_pressure(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_temperature(:^D&)
  real(kind=gr_rpknd), allocatable :: gr_dust_temperature(:^D&)
  !     Grid size and dimension
  !     grid_start and grid_end are used to ignore ghost zones.
  !     grid_dx is used in H2 self-shielding approximation only
  INTEGER                          :: grid_rank
  INTEGER                          :: grid_dimension(3)
  INTEGER                          :: grid_start(3)
  INTEGER                          :: grid_end(3)
end type gr_fields

end module mod_grackle_parameters
