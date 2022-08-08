module mod_grackle_parameters
  use mod_constants
  use mod_global_parameters
  use mod_obj_global_parameters
  use grackle_header

!=================================================================================

! Define constants

real*8,parameter    :: mh_gr=mp_cgs, kboltz_gr=kB_cgs, me_gr=const_me
real(kind=gr_rpknd) :: tiny_number= 1.0e-20_gr_RKIND,huge_number= 1.0e+20_gr_RKIND
integer,parameter :: max_num_parameters = 5
real(kind=gr_rpknd),parameter :: abundance_tolerance = 1.0e-2_gr_RKIND
real(kind=gr_rpknd),parameter :: abundance_density = 1.0e-2_gr_RKIND
real(kind=gr_rpknd),parameter :: abundance_relative_density = 5.0e-2_gr_RKIND
real(kind=gr_rpknd),parameter :: abundance_sum = 1.0e-1_gr_RKIND
real(kind=gr_rpknd),parameter :: chiD_default = 0.0000340249732407893951507840462224
real(kind=gr_rpknd),parameter :: xHI_default = 0.999965975041198345221501



! species weight
real(kind=gr_rpknd),parameter :: w_HI = 1
real(kind=gr_rpknd),parameter :: w_HII = 1
real(kind=gr_rpknd),parameter :: w_HM = 1
real(kind=gr_rpknd),parameter :: w_HeI = 3.971525936989
real(kind=gr_rpknd),parameter :: w_HeII = 3.971525936989
real(kind=gr_rpknd),parameter :: w_HeIII = 3.971525936989
real(kind=gr_rpknd),parameter :: w_H2I = 2.000000000004
real(kind=gr_rpknd),parameter :: w_H2II = 2.000000000004
real(kind=gr_rpknd),parameter :: w_DI = 1.998463735368
real(kind=gr_rpknd),parameter :: w_DII = 1.998463735368
real(kind=gr_rpknd),parameter :: w_HDI = 2.998463735326
real(kind=gr_rpknd),parameter :: w_e = 5.443205831355435e-4

! Grackle parameters for this solver time step


! Grackle fields
type gr_fields
  character(len=257):: filename
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
