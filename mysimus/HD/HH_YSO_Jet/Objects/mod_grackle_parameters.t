module mod_grackle_parameters
  use mod_constants
  use mod_global_parameters
  use mod_obj_global_parameters
  use grackle_header

!=================================================================================

! Define constants

real*8,parameter    :: mh_gr=mp_cgs, kboltz_gr=kB_cgs, me_gr=const_me
real(kind=dp) :: tiny_number= 1.0d-20,huge_number= 1.0d20
integer,parameter :: max_num_parameters = 5
real(kind=dp),parameter :: abundance_tolerance = 1.0e-2_gr_RKIND
real(kind=dp),parameter :: abundance_density = 1.0e-2_gr_RKIND
real(kind=dp),parameter :: abundance_relative_density = 5.0e-2_gr_RKIND
real(kind=dp),parameter :: abundance_sum = 1.0e-1_gr_RKIND
real(kind=dp),parameter :: chiD_default = 0.0000340249732407893951507840462224
real(kind=dp),parameter :: xHI_default = 0.999965975041198345221501



! species weight
real(kind=dp),parameter :: w_HI = 1.0d0 ! 1
real(kind=dp),parameter :: w_HII = 1.0d0 ! 1
real(kind=dp),parameter :: w_HM = 1.0d0 ! 1
real(kind=dp),parameter :: w_HeI = 4.0d0 ! 3.971525936989
real(kind=dp),parameter :: w_HeII = 4.0d0 ! 3.971525936989
real(kind=dp),parameter :: w_HeIII = 4.0d0 ! 3.971525936989
real(kind=dp),parameter :: w_H2I = 2.0d0 ! 2.000000000004
real(kind=dp),parameter :: w_H2II = 2.0d0 ! 2.000000000004
real(kind=dp),parameter :: w_DI = 2.0d0 ! 1.998463735368
real(kind=dp),parameter :: w_DII = 2.0d0 ! 1.998463735368
real(kind=dp),parameter :: w_HDI = 3.0d0 ! 2.998463735326
real(kind=dp),parameter :: w_e = 5.443205831355435e-4

! Grackle parameters for this solver time step


end module mod_grackle_parameters
