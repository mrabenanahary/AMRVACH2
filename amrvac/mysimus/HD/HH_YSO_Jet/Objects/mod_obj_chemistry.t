      module mod_obj_chemistry

      use mod_global_parameters
      use mod_obj_global_parameters

      USE ISO_C_BINDING
      implicit None
#include "grackle.def"

      integer iresult, i

c     Define constants

      real*8 mh, kboltz, fH
      R_PREC tiny_number
      parameter (tiny_number = 1.0e-20_RKIND)
      parameter (mh = 1.67262171d-24)
      parameter (kboltz = 1.3806504d-16)
      parameter (fH = 0.76)

c     Initialization parameters

      real*8 initial_redshift
      character(len=80), TARGET :: filename

c     Field data arrays

      integer field_size
      parameter (field_size = 5)

      real*8 temperature_units, pressure_units, dtchem

      R_PREC, TARGET :: density(field_size), energy(field_size),
     &     x_velocity(field_size), y_velocity(field_size),
     &     z_velocity(field_size),
     &     HI_density(field_size), HII_density(field_size),
     &     HM_density(field_size),
     &     HeI_density(field_size), HeII_density(field_size),
     &     HeIII_density(field_size),
     &     H2I_density(field_size), H2II_density(field_size),
     &     DI_density(field_size), DII_density(field_size),
     &     HDI_density(field_size),
     &     e_density(field_size), metal_density(field_size),
     &     dust_density_ch(field_size),
     &     volumetric_heating_rate(field_size),
     &     specific_heating_rate(field_size),
     &     RT_HI_ionization_rate(field_size),
     &     RT_HeI_ionization_rate(field_size),
     &     RT_HeII_ionization_rate(field_size),
     &     RT_H2_dissociation_rate(field_size),
     &     RT_heating_rate(field_size)

      R_PREC, TARGET :: cooling_time(field_size), gamma(field_size),
     &     pressure(field_size), temperature(field_Size),
     &     dust_temperature(field_Size)

c     Grid size and dimension
c     grid_start and grid_end are used to ignore ghost zones.
c     grid_dx is used in H2 self-shielding approximation only

      INTEGER, TARGET :: grid_rank, grid_dimension(3),
     &     grid_start(3), grid_end(3)

      real*8 grid_dx

c     Define storage for grackle units, fields and parameter data

      TYPE (grackle_units) :: my_units
      TYPE (grackle_field_data) :: my_fields
      TYPE (grackle_chemistry_data) :: grackle_data

      contains

      subroutine test_chemistry()

      PRINT *, 'Linking with obj_global_parameters'
      PRINT *, 'ndim='
      PRINT *, ndim
      PRINT *, 'mod_global_parameters:w_refine_weight='
      PRINT *, w_refine_weight
      PRINT *, 'mod_global_parameters:clight='
      PRINT *, constusr%clight


      end subroutine test_chemistry

      end module mod_obj_chemistry
