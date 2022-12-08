module mod_grackle_parameters
  use mod_global_parameters
  use mod_obj_global_parameters
  use grackle_header


  ! Definitions of the parameters necessary for Grackle and to link with AMRVAC
  ! parameters space

  integer :: iresult, i
  real*8,parameter :: mh=1.67262171d-24, kboltz=1.3806504d-16, fH=0.76
  real(kind=gr_rpknd) :: tiny_number= 1.0e-20_gr_RKIND,huge_number= 1.0e+20_gr_RKIND


  !    Initialization parameters

  real*8           :: initial_redshift
  character(len=80), TARGET :: filename

  !     Field data arrays

  integer          :: field_size

  real*8                     :: temperature_units, pressure_units, dtchem

  ! real(8), allocatable, target

  type amrvac_grackle_arrays

    !     Field data arrays

    integer, allocatable ::                  :: field_size(:)

    real*8                                   :: temperature_units, pressure_units, dtchem

    real(kind=gr_rpknd), allocatable :: density(:)
    real(kind=gr_rpknd), allocatable :: energy(:)
    real(kind=gr_rpknd), allocatable :: x_velocity(:)
    real(kind=gr_rpknd), allocatable :: y_velocity(:)
    real(kind=gr_rpknd), allocatable :: z_velocity(:)
    real(kind=gr_rpknd), allocatable :: HI_density(:)
    real(kind=gr_rpknd), allocatable :: HII_density(:)
    real(kind=gr_rpknd), allocatable :: HM_density(:)
    real(kind=gr_rpknd), allocatable :: HeI_density(:)
    real(kind=gr_rpknd), allocatable :: HeII_density(:)
    real(kind=gr_rpknd), allocatable :: HeIII_density(:)
    real(kind=gr_rpknd), allocatable :: H2I_density(:)
    real(kind=gr_rpknd), allocatable :: H2II_density(:)
    real(kind=gr_rpknd), allocatable :: DI_density(:)
    real(kind=gr_rpknd), allocatable :: DII_density(:)
    real(kind=gr_rpknd), allocatable :: HDI_density(:)
    real(kind=gr_rpknd), allocatable :: e_density(:)
    real(kind=gr_rpknd), allocatable :: metal_density(:)
    real(kind=gr_rpknd), allocatable :: dust_density(:)
    real(kind=gr_rpknd), allocatable :: volumetric_heating_rate(:)
    real(kind=gr_rpknd), allocatable :: specific_heating_rate(:)
    real(kind=gr_rpknd), allocatable :: RT_HI_ionization_rate(:)
    real(kind=gr_rpknd), allocatable :: RT_HeI_ionization_rate(:)
    real(kind=gr_rpknd), allocatable :: RT_HeII_ionization_rate(:)
    real(kind=gr_rpknd), allocatable :: RT_H2_dissociation_rate(:)
    real(kind=gr_rpknd), allocatable :: RT_heating_rate(:)

    !output parameters
    real(kind=gr_rpknd), allocatable :: cooling_time(:)
    real(kind=gr_rpknd), allocatable :: gamma(:)
    real(kind=gr_rpknd), allocatable :: pressure(:)
    real(kind=gr_rpknd), allocatable :: temperature(:)
    real(kind=gr_rpknd), allocatable :: dust_temperature(:)

    !     Grid size and dimension
    !     grid_start and grid_end are used to ignore ghost zones.
    !     grid_dx is used in H2 self-shielding approximation only

    integer             :: grid_rank
    integer             :: grid_dimension(ndim)
    integer             :: grid_start(ndim)
    integer             :: grid_end(ndim)
    real*8                           :: grid_dx

  end type amrvac_grackle_arrays

  real(kind=gr_rpknd), allocatable, target :: density(:), energy(:), &
      x_velocity(:), y_velocity(:), &
      z_velocity(:), &
      HI_density(:), HII_density(:), &
      HM_density(:), &
      HeI_density(:), HeII_density(:), &
      HeIII_density(:), &
      H2I_density(:), H2II_density(:), &
      DI_density(:), DII_density(:), &
      HDI_density(:), &
      e_density(:), metal_density(:), &
      dust_density(:), &
      volumetric_heating_rate(:), &
      specific_heating_rate(:), &
      RT_HI_ionization_rate(:), &
      RT_HeI_ionization_rate(:), &
      RT_HeII_ionization_rate(:), &
      RT_H2_dissociation_rate(:), &
      RT_heating_rate(:)


  real(kind=gr_rpknd), allocatable, target :: cooling_time(:), gamma(:), &
      pressure(:), temperature(:), &
      dust_temperature(:)

  !     Grid size and dimension
  !     grid_start and grid_end are used to ignore ghost zones.
  !     grid_dx is used in H2 self-shielding approximation only

  INTEGER, allocatable, target :: grid_rank, grid_dimension(ndim), &
      grid_start(ndim), grid_end(ndim)

  real*8 :: grid_dx



contains

subroutine check_indexes_inconsistency(ixI^L,ixO^L)
  implicit none
  integer, intent(in)                     :: ixI^L,ixO^L



  {^D&
  if(ixImax^D<ixImin^D)

    write(*,*) 'We get ixImin^D=',ixImin^D
    write(*,*) 'We get ixImax^D=',ixImax^D
    call mpistop('Inconsistency with ixImax^D<ixImin^D')

  end if

  if(ixOmax^D<ixOmin^D)

    write(*,*) 'We get ixOmin^D=',ixOmin^D
    write(*,*) 'We get ixOmax^D=',ixOmax^D
    call mpistop('Inconsistency with ixOmax^D<ixOmin^D')
  end if
  }


end check_indexes_inconsistency

subroutine deallocate_grackle_grids(ixI^L,ixO^L,qt,x,self)
  implicit none
  integer, intent(in)                     :: ixI^L,ixO^L
  real(kind=dp), intent(in)               :: qt
  real(kind=dp), intent(in)               :: x(ixI^S,1:ndim)
  class(amrvac_grackle_arrays)            :: self

  call check_indexes_inconsistency(ixI^L,ixO^L)

  if(allocated(self%field_size)) deallocate(self%field_size)
  if(allocated(self%density)) deallocate(self%density)
  if(allocated(self%energy)) deallocate(self%energy)
  if(allocated(self%x_velocity)) deallocate(self%x_velocity)
  if(allocated(self%y_velocity)) deallocate(self%y_velocity)
  if(allocated(self%z_velocity)) deallocate(self%z_velocity)
  if(allocated(self%HI_density)) deallocate(self%HI_density)
  if(allocated(self%HII_density)) deallocate(self%HII_density)
  if(allocated(self%HM_density)) deallocate(self%HM_density)
  if(allocated(self%HeI_density)) deallocate(self%HeI_density)
  if(allocated(self%HeII_density)) deallocate(self%HeII_density)
  if(allocated(self%HeIII_density)) deallocate(self%HeIII_density)
  if(allocated(self%H2I_density)) deallocate(self%H2I_density)
  if(allocated(self%H2II_density)) deallocate(self%H2II_density)
  if(allocated(self%DI_density)) deallocate(self%DI_density)
  if(allocated(self%DII_density)) deallocate(self%DII_density)
  if(allocated(self%HDI_density)) deallocate(self%HDI_density)
  if(allocated(self%e_density)) deallocate(self%e_density)
  if(allocated(self%metal_density)) deallocate(self%metal_density)
  if(allocated(self%dust_density)) deallocate(self%dust_density)
  if(allocated(self%volumetric_heating_rate)) deallocate(self%volumetric_heating_rate)
  if(allocated(self%specific_heating_rate)) deallocate(self%specific_heating_rate)
  if(allocated(self%RT_HI_ionization_rate)) deallocate(self%RT_HI_ionization_rate)
  if(allocated(self%RT_HeI_ionization_rate)) deallocate(self%RT_HeI_ionization_rate)
  if(allocated(self%RT_HeII_ionization_rate)) deallocate(self%RT_HeII_ionization_rate)
  if(allocated(self%RT_H2_dissociation_rate)) deallocate(self%RT_H2_dissociation_rate)
  if(allocated(self%RT_heating_rate)) deallocate(self%RT_heating_rate)
  if(allocated(self%cooling_time)) deallocate(self%cooling_time)
  if(allocated(self%gamma)) deallocate(self%gamma)
  if(allocated(self%pressure)) deallocate(self%pressure)
  if(allocated(self%temperature)) deallocate(self%temperature)
  if(allocated(self%dust_temperature)) deallocate(self%dust_temperature)

end subroutine deallocate_grackle_grids

subroutine allocate_grackle_grids(ixI^L,ixO^L,qt,x,self)
  implicit none
  integer, intent(in)                     :: ixI^L,ixO^L
  real(kind=dp), intent(in)               :: qt
  real(kind=dp), intent(in)               :: x(ixI^S,1:ndim)
  class(amrvac_grackle_arrays)            :: self

  call check_indexes_inconsistency(ixI^L,ixO^L)
  call self%deallocate_grackle_grids(ixI^L,ixO^L,qt,x)

  !ultimately, field size must be computed based on ixI^L:
  ! 1D: = (ixImin1),...,(ixI^max1)
  ! 2D: = (ixImin1,ixImin2),(ixImin1,ixImin2+1),...,(ixImax1,ixImax2)
  ! 3D: = (ixImin1,ixImin2,ixImin3),(ixImin1,ixImin2,ixImin3+1),...,(ixImax1,ixImax2,ixImax3)

  ! Field data arrays

  allocate(self%field_size(ndim))
  allocate(self%density(ixI^S))
  allocate(self%energy(ixI^S))
  allocate(self%x_velocity(ixI^S))
  allocate(self%y_velocity(ixI^S))
  allocate(self%z_velocity(ixI^S))
  allocate(self%HI_density(ixI^S))
  allocate(self%HII_density(ixI^S))
  allocate(self%HM_density(ixI^S))
  allocate(self%HeI_density(ixI^S))
  allocate(self%HeII_density(ixI^S))
  allocate(self%HeIII_density(ixI^S))
  allocate(self%H2I_density(ixI^S))
  allocate(self%H2II_density(ixI^S))
  allocate(self%DI_density(ixI^S))
  allocate(self%DII_density(ixI^S))
  allocate(self%HDI_density(ixI^S))
  allocate(self%e_density(ixI^S))
  allocate(self%metal_density(ixI^S))
  allocate(self%dust_density(ixI^S))
  allocate(self%volumetric_heating_rate(ixI^S))
  allocate(self%specific_heating_rate(ixI^S))
  allocate(self%RT_HI_ionization_rate(ixI^S))
  allocate(self%RT_HeI_ionization_rate(ixI^S))
  allocate(self%RT_HeII_ionization_rate(ixI^S))
  allocate(self%RT_H2_dissociation_rate(ixI^S))
  allocate(self%RT_heating_rate(ixI^S))
  allocate(self%cooling_time(ixI^S))
  allocate(self%gamma(ixI^S))
  allocate(self%pressure(ixI^S))
  allocate(self%temperature(ixI^S))
  allocate(self%dust_temperature(ixI^S))


end allocate_grackle_grids

! Arrays total size
self%field_size = {(ixImax^D-ixImin^D+1)|*}

!     Grid size and dimension
!     grid_start and grid_end are used to ignore ghost zones.
!     grid_dx is used in H2 self-shielding approximation only

integer             :: grid_rank
integer             :: grid_dimension(ndim)
integer             :: grid_start(ndim)
integer             :: grid_end(ndim)
real*8                           :: grid_dx

subroutine grackle_fields_set_w(ixI^L,ixO^L,qt,x,self)
  implicit none
  integer, intent(in)                     :: ixI^L,ixO^L
  real(kind=dp), intent(in)               :: qt
  real(kind=dp), intent(in)               :: x(ixI^S,1:ndim)
  class(amrvac_grackle_arrays)            :: self
  !-----local----------------------------------
  integer                                 :: ifield

  call check_indexes_inconsistency(ixI^L,ixO^L)

  ! Test case:

  !     Create a grackle chemistry object for parameters and set defaults

  iresult = set_default_chemistry_parameters(grackle_data)

  !     Set parameters

  grackle_data%use_grackle = 1            ! chemistry on

  PRINT*, 'on est arrive jusqu ici 1 '

  !     Set parameters

  grackle_data%use_grackle = 1            ! chemistry on
  grackle_data%with_radiative_cooling = 1 ! cooling on
  grackle_data%primordial_chemistry = 3   ! network with H, He, D
  grackle_data%dust_chemistry = 1         ! dust processes
  grackle_data%metal_cooling = 1          ! metal cooling on
  grackle_data%UVbackground = 1           ! UV background on
  !     cooling data for Haardt & Madau 2012 background
  filename = "input/CloudyData_UVB=HM2012.h5"//C_NULL_CHAR
  grackle_data%grackle_data_file = C_LOC(filename(1:1))
  grackle_data%h2_on_dust = 0             ! no dust
  grackle_data%cmb_temperature_floor = 1  ! include CMB cooling floor
  grackle_data%Gamma = 5./3.;          ! monoatomic gas

  !     Set units

  my_units%comoving_coordinates = 0
  my_units%density_units = 1.67d-24
  my_units%length_units = 1.0d0
  my_units%time_units = 1.0d12
  my_units%a_units = 1.0d0

  !     Set initial expansion factor (for internal units).
  !     Set expansion factor to 1 for non-cosmological simulation.
  initial_redshift = 0.;
  my_units%a_value = 1. / (1. + initial_redshift);

  PRINT*, 'on est arrive jusqu ici 2 '

  call set_velocity_units(my_units)

  !     Initialize the Grackle

  write(6,*) "primordial_chemistry:", &
      grackle_data%primordial_chemistry
  write(6,*) "metal_cooling:", &
      grackle_data%metal_cooling
  iresult = initialize_chemistry_data(my_units)

  !     Set field arrays

  !     If grid rank is less than 3, set the other dimensions,
  !     start indices, and end indices to 0.
  self%grid_rank = 3
  do i = 1, grid_rank
     self%grid_dimension(i) = 1
     self%grid_start(i) = 0
     self%grid_end(i) = 0
  enddo
  grid_dx = 0.0
  grid_dimension(1) = self%field_size
  !     0-based
  grid_end(1) = self%field_size - 1

  temperature_units = get_temperature_units(my_units)

  do ifield=1,ndim
    self%field_size(ndim)=1 !one cell only
  end do

  self%density(ixO^S)= 1.0
  ! initilize internal energy (here 1000 K for no reason)
  self%energy(ixO^S)= 1000. / temperature_units
  self%x_velocity(ixO^S)= 0.0
  self%y_velocity(ixO^S)= 0.0
  self%z_velocity(ixO^S)= 0.0
  self%HI_density(ixO^S)= fH * self%density(ixO^S)
  self%HII_density(ixO^S)= tiny_number * self%density(ixO^S)
  self%HM_density(ixO^S)= tiny_number * self%density(ixO^S)
  self%HeI_density(ixO^S)= (1.0 - fH) * self%density(ixO^S)
  self%HeII_density(ixO^S)= tiny_number * self%density(ixO^S)
  self%HeIII_density(ixO^S)= tiny_number * self%density(ixO^S)
  self%H2I_density(ixO^S)= tiny_number * self%density(ixO^S)
  self%H2II_density(ixO^S)= tiny_number * self%density(ixO^S)
  self%DI_density(ixO^S)= 2.0 * 3.4e-5 * self%density(ixO^S)
  self%DII_density(ixO^S)= tiny_number * self%density(ixO^S)
  self%HDI_density(ixO^S)= tiny_number * self%density(ixO^S)
  self%e_density(ixO^S)= tiny_number * self%density(ixO^S)
  !solar metallicity
  self%metal_density(ixO^S)= grackle_data%SolarMetalFractionByMass * &
                             self%density(ixO^S)
  self%dust_density(ixO^S)= grackle_data%local_dust_to_gas_ratio * &
                            self%density(ixO^S)
  self%volumetric_heating_rate(ixO^S)= 0.0
  self%specific_heating_rate(ixO^S)= 0.0
  self%RT_HI_ionization_rate(ixO^S)= 0.0
  self%RT_HeI_ionization_rate(ixO^S)= 0.0
  self%RT_HeII_ionization_rate(ixO^S)= 0.0
  self%RT_H2_dissociation_rate(ixO^S)= 0.0
  self%RT_heating_rate(ixO^S)= 0.0

end subroutine grackle_fields_set_w

subroutine test_chemistry()

  PRINT *, 'Linking with obj_global_parameters'
  PRINT *, 'ndim='
  PRINT *, ndim
  PRINT *, 'mod_global_parameters:w_refine_weight='
  PRINT *, w_refine_weight
  PRINT *, 'mod_global_parameters:clight='
  PRINT *, constusr%clight

  call test_grackle_header()

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  !     Create a grackle chemistry object for parameters and set defaults

        iresult = set_default_chemistry_parameters(grackle_data)

  !     Set parameters

        grackle_data%use_grackle = 1            ! chemistry on

        PRINT*, 'on est arrive jusqu ici 1 '

  !     Set parameters

        grackle_data%use_grackle = 1            ! chemistry on
        grackle_data%with_radiative_cooling = 1 ! cooling on
        grackle_data%primordial_chemistry = 3   ! network with H, He, D
        grackle_data%dust_chemistry = 1         ! dust processes
        grackle_data%metal_cooling = 1          ! metal cooling on
        grackle_data%UVbackground = 1           ! UV background on
  !     cooling data for Haardt & Madau 2012 background
        filename = "input/CloudyData_UVB=HM2012.h5"//C_NULL_CHAR
        grackle_data%grackle_data_file = C_LOC(filename(1:1))
        grackle_data%h2_on_dust = 0             ! no dust
        grackle_data%cmb_temperature_floor = 1  ! include CMB cooling floor
        grackle_data%Gamma = 5./3.;          ! monoatomic gas

  !     Set units

        my_units%comoving_coordinates = 0
        my_units%density_units = 1.67d-24
        my_units%length_units = 1.0d0
        my_units%time_units = 1.0d12
        my_units%a_units = 1.0d0

  !     Set initial expansion factor (for internal units).
  !     Set expansion factor to 1 for non-cosmological simulation.
        initial_redshift = 0.;
        my_units%a_value = 1. / (1. + initial_redshift);

        PRINT*, 'on est arrive jusqu ici 2 '

        call set_velocity_units(my_units)

  !     Initialize the Grackle

        write(6,*) "primordial_chemistry:", &
            grackle_data%primordial_chemistry
        write(6,*) "metal_cooling:", &
            grackle_data%metal_cooling
        iresult = initialize_chemistry_data(my_units)

  !     Set field arrays

  !     If grid rank is less than 3, set the other dimensions,
  !     start indices, and end indices to 0.
        grid_rank = 3
        do i = 1, grid_rank
           grid_dimension(i) = 1
           grid_start(i) = 0
           grid_end(i) = 0
        enddo
        grid_dx = 0.0
        grid_dimension(1) = field_size
  !     0-based
        grid_end(1) = field_size - 1

        temperature_units = get_temperature_units(my_units)

        do i = 1,field_size
           density(i) = 1.0
           HI_density(i) = fH * density(i)
           HII_density(i) = tiny_number * density(i)
           HM_density(i) = tiny_number * density(i)
           HeI_density(i) = (1.0 - fH) * density(i)
           HeII_density(i) = tiny_number * density(i)
           HeIII_density(i) = tiny_number * density(i)
           H2I_density(i) = tiny_number * density(i)
           H2II_density(i) = tiny_number * density(i)
           DI_density(i) = 2.0 * 3.4e-5 * density(i)
           DII_density(i) = tiny_number * density(i)
           HDI_density(i) = tiny_number * density(i)
           e_density(i) = tiny_number * density(i)
  !        solar metallicity
           metal_density(i) = grackle_data%SolarMetalFractionByMass * &
               density(i)
           dust_density(i) = grackle_data%local_dust_to_gas_ratio * &
               density(i)

           x_velocity(i) = 0.0
           y_velocity(i) = 0.0
           z_velocity(i) = 0.0

  !        initilize internal energy (here 1000 K for no reason)
           energy(i) = 1000. / temperature_units

           volumetric_heating_rate(i) = 0.0
           specific_heating_rate(i) = 0.0
           RT_HI_ionization_rate(i) = 0.0
           RT_HeI_ionization_rate(i) = 0.0
           RT_HeII_ionization_rate(i) = 0.0
           RT_H2_dissociation_rate(i) = 0.0
           RT_heating_rate(i) = 0.0
        enddo
  !
  !     Fill in structure to be passed to Grackle
  !
        my_fields%grid_rank = 1
        my_fields%grid_dimension = C_LOC(grid_dimension)
        my_fields%grid_start = C_LOC(grid_start)
        my_fields%grid_end = C_LOC(grid_end)
        my_fields%grid_dx  = grid_dx

        my_fields%density = C_LOC(density)
        my_fields%HI_density = C_LOC(HI_density)
        my_fields%HII_density = C_LOC(HII_density)
        my_fields%HM_density = C_LOC(HM_density)
        my_fields%HeI_density = C_LOC(HeI_density)
        my_fields%HeII_density = C_LOC(HeII_density)
        my_fields%HeIII_density = C_LOC(HeIII_density)
        my_fields%H2I_density = C_LOC(H2I_density)
        my_fields%H2II_density = C_LOC(H2II_density)
        my_fields%DI_density = C_LOC(DI_density)
        my_fields%DII_density = C_LOC(DII_density)
        my_fields%HDI_density = C_LOC(HDI_density)
        my_fields%e_density = C_LOC(e_density)
        my_fields%metal_density = C_LOC(metal_density)
        my_fields%internal_energy = C_LOC(energy)
        my_fields%x_velocity = C_LOC(x_velocity)
        my_fields%y_velocity = C_LOC(y_velocity)
        my_fields%z_velocity = C_LOC(z_velocity)
        my_fields%volumetric_heating_rate = &
                                       C_LOC(volumetric_heating_rate)
        my_fields%specific_heating_rate = C_LOC(specific_heating_rate)
        my_fields%RT_HI_ionization_rate = C_LOC(RT_HI_ionization_rate)
        my_fields%RT_HeI_ionization_rate = C_LOC(RT_HeI_ionization_rate)
        my_fields%RT_HeII_ionization_rate = &
                                         C_LOC(RT_HeII_ionization_rate)
        my_fields%RT_H2_dissociation_rate = C_LOC(RT_H2_dissociation_rate)
        my_fields%RT_heating_rate = C_LOC(RT_heating_rate)

  !
  !     Calling the chemistry solver
  !     These routines can now be called during the simulation.

  !     Evolving the chemistry.

        dtchem = 3.15e7 * 1e6 / my_units%time_units    ! some timestep

        !write(*,*) 'HI density before (local,global):',my_fields%density, density
        iresult = solve_chemistry(my_units, my_fields, dtchem)
        !write(*,*) 'HI density after (local,global):' ,my_fields%density, density

  !     Calculate cooling time.

        iresult = calculate_cooling_time(my_units, my_fields,cooling_time)
        write(*,*) "Cooling time = ", (cooling_time(1) * &
            my_units%time_units), "s."

  !     Calculate temperature.

        iresult = calculate_temperature(my_units, my_fields, temperature)
        write(*,*) "Temperature = ", temperature(1), "K."

  !     Calcualte pressure.

        pressure_units = my_units%density_units * &
            my_units%velocity_units**2
        iresult = calculate_pressure(my_units, my_fields, pressure)
        write(*,*) "Pressure = ", pressure(1)*pressure_units, "dyne/cm^2."

  !     Calculate gamma.

        iresult = calculate_gamma(my_units, my_fields, gamma)
        write(*,*) "Gamma = ", gamma(1)

  !     Calculate dust temperature.

        iresult = calculate_dust_temperature(my_units, my_fields,&
            dust_temperature)
        write(*,*) "Dust temperature = ", dust_temperature(1), "K."

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end subroutine test_chemistry

end module mod_grackle_parameters
