module mod_grackle_chemistry
  use mod_constants
  use mod_global_parameters
  use mod_physics
  use mod_obj_dust, only : dust
  use mod_obj_global_parameters
  use mod_obj_usr_unit
  use mod_obj_mat
  use grackle_header
  use mod_grackle_parameters

character(len=*),parameter :: mod_grackle_chemistry_name='mod_grackle_chemistry'


! Grackle AMRVAC-defined parameters for this solver time step
type gr_config
  character(len=78)    :: obj_name       !> Obj name that call it
  character(len=20)    :: unit           !> physical unit at parameter file
  integer              :: myindice       !> ism indices associated with ism in use

  ! Parameters as taken from RAMSES implimentation with Grackle
  ! Chemistry with Grackle : on=1 ; off=0
  INTEGER :: use_grackle

  ! Cooling with Grackle : on=1 ; off=0
  INTEGER :: gr_with_radiative_cooling

  ! Molecular/atomic network solved by Grackle :
  ! 0: no chemistry network. Radiative cooling for primordial species is solved by interpolating from lookup tables calculated with Cloudy.
  ! 1: 6-species atomic H and He.
  ! Active species: H, H+, He, He+, ++, e-.
  ! 2: 9-species network including atomic species above and species
  ! for molecular hydrogen formation.
  ! This network includes formation from the H- and H2+ channels,
  ! three-body formation (H+H+H and H+H+H2),
  ! H2 rotational transitions, chemical heating,
  ! and collision-induced emission (optional).
  ! Active species: above + H-, H2, H2+.
  ! 3: 12-species network include all above plus HD rotation cooling.
  ! Active species: above + D, D+, HD.
  ! Reactions listed in Table 3 of Smith et al. 2017
  INTEGER :: gr_primordial_chemistry

  ! Include metal cooling (Smith et al. 2017) :
  ! on = 1 ; off = 0
  ! If enabled, the cooling table to be
  ! used must be specified with the grackle_data_file
  INTEGER :: gr_metal_cooling

  ! Include UV background :
  ! on = 1 ; off = 0
  ! If enabled, the cooling table
  ! to be used must be specified with the grackle_data_file parameter
  INTEGER :: gr_UVbackground

  ! Flag to enable an effective CMB temperature floor.
  ! on = 1 ; off = 0
  ! This is implemented by subtracting the value
  ! of the cooling rate at TCMB from the total cooling rate.
  ! If enabled, the cooling table
  ! to be used must be specified with the grackle_data_file parameter
  INTEGER :: gr_cmb_temperature_floor

  ! Flag to enable H2 formation on dust grains,
  ! dust cooling, and dust-gas heat transfer
  ! follow Omukai (2000). This assumes that the dust
  ! to gas ratio scales with the metallicity.
  INTEGER :: gr_h2_on_dust


  ! Flag to control additional dust cooling
  ! and chemistry processes.
  ! 0: no dust-related processes included.
  ! 1: adds the following processes:
  !   1. photo-electric heating (sets photoelectric_heating to 2).
  !   2. cooling from electron recombination onto dust (equation 9
  !      from Wolfire et al. 1995). Both the photo-electric heating
  !      and recombination cooling are scaled by the value of the
  !      interstellar_radiation_field.
  !   3. H2 formation on dust (sets h2_on_dust to 1 if primordial_chemistry > 1).
  ! Setting dust_chemistry greater than 0 requires metal_cooling to be enabled.
  INTEGER :: gr_dust_chemistry

  ! Flag to provide the dust density as a field using the dust_density
  ! pointer in the grackle_field_data struct. If set to 0,
  ! the dust density takes the value of local_dust_to_gas_ratio
  ! multiplied by the metallicity
  INTEGER :: gr_use_dust_density_field

  ! Flag to enable photo-electric heating from irradiated dust grains. Default: 0.
  ! 0: no photo-electric heating.
  ! 1: a spatially uniform heating term
  ! from Tasker & Bryan (2008). The exact
  ! heating rate used must be specified
  ! with the photoelectric_heating_rate
  ! parameter. For temperatures above
  ! 20,000 K, the photo-electric heating
  ! rate is set to 0.
  ! 2: similar to option 1, except the
  ! heating rate is calculated using
  ! equation 1 of Wolfire et al. (1995)
  ! and the user must supply the intensity
  ! of the interstellar radiation field with the
  ! interstellar_radiation_field parameter.
  ! The value of epsilon is taken as a constant equal to 0.05
  ! for gas below 20,000 K and 0 otherwise.
  ! 3: similar to option 1, except the value of
  ! epsilon is calculated directly
  ! from equation 2 of Wolfire et al. (1995).
  INTEGER :: gr_photoelectric_heating

  ! Flag to signal that an array of volumetric
  ! heating rates is being provided in the
  ! volumetric_heating_rate field of the
  ! grackle_field_data struct.
  ! on = 1 ; off = 0
  INTEGER :: gr_use_volumetric_heating_rate

  ! Flag to signal that an array of
  ! specific heating rates is being
  ! provided in the specific_heating_rate
  ! field of the grackle_field_data struct.
  ! on = 1 ; off = 0
  INTEGER :: gr_use_specific_heating_rate

  ! Flag to control which three-body H2 formation rate is used.
  !0: Abel, Bryan & Norman (2002)
  !1: Palla, Salpeter & Stahler (1983)
  !2: Cohen & Westberg (1983)
  !3: Flower & Harris (2007)
  !4: Glover (2008)
  !5: Forrey (2013).
  !The first five options are discussed in Turk et. al. (2011).
  INTEGER :: gr_three_body_rate

  ! Flag to enable H2 collision-induced
  ! emission cooling from Ripamonti & Abel (2004).
  ! on = 1 ; off = 0
  INTEGER :: gr_cie_cooling

  ! Flag to enable H2 cooling attenuation
  ! from Ripamonti & Abel (2004).
  ! on = 1 ; off = 0
  INTEGER :: gr_h2_optical_depth_approximation = 0

  INTEGER :: gr_ih2co
  INTEGER :: gr_ipiht
  INTEGER :: gr_NumberOfTemperatureBins

  ! 0 : The recombination of H + , He + and He ++
  ! is modelled using the case A recombination rate coefficients
  ! (the optically-thin approximation in which recombination
  ! photons above 1 Ryd escape).
  ! 1 : case B rate coefficients (in which recombination photons above 1 Ryd
  ! are locally re-absorbed, Osterbrock 1989) can instead be se-
  ! lected by setting CaseBRecombination = 1.
  INTEGER :: gr_CaseBRecombination

  ! Flag to enable Compton heating
  ! from an X-ray background following
  ! Madau & Efstathiou (1999).
  ! on = 1 ; off = 0
  INTEGER :: gr_Compton_xray_heating

  ! Flag to enable suppression of Lyman-Werner flux due to Lyman-series
  ! absorption (giving a sawtooth pattern), taken from Haiman & Abel, & Rees (2000)
  ! on = 1 ; off = 0
  INTEGER :: gr_LWbackground_sawtooth_suppression

  INTEGER :: gr_NumberOfDustTemperatureBins


  ! Flag to signal that arrays of ionization and heating rates from
  ! radiative transfer solutions are being provided.
  ! Only available if primordial_chemistry is greater than 0.
  ! HI, HeI, and HeII ionization arrays are provided in
  ! RT_HI_ionization_rate, RT_HeI_ionization_rate, and RT_HeII_ionization_rate fields, respectively,
  ! of the grackle_field_data struct. Associated heating rate is
  ! provided in the RT_heating_rate field, and H2photodissociation
  ! rate can also be provided in the RT_H2_dissociation_rate field
  ! when primordial_chemistry is set to either 2 or 3.
  ! on = 1 ; off = 0
  INTEGER :: gr_use_radiative_transfer


  ! When used with use_radiative_transfer set to 1,
  ! this flag makes it possible to solve the chemistry and cooling
  ! of the computational elements for which the radiation field is non-zero
  ! separately from those with no incident radiation. This allows radiation transfer calculations to
  ! be performed on a smaller timestep than the global timestep. The parameter, radiative_transfer_intermediate_step,
  ! is then used to toggle between updating the cells/particles receiving radiative input and those that
  ! are not.
  ! on = 1 ; off = 0
  INTEGER :: gr_radiative_transfer_coupled_rate_solver

  ! Used in conjunction with radiative_transfer_coupled_rate_solver set to 1, setting this parameter to 1 tells the solver
  ! to only update cells/particles where the radiation field is non-zero. Setting this to 0 updates only those elements with
  ! no incident radiation. When radiative_transfer_coupled_rate_solver is set to 0, changing this parameter
  ! will have no effect.
  INTEGER :: gr_radiative_transfer_intermediate_step

  ! Flag to only use hydrogen ionization and heating rates from the radiative transfer solutions.
  ! on = 1 ; off = 0
  INTEGER :: gr_radiative_transfer_hydrogen_only

  ! Switch to enable approximate self-shielding from the UV background. All three of the below methods incorporate Eq. 13 and
  ! 14 from Rahmati et. al. 2013. These equations involve using the spectrum averaged photoabsorption cross for the given species (HI or HeI).
  ! These redshift dependent values are pre-computed for the HM2012 and FG2011 UV backgrounds and included in their respective cooling data tables.
  ! Care is advised in using any of these methods. The default behavior is to apply no self-shielding, but this is not necessarily the proper assumption,
  ! depending on the use case. If the user desires to turn on self-shielding, we strongly advise using option 3. All options include HI self-shielding,
  ! and vary only in treatment of HeI and HeII. In options 2 and 3, we approximately account for HeI self-shielding by applying the Rahmati et. al. 2013 relations,
  ! which are only strictly valid for HI, !to HeI under the assumption that it behaves similarly to HI. None of these options are completely correct in practice,
  ! but option 3 has produced the most reasonable results in test simulations. Repeating the analysis of Rahmati et. al. 2013 to directly parameterize HeI and HeII
  ! self-shielding behavior would be a valuable avenue of future research in developing a more complete self-shielding model. Each self-shielding option is described below.
  ! 0: No self shielding. Elements are optically thin to the UV background.
  ! 1: Not Recommended. Approximate self-shielding in HI only.
  ! HeI and HeII are left as optically thin.
  ! 2: Approximate self-shielding in both HI and HeI. HeII remains
  ! optically thin.
  ! 3: Approximate self-shielding in both HI and HeI, but ignoring
  ! HeII ionization and heating from the UV background entirely (HeII ionization and heating rates are set to zero).
  ! These methods only work in conjunction with using updated Cloudy cooling tables, denoted with “_shielding”. These tables properly account for the decrease
  ! in metal line cooling rates in self-shielded regions, which can be significant.
  !For consistency, when primordial_chemistry > 2, the self-shielding attenutation factors calculated for HI and HeI are applied to the H2ionization (15.4 eV) and H2+
  ! dissociation rates (30 eV) respectively. These reaction rates are distinct from the H2self-shielding computed using the H2_self_shielding flag.
  ! on = 1 ; off = 0
  INTEGER :: gr_self_shielding_method

  ! The ratio of specific heats for an ideal gas.
  ! A direct calculation for the molecular component
  ! is used if primordial_chemistry > 1.
  REAL(kind=gr_rpknd)    :: gr_Gamma

  ! If photoelectric_heating is enabled, the heating rate in units of (erg cm-3/s)
  ! n-1, where n is the total hydrogen number density. In other words, this is the
  ! volumetric heating rate at a hydrogen number density of n = 1 cm-3.
  REAL(kind=gr_rpknd)    :: gr_photoelectric_heating_rate

  ! The fraction by mass of Hydrogen in the metal-free portion of the gas (i.e., just the H and He).
  ! In the non-equilibrium solver, this is used to ensure consistency in the densities of the individual species.
  ! In tabulated mode, this is used to calculate the H number density from the total gas density,
  ! which is a parameter of the heating/cooling tables. When using the non-equilibrium solver,
  ! a sensible default is 0.76. However, the tables for tabulated mode were created assuming nHe/nH = 0.1,
  ! which corresponds to an H mass fraction of about 0.716.
  ! When running in tabulated mode, this parameter will automatically be changed to this value.
  REAL(kind=gr_rpknd)    :: gr_HydrogenFractionByMass

  ! The ratio by mass of Deuterium to Hydrogen.
  ! Default: 6.8e-5 (the value from Burles & Tytler (1998)
  ! multiplied by 2 for the mass of Deuterium).
  REAL(kind=gr_rpknd)    :: gr_DeuteriumToHydrogenRatio

  ! The fraction of total gas mass in metals for a solar composition.
  ! Default: 0.01295 (consistent with the default abundances in the Cloudy code).
  REAL(kind=gr_rpknd)    :: gr_SolarMetalFractionByMass

  ! Temperature limits
  REAL(kind=gr_rpknd)    :: TemperatureStart
  REAL(kind=gr_rpknd)    :: TemperatureEnd

  ! Dust temperature limits
  REAL(kind=gr_rpknd)    :: DustTemperatureStart
  REAL(kind=gr_rpknd)    :: DustTemperatureEnd

  ! Intensity of a constant Lyman-Werner H2 photo-dissociating
  ! radiation field in units of 10-21 erg /s /cm2 Hz-1 sr-1. Default: 0.
  REAL(kind=gr_rpknd)    :: gr_LWbackground_intensity

  ! Used in combination with UVbackground_redshift_fullon, UVbackground_redshift_drop,
  ! and UVbackground_redshift_off to set an attenuation factor
  ! for the photo-heating and photo-ionization rates of the UV background model.
  ! See the figure below for an illustration its behavior. If not set, this parameter will
  ! be set to the highest redshift of the UV background data being used.
  REAL(kind=gr_rpknd)    :: gr_UVbackground_redshift_on

  ! Used in combination with UVbackground_redshift_on, UVbackground_redshift_fullon, and
  ! UVbackground_redshift_drop to set an attenuation factor for the photo-heating and
  ! photo-ionization rates of the UV background model. See the figure below for an illustration its behavior.
  ! If not set, this parameter will be set to the lowest redshift of the UV background data being used.
  REAL(kind=gr_rpknd)    :: gr_UVbackground_redshift_off

  ! Used in combination with UVbackground_redshift_on, UVbackground_redshift_drop, and UVbackground_redshift_off
  ! to set an attenuation factor for the photo-heating and photo-ionization rates of the UV background model.
  ! See the figure below for an illustration its behavior. If not set, this parameter will be set to the highest
  ! redshift of the UV background data being used.
  REAL(kind=gr_rpknd)    :: gr_UVbackground_redshift_fullon

  ! Used in combination with UVbackground_redshift_on, UVbackground_redshift_fullon, and UVbackground_redshift_off
  ! to set an attenuation factor for the photo-heating and photo-ionization rates of the UV background model.
  ! See the figure below for an illustration its behavior. If not set, this parameter will be set to the lowest
  ! redshift of the UV background data being used.
  REAL(kind=gr_rpknd)    :: gr_UVbackground_redshift_drop

  ! Cloudy 07.02 abundances :
  ! A float value to account for additional electrons contributed by metals. This is only used with Cloudy datasets
  ! with dimension greater than or equal to 4. The value of this factor is calculated as the sum of (Ai * i) over all
  ! elements i heavier than He, where Ai is the solar number abundance relative to H. For the solar abundance pattern
  ! from the latest version of Cloudy, using all metals through Zn, this value is 9.153959e-3. Default: 9.153959e-3.
  REAL(kind=gr_rpknd)    :: cloudy_electron_fraction_factor
  CHARACTER(LEN=128) :: data_filename
  CHARACTER(LEN=128) :: data_dir

  logical :: normalize_done
  INTEGER :: gr_comoving_coordinates
  REAL(kind=gr_rpknd)    :: gr_a_units
  REAL(kind=gr_rpknd)    :: gr_current_redshift
end type gr_config

integer,allocatable              :: tstst(:)

! Fields parameters
! parts to be initialized with chemical species
! 'ism' : ISM
! 'jet' : jet (on top of the ISM)
integer :: i_zero_ism,i_zero_jet,i_zero_cloud
integer :: i_end_ism,i_end_jet,i_end_cloud

!for namelist:
CHARACTER(LEN=30), allocatable :: gr_patches_name(:)

! indices of the objects to be initialized with chemical species
INTEGER, allocatable :: gr_patches_indices_global(:)
INTEGER, allocatable :: gr_patches_indices_local(:)
CHARACTER(LEN=64), allocatable :: gr_profiles(:)

! fields config in the fraction of ISM/JET n_H
real(kind=gr_rpknd), allocatable :: gr_epsilon_tol(:) !density fraction tolerance
real(kind=gr_rpknd), allocatable :: gr_x_velocity(:)
real(kind=gr_rpknd), allocatable :: gr_y_velocity(:)
real(kind=gr_rpknd), allocatable :: gr_z_velocity(:)
CHARACTER(LEN=30), allocatable :: gr_density_method(:)
!real(kind=gr_rpknd), allocatable :: gr_density(:)!> this is furnished by AMRVAC
real(kind=gr_rpknd), allocatable :: gr_HI_density(:)
real(kind=gr_rpknd), allocatable :: gr_HII_density(:)
real(kind=gr_rpknd), allocatable :: gr_HM_density(:)
real(kind=gr_rpknd), allocatable :: gr_HeI_density(:)
real(kind=gr_rpknd), allocatable :: gr_HeII_density(:)
real(kind=gr_rpknd), allocatable :: gr_HeIII_density(:)
real(kind=gr_rpknd), allocatable :: gr_H2I_density(:)
real(kind=gr_rpknd), allocatable :: gr_H2II_density(:)
real(kind=gr_rpknd), allocatable :: gr_DI_density(:)
real(kind=gr_rpknd), allocatable :: gr_DII_density(:)
real(kind=gr_rpknd), allocatable :: gr_HDI_density(:)
real(kind=gr_rpknd), allocatable :: gr_e_density(:)
real(kind=gr_rpknd), allocatable :: gr_metal_density(:)
real(kind=gr_rpknd), allocatable :: gr_dust_density(:)
real(kind=gr_rpknd), allocatable :: gr_volumetric_heating_rate(:)
real(kind=gr_rpknd), allocatable :: gr_specific_heating_rate(:)
real(kind=gr_rpknd), allocatable :: gr_RT_HI_ionization_rate(:)
real(kind=gr_rpknd), allocatable :: gr_RT_HeI_ionization_rate(:)
real(kind=gr_rpknd), allocatable :: gr_RT_HeII_ionization_rate(:)
real(kind=gr_rpknd), allocatable :: gr_RT_H2_dissociation_rate(:)
real(kind=gr_rpknd), allocatable :: gr_RT_heating_rate(:)






type gr_objects
  logical, allocatable            :: patch(:,:^D&)           !> spatial patch
  character(len=78)               :: subname                 !> subroutine name that call it
  type(gr_config) :: myconfig
  type(gr_params) :: myparams
  type(gr_fields) :: myfields
  type(usrphysical_unit), pointer :: myphysunit
  contains
  !PRIVATE
!=BEGIN PREINITIALIZATION PROCESS SUBROUTINES=====================================================================================
  PROCEDURE, PASS(self) :: check_indexes          => check_indexes_consistency
  PROCEDURE, PASS(self) :: dealloc_grids          => deallocate_grids
  PROCEDURE, PASS(self) :: alloc_grids            => allocate_grids
  PROCEDURE, PASS(self) :: reset_grid             => reset_field_size
  !     Create a grackle chemistry object for parameters and set defaults
  PROCEDURE, PASS(self) :: default_filename       => grackle_object_default_filename
  PROCEDURE, PASS(self) :: default_data           => grackle_data_set_default
  !     Set default units
  PROCEDURE, PASS(self) :: default_units          => grackle_units_set_default
  !     Initialize the Grackle and set default field arrays
  PROCEDURE, PASS(self) :: init_chem              => grackle_init_chem
  PROCEDURE, PASS(self) :: default_fields         => grackle_fields_set_default
  !     Defaulting every preinitialization default parameters, data and field
  PROCEDURE, PASS(self) :: default_all_one            => grackle_set_all_default_part_one
  PROCEDURE, PASS(self) :: default_all_two            => grackle_set_all_default_part_two
  PROCEDURE, PASS(self) :: link_to_grackle        => link_gr_obj_to_gr_struct
!=END PREINITIALIZATION PROCESS SUBROUTINES=====================================================================================

  PROCEDURE, PASS(self) :: set_default_config => grackle_set_default
  PROCEDURE, PASS(self) :: read_parameters    => grackle_params_read
  PROCEDURE, PASS(self) :: allocate_fields_parameters => grackle_fields_config_allocate
  PROCEDURE, PASS(self) :: set_fields_default => grackle_fields_config_set_default
  PROCEDURE, PASS(self) :: read_fields_parameters => grackle_fields_config_read
  PROCEDURE, PASS(self) :: set_complet        => grackle_set_complet
  PROCEDURE, PASS(self) :: normalize          => grackle_normalize
  PROCEDURE, PASS(self) :: test_chemistry         => use_chemistry
  PROCEDURE, PASS(self) :: write_setting        => grackle_write_setting
end type gr_objects


contains



!=BEGIN PREINITIALIZATION PROCESS SUBROUTINES=====================================================================================
subroutine check_indexes_consistency(ixI^L,ixO^L,subroutinename,modulename,parentname,self)
  implicit none
  integer, intent(in)           :: ixI^L,ixO^L
  character(len=*), intent(in) :: subroutinename,modulename
  character(len=*), intent(in) :: parentname
  class(gr_objects), TARGET :: self

  {^D&

  if(ixImax^D<ixImin^D)then
    write(*,*) '=================================================='
    write(*,*) '0 Called from ', parentname
    write(*,*) '> In module : ', modulename
    write(*,*) '>> In subroutine : ', subroutinename
    write(*,*) '>>> we get ixImin',^D,' = ', ixImin^D
    write(*,*) '>>> we get ixImax',^D,' = ', ixImax^D
    call mpistop('>>> Error : Inconsistency with ixImax^D<ixImin^D')
    write(*,*) '=================================================='
  end if

  if(ixOmax^D<ixOmin^D)then
    write(*,*) '=================================================='
    write(*,*) '0 Called from ', parentname
    write(*,*) '> In module : ', modulename
    write(*,*) '>> In subroutine : ', subroutinename
    write(*,*) '>>> we get ixImin',^D,' = ', ixOmin^D
    write(*,*) '>>> we get ixImax',^D,' = ', ixOmax^D
    call mpistop('>>> Error : Inconsistency with ixOmax^D<ixOmin^D')
    write(*,*) '=================================================='
  end if

  }

end

subroutine deallocate_grids(self)
  implicit none
  class(gr_objects), TARGET :: self
  !=================================================================================

  if(allocated(self%myfields%gr_density))deallocate(self%myfields%gr_density)
  if(allocated(self%myfields%gr_energy))deallocate(self%myfields%gr_energy)
  if(allocated(self%myfields%gr_x_velocity))deallocate(self%myfields%gr_x_velocity)
  if(allocated(self%myfields%gr_y_velocity))deallocate(self%myfields%gr_y_velocity)
  if(allocated(self%myfields%gr_z_velocity))deallocate(self%myfields%gr_z_velocity)
  if(allocated(self%myfields%gr_HI_density))deallocate(self%myfields%gr_HI_density)
  if(allocated(self%myfields%gr_HII_density))deallocate(self%myfields%gr_HII_density)
  if(allocated(self%myfields%gr_HM_density))deallocate(self%myfields%gr_HM_density)
  if(allocated(self%myfields%gr_HeI_density))deallocate(self%myfields%gr_HeI_density)
  if(allocated(self%myfields%gr_HeII_density))deallocate(self%myfields%gr_HeII_density)
  if(allocated(self%myfields%gr_HeIII_density))deallocate(self%myfields%gr_HeIII_density)
  if(allocated(self%myfields%gr_H2I_density))deallocate(self%myfields%gr_H2I_density)
  if(allocated(self%myfields%gr_H2II_density))deallocate(self%myfields%gr_H2II_density)
  if(allocated(self%myfields%gr_DI_density))deallocate(self%myfields%gr_DI_density)
  if(allocated(self%myfields%gr_DII_density))deallocate(self%myfields%gr_DII_density)
  if(allocated(self%myfields%gr_HDI_density))deallocate(self%myfields%gr_HDI_density)
  if(allocated(self%myfields%gr_e_density))deallocate(self%myfields%gr_e_density)
  if(allocated(self%myfields%gr_metal_density))deallocate(self%myfields%gr_metal_density)
  if(allocated(self%myfields%gr_dust_density))deallocate(self%myfields%gr_dust_density)
  if(allocated(self%myfields%gr_volumetric_heating_rate))deallocate(self%myfields%gr_volumetric_heating_rate)
  if(allocated(self%myfields%gr_specific_heating_rate))deallocate(self%myfields%gr_specific_heating_rate)
  if(allocated(self%myfields%gr_RT_HI_ionization_rate))deallocate(self%myfields%gr_RT_HI_ionization_rate)
  if(allocated(self%myfields%gr_RT_HeI_ionization_rate))deallocate(self%myfields%gr_RT_HeI_ionization_rate)
  if(allocated(self%myfields%gr_RT_HeII_ionization_rate))deallocate(self%myfields%gr_RT_HeII_ionization_rate)
  if(allocated(self%myfields%gr_RT_H2_dissociation_rate))deallocate(self%myfields%gr_RT_H2_dissociation_rate)
  if(allocated(self%myfields%gr_RT_heating_rate))deallocate(self%myfields%gr_RT_heating_rate)
  if(allocated(self%myfields%gr_cooling_time))deallocate(self%myfields%gr_cooling_time)
  if(allocated(self%myfields%gr_gamma))deallocate(self%myfields%gr_gamma)
  if(allocated(self%myfields%gr_pressure))deallocate(self%myfields%gr_pressure)
  if(allocated(self%myfields%gr_temperature))deallocate(self%myfields%gr_temperature)
  if(allocated(self%myfields%gr_dust_temperature))deallocate(self%myfields%gr_dust_temperature)

  !=================================================================================
end subroutine deallocate_grids

subroutine allocate_grids(self)
  implicit none
  class(gr_objects), TARGET :: self
  !=================================================================================

  if(.not.allocated(self%myfields%gr_density))allocate(self%myfields%gr_density({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_energy))allocate(self%myfields%gr_energy({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_x_velocity))allocate(self%myfields%gr_x_velocity({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_y_velocity))allocate(self%myfields%gr_y_velocity({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_z_velocity))allocate(self%myfields%gr_z_velocity({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_HI_density))allocate(self%myfields%gr_HI_density({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_HII_density))allocate(self%myfields%gr_HII_density({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_HM_density))allocate(self%myfields%gr_HM_density({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_HeI_density))allocate(self%myfields%gr_HeI_density({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_HeII_density))allocate(self%myfields%gr_HeII_density({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_HeIII_density))allocate(self%myfields%gr_HeIII_density({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_H2I_density))allocate(self%myfields%gr_H2I_density({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_H2II_density))allocate(self%myfields%gr_H2II_density({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_DI_density))allocate(self%myfields%gr_DI_density({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_DII_density))allocate(self%myfields%gr_DII_density({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_HDI_density))allocate(self%myfields%gr_HDI_density({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_e_density))allocate(self%myfields%gr_e_density({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_metal_density))allocate(self%myfields%gr_metal_density({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_dust_density))allocate(self%myfields%gr_dust_density({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_volumetric_heating_rate))allocate(self%myfields%gr_volumetric_heating_rate({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_specific_heating_rate))allocate(self%myfields%gr_specific_heating_rate({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_RT_HI_ionization_rate))allocate(self%myfields%gr_RT_HI_ionization_rate({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_RT_HeI_ionization_rate))allocate(self%myfields%gr_RT_HeI_ionization_rate({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_RT_HeII_ionization_rate))allocate(self%myfields%gr_RT_HeII_ionization_rate({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_RT_H2_dissociation_rate))allocate(self%myfields%gr_RT_H2_dissociation_rate({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_RT_heating_rate))allocate(self%myfields%gr_RT_heating_rate({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_cooling_time))allocate(self%myfields%gr_cooling_time({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_gamma))allocate(self%myfields%gr_gamma({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_pressure))allocate(self%myfields%gr_pressure({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_temperature))allocate(self%myfields%gr_temperature({self%myfields%field_size(^D)|,}))
  if(.not.allocated(self%myfields%gr_dust_temperature))allocate(self%myfields%gr_dust_temperature({self%myfields%field_size(^D)|,}))

  !=================================================================================
end subroutine allocate_grids

subroutine reset_field_size(self)
  implicit none
  class(gr_objects), TARGET :: self

  call self%dealloc_grids()
  call self%alloc_grids()

end subroutine reset_field_size

subroutine grackle_object_default_filename(self)
  implicit none
  integer                        :: iresult
  class(gr_objects), TARGET      :: self

  !     cooling data for Haardt & Madau 2012 background
  ! must adapt and change if grackle src directory location
  ! is moved !!!
  self%myfields%filename = "../../../src/grackle/input/"//"CloudyData_UVB=HM2012.h5"//C_NULL_CHAR

end subroutine grackle_object_default_filename

subroutine grackle_data_set_default(gr_struct,self)
  implicit none
  integer                        :: iresult
  class(grackle_type)            :: gr_struct
  class(gr_objects), TARGET      :: self
  !----------------------------------

  !     Create a grackle chemistry object for parameters and set defaults

  iresult = set_default_chemistry_parameters(gr_struct%grackle_data)

  !     Set parameters

  gr_struct%grackle_data%use_grackle = 1            ! chemistry on
  gr_struct%grackle_data%with_radiative_cooling = 1 ! cooling on
  gr_struct%grackle_data%primordial_chemistry = 3   ! network with H, He, D
  gr_struct%grackle_data%dust_chemistry = 1         ! dust processes
  gr_struct%grackle_data%metal_cooling = 1          ! metal cooling on
  gr_struct%grackle_data%UVbackground = 1           ! UV background on

  gr_struct%grackle_data%grackle_data_file = C_LOC(self%myfields%filename(1:1))
  gr_struct%grackle_data%h2_on_dust = 0             ! no dust
  gr_struct%grackle_data%cmb_temperature_floor = 1  ! include CMB cooling floor
  gr_struct%grackle_data%Gamma = 5./3.;          ! monoatomic gas

end subroutine grackle_data_set_default

subroutine grackle_units_set_default(gr_struct,self)
  implicit none
  class(grackle_type)            :: gr_struct
  class(gr_objects), TARGET      :: self
  !----------------------------------

   gr_struct%units%comoving_coordinates = 0
   gr_struct%units%density_units = 1.67d-24
   gr_struct%units%length_units = 1.0d0
   gr_struct%units%time_units = 1.0d12
   gr_struct%units%a_units = 1.0d0

   !     Set initial expansion factor (for internal units).
   !     Set expansion factor to 1 for non-cosmological simulation.
   self%myparams%initial_redshift = 0.;
   gr_struct%units%a_value = 1. / (1. + self%myparams%initial_redshift);
   call set_velocity_units(gr_struct%units)

end subroutine grackle_units_set_default

subroutine grackle_init_chem(gr_struct,self)
  implicit none
  integer                        :: iresult
  class(grackle_type)            :: gr_struct
  class(gr_objects), TARGET      :: self

  !     Initialize the Grackle

  iresult = initialize_chemistry_data(gr_struct%units)

end subroutine grackle_init_chem

subroutine grackle_fields_set_default(gr_struct,self)
  implicit none
  integer                        :: iresult,ifield, irank^D
  class(grackle_type)            :: gr_struct
  class(gr_objects), TARGET      :: self
  !----------------------------------


  !     Set field arrays

  !     If grid rank is less than 3, set the other dimensions,
  !     start indices, and end indices to 0.
  self%myfields%grid_rank = 3
  do ifield = 1, self%myfields%grid_rank
   self%myfields%grid_dimension(ifield) = 1
   self%myfields%grid_start(ifield) = 0
   self%myfields%grid_end(ifield) = 0
  enddo
  self%myparams%grid_dx = 0.0
  {^D&
  self%myfields%grid_dimension(^D) = self%myfields%field_size(^D)
  !     0-based
  self%myfields%grid_end(^D) = self%myfields%field_size(^D) - 1
  }

  self%myparams%temperature_units = get_temperature_units(gr_struct%units)

  {^D& do irank^D = 1,self%myfields%field_size(^D)\}
   self%myfields%gr_density(irank^D) = 1.0
   self%myfields%gr_HI_density(irank^D) = gr_struct%grackle_data%HydrogenFractionByMass * &
   self%myfields%gr_density(irank^D)
   self%myfields%gr_HII_density(irank^D) = tiny_number * self%myfields%gr_density(irank^D)
   self%myfields%gr_HM_density(irank^D) = tiny_number * self%myfields%gr_density(irank^D)
   self%myfields%gr_HeI_density(irank^D) = (1.0 - gr_struct%grackle_data%HydrogenFractionByMass) * &
   self%myfields%gr_density(irank^D)
   self%myfields%gr_HeII_density(irank^D) = tiny_number * self%myfields%gr_density(irank^D)
   self%myfields%gr_HeIII_density(irank^D) = tiny_number * self%myfields%gr_density(irank^D)
   self%myfields%gr_H2I_density(irank^D) = tiny_number * self%myfields%gr_density(irank^D)
   self%myfields%gr_H2II_density(irank^D) = tiny_number * self%myfields%gr_density(irank^D)
   self%myfields%gr_DI_density(irank^D) = 2.0 * 3.4e-5 * self%myfields%gr_density(irank^D)
   self%myfields%gr_DII_density(irank^D) = tiny_number * self%myfields%gr_density(irank^D)
   self%myfields%gr_HDI_density(irank^D) = tiny_number * self%myfields%gr_density(irank^D)
   self%myfields%gr_e_density(irank^D) = tiny_number * self%myfields%gr_density(irank^D)
   !        solar metallicity
   self%myfields%gr_metal_density(irank^D) = gr_struct%grackle_data%SolarMetalFractionByMass * &
           self%myfields%gr_density(irank^D)
   self%myfields%gr_dust_density(irank^D) = gr_struct%grackle_data%local_dust_to_gas_ratio * &
           self%myfields%gr_density(irank^D)

   self%myfields%gr_x_velocity(irank^D) = 0.0
   self%myfields%gr_y_velocity(irank^D) = 0.0
   self%myfields%gr_z_velocity(irank^D) = 0.0

   !        initilize internal energy (here 1000 K for no reason)
   self%myfields%gr_energy(irank^D) = 1000. / self%myparams%temperature_units

   self%myfields%gr_volumetric_heating_rate(irank^D) = 0.0
   self%myfields%gr_specific_heating_rate(irank^D) = 0.0
   self%myfields%gr_RT_HI_ionization_rate(irank^D) = 0.0
   self%myfields%gr_RT_HeI_ionization_rate(irank^D) = 0.0
   self%myfields%gr_RT_HeII_ionization_rate(irank^D) = 0.0
   self%myfields%gr_RT_H2_dissociation_rate(irank^D) = 0.0
   self%myfields%gr_RT_heating_rate(irank^D) = 0.0
  {^D& end do\}

end subroutine grackle_fields_set_default

subroutine grackle_set_all_default_part_one(gr_struct,self)
  implicit none
  class(grackle_type)            :: gr_struct
  class(gr_objects), TARGET      :: self

  !     Create a grackle chemistry object for parameters and set defaults
  call self%default_filename()
  call self%default_data(gr_struct)
  !     Set units
  call self%default_units(gr_struct)
end subroutine grackle_set_all_default_part_one

subroutine grackle_set_all_default_part_two(gr_struct,self)
  implicit none
  class(grackle_type)            :: gr_struct
  class(gr_objects), TARGET      :: self

  call self%default_fields(gr_struct)

end subroutine grackle_set_all_default_part_two

subroutine link_gr_obj_to_gr_struct(gr_struct,self)
  implicit none
  class(grackle_type)            :: gr_struct
  class(gr_objects), TARGET      :: self

  !
  !     Fill in structure to be passed to Grackle
  !

  gr_struct%fields%grid_rank = ndim
  gr_struct%fields%grid_dimension = C_LOC(self%myfields%grid_dimension)
  gr_struct%fields%grid_start = C_LOC(self%myfields%grid_start)
  gr_struct%fields%grid_end = C_LOC(self%myfields%grid_end)
  gr_struct%fields%grid_dx  = self%myparams%grid_dx

  gr_struct%fields%density = C_LOC(self%myfields%gr_density)
  gr_struct%fields%HI_density = C_LOC(self%myfields%gr_HI_density)
  gr_struct%fields%HII_density = C_LOC(self%myfields%gr_HII_density)
  gr_struct%fields%HM_density = C_LOC(self%myfields%gr_HM_density)
  gr_struct%fields%HeI_density = C_LOC(self%myfields%gr_HeI_density)
  gr_struct%fields%HeII_density = C_LOC(self%myfields%gr_HeII_density)
  gr_struct%fields%HeIII_density = C_LOC(self%myfields%gr_HeIII_density)
  gr_struct%fields%H2I_density = C_LOC(self%myfields%gr_H2I_density)
  gr_struct%fields%H2II_density = C_LOC(self%myfields%gr_H2II_density)
  gr_struct%fields%DI_density = C_LOC(self%myfields%gr_DI_density)
  gr_struct%fields%DII_density = C_LOC(self%myfields%gr_DII_density)
  gr_struct%fields%HDI_density = C_LOC(self%myfields%gr_HDI_density)
  gr_struct%fields%e_density = C_LOC(self%myfields%gr_e_density)
  gr_struct%fields%metal_density = C_LOC(self%myfields%gr_metal_density)
  gr_struct%fields%internal_energy = C_LOC(self%myfields%gr_energy)
  gr_struct%fields%x_velocity = C_LOC(self%myfields%gr_x_velocity)
  gr_struct%fields%y_velocity = C_LOC(self%myfields%gr_y_velocity)
  gr_struct%fields%z_velocity = C_LOC(self%myfields%gr_z_velocity)
  gr_struct%fields%volumetric_heating_rate = &
                                 C_LOC(self%myfields%gr_volumetric_heating_rate)
  gr_struct%fields%specific_heating_rate = C_LOC(self%myfields%gr_specific_heating_rate)
  gr_struct%fields%RT_HI_ionization_rate = C_LOC(self%myfields%gr_RT_HI_ionization_rate)
  gr_struct%fields%RT_HeI_ionization_rate = C_LOC(self%myfields%gr_RT_HeI_ionization_rate)
  gr_struct%fields%RT_HeII_ionization_rate = &
                                   C_LOC(self%myfields%gr_RT_HeII_ionization_rate)
  gr_struct%fields%RT_H2_dissociation_rate = C_LOC(self%myfields%gr_RT_H2_dissociation_rate)
  gr_struct%fields%RT_heating_rate = C_LOC(self%myfields%gr_RT_heating_rate)

end subroutine link_gr_obj_to_gr_struct

!=END PREINITIALIZATION PROCESS SUBROUTINES=====================================================================================






































!-------------------------------------------------------------------------
!> subroutine default setting for Grackle
subroutine grackle_set_default(self)
  implicit none
  class(gr_objects),TARGET            :: self
  !----------------------------------
  self%myconfig%obj_name              = 'grackle_chemistry_config'
  self%myconfig%unit                  = 'code'
  self%myconfig%myindice              = 0

  self%myconfig%use_grackle = 1
  self%myconfig%gr_with_radiative_cooling = 1
  self%myconfig%gr_primordial_chemistry = 0
  self%myconfig%gr_metal_cooling = 1
  self%myconfig%gr_UVbackground = 0
  self%myconfig%gr_cmb_temperature_floor = 0
  self%myconfig%gr_h2_on_dust = 0
  self%myconfig%gr_dust_chemistry = 0
  self%myconfig%gr_use_dust_density_field = 0
  self%myconfig%gr_photoelectric_heating = 0
  self%myconfig%gr_use_volumetric_heating_rate = 0
  self%myconfig%gr_use_specific_heating_rate = 0
  self%myconfig%gr_three_body_rate = 0
  self%myconfig%gr_cie_cooling = 0
  self%myconfig%gr_h2_optical_depth_approximation = 0
  self%myconfig%gr_ih2co = 1
  self%myconfig%gr_ipiht = 1
  self%myconfig%gr_NumberOfTemperatureBins = 600
  self%myconfig%gr_CaseBRecombination = 0
  self%myconfig%gr_Compton_xray_heating = 0
  self%myconfig%gr_LWbackground_sawtooth_suppression = 0
  self%myconfig%gr_NumberOfDustTemperatureBins = 250
  self%myconfig%gr_use_radiative_transfer = 0
  self%myconfig%gr_radiative_transfer_coupled_rate_solver = 0
  self%myconfig%gr_radiative_transfer_intermediate_step = 0
  self%myconfig%gr_radiative_transfer_hydrogen_only = 0
  self%myconfig%gr_self_shielding_method = 0
  self%myconfig%gr_Gamma = 5.d0/3.d0
  self%myconfig%gr_photoelectric_heating_rate = 8.5D-26
  self%myconfig%gr_HydrogenFractionByMass = 0.76d0
  self%myconfig%gr_DeuteriumToHydrogenRatio = 2.0d0*3.4d-5
  self%myconfig%gr_SolarMetalFractionByMass = 0.01295d0
  self%myconfig%TemperatureStart = 1.0d0
  self%myconfig%TemperatureEnd = 1.0D9
  self%myconfig%DustTemperatureStart = 1.0d0
  self%myconfig%DustTemperatureEnd = 1500.0d0
  self%myconfig%gr_LWbackground_intensity = 0.0d0
  self%myconfig%gr_UVbackground_redshift_on = 7.0d0
  self%myconfig%gr_UVbackground_redshift_off = 0.0d0
  self%myconfig%gr_UVbackground_redshift_fullon = 6.0d0
  self%myconfig%gr_UVbackground_redshift_drop = 0.0d0
  self%myconfig%cloudy_electron_fraction_factor = 9.153959D-3
  self%myconfig%data_dir = "../../../src/grackle/input/"
  self%myconfig%data_filename = "CloudyData_UVB=HM2012.h5"
  self%myparams%data_file = self%myconfig%data_dir//self%myconfig%data_filename//C_NULL_CHAR
  self%myconfig%normalize_done = .false.
  self%myconfig%gr_comoving_coordinates = 0
  self%myconfig%gr_a_units = 0.0d0
  self%myconfig%gr_current_redshift = 0.
  self%myparams%number_of_objects = 0

  write(*,*) 'Grackle configuration defaulting successfully done !'



end subroutine grackle_set_default



!-------------------------------------------------------------------------
 !> Read grackle parameters  from a parfile
 subroutine grackle_params_read(self,grackle_config,files)
   implicit none
   class(gr_objects),TARGET                         :: self
   character(len=*),intent(in)        :: files(:)
   type(gr_config), intent(out)  :: grackle_config

   ! .. local ..
   integer                            :: i_file,i_error_read
   integer                  :: idim,iside
   character(len=70)            :: error_message
   !-------------------------------------------------------------------------
   namelist /grackle_par_list/  grackle_config,tstst
   namelist /grackle_par1_list/ grackle_config,tstst
   namelist /grackle_par2_list/ grackle_config,tstst
   namelist /grackle_par3_list/ grackle_config,tstst


   error_message = 'In the procedure : grackle_params_read'

   if(mype==0)write(*,*)'Reading grackle_par_list'
   do i_file = 1, size(files)
      open(unitpar, file=trim(files(i_file)), status="old")
      select case(grackle_config%myindice)
      case(1)
        read(unitpar, grackle_par1_list, iostat=i_error_read)
      case(2)
        read(unitpar, grackle_par2_list, iostat=i_error_read)
      case(3)
        read(unitpar, grackle_par3_list, iostat=i_error_read)
      case default
        read(unitpar, grackle_par_list, iostat=i_error_read)
        write(*,*) 'use_grackle : ', grackle_config%use_grackle
        write(*,*) 'gr_primordial_chemistry : ', grackle_config%gr_primordial_chemistry
        write(*,*) 'gr_metal_cooling : ', grackle_config%gr_metal_cooling
        write(*,*) 'gr_dust_chemistry : ', grackle_config%gr_dust_chemistry
      end select
      call usr_mat_read_error_message(i_error_read,grackle_config%myindice,&
                                      self%myconfig%obj_name)
      close(unitpar)

      !add a routine to check that the sum of each density is equal to one times
      ! the total density in ISM AND JET
   end do

   write(*,*) 'tstst = ', tstst
   write(*,*) 'Grackle usr configuration reading successfully done !'
 end subroutine grackle_params_read


 subroutine grackle_fields_config_allocate(self,number_of_isms,number_of_jets,number_of_clouds)
  implicit none
  class(gr_objects),TARGET           :: self
  integer,intent(inout),optional        :: number_of_isms
  integer,intent(inout),optional        :: number_of_jets
  integer,intent(inout),optional        :: number_of_clouds
  integer        :: default_number_of_isms
  integer        :: default_number_of_jets
  integer        :: default_number_of_clouds
  ! .. local ..
  integer                  :: n_total_objects
  !-------------------------------------------------------------------------

  ! Defaults
  default_number_of_isms = 1
  default_number_of_jets = 0
  default_number_of_clouds = 0
  if(.not.present(number_of_isms))number_of_isms = default_number_of_isms
  if(.not.present(number_of_jets))number_of_jets = default_number_of_jets
  if(.not.present(number_of_clouds))number_of_clouds = default_number_of_clouds


  n_total_objects = 0
  i_zero_ism = 1
  i_zero_jet = 1
  i_zero_cloud = 1
  i_end_ism = 1
  i_end_jet = 1
  i_end_cloud = 1

  ! ex of slicing :
  ! 1 -- ISM --> 4; 5-- JET-->7; 8-- CLOUDS--> 9
  if(number_of_isms>0)then
    i_zero_ism = i_zero_ism ! ex: i_zero_ism = 1
    n_total_objects = n_total_objects + number_of_isms ! ex : n_tot = 0 + 4 = 4
    i_end_ism = n_total_objects
  end if
  if(number_of_jets>0)then
    i_zero_jet = i_zero_ism + number_of_isms ! ex: i_zero_jet = 1 + 4 = 5
    n_total_objects = n_total_objects + number_of_jets ! ex : n_tot = 0 + 4 + 3 = 7
    i_end_jet = n_total_objects
  end if
  if(number_of_clouds>0)then
    i_zero_cloud = i_zero_ism + number_of_isms + number_of_jets
    ! ex: i_zero_cloud = 1 + 4 + 3 = 8
    n_total_objects = n_total_objects + number_of_clouds ! ex : n_tot = 0 + 4 + 3 + 2 = 9
    i_end_cloud = n_total_objects
  end if
  if(n_total_objects==0)then
   WRITE(*,*) 'Error, the code found 0 objects ism+jets+clouds, n_total_objects = ', n_total_objects
   call mpistop('Incoherence in the number of objects. The code stops.')
  end if

  ! Store the total number of objects
  self%myparams%number_of_objects = n_total_objects

  ! First, deallocate already allocated arrays


  if(allocated(gr_patches_indices_global))deallocate(gr_patches_indices_global)
  if(allocated(gr_patches_indices_local))deallocate(gr_patches_indices_local)
  if(allocated(gr_patches_name))deallocate(gr_patches_name)
  if(allocated(gr_profiles))deallocate(gr_profiles)
  if(allocated(gr_epsilon_tol))deallocate(gr_epsilon_tol)
  if(allocated(gr_density_method))deallocate(gr_density_method)
  if(allocated(gr_x_velocity))deallocate(gr_x_velocity)
  if(allocated(gr_y_velocity))deallocate(gr_y_velocity)
  if(allocated(gr_z_velocity))deallocate(gr_z_velocity)
  if(allocated(gr_HI_density))deallocate(gr_HI_density)
  if(allocated(gr_HII_density))deallocate(gr_HII_density)
  if(allocated(gr_HM_density))deallocate(gr_HM_density)
  if(allocated(gr_HeI_density))deallocate(gr_HeI_density)
  if(allocated(gr_HeII_density))deallocate(gr_HeII_density)
  if(allocated(gr_HeIII_density))deallocate(gr_HeIII_density)
  if(allocated(gr_H2I_density))deallocate(gr_H2I_density)
  if(allocated(gr_H2II_density))deallocate(gr_H2II_density)
  if(allocated(gr_DI_density))deallocate(gr_DI_density)
  if(allocated(gr_DII_density))deallocate(gr_DII_density)
  if(allocated(gr_HDI_density))deallocate(gr_HDI_density)
  if(allocated(gr_e_density))deallocate(gr_e_density)
  if(allocated(gr_metal_density))deallocate(gr_metal_density)
  if(allocated(gr_dust_density))deallocate(gr_dust_density)
  if(allocated(gr_volumetric_heating_rate))deallocate(gr_volumetric_heating_rate)
  if(allocated(gr_specific_heating_rate))deallocate(gr_specific_heating_rate)
  if(allocated(gr_RT_HI_ionization_rate))deallocate(gr_RT_HI_ionization_rate)
  if(allocated(gr_RT_HeI_ionization_rate))deallocate(gr_RT_HeI_ionization_rate)
  if(allocated(gr_RT_HeII_ionization_rate))deallocate(gr_RT_HeII_ionization_rate)
  if(allocated(gr_RT_H2_dissociation_rate))deallocate(gr_RT_H2_dissociation_rate)
  if(allocated(gr_RT_heating_rate))deallocate(gr_RT_heating_rate)

  ! Second, allocate the size of the allocatable arrays


  if(.not.allocated(gr_patches_indices_global))allocate(gr_patches_indices_global(1:n_total_objects))
  if(.not.allocated(gr_patches_indices_local))allocate(gr_patches_indices_local(1:n_total_objects))
  if(.not.allocated(gr_patches_name))allocate(gr_patches_name(1:n_total_objects))
  if(.not.allocated(gr_profiles))allocate(gr_profiles(1:n_total_objects))
  if(.not.allocated(gr_epsilon_tol))allocate(gr_epsilon_tol(1:n_total_objects))
  if(.not.allocated(gr_density_method))allocate(gr_density_method(1:n_total_objects))
  if(.not.allocated(gr_x_velocity))allocate(gr_x_velocity(1:n_total_objects))
  if(.not.allocated(gr_y_velocity))allocate(gr_y_velocity(1:n_total_objects))
  if(.not.allocated(gr_z_velocity))allocate(gr_z_velocity(1:n_total_objects))
  if(.not.allocated(gr_HI_density))allocate(gr_HI_density(1:n_total_objects))
  if(.not.allocated(gr_HII_density))allocate(gr_HII_density(1:n_total_objects))
  if(.not.allocated(gr_HM_density))allocate(gr_HM_density(1:n_total_objects))
  if(.not.allocated(gr_HeI_density))allocate(gr_HeI_density(1:n_total_objects))
  if(.not.allocated(gr_HeII_density))allocate(gr_HeII_density(1:n_total_objects))
  if(.not.allocated(gr_HeIII_density))allocate(gr_HeIII_density(1:n_total_objects))
  if(.not.allocated(gr_H2I_density))allocate(gr_H2I_density(1:n_total_objects))
  if(.not.allocated(gr_H2II_density))allocate(gr_H2II_density(1:n_total_objects))
  if(.not.allocated(gr_DI_density))allocate(gr_DI_density(1:n_total_objects))
  if(.not.allocated(gr_DII_density))allocate(gr_DII_density(1:n_total_objects))
  if(.not.allocated(gr_HDI_density))allocate(gr_HDI_density(1:n_total_objects))
  if(.not.allocated(gr_e_density))allocate(gr_e_density(1:n_total_objects))
  if(.not.allocated(gr_metal_density))allocate(gr_metal_density(1:n_total_objects))
  if(.not.allocated(gr_dust_density))allocate(gr_dust_density(1:n_total_objects))
  if(.not.allocated(gr_volumetric_heating_rate))allocate(gr_volumetric_heating_rate(1:n_total_objects))
  if(.not.allocated(gr_specific_heating_rate))allocate(gr_specific_heating_rate(1:n_total_objects))
  if(.not.allocated(gr_RT_HI_ionization_rate))allocate(gr_RT_HI_ionization_rate(1:n_total_objects))
  if(.not.allocated(gr_RT_HeI_ionization_rate))allocate(gr_RT_HeI_ionization_rate(1:n_total_objects))
  if(.not.allocated(gr_RT_HeII_ionization_rate))allocate(gr_RT_HeII_ionization_rate(1:n_total_objects))
  if(.not.allocated(gr_RT_H2_dissociation_rate))allocate(gr_RT_H2_dissociation_rate(1:n_total_objects))
  if(.not.allocated(gr_RT_heating_rate))allocate(gr_RT_heating_rate(1:n_total_objects))


 end subroutine grackle_fields_config_allocate

 subroutine grackle_fields_config_set_default(self,ism_present,jet_present,cloud_present)
  implicit none
  class(gr_objects),TARGET           :: self
  logical,intent(in)                 :: ism_present,jet_present,cloud_present
  ! .. local ..
  integer                            :: i_file,i_error_read,i_ism,i_jet,i_cloud,i_all
  integer                            :: idim,iside,n_total_objects,counter
  character(len=70)                  :: error_message
  !-------------------------------------------------------------------------




  ! First, fill in specifically for ism

  counter = 0
  ism_is_present :if(ism_present)then
    Loop_fill_ism : do i_ism=i_zero_ism,i_end_ism

      counter=counter+1
      gr_patches_name(i_ism) = 'ism'
      gr_patches_indices_global(i_ism) = i_ism
      gr_patches_indices_local(i_ism) = counter

    end do Loop_fill_ism
  end if ism_is_present

  ! do the same  specifically for jet
  counter = 0
  jet_is_present :if(jet_present)then
    Loop_fill_jet : do i_jet=i_zero_jet,i_end_jet

      counter=counter+1
      gr_patches_name(i_jet) = 'jet'
      gr_patches_indices_global(i_jet) = i_jet
      gr_patches_indices_local(i_jet) = counter

    end do Loop_fill_jet
  end if jet_is_present


  ! do the same  specifically for cloud
  counter = 0
  cloud_is_present :if(cloud_present)then
    Loop_fill_cloud : do i_cloud=i_zero_cloud,i_end_cloud

      counter=counter+1
      gr_patches_name(i_cloud) = 'cloud'
      gr_patches_indices_global(i_cloud) = i_cloud
      gr_patches_indices_local(i_cloud) = counter

    end do Loop_fill_cloud
  end if cloud_is_present


   ! Second, fill in default values for everything, independently

   ! Fields parameters

   Loop_fill_everything : do i_all=1,size(gr_patches_indices_global)

     gr_profiles(gr_patches_indices_global(i_all)) = 'uniform'
     ! fields config
     gr_epsilon_tol(gr_patches_indices_global(i_all)) = 14.01*tiny_number !density fraction tolerance
     gr_x_velocity(gr_patches_indices_global(i_all)) = 0.
     gr_y_velocity(gr_patches_indices_global(i_all)) = 0.
     gr_z_velocity(gr_patches_indices_global(i_all)) = 0.

     ! Densities
     gr_density_method(gr_patches_indices_global(i_all)) = 'chemical_coefficient' !method to estimate following densities
     gr_HI_density(gr_patches_indices_global(i_all)) = 0.
     gr_HII_density(gr_patches_indices_global(i_all)) = 0.
     gr_HM_density(gr_patches_indices_global(i_all)) = 0.
     gr_H2I_density(gr_patches_indices_global(i_all)) = 0.
     gr_H2II_density(gr_patches_indices_global(i_all)) = 0.
     gr_DI_density(gr_patches_indices_global(i_all)) = 0.
     gr_DII_density(gr_patches_indices_global(i_all)) = 0.
     gr_HDI_density(gr_patches_indices_global(i_all)) = 0.
     gr_e_density(gr_patches_indices_global(i_all)) = 0.
     ! Densities as fraction of n_H=n(H+H^+)+2n(H_2)
     gr_HeI_density(gr_patches_indices_global(i_all)) = 0.
     gr_HeII_density(gr_patches_indices_global(i_all)) = 0.
     gr_HeIII_density(gr_patches_indices_global(i_all)) = 0.
     ! metal density : treat separately ==> density = chem species + metals + dust
     ! Densities in the fraction of ISM/JET in fraction of rho_tot = mean_mass * n_H
     gr_metal_density(gr_patches_indices_global(i_all)) = 0.
     ! dust density : treat separately ==> density = chem species + metals + dust
     ! Densities in the fraction of ISM/JET in fraction of rho_tot = mean_mass * n_H
     gr_dust_density(gr_patches_indices_global(i_all)) = 0.

     gr_volumetric_heating_rate(gr_patches_indices_global(i_all)) = 0.
     gr_specific_heating_rate(gr_patches_indices_global(i_all)) = 0.
     gr_RT_HI_ionization_rate(gr_patches_indices_global(i_all)) = 0.
     gr_RT_HeI_ionization_rate(gr_patches_indices_global(i_all)) = 0.
     gr_RT_HeII_ionization_rate(gr_patches_indices_global(i_all)) = 0.
     gr_RT_H2_dissociation_rate(gr_patches_indices_global(i_all)) = 0.
     gr_RT_heating_rate(gr_patches_indices_global(i_all)) = 0.

   end do Loop_fill_everything

   write(*,*) 'Grackle fields configuration defaulting successfully done!'

 end subroutine grackle_fields_config_set_default

 subroutine grackle_fields_config_read(self,files)
  implicit none
  class(gr_objects),TARGET                         :: self
  character(len=*),intent(in)                      :: files(:)
  ! .. local ..
  integer                            :: i_file,i_error_read,i_ism,i_jet
  integer                  :: idim,iside,n_total_objects
  character(len=70)            :: error_message
  !--------------------------------------------------------------------
  namelist /grackle_fields_par_list/gr_patches_name,&
  gr_profiles,&
  gr_epsilon_tol,&
  gr_x_velocity,&
  gr_y_velocity,&
  gr_z_velocity,&
  gr_density_method,&
  gr_HI_density,&
  gr_HII_density,&
  gr_HM_density,&
  gr_HeI_density,&
  gr_HeII_density,&
  gr_HeIII_density,&
  gr_H2I_density,&
  gr_H2II_density,&
  gr_DI_density,&
  gr_DII_density,&
  gr_HDI_density,&
  gr_e_density,&
  gr_metal_density,&
  gr_dust_density,&
  gr_volumetric_heating_rate,&
  gr_specific_heating_rate,&
  gr_RT_HI_ionization_rate,&
  gr_RT_HeI_ionization_rate,&
  gr_RT_HeII_ionization_rate,&
  gr_RT_H2_dissociation_rate,&
  gr_RT_heating_rate

  Loop_read_fields_par : do i_file = 1, size(files)
    open(unitpar, file=trim(files(i_file)), status="old")
    read(unitpar, grackle_fields_par_list, iostat=i_error_read)
    cond_ierror : if(i_error_read>0)then
     write(*,*)' Error in reading the parameters file : ',trim(files(i_file))
     write(*,*)' Error at namelist: ', 'grackle_fields_par_list'
     write(*,*)' Error number = ',i_error_read
     write(*,*)' The code stops now '
     call mpistop(trim(error_message))
    elseif(i_error_read<0)then cond_ierror
     write(*,*)' Reache the end of the file  : ',trim(files(i_file))
     write(*,*)' Error at namelist: grackle_fields_par_list'
     write(*,*)' Error number = ',i_error_read
     write(*,*)' The code stops now '
     call mpistop(trim(error_message))
    else cond_ierror
     write(*,*)' End of reading of the grackle_fields_par_list'
    end if cond_ierror
    close(unitpar)

  end do Loop_read_fields_par


  write(*,*) 'Grackle fields configuration reading successfully done!'

 end subroutine grackle_fields_config_read

 !--------------------------------------------------------------------
 !> subroutine check the parfile setting for ism
 subroutine grackle_set_complet(self)
   implicit none
   class(gr_objects), TARGET                                :: self
   ! .. local ..
   integer :: i_all
   !-----------------------------------

   ! TO DO: nullify and ajust the densities according to which components
   ! (composition,coolings,rates,etc...) is not enabled

   Loop_through_all : do i_all = 1,size(gr_density_method)

      select case(gr_density_method(gr_patches_indices_global(i_all)))
      !TO DO: other cases than chemical_coefficient

      case('chemical_coefficient')
        !DO NOTHING TO COMPLETE
      case default
        !DO NOTHING TO COMPLETE
      end select

   end do Loop_through_all

 end subroutine grackle_set_complet

 !--------------------------------------------------------------------
 !> subroutine normalize setting for Chemistry configuration
  subroutine grackle_normalize(self,gr_struct,physunit_inuse)
   use mod_obj_usr_unit
   implicit none
   class(gr_objects), TARGET                              :: self
   class(grackle_type)                                    :: gr_struct
   type(usrphysical_unit), target,intent(in)      :: physunit_inuse
   !----------------------------------
   self%myphysunit => physunit_inuse

   if(trim(self%myconfig%unit)=='code')then
      return
   end if

   if(SI_unit) then
     gr_struct%units%density_units = self%myphysunit%myconfig%density*0.001
     gr_struct%units%length_units = self%myphysunit%myconfig%length*0.01
   else
     gr_struct%units%density_units = self%myphysunit%myconfig%density
     gr_struct%units%length_units = self%myphysunit%myconfig%length
   end if
   gr_struct%units%time_units = self%myphysunit%myconfig%time

   if(self%myconfig%gr_comoving_coordinates==0)then
    gr_struct%units%a_units = 1.0d0
    self%myconfig%gr_current_redshift = 0.
    !     Set initial expansion factor (for internal units)
    !     Set expansion factor to 1 for non-cosmological simulation.
    gr_struct%units%a_value = 1. / (1. + self%myconfig%gr_current_redshift)
   elseif(self%myconfig%gr_comoving_coordinates==1)then
    if(DABS(self%myconfig%gr_a_units)<smalldouble)then
      if(SI_unit) then
        gr_struct%units%a_units = self%myphysunit%myconfig%length*0.01/(self%myphysunit%myconfig%time*&
        self%myphysunit%myconfig%time)
      else
        gr_struct%units%a_units = self%myphysunit%myconfig%length/(self%myphysunit%myconfig%time*&
        self%myphysunit%myconfig%time)
      end if
    else
      gr_struct%units%a_units = self%myconfig%gr_a_units
    end if
    !     Set initial expansion factor (for internal units) for cosmological simulation
    gr_struct%units%a_value = 1. / (1. + self%myconfig%gr_current_redshift)
   end if
   call set_velocity_units(gr_struct%units)

  end subroutine grackle_normalize



  subroutine grackle_write_setting(self,unit_config)
    implicit none
    class(gr_objects),TARGET            :: self
    integer,intent(in)                  :: unit_config
    integer                             :: idims2,iside2,iB2
    real(kind=dp)                       :: rto_print
    character(len=64)                   :: sto_print
    character(len=128)                   :: wto_print
    integer                             :: idim,iside,idims,iw2,i_all
    ! .. local ..

    !-----------------------------------

    write(unit_config,*)'************************************'
    write(unit_config,*)'************Grackle settings ************'
    write(unit_config,*)'************************************'

    write(unit_config,*)'      ****** Grackle general parameters  *******      '
    write(unit_config,*) 'use_grackle =',  self%myconfig%use_grackle
    if(self%myconfig%use_grackle==1)then
      write(unit_config,*) 'Unit :',  self%myconfig%unit
      write(unit_config,*) 'gr_with_radiative_cooling = ',  self%myconfig%gr_with_radiative_cooling
      write(unit_config,*) 'gr_primordial_chemistry = ',  self%myconfig%gr_primordial_chemistry
      write(unit_config,*) 'gr_metal_cooling = ',  self%myconfig%gr_metal_cooling
      write(unit_config,*) 'gr_UVbackground = ',  self%myconfig%gr_UVbackground
      write(unit_config,*) 'gr_cmb_temperature_floor = ',  self%myconfig%gr_cmb_temperature_floor
      write(unit_config,*) 'gr_h2_on_dust = ',  self%myconfig%gr_h2_on_dust
      write(unit_config,*) 'gr_dust_chemistry = ',  self%myconfig%gr_dust_chemistry
      write(unit_config,*) 'gr_use_dust_density_field = ',  self%myconfig%gr_use_dust_density_field
      write(unit_config,*) 'gr_photoelectric_heating = ',  self%myconfig%gr_photoelectric_heating
      write(unit_config,*) 'gr_use_volumetric_heating_rate = ',  self%myconfig%gr_use_volumetric_heating_rate
      write(unit_config,*) 'gr_use_specific_heating_rate = ',  self%myconfig%gr_use_specific_heating_rate
      write(unit_config,*) 'gr_three_body_rate = ',  self%myconfig%gr_three_body_rate
      write(unit_config,*) 'gr_cie_cooling = ',  self%myconfig%gr_cie_cooling
      write(unit_config,*) 'gr_h2_optical_depth_approximation = ',  self%myconfig%gr_h2_optical_depth_approximation
      write(unit_config,*) 'gr_ih2co = ',  self%myconfig%gr_ih2co
      write(unit_config,*) 'gr_ipiht = ',  self%myconfig%gr_ipiht
      write(unit_config,*) 'gr_NumberOfTemperatureBins = ',  self%myconfig%gr_NumberOfTemperatureBins
      write(unit_config,*) 'gr_CaseBRecombination = ',  self%myconfig%gr_CaseBRecombination
      write(unit_config,*) 'gr_Compton_xray_heating = ',  self%myconfig%gr_Compton_xray_heating
      write(unit_config,*) 'gr_LWbackground_sawtooth_suppression = ',  self%myconfig%gr_LWbackground_sawtooth_suppression
      write(unit_config,*) 'gr_NumberOfDustTemperatureBins = ',  self%myconfig%gr_NumberOfDustTemperatureBins
      write(unit_config,*) 'gr_use_radiative_transfer = ',  self%myconfig%gr_use_radiative_transfer
      write(unit_config,*) 'gr_radiative_transfer_coupled_rate_solver = ',  self%myconfig%gr_radiative_transfer_coupled_rate_solver
      write(unit_config,*) 'gr_radiative_transfer_intermediate_step = ',  self%myconfig%gr_radiative_transfer_intermediate_step
      write(unit_config,*) 'gr_radiative_transfer_hydrogen_only = ',  self%myconfig%gr_radiative_transfer_hydrogen_only
      write(unit_config,*) 'gr_self_shielding_method = ',  self%myconfig%gr_self_shielding_method
      write(unit_config,*) 'gr_Gamma = ',  self%myconfig%gr_Gamma
      write(unit_config,*) 'gr_photoelectric_heating_rate = ',  self%myconfig%gr_photoelectric_heating_rate
      write(unit_config,*) 'gr_HydrogenFractionByMass = ',  self%myconfig%gr_HydrogenFractionByMass
      write(unit_config,*) 'gr_DeuteriumToHydrogenRatio = ',  self%myconfig%gr_DeuteriumToHydrogenRatio
      write(unit_config,*) 'gr_SolarMetalFractionByMass = ',  self%myconfig%gr_SolarMetalFractionByMass
      write(unit_config,*) 'TemperatureStart = ',  self%myconfig%TemperatureStart
      write(unit_config,*) 'TemperatureEnd = ',  self%myconfig%TemperatureEnd
      write(unit_config,*) 'DustTemperatureStart = ',  self%myconfig%DustTemperatureStart
      write(unit_config,*) 'DustTemperatureEnd = ',  self%myconfig%DustTemperatureEnd
      write(unit_config,*) 'gr_LWbackground_intensity = ',  self%myconfig%gr_LWbackground_intensity
      write(unit_config,*) 'gr_UVbackground_redshift_on = ',  self%myconfig%gr_UVbackground_redshift_on
      write(unit_config,*) 'gr_UVbackground_redshift_off = ',  self%myconfig%gr_UVbackground_redshift_off
      write(unit_config,*) 'gr_UVbackground_redshift_fullon = ',  self%myconfig%gr_UVbackground_redshift_fullon
      write(unit_config,*) 'gr_UVbackground_redshift_drop = ',  self%myconfig%gr_UVbackground_redshift_drop
      write(unit_config,*) 'cloudy_electron_fraction_factor = ',  self%myconfig%cloudy_electron_fraction_factor
      write(unit_config,*) 'data_dir = ',  self%myconfig%data_dir
      write(unit_config,*) 'data_filename = ',  self%myconfig%data_filename
      write(unit_config,*) 'data_file = ',  trim(self%myconfig%data_dir//self%myconfig%data_filename)
      write(unit_config,*) 'normalize_done = ',  self%myconfig%normalize_done
      write(unit_config,*) 'gr_comoving_coordinates = ',  self%myconfig%gr_comoving_coordinates
      write(unit_config,*) 'gr_a_units = ',  self%myconfig%gr_a_units
      write(unit_config,*) 'gr_current_redshift = ',  self%myconfig%gr_current_redshift

      Loop_fill_everything : do i_all=1,size(gr_patches_indices_global)
        write(unit_config,*)'      ****** Parameters for ',gr_patches_name(i_all),&
        '(',gr_patches_indices_local(i_all),') ******      '
        write(unit_config,*)'      ****** Physical Unit *******   '
        ! fields config
        write(unit_config,*) 'gr_epsilon_tol = ',gr_epsilon_tol(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_x_velocity = ',gr_x_velocity(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_y_velocity = ',gr_y_velocity(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_z_velocity = ',gr_z_velocity(gr_patches_indices_global(i_all))

        write(unit_config,*) 'gr_profiles = ',gr_profiles(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_density_method = ',gr_density_method(gr_patches_indices_global(i_all))
        ! Densities
        write(unit_config,*) 'gr_HI_density = ',gr_HI_density(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_HII_density = ',gr_HII_density(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_HM_density = ',gr_HM_density(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_H2I_density = ',gr_H2I_density(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_H2II_density = ',gr_H2II_density(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_DI_density = ',gr_DI_density(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_DII_density = ',gr_DII_density(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_HDI_density = ',gr_HDI_density(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_e_density = ',gr_e_density(gr_patches_indices_global(i_all))
        ! Densities as fraction of n_H=n(H+H^+)+2n(H_2)
        write(unit_config,*) 'gr_HeI_density = ',gr_HeI_density(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_HeII_density = ',gr_HeII_density(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_HeIII_density = ',gr_HeIII_density(gr_patches_indices_global(i_all))
        ! metal density : treat separately ==> density = chem species + metals + dust
        ! Densities in the fraction of ISM/JET in fraction of rho_tot = mean_mass * n_H
        write(unit_config,*) 'gr_metal_density = ',gr_metal_density(gr_patches_indices_global(i_all))
        ! dust density : treat separately ==> density = chem species + metals + dust
        ! Densities in the fraction of ISM/JET in fraction of rho_tot = mean_mass * n_H
        write(unit_config,*) 'gr_dust_density = ',gr_dust_density(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_volumetric_heating_rate = ',gr_volumetric_heating_rate(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_specific_heating_rate = ',gr_specific_heating_rate(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_RT_HI_ionization_rate = ',gr_RT_HI_ionization_rate(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_RT_HeI_ionization_rate = ',gr_RT_HeI_ionization_rate(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_RT_HeII_ionization_rate = ',gr_RT_HeII_ionization_rate(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_RT_H2_dissociation_rate = ',gr_RT_H2_dissociation_rate(gr_patches_indices_global(i_all))
        write(unit_config,*) 'gr_RT_heating_rate = ',gr_RT_heating_rate(gr_patches_indices_global(i_all))

      end do Loop_fill_everything

      write(unit_config,*) '======================================================='
      write(unit_config,*)'************************************'
      write(unit_config,*)'******** END Grackle setting **********'
      write(unit_config,*)'************************************'
    end if
    write(*,*)'************************************'
    write(*,*)'Finished writing grackle settings'
    write(*,*)'************************************'

  end    subroutine grackle_write_setting


























subroutine use_chemistry(parentname,self)
  implicit none
  character(len=*), intent(in) :: parentname

  !=================================================================================
  character(len=*),parameter :: subname='use_chemistry=>test_chemistry'
  integer :: ifield, irank^D
  integer :: iresult
  ! Define constants
  ! Grackle Fortran-to-C targets for field grackle data arrays
  class(gr_objects), TARGET :: self ! TARGET is very critically necessary !!
  ! Grackle C-to-Fortran pointers for field grackle data arrays
  type(grackle_type) gr_struct
  !=================================================================================

  write(*,*) '=================================================='
  call test_grackle_header()

  self%myfields%field_size(1)=3
  self%myfields%field_size(2)=5
  call self%reset_grid()

  ! set the length of the fields arrays along each dimension
  self%myfields%field_size(1)=2
  self%myfields%field_size(2)=3
  ! make sure to reset all fields arrays according to self%myfields%field_size attribute
  call self%reset_grid()

  !=================================================================================

  !     Create a grackle chemistry object for parameters and set defaults
  !     , set units, initialize the Grackle and set default field arrays

  call self%default_all_one(gr_struct)
  ! ==================================== !
  !> HERE: can add user-defined parameters and units

  ! ==================================== !
  call self%init_chem(gr_struct)

  call self%default_all_two(gr_struct)
  ! ==================================== !
  !> HERE: can add user-defined field arrays


  ! ==================================== !

  write(6,*) "primordial_chemistry:", &
  gr_struct%grackle_data%primordial_chemistry
  write(6,*) "metal_cooling:", &
  gr_struct%grackle_data%metal_cooling


  !=================================================================================
  !
  !     Fill in structure to be passed to Grackle
  !

  call self%link_to_grackle(gr_struct)

  !=================================================================================
  !
  !     Calling the chemistry solver
  !     These routines can now be called during the simulation.

  !     Evolving the chemistry.

  !write(*,*) 'Before solver >>> HI_density = '
  !write(*,*) 'global = ', self%myfields%gr_HI_density({1:self%myfields%field_size(^D)|,}) !<-defines the value of gr_HI_density
  !write(*,*) "local = ", gr_struct%fields%HI_density !<-defines the C-address to gr_HI_density
  !write(*,*) 'Before solver >>> internal energy = '
  !write(*,*) 'global = ', self%myfields%gr_energy({1:self%myfields%field_size(^D)|,}) !<-defines the value of gr_HI_density
  !write(*,*) "local = ", gr_struct%fields%internal_energy !<-defines the C-address to gr_HI_density
  self%myparams%dtchem = 3.15e7 * 1e6 / gr_struct%units%time_units    ! some timestep
  iresult = solve_chemistry(gr_struct%units, gr_struct%fields, self%myparams%dtchem)
  !write(*,*) 'After solver >>> HI_density = '
  !write(*,*) 'global = ', self%myfields%gr_HI_density({1:self%myfields%field_size(^D)|,}) !<-defines the value of gr_HI_density
  !write(*,*) "local = ", gr_struct%fields%HI_density !<-defines the C-address to gr_HI_density
  !write(*,*) 'After solver >>> internal energy = '
  !write(*,*) 'global = ', self%myfields%gr_energy({1:self%myfields%field_size(^D)|,}) !<-defines the value of gr_HI_density
  !write(*,*) "local = ", gr_struct%fields%internal_energy !<-defines the C-address to gr_HI_density

  !     Calculate cooling time.

  iresult = calculate_cooling_time(gr_struct%units, gr_struct%fields, &
    self%myfields%gr_cooling_time)
  write(*,*) "Cooling time = ", (self%myfields%gr_cooling_time({1:self%myfields%field_size(^D)|,}) * &
       gr_struct%units%time_units), "s."

  !     Calculate temperature.

  iresult = calculate_temperature(gr_struct%units, gr_struct%fields, self%myfields%gr_temperature)
  write(*,*) "Temperature = ", self%myfields%gr_temperature({1:self%myfields%field_size(^D)|,}), "K."

  !     Calcualte pressure.

  self%myparams%pressure_units = gr_struct%units%density_units * &
       gr_struct%units%velocity_units**2
  iresult = calculate_pressure(gr_struct%units, gr_struct%fields, self%myfields%gr_pressure)
  write(*,*) "Pressure = ", self%myfields%gr_pressure({1:self%myfields%field_size(^D)|,})* &
  self%myparams%pressure_units, "dyne/cm^2."

  !     Calculate gamma.

  iresult = calculate_gamma(gr_struct%units, gr_struct%fields, self%myfields%gr_gamma)
  write(*,*) "Gamma = ", self%myfields%gr_gamma({1:self%myfields%field_size(^D)|,})

  !     Calculate dust temperature.

  iresult = calculate_dust_temperature(gr_struct%units, gr_struct%fields, &
       self%myfields%gr_dust_temperature)
  write(*,*) "Dust temperature = ", self%myfields%gr_dust_temperature({1:self%myfields%field_size(^D)|,}), "K."

  ! Indexes inconsistency routine check with correct cases :
  ! call self%check_indexes(0,0,0,1,1,1,1,1,subname,mod_grackle_chemistry_name,parentname)
  ! Indexes inconsistency routine check with wrong cases :
  ! call self%check_indexes(0,0,0,1,1,1,1,0,subname,mod_grackle_chemistry_name,parentname)
  !call self%check_indexes({0,0,0,-1^D&|,},subname,mod_grackle_chemistry_name,parentname)

  write(*,*) '=================================================='

end subroutine use_chemistry

end module mod_grackle_chemistry
