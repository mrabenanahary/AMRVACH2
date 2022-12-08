module mod_obj_cla_jet
  use mod_constants
  use mod_global_parameters
  use mod_physics
  use mod_obj_dust, only : dust
  use mod_obj_global_parameters
  use mod_obj_mat
  use mod_obj_usr_unit
  implicit none
  type cla_jet_parameters
    character(len=20)    :: unit                   !> physical unit at parameter file
    logical              :: normalize_done         !> cla_jet is the normalisation is already done
    real(dp)             :: time_cla_jet_on        !> initial time the cla_jet is set in the simulation box

    integer              :: myindice               !> cla_jet associated indice
    real(dp)             :: density                !> cla_jet density  (g/cm^3)
    real(dp)             :: number_density         !> cla_jet number density (1/cm^3)
    real(dp)             :: temperature            !> cla_jet temperature  (K)
    real(dp)             :: pressure               !> cla_jet pressure
    real(dp)             :: pressure_toism         !> cla_jet pressure relative to ism pressure
    real(dp)             :: pressure_associate_ism !> cla_jet pressure of associated ism pressure
    real(dp)             :: density_toism         !> cla_jet density relative to ism pressure
    real(dp)             :: density_associate_ism !> cla_jet density of associated ism pressure
    real(dp)             :: magnetic(1:3)          !> cla_jet magnetic field components
    real(dp)             :: Magnetic_poloidal
    real(dp)             :: Magnetic_toroidal
    real(dp)             :: xisigma                !> cla_jet magnetisation
    real(dp)             :: magn_anglePHItoPol     !> cla_jet magnetisation field incl

    character(len=30)    :: chemical_gas_type
    real(kind=dp)        :: He_abundance
    real(kind=dp)        :: mean_mass
    real(kind=dp)        :: mean_mup
    real(kind=dp)        :: mean_ne_to_nH
    real(kind=dp)        :: mean_nall_to_nH

    real(dp)             :: mach_number            !> cla_jet  mach number
    real(dp)             :: c_sound                !> cla_jet  sound speed (cm/s)
    real(dp)             :: velocity(1:3)          !> cla_jet  velocity (cm/s)
    real(dp)             :: velocity_toroidal      !> cla_jet  toroidal velocity (cm/s)
    real(dp)             :: velocity_poloidal      !> cla_jet  poloidal velocity (cm/s)
    real(dp)             :: power                  !> wind power flux
    real(dp)             :: mass_flux              !> wind mass flux


    real(dp)             :: open_angle             !> cla_jet initial open angle  (degre)
    real(dp)             :: z_in                   !> cla_jet inner boundary position
    real(dp)             :: z_impos                !> cla_jet  impose r
    real(dp)             :: z_out_init             !> cla_jetinitial wind region
    real(dp)             :: r_out_init             !> cla_jet inner boundary wind position
    real(dp)             :: r_out_impos                !> cla_jet  impose radius
    real(dp)             :: r_in_init              !> initial wind region
    real(dp)             :: center(1:3)            !> cla_jet center position (cm)
    real(dp)             :: extend(1:3)            !> cla_jet region in space (cm)

    logical              :: tracer_on              !> cla_jet logical to set tracer
    integer              :: itr                    !> cla_jet tracer indice
    real(dp)             :: tracer_init_density    !> cla_jet tracer initial density
    real(dp)             :: tracer_small_density   !> cla_jet tracer small density cut
    character(len=20)    :: shape                  !> cla_jet shape
    character(len=20)    :: head_shape                  !> cla_jet head shape
    character(len=20)    :: profile                !> cla_jet profile

    integer              :: refine_min_level       !> jet minimum refinent level
    integer              :: refine_max_level       !> jet maximum refinent level
    real(dp)             :: coarsen_distance       !> jet  distance to start coarsen
    real(dp)             :: coarsen_var_distance   !> jetscaling distance for coarsen

    logical              :: variation_on           !>  jet variation  condition
    character(len=40)    :: variation_type         !> jet type of the power variation
    logical              :: density_variation_on   !>  jet density variation switch
    logical              :: pressure_variation_on  !>  jet pressure variation switch
    integer              :: variation_n_cells      !> jet ejecta variation cell number width
    real(dp)             :: variation_time         !> jet variation characteristic time
    integer              :: variation_phys_nvariable!> jet variation number of variable will change in time
    integer              :: variation_phys_variable(5)!> jet variation  variable will change in time
    real(dp)             :: variation_velocity(3)  !> jet variation bottom velocity
    real(dp)             :: variation_velocity_poloidal  !> jet variation bottom velocity  poloidal
    real(dp)             :: variation_velocity_toroidal  !> jet variation bottom velocity  toroidal
    real(dp)             :: variation_velocity_amplitude !> jet variation characteristic time

    real(dp)             :: variation_start_time  !> jet variation starting time
    real(dp)             :: variation_end_time    !> jet variation end  time
    real(dp)             :: variation_position(3) !> jet variation space position
    logical              :: dust_on                !> cla_jet with dust in is true
    real(dp)             :: dust_frac              !> cla_jet dust fraction
    character(len=20)    :: dust_profile           !> cla_jet dust inside profile
  end type cla_jet_parameters
  ! cla_jet features
  type cla_jet
    type(cla_jet_parameters)          :: myconfig           !> cla_jet configuration parameters
    type (dust)                       :: mydust             !> cla_jet dust
    type(usrphysical_unit), pointer   :: myphysunit         !> cla_jet physics unit in use
    logical, allocatable              :: patch(:^D&)        !> spatial patch
    logical, allocatable              :: ejecta_patch(:^D&)        !> spatial patch
    logical, allocatable              :: escape_patch(:^D&) !> spatial patch
    character(len=78)                 :: subname            !> subroutine name that call it

    contains
     !PRIVATE
     PROCEDURE, PASS(self) :: set_default     => usr_cla_jet_set_default
     PROCEDURE, PASS(self) :: set_complet     => usr_cla_jet_set_complet
     PROCEDURE, PASS(self) :: normalize       => usr_cla_jet_normalize
     PROCEDURE, PASS(self) :: set_w           => usr_cla_jet_set_w
     PROCEDURE, PASS(self) :: process_grid    => usr_cla_process_grid
     PROCEDURE, PASS(self) :: add_ejecta      => usr_cla_add_ejecta
     PROCEDURE, PASS(self) :: ejecta_set_patch=> usr_cla_ejecta_set_patch
     PROCEDURE, PASS(self) :: read_parameters => usr_cla_jet_read_p
     PROCEDURE, PASS(self) :: write_setting   => usr_cla_jet_write_setting
     PROCEDURE, PASS(self) :: alloc_set_patch => usr_cla_jet_alloc_set_patch
     PROCEDURE, PASS(self) :: set_patch       => usr_cla_jet_set_patch
     PROCEDURE, PASS(self) :: clean_memory    => usr_cla_jet_clean_memory
  end type
  integer, save            ::  zjet_,thetajet_,rjet_
contains

  !-------------------------------------------------------------------------
   !> Read the ism parameters  from a parfile
    subroutine usr_cla_jet_read_p(files,cla_jet_config,self)
      class(cla_jet)                          :: self
      character(len=*), intent(in)            :: files(:)
      type(cla_jet_parameters), intent(inout) :: cla_jet_config
      ! .. local ..
      integer                                 :: i_file
      !----------------------------------------------------------------
      namelist /usr_cla_jet_list/ cla_jet_config
      namelist /usr_cla_jet_1_list/ cla_jet_config
      namelist /usr_cla_jet_2_list/ cla_jet_config
      namelist /usr_cla_jet_3_list/ cla_jet_config



      select case(cla_jet_config%myindice)
      case(1)
       do i_file = 1, size(files)
          open(unitpar, file=trim(files(i_file)), status="old")
          read(unitpar, usr_cla_jet_1_list)
          close(unitpar)
       end do
      case(2)
       do i_file = 1, size(files)
          open(unitpar, file=trim(files(i_file)), status="old")
          read(unitpar, usr_cla_jet_2_list)
          close(unitpar)
       end do
      case(3)
       do i_file = 1, size(files)
          open(unitpar, file=trim(files(i_file)), status="old")
          read(unitpar, usr_cla_jet_3_list)
          close(unitpar)
       end do
      case default
       if(mype==0)write(*,*)'Reading usr_cla_jet_list'
       do i_file = 1, size(files)
          open(unitpar, file=trim(files(i_file)), status="old")
          read(unitpar, usr_cla_jet_list)
          close(unitpar)
       end do


     end select
     if(self%myconfig%dust_on)then
        self%mydust%myconfig%associated_medium = 'cloud'
        call self%mydust%read_parameters(self%mydust%myconfig,files)
    end if


    end subroutine usr_cla_jet_read_p

  !------------------------------------------------------------------------
  !> write the cla_jet setting
  subroutine usr_cla_jet_write_setting(self,unit_config)
    implicit none
    class(cla_jet)                        :: self
    integer,intent(in)                    :: unit_config
    ! .. local ..

    !-----------------------------------

    write(unit_config,*)'************************************'
    write(unit_config,*)'**********cla_jet setting **********'
    write(unit_config,*)'************************************'
    write(unit_config,*) 'power       = ', self%myconfig%power
    write(unit_config,*) 'mass flux   = ', self%myconfig%mass_flux
    write(unit_config,*) 'Density     = ', self%myconfig%density
    write(unit_config,*) 'Pressure    = ', self%myconfig%pressure
    write(unit_config,*) 'Temperature = ', self%myconfig%temperature
    write(unit_config,*) 'Speed       = ', self%myconfig%velocity
    write(unit_config,*) 'magnetic    = ', self%myconfig%magnetic
    write(unit_config,*) 'sound speed = ', self%myconfig%c_sound
    write(unit_config,*) 'open angle  = ', self%myconfig%open_angle
    write(unit_config,*) 'radius      = ', self%myconfig%r_out_impos
    write(unit_config,*) 'z_impos     = ', self%myconfig%z_impos
    write(unit_config,*) 'jet variation= ', self%myconfig%variation_on
    write(unit_config,*) 'mean_mass          = ', self%myconfig%mean_mass

    if(self%myconfig%variation_on)then
     write(unit_config,*) 'jet variation type = ', &
                              self%myconfig%variation_type
    write(unit_config,*) 'jet variation type = ', &
                              self%myconfig%variation_time
    write(unit_config,*) 'jet variation position = ', &
                              self%myconfig%variation_position

    write(unit_config,*) 'jet variation start time = ', &
                              self%myconfig%variation_start_time

    write(unit_config,*) 'jet variation end time = ', &
                              self%myconfig%variation_end_time

    write(unit_config,*) 'jet variation velocity = ', &
                              self%myconfig%variation_velocity

    write(unit_config,*) 'jet variation velocity amplitude = ', &
                              self%myconfig%variation_velocity_amplitude
    end if

    if(self%myconfig%dust_on)  call self%mydust%write_setting(unit_config)
    write(unit_config,*)'************************************'
    write(unit_config,*)'******* END cla_jet setting ********'
    write(unit_config,*)'************************************'
  end    subroutine usr_cla_jet_write_setting
  !--------------------------------------------------------------------
!> subroutine default setting for cla_jet
!M> sets the default values for the cla_jet class user parameters to ensure they all
!M> get an initialized value
 subroutine usr_cla_jet_set_default(self)
  implicit none
  class(cla_jet)          :: self
  !----------------------------------
  self%myconfig%unit                   = 'code'
  self%myconfig%myindice               = 0

  self%myconfig%density                = 0.0_dp
  self%myconfig%number_density         = 0.0_dp
  self%myconfig%temperature            = 0.0_dp
  self%myconfig%pressure               = 0.0_dp
  self%myconfig%pressure_toism         = 0.0_dp
  self%myconfig%pressure_associate_ism = 0.0_dp
  self%myconfig%density_toism          = 0.0_dp
  self%myconfig%density_associate_ism  = 0.0_dp
  self%myconfig%velocity(:)            = 0.0_dp
  self%myconfig%velocity_poloidal      = 0.0_dp
  self%myconfig%velocity_toroidal      = 0.0_dp
  self%myconfig%magnetic(:)            = 0.0_dp
  self%myconfig%magnetic_poloidal      = 0.0_dp
  self%myconfig%magnetic_toroidal      = 0.0_dp
  self%myconfig%mach_number            = 0.0_dp
  self%myconfig%c_sound                = 0.0_dp

  self%myconfig%power                  = 0.0_dp
  self%myconfig%mass_flux              = 0.0_dp

  self%myconfig%He_abundance           = 0.1_dp
  self%myconfig%chemical_gas_type      = 'fullyionised'
  call unit_chemical_ionisation( self%myconfig%He_abundance, &
                                      self%myconfig%chemical_gas_type,     &
                                      self%myconfig%mean_nall_to_nH, &
                                      self%myconfig%mean_mass,       &
                                      self%myconfig%mean_mup,        &
                                      self%myconfig%mean_ne_to_nH)




  self%myconfig%z_in                   = 0.0_dp
  self%myconfig%z_out_init             = 0.0_dp
  self%myconfig%z_impos                = 0.0_dp
  self%myconfig%r_out_init             = 0.0_dp
  self%myconfig%center(:)              = 0.0_dp
  self%myconfig%extend(:)              = 0.0_dp


  self%myconfig%shape                  = 'cylinder'
  self%myconfig%head_shape             = 'cartesian'
  self%myconfig%profile                = 'uniform'
  self%myconfig%open_angle             =  0.0_dp
  self%myconfig%tracer_on              = .false.
  self%myconfig%itr                    = 0
  self%myconfig%tracer_init_density    = 0.0_dp
  self%myconfig%tracer_small_density   = 0.0_dp

  self%myconfig%refine_min_level       = 0
  self%myconfig%refine_max_level       = 0
  self%myconfig%coarsen_var_distance   = 0.0_dp
  self%myconfig%coarsen_distance       = 0.0_dp



  self%myconfig%variation_on           = .false.
  self%myconfig%density_variation_on   = .false.
  self%myconfig%pressure_variation_on  = .false.
  self%myconfig%variation_n_cells      = 2
  self%myconfig%variation_type         = 'none'
  self%myconfig%variation_time         = 0.0_dp
  self%myconfig%variation_position     = 0.0_dp
  self%myconfig%variation_start_time   = 0.0_dp
  self%myconfig%variation_end_time     = 0.0_dp
  self%myconfig%variation_velocity     = 0.0_dp
  self%myconfig%variation_velocity_poloidal = 0.0_dp
  self%myconfig%variation_velocity_toroidal = 0.0_dp
  self%myconfig%variation_velocity_amplitude     = 0.0_dp
  self%myconfig%variation_phys_nvariable= 0
  self%myconfig%variation_phys_variable = 1
  self%myconfig%normalize_done         = .false.
  zjet_                                = 2
  thetajet_                            = 2
  rjet_                                = 1
  call self%mydust%set_default()
  !M> List of user-defined (then used/manipulated in this .t fortran module) jet parameters:
  !M>
  !M> self%myconfig%unit = physical unit of parameter file
  !M> self%myconfig%time_cla_jet_on = jet launching initial time inside the simulation box
  !M> self%myconfig%density = rho = jet total density [n(H+H^+)+2n(H_2)+n(He)]*m_H (g/cm^3)
  !M> self%myconfig%number_density = n_H = jet hydrogen-formed species number density [n(H+H^2)+2*n(H_2)] (1/cm^3)
  !M> self%myconfig%temperature = T = jet temperature (K)
  !M> self%myconfig%pressure = p = jet pressure (Ba)
  !M> self%myconfig%chemical_gas_type = f_gas = form of the simulated fluid (molecular, atomic, fully ionised, ionised) : determines the mean H-formed species mass factor mean_mass = rho/(m_H*n_H)
  !M>  and mean_mup = mean_mass/(n())
  !M> self%myconfig%He_abundance = x(He) = He abundance relative to the density of hydrogen-formed species proton density n_H=n(H+H^)+2*n(H_2) : x(He)=[n(He)*n_H*m_H]/[n_H*m_H]=n(He)/n_H
  !real(kind=dp)        :: mean_mass
  !real(kind=dp)        :: mean_mup
  !real(kind=dp)        :: mean_ne_to_nH
  !real(kind=dp)        :: mean_nall_to_nH

  !real(dp)             :: mach_number            !> cla_jet  mach number
  !real(dp)             :: c_sound                !> cla_jet  sound speed (cm/s)
  !real(dp)             :: velocity(1:3)          !> cla_jet  velocity (cm/s)
  !real(dp)             :: velocity_toroidal      !> cla_jet  toroidal velocity (cm/s)
  !real(dp)             :: velocity_poloidal      !> cla_jet  poloidal velocity (cm/s)
  !real(dp)             :: power                  !> wind power flux
  !real(dp)             :: mass_flux              !> wind mass flux


  !real(dp)             :: open_angle             !> cla_jet initial open angle  (degre)
  !real(dp)             :: z_in                   !> cla_jet inner boundary position
  !real(dp)             :: z_impos                !> cla_jet  impose r
  !real(dp)             :: r_out_init             !> cla_jet inner boundary wind position
  !real(dp)             :: z_out_init             !> cla_jetinitial wind region
  !real(dp)             :: r_out_impos                !> cla_jet  impose radius
  !real(dp)             :: r_in_init              !> initial wind region
  !real(dp)             :: center(1:3)            !> cla_jet center position (cm)
  !real(dp)             :: extend(1:3)            !> cla_jet region in space (cm)

  !logical              :: tracer_on              !> cla_jet logical to set tracer
  !integer              :: itr                    !> cla_jet tracer indice
  !real(dp)             :: tracer_init_density    !> cla_jet tracer initial density
  !real(dp)             :: tracer_small_density   !> cla_jet tracer small density cut
  !character(len=20)    :: shape                  !> cla_jet shape
  !character(len=20)    :: profile                !> cla_jet profile

  !integer              :: refine_min_level       !> jet minimum refinent level
  !integer              :: refine_max_level       !> jet maximum refinent level
  !real(dp)             :: coarsen_distance       !> jet  distance to start coarsen
  !real(dp)             :: coarsen_var_distance   !> jetscaling distance for coarsen

  !logical              :: variation_on           !>  jet variation  condition
  !character(len=40)    :: variation_type         !> jet type of the power variation
  !real(dp)             :: variation_time         !> jet variation characteristic time
  !integer              :: variation_phys_nvariable!> jet variation number of variable will change in time
  !integer              :: variation_phys_variable(5)!> jet variation  variable will change in time
  !real(dp)             :: variation_velocity(3)  !> jet variation bottom velocity
  !real(dp)             :: variation_velocity_poloidal  !> jet variation bottom velocity  poloidal
  !real(dp)             :: variation_velocity_toroidal  !> jet variation bottom velocity  toroidal
  !real(dp)             :: variation_velocity_amplitude !> jet variation characteristic time

  !real(dp)             :: variation_start_time  !> jet variation starting time
  !real(dp)             :: variation_end_time    !> jet variation end  time
  !real(dp)             :: variation_position(3) !> jet variation space position
  !logical              :: dust_on                !> cla_jet with dust in is true
  !real(dp)             :: dust_frac              !> cla_jet dust fraction
  !character(len=20)    :: dust_profile           !> cla_jet dust inside profile
 end subroutine usr_cla_jet_set_default
 !--------------------------------------------------------------------
 !> subroutine check the parfile setting for cla_jet
 subroutine usr_cla_jet_set_complet(self)
   implicit none
   class(cla_jet)             :: self
   ! .. local ..
   logical                    :: dust_is_frac
   real(dp)                   :: mp,kb, jet_surface_init,&
                                 Magnetic_poloidal
   real(kind=dp)              :: jet_radius_out_init, jet_radius_in_init
   real(kind=dp)              :: open_angle_in,z0, ztot
   integer                    :: idims
   !-----------------------------------
   zjet_     = min(zjet_,ndim)
   thetajet_ = zjet_

   if(SI_unit) then
     mp=mp_SI
     kB=kB_SI
   else
     mp=mp_cgs
     kB=kB_cgs
   end if

    call phys_fill_chemical_ionisation(self%myconfig%He_abundance,self%myconfig%chemical_gas_type, &
       self%myconfig%mean_nall_to_nH,self%myconfig%mean_mass,&
      self%myconfig%mean_mup,self%myconfig%mean_ne_to_nH)


   if(ndim>1) then
   self%myconfig%r_in_init   = nint((self%myconfig%r_in_init/unit_length-box_limit(1,r_)) &
                                  /dx(1,refine_max_level))  &
                               *dx(1,refine_max_level)*unit_length&
                               -0.5_dp*dx(1,refine_max_level)*unit_length
   print*,'1)the jet radius = ',self%myconfig%r_out_init
   self%myconfig%r_out_init  = nint((self%myconfig%r_out_init/unit_length-box_limit(1,r_)) &
                                    /dx(1,refine_max_level))  &
                               *dx(1,refine_max_level)*unit_length&
                               +0.5_dp*dx(1,refine_max_level)*unit_length
   print*,'2)the jet radius = ',self%myconfig%r_out_init

   self%myconfig%r_out_impos = nint((self%myconfig%r_out_impos/unit_length-box_limit(1,r_)) &
                                      /dx(1,refine_max_level))  &
                                 *dx(1,refine_max_level)*unit_length &
                                 +0.5_dp*dx(1,refine_max_level)*unit_length
  end if
  self%myconfig%z_in        = nint((self%myconfig%z_in/unit_length-box_limit(1,zjet_)) &
                                    /dx(zjet_,refine_max_level))  &
                               *dx(zjet_,refine_max_level)*unit_length&
                              +0.5_dp*dx(zjet_,refine_max_level)*unit_length
  self%myconfig%z_out_init  = nint((self%myconfig%z_out_init/unit_length-box_limit(1,zjet_)) &
                                    /dx(zjet_,refine_max_level))  &
                               *dx(zjet_,refine_max_level)*unit_length+0.5_dp*dx(zjet_,refine_max_level)*unit_length
  self%myconfig%z_impos     = nint((self%myconfig%z_impos/unit_length-box_limit(1,zjet_)) &
                                    /dx(zjet_,refine_max_level))   &
                               *dx(zjet_,refine_max_level)*unit_length&
                              +0.5_dp*dx(zjet_,refine_max_level)*unit_length



  Loop_dim : do idims = 1,ndim
   self%myconfig%variation_position(idims) = nint((self%myconfig%variation_position(idims)/unit_length-box_limit(1,idims)) &
                                      /dx(idims,refine_max_level))   &
                                 *dx(idims,refine_max_level)*unit_length&
                                +0.5_dp*dx(idims,refine_max_level)*unit_length
  end do  Loop_dim

   cond_vtor_set : if (self%myconfig%velocity_toroidal>0.0_dp) then
     self%myconfig%velocity(phi_)   = self%myconfig%velocity_toroidal
   else  cond_vtor_set
     self%myconfig%velocity_toroidal     = self%myconfig%velocity(phi_)
   end if cond_vtor_set

   cond_vpol_set : if (self%myconfig%velocity_poloidal>0.0_dp) then
      self%myconfig%velocity(r_)     = self%myconfig%velocity_poloidal*&
            dsin(self%myconfig%open_angle*dpi/180.0_dp)
      self%myconfig%velocity(zjet_) = self%myconfig%velocity_poloidal*&
             dcos(self%myconfig%open_angle*dpi/180.0_dp)
      self%myconfig%velocity(phi_)   = 0.0_dp

   else cond_vpol_set
      self%myconfig%velocity_poloidal=dsqrt(self%myconfig%velocity(r_)**2.0_dp&
                                           +self%myconfig%velocity(zjet_)**2.0_dp)
   end if cond_vpol_set

   ! if(self%myconfig%r_out_init>0.0_dp) then
   !   if(self%myconfig%open_angle>0.0_dp)then
   !    self%myconfig%z_in  = self%myconfig%r_out_init/dtan(self%myconfig%open_angle)
   !   else
   !    self%myconfig%z_in=min(box_limit(1,zjet_)/2.0_dp,self%myconfig%z_in)
   !   end if
   ! end if

print*,'the jet self%myconfig%r_out_init = ',self%myconfig%r_out_init
   print*, 'jet shape ', self%myconfig%shape
   cond_theta_O : if(self%myconfig%open_angle>0.0_dp)then
      self%myconfig%shape='conical'
      ztot                  = self%myconfig%r_out_init/dtan(self%myconfig%open_angle*dpi/180.0_dp)
      z0 = ztot - self%myconfig%z_out_init
      jet_radius_out_init = dsqrt(self%myconfig%r_out_init**2.0_dp&
                          +(z0+self%myconfig%z_out_init)**2.0_dp)
      jet_radius_in_init  = jet_radius_out_init
      if(self%myconfig%r_in_init>0.0_dp)then
        open_angle_in       = self%myconfig%open_angle*self%myconfig%r_out_init/self%myconfig%r_in_init
      else
        open_angle_in       =  0.0_dp
      end if
   else cond_theta_O
     jet_radius_out_init = self%myconfig%r_out_init
     jet_radius_in_init  = self%myconfig%r_in_init
   end if cond_theta_O
   print*, 'jet shape ', self%myconfig%shape
   print*, 'jet typeaxial ', typeaxial
   print*,'the jet Surface = ',jet_surface_init
   print*,'the jet open_angle = ',self%myconfig%open_angle
   print*,'the jet open_angle_in = ',open_angle_in
   print*,'the jet z0 = ', z0
   print*,'the jet jet_radius_out_init = ',jet_radius_out_init

   select case(typeaxial)
   case('cylindrical','spherical')
    cond_theta_1 : if(self%myconfig%open_angle>0.0_dp)then
      jet_surface_init = 2.0_dp*dpi *(jet_radius_out_init**2.0_dp*&
                             (dcos(open_angle_in*dpi/180.0_dp)- &
                             dcos(self%myconfig%open_angle*dpi/180.0_dp)))
    else  cond_theta_1
      jet_surface_init = dpi *(jet_radius_out_init**2.0_dp&
                             -max(self%myconfig%r_in_init,0.0_dp)**2.0_dp)
    end if cond_theta_1
   case('slab')
      if(self%myconfig%r_out_init*self%myconfig%r_in_init>0.0_dp&
       {^NOONED.or.self%myconfig%r_in_init>xprobmin1}) then
        jet_surface_init = dpi *(self%myconfig%r_out_init**2.0_dp&
                             -self%myconfig%r_in_init**2.0_dp)
      else
        jet_surface_init = dpi *(self%myconfig%r_out_init-self%myconfig%r_in_init)**2.0_dp/2.0_dp
      end if
   end select
   print*, 'jet typeaxial ', typeaxial
   print*,'the jet Surface = ',jet_surface_init









   if(dabs(self%myconfig%power)>smalldouble) then
     self%myconfig%mass_flux           = self%myconfig%power&
                                      /self%myconfig%velocity_poloidal**2.0_dp
   end if


  cond_massflux : if (dabs(self%myconfig%mass_flux)>smalldouble)then
    self%myconfig%density              = self%myconfig%mass_flux/(jet_surface_init*  &
                                         dabs(self%myconfig%velocity_poloidal))

    self%myconfig%number_density = self%myconfig%density/(mp*phys_config%mean_mass)
    self%myconfig%power                   = self%myconfig%mass_flux*unit_velocity
  else cond_massflux
    if (dabs(self%myconfig%density)<smalldouble*mp)then
      self%myconfig%density               = self%myconfig%number_density*mp*phys_config%mean_mass
    else if (dabs(self%myconfig%number_density)<smalldouble)then
      self%myconfig%number_density        = self%myconfig%density/(mp*phys_config%mean_mass)
    end if
    self%myconfig%mass_flux      = self%myconfig%density*jet_surface_init *  &
                                    dabs(self%myconfig%velocity_poloidal)
    self%myconfig%power          = self%myconfig%mass_flux*self%myconfig%velocity_poloidal**2.0_dp
  end if cond_massflux




   if(self%myconfig%pressure_toism>0.0_dp.and.self%myconfig%pressure_associate_ism>0.0_dp) then
     self%myconfig%pressure    = self%myconfig%pressure_toism*self%myconfig%pressure_associate_ism
   end if

   if(self%myconfig%density_toism>0.0_dp.and.self%myconfig%density_associate_ism>0.0_dp) then
    self%myconfig%density    = self%myconfig%density_toism*self%myconfig%density_associate_ism
   end if
   cond_Mach_set : if(self%myconfig%mach_number>0.0_dp) then
     self%myconfig%c_sound = dsqrt(sum(self%myconfig%velocity**2.0_dp))/&
       self%myconfig%mach_number
   else if(self%myconfig%c_sound>0.0_dp) then
     self%myconfig%mach_number = dsqrt(sum(self%myconfig%velocity**2.0_dp))/self%myconfig%c_sound
   end if cond_Mach_set


   cond_csound_set : if(self%myconfig%c_sound>0.0_dp) then
      self%myconfig%pressure = self%myconfig%c_sound**2.0_dp * self%myconfig%density /&
                               phys_config%gamma
      self%myconfig%temperature = self%myconfig%mean_mup*self%myconfig%pressure &
       /(kB*self%myconfig%number_density*self%myconfig%mean_mass)
   else  cond_csound_set

    if(dabs(self%myconfig%pressure)<=0.0_dp) then
    self%myconfig%pressure = (self%myconfig%number_density*self%myconfig%mean_mass/self%myconfig%mean_mup)*&
                            kB* self%myconfig%temperature
    else
     self%myconfig%temperature = self%myconfig%mean_mup*self%myconfig%pressure &
       /(kB*self%myconfig%number_density*self%myconfig%mean_mass)
    end if
    self%myconfig%c_sound = sqrt(phys_config%gamma*self%myconfig%pressure/self%myconfig%density)
   end if cond_csound_set




  if(self%myconfig%xisigma>0.0_dp)then
    if(dabs(dabs(dsin(self%myconfig%magn_anglePHItoPol*dpi/180.0_dp))-1.0_dp)>smalldouble)then
     self%myconfig%Magnetic_poloidal =  dsqrt( self%myconfig%density*&
                          self%myconfig%xisigma&
                         /(1.0_dp+(dtan(self%myconfig%magn_anglePHItoPol*dpi/180.0_dp))**2.0_dp))

     self%myconfig%Magnetic_toroidal = Magnetic_poloidal * dtan(self%myconfig%magn_anglePHItoPol*dpi/180.0d0)
     self%myconfig%magnetic(r_)   = Magnetic_poloidal * dsin(self%myconfig%open_angle*dpi/180.0d0)
     self%myconfig%magnetic(zjet_)   = Magnetic_poloidal * dcos(self%myconfig%open_angle*dpi/180.0d0)
     self%myconfig%magnetic(phi_) = Magnetic_poloidal * dtan(self%myconfig%magn_anglePHItoPol*dpi/180.0d0)
    else
     self%myconfig%magnetic(r_)    =  0.0_dp
     self%myconfig%magnetic(zjet_)    =  0.0_dp
     self%myconfig%magnetic(phi_)  = sign(dsqrt( self%myconfig%density*&
                          self%myconfig%xisigma),self%myconfig%magn_anglePHItoPol)
     self%myconfig%Magnetic_poloidal = 0.0_dp
    end if
  end if

   if(any(self%myconfig%extend(:)<=0.0_dp))then
     self%myconfig%extend(r_)                                 = self%myconfig%r_out_init
     if(z_in.or.zjet_<=ndim) self%myconfig%extend(thetajet_)  = self%myconfig%open_angle
     if(phi_in) self%myconfig%extend(phi_)                    = 360.0_dp
   end if





 if(self%myconfig%refine_min_level==0)then
    self%myconfig%refine_min_level=1
  end if
  if(self%myconfig%refine_max_level==0)then
    self%myconfig%refine_max_level=refine_max_level
  end if
  if(dabs(self%myconfig%coarsen_var_distance )<= 0.0_dp)then
   self%myconfig%coarsen_var_distance=HUGE(0.0_dp)
  end if
  if(dabs(self%myconfig%coarsen_distance )<= 0.0_dp)then
   self%myconfig%coarsen_distance=HUGE(0.0_dp)
  end if




   cond_duston : if(self%myconfig%dust_on)then
     dust_is_frac = .false.
     call self%mydust%set_complet(dust_is_frac,self%myconfig%dust_frac,self%myconfig%density&
               ,self%myconfig%velocity)
   end if cond_duston

   cond_traceron : if(self%myconfig%tracer_on)then
     prim_wnames(self%myconfig%itr+(phys_ind%tracer(1)-1)) = 'tracer_jet'
     cons_wnames(self%myconfig%itr+(phys_ind%tracer(1)-1)) = 'tracer_jet'
   else cond_traceron ! OFF NOW
     self%myconfig%itr=0
   end if cond_traceron



  cond_var0 : if(trim(self%myconfig%variation_type)=='none'.or.&
          dabs(self%myconfig%variation_time)<smalldouble)then
    self%myconfig%variation_on=.false.
    self%myconfig%density_variation_on =.false.
    self%myconfig%pressure_variation_on = .false.
  end if cond_var0
  cond_var1 : if(self%myconfig%variation_on)then
    self%myconfig%variation_phys_nvariable = max(self%myconfig%variation_phys_nvariable,1)
  end if cond_var1


   cond_vtor_var : if (self%myconfig%variation_velocity_toroidal>0.0_dp) then
     self%myconfig%variation_velocity(phi_)   = self%myconfig%variation_velocity_toroidal
   else  cond_vtor_var
     self%myconfig%variation_velocity_toroidal     = self%myconfig%variation_velocity(phi_)
   end if cond_vtor_var

   cond_vpol_var : if (self%myconfig%variation_velocity_poloidal>0.0_dp) then
      self%myconfig%variation_velocity(r_)     = self%myconfig%variation_velocity_poloidal*&
            dsin(self%myconfig%open_angle*dpi/180.0_dp)
      self%myconfig%variation_velocity(zjet_) = self%myconfig%variation_velocity_poloidal*&
             dcos(self%myconfig%open_angle*dpi/180.0_dp)
      self%myconfig%variation_velocity(phi_)   = 0.0_dp

   else cond_vpol_var
      self%myconfig%variation_velocity_poloidal=dsqrt(self%myconfig%variation_velocity(r_)**2.0_dp&
                                           +self%myconfig%variation_velocity(zjet_)**2.0_dp)
   end if cond_vpol_var
  if(ndim==1)then
    write(*,*)'With ndim=1 : The shape should be cartesian'
   self%myconfig%shape='cartesian'
  end if
  if(mype==0)then
    print*,'the jet temperature = ',self%myconfig%temperature
    print*, 'jet density', self%myconfig%density
    print*, 'jet ndensity', self%myconfig%number_density
    print*, 'set_complet jet shape', self%myconfig%shape
    print*, 'jet v_j', self%myconfig%velocity
    print*,'the jet mflux = ',self%myconfig%mass_flux
    print*,'3)the jet radius = ',self%myconfig%r_out_init
    print*,'3)the jet inner radius = ',self%myconfig%r_in_init
    print*,'the jet Surface = ',jet_surface_init

  end if

 end subroutine usr_cla_jet_set_complet
!--------------------------------------------------------------------
 subroutine usr_cla_jet_normalize(self,physunit_inuse)
  use mod_obj_usr_unit
  implicit none
  class(cla_jet)                                 :: self
  type(usrphysical_unit), target, intent(in)     :: physunit_inuse
  !----------------------------------

  self%myphysunit =>physunit_inuse
  if(trim(self%myconfig%unit)=='code'.or.self%myconfig%normalize_done)then
     if(self%myconfig%normalize_done)then
      write(*,*) 'WARNING: Second call for cloud normalisation', &
                   'no new normalisation will be done'
     end if
     return
  end if


  self%myconfig%density          = self%myconfig%density       /physunit_inuse%myconfig%density
  self%myconfig%number_density   = self%myconfig%number_density/physunit_inuse%myconfig%number_density
  self%myconfig%temperature      = self%myconfig%temperature   /physunit_inuse%myconfig%temperature
  self%myconfig%pressure         = self%myconfig%pressure      /physunit_inuse%myconfig%pressure
  self%myconfig%pressure_associate_ism  =self%myconfig%pressure_associate_ism/physunit_inuse%myconfig%pressure

  self%myconfig%mean_mup         =  self%myconfig%mean_mup      / physunit_inuse%myconfig%mean_mup

  self%myconfig%velocity         = self%myconfig%velocity /physunit_inuse%myconfig%velocity
  self%myconfig%velocity_poloidal= self%myconfig%velocity_poloidal/physunit_inuse%myconfig%velocity

  self%myconfig%velocity_toroidal= self%myconfig%velocity_toroidal/physunit_inuse%myconfig%velocity


  self%myconfig%magnetic         = self%myconfig%magnetic      /physunit_inuse%myconfig%magnetic
  self%myconfig%magnetic_poloidal= self%myconfig%magnetic_poloidal /physunit_inuse%myconfig%magnetic
  self%myconfig%magnetic_toroidal= self%myconfig%magnetic_toroidal /physunit_inuse%myconfig%magnetic
  self%myconfig%c_sound          =  self%myconfig%c_sound      /physunit_inuse%myconfig%velocity



  self%myconfig%open_angle       = self%myconfig%open_angle   *(dpi/180._dp)
  self%myconfig%magn_anglePHItoPol = self%myconfig%magn_anglePHItoPol *(dpi/180._dp)
  self%myconfig%center           = self%myconfig%center        /physunit_inuse%myconfig%length



  self%myconfig%z_in             = self%myconfig%z_in          /physunit_inuse%myconfig%length
  self%myconfig%z_impos          = self%myconfig%z_impos       /physunit_inuse%myconfig%length
  self%myconfig%z_out_init       = self%myconfig%z_out_init    /physunit_inuse%myconfig%length

  self%myconfig%r_out_init       = self%myconfig%r_out_init     /physunit_inuse%myconfig%length
  self%myconfig%r_out_impos      = self%myconfig%r_out_impos    /physunit_inuse%myconfig%length
  self%myconfig%r_in_init        = self%myconfig%r_in_init      /physunit_inuse%myconfig%length

  self%myconfig%extend(r_)       = self%myconfig%extend(r_)     /physunit_inuse%myconfig%length

  if(z_in)then
    select case(typeaxial)
    case('spherical')
     self%myconfig%extend(theta_) = self%myconfig%extend(theta_) *(dpi/180._dp)
    case('slab','cylindrical')
      self%myconfig%extend(zjet_) = self%myconfig%extend(zjet_)/physunit_inuse%myconfig%length
    end select
  end if
  if(phi_in)then
    select case(typeaxial)
    case('spherical','cylindrical')
       self%myconfig%extend(phi_) = self%myconfig%extend(phi_) *(dpi/180._dp)
    case('slab')
      self%myconfig%extend(phi_) = self%myconfig%extend(phi_)/physunit_inuse%myconfig%length
    end select
  end if

  self%myconfig%power            = self%myconfig%power            /unit_user%luminosity
  self%myconfig%mass_flux        = self%myconfig%mass_flux        /unit_user%mass_flux
  self%myconfig%time_cla_jet_on  = self%myconfig%time_cla_jet_on  / unit_time


  self%myconfig%coarsen_distance     = self%myconfig%coarsen_distance /physunit_inuse%myconfig%length
  self%myconfig%coarsen_var_distance = self%myconfig%coarsen_var_distance /physunit_inuse%myconfig%length

  self%myconfig%variation_time       = self%myconfig%variation_time/physunit_inuse%myconfig%time
  self%myconfig%variation_position   = self%myconfig%variation_position/physunit_inuse%myconfig%length
  self%myconfig%variation_start_time = self%myconfig%variation_start_time/physunit_inuse%myconfig%time
  self%myconfig%variation_end_time   = self%myconfig%variation_end_time/physunit_inuse%myconfig%time
  self%myconfig%variation_velocity   = self%myconfig%variation_velocity/physunit_inuse%myconfig%velocity
  self%myconfig%variation_velocity_poloidal = self%myconfig%variation_velocity_poloidal/physunit_inuse%myconfig%velocity
  self%myconfig%variation_velocity_toroidal = self%myconfig%variation_velocity_toroidal/physunit_inuse%myconfig%velocity

  self%myconfig%variation_velocity_amplitude   = self%myconfig%variation_velocity_amplitude/&
                                                 physunit_inuse%myconfig%velocity

  if(self%myconfig%dust_on)then
    call self%mydust%normalize(physunit_inuse)
    call self%mydust%to_phys()
  end if
  self%myconfig%normalize_done=.true.


 end subroutine usr_cla_jet_normalize


!--------------------------------------------------------------------

 !--------------------------------------------------------------------
 !> Subroutine to clean array memory of associated with cla_jet object
 subroutine usr_cla_jet_alloc_set_patch(ixI^L,ixO^L,qt,x,use_tracer,w,&
                                        escape_patch,self)
   implicit none
   integer, intent(in)                     :: ixI^L,ixO^L
   real(kind=dp), intent(in)               :: qt
   real(kind=dp), intent(in)               :: x(ixI^S,1:ndim)
   real(kind=dp), intent(in), optional     :: w(ixI^S,1:nw)
   logical, intent(in), optional           :: use_tracer
   logical, intent(in),optional            :: escape_patch(ixI^S)
   class(cla_jet)                          :: self
   ! .. local ..
   logical                                 :: usr_tracer_loc
   !---------------------------------------------------------
     if(present(use_tracer).and. .not.present(w)) then
       usr_tracer_loc = use_tracer
     else
       usr_tracer_loc = .false.
     end if
     cond_tracer : if(usr_tracer_loc) then
       if(allocated(self%patch))deallocate(self%patch)
       allocate(self%patch(ixI^S))
       self%patch(ixO^S) =.true.

     else cond_tracer
       cond_jettracer_on : if(self%myconfig%tracer_on)then
        if(allocated(self%patch))deallocate(self%patch)
        allocate(self%patch(ixI^S))
        where(w(ixO^S,phys_ind%tracer(self%myconfig%itr))>small_density)
          self%patch(ixO^S)=.true.
        else where
          self%patch(ixO^S)=.false.
        end where
      end if cond_jettracer_on

     end if cond_tracer

     if(allocated(self%escape_patch))deallocate(self%escape_patch)
     allocate(self%escape_patch(ixI^S))
     if(present(escape_patch))then
      self%escape_patch(ixO^S) = escape_patch(ixO^S)
     else
      self%escape_patch(ixO^S) =.false.
     end if
     self%patch(ixO^S)        = self%patch(ixO^S) .and. (.not.self%escape_patch(ixO^S))
 end subroutine usr_cla_jet_alloc_set_patch
!---------------------------------------------------------------------
!--------------------------------------------------------------------
 !> subroutine patch for the cla_jet
 subroutine usr_cla_jet_set_patch(ixI^L,ixO^L,qt,x,self,force_refine,dx_loc,&
                                  r_limit,z_out_force)
  implicit none
  integer, intent(in)       :: ixI^L,ixO^L
  real(kind=dp), intent(in) :: x(ixI^S,1:ndim)
  real(kind=dp), intent(in) :: qt
  class(cla_jet)            :: self
  integer, optional         :: force_refine
  real(kind=dp),optional    :: dx_loc(1:ndim)
  real(kind=dp),optional    :: r_limit(2)
  real(kind=dp),optional    :: z_out_force
  ! .. local ..
  integer                    :: idims
  real(dp)                   :: x_edge(ixI^S,1:ndim)
  real(dp), dimension(ixI^S) :: dist_edge,jet_radius, jet_h, jet_ortho_r
  real(dp)                   :: r_in,r_out
  real(dp)                   :: z_down,z_out,min_dist,ztot,z0,r0
  !----------------------------------
  ! check distance

  r_in =self%myconfig%r_in_init
  z_down =self%myconfig%z_in
  if(dabs(qt-self%myconfig%time_cla_jet_on)<=smalldouble)then
    r_out=self%myconfig%r_out_init
    z_out=self%myconfig%z_out_init
  else
    r_out=self%myconfig%r_out_impos
    z_out=self%myconfig%z_impos
  end if
  if(present(z_out_force)) then
    z_out = nint(z_out_force/dxlevel(zjet_))*dxlevel(zjet_)+0.5_dp*dxlevel(zjet_)
  end if
  if(present(r_limit))then
   r_limit(1) = r_in
   r_limit(2) = r_out
  end if

  cond_from_refine : if(present(force_refine)) then
    !if(force_refine==-1)r_in=r_out
   Loop_idim : do idims=1,ndim
    x_edge(ixO^S,idims) = x(ixO^S,idims) - dx_loc(idims)/2.0_dp
   end do Loop_idim
    Dist_edge(ixO^S) = x_edge(ixO^S,1)

   need_refine: if(minval(Dist_edge(ixO^S)-r_out)&
                   *maxval(x(ixO^S,1)-r_out)<=smalldouble.or. &
                   minval(Dist_edge(ixO^S)-r_in)&
                   *maxval(x(ixO^S,1)-r_in)<=smalldouble) then
     min_dist=minval(Dist_edge(ixO^S),Dist_edge(ixO^S)>smalldouble)
     r_out=r_out+min_dist
     r_in=r_in-min_dist
   end if  need_refine
  end if  cond_from_refine

  if(.not.allocated(self%patch))allocate(self%patch(ixI^S))

  self%patch              = .false.
  !print*, 'jet shape ', self%myconfig%shape
  !print*, 'jet typeaxial ', typeaxial
  !------------------------------------------------
  ! * Select the actual jet shape :
  ! * <!> : it is 'conical' IF user chose jet open_angle>0°
  !------------------------------------------------
  select case(trim(self%myconfig%shape))
  !------------------------------------------------
  ! * Case: conical jet shape:
  !------------------------------------------------
  case('conical')
    !------------------------------------------------
    ! ** Select the actual axial geometry frame used:
    !------------------------------------------------
    select case(typeaxial)
    !------------------------------------------------
    ! ** If spherical frame:
    !------------------------------------------------
    case('spherical')
      self%patch(ixO^S)  = dabs(x(ixO^S,theta_))<=self%myconfig%open_angle
      !------------------------------------------------
      ! *** If any dimension (since 1D==>zjet=1) :
      !------------------------------------------------
      if(zjet_<=ndim) then
        self%patch(ixO^S)  = self%patch(ixO^S).and. x(ixO^S,r_)<=z_out .and. &
                         x(ixO^S,r_)>=z_down

      end if
      !------------------------------------------------
      ! ***End if any dimension
      !------------------------------------------------
    !------------------------------------------------
    ! **If slab or cylindrical frame:
    !------------------------------------------------
    case('slab','cylindrical')
      !------------------------------------------------
      ! ***If we're >= 2D:
      !------------------------------------------------
      if(ndim>1)then
        !self%patch(ixO^S)  = dabs(x(ixO^S,rjet_))<=r_out&
        !          +dabs(x(ixO^S,zjet_)-box_limit(1,zjet_))&
        !          *dtan(self%myconfig%open_angle) .and. &
        !           dabs(x(ixO^S,r_))>=r_in&
        !          -dabs(x(ixO^S,zjet_)-box_limit(1,zjet_))&
        !          *dtan(self%myconfig%open_angle)
        ! Make the domaine in cylindrical frame for jet
        ! between r_out at z=z_impos and =r_out-self%myconfig%z_impos
        ! *dtan(self%myconfig%open_angle) at z=0
        !
        ! at any z and
        !
        ! r=r_in
        ! -dabs(self%myconfig%z_impos+
        ! at z=z_impos and =r_in at z=0
        !
        ! at any z
        jet_ortho_r(ixO^S) = dabs(z_out-x(ixO^S,zjet_))&
        *dtan(self%myconfig%open_angle)
        self%patch(ixO^S)  = dabs(x(ixO^S,rjet_))<=r_out&
                  -jet_ortho_r(ixO^S) .and. &
                   dabs(x(ixO^S,rjet_))>=r_in&
                  -jet_ortho_r(ixO^S)

        if(zjet_<=ndim) then
          !------------------------------------------------
          ! **** Select the jet head shape
          !------------------------------------------------
          select case(trim(self%myconfig%head_shape))
          !------------------------------------------------
          ! **** If you want the head to be spherical shaped:
          !------------------------------------------------
          case('spherical')
            ztot                = r_out/dtan(self%myconfig%open_angle) !
            z0                  = ztot-z_out !ztot = z_impos + z0
            !radius relative to (-z0,0) at cylindrical r=r_out_init:
            r0                = dsqrt((z_out+z0)**2.0_dp&
                                +r_out**2.0_dp)
            !radius relative to (-z0,0) for every cylindrical (r,z):
            jet_radius(ixO^S)      = dsqrt((x(ixO^S,zjet_)+z0)**2.0_dp&
                                    +x(ixO^S,rjet_)**2.0_dp)
            !for every r and z, select cylindrical r that are below r0
            self%patch(ixO^S)  = self%patch(ixO^S).and. jet_radius(ixO^S)<=r0&
                                .and. x(ixO^S,zjet_)>=z_down
          !------------------------------------------------
          ! **** Default for head is cartesian shaped:
          !------------------------------------------------
          case default
            self%patch(ixO^S)  = self%patch(ixO^S).and.  x(ixO^S,zjet_)<z_out .and. &
                                    x(ixO^S,zjet_)>=z_down
          end select
          !------------------------------------------------
          ! **** End select the jet head shape
          !------------------------------------------------
        end if
      !------------------------------------------------
      ! *** Else, if we're 1D:
      !------------------------------------------------
      else
        self%patch(ixO^S)  = self%patch(ixO^S).and. x(ixO^S,zjet_)<=z_out .and. &
                                x(ixO^S,zjet_)>=z_down
      end if
      !------------------------------------------------
      ! *** End if we're >=2D or if we're 1D
      !------------------------------------------------
    end select
    !------------------------------------------------
    ! ** End select the actual axial geometry frame used
    !------------------------------------------------
  !------------------------------------------------
  ! * End case conical jet shape
  !------------------------------------------------
  !------------------------------------------------
  ! * Case: cylindrical jet shape
  ! * (assuming cylindrical geometry frame):
  !------------------------------------------------
  case('cylindrical')
    self%patch(ixO^S)  = dabs(x(ixO^S,rjet_))<=r_out .and. &
                         dabs(x(ixO^S,rjet_))>=r_in
    !------------------------------------------------
    ! ** If we're >= 2D:
    !------------------------------------------------
    !if(zjet_<=ndim) then
    if(ndim>1) then
        !------------------------------------------------
        ! *** Select the jet head shape
        !------------------------------------------------
        select case(trim(self%myconfig%head_shape))
        !------------------------------------------------
        ! *** If you want the head to be spherical shaped:
        !------------------------------------------------
        case('spherical')
          !height z relative to (0,0) at any cylindrical r:
          jet_h(ixO^S)      = z_out&
          +dsqrt(dabs(r_out**2.0_dp-x(ixO^S,rjet_)**2.0_dp))
          !for every r and z, select z that are below jet_h
          self%patch(ixO^S)  = self%patch(ixO^S).and. x(ixO^S,zjet_)&
                              <=jet_h(ixO^S)&
                              .and. x(ixO^S,zjet_)>=z_down
        !------------------------------------------------
        ! *** Default for head is cartesian shaped:
        !------------------------------------------------
        case default
          self%patch(ixO^S)  = self%patch(ixO^S).and. x(ixO^S,zjet_)<=z_out .and. &
                           x(ixO^S,zjet_)>=z_down
        end select
        !------------------------------------------------
        ! *** End select the jet head shape
        !------------------------------------------------
        !self%patch(ixO^S)  = self%patch(ixO^S).and. x(ixO^S,zjet_)<=z_out .and. &
        !                  x(ixO^S,zjet_)>=z_down
    !------------------------------------------------
    ! ** Else, if we're 1D:
    !------------------------------------------------
    else
        self%patch(ixO^S)  = self%patch(ixO^S).and. x(ixO^S,zjet_)<=z_out .and. &
                          x(ixO^S,zjet_)>=z_down
    end if
    !------------------------------------------------
    ! ** End if we're >=2D or if we're 1D
    !------------------------------------------------
  !------------------------------------------------
  ! * End case cylindrical jet shape
  !------------------------------------------------
  !------------------------------------------------
  ! * Case: cartesian jet shape
  ! * (assuming cylindrical geometry frame):
  !------------------------------------------------
  case('cartesian')
    !------------------------------------------------
    ! ** If we're >= 2D:
    !------------------------------------------------
    if(ndim>1) then
     self%patch(ixO^S)  = x(ixO^S,rjet_)<=r_out .and. &
                         x(ixO^S,rjet_)>=r_in
    !------------------------------------------------
    ! ** Else, if we're 1D:
    !------------------------------------------------
    else
       self%patch(ixO^S)  = .true.
    end if
    !------------------------------------------------
    ! ** End if we're >= 2D or if we're 1D
    !------------------------------------------------


    if(zjet_<=ndim) then
      !------------------------------------------------
      ! *** Select the jet head shape
      !------------------------------------------------
      select case(trim(self%myconfig%head_shape))
      !------------------------------------------------
      ! *** If you want the head to be spherical shaped
      ! in 2D :
      !------------------------------------------------
      case('spherical')
        !------------------------------------------------
        ! **** If we're >= 2D:
        !------------------------------------------------
        if(ndim>1) then
          !height z relative to (0,0) at any cylindrical r:
          jet_h(ixO^S)      = z_out&
          +dsqrt(dabs(r_out**2.0_dp-x(ixO^S,rjet_)**2.0_dp))
          !for every r and z, select z that are below jet_h
          self%patch(ixO^S)  = self%patch(ixO^S).and. x(ixO^S,zjet_)&
                              <=jet_h(ixO^S)&
                              .and. x(ixO^S,zjet_)>=z_down
        !------------------------------------------------
        ! **** Else, if we're 1D:
        !------------------------------------------------
        else
          self%patch(ixO^S)  = self%patch(ixO^S).and. x(ixO^S,zjet_)<=z_out .and. &
                            x(ixO^S,zjet_)>=z_down
        end if
        !------------------------------------------------
        ! **** End if we're >= 2D or if we're 1D
        !------------------------------------------------
      !------------------------------------------------
      ! *** Default for head is cartesian shaped for any
      ! dimension :
      !------------------------------------------------
      case default
        self%patch(ixO^S)  = self%patch(ixO^S).and. x(ixO^S,zjet_)<=z_out .and. &
                          x(ixO^S,zjet_)>=z_down
      end select
      !------------------------------------------------
      ! *** End select the jet head shape
      !------------------------------------------------
    end if
  !------------------------------------------------
  ! * End case cartesian jet shape
  !------------------------------------------------
  !------------------------------------------------
  ! * Case: any other unknown jet shape is entered,
  ! * the code stops and raise an exception
  !------------------------------------------------
  case default
     write(*,*)'this cla_jet shape ',trim(self%myconfig%shape),' is not implimented'
     call mpistop('This cla_jet shape is not implimented in usr_cla_jet_patch at mod_usr.t')
  !------------------------------------------------
  ! * End case any other unknown jet shape
  !----------------------------------------------
  end select
  !------------------------------------------------
  ! * End select the actual jet shape
  !------------------------------------------------
end subroutine usr_cla_jet_set_patch
!-------------------------------------------------------------------
 !> subroutine setting for cla_jet
 subroutine usr_cla_jet_set_w(ixI^L,ixO^L,qt,x,w,self)
  implicit none
  integer, intent(in)          :: ixI^L,ixO^L
  real(kind=dp), intent(in)    :: qt
  real(kind=dp)                :: x(ixI^S,1:ndim)
  real(kind=dp)                :: w(ixI^S,1:nw)
  class(cla_jet)               :: self
  ! .. local..
  integer                         :: idir
  logical                         :: dust_is_frac
  real(kind=dp), dimension(ixI^S) :: fprofile,angle_theta,&
                                    jet_surface,&
                                    jet_radius_in,jet_radius
  real(kind=dp)                   :: z0,ztot,open_angle_in
  real(kind=dp), dimension(ixI^S,1:ndim)::project_speed
  !integer::ix2
  !----------------------------------
  if(.not.allocated(self%patch))then
    call self%set_patch(ixI^L,ixO^L,qt,x)
  end if


  cond_inside_cla_jet: if(any(self%patch(ixO^S))) then

  cond_conical : if(trim(self%myconfig%shape)=='conical')then
   select case(typeaxial)
    case('spherical')
      where(self%patch(ixO^S))
        jet_radius(ixO^S)  = x(ixO^S,r_)

        jet_surface(ixO^S) = 2.0_dp*dpi *x(ixO^S,r_)**2.0_dp*&
                            (1.0_dp- dcos(self%myconfig%open_angle))

        angle_theta(ixO^S) =  self%myconfig%open_angle
        project_speed(ixO^S,rjet_)  = dsin(self%myconfig%open_angle)
        project_speed(ixO^S,zjet_)  = dcos(self%myconfig%open_angle)

      end where
    case('slab','cylindrical')
      if(ndim>1)then
        open_angle_in       = self%myconfig%open_angle*self%myconfig%r_in_init/self%myconfig%r_out_init
        ztot                = self%myconfig%r_out_init/dtan(self%myconfig%open_angle)
        z0 = ztot-self%myconfig%z_impos
        where(self%patch(ixO^S))
          jet_radius(ixO^S) = dsqrt((z0 + x(ixO^S,zjet_))**2.0_dp+x(ixO^S,rjet_)**2.0_dp)


          jet_surface(ixO^S)    = 2.0_dp*dpi *(jet_radius(ixO^S)**2.0_dp*&
                                   (dcos(open_angle_in)- &
                                 dcos(self%myconfig%open_angle)))

          angle_theta(ixO^S)  =  datan(x(ixO^S,r_)/(z0+x(ixO^S,zjet_)))
          project_speed(ixO^S,rjet_)  = x(ixO^S,r_)/jet_radius(ixO^S)
          project_speed(ixO^S,zjet_)  = (z0 +x(ixO^S,z_))/jet_radius(ixO^S)
        end where
      else
        z0                          = 0.0_dp
        jet_surface(ixO^S)          = dpi *(self%myconfig%r_out_init**2.0_dp&
                                  -max(self%myconfig%r_in_init,0.0_dp)**2.0_dp)
        angle_theta(ixO^S)          =  0.0_dp
        project_speed(ixO^S,rjet_)  =0.0_dp
        project_speed(ixO^S,zjet_)  = 1.0_dp
      end if
    end select
  else  cond_conical
    jet_surface(ixO^S) =1.0_dp
    where(self%patch(ixO^S))
      jet_surface(ixO^S) = dpi *(self%myconfig%r_out_init**2.0_dp&
                            -max(self%myconfig%r_in_init,0.0_dp)**2.0_dp)

      angle_theta(ixO^S)  = 0.0_dp
      project_speed(ixO^S,rjet_) = 0.0_dp
      project_speed(ixO^S,zjet_) = 1.0_dp
    end where
  end if cond_conical




   where(self%patch(ixO^S))
    w(ixO^S,phys_ind%mom(rjet_)) = self%myconfig%velocity_poloidal &
                    * project_speed(ixO^S,rjet_)

    w(ixO^S,phys_ind%mom(zjet_)) = self%myconfig%velocity_poloidal &
                    * project_speed(ixO^S,zjet_)

    w(ixO^S,phys_ind%mom(phi_))  = self%myconfig%velocity(phi_)
   end where

   where(self%patch(ixO^S))
    w(ixO^S,phys_ind%rho_)        = self%myconfig%mass_flux/(jet_surface(ixO^S)*  &
         dabs(self%myconfig%velocity_poloidal))
  end where
  if(phys_config%energy)then
    where(self%patch(ixO^S))
      w(ixO^S,phys_ind%pressure_)   = self%myconfig%pressure
    end where
  end if
   cond_mhd0 : if(phys_config%ismhd)then
    where(self%patch(ixO^S))
     w(ixO^S,phys_ind%mag(r_))    = self%myconfig%Magnetic_poloidal*dsin(angle_theta(ixO^S))
     w(ixO^S,phys_ind%mag(zjet_)) = self%myconfig%Magnetic_poloidal*dcos(angle_theta(ixO^S))
     w(ixO^S,phys_ind%mag(phi_))  = self%myconfig%magnetic(phi_)
    end where
   end if cond_mhd0




  cond_profile : if(self%myconfig%profile/='none') then

      call usr_mat_profile(ixI^L,ixO^L,typeaxial,self%myconfig%profile,&
                           self%myconfig%center,self%myconfig%extend,&
                            x,fprofile)

      where(self%patch(ixO^S))
       w(ixO^S,phys_ind%rho_)       = w(ixO^S,phys_ind%rho_) * fprofile(ixO^S)
       w(ixO^S,phys_ind%mom(zjet_)) = w(ixO^S,phys_ind%mom(zjet_)) * fprofile(ixO^S)
      end where
      if(phys_config%energy)then
        where(self%patch(ixO^S))
          w(ixO^S,phys_ind%pressure_)  = w(ixO^S,phys_ind%pressure_) * fprofile(ixO^S)
        end where
      end if
  end if cond_profile


  cond_mhd : if(phys_config%ismhd.and.phys_config%energy)then
    where(self%patch(ixO^S))
     w(ixO^S,phys_ind%pressure_)   = w(ixO^S,phys_ind%pressure_) -&
                   ((1.0_dp-w(ixO^S,phys_ind%mom(zjet_))**2.0_dp)*&
                    0.5_dp*w(ixO^S,phys_ind%mag(phi_))**2.0_dp&
                    + w(ixO^S,phys_ind%mag(zjet_))**2.0_dp/2.0_dp)
    end where
  end if cond_mhd


   cond_var_0 : if(self%myconfig%variation_on.and.&
                   self%myconfig%variation_position(zjet_)<=self%myconfig%z_impos)then
    !  self%myconfig%variation_position(zjet_)=&
    !     (floor((self%myconfig%z_impos-box_limit(1,zjet_))/dxlevel(zjet_))+0.5_dp)&
    !     *dxlevel(zjet_)

     call self%add_ejecta(ixI^L,ixO^L,qt,.false.,x,w)

   end if cond_var_0


   if(phys_config%mean_mup_on) then
     where(self%patch(ixO^S))
        w(ixO^S,phys_ind%mup_) = self%myconfig%mean_mup
     end where
   end if
   cond_tracer_on :if(self%myconfig%tracer_on.and.phys_config%n_tracer>0&
                      .and.self%myconfig%itr<=phys_config%n_tracer)then
     cond_jet_on : if(qt< self%myconfig%time_cla_jet_on)then
      where(.not.self%patch(ixO^S))
       w(ixO^S,phys_ind%tracer(self%myconfig%itr)) = 0.0_dp
      end where
     else cond_jet_on
      if(self%myconfig%tracer_init_density>0.0_dp) then
       where(self%patch(ixO^S))
        w(ixO^S,phys_ind%tracer(self%myconfig%itr)) = self%myconfig%tracer_init_density
       end where
      else
       where(self%patch(ixO^S))
        w(ixO^S,phys_ind%tracer(self%myconfig%itr)) =  w(ixO^S,phys_ind%rho_)
       end where
      end if
     end if cond_jet_on
     itr=itr+1
   end if cond_tracer_on

   cond_dust_on : if(self%myconfig%dust_on)then
       if(self%myconfig%dust_profile/='none') then
         call usr_mat_profile(ixI^L,ixO^L,typeaxial,self%myconfig%dust_profile,&
                              self%myconfig%center,self%myconfig%extend,&
                               x,fprofile)
       else
         fprofile = 1.0_dp
       end if
       call self%mydust%set_patch(ixI^L,ixO^L,self%patch)
       self%mydust%myconfig%velocity=self%myconfig%velocity
       dust_is_frac=.false.
       call self%mydust%set_w(ixI^L,ixO^L,qt,dust_is_frac,self%myconfig%dust_frac,fprofile,x,w)
    end if cond_dust_on
   end if cond_inside_cla_jet



 end subroutine usr_cla_jet_set_w


 !--------------------------------------------------------------------
 subroutine usr_cla_add_source(ixI^L,ixO^L,iw^LIM,x,qdt,qtC,wCT,qt,w,self,&
                               use_tracer,escape_patch,source_filter)
   implicit none
   integer, intent(in)                     :: ixI^L,ixO^L,iw^LIM
   real(kind=dp), intent(in)               :: qdt,qtC,qt
   real(kind=dp), intent(in)               :: x(ixI^S,1:ndim)
   real(kind=dp), intent(in)               :: wCT(ixI^S,1:nw)
   real(kind=dp), intent(inout)            :: w(ixI^S,1:nw)
   logical, intent(in), optional           :: use_tracer
   logical, intent(in),optional            :: escape_patch(ixI^S)
   real(kind=dp), intent(in),optional      :: source_filter(ixI^S)

   class(cla_jet)                          :: self
   ! .. local ..
   real(kind=dp)                           :: source_filter_loc(ixI^S)
   real(kind=dp)                           :: f_profile(ixI^S,1:ndim)
   real(kind=dp)                           :: w_init(ixI^S,1:nw)
   integer                                 :: idir
   !---------------------------------------------------------



end subroutine usr_cla_add_source

 !--------------------------------------------------------------------
 !> Subroutine to process variables in cloud object
  subroutine usr_cla_process_grid(ixI^L,ixO^L,qt,x,w,self)
   implicit none
   integer, intent(in)            :: ixI^L,ixO^L
   real(kind=dp), intent(in)      :: qt
   real(kind=dp)                  :: x(ixI^S,1:ndim)
   real(kind=dp)                  :: w(ixI^S,1:nw)
   class(cla_jet)                 :: self
   ! .. local ..

   !----------------------------------------------------------

   cond_var_0 : if(self%myconfig%variation_on.and.&
                   self%myconfig%variation_position(zjet_)>self%myconfig%z_impos)then
     call self%add_ejecta(ixI^L,ixO^L,qt,.true.,x,w)
   end if cond_var_0


   cond_dust_on : if(self%myconfig%dust_on)then
     call self%mydust%handel_small_val(ixI^L,ixO^L,qt,x,w)
   end if cond_dust_on
 end subroutine usr_cla_process_grid


 !--------------------------------------------------------------------
  subroutine usr_cla_ejecta_set_patch(ixI^L,ixO^L,qt,x,w,self)
     integer, intent(in)            :: ixI^L,ixO^L
     real(kind=dp), intent(in)      :: qt
     real(kind=dp), intent(in)      :: x(ixI^S,1:ndim)
     real(kind=dp), intent(in)      :: w(ixI^S,1:nw)
     class(cla_jet)                 :: self
     ! .. local ..
     integer                        :: iw,n_cells
     real(kind=dp)                  :: z0,ztot,r0
     real(kind=dp), dimension(ixI^S):: jet_radius, jet_h
     !------------------------------------------------------
     if(allocated(self%ejecta_patch))deallocate(self%ejecta_patch)
     allocate(self%ejecta_patch(ixI^S))
      n_cells= self%myconfig%variation_n_cells !1!2
      cond_time_var : if(qt>=self%myconfig%variation_start_time&
                         .and.qt<=self%myconfig%variation_end_time) then
       !------------------------------------------------
       ! * If the actual jet shape is conical :
       !------------------------------------------------
        cond_conical : if(trim(self%myconfig%shape)=='conical')then
          !------------------------------------------------
          ! ** Select the actual axial geometry frame used:
          !------------------------------------------------
          select case(typeaxial)
          !------------------------------------------------
          ! ** Case: spherical frame:
          !------------------------------------------------
          case('spherical')
            self%ejecta_patch(ixO^S) = (dabs(x(ixO^S,rjet_)&
                                          -self%myconfig%variation_position(zjet_))&
                                          <dxlevel(zjet_)*n_cells)
!            angle_theta(ixO^S)          =  self%myconfig%open_angle
          !------------------------------------------------
          ! ** End case: spherical frame
          !------------------------------------------------
          !------------------------------------------------
          ! ** Case: slab or cylindrical frame:
          !------------------------------------------------
          case('slab','cylindrical')
            !------------------------------------------------
            ! *** Select the jet head shape
            !------------------------------------------------
            select case(trim(self%myconfig%head_shape))
            !------------------------------------------------
            ! *** If you want the head to be spherical shaped
            ! in 2D :
            !------------------------------------------------
            case('spherical')
              ztot                = self%myconfig%r_out_init/dtan(self%myconfig%open_angle) !
              z0                  = ztot-self%myconfig%variation_position(zjet_) !ztot = z_impos + z0
            !  !radius relative to (-z0,0) at cylindrical r=r_out_init:
              r0                = dsqrt((self%myconfig%variation_position(zjet_)+z0)**2.0_dp&
              +self%myconfig%r_out_init**2.0_dp)
            !  !radius relative to (-z0,0) for every cylindrical (r,z):
              jet_radius(ixO^S)      = dsqrt((x(ixO^S,zjet_)+z0)**2.0_dp&
                                      +x(ixO^S,rjet_)**2.0_dp)
            !  !for every r and z, select cylindrical r that are between
            !  ! r0 +/- dz(maxlevel,z)
              self%ejecta_patch(ixO^S)  = (jet_radius(ixO^S)-r0)<=dxlevel(zjet_)*n_cells



    !           angle_theta(ixO^S)          =  datan(x(ixO^S,r_)/(z0+x(ixO^S,zjet_)))
            !------------------------------------------------
            ! *** Default for head is cartesian shaped:
            !------------------------------------------------
            case default
              self%ejecta_patch(ixO^S)  = dabs(x(ixO^S,zjet_)                           &
                  -self%myconfig%variation_position(zjet_))&
                  <=dxlevel(zjet_)*n_cells
            end select
            !------------------------------------------------
            ! *** End select the jet head shape
            !------------------------------------------------
          !------------------------------------------------
          ! ** End case: slab or cylindrical frame
          !------------------------------------------------
          end select
          !------------------------------------------------
          ! ** End select the actual axial geometry frame used
          !------------------------------------------------
        !------------------------------------------------
        ! * If the actual jet shape is anything else :
        ! * 'cylindrical', 'cartesian'
        !------------------------------------------------
        else  cond_conical
          !------------------------------------------------
          ! *** Select the jet head shape
          !------------------------------------------------
          select case(trim(self%myconfig%head_shape))
          !------------------------------------------------
          ! *** If you want the head to be spherical shaped
          ! in 2D :
          !------------------------------------------------
          case('spherical')
            !height z relative to (0,0) at any cylindrical r:
            jet_h(ixO^S)      = self%myconfig%variation_position(zjet_)&
            +dsqrt(dabs(self%myconfig%r_out_init**2.0_dp-x(ixO^S,rjet_)**2.0_dp))
            !for every r and z, select z that are below jet_h
            self%ejecta_patch(ixO^S)  = (x(ixO^S,zjet_)-jet_h(ixO^S))&
                                      <=dxlevel(zjet_)*n_cells



          !           angle_theta(ixO^S)          =  datan(x(ixO^S,r_)/(z0+x(ixO^S,zjet_)))
          !------------------------------------------------
          ! *** Default for head is cartesian shaped:
          !------------------------------------------------
          case default
            self%ejecta_patch(ixO^S) = (x(ixO^S,zjet_)                           &
                                           -self%myconfig%variation_position(zjet_))&
                                           <=dxlevel(zjet_)*n_cells
          end select
          !------------------------------------------------
          ! *** End select the jet head shape
          !------------------------------------------------


    !      angle_theta(ixO^S)  = 0.0_dp

        end if cond_conical
        !------------------------------------------------
        ! * End if the actual jet shape is conical
        ! * or if the actual jet shape is anything else
        !------------------------------------------------


        cond_inside : if(any(self%ejecta_patch(ixO^S)))then

          if(.not.allocated(self%patch))then
            call self%set_patch(ixI^L,ixO^L,qt,x)
          end if

          self%ejecta_patch(ixO^S) = self%ejecta_patch(ixO^S).and.self%patch(ixO^S)

        end if  cond_inside
      end if cond_time_var
  end subroutine usr_cla_ejecta_set_patch

 subroutine usr_cla_add_ejecta(ixI^L,ixO^L,qt,is_conserved,x,w,self)
   integer, intent(in)            :: ixI^L,ixO^L
   real(kind=dp), intent(in)      :: qt
   logical, intent(in)            :: is_conserved
   real(kind=dp)                  :: x(ixI^S,1:ndim)
   real(kind=dp)                  :: w(ixI^S,1:nw)
   class(cla_jet)                 :: self
   ! .. local ..
   integer                        :: iw,n_cells
   real(kind=dp)                  :: Vprofile,z0,ztot,amplitude_profile
  ! logical, dimension(ixI^S)      :: patch_ejecta_surface
   real(kind=dp), dimension(ixI^S)       :: angle_theta,jet_radius
   real(kind=dp), dimension(ixI^S,1:ndim)::project_speed
   !------------------------------------------------------
    n_cells= self%myconfig%variation_n_cells
    cond_time_var : if(qt>=self%myconfig%variation_start_time&
                       .and.qt<=self%myconfig%variation_end_time) then

      cond_conical : if(trim(self%myconfig%shape)=='conical')then
        select case(typeaxial)
          case('spherical')
            angle_theta(ixO^S)          =  self%myconfig%open_angle
            project_speed(ixO^S,rjet_)  = dsin(self%myconfig%open_angle)
            project_speed(ixO^S,zjet_)  = dcos(self%myconfig%open_angle)
          case('slab','cylindrical')
            if(ndim>1)then
              ztot                = self%myconfig%r_out_init/dtan(self%myconfig%open_angle)
              z0 = ztot -self%myconfig%z_out_init
              jet_radius(ixO^S) = dsqrt((z0 + x(ixO^S,zjet_))**2.0_dp+x(ixO^S,rjet_)**2.0_dp)
              angle_theta(ixO^S)          =  datan(x(ixO^S,r_)/(z0+x(ixO^S,zjet_)))
              project_speed(ixO^S,rjet_)  = x(ixO^S,r_)/jet_radius(ixO^S)
              project_speed(ixO^S,zjet_)  = (z0 +x(ixO^S,z_))/jet_radius(ixO^S)
            else
              z0                          = 0.0_dp
              angle_theta(ixO^S)          = 0.0_dp
              project_speed(ixO^S,rjet_)  = 0.0_dp
              project_speed(ixO^S,zjet_)  = 1.0_dp
            end if
          end select
      else  cond_conical
        angle_theta(ixO^S)  = 0.0_dp
        project_speed(ixO^S,rjet_)  = 0.0_dp
        project_speed(ixO^S,zjet_)  = 1.0_dp
      end if cond_conical

     if(.not.allocated(self%ejecta_patch))then
        call self%ejecta_set_patch(ixI^L,ixO^L,qt,x,w)
     end if


     cond_inside : if(any(self%ejecta_patch(ixO^S)))then


        call usr_mat_profile_scalar(qt,self%myconfig%variation_time,&
                                    self%myconfig%variation_type,Vprofile)

        if(is_conserved)call phys_to_primitive(ixI^L,ixO^L,w,x)



      !   Loop_var: do iw=1,self%myconfig%variation_phys_nvariable
      !    if(self%myconfig%variation_phys_variable(iw)==phys_ind%mom(zjet_))then
           where(self%ejecta_patch(ixO^S))

            w(ixO^S,phys_ind%mom(zjet_))= (self%myconfig%velocity_poloidal+&
                                      self%myconfig%variation_velocity_poloidal*&
                                      Vprofile)*project_speed(ixO^S,zjet_)
!*(Z0 +x(ixO^S,z_))/jet_radius(ixO^S)
                                      !dcos(angle_theta(ixO^S))
            w(ixO^S,phys_ind%mom(rjet_))= (self%myconfig%velocity_poloidal+&
                                       self%myconfig%variation_velocity_poloidal*&
                                       Vprofile)*project_speed(ixO^S,rjet_)
!*x(ixO^S,r_)/jet_radius(ixO^S)
                                       !dsin(angle_theta(ixO^S))

           end where
           cond_density_variation_on : if(self%myconfig%density_variation_on)then
            amplitude_profile =  self%myconfig%variation_velocity_poloidal/&
                                 self%myconfig%velocity_poloidal
            where(self%ejecta_patch(ixO^S))
             w(ixO^S,phys_ind%rho_) = w(ixO^S,phys_ind%rho_)/&
                                      (1.0_dp+amplitude_profile*Vprofile)
            end where
            !----------------------------------------------------
            ! If switched on, enabling constant temperature in the jet inlet by changing the pressure
            ! profile accordingly
            !----------------------------------------------------
            cond_pressure_variation_on : if(self%myconfig%pressure_variation_on)then
              where(self%ejecta_patch(ixO^S))
               w(ixO^S,phys_ind%pressure_) = w(ixO^S,phys_ind%pressure_)/&
                                        (1.0_dp+amplitude_profile*Vprofile)
              end where
            end if cond_pressure_variation_on
           end if cond_density_variation_on
        if(is_conserved)call phys_to_conserved(ixI^L,ixO^L,w,x)
      end if cond_inside
    end if  cond_time_var


 end subroutine usr_cla_add_ejecta
!--------------------------------------------------------------------
!> Subroutine to clean array memory of associated with cla_jet object
subroutine usr_cla_jet_clean_memory(self)
  class(cla_jet)    :: self
  if(allocated(self%patch))deallocate(self%patch)
  if(allocated(self%ejecta_patch))deallocate(self%ejecta_patch)
  if(self%myconfig%dust_on)call self%mydust%clean_memory()

end subroutine usr_cla_jet_clean_memory
end module mod_obj_cla_jet
