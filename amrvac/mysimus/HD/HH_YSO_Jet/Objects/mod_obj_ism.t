module mod_obj_ism
  use mod_constants
  use mod_global_parameters
  use mod_physics
  use mod_obj_dust, only : dust
  use mod_obj_global_parameters
  use mod_obj_mat
  use mod_obj_usr_unit
  implicit none

    ! ISM features
    type ism_parameters
      character(len=20)    :: unit           !> physical unit at parameter file
      character(len=78)    :: obj_name       !> Obj name that call it
      logical              :: normalize_done !> ism is the normalisation is already done
      real(dp)             :: density        !> ISM density  (g/cm^3)
      real(dp)             :: number_density !> ISM number density (1/cm^3)
      real(dp)             :: temperature    !> ISM temperature  (K)
      real(dp)             :: pressure       !> ISM pressure
      real(dp)             :: Mach_number    !> ISM mach number
      real(dp)             :: Mach_number_tomov_obj    !> ISM mach number associate with moving object in it
      real(dp)             :: velocity_ofmov_obj(3)  !> ISM velocity of the moving object
      logical              :: temperature_inuse !> ISM the temerature to be use and not pressure
      real(dp)             :: extend(2,3)    !> region in space (cm)
      integer              :: myindice       !> ism indices associated with ism in use

      real(dp)             :: profile_kappa  !> ISM index power in pressure
      real(dp)             :: profile_rw  !> ISM reference radius for Lee 2001 profile
      real(dp)             :: profile_zc     !> ISM typical value of height  (cm)
      real(dp)             :: profile_shiftstart !> ISM typical value of height  (cm)
      real(dp)             :: profile_shift_density(3,2) !> ISM profile shift in density
      real(dp)             :: profile_shift_number_density(3,2) !> ISM profile shift in number density
      real(dp)             :: profile_shift_pressure(3,2) !> ISM profile shift in pressure
      real(dp)             :: profile_shift_Temperature(3,2) !> ISM profile shift in Temperature

      real(dp)             :: profile_center(3)    !> ISM  profile center for variation distance
      character(len=30)    :: profile_pressure     !> ism profile pressure
      logical              :: profile_pressure_on  !> ISM pressure profile on
      logical              :: profile_force_on     !> ISM force profile on
      logical              :: profile_force_gradP_on !> ISM profile numerical gradient of p
      logical              :: self_gravity_force_on     !> ISM self-gravity force on
      character(len=30)    :: self_gravity_force_profile !> ISM self-gravity force profile name

      character(len=30)    :: profile_density     !> ism profile density
      logical              :: profile_density_on  !> ISM density profile on
      logical              :: profile_density_keep_pressure !> ISM density profile but keep pressure
      real(dp)             :: theta_floored

      logical              :: ism_display_parameters     !> ISM display of parameters at the end of the set_profile subroutine
      logical              :: profile_on        !> ISM profile set
      real(dp)             :: profile_idir(3)   !> IMS profile direction
      character(len=30)    :: profile_typeaxial !> ISM profile axial type
      real(dp)             :: velocity(3)    !> ISM velocity (cm/s)
      real(dp)             :: magnetic(3)    !> ISM magnetic field (gauss)
      real(dp)             :: c_sound        !> ISM sound speed
      real(dp)             :: xisigma        !> ISM magentisation
      logical              :: tracer_on      !> logical to set tracer
      integer              :: itr            !> ISM tracer indice
      real(dp)             :: tracer_init_density    !> ISM tracer initial density
      real(dp)             :: tracer_small_density   !> ISM tracer small density cut
      logical              :: reset_on       !> ISM reset
      real(dp)             :: reset_coef     !> ISM relaxation coefficient
      real(dp)             :: reset_dTemperature
      real(dp)             :: reset_distance(3)
      real(dp)             :: reset_scale(3)
      logical              :: boundary_on   !> ISM logical to check if it will use implimented boundary
      character(len=30)    :: boundary_cond(3,2)!> ism boundary condition
      logical              :: dust_on        !> logical to set dust
      real(dp)             :: dust_frac      !> dust fraction

      character(len=30)    :: chemical_gas_type
      real(kind=dp)        :: He_abundance
      real(kind=dp)        :: mean_mass
      real(kind=dp)        :: mean_mup
      real(kind=dp)        :: mean_ne_to_nH
      real(kind=dp)        :: mean_nall_to_nH
    end type ism_parameters

    type ISM
      logical, allocatable            :: patch(:^D&)           !> spatial patch
      logical, allocatable            :: escape_patch(:^D&)    !> spatial patch
      character(len=78)               :: subname               !> subroutine name that call it
      type(ism_parameters)            :: myconfig              !> ISM configuation parameters
      type(usrphysical_unit), pointer :: myphysunit            !> ISM physic unity in use
      type (dust)                     :: mydust                !> ISM dust
      type(usrboundary_type)          :: myboundaries          !> ISM boundary condition
     contains
     !PRIVATE
     PROCEDURE, PASS(self) :: set_default          => usr_ism_set_default
     PROCEDURE, PASS(self) :: set_complet          => usr_ism_set_complet
     PROCEDURE, PASS(self) :: normalize            => usr_ism_normalize
     PROCEDURE, PASS(self) :: set_w                => usr_ism_set_w
     PROCEDURE, PASS(self) :: process_grid         => usr_ism_process_grid
     PROCEDURE, PASS(self) :: read_parameters      => usr_ism_read_p
     PROCEDURE, PASS(self) :: write_setting        => usr_ism_write_setting
     PROCEDURE, PASS(self) :: alloc_set_patch      => usr_ism_alloc_set_patch
     PROCEDURE, PASS(self) :: clean_memory         => usr_ism_clean_memory
     PROCEDURE, PASS(self) :: get_patch_escape     => usr_ism_get_patch_escape

     PROCEDURE, PASS(self) :: set_profile          => usr_ism_set_profile
     PROCEDURE, PASS(self) :: get_pforce_profile   => usr_ism_get_pforce_profile
     PROCEDURE, PASS(self) :: get_fgravity_profile  => usr_ism_get_fgravity_profile
     PROCEDURE, PASS(self) :: set_profile_distance => usr_ism_set_profile_distance
     PROCEDURE, PASS(self) :: add_source           => usr_ism_add_source

    end type



contains


  !--------------------------------------------------------------------
   !> Read the ism parameters  from a parfile
    subroutine usr_ism_read_p(self,ism_config,files)
      use mod_obj_mat
      implicit none
      class(ism)                         :: self
      character(len=*),intent(in)        :: files(:)
      type(ism_parameters), intent(out)  :: ism_config
      ! .. local ..
      integer                            :: i_file,i_error_read

      namelist /usr_ism_list/  ism_config
      namelist /usr_ism1_list/ ism_config
      namelist /usr_ism2_list/ ism_config
      namelist /usr_ism3_list/ ism_config

      if(mype==0)write(*,*)'Reading usr_ism_list'
      do i_file = 1, size(files)
         open(unitpar, file=trim(files(i_file)), status="old")
         select case(ism_config%myindice)
         case(1)
           read(unitpar, usr_ism1_list, iostat=i_error_read)
         case(2)
           read(unitpar, usr_ism2_list, iostat=i_error_read)
         case(3)
           read(unitpar, usr_ism3_list, iostat=i_error_read)
         case default
           read(unitpar, usr_ism_list, iostat=i_error_read)
         end select
         call usr_mat_read_error_message(i_error_read,ism_config%myindice,&
                                         self%myconfig%obj_name)
         close(unitpar)
      end do



      if(ism_config%boundary_on)then
        self%myboundaries%myconfig%myindice =1
        call self%myboundaries%read_parameters(self%myboundaries%myconfig,files)
      end if
      if(ism_config%dust_on)then
        self%mydust%myconfig%associated_medium = 'ism'
        call self%mydust%read_parameters(self%mydust%myconfig,files)
      end if


    end subroutine usr_ism_read_p
   !------------------------------------------------------------------------
   subroutine usr_ism_write_setting(self,unit_config)
     implicit none
     class(ism)                          :: self
     integer,intent(in)                  :: unit_config
     ! .. local ..

     !-----------------------------------

     write(unit_config,*)'************************************'
     write(unit_config,*)'************ISM setting ************'
     write(unit_config,*)'************************************'
     write(unit_config,*)'      ****** Code Unit *******      '
     write(unit_config,*) 'Density     = ',  self%myconfig%density, '  code unit'
     write(unit_config,*) 'Pressure    = ',  self%myconfig%pressure, '  code unit'
     write(unit_config,*) 'Temperature = ',  self%myconfig%temperature, '  code unit'
     write(unit_config,*) 'Speed       = ',  self%myconfig%velocity, '  code unit'
     write(unit_config,*) 'Sound speed =' ,  self%myconfig%c_sound, '  code unit'
     write(unit_config,*)'      ****** Physical Unit *******   '
     write(unit_config,*) 'Density     = ',  self%myconfig%density*self%myphysunit%myconfig%density,&
                                             '  ',self%myphysunit%myunit%density
     write(unit_config,*) 'Pressure    = ',  self%myconfig%pressure*self%myphysunit%myconfig%pressure,&
                                             '  ',self%myphysunit%myunit%pressure
     write(unit_config,*) 'Temperature = ',  self%myconfig%temperature*self%myphysunit%myconfig%temperature,&
                                             '  ',self%myphysunit%myunit%temperature
     write(unit_config,*) 'Speed       = ',  self%myconfig%velocity*self%myphysunit%myconfig%velocity,&
                                             '  ',self%myphysunit%myunit%velocity
     write(unit_config,*) 'Sound speed =' ,  self%myconfig%c_sound*self%myphysunit%myconfig%velocity,&
                                             '  ',self%myphysunit%myunit%velocity
     if(self%myconfig%dust_on) then
      call self%mydust%write_setting(unit_config)
     end if
     write(unit_config,*)'************************************'
     write(unit_config,*)'******** END ISM setting **********'
     write(unit_config,*)'************************************'
   end    subroutine usr_ism_write_setting
   !-------------------------------------------------------------------------
   !> subroutine default setting for ISM
   subroutine usr_ism_set_default(self)
    implicit none
    class(ism)            :: self
    !----------------------------------
     self%myconfig%obj_name              = 'ism'
     self%myconfig%unit                  = 'code'
     self%myconfig%myindice              = 0
     self%myconfig%density               = 0.0_dp
     self%myconfig%number_density        = 0.0_dp
     self%myconfig%temperature           = 0.0_dp
     self%myconfig%pressure              = 0.0_dp
     self%myconfig%Mach_number_tomov_obj = 0.0_dp
     self%myconfig%velocity_ofmov_obj    = 0.0_dp
     self%myconfig%Mach_number           = 0.0_dp
     self%myconfig%temperature_inuse     = .false.
     self%myconfig%extend(1:2,1:3)       = 0.0_dp!box_limit(1:2,1:ndim)
     self%myconfig%velocity              = 0.0_dp
     self%myconfig%magnetic              = 0.0_dp
     self%myconfig%xisigma               = 0.0_dp



     self%myconfig%profile_kappa         = 0.0_dp
     self%myconfig%profile_rw            = 2.5d15  !from Lee 2001 profile
     self%myconfig%profile_zc            = 0.0_dp
     self%myconfig%profile_shiftstart    = 0.0_dp

     self%myconfig%profile_shift_number_density = 0.0_dp
     self%myconfig%profile_shift_density        = 0.0_dp
     self%myconfig%profile_shift_pressure       = 0.0_dp
     self%myconfig%profile_shift_Temperature    = 0.0_dp

     self%myconfig%profile_pressure       = 'none'
     self%myconfig%self_gravity_force_profile       = 'none'
     self%myconfig%profile_pressure_on    = .false.
     self%myconfig%profile_force_on       = .false.
     self%myconfig%self_gravity_force_on       = .false.
     self%myconfig%profile_force_gradP_on = .false.

     self%myconfig%ism_display_parameters    = .false.
     self%myconfig%profile_on            = .false.
     self%myconfig%profile_idir          = 0

     self%myconfig%profile_density       = 'none'
     self%myconfig%profile_density_on    = .false.
     self%myconfig%theta_floored         = 0.0_dp
     self%myconfig%profile_typeaxial     = 'slab'
     self%myconfig%profile_center(3)     = 0.0_dp
     self%myconfig%profile_density_keep_pressure =.false.


     self%myconfig%reset_coef            = 0.0_dp
     self%myconfig%reset_dtemperature    = 0.0_dp
     self%myconfig%reset_distance        = 0.0_dp
     self%myconfig%reset_scale           = 0.0_dp
     self%myconfig%reset_on              = .false.
     self%myconfig%boundary_on          = .false.
     self%myconfig%boundary_cond         = 'fix'
     self%myconfig%c_sound               = 0.0_dp


     self%myconfig%tracer_on             = .false.
     self%myconfig%itr                   = 0
     self%myconfig%tracer_init_density   = 0.0_dp
     self%myconfig%tracer_small_density  = 0.0_dp

     self%myconfig%dust_on               = .false.
     self%myconfig%dust_frac             = 0.0_dp

     self%myconfig%normalize_done        = .false.


  self%myconfig%He_abundance   = 0.1_dp
  self%myconfig%chemical_gas_type       = 'fullyionised'
  call unit_chemical_ionisation( self%myconfig%He_abundance, &
                                      self%myconfig%chemical_gas_type,     &
                                      self%myconfig%mean_nall_to_nH, &
                                      self%myconfig%mean_mass,       &
                                      self%myconfig%mean_mup,        &
                                      self%myconfig%mean_ne_to_nH,&
                                      'usr_ism_set_default')

     !if(phys_config%dust_on)then
     call self%mydust%set_default

     call self%myboundaries%set_default

     !end if
   end subroutine usr_ism_set_default
  !--------------------------------------------------------------------
  !> subroutine check the parfile setting for ism
  subroutine usr_ism_set_complet(self)

    implicit none
    class(ism)                                :: self
    ! .. local ..
    logical                  :: dust_is_frac
    real(dp)                 :: mp,kb
    integer                  :: idim,iside
    !-----------------------------------

    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
    else
      mp=mp_cgs
      kB=kB_cgs
    end if
    !----------------------------------------------------
    ! set the parameters :
    !   - mean_nall_to_nH = n_e/n_i
    !   - mean_mass = rho_tot/(nH*mH), where nH=n(H+H^+)+2n(H_2)
    !   - mean_mup,mean_ne_to_nH = n(e^-)/nH
    ! from calling the mod_hd_phys/mod_mhd_phys module
    !----------------------------------------------------

    write(*,*) ' mod_obj_ism.t--> usr_ism_set_complet'
    write(*,*) 'self (before change): chemical_composition | ',&
    'xHe | mean_nall_to_nH |',&
    ' mean_mass | mean_ne_to_nH | mean_mup '
    write(*,*) self%myconfig%chemical_gas_type,&
    self%myconfig%He_abundance,&
    self%myconfig%mean_nall_to_nH,&
    self%myconfig%mean_mass,&
    self%myconfig%mean_ne_to_nH,&
    self%myconfig%mean_mup

    call phys_fill_chemical_ionisation(self%myconfig%He_abundance,  &
       self%myconfig%chemical_gas_type,                             &
       self%myconfig%mean_nall_to_nH,self%myconfig%mean_mass,       &
       self%myconfig%mean_mup,self%myconfig%mean_ne_to_nH)

     write(*,*) ' mod_obj_ism.t--> usr_ism_set_complet'
     write(*,*) 'self (after change): chemical_composition | ',&
     'xHe | mean_nall_to_nH |',&
     ' mean_mass | mean_ne_to_nH | mean_mup '
     write(*,*) self%myconfig%chemical_gas_type,&
     self%myconfig%He_abundance,&
     self%myconfig%mean_nall_to_nH,&
     self%myconfig%mean_mass,&
     self%myconfig%mean_ne_to_nH,&
     self%myconfig%mean_mup

   !----------------------------------------------------
   ! set the density and/or number density parameters accordingly
   ! (provided density takes priority)
   !----------------------------------------------------
    if (dabs( self%myconfig%density)<smalldouble*mp)then
         self%myconfig%density= self%myconfig%number_density*mp*self%myconfig%mean_mass
    elseif(dabs( self%myconfig%number_density)<smalldouble)then
         self%myconfig%number_density= self%myconfig%density/(mp*self%myconfig%mean_mass)
    end if


    !----------------------------------------------------
    ! According to what is input, set the ism
    ! speed of sound c_sound,ism or its Mach Number Ma
    ! (Priority order is :
    !   1) Ma
    !   2) c_sound,ism
    !   3) Mach number relat. to moving objet = ||v_jet||/c_sound,ism)
    !----------------------------------------------------
    cond_Mach_set : if(self%myconfig%mach_number>0.0_dp) then
      self%myconfig%c_sound = dsqrt(sum(self%myconfig%velocity**2.0_dp))/&
        self%myconfig%mach_number
    else if(self%myconfig%c_sound>0.0_dp) then cond_Mach_set
      self%myconfig%mach_number = dsqrt(sum(self%myconfig%velocity**2.0_dp))/self%myconfig%c_sound
    else if(self%myconfig%Mach_number_tomov_obj>0.0_dp.and.any(self%myconfig%velocity_ofmov_obj>0.0_dp)) then cond_Mach_set
      self%myconfig%c_sound = dsqrt(sum(self%myconfig%velocity_ofmov_obj**2.0_dp))/&
        self%myconfig%mach_number_tomov_obj
    end if cond_Mach_set

    !----------------------------------------------------
    ! According to what is input and treated previously,
    ! set the ism pressure from c_sound
    ! 19-03-21 : added the missing setting of temperature accordingly
    ! Priority is for c_sound
    !----------------------------------------------------
    cond_csound_set : if(self%myconfig%c_sound>0.0_dp) then
       self%myconfig%pressure = self%myconfig%c_sound**2.0_dp * self%myconfig%density /&
                                phys_config%gamma
       self%myconfig%temperature = self%myconfig%c_sound**2.0_dp*mp*&
                                self%myconfig%mean_mup/(kB*phys_config%gamma)
    !----------------------------------------------------
    ! *If not yet provided or input,
    ! *set the pressure parameter of the ambient medium
    ! *or the ism temperature.
    ! *At the end, we can
    ! *assign the sound speed accordingly if not yet set
    ! *19-03-21 : made sure that temperature isn't changed if provided
    ! * and takes priority over pressure, which can be changed
    !----------------------------------------------------
    else  cond_csound_set
     if(self%myconfig%temperature>0.0_dp) then
      self%myconfig%pressure = self%myconfig%density*&
                            kB* self%myconfig%temperature/(self%myconfig%mean_mup*mp)
     else if(dabs(self%myconfig%pressure)>0.0_dp) then
      self%myconfig%temperature = self%myconfig%mean_mup*mp*self%myconfig%pressure&
                            /(kB*self%myconfig%density)
     end if
     self%myconfig%c_sound = sqrt(phys_config%gamma*self%myconfig%pressure/self%myconfig%density)
    end if cond_csound_set

    !----------------------------------------------------
    ! End of the ISM physical parameters initialization
    !----------------------------------------------------



    !self%myconfig%profile_idir=min(self%myconfig%profile_idir,ndim)

    self%myconfig%profile_idir=self%myconfig%profile_idir &
                             /dsqrt(sum(self%myconfig%profile_idir**2.0_dp))
    select case(self%myconfig%profile_pressure)
    case('none')
     self%myconfig%profile_pressure_on = .false.
    case default
     if(.not.(dabs(self%myconfig%profile_kappa)>smalldouble.and.self%myconfig%profile_zc>smalldouble)&
        .or.all(self%myconfig%profile_idir<smalldouble))then
       self%myconfig%profile_pressure_on = .false.
       self%myconfig%profile_force_on    = .false.
     end if
    end select

    select case(self%myconfig%profile_density)
    case('none')
     self%myconfig%profile_density_on = .false.
    case default
      if(self%myconfig%profile_zc<smalldouble.or.all(self%myconfig%profile_idir<=smalldouble))then
       self%myconfig%profile_density_on = .false.
       self%myconfig%profile_force_on   = .false.
      end if
    !  if(self%myconfig%profile_force_on)self%myconfig%profile_pressure_on =.true.
    end select

    Loop_idim_shift :  do idim=1,ndim
      Loop_iside_shift :  do iside = 1,2
        if (dabs( self%myconfig%profile_shift_density(idim,iside))<smalldouble*mp.and. &
         dabs( self%myconfig%profile_shift_number_density(idim,iside))<smalldouble)then
          self%myconfig%profile_shift_density(idim,iside) = self%myconfig%density
          self%myconfig%profile_shift_number_density(idim,iside) = self%myconfig%number_density
        else
          if(dabs( self%myconfig%profile_shift_density(idim,iside))<smalldouble*mp)then
             self%myconfig%profile_shift_density(idim,iside) = &
                 self%myconfig%profile_shift_number_density(idim,iside)*mp*self%myconfig%mean_mass
          elseif(dabs( self%myconfig%profile_shift_number_density(idim,iside))<smalldouble)then
             self%myconfig%profile_shift_number_density(idim,iside) = &
                self%myconfig%profile_shift_density(idim,iside)/(mp*self%myconfig%mean_mass)
          end if
        end if
        if(dabs( self%myconfig%profile_shift_pressure(idim,iside))<smalldouble.and. &
          self%myconfig%profile_shift_temperature(idim,iside)<smalldouble) then
            self%myconfig%profile_shift_pressure(idim,iside) = self%myconfig%pressure
            self%myconfig%profile_shift_temperature(idim,iside) = self%myconfig%temperature
        else
          if(dabs( self%myconfig%profile_shift_pressure(idim,iside))<smalldouble)then
           self%myconfig%profile_shift_pressure(idim,iside) = &
              self%myconfig%profile_shift_density(idim,iside) &
              /(self%myconfig%mean_mup*mp)*kB* &
              self%myconfig%profile_shift_temperature(idim,iside)
          else
            self%myconfig%profile_shift_temperature(idim,iside) =  &
              self%myconfig%mean_mup*mp*self%myconfig%profile_shift_pressure(idim,iside)  / &
              (self%myconfig%profile_shift_density(idim,iside)*kB)
          end if
        end if
      end do Loop_iside_shift
    end do Loop_idim_shift
    if( self%myconfig%dust_on)then
      dust_is_frac=.false.
      call  self%mydust%set_complet(dust_is_frac, self%myconfig%dust_frac,&
                              self%myconfig%density, self%myconfig%velocity)
    end if


   cond_traceron : if(self%myconfig%tracer_on)then
     prim_wnames(self%myconfig%itr+(phys_ind%tracer(1)-1)) = 'tracer_ism'
     cons_wnames(self%myconfig%itr+(phys_ind%tracer(1)-1)) = 'tracer_ism'
   else cond_traceron ! OFF NOW
     self%myconfig%itr=0
   end if cond_traceron

   if(self%myconfig%reset_on)then
      Loop_dim : do idim=1,ndim
       if(self%myconfig%reset_scale(idim)<smalldouble)then
         self%myconfig%reset_scale(idim) =box_limit(2,idim)
        end if
      end do Loop_dim
   end if
   if(self%myconfig%profile_density_on.or.self%myconfig%profile_pressure_on)then
     self%myconfig%profile_on   = .true.
   end if
    if(.not.self%myconfig%boundary_on)                                       &
        self%myboundaries%myconfig%boundary_type=self%myconfig%boundary_cond
    call self%myboundaries%set_complet

    if(mype==0)then
      print*,'the ism temperature = ',self%myconfig%temperature
      print*, 'ism density', self%myconfig%density
      print*, 'ism nH', self%myconfig%number_density
    end if
  end subroutine usr_ism_set_complet
  !--------------------------------------------------------------------
  !> subroutine normalize setting for ISM
   subroutine usr_ism_normalize(self,physunit_inuse)
    use mod_obj_usr_unit
    implicit none
    class(ism)                                     :: self
    type(usrphysical_unit), target,intent(in)      :: physunit_inuse
    !----------------------------------
    self%myphysunit =>physunit_inuse

    write(*,*) ' mod_obj_ism.t--> usr_ism_normalize'
    write(*,*) 'self (before normalization) : chemical_composition | ',&
    'xHe | mean_nall_to_nH |',&
    ' mean_mass | mean_ne_to_nH | mean_mup '
    write(*,*) self%myconfig%chemical_gas_type,&
    self%myconfig%He_abundance,&
    self%myconfig%mean_nall_to_nH,&
    self%myconfig%mean_mass,&
    self%myconfig%mean_ne_to_nH,&
    self%myconfig%mean_mup

    if(trim(self%myconfig%unit)=='code'.or.self%myconfig%normalize_done)then
       if(self%myconfig%normalize_done)then
        write(*,*) 'WARNING: Second call for ISM normalisation', &
                     'no new normalisation will be done'
       end if
       return
    end if
     !self%myconfig%mup                =  self%myconfig%mup            /physunit_inuse%myconfig%mup
     self%myconfig%density            =  self%myconfig%density       /physunit_inuse%myconfig%density
     self%myconfig%number_density     =  self%myconfig%number_density/physunit_inuse%myconfig%number_density
     self%myconfig%temperature        =  self%myconfig%temperature   /physunit_inuse%myconfig%temperature
     self%myconfig%pressure           =  self%myconfig%pressure      /physunit_inuse%myconfig%pressure
     self%myconfig%velocity           =  self%myconfig%velocity      /physunit_inuse%myconfig%velocity
     self%myconfig%velocity_ofmov_obj =  self%myconfig%velocity_ofmov_obj  /physunit_inuse%myconfig%velocity
     self%myconfig%extend             =  self%myconfig%extend        /physunit_inuse%myconfig%length
     self%myconfig%c_sound            =  self%myconfig%c_sound       /physunit_inuse%myconfig%velocity
     self%myconfig%mean_mup           =  self%myconfig%mean_mup      / physunit_inuse%myconfig%mean_mup

     self%myconfig%theta_floored = self%myconfig%theta_floored *(dpi/180._dp)
     self%myconfig%profile_rw            = self%myconfig%profile_rw/physunit_inuse%myconfig%length
     self%myconfig%profile_zc              =  self%myconfig%profile_zc           /physunit_inuse%myconfig%length
     self%myconfig%profile_center              =  self%myconfig%profile_center           /physunit_inuse%myconfig%length
     self%myconfig%profile_shiftstart      =  self%myconfig%profile_shiftstart   /physunit_inuse%myconfig%length
     self%myconfig%profile_shift_density          =  self%myconfig%profile_shift_density        /physunit_inuse%myconfig%density
     self%myconfig%profile_shift_number_density   =  self%myconfig%profile_shift_number_density /physunit_inuse%myconfig%number_density
     self%myconfig%profile_shift_temperature      =  self%myconfig%profile_shift_temperature    /physunit_inuse%myconfig%temperature
     self%myconfig%profile_shift_pressure         =  self%myconfig%profile_shift_pressure       /physunit_inuse%myconfig%pressure
    if( self%myconfig%dust_on)then
      call self%mydust%normalize(physunit_inuse)
      call self%mydust%to_phys
    end if
    self%myconfig%normalize_done =.true.

    write(*,*) ' mod_obj_ism.t--> usr_ism_normalize'
    write(*,*) 'self (after normalization): chemical_composition | ',&
    'xHe | mean_nall_to_nH |',&
    ' mean_mass | mean_ne_to_nH | mean_mup '
    write(*,*) self%myconfig%chemical_gas_type,&
    self%myconfig%He_abundance,&
    self%myconfig%mean_nall_to_nH,&
    self%myconfig%mean_mass,&
    self%myconfig%mean_ne_to_nH,&
    self%myconfig%mean_mup

   end subroutine usr_ism_normalize
  !--------------------------------------------------------------------
   !> subroutine setting for ISM
   subroutine usr_ism_set_w(ixI^L,ixO^L,qt,x,w,self,isboundary_iB)
    implicit none
    integer, intent(in)           :: ixI^L,ixO^L
    real(kind=dp), intent(in)     :: qt
    real(kind=dp), intent(in)     :: x(ixI^S,1:ndim)
    real(kind=dp)                 :: w(ixI^S,1:nw)
    integer,             optional :: isboundary_iB(2)
    class(ism)                    :: self
    ! .. local..
    integer                    :: idir,IB,idims,idims_bound,iw
    real(kind=dp)              :: fprofile(ixI^S)
    logical                    :: isboundary
    character(len=30)          :: myboundary_cond
    !----------------------------------

    if(.not.allocated(self%patch)) then
     allocate(self%patch(ixI^S))
     self%patch              = .true.
    end if
    cond_B_present : if(present(isboundary_iB))then
      {^D&
        idims=^D
        if(isboundary_iB(1)==idims.and.isboundary_iB(2)==1) then
         if(all(x(ixO^S,idims)<=xprobmin^D))then
           isboundary=.true.
           IB = 2*(idims-1)+1
           idims_bound = idims
         end if
        end if
        if(isboundary_iB(1)==idims.and.isboundary_iB(2)==2) then
         if(all(x(ixO^S,idims)>=xprobmax^D))then
           isboundary=.true.
           IB = 2*idims
           idims_bound = idims
         end if
        end if
      \}
      myboundary_cond = self%myconfig%boundary_cond(isboundary_iB(1),isboundary_iB(2))
    else cond_B_present
      isboundary=.false.
      myboundary_cond = 'fix'
    end if cond_B_present

    boundary_cond : if(.not.isboundary.or.&
                       trim(myboundary_cond)=='fix') then


      where(self%patch(ixO^S))
        w(ixO^S,phys_ind%rho_)        =  self%myconfig%density
      end where
      if(phys_config%energy)then
        where(self%patch(ixO^S))
          w(ixO^S,phys_ind%pressure_)   =  self%myconfig%pressure
        end where
      end if

      Loop_idir : do idir=1,ndir
       where(  self%patch(ixO^S))
       !----------------------------------------------------
       !line to modify if we want to change the velocity
       !put for each time step at any boundary
       !set to 'fix'
       !----------------------------------------------------
        w(ixO^S,phys_ind%mom(idir)) = self%myconfig%velocity(idir)
       end where
      end do Loop_idir


      ! add profile to ISM density and pressure
      if(self%myconfig%profile_on) then
         call self%set_profile(ixI^L,ixO^L,x,w)
      end if




      if(phys_config%mean_mup_on) then
        where(self%patch(ixO^S))
         w(ixO^S,phys_ind%mup_) = self%myconfig%mean_mup
        end where
      end if

      cond_tracer_on :if(self%myconfig%tracer_on.and.phys_config%n_tracer>0&
                         .and.self%myconfig%itr<=phys_config%n_tracer)then
        if(self%myconfig%tracer_init_density>0.0_dp) then
        where(self%patch(ixO^S))
         w(ixO^S,phys_ind%tracer(self%myconfig%itr)) = self%myconfig%tracer_init_density
        end where
        else
        where(self%patch(ixO^S))
         w(ixO^S,phys_ind%tracer(self%myconfig%itr)) =  w(ixO^S,phys_ind%rho_)
        end where
        end if
        itr=itr+1
      end if cond_tracer_on



      cond_dust_on : if( self%myconfig%dust_on)then
        call self%mydust%set_patch(ixI^L,ixO^L,self%patch)
        self%mydust%myconfig%velocity= self%myconfig%velocity
        fprofile = 1.0_dp
        call   self%mydust%set_w(ixI^L,ixO^L,qt,.false., self%myconfig%dust_frac,fprofile,x,w)
      end if cond_dust_on
    else boundary_cond
     if(any(self%patch(ixO^S))) then
       call self%myboundaries%set_w(ixI^L,ixO^L,iB,isboundary_iB(1),isboundary_iB(2),&
                                self%patch,x,w)
     end if
    end if boundary_cond

   end subroutine usr_ism_set_w

   !--------------------------------------------------------------------
  subroutine usr_ism_set_profile_distance(ixI^L,ixO^L,x,d_profile,self)
    implicit none
    integer, intent(in)             :: ixI^L,ixO^L
    real(kind=dp), intent(in)       :: x(ixI^S,1:ndim)
    class(ism)                      :: self
    ! ..local ..
    integer                         :: idims
    real(kind=dp), dimension(ixI^S) :: d_profile
    !----------------------------------------------------
    select case(self%myconfig%profile_typeaxial)
    case('projection_direction')
      Loop_idims0: do idims = 1,ndim
       d_profile(ixO^S) = self%myconfig%profile_idir(idims) * dabs(x(ixO^S,idims))
      end do Loop_idims0
    case default
      call usr_distance(ixI^L,ixO^L,typeaxial,&
                        self%myconfig%profile_center,x,d_profile)
    end select
  end subroutine usr_ism_set_profile_distance

   subroutine usr_ism_set_profile(ixI^L, ixO^L,x,w,self)
    implicit none
    integer, intent(in)             :: ixI^L,ixO^L
    real(kind=dp), intent(in)       :: x(ixI^S,1:ndim)
    real(kind=dp), intent(inout)    :: w(ixI^S,1:nw)
    class(ism)                      :: self
    ! ..local ..
    real(kind=dp), dimension(ixI^S) :: p_profile,d_profile,theta_profile
    !----------------------------------------------------

    call self%set_profile_distance(ixI^L,ixO^L,x,d_profile)
    cond_pressure_profile : if(self%myconfig%profile_pressure_on) then
     select case(trim(self%myconfig%profile_pressure))
      case('king')
        if(dabs(self%myconfig%profile_kappa)>smalldouble.and.self%myconfig%profile_zc>smalldouble)then
          p_profile(ixO^S) = 1.0_dp/(1.0_dp+(d_profile(ixO^S)&
                            /self%myconfig%profile_zc)**(-self%myconfig%profile_kappa) ) !/
        else
          p_profile(ixO^S) = 1.0_dp
        end if
      case('komissarov')
        if(dabs(self%myconfig%profile_kappa)>smalldouble.and.self%myconfig%profile_zc>smalldouble)then
          p_profile(ixO^S) = (d_profile(ixO^S)&
                              /self%myconfig%profile_zc)**(-self%myconfig%profile_kappa) !/
        else
          p_profile(ixO^S) = 1.0_dp
        end if
      case default
        p_profile(ixO^S) = 1.0_dp
      end select
    end if   cond_pressure_profile

    cond_density_profile : if(self%myconfig%profile_density_on) then
     select case(trim(self%myconfig%profile_density))
        case('cabrit1997')
        if(self%myconfig%profile_zc>smalldouble)then
          where(d_profile(ixO^S)>self%myconfig%profile_shiftstart)
           p_profile(ixO^S) = (1.0_dp+&
              (d_profile(ixO^S)-self%myconfig%profile_shiftstart)&
              /self%myconfig%profile_zc)**(-self%myconfig%profile_kappa) !/
          else where
            p_profile(ixO^S) = 1.0_dp
          end where
        else
          p_profile(ixO^S) = 1.0_dp
        end if
        case('HSAM')
        if(self%myconfig%profile_zc>smalldouble)then
          where(d_profile(ixO^S)>self%myconfig%profile_shiftstart)
           p_profile(ixO^S) = 1.0_dp/(1.0_dp+&
              ((d_profile(ixO^S)-self%myconfig%profile_shiftstart)&
              /self%myconfig%profile_zc)**(self%myconfig%profile_kappa)) !/
          else where
            p_profile(ixO^S) = 1.0_dp
          end where
        else
          p_profile(ixO^S) = 1.0_dp
        end if
        case('Lee2001','Lee2001_floored')
          select case(typeaxial)
            case('spherical')
              theta_profile(ixO^S)          = x(ixO^S,theta_)
            case('slab','cylindrical')
              if(ndim>1)then
                theta_profile(ixO^S)        = datan(x(ixO^S,r_)/x(ixO^S,z_))
              else
                theta_profile(ixO^S)        = 0.0_dp
              end if
           end select
              !p_profile(ixO^S) = d_profile(ixO^S)*self%myphysunit%myconfig%length!1.0_dp/(d_profile(ixO^S)**(2.0_dp)) !((dsin(theta_profile(ixO^S)))**2.0_dp)
              ! 18-10-20 : les tests montrent qu'ici d_profile(ixO^S) est normalisé alors que
              ! il ne devrait pas l'être, la normalisation des grandeurs intervenant bien plus en avant
              p_profile(ixO^S) = ((dsin(theta_profile(ixO^S)))**2.0_dp)/&
                                 ((d_profile(ixO^S)/self%myconfig%profile_rw)**2.0_dp)
              !print*,p_profile(ixO^S)
          if(trim(self%myconfig%profile_density)=='Lee2001_floored')then
            where(theta_profile(ixO^S)<=self%myconfig%theta_floored)
            p_profile(ixO^S) = ((dsin(self%myconfig%theta_floored))**2.0_dp)/&
                               ((d_profile(ixO^S)/self%myconfig%profile_rw)**2.0_dp)
            end where
          end if
      case default
        p_profile(ixO^S) = 1.0_dp
      end select
    end if   cond_density_profile

    if(self%myconfig%profile_pressure_on) then
      if(phys_config%energy)then
          where(self%patch(ixO^S))
            where(d_profile(ixO^S)>self%myconfig%profile_shiftstart)
              w(ixO^S,phys_ind%pressure_) = w(ixO^S,phys_ind%pressure_) * p_profile(ixO^S)
            elsewhere
                w(ixO^S,phys_ind%pressure_) = &
                   self%myconfig%profile_shift_pressure(1,1)
            end where
          endwhere
      end if
      where(self%patch(ixO^S))
        where(d_profile(ixO^S)>self%myconfig%profile_shiftstart)
          w(ixO^S,phys_ind%rho_)      = w(ixO^S,phys_ind%rho_) * p_profile(ixO^S)
        elsewhere
            w(ixO^S,phys_ind%rho_)      = &
               self%myconfig%profile_shift_density(1,1)
        end where
      endwhere
    elseif(self%myconfig%profile_density_on)then
      where(self%patch(ixO^S))
        where(d_profile(ixO^S)>self%myconfig%profile_shiftstart)
          w(ixO^S,phys_ind%rho_)      = w(ixO^S,phys_ind%rho_)* p_profile(ixO^S)
        elsewhere
          w(ixO^S,phys_ind%rho_)      = &
             self%myconfig%profile_shift_density(1,1)

        end where
      end where
      if(phys_config%energy)then
        where(self%patch(ixO^S))
          where(d_profile(ixO^S)>self%myconfig%profile_shiftstart)
            w(ixO^S,phys_ind%pressure_) =w(ixO^S,phys_ind%pressure_)*p_profile(ixO^S)
          elsewhere
            w(ixO^S,phys_ind%pressure_) = &
               self%myconfig%profile_shift_pressure(1,1)

          end where
        end where
      end if
    end if

if(self%myconfig%ism_display_parameters) then
   print*, 'in mod_obj_ism.t : '
   print*, 'ism density', self%myconfig%density*self%myphysunit%myconfig%density
   print*, 'ism ndensity', self%myconfig%number_density*self%myphysunit%myconfig%number_density
   print*, 'ism pressure', self%myconfig%pressure*self%myphysunit%myconfig%pressure
   print*, 'ism temperature', self%myconfig%temperature*self%myphysunit%myconfig%temperature
   print*, 'ism velocity', self%myconfig%velocity*self%myphysunit%myconfig%velocity
end if
   end subroutine usr_ism_set_profile


   !--------------------------------------------------------------------

   subroutine usr_ism_get_pforce_profile(ixI^L,ixO^L,qt,qdt,x,w,f_profile,divpv,self)
    use mod_radiative_cooling
    implicit none
    integer, intent(in)                  :: ixI^L,ixO^L
    real(kind=dp), intent(in)            :: qt,qdt
    real(kind=dp), intent(in)            :: x(ixI^S,1:ndim)
    real(kind=dp), intent(in)            :: w(ixI^S,1:nw)
    class(ism)                           :: self
    real(kind=dp), intent(inout)         :: f_profile(ixI^S,1:ndim)
    real(kind=dp), intent(inout)         :: divpv(ixI^S)
    ! ..local ..
    integer                              :: idims
    real(kind=dp), dimension(ixI^S,1:ndir) :: pv
    real(kind=dp), dimension(ixI^S,1:nw) :: w_init, w_tmp
    real(kind=dp), dimension(ixI^S)      :: ptherm,gradp,f_projection,d_profile
    logical                              :: tmp_active,src_active
    logical, dimension(ixI^S)            :: patch_tosave
    !----------------------------------------------------
    f_profile(ixO^S,:) = 0.0
    call self%set_profile_distance(ixI^L,ixO^L,x,d_profile)
    cond_pressure_fprofile : if(self%myconfig%profile_pressure_on) then
      select case(trim(self%myconfig%profile_pressure))
      case('komissarov')
       cond_force1: if(dabs(self%myconfig%profile_kappa)&
            >smalldouble.and.self%myconfig%profile_zc>smalldouble)then
        where(self%patch(ixO^S))
         f_projection(ixO^S) =&
            - self%myconfig%profile_kappa/self%myconfig%profile_zc*&
           (d_profile(ixO^S)/self%myconfig%profile_zc)**(-(self%myconfig%profile_kappa+1))
        end  where
       else cond_force1
        where(self%patch(ixO^S))
         f_projection(ixO^S) = 0.0_dp
        end where
       end if cond_force1
      case default
       where(self%patch(ixO^S))
        f_projection(ixO^S) = 0.0_dp
       end where
      end select
    end if cond_pressure_fprofile

    cond_density_fprofile : if(self%myconfig%profile_density_on) then


     cond_fprofile_numerical : if(self%myconfig%profile_force_gradP_on)then
      patch_tosave(ixI^S)=self%patch(ixI^S)
      self%patch(ixI^S)  =.true.
      call self%set_w(ixI^L, ixI^L,qt,x,w_init)
      self%patch(ixI^S)=patch_tosave(ixI^S)

      call phys_to_conserved(ixI^L,ixI^L,w_init,x)

      w_tmp(ixI^S,1:nw) = w_init(ixI^S,1:nw)
!

      if(phys_config%radiative_cooling)then
        tmp_active=.false.
        src_active=.false.
        call radiative_cooling_add_source(qdt,ixI^L,ixI^L,w_init,w_tmp,x,&
                                           src_active,tmp_active)
      end if
      call phys_get_pthermal(w_tmp,x,ixI^L,ixI^L,ptherm)
      Loop_idir_force : do idims=1,ndim
        call gradient(ptherm,ixI^L,ixO^L,idims,gradp)

        where(self%patch(ixO^S))
          f_profile(ixO^S,idims) = gradp(ixO^S)
        end where
        pv(ixI^S,idims) =  ptherm(ixI^S)*w(ixI^S,phys_ind%mom(idims))/w(ixI^S,phys_ind%rho_)
      end do  Loop_idir_force
      if(ndir>ndim)pv(ixI^S,ndim+1:ndir)=0.0_dp
      call  divvector(pv,ixI^L,ixO^L,divpv)


     else cond_fprofile_numerical


      select case(trim(self%myconfig%profile_density))
      case('cabrit1997')
       cond_force_rho1: if(self%myconfig%profile_zc>smalldouble)then

        where(self%patch(ixO^S))
          where(d_profile(ixO^S)>self%myconfig%profile_shiftstart)

            f_projection(ixO^S) = -phys_config%gamma/&
                   self%myconfig%profile_zc*self%myconfig%profile_kappa*&
              (1.0_dp+(d_profile(ixO^S)&
              -self%myconfig%profile_shiftstart)&
              /self%myconfig%profile_zc)**(-phys_config%gamma*&
                                              self%myconfig%profile_kappa-1)

          else where
            f_projection(ixO^S) = 0.0_dp
          end where
        end  where

        divpv(ixO^S) = 0.0_dp ! A REMPLIR !
       else cond_force_rho1
        where(self%patch(ixO^S))
         f_projection(ixO^S) = 0.0_dp
        end where
      end if cond_force_rho1
      case default
       where(self%patch(ixO^S))
        f_projection(ixO^S) = 0.0_dp
       end where
      end select
      select case(self%myconfig%profile_typeaxial)
      case('projection_direction')
        Loop_idim0 : do idims=1,ndim
          where(self%patch(ixO^S))
            f_profile(ixO^S,idims) = f_projection(ixO^S)*self%myconfig%profile_idir(idims)
          end where
        end do Loop_idim0
      case default
        Loop_idim : do idims=1,ndim
          where(self%patch(ixO^S))
            f_profile(ixO^S,idims) = f_projection(ixO^S) &
                          *(x(ixO^S,idims)-self%myconfig%profile_center(idims))/d_profile(ixO^S)
          end where
        end do Loop_idim

      end select
     end if cond_fprofile_numerical
    end if cond_density_fprofile



   end subroutine usr_ism_get_pforce_profile

   !--------------------------------------------------------------------
   ! This subroutine adds gravity potential gradient
   subroutine usr_ism_get_fgravity_profile(ixI^L,ixO^L,qt,qdt,x,w,f_profile,divpv,self)
    use mod_radiative_cooling
    implicit none
    integer, intent(in)                  :: ixI^L,ixO^L
    real(kind=dp), intent(in)            :: qt,qdt
    real(kind=dp), intent(in)            :: x(ixI^S,1:ndim)
    real(kind=dp), intent(in)            :: w(ixI^S,1:nw)
    class(ism)                           :: self
    real(kind=dp), intent(inout)         :: f_profile(ixI^S,1:ndim)
    real(kind=dp), intent(inout)         :: divpv(ixI^S)
    ! ..local ..
    integer                              :: idims
    real(kind=dp), dimension(ixI^S,1:ndir) :: pv
    real(kind=dp), dimension(ixI^S,1:nw) :: w_init, w_tmp
    real(kind=dp), dimension(ixI^S)      :: ptherm,gradp,f_projection,d_profile
    logical                              :: tmp_active,src_active
    logical, dimension(ixI^S)            :: patch_tosave
    !----------------------------------------------------
    !TODO

   end subroutine usr_ism_get_fgravity_profile

   !--------------------------------------------------------------------
   subroutine usr_ism_add_source(ixI^L,ixO^L,iw^LIM,x,qdt,qtC,wCT,qt,w,self,&
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

     class(ism)                              :: self
     ! .. local ..
     real(kind=dp), dimension(ixI^S)         :: source_filter_loc,divpv
     real(kind=dp)                           :: f_profile(ixI^S,1:ndim)
     real(kind=dp), dimension(ixI^S)         :: temperature,temperature_init
     real(kind=dp)                           :: w_init(ixI^S,1:nw)
     real(kind=dp)                           :: coef_eos
     integer                                 :: idir,idims
     logical, dimension(ixI^S)               :: patch_temperature
     !---------------------------------------------------------
     ! This subroutine implicitly uses conservative variables


      call self%alloc_set_patch(ixI^L,ixI^L,qt,x,&
                   use_tracer=use_tracer,w=w,escape_patch=escape_patch)

      cond_add_force : if(self%myconfig%profile_force_on) then

        cond_inside_prof: if(any(self%patch(ixO^S)))then
          call self%get_pforce_profile(ixI^L,ixO^L,qtC,qdt,x,wCT,f_profile,divpv)
          Loop_idim_force_m : do idims = 1,ndim
          where(self%patch(ixO^S))
           ! This subroutine implicitly uses conservative variables
            w(ixO^S,phys_ind%mom(idims)) = w(ixO^S,phys_ind%mom(idims))&
                +qdt*f_profile(ixO^S,idims)

          end where
        end do Loop_idim_force_m
          if(phys_config%energy) then
            coef_eos = phys_config%gamma/(phys_config%gamma-1.0_dp)

            where(self%patch(ixO^S))
            ! This subroutine implicitly uses conservative variables
             w(ixO^S,phys_ind%e_)=w(ixO^S,phys_ind%e_)  &
                + qdt*coef_eos*divpv(ixO^S)
            end where

          end if
        end if  cond_inside_prof
      end if cond_add_force


     cond_reset : if(self%myconfig%reset_coef>0.0_dp.and.phys_config%energy)then
      cond_inside : if(any(self%patch(ixO^S)))then

       cond_filter : if(present(source_filter))then
         source_filter_loc(ixO^S) = source_filter(ixO^S)
       else cond_filter
         source_filter_loc(ixO^S) = 1.0_dp
       end if cond_filter


       call self%set_w(ixI^L,ixO^L,qt,x,w_init)
       !Here the codes needs to pass from conservative to primitive variables w
       call phys_to_primitive(ixI^L,ixO^L,w,x)
       where(self%patch(ixO^S))
         temperature(ixO^S)      = w(ixO^S,phys_ind%pressure_)/w(ixO^S,phys_ind%rho_)
         temperature_init(ixO^S) =  w_init(ixO^S,phys_ind%pressure_)&
                                   / w_init(ixO^S,phys_ind%rho_)
        patch_temperature(ixO^S) = dabs(temperature(ixO^S)-&
                                    temperature_init(ixO^S))&
                                    >=self%myconfig%reset_dtemperature
       elsewhere
          patch_temperature(ixO^S) = .false.
       end where
        cond_rest_pressure: if(any(patch_temperature(ixO^S)))then
          where(self%patch(ixO^S))
            self%patch(ixO^S) =.false.
          end where
        else cond_rest_pressure
          self%patch(ixO^S) =.false.
        end if  cond_rest_pressure
       where(self%patch(ixO^S))
         source_filter_loc(ixO^S) = max(dabs(source_filter_loc(ixO^S)*self%myconfig%reset_coef),1.0_dp)

         w(ixO^S,phys_ind%rho_) =   w(ixO^S,phys_ind%rho_)*(1.0_dp-source_filter_loc(ixO^S)) +&
                          w_init(ixO^S,phys_ind%rho_)*source_filter_loc(ixO^S)
         w(ixO^S,phys_ind%pressure_) =  w(ixO^S,phys_ind%pressure_)*(1.0_dp-source_filter_loc(ixO^S)) +&
                          w_init(ixO^S,phys_ind%pressure_)*source_filter_loc(ixO^S)
       end where

       loop_idir :  do idir = 1,ndir
        where(patch_temperature(ixO^S))
          w(ixO^S,phys_ind%mom(idir)) = w(ixO^S,phys_ind%mom(idir))*(1.0_dp-source_filter_loc(ixO^S)) +&
                          w_init(ixO^S,phys_ind%mom(idir))*source_filter_loc(ixO^S)
        end where
       end do loop_idir



       cond_dust_on : if( self%myconfig%dust_on)then
          call self%mydust%set_patch(ixI^L,ixO^L,self%patch)
          self%mydust%myconfig%velocity= self%myconfig%velocity
          f_profile = 1.0_dp
          call   self%mydust%set_w(ixI^L,ixO^L,qt,.false., &
                                   self%myconfig%dust_frac,f_profile,x,w)
       end if cond_dust_on

       !Don't forget to switch back to conservative variables w
       call phys_to_conserved(ixI^L,ixO^L,w,x)

      end if cond_inside
    end if cond_reset
    !if(any(self%patch(ixO^S)))print*,' test ism force',maxval(dabs(w(ixO^S,phys_ind%mom(z_))),mask=self%patch(ixO^S))
   end subroutine usr_ism_add_source

   !--------------------------------------------------------------------
   !> Subroutine to process variables in cloud object
    subroutine usr_ism_process_grid(ixI^L,ixO^L,qt,x,w,self)
     implicit none
     integer, intent(in)        :: ixI^L,ixO^L
     real(kind=dp), intent(in)  :: qt
     real(kind=dp)              :: x(ixI^S,1:ndim)
     real(kind=dp)              :: w(ixI^S,1:nw)
     class(ism)                 :: self
     ! .. local..
     !----------------------------------------------------------
     cond_dust_on : if(self%myconfig%dust_on)then
       call self%mydust%handel_small_val(ixI^L,ixO^L,qt,x,w)
     end if cond_dust_on
   end subroutine usr_ism_process_grid
   !--------------------------------------------------------------------
   !> Subroutine to clean array memory of associated with cloud object
   subroutine usr_ism_alloc_set_patch(ixI^L,ixO^L,qt,x,self,use_tracer,w,escape_patch)
     implicit none
     integer, intent(in)                     :: ixI^L,ixO^L
     real(kind=dp), intent(in)               :: qt
     real(kind=dp), intent(in)               :: x(ixI^S,1:ndim)
     real(kind=dp), intent(in), optional     :: w(ixI^S,1:nw)
     logical, intent(in), optional           :: use_tracer
     logical, intent(in),optional            :: escape_patch(ixI^S)
     class(ism)                              :: self
     !---------------------------------------------------------
     cond_tracer : if(.not.present(use_tracer).and. .not.present(w)) then
       if(allocated(self%patch))deallocate(self%patch)
       allocate(self%patch(ixI^S))
       self%patch(ixO^S) =.true.

     else cond_tracer

       cond_ismtracer_on : if(self%myconfig%tracer_on)then
        if(allocated(self%patch))deallocate(self%patch)
        allocate(self%patch(ixI^S))
        where(w(ixO^S,phys_ind%tracer(self%myconfig%itr))>small_density)
          self%patch(ixO^S)=.true.
        else where
          self%patch(ixO^S)=.false.
        end where
       end if cond_ismtracer_on

     end if cond_tracer

     if(allocated(self%escape_patch))deallocate(self%escape_patch)
     allocate(self%escape_patch(ixI^S))
     if(present(escape_patch))then
      self%escape_patch(ixO^S) = escape_patch(ixO^S)
     else
      self%escape_patch(ixO^S) =.false.
     end if

     where(self%patch(ixO^S))self%patch(ixO^S)         = .not.self%escape_patch(ixO^S)
   end subroutine usr_ism_alloc_set_patch
  !---------------------------------------------------------------------
   subroutine  usr_ism_get_patch_escape(ixI^L,ixO^L,need_dealloc,escape_patch,self)
     implicit none
     integer, intent(in)           :: ixI^L,ixO^L
     logical, intent(in)           :: need_dealloc
     logical, intent(in)           :: escape_patch(ixI^S)
     class(ism)                    :: self
     !----------------------------------------------------------
     if(allocated(self%escape_patch))deallocate(self%escape_patch)
      allocate(self%escape_patch(ixI^S))
      self%escape_patch(ixO^S) = .false.

     self%escape_patch(ixO^S)=self%escape_patch(ixO^S).or.escape_patch(ixO^S)
   end subroutine  usr_ism_get_patch_escape
   !--------------------------------------------------------------------
   !> Subroutine to clean array memory of associated with cloud object
   subroutine usr_ism_clean_memory(self)
     class(ism)    :: self



     if(allocated(self%patch))deallocate(self%patch)

     if(allocated(self%escape_patch))deallocate(self%escape_patch)
     if(self%myconfig%dust_on)call self%mydust%clean_memory
   end subroutine usr_ism_clean_memory
  !---------------------------------------------------------------------
end module mod_obj_ism
