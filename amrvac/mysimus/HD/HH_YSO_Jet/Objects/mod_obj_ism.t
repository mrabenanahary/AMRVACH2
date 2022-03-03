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
      real(dp)             :: pressure_Ulrich !> ISM pressure in case of Ulrich model
      real(dp)             :: Mach_number    !> ISM mach number
      real(dp)             :: Mach_number_tomov_obj    !> ISM mach number associate with moving object in it
      real(dp)             :: velocity_ofmov_obj(3)  !> ISM velocity of the moving object
      logical              :: temperature_inuse !> ISM the temerature to be use and not pressure
      real(dp)             :: extend(2,3)    !> region in space (cm)
      integer              :: myindice       !> ism indices associated with ism in use

      real(dp)             :: profile_kappa  !> ISM index power in pressure
      real(dp)             :: profile_rw  !> ISM reference radius for Lee 2001 profile
      real(dp)             :: profile_rw_force  !> ISM reference radius for balancing force amplification
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
      !> Ulrich model ambient medium parameters
      real(dp)             :: profile_jspec !> ISM envelope specific angular momentum
      real(dp)             :: profile_Minf  !> ISM envelope infalling mass-loss rate on the protostar
      real(dp)             :: profile_Mstar !> central objet or protostar mass
      real(dp)             :: profile_rd    !> ISM envelope lower ﬁducial radius
      real(dp)             :: profile_vd    !> ISM envelope kepler velocity at r=profile_rd
      real(dp)             :: profile_rho_0 !> ISM envelope fiducial density


      logical              :: ism_display_parameters     !> ISM display of parameters at the end of the set_profile subroutine
      logical              :: profile_on        !> ISM profile set
      logical              :: bc_profile_on        !> ISM profile set
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
      character(len=30)    :: weight_mean !> ISM weighting mean type
      character(len=30)    :: weight_mean_variable !> ISM weighting mean type
      real(dp)             :: weight_analy     !> ISM analytical force weight
      real(dp)             :: weight_num     !> ISM numerical force weight
      real(dp)             :: weight_analy_index !> ISM power law index for amplification
      real(dp)             :: weight_num_index !> ISM power law index for amplification
      real(dp)             :: reset_dTemperature
      real(dp)             :: reset_distance(3)
      real(dp)             :: reset_scale(3)
      logical              :: boundary_on   !> ISM logical to check if it will use implimented boundary
      logical              :: debug
      character(len=30)    :: boundary_cond(3,2)!> ism boundary condition
      real(dp)             :: flux_frac
      logical              :: mixed_fixed_bound(3,2,7) !> ism 'fix' flux variables
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
     PROCEDURE, PASS(self) :: set_ulrich_profile   => usr_set_ulrich_profile
     PROCEDURE, PASS(self) :: set_profile          => usr_ism_set_profile
     PROCEDURE, PASS(self) :: get_pforce_profile   => usr_ism_get_pforce_profile
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
     integer                             :: idims2,iside2
     real(kind=dp)                       :: rto_print
     character(len=64)                   :: sto_print
     ! .. local ..

     !-----------------------------------

     write(unit_config,*)'************************************'
     write(unit_config,*)'************ISM setting ************'
     write(unit_config,*)'************************************'
     write(unit_config,*)'      ****** Code Unit *******      '
     write(unit_config,*) 'Density  at (R,z)=(R_j,0) = ',  self%myconfig%density, '  code unit'
     write(unit_config,*) 'Pressure    = ',  self%myconfig%pressure, '  code unit'
     write(unit_config,*) 'Temperature = ',  self%myconfig%temperature, '  code unit'
     if(trim(self%myconfig%profile_density)=='Ulrich1976')then
      write(unit_config,*) 'Ulrich1976 model used !'
      write(unit_config,*) 'jspec = ',  self%myconfig%profile_jspec, '  code unit'
      write(unit_config,*) 'Mstar = ',  self%myconfig%profile_Mstar, '  code unit'
      write(unit_config,*) 'Mdot_inf = ',  self%myconfig%profile_Minf, '  code unit'
      write(unit_config,*) 'rD = jspec^2/(G*Mstar) = ',  self%myconfig%profile_rd, '  code unit'
      write(unit_config,*) 'vK = sqrt((G*Mstar)/rD) = ',  self%myconfig%profile_vd, '  code unit'
      write(unit_config,*) 'Density rho_0 ',  self%myconfig%profile_rho_0, '  code unit'
     end if
     write(unit_config,*) 'Speed       = ',  self%myconfig%velocity, '  code unit'
     write(unit_config,*) 'Sound speed =' ,  self%myconfig%c_sound, '  code unit'
     write(unit_config,*)'      ****** Physical Unit *******   '
     write(unit_config,*) 'Density     = ',  self%myconfig%density*self%myphysunit%myconfig%density,&
                                             '  ',self%myphysunit%myunit%density
     write(unit_config,*) 'Pressure    = ',  self%myconfig%pressure*self%myphysunit%myconfig%pressure,&
                                             '  ',self%myphysunit%myunit%pressure
     write(unit_config,*) 'Temperature = ',  self%myconfig%temperature*self%myphysunit%myconfig%temperature,&
                                             '  ',self%myphysunit%myunit%temperature
     if(trim(self%myconfig%profile_density)=='Ulrich1976')then

      write(unit_config,*) 'Ulrich1976 model used !'

      rto_print = self%myconfig%profile_jspec*self%myphysunit%myconfig%velocity
      rto_print = rto_print*self%myphysunit%myconfig%length
      sto_print = (self%myphysunit%myunit%velocity // ' ') // self%myphysunit%myunit%length
      write(unit_config,*) 'jspec = ', rto_print, ' ', sto_print

      rto_print = self%myconfig%profile_Mstar*self%myphysunit%myconfig%mass
      sto_print = self%myphysunit%myunit%mass
      write(unit_config,*) 'Mstar = ', rto_print, ' ', sto_print

      rto_print = self%myconfig%profile_Minf*self%myphysunit%myconfig%mass
      rto_print = rto_print/self%myphysunit%myconfig%time
      sto_print = (self%myphysunit%myunit%mass // '/') // self%myphysunit%myunit%time
      write(unit_config,*) 'Mdot_inf = ', rto_print, ' ', sto_print

      rto_print = self%myconfig%profile_rd*self%myphysunit%myconfig%length
      sto_print = self%myphysunit%myunit%length
      write(unit_config,*) 'rD = jspec^2/(G*Mstar) = ', rto_print, ' ', sto_print

      rto_print = self%myconfig%profile_vd*self%myphysunit%myconfig%velocity
      sto_print = self%myphysunit%myunit%velocity
      write(unit_config,*) 'vK = sqrt((G*Mstar)/rD) = ', rto_print, ' ', sto_print

      rto_print = self%myconfig%profile_rho_0*self%myphysunit%myconfig%density
      sto_print = self%myphysunit%myunit%density
      write(unit_config,*) 'Density rho_0 ', rto_print, ' ', sto_print
     end if
     write(unit_config,*) 'Speed       = ',  self%myconfig%velocity*self%myphysunit%myconfig%velocity,&
                                             '  ',self%myphysunit%myunit%velocity
     write(unit_config,*) 'Sound speed =' ,  self%myconfig%c_sound*self%myphysunit%myconfig%velocity,&
                                             '  ',self%myphysunit%myunit%velocity
     if(self%myconfig%dust_on) then
      call self%mydust%write_setting(unit_config)
     end if

     write(unit_config,*)'** Boundary conditions'

     do idims2=1,ndim
      do iside2=1,2
        write(unit_config,*)' *** (idims,iside)=', idims2,iside2
        write(unit_config,*)' **** BC = ', self%myconfig%boundary_cond(idims2,iside2)
        write(unit_config,*) ' **** ... is this BC fixed ? : '
        write(unit_config,*) '--------------------------------'
        if(phys_config%energy)then
          write(unit_config,*) ' rho v1 v2 v3 p(or e) tracer1 tracer2'
          write(unit_config,*) self%myconfig%mixed_fixed_bound(idims2,iside2,1:7)
        else
          write(unit_config,*) ' rho v1 v2 v3 tracer1 tracer2'
          write(unit_config,*) self%myconfig%mixed_fixed_bound(idims2,iside2,1:6)
        end if
        write(unit_config,*) '--------------------------------'

      end do
     end do
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
     self%myconfig%pressure_Ulrich       = 0.0_dp
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
     self%myconfig%profile_rw_force      = 2.5d16  !from Lee 2001 profile
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
     self%myconfig%bc_profile_on            = .true. !< when profile_on = .true.,
     !bc_profile_on should be true by default by convenience
     self%myconfig%profile_idir          = 0

     self%myconfig%profile_density       = 'none'
     self%myconfig%profile_density_on    = .false.
     self%myconfig%theta_floored         = 0.0_dp
     self%myconfig%profile_typeaxial     = 'slab'
     self%myconfig%profile_center(3)     = 0.0_dp
     self%myconfig%profile_density_keep_pressure =.false.
     !Ulrich model ambient medium parameters
     self%myconfig%profile_jspec = 0.0_dp       !> ISM envelope specific angular momentum
     self%myconfig%profile_Minf  = 0.0_dp !> ISM envelope infalling mass-loss rate on the protostar
     self%myconfig%profile_Mstar = 0.0_dp   !> central objet or protostar mass
     self%myconfig%profile_rd    = 0.0_dp !> ISM envelope lower ﬁducial radius
     self%myconfig%profile_vd    = 0.0_dp   !> ISM envelope kepler velocity at r=profile_rd
     self%myconfig%profile_rho_0 = 0.0_dp  !> ISM envelope fiducial density

     self%myconfig%reset_coef            = 0.0_dp
     self%myconfig%weight_mean           = 'arithmetic'
     self%myconfig%weight_mean_variable  = 'none'
     self%myconfig%weight_analy          = 0.0_dp
     self%myconfig%weight_num            = 1.0_dp
     self%myconfig%weight_analy_index    = 1.0_dp
     self%myconfig%weight_num_index      = 1.0_dp
     self%myconfig%reset_dtemperature    = 0.0_dp
     self%myconfig%reset_distance        = 0.0_dp
     self%myconfig%reset_scale           = 0.0_dp
     self%myconfig%reset_on              = .false.
     self%myconfig%boundary_on           = .false.
     self%myconfig%boundary_cond         = 'open'
     self%myconfig%flux_frac             = 1.0d-2
     self%myconfig%debug                 = .false.
     self%myconfig%mixed_fixed_bound     = .false.
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
    real(dp)                 :: mp,kb,Ggrav,r_normalized,costhetazero
    integer                  :: idim,iside
    !-----------------------------------

    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
    else
      mp=mp_cgs
      kB=kB_cgs
    end if
    Ggrav = constusr%G


    !----------------------------------------------------
    ! set the parameters :
    !   - mean_nall_to_nH = n_e/n_i
    !   - mean_mass = rho_tot/(nH*mH), where nH=n(H+H^+)+2n(H_2)
    !   - mean_mup,mean_ne_to_nH = n(e^-)/nH
    ! from calling the mod_hd_phys/mod_mhd_phys module
    !----------------------------------------------------

    !write(*,*) ' mod_obj_ism.t--> usr_ism_set_complet'
    !write(*,*) 'self (before change): chemical_composition | ',&
    !'xHe | mean_nall_to_nH |',&
    !' mean_mass | mean_ne_to_nH | mean_mup '
    !write(*,*) self%myconfig%chemical_gas_type,&
    !self%myconfig%He_abundance,&
    !self%myconfig%mean_nall_to_nH,&
    !self%myconfig%mean_mass,&
    !self%myconfig%mean_ne_to_nH,&
    !self%myconfig%mean_mup

    call phys_fill_chemical_ionisation(self%myconfig%He_abundance,  &
       self%myconfig%chemical_gas_type,                             &
       self%myconfig%mean_nall_to_nH,self%myconfig%mean_mass,       &
       self%myconfig%mean_mup,self%myconfig%mean_ne_to_nH)

     !write(*,*) ' mod_obj_ism.t--> usr_ism_set_complet'
     !write(*,*) 'self (after change): chemical_composition | ',&
     !'xHe | mean_nall_to_nH |',&
     !' mean_mass | mean_ne_to_nH | mean_mup '
     !write(*,*) self%myconfig%chemical_gas_type,&
     !self%myconfig%He_abundance,&
     !self%myconfig%mean_nall_to_nH,&
     !self%myconfig%mean_mass,&
     !self%myconfig%mean_ne_to_nH,&
     !self%myconfig%mean_mup

   !----------------------------------------------------
   ! set the density and/or number density parameters accordingly
   ! (provided density takes priority)
   ! 31-01-22 : density successfully implemented
   !----------------------------------------------------

    !if we do Ulrich model
    if(trim(self%myconfig%profile_density)=='Ulrich1976')then
      if(mype==0)then
        write(*,*) '====================Ulrich1976s model===================='
        write(*,*) '-Setting completely all parameters:'
      end if
      ! rd
      self%myconfig%profile_rd = (self%myconfig%profile_jspec**2.0_dp)/&
      (Ggrav*self%myconfig%profile_Mstar)
      if(mype==0)then
        write(*,*) '* rW = ', self%myconfig%profile_rw, ' cm'
        write(*,*) '* rD = ', self%myconfig%profile_rd, ' cm'
      end if
      ! v_Kepler
      self%myconfig%profile_vd = dsqrt((Ggrav*self%myconfig%profile_Mstar)/&
      self%myconfig%profile_rd)
      if(mype==0)then
        write(*,*) '* vK = ', self%myconfig%profile_vd, ' cm/s'
      end if
      ! rho_0
      self%myconfig%profile_rho_0 = self%myconfig%profile_Minf/(4.0_dp*&
      dpi*(self%myconfig%profile_rd**2.0_dp)*self%myconfig%profile_vd)
      if(mype==0)then
        write(*,*) '* rho_0 = ', self%myconfig%profile_rho_0, ' g/cm3'
        r_normalized = (self%myconfig%profile_rw/self%myconfig%profile_rd)
        write(*,*) '* r = (rD/rw) = ', r_normalized
      end if
      if(dabs(r_normalized-1.0_dp)<smalldouble)then
        costhetazero = 0.0_dp
      else if (r_normalized>1.0_dp)then
        costhetazero = 0.0_dp
      else if (r_normalized<1.0_dp&
      .and.(-((1.0_dp-r_normalized)/3.0_dp)**3.0_dp)>0.0_dp)then
        costhetazero = dsqrt(1.0_dp-r_normalized)
      else if (r_normalized<1.0_dp&
      .and.(-((1.0_dp-r_normalized)/3.0_dp)**3.0_dp)<0.0_dp)then
        costhetazero = dsqrt(1.0_dp-r_normalized)
      end if
      ! rho_amb(r=rw,theta=0.5*pi)=rho_a0
      self%myconfig%density = self%myconfig%profile_rho_0*&
      (r_normalized**(-1.5_dp))*(1.0_dp+(2.0_dp/r_normalized)*&
      0.5_dp*(3*costhetazero*costhetazero-1.0_dp))**(-1.0_dp)

      if(mype==0)then
        write(*,*) '* rho_a0 = ', self%myconfig%density, ' g/cm3'
      end if
      self%myconfig%number_density= self%myconfig%density/(mp*self%myconfig%mean_mass)
      if(mype==0)then
        write(*,*) '* nH_a0 = ', self%myconfig%number_density, ' cm-3'
        write(*,*) '========================================================='
      end if
    else
      if (dabs( self%myconfig%density)<smalldouble*mp)then
           self%myconfig%density= self%myconfig%number_density*mp*self%myconfig%mean_mass
      elseif(dabs( self%myconfig%number_density)<smalldouble)then
           self%myconfig%number_density= self%myconfig%density/(mp*self%myconfig%mean_mass)
      end if
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
       if(trim(self%myconfig%profile_density)=='Ulrich1976')then
         self%myconfig%pressure_Ulrich = self%myconfig%c_sound**2.0_dp *&
                        self%myconfig%profile_rho_0 / phys_config%gamma
       end if
       self%myconfig%temperature = self%myconfig%c_sound**2.0_dp*mp*&
                                self%myconfig%mean_mup/(kB*phys_config%gamma)
    !----------------------------------------------------
    ! *If not yet provided or input,
    ! *set the pressure parameter of the ambient medium
    ! *or the ism temperature.
    ! *At the end, we can
    ! *assign the sound speed accordingly if not yet set
    ! *19-03-21 : made sure that temperature isn t changed if provided
    ! * and takes priority over pressure, which can be changed
    !----------------------------------------------------
    else  cond_csound_set
     if(self%myconfig%temperature>0.0_dp) then
      self%myconfig%pressure = self%myconfig%density*&
                            kB* self%myconfig%temperature/(self%myconfig%mean_mup*mp)
      if(trim(self%myconfig%profile_density)=='Ulrich1976')then
        self%myconfig%pressure_Ulrich = self%myconfig%profile_rho_0*&
                              kB* self%myconfig%temperature/(self%myconfig%mean_mup*mp)
      end if
     else if(dabs(self%myconfig%pressure)>0.0_dp) then
      self%myconfig%temperature = self%myconfig%mean_mup*mp*self%myconfig%pressure&
                            /(kB*self%myconfig%density)
     end if
     if(trim(self%myconfig%profile_density)=='Ulrich1976')then
       self%myconfig%c_sound = sqrt(phys_config%gamma*self%myconfig%pressure/&
            self%myconfig%profile_rho_0)
     else
        self%myconfig%c_sound = sqrt(phys_config%gamma*self%myconfig%pressure/self%myconfig%density)
     end if
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
         if(trim(self%myconfig%profile_density)=='Ulrich1976')then
           self%myconfig%profile_shift_density(idim,iside) = self%myconfig%profile_rho_0
         else
           self%myconfig%profile_shift_density(idim,iside) = self%myconfig%density
         end if

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
    ! original :
    ! if(.not.self%myconfig%boundary_on)                                       &
    !    self%myboundaries%myconfig%boundary_type=self%myconfig%boundary_cond


    ! Here the boundary conditions input from ism_config takes priority from what is
    ! input from usrboundary_config
    if(self%myconfig%boundary_on)then
        !write(*,*) 'self%myconfig%boundary_on is .true.'
        !write(*,*) '- Before the prioritization of boundary conditions :'
        !write(*,*) 'we have self%myboundaries%myconfig%boundary_type = ', self%myboundaries%myconfig%boundary_type
        !write(*,*) 'we have self%myconfig%boundary_cond = ', self%myconfig%boundary_cond
        self%myboundaries%myconfig%boundary_type=self%myconfig%boundary_cond
        self%myboundaries%myconfig%flux_frac=self%myconfig%flux_frac
        !write(*,*) '- After the prioritization of boundary conditions :'
        !write(*,*) 'we have self%myboundaries%myconfig%boundary_type = ', self%myboundaries%myconfig%boundary_type
        !write(*,*) 'we have self%myconfig%boundary_cond = ', self%myconfig%boundary_cond
    else
      !write(*,*) 'self%myconfig%boundary_on is .false.'
      !write(*,*) 'we have self%myboundaries%myconfig%boundary_type = ', self%myboundaries%myconfig%boundary_type
      !write(*,*) 'we have self%myconfig%boundary_cond = ', self%myconfig%boundary_cond
    end if

    call self%myboundaries%set_complet

    if(mype==0)then
      write(*,*) '======================ISM settings======================='
      print*, '* at (R,z) = (r_j,0) : '
      print*, '** The ism base temperature is = ',self%myconfig%temperature
      print*, '** The ism base density is =', self%myconfig%density
      print*, '** The ism base number density nH is ', self%myconfig%number_density
      write(*,*) '========================================================='
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

    !write(*,*) ' mod_obj_ism.t--> usr_ism_normalize'
    !write(*,*) 'self (before normalization) : chemical_composition | ',&
    !'xHe | mean_nall_to_nH |',&
    !' mean_mass | mean_ne_to_nH | mean_mup '
    !write(*,*) self%myconfig%chemical_gas_type,&
    !self%myconfig%He_abundance,&
    !self%myconfig%mean_nall_to_nH,&
    !self%myconfig%mean_mass,&
    !self%myconfig%mean_ne_to_nH,&
    !self%myconfig%mean_mup

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
     self%myconfig%pressure_Ulrich    =  self%myconfig%pressure_Ulrich /physunit_inuse%myconfig%pressure
     self%myconfig%velocity           =  self%myconfig%velocity      /physunit_inuse%myconfig%velocity
     self%myconfig%velocity_ofmov_obj =  self%myconfig%velocity_ofmov_obj  /physunit_inuse%myconfig%velocity
     self%myconfig%extend             =  self%myconfig%extend        /physunit_inuse%myconfig%length
     self%myconfig%c_sound            =  self%myconfig%c_sound       /physunit_inuse%myconfig%velocity
     self%myconfig%mean_mup           =  self%myconfig%mean_mup      / physunit_inuse%myconfig%mean_mup

     self%myconfig%theta_floored = self%myconfig%theta_floored *(dpi/180._dp)
     self%myconfig%profile_rw            = self%myconfig%profile_rw/physunit_inuse%myconfig%length
     self%myconfig%profile_rw_force      = self%myconfig%profile_rw_force/physunit_inuse%myconfig%length
     self%myconfig%profile_zc              =  self%myconfig%profile_zc           /physunit_inuse%myconfig%length
     self%myconfig%profile_center              =  self%myconfig%profile_center           /physunit_inuse%myconfig%length
     !> Ulrich model ambient medium parameters
     self%myconfig%profile_jspec = self%myconfig%profile_jspec / (physunit_inuse%myconfig%velocity*physunit_inuse%myconfig%length)
     self%myconfig%profile_Minf  = self%myconfig%profile_Minf  / (physunit_inuse%myconfig%mass/physunit_inuse%myconfig%time)
     self%myconfig%profile_Mstar = self%myconfig%profile_Mstar / physunit_inuse%myconfig%mass
     self%myconfig%profile_rd    = self%myconfig%profile_rd    / physunit_inuse%myconfig%length
     self%myconfig%profile_vd    = self%myconfig%profile_vd    / physunit_inuse%myconfig%velocity
     self%myconfig%profile_rho_0 = self%myconfig%profile_rho_0 / physunit_inuse%myconfig%density

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

    !write(*,*) ' mod_obj_ism.t--> usr_ism_normalize'
    !write(*,*) 'self (after normalization): chemical_composition | ',&
    !'xHe | mean_nall_to_nH |',&
    !' mean_mass | mean_ne_to_nH | mean_mup '
    !write(*,*) self%myconfig%chemical_gas_type,&
    !self%myconfig%He_abundance,&
    !self%myconfig%mean_nall_to_nH,&
    !self%myconfig%mean_mass,&
    !self%myconfig%mean_ne_to_nH,&
    !self%myconfig%mean_mup

   end subroutine usr_ism_normalize
  !--------------------------------------------------------------------
   !> subroutine setting for ISM
   subroutine usr_ism_set_w(ixI^L,ixO^L,qt,x,w,self,isboundary_iB)
    implicit none
    integer, intent(in)           :: ixI^L,ixO^L
    real(kind=dp), intent(in)     :: qt
    real(kind=dp), intent(in)     :: x(ixI^S,1:ndim)
    real(kind=dp)                 :: w(ixI^S,1:nw)
    real(kind=dp),allocatable     :: w_tmp(:^D&,:)
    integer,             optional :: isboundary_iB(2)
    class(ism)                    :: self
    ! .. local..
    integer                    :: idir,IB,idims,idims_bound,iw,iwfluxbc,idims2,iside2
    real(kind=dp)              :: fprofile(ixI^S)
    logical                    :: isboundary,to_fix,bc_to_fix,some_unfixed
    character(len=30)          :: myboundary_cond
    character(len=64)          :: aString
    !----------------------------------



    if(.not.allocated(self%patch)) then
     allocate(self%patch(ixI^S))
     self%patch              = .true.
    end if
    cond_B_present : if(present(isboundary_iB))then
    !write(*,*) '>1)isiB1= isib2=', isboundary_iB(1),' ', isboundary_iB(2)
    !write(*,*) '>boundary_cond = ', self%myconfig%boundary_cond(isboundary_iB(1),isboundary_iB(2))



      {^D&
        idims=^D
        if(isboundary_iB(1)==idims.and.isboundary_iB(2)==1) then
         if(all(x(ixO^S,idims)<=xprobmin^D))then
           isboundary=.true.
           ! original : IB = 2*(idims-1)+1
           ! 03-03-2022 : we defined idims=(iB+1)/2 in mod_usr_yso_jet.t
           ! ==> 2*idims = iB +1 ==> iB = 2*idims-1= 2*idims - 2 +1
           ! ==> iB = 2*(idims-1)+1 which is identical to the original
           IB = 2*(idims-1)+1
           idims_bound = idims
         end if
        end if
        if(isboundary_iB(1)==idims.and.isboundary_iB(2)==2) then
         if(all(x(ixO^S,idims)>=xprobmax^D))then
           isboundary=.true.
           ! original : IB = 2*idims
           ! 03-03-2022 : we defined idims=iB/2 in mod_usr_yso_jet.t
           ! ==> iB = 2*idims which is identical to the original
           IB = 2*idims
           idims_bound = idims
         end if
        end if
      \}

      !if(mype==0)then
      !  write(*,*) '*****Conditions for boundary number iB=', iB,'*****'
      !  write(*,*) '* ...corresponding to (idims,iside)=', isboundary_iB(1),isboundary_iB(2)
      !end if

      myboundary_cond = self%myconfig%boundary_cond(isboundary_iB(1),isboundary_iB(2))

      !if(mype==0)then
      !  write(*,*) '* ...and which boundary condition is =', myboundary_cond
      !  write(*,*) '************************************************************'
        !write(*,*) 'present myboundary_cond=', myboundary_cond
      !end if
    else cond_B_present
      isboundary=.false.
      myboundary_cond = 'notbc_but_fix'
      !write(*,*) 'not present myboundary_cond=', myboundary_cond
    end if cond_B_present



    debug_mode : if(.not.self%myconfig%debug)then
    !write(*,*) 'Not debug mode'
    if(.not.isboundary.or.&
                       trim(myboundary_cond)=='fix') then

      !write(*,*) 'call it toooo here !!!!!!!!!!!!!!!!!!!!!!!!'
      if(trim(self%myconfig%profile_density)=='Ulrich1976')then
        where(self%patch(ixO^S))
            w(ixO^S,phys_ind%rho_)        =  self%myconfig%profile_rho_0
        end where
      else
        where(self%patch(ixO^S))
            w(ixO^S,phys_ind%rho_)        =  self%myconfig%density
        end where
      end if
      if(phys_config%energy)then
        if(trim(self%myconfig%profile_density)=='Ulrich1976')then
          where(self%patch(ixO^S))
              w(ixO^S,phys_ind%pressure_)   =  self%myconfig%pressure_Ulrich
          end where
        else
          where(self%patch(ixO^S))
              w(ixO^S,phys_ind%pressure_)   =  self%myconfig%pressure
          end where
        end if
      end if


      if(trim(self%myconfig%profile_density)=='Ulrich1976')then
         !----------------------------------------------------
         !line to modify if we want to change the velocity
         !put for each time step at any boundary
         !set to 'fix'
         !----------------------------------------------------
         do idir=1,ndir
           where(  self%patch(ixO^S))
             !----------------------------------------------------
             !line to modify if we want to change the velocity
             !put for each time step at any boundary
             !set to 'fix'
             !----------------------------------------------------
              w(ixO^S,phys_ind%mom(idir)) = self%myconfig%profile_vd
           end where
         end do
      else
         do idir=1,ndir
         where(  self%patch(ixO^S))
         !----------------------------------------------------
         !line to modify if we want to change the velocity
         !put for each time step at any boundary
         !set to 'fix'
         !----------------------------------------------------
          w(ixO^S,phys_ind%mom(idir)) = self%myconfig%velocity(idir)
         end where
        end do
      end if


      ! add profile to ISM density and pressure
      if(self%myconfig%profile_on) then
         call self%set_profile(ixI^L,ixO^L,x,w)
      end if




      if(phys_config%mean_mup_on) then
        where(self%patch(ixO^S))
         w(ixO^S,phys_ind%mup_) = self%myconfig%mean_mup
        end where
      end if

      if(self%myconfig%tracer_on.and.phys_config%n_tracer>0&
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
      end if



      if( self%myconfig%dust_on)then
        call self%mydust%set_patch(ixI^L,ixO^L,self%patch)
        self%mydust%myconfig%velocity= self%myconfig%velocity
        fprofile = 1.0_dp
        call   self%mydust%set_w(ixI^L,ixO^L,qt,.false., self%myconfig%dust_frac,fprofile,x,w)
      end if
    else
     if(any(self%patch(ixO^S))) then
       call self%myboundaries%set_w(ixI^L,ixO^L,iB,isboundary_iB(1),isboundary_iB(2),&
                                self%patch,x,w)
     end if
    end if

    else debug_mode
    ! write(*,*) 'Debug mode'
    ! Start of the conditional loop when we treat the domain without the boundaries
    ! or when my_boundary=='fix' : to modify in order
    ! to change what 'fix' do
    ! + Start of  the loop for any my_boundary is other than 'fix' and
    ! mixed_fixed_bound comport a flux variable which is set to true : similar to what
    ! my_boundary=='fix' does but only along one flux index (e.g. v_phi)


    if(isboundary.and.(trim(myboundary_cond)=='fix')) then
      self%myconfig%mixed_fixed_bound(isboundary_iB(1),isboundary_iB(2),1:7)=.true.
    end if
    !if(present(isboundary_iB))then
    !write(*,*) '>2)isiB1= isib2=', isboundary_iB(1),' ', isboundary_iB(2)
    !write(*,*) '>mixed_fixed_bound = ', self%myconfig%mixed_fixed_bound(isboundary_iB(1),isboundary_iB(2),1:7)
    !end if

    !if we want to fix the inner domain or if we need to fix all variables in the boundaries...
    to_fix = ((.not.isboundary).or.(trim(myboundary_cond)=='fix'))

    if(isboundary) then
      !...or if we have boundary conditions such that one variable must be fixed in those boundaries
      !write(*,*) '>mixed_fixed_bound = ', self%myconfig%mixed_fixed_bound(isboundary_iB(1),isboundary_iB(2),1:7)
      !write(*,*) 'any(self%myconfig%mixed_fixed_bound(isboundary_iB(1),isboundary_iB(2),1:7))'
      !write(*,*) '=', any(self%myconfig%mixed_fixed_bound(isboundary_iB(1),isboundary_iB(2),1:7))
      to_fix = to_fix.or.((isboundary).and.any(self%myconfig%mixed_fixed_bound(isboundary_iB(1),isboundary_iB(2),1:7)))
      !write(*,*) '>mixed_fixed_bound = ', self%myconfig%mixed_fixed_bound(isboundary_iB(1),isboundary_iB(2),1:7)
      !write(*,*) 'any(.not.self%myconfig%mixed_fixed_bound(isboundary_iB(1),'
      !write(*,*) 'isboundary_iB(2),1:7))'
      !write(*,*) '=', any(.not.self%myconfig%mixed_fixed_bound(isboundary_iB(1),&
      !isboundary_iB(2),1:7))
      some_unfixed = any(.not.self%myconfig%mixed_fixed_bound(isboundary_iB(1),&
      isboundary_iB(2),1:7))
    end if



    boundary_cond : if(to_fix) then

      !1) First, set all boundary conditions related to
      ! (isboundary_iB(1),isboundary_iB(2))
      ! = (/idims,iside/) according to what
      ! my_boundary = self%myconfig%boundary_cond(isboundary_iB(1),isboundary_iB(2))
      ! tells
      ! <=> If some boundary conditions are not 'fixed'
      if(isboundary.and.some_unfixed)then
        !write(*,*) '>>some bc are unfixed'
        !write(*,*) '> mixed_fixed_bound = ', self%myconfig%mixed_fixed_bound(isboundary_iB(1),isboundary_iB(2),1:7)
        do idims2=1,ndim
          do iside2=1,2
            !write(*,*) '>idims=', idims2
            !write(*,*) '>iside=', iside2
            !write(*,*) 'self%myconfig%boundary_cond = ',self%myconfig%boundary_cond(idims2,iside2)
            !write(*,*) 'self%myboundaries%myconfig%boundary_type = ',self%myboundaries%myconfig%boundary_type(idims2,iside2)
          end do
        end do
        if(any(self%patch(ixO^S))) then
          call self%myboundaries%set_w(ixI^L,ixO^L,iB,isboundary_iB(1),isboundary_iB(2),&
                                  self%patch,x,w)
        end if
      end if

      bc_to_fix = .false.

      !2) Then, in any cases,
      ! make sure that the variable w(ixO^S,1:7) are fixed accordingly:

      if(isboundary)bc_to_fix = self%myconfig%mixed_fixed_bound(isboundary_iB(1),isboundary_iB(2),phys_ind%rho_)
      ! initialize uniform inner grid or fix density in boundaries
      update_rho : if((.not.isboundary).or.bc_to_fix)then
        if(trim(self%myconfig%profile_density)=='Ulrich1976')then
          where(self%patch(ixO^S))
              w(ixO^S,phys_ind%rho_)        =  self%myconfig%profile_rho_0
          end where
        else
          where(self%patch(ixO^S))
              w(ixO^S,phys_ind%rho_)        =  self%myconfig%density
          end where
        end if
      end if update_rho

      if(isboundary)bc_to_fix = self%myconfig%mixed_fixed_bound(isboundary_iB(1),isboundary_iB(2),phys_ind%pressure_)
      ! if energy is evolved: initialize uniform inner grid or fix pressure in boundaries
      if(phys_config%energy)then
        update_p : if((.not.isboundary).or.bc_to_fix)then
          if(trim(self%myconfig%profile_density)=='Ulrich1976')then
            where(self%patch(ixO^S))
                w(ixO^S,phys_ind%pressure_)   =  self%myconfig%pressure_Ulrich
            end where
          else
            where(self%patch(ixO^S))
                w(ixO^S,phys_ind%pressure_)   =  self%myconfig%pressure
            end where
          end if
        end if update_p
      end if

      ! initialize uniform inner grid or fix velocities in boundaries
      Loop_idir : do idir=1,ndir
        if(isboundary)bc_to_fix = self%myconfig%mixed_fixed_bound(isboundary_iB(1),isboundary_iB(2),phys_ind%mom(idir))
        update_mom : if((.not.isboundary).or.bc_to_fix)then
          if(trim(self%myconfig%profile_density)=='Ulrich1976')then

               where(  self%patch(ixO^S))
                 !----------------------------------------------------
                 !lines to modify if we want to change the velocity
                 !put for each time step at any boundary
                 !set to 'fix' with Ulrich1976 model
                 !----------------------------------------------------
                  w(ixO^S,phys_ind%mom(idir)) = self%myconfig%profile_vd
               end where

          else

             where(  self%patch(ixO^S))
               !----------------------------------------------------
               !lines to modify if we want to change the velocity
               !put for each time step at any boundary
               !set to 'fix'
               !----------------------------------------------------
                w(ixO^S,phys_ind%mom(idir)) = self%myconfig%velocity(idir)
             end where

          end if
        end if update_mom
      end do Loop_idir


      ! add profile to ISM density, pressure and eventually velocity
      ! in the inner grid or to boundaries
      if(isboundary)bc_to_fix = any(self%myconfig%mixed_fixed_bound(isboundary_iB(1),isboundary_iB(2),1:7))
      treat_profile: if((.not.isboundary).or.bc_to_fix)then
        add_profile : if(self%myconfig%profile_on) then
          !when treating the inner domain
          !or in case 'fix' conditions is called at every bc flux variables
          cond_fix : if((.not.isboundary).or.trim(myboundary_cond)=='fix')then
            !apply the profile to rho,p and mom(:)
            call self%set_profile(ixI^L,ixO^L,x,w)

            ! when treating the boundaries of the domain
            ! in case any other bc conditions is called and that one of the
            ! nwfluxbc flux bc variables must be fixed
          else if(isboundary.and.&
          any(self%myconfig%mixed_fixed_bound(isboundary_iB(1),isboundary_iB(2),1:7)))then cond_fix

            !make sure to work with a fresh untouched wtmp version
            if(.not.allocated(w_tmp)) then
              allocate(w_tmp(ixI^S,1:nw))
            else
              deallocate(w_tmp)
              allocate(w_tmp(ixI^S,1:nw))
            end if

            !store the uniformly initialized w into a separate w_tmp
            w_tmp(ixO^S,1:nw) = w(ixO^S,1:nw)
            !apply the profile on w_tmp only
            call self%set_profile(ixI^L,ixO^L,x,w_tmp)

            ! modify only the flux variable to be fixed in w

            do iwfluxbc=1,nwfluxbc
              if(self%myconfig%mixed_fixed_bound(isboundary_iB(1),isboundary_iB(2),iwfluxbc))then
                where(self%patch(ixO^S))
                  w(ixO^S,iwfluxbc)=w_tmp(ixO^S,iwfluxbc)
                end where
              end if
            end do

            !don t forget to free memory :
            deallocate(w_tmp)
          else cond_fix
            !write(*,*) 'Unknown case of bc conditions or inner domain ...'
            call mpistop(' ... initialization at cond_fix')
          end if cond_fix
        end if add_profile
      end if treat_profile

      ! End of the conditional loop for initialization of the inner domaine
      ! or when my_boundary=='fix' or when mixed 'fix' boundary conditions are input
      ! Start of  the loop for any other my_boundary : call to set_w which
      ! set the conditions in the boundaries from mod_obj_mat.t
      else if(isboundary.and.(trim(myboundary_cond)/='fix'))then boundary_cond
        !write(*,*) '>> here is the boundary (idims,iside) = ',isboundary_iB(1),isboundary_iB(2),'...'
        !write(*,*) '>> ...where no bc are fixed but we get the condition : ', myboundary_cond
        if(any(self%patch(ixO^S))) then
          call self%myboundaries%set_w(ixI^L,ixO^L,iB,isboundary_iB(1),isboundary_iB(2),&
                                  self%patch,x,w)
        end if
      else boundary_cond
        write(*,*) 'Unknown case of bc conditions or inner domain...'
        call mpistop('...initialization at boundary_cond')
    end if boundary_cond


    !Treat mean_mup the same way everywhere :
    if(phys_config%mean_mup_on) then
      where(self%patch(ixO^S))
       w(ixO^S,phys_ind%mup_) = self%myconfig%mean_mup
      end where
    end if

    !Treat tracers the same way everywhere:
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

    !Treat dust the same way  everywhere
    cond_dust_on : if( self%myconfig%dust_on)then
      call self%mydust%set_patch(ixI^L,ixO^L,self%patch)
      self%mydust%myconfig%velocity= self%myconfig%velocity
      fprofile = 1.0_dp
      call   self%mydust%set_w(ixI^L,ixO^L,qt,.false., self%myconfig%dust_frac,fprofile,x,w)
    end if cond_dust_on

    end if debug_mode
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







   !--------------------------------------------------------------------
   subroutine usr_ism_set_profile(ixI^L, ixO^L,x,w,self)
    implicit none
    integer, intent(in)             :: ixI^L,ixO^L
    real(kind=dp), intent(in)       :: x(ixI^S,1:ndim)
    real(kind=dp), intent(inout)    :: w(ixI^S,1:nw)
    class(ism)                      :: self
    ! ..local ..
    real(kind=dp), dimension(ixI^S) :: p_profile,d_profile,theta_profile
    real(kind=dp), dimension(ixI^S) :: vr_profile,vt_profile,vp_profile
    real(kind=dp), dimension(ixI^S) :: cos_theta_zero,sin_theta_zero,r_normalized
    real(kind=dp), dimension(ixI^S) :: theta_zero
    real(kind=dp), dimension(ixI^S,1:ndir) :: project_speed
    character(len=30)               :: cos_or_sin
    integer                         :: idir
    !----------------------------------------------------

    call self%set_profile_distance(ixI^L,ixO^L,x,d_profile)

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

              call usr_get_theta(ixI^L,ixO^L,x,theta_profile)

              p_profile(ixO^S) = ((dsin(theta_profile(ixO^S)))**2.0_dp)/&
                                 ((d_profile(ixO^S)/self%myconfig%profile_rw)**2.0_dp)
              !print*,p_profile(ixO^S)
          if(trim(self%myconfig%profile_density)=='Lee2001_floored')then
            where(theta_profile(ixO^S)<=self%myconfig%theta_floored)
            p_profile(ixO^S) = ((dsin(self%myconfig%theta_floored))**2.0_dp)/&
                               ((d_profile(ixO^S)/self%myconfig%profile_rw)**2.0_dp)
            end where
          end if

        case('Ulrich1976')

          call self%set_ulrich_profile(ixI^L, ixO^L,x,w,p_profile)

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


   end subroutine usr_ism_set_profile

   !--------------------------------------------------------------------
   subroutine usr_set_ulrich_profile(ixI^L, ixO^L,x,w,p_profile,self)
   implicit none
   integer, intent(in)             :: ixI^L,ixO^L
   real(kind=dp), intent(in)       :: x(ixI^S,1:ndim)
   real(kind=dp), intent(inout)    :: w(ixI^S,1:nw)
   real(kind=dp), intent(inout)    :: p_profile(ixI^S)
   class(ism)                      :: self
   ! ..local ..
   real(kind=dp), dimension(ixI^S) :: d_profile,theta_profile
   real(kind=dp), dimension(ixI^S) :: vr_profile,vt_profile,vp_profile
   real(kind=dp), dimension(ixI^S) :: cos_theta_zero,sin_theta_zero,r_normalized
   real(kind=dp), dimension(ixI^S) :: theta_zero
   real(kind=dp), dimension(ixI^S,1:ndir) :: project_speed
   integer                         :: idir
   real(kind=dp), dimension(1:ndim) :: zero_dim
   !--------------------------------------------------------------------


   {zero_dim(^D)=0.0_dp\}!-dx(1,1)

   p_profile(ixO^S) = 1.0_dp
   vr_profile(ixO^S) = 1.0_dp
   vt_profile(ixO^S) = 1.0_dp
   vp_profile(ixO^S) = 1.0_dp

   call self%set_profile_distance(ixI^L,ixO^L,x,d_profile)

   call usr_get_theta(ixI^L,ixO^L,x,theta_profile)

   r_normalized(ixI^S) = (d_profile(ixI^S)/self%myconfig%profile_rd)

   call usr_ulrich1976_costheta_zero(ixI^L, ixO^L,x,r_normalized,&
   theta_profile,cos_theta_zero)

   theta_zero(ixI^S) = DACOS(cos_theta_zero(ixI^S))

   p_profile(ixO^S) = (r_normalized(ixO^S)**(-1.5_dp))*&
   ((1.0_dp + (DCOS(theta_profile(ixO^S))/&
   cos_theta_zero(ixO^S)))**(-0.5_dp))/&
   ((1.0_dp + (3.0_dp * cos_theta_zero(ixO^S) * &
   cos_theta_zero(ixO^S) - 1.0_dp)/r_normalized(ixO^S)))

   vr_profile(ixO^S) = -(r_normalized(ixO^S)**(-0.5_dp))*&
   ((1.0_dp + (DCOS(theta_profile(ixO^S))/cos_theta_zero(ixO^S)))**(0.5_dp))

   vt_profile(ixO^S) = (r_normalized(ixO^S)**(-0.5_dp))*&
   ((1.0_dp + (DCOS(theta_profile(ixO^S))/cos_theta_zero(ixO^S)))**(0.5_dp))*&
   ((cos_theta_zero(ixO^S)-DCOS(theta_profile(ixO^S)))/DSIN(theta_profile(ixO^S)))

   vp_profile(ixO^S) = (r_normalized(ixO^S)**(-0.5_dp))*&
   ((1.0_dp - (DCOS(theta_profile(ixO^S))/cos_theta_zero(ixO^S)))**(0.5_dp))*&
   ((DSIN(theta_zero(ixO^S)))/DSIN(theta_profile(ixO^S)))


   select case(typeaxial)
     case('spherical')
        where(self%patch(ixO^S))
          project_speed(ixO^S,r_) = vr_profile(ixO^S)
        end where
        if(ndir>1)then
          where(self%patch(ixO^S))
            project_speed(ixO^S,theta_) = vt_profile(ixO^S)
          end where
        end if
        if(ndir>2)then
          where(self%patch(ixO^S))
            project_speed(ixO^S,phi_) = vp_profile(ixO^S)
          end where
        end if
     case('cylindrical')
        if(ndim<3)then
          where(self%patch(ixO^S))
            project_speed(ixO^S,r_) = (vr_profile(ixO^S)*&
            DSIN(theta_profile(ixO^S)))+(vt_profile(ixO^S)*&
            DCOS(theta_profile(ixO^S)))
            project_speed(ixO^S,z_) = (vr_profile(ixO^S)*&
            DCOS(theta_profile(ixO^S)))-(vt_profile(ixO^S)*&
            DSIN(theta_profile(ixO^S)))
          end where
        else
          where(self%patch(ixO^S))
            project_speed(ixO^S,r_) = 1.0_dp ! TO DO
            project_speed(ixO^S,z_) = 1.0_dp ! TO DO
          end where
        end if
        if(ndir>2)then
           where(self%patch(ixO^S))
             project_speed(ixO^S,phi_) = vp_profile(ixO^S)
           end where
        end if
     case('slab','slabstretch')
        if(ndim<3)then
          where(self%patch(ixO^S))
            project_speed(ixO^S,x_) = (vr_profile(ixO^S)*&
            DSIN(theta_profile(ixO^S)))+(vt_profile(ixO^S)*&
            DCOS(theta_profile(ixO^S)))
            project_speed(ixO^S,y_) = (vr_profile(ixO^S)*&
            DCOS(theta_profile(ixO^S)))-(vt_profile(ixO^S)*&
            DSIN(theta_profile(ixO^S)))
          end where
          if(ndir>2)then
             where(self%patch(ixO^S))
               project_speed(ixO^S,z_) = vp_profile(ixO^S)
             end where
          end if
        else
          where(self%patch(ixO^S))
            project_speed(ixO^S,x_) = 1.0_dp ! TO DO
            project_speed(ixO^S,y_) = 1.0_dp ! TO DO
          end where
          if(ndir>2)then
             where(self%patch(ixO^S))
               project_speed(ixO^S,z_) = 1.0_dp ! TO DO
             end where
          end if
        end if
     case default
        call mpistop('Unknown typeaxial')
   end select



       Loop_idir_profile : do idir=1,ndir
        where(self%patch(ixO^S))
        !----------------------------------------------------
        !line to modify if Ulrich model and we want to change the velocity
        !put for each time step at any boundary
        !set to 'fix'
        !----------------------------------------------------
         w(ixO^S,phys_ind%mom(idir)) = w(ixO^S,phys_ind%mom(idir)) * &
         project_speed(ixO^S,idir)
        end where
       end do Loop_idir_profile

       end  subroutine usr_set_ulrich_profile


   !--------------------------------------------------------------------

   subroutine usr_ism_get_pforce_profile(ixI^L,ixO^L,qt,qdt,x,w,f_profile,divpv,self)
    use mod_radiative_cooling
    implicit none
    integer, intent(in)                  :: ixI^L,ixO^L
    real(kind=dp), intent(in)            :: qt,qdt
    real(kind=dp)                        :: maxv,minv
    real(kind=dp), intent(in)            :: x(ixI^S,1:ndim)
    real(kind=dp), intent(in)            :: w(ixI^S,1:nw)
    class(ism)                           :: self
    real(kind=dp), intent(inout)         :: f_profile(ixI^S,1:ndim)
    real(kind=dp)                        :: f_profile1(ixI^S,1:ndim)
    real(kind=dp)                        :: f_profile2(ixI^S,1:ndim)
    real(kind=dp), intent(inout)         :: divpv(ixI^S)
    ! ..local ..
    integer                              :: idims
    real(kind=dp), dimension(ixI^S,1:ndir) :: pv
    real(kind=dp), dimension(ixI^S,1:nw) :: w_init, w_tmp
    real(kind=dp), dimension(ixI^S,1:nw) :: w_init2, w_tmp2
    real(kind=dp), dimension(ixI^S)      :: ptherm,gradp,d_profile,theta_profile
    real(kind=dp), dimension(ixI^S)      :: alpha_analy,alpha_num,fsign
    logical                              :: tmp_active,src_active
    logical, dimension(ixI^S)            :: patch_tosave
    !----------------------------------------------------
    f_profile(ixO^S,:) = 0.0
    cond_fprofile : if(self%myconfig%profile_force_gradP_on)then
      cond_density_fprofile : if(self%myconfig%profile_density_on) then

         !1> Compute analytical force
         call self%set_profile_distance(ixI^L,ixO^L,x,d_profile)
         !initiate the ism positions at current time step
         !and the primitive state vector w_init, all previously to
         !radiative cooling treatment

         patch_tosave(ixI^S)=self%patch(ixI^S)
         ! allow the following lines to be applied to the whole input
         ! part of the grid in self%patch
         self%patch(ixI^S)  =.true.
         call self%set_w(ixI^L, ixI^L,qt,x,w_init)


         w_tmp(ixI^S,1:nw) = w_init(ixI^S,1:nw)




         if(ndim>1)then

           select case(trim(self%myconfig%profile_density))
               case('Lee2001','Lee2001_floored')
                  call usr_get_theta(ixI^L,ixO^L,x,theta_profile)

                  select case(typeaxial)
                   case('spherical')
                     where(self%patch(ixO^S))
                       where(d_profile(ixO^S)>self%myconfig%profile_shiftstart)
                        f_profile1(ixO^S,r_) = -w(ixO^S,phys_ind%pressure_)/&
                        d_profile(ixO^S)
                        f_profile1(ixO^S,theta_) = 2.0_dp*w(ixO^S,phys_ind%pressure_)/&
                        (d_profile(ixO^S)*dtan(theta_profile(ixO^S)))
                       end where
                     end where
                   case('cylindrical')
                     where(self%patch(ixO^S))
                       where(d_profile(ixO^S)>self%myconfig%profile_shiftstart)
                         f_profile1(ixO^S,r_) = -4.0_dp*w(ixO^S,phys_ind%pressure_)*&
                         x(ixO^S,r_)/(x(ixO^S,r_)**2.0_dp+x(ixO^S,z_)**2.0_dp)
                         f_profile1(ixO^S,z_) = 2.0_dp*w(ixO^S,phys_ind%pressure_)/&
                         x(ixO^S,z_)
                       end where
                     end where
                   case('slab','slabstretch')
                     where(self%patch(ixO^S))
                       where(d_profile(ixO^S)>self%myconfig%profile_shiftstart)
                         f_profile1(ixO^S,x_) = -4.0_dp*w(ixO^S,phys_ind%pressure_)*&
                         x(ixO^S,x_)/(x(ixO^S,x_)**2.0_dp+x(ixO^S,y_)**2.0_dp)
                         f_profile1(ixO^S,y_) = 2.0_dp*w(ixO^S,phys_ind%pressure_)/&
                         x(ixO^S,y_)
                       end where
                     end where
                   case default
                    call mpistop('Unknown geometry')
                  end select

                   if(trim(self%myconfig%profile_density)=='Lee2001_floored')then
                    select case(typeaxial)
                     case('spherical')
                      where(self%patch(ixO^S))
                       where(theta_profile(ixO^S)<=self%myconfig%theta_floored)
                         f_profile1(ixO^S,theta_) = 2.0_dp*w(ixO^S,phys_ind%pressure_)/&
                         (d_profile(ixO^S)*dtan(self%myconfig%theta_floored))
                       end where
                      end where
                     case('cylindrical')
                      where(self%patch(ixO^S))
                        where(theta_profile(ixO^S)<=self%myconfig%theta_floored)
                          f_profile1(ixO^S,r_) = -4.0_dp*w(ixO^S,phys_ind%pressure_)*&
                          dcos(self%myconfig%theta_floored)/d_profile(ixO^S)
                          f_profile1(ixO^S,z_) = 2.0_dp*w(ixO^S,phys_ind%pressure_)/&
                          (d_profile(ixO^S)*dsin(self%myconfig%theta_floored))
                        end where
                      end where
                     case('slab','slabstretch')
                       where(self%patch(ixO^S))
                         where(theta_profile(ixO^S)<=self%myconfig%theta_floored)
                           f_profile1(ixO^S,x_) = -4.0_dp*w(ixO^S,phys_ind%pressure_)*&
                           dcos(self%myconfig%theta_floored)/d_profile(ixO^S)
                           f_profile1(ixO^S,y_) = 2.0_dp*w(ixO^S,phys_ind%pressure_)/&
                           (d_profile(ixO^S)*dsin(self%myconfig%theta_floored))
                         end where
                        end where
                     case default
                      call mpistop('Unknown geometry')
                    end select
                   end if

               case default
                 Loop_idim_default : do idims=1,ndim
                  where(self%patch(ixO^S))
                    where(d_profile(ixO^S)>self%myconfig%profile_shiftstart)
                     f_profile1(ixO^S,idims) = 0.0_dp
                    end where
                  end where
                 end do Loop_idim_default
           end select

           !if ndim==1:
           !else
           !TO DO
         end if
         ! restore the self%patch to its state before the call of this subroutine
         self%patch(ixI^S)=patch_tosave(ixI^S)

         !2> Compute numerical force separately
         patch_tosave(ixI^S)=self%patch(ixI^S)
         self%patch(ixI^S)  =.true.
         call self%set_w(ixI^L, ixI^L,qt,x,w_init2)
         self%patch(ixI^S)=patch_tosave(ixI^S)
         call phys_to_conserved(ixI^L,ixI^L,w_init2,x)

         w_tmp2(ixI^S,1:nw) = w_init2(ixI^S,1:nw)

         if(phys_config%radiative_cooling)then
           tmp_active=.false.
           src_active=.false.
           call radiative_cooling_add_source(qdt,ixI^L,ixI^L,w_init2,w_tmp2,x,&
                                              src_active,tmp_active)
         end if
         call phys_get_pthermal(w_tmp2,x,ixI^L,ixI^L,ptherm)
         Loop_idir_force : do idims=1,ndim
           call gradient(ptherm,ixI^L,ixO^L,idims,gradp)

           where(self%patch(ixO^S))
             f_profile2(ixO^S,idims) = gradp(ixO^S)
           end where
           pv(ixI^S,idims) =  ptherm(ixI^S)*w(ixI^S,phys_ind%mom(idims))/w(ixI^S,phys_ind%rho_)
         end do  Loop_idir_force
         if(ndir>ndim)pv(ixI^S,ndim+1:ndir)=0.0_dp
         call  divvector(pv,ixI^L,ixO^L,divpv)

         !3> Compute the mean
         select case(trim(self%myconfig%weight_mean_variable))
           case('density_grad')
             maxv = maxval(w_tmp(ixO^S,phys_ind%rho_))
             where(self%patch(ixO^S))
               alpha_analy(ixO^S) = self%myconfig%weight_analy*&
               (w_tmp(ixO^S,phys_ind%rho_)/maxv)
               alpha_num(ixO^S) = self%myconfig%weight_num*&
               ((maxv-w_tmp(ixO^S,phys_ind%rho_))/maxv)
             end where
           case('minus_density_grad')
             maxv = maxval(w_tmp(ixO^S,phys_ind%rho_))
             where(self%patch(ixO^S))
               alpha_analy(ixO^S) = self%myconfig%weight_analy*&
               (maxv/w_tmp(ixO^S,phys_ind%rho_))
               alpha_num(ixO^S) = self%myconfig%weight_num*&
               (maxv/(maxv-w_tmp(ixO^S,phys_ind%rho_)))
             end where
           case default
             where(self%patch(ixO^S))
               alpha_analy(ixO^S) = self%myconfig%weight_analy
               alpha_num(ixO^S) = self%myconfig%weight_num
             end where
         end select
         Loop_idir_force2 : do idims=1,ndim
           select case(trim(self%myconfig%weight_mean))
             case('arithmetic')
                 Where_amplify_force : where(self%patch(ixO^S))
                   f_profile(ixO^S,idims) = ((alpha_analy(ixO^S)**self%myconfig%weight_analy_index)*&
                   f_profile1(ixO^S,idims)+(alpha_num(ixO^S)**self%myconfig%weight_num_index)*&
                   f_profile2(ixO^S,idims))/&
                   (alpha_analy(ixO^S)**self%myconfig%weight_analy_index+&
                   alpha_num(ixO^S)**self%myconfig%weight_num_index)
                 end where Where_amplify_force
             !Not working:
             !case('geometric')
                 !where(self%patch(ixO^S))
                   !fsign(ixO^S) = ((alpha_analy(ixO^S)**self%myconfig%weight_analy_index)*&
                   !f_profile1(ixO^S,idims)+(alpha_num(ixO^S)**self%myconfig%weight_num_index)*&
                   !f_profile2(ixO^S,idims))/&
                   !(alpha_analy(ixO^S)**self%myconfig%weight_analy_index+&
                   !alpha_num(ixO^S)**self%myconfig%weight_num_index)
                   !where(fsign(ixO^S)<smalldouble)
                    !fsign(ixO^S) = 0.0_dp
                   !elsewhere
                    !fsign(ixO^S) = fsign(ixO^S) / dabs(fsign(ixO^S))
                   !end where
                   !f_profile(ixO^S,idims) = fsign(ixO^S)*&
                   !((dabs(f_profile1(ixO^S,idims))**(alpha_analy(ixO^S)**&
                   !self%myconfig%weight_analy_index))*&
                   !(dabs(f_profile2(ixO^S,idims))**(alpha_num(ixO^S)**&
                   !self%myconfig%weight_num_index)))**(1.0_dp/&
                   !((alpha_analy(ixO^S)**self%myconfig%weight_analy_index)+&
                   !(alpha_num(ixO^S)**self%myconfig%weight_num_index)))
                 !end where
             case default
              call mpistop('Unknown weighting mean type')
           end select
         end do  Loop_idir_force2
      end if cond_density_fprofile
    end if cond_fprofile




   end subroutine usr_ism_get_pforce_profile



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
            coef_eos = 1.0_dp!phys_config%gamma/(phys_config%gamma-1.0_dp)

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
