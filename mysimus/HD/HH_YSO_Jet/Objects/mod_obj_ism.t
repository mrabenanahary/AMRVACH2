module mod_obj_ism
  use mod_constants
  use mod_global_parameters
  use mod_physics
  use mod_obj_dust, only : dust
  use mod_obj_global_parameters
  use mod_obj_mat
  use mod_obj_usr_unit
  use mod_grackle_parameters
  use mod_grackle_chemistry
  implicit none

    type ism_parameters
      character(len=20)    :: unit           !> physical unit at parameter file
      character(len=78)    :: obj_name       !> Obj name that call it
      logical              :: normalize_done !> ism is the normalisation is already done
      real(dp)             :: density        !> ISM density  (g/cm^3)
      real(dp)             :: number_density !> ISM number density (1/cm^3)
      real(dp)             :: temperature    !> ISM temperature  (K)
      real(dp)             :: gamma
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
      real(dp)             :: profile_p  !> ISM reference pressire for Lee2001 profile
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
      logical              :: profile_force_with_grav
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
      real(dp)             :: profile_Beta !> ISM envelope fiducial density
      real(dp)             :: profile_rhomax !> ISM envelope fiducial density
      integer              :: patching_idim !> ISM tracer indice
      logical              :: profile_add_pmgrav_to_hse
      logical              :: profile_addvphi
      logical              :: profile_zero_velocities


      logical              :: ism_display_parameters     !> ISM display of parameters at the end of the set_profile subroutine
      logical              :: w_already_initialized
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
      logical              :: useprimitive
      integer              :: nghostcells(3)
      integer              :: escapencells(3)
      integer              :: escapencellsglobal(3)
      logical              :: mixed_fixed_bound(3,2,100) !> ism 'fix' flux variables
      !TO DO    : f95 does not allow to namelist an allocatable array, but fortran 2003+ does
      !So, for now, take a big enough useful array and the same size as in mod_obj_mat.t
      logical              :: dust_on        !> logical to set dust
      real(dp)             :: dust_frac      !> dust fraction

      character(len=30)    :: chemical_gas_type
      real(kind=dp)        :: He_abundance
      real(kind=dp)        :: mean_mass
      real(kind=dp)        :: mean_mup
      real(kind=dp)        :: mean_ne_to_nH
      real(kind=dp)        :: mean_nall_to_nH
    end type ism_parameters


    ! ISM features

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
      integer                  :: idim,iside
      logical                  :: ism_and_usr_boundary

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


      ! original :
      ! if(.not.self%myconfig%boundary_on)                                       &
      !    self%myboundaries%myconfig%boundary_type=self%myconfig%boundary_cond


      if(ism_config%boundary_on)then
        self%myboundaries%myconfig%myindice =1
        call self%myboundaries%read_parameters(self%myboundaries%myconfig,files)
      end if
      if(ism_config%dust_on)then
        self%mydust%myconfig%associated_medium = 'ism'
        call self%mydust%read_parameters(self%mydust%myconfig,files)
      end if

      ! Here the boundary conditions input from ism_config takes priority from what is
      ! input from usrboundary_config
      if(ism_config%boundary_on)then
          do idim=1,ndim
            do iside=1,2
              ism_and_usr_boundary = (trim(ism_config%boundary_cond(idim,iside))/=&
              trim(self%myboundaries%myconfig%boundary_type(idim,iside)))
              if(ism_and_usr_boundary) then
                 self%myboundaries%myconfig%boundary_type(idim,iside)=&
                 ism_config%boundary_cond(idim,iside)
              end if
            end do
          end do
        end if

    end subroutine usr_ism_read_p
   !------------------------------------------------------------------------
   subroutine usr_ism_write_setting(self,unit_config)
     implicit none
     class(ism)                          :: self
     integer,intent(in)                  :: unit_config
     integer                             :: idims2,iside2,iB2
     real(kind=dp)                       :: rto_print
     character(len=64)                   :: sto_print
     character(len=128)                   :: wto_print
     integer                             :: idim,iside,idims,iw2
     ! .. local ..

     !-----------------------------------

     write(unit_config,*)'************************************'
     write(unit_config,*)'************ISM setting ************'
     write(unit_config,*)'************************************'
     write(unit_config,*)'      ****** Code Unit *******      '
     write(unit_config,*) 'Density  at (R,z)=(R_j,0) = ',  self%myconfig%density, '  code unit'
     write(unit_config,*) 'Pressure    = ',  self%myconfig%pressure, '  code unit'
     write(unit_config,*) 'Temperature = ',  self%myconfig%temperature, '  code unit'
     write(unit_config,*) 'Mean molecular weight = ', self%myconfig%mean_mup
     write(unit_config,*) 'Gamma = ',  self%myconfig%gamma
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
     write(unit_config,*) 'Mean molecular weight = ', self%myconfig%mean_mup*&
     self%myphysunit%myconfig%mean_mup
     write(unit_config,*) 'Gamma = ',  self%myconfig%gamma
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
     write(unit_config,*) '======================================================='
     write(unit_config,*)'** Boundary conditions'

     do idims2=1,ndim
      do iside2=1,2
        write(unit_config,*)' *** (idims,iside)=', idims2,iside2

        if(iside2==1) then
           iB2 = 2*(idims2-1)+1
        end if
        if(iside2==2) then
           iB2 = 2*idims2
        end if


        if (any(typeboundary(1:nwfluxbc,iB2)=="special")) then
          if (.not. associated(usr_special_bc)) &
               call mpistop("usr_special_bc not defined")
             write(unit_config,*)' **** BC = ', self%myconfig%boundary_cond(idims2,iside2)
             write(unit_config,*) ' **** ... is this BC fixed ? : '
             if(phys_config%energy)then
               write(unit_config,*) ' **** rho v1 v2 v3 p(or e) tracer1 tracer2'
             else
               write(unit_config,*) ' **** rho v1 v2 v3 tracer1 tracer2'
             end if
             write(unit_config,*) ' **** ',self%myconfig%mixed_fixed_bound(idims2,iside2,1:nwfluxbc),'(ism_config%mixed_fixed_bound)'
             write(unit_config,*) ' **** ',self%myboundaries%myconfig%mixed_fixed_bound(idims2,iside2,1:nwfluxbc),'(ism_config%myboundaries%myconfig%mixed_fixed_bound)'
             write(unit_config,*) ' **** ',self%myboundaries%mixed_fixed_bound(idims2,iside2,1:nwfluxbc),'(ism_config%myboundaries%mixed_fixed_bound)'
             write(unit_config,*) ' **** BC : '
             wto_print=''
             do iw2=1,nwfluxbc
                if(iw2/=1)wto_print = trim(trim(wto_print) // trim(','))
                wto_print = trim(trim(wto_print) // trim(self%myboundaries%variable_typebound(idims2,iside2,iw2)))
             end do
             write(unit_config,*) ' **** ', wto_print
             !write(unit_config,*) ' **** initial BC before eventual fix : '
             !wto_print=trim('')
             !do iw2=1,nwfluxbc
                !if(iw2/=1)wto_print = trim(trim(wto_print) // trim(','))
                !wto_print = trim(trim(wto_print) // trim(typeboundary(iw2,iB2)))
             !end do
             !write(unit_config,*) ' **** ', wto_print
             write(unit_config,*) '--------------------------------'
        else
           wto_print=''
           do iw2=1,nwfluxbc
              if(iw2/=1)wto_print = trim(trim(wto_print) // trim(','))
              wto_print = trim(trim(wto_print) // trim(typeboundary(iw2,iB2)))
           end do
           if(phys_config%energy)then
             write(unit_config,*) ' **** w_BC = rho v1 v2 v3 p(or e) tracer1 tracer2'
           else
             write(unit_config,*) ' **** w_BC = rho v1 v2 v3 tracer1 tracer2'
           end if
           write(unit_config,*) ' **** BC = ', wto_print
           write(unit_config,*) ' **** no need to be fixed '
           write(unit_config,*) '--------------------------------'
        end if
      end do
     end do
     write(unit_config,*) '======================================================='
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
     self%myconfig%gamma                 = 0.0_dp
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
     self%myconfig%profile_p             = 0.0_dp !from Lee 2001 profile
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
     self%myconfig%profile_force_with_grav = .false.

     self%myconfig%ism_display_parameters    = .false.
     self%myconfig%w_already_initialized   = .false.
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
     self%myconfig%profile_Beta  = 0.0_dp
     self%myconfig%patching_idim = 1
     self%myconfig%profile_rhomax= 0.0_dp
     self%myconfig%profile_add_pmgrav_to_hse = .false.
     self%myconfig%profile_addvphi           = .false.
     self%myconfig%profile_zero_velocities   = .false.

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
     self%myconfig%boundary_cond         = 'fix'!'open'
     self%myconfig%flux_frac             = 1.0d-2
     self%myconfig%nghostcells(1:3)           = 2
     self%myconfig%escapencells(1:3)     = 4
     self%myconfig%escapencellsglobal(1:3) = 2
     self%myconfig%useprimitive         = .false.
     self%myconfig%debug                 = .false.
     self%myconfig%mixed_fixed_bound(1:3,1:2,1:nwfluxbc)     = .false.
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
  subroutine usr_ism_set_complet(i_obj,gr_solv,self)
    implicit none
    class(ism)                                :: self
    type(gr_solver), optional :: gr_solv
    integer                   :: i_obj
    ! .. local ..
    logical                  :: dust_is_frac
    real(dp)                 :: Ggrav,r_normalized,costhetazero
    integer                  :: idim,iside,idims,iB2,iw2
    real(dp)                 :: temperature_intermediate, temperature_ratio
    real(dp)                 :: ism_temper, ism_pressu
    real(kind=dp):: nH2, nother
    real(kind=dp):: gammaoth,gammaeff,tvar
    real(kind=dp):: inv_gammaH2_1, inv_gammaoth_1
    real(dp)                 :: mp,kB,me,uenergy, nden
    !-------------------------------------------------------


    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      me = const_me*1.0d-3
    else
      mp=mp_cgs
      kB=kB_cgs
      me = const_me
    end if

    write(*,*) '======================Beginning of ism_set_w settings======================='
    print*, '* at (R,z) = (r_j,0) : '
    print*, '** The ism base temperature is = ',self%myconfig%temperature
    write(*,*) '========================================================'


    ! original :
    ! if(.not.self%myconfig%boundary_on)                                       &
    !    self%myboundaries%myconfig%boundary_type=self%myconfig%boundary_cond


    ! Here the boundary conditions input from ism_config takes priority from what is
    ! input from usrboundary_config





    if(self%myconfig%boundary_on)then
        do idim=1,ndim
          do iside=1,2
            if(trim(self%myconfig%boundary_cond(idim,iside))=='fix') then
               self%myconfig%mixed_fixed_bound(idim,iside,1:nwfluxbc)=.true.
            end if
          end do
        end do
        self%myboundaries%mixed_fixed_bound=self%myconfig%mixed_fixed_bound
        self%myboundaries%myconfig%boundary_type=self%myconfig%boundary_cond
        self%myboundaries%myconfig%flux_frac=self%myconfig%flux_frac
        self%myboundaries%myconfig%useprimitive=self%myconfig%useprimitive
        self%myboundaries%myconfig%nghostcells(1:3)=self%myconfig%nghostcells(1:3)
        self%myboundaries%myconfig%escapencells(1:3)=self%myconfig%escapencells(1:3)
    end if



    Ggrav = constusr%G


    !----------------------------------------------------
    ! set the parameters :
    !   - mean_nall_to_nH = n_e/n_i
    !   - mean_mass = rho_tot/(nH*mH), where nH=n(H+H^+)+2n(H_2)
    !   - mean_mup,mean_ne_to_nH = n(e^-)/nH
    ! from calling the mod_hd_grackle_phys/mod_mhd_phys module
    !----------------------------------------------------

    call phys_fill_chemical_ionisation(self%myconfig%He_abundance,  &
       self%myconfig%chemical_gas_type,                             &
       self%myconfig%mean_nall_to_nH,self%myconfig%mean_mass,       &
       self%myconfig%mean_mup,self%myconfig%mean_ne_to_nH)

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
    !   3) Mach number relat. to moving objet = ||v_jet|| \ c_sound,ism)
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

    if(phys_config%use_grackle)then
      ! gr_solv%myconfig%density(i_obj) = rhoH + rhoHe = rhoXY, then :
      ! self%myconfig%density \simeq rhotot (gaz so no Dust)
      ! = rhoH + rhoHe + rhoZ (ignores deuterium and electron mass densities)
      ! = rhoXY
      ! + self%myconfig%MetalFractionByMass*rhoXY
      ! = rhoXY*(1 + self%myconfig%MetalFractionByMass)
      ! So gr_solv%myconfig%density(i_obj = rhoXY
      ! = self%myconfig%density/(1+ self%myconfig%MetalFractionByMass)
      if(trim(self%myconfig%profile_density)=='Ulrich1976')then
        gr_solv%myconfig%density(i_obj) = self%myconfig%profile_rho_0
      else
        gr_solv%myconfig%density(i_obj) = self%myconfig%density
      end if



      ! here : gr_solv%myconfig%density(i_obj) = rhotot
      call gr_solv%set_complet(gr_solv%myconfig%density(i_obj),&
      i_obj,.false.)
      ! here : gr_solv%myconfig%density(i_obj) = rhoXY
      

write(*,*) 'HI dens 1 =', gr_solv%myconfig%densityHI(i_obj)

      if(self%myconfig%temperature<smalldouble)then
        call mpistop('Need temperature par for Grackle !')
      else
        if(phys_config%primordial_chemistry==0)then
          self%myconfig%mean_mup = 1.0_dp !by default
          self%myconfig%gamma = phys_config%gamma
        else if(phys_config%primordial_chemistry>0)then
          nden = 0.0_dp
          if(phys_config%use_metal_field==1)then
            nden = nden +&
            (gr_solv%myconfig%density_Z(i_obj)/(metal_mu*mp))
          end if
          nden = nden + &
          (gr_solv%myconfig%densityHI(i_obj)/mp)+&
          (gr_solv%myconfig%densityHII(i_obj)/mp)+&
          (gr_solv%myconfig%densityElectrons(i_obj)/me)+&
          ((gr_solv%myconfig%densityHeI(i_obj)+&
          gr_solv%myconfig%densityHeII(i_obj)+&
          gr_solv%myconfig%densityHeIII(i_obj))/(4.0_dp*mp))
          if(phys_config%primordial_chemistry>1)then
            nden = nden + &
            (gr_solv%myconfig%densityHM(i_obj)/mp)+&
            ((gr_solv%myconfig%densityH2I(i_obj)+&
            gr_solv%myconfig%densityH2II(i_obj))/(2.0_dp*mp))           
          end if
          self%myconfig%mean_mup = &
          (gr_solv%myconfig%density(i_obj)/mp)/nden

          if(phys_config%primordial_chemistry==1)then
            self%myconfig%gamma = phys_config%gamma
          else if(phys_config%primordial_chemistry>1)then
            ! from given temperature
            nother = ((gr_solv%myconfig%densityHeI(i_obj)+&
            gr_solv%myconfig%densityHeII(i_obj)+&
            gr_solv%myconfig%densityHeIII(i_obj))/(4.0d0*mp)) +&
            (gr_solv%myconfig%densityHI(i_obj)/mp) + &
            (gr_solv%myconfig%densityHII(i_obj)/mp) +&
            (gr_solv%myconfig%densityHM(i_obj)/mp) + &
            (gr_solv%myconfig%densityElectrons(i_obj)/me)
            
            nH2 = (gr_solv%myconfig%densityH2I(i_obj)+&
            gr_solv%myconfig%densityH2II(i_obj))/(2.0_dp*mp)
              
            inv_gammaH2_1 = 2.5_dp
            if(nH2/nother > 1.0d-3)then
              tvar = 6100.0d0/self%myconfig%temperature
            if (tvar < 10.0)then
          
              inv_gammaH2_1 = 0.5d0*(5.0d0 +&
              2.0d0 * tvar**2.0d0 *&
              DEXP(tvar)/(DEXP(tvar)-1.0d0)**2.0d0)
        
            end if
            end if

            gammaoth = phys_config%gamma
            inv_gammaoth_1 = 1.0d0/(gammaoth-1.0d0)
            gammaeff = 1.0d0 + ((nH2 + nother)/&
            (nH2*inv_gammaH2_1 +&
            nother*inv_gammaoth_1))
write(*,*) 'ISM set_complet nH2 =', nH2
write(*,*) 'ISM set_complet nother =', nother
write(*,*) 'ISM set_complet inv_gammaH2_1 = ', inv_gammaH2_1
write(*,*) 'ISM set_complet gammaoth = ', gammaoth
write(*,*) 'ISM set_complet inv_gammaoth_1 = ', inv_gammaoth_1
write(*,*) 'ISM set_complet gammaeff = ', gammaeff

            self%myconfig%gamma = gammaeff
          end if
        end if
      end if

      ! gr_solv%myconfig%density(i_obj) 
      ! = density for Grackle = rhoH + rhoHe
      ! Eqs. (14) & (27) of Smith et al. 2017, MNRAS 466, 2217–2234      
      ism_pressu = gr_solv%myconfig%density(i_obj)*&
      kB*self%myconfig%temperature/&
      (self%myconfig%mean_mup*mp)

      ism_temper = self%myconfig%temperature

      if(trim(self%myconfig%profile_density)=='Ulrich1976')then
        self%myconfig%pressure_Ulrich = ism_pressu
      else
        self%myconfig%pressure = ism_pressu
      end if
      write(*,*) 'ISM set_complet rhoXY',gr_solv%myconfig%density(i_obj)
      write(*,*) 'ISM set_complet rhogas',gr_solv%myconfig%density_gas(i_obj)
      write(*,*) 'ISM set_complet density',self%myconfig%density
      write(*,*) 'ISM set_complet rhotot',gr_solv%myconfig%density_tot(i_obj)
      write(*,*) 'ISM set_complet HIdensity',gr_solv%myconfig%densityHI(i_obj)
      write(*,*) 'ISM set_complet H2Idensity',gr_solv%myconfig%densityH2I(i_obj)
      write(*,*) 'ISM set_complet HeIdensity',gr_solv%myconfig%densityHeI(i_obj)
      write(*,*) 'ISM set_complet e_density',gr_solv%myconfig%densityElectrons(i_obj)
      write(*,*) 'ISM set_complet temperature',self%myconfig%temperature
      write(*,*) 'ISM set_complet pressure',self%myconfig%pressure
      write(*,*) 'ISM set_complet gamma',self%myconfig%gamma
      write(*,*) 'ISM set_complet mean_mup',self%myconfig%mean_mup


      self%myconfig%c_sound =&
      sqrt(self%myconfig%gamma*ism_pressu/&
      gr_solv%myconfig%density(i_obj))

    else

      !----------------------------------------------------
      ! According to what is input and treated previously,
      ! set the ism pressure from c_sound
      ! 19-03-21 : added the missing setting of temperature accordingly
      ! Priority is for c_sound
      !----------------------------------------------------

      self%myconfig%gamma = phys_config%gamma
      cond_csound_set : if(self%myconfig%c_sound>0.0_dp) then
         self%myconfig%pressure = self%myconfig%c_sound**2.0_dp * self%myconfig%density /&
                                  self%myconfig%gamma
         if(trim(self%myconfig%profile_density)=='Ulrich1976')then
           self%myconfig%pressure_Ulrich = self%myconfig%c_sound**2.0_dp *&
                          self%myconfig%profile_rho_0 / self%myconfig%gamma
         end if
         self%myconfig%temperature = self%myconfig%c_sound**2.0_dp*mp*&
                                  self%myconfig%mean_mup/(kB*self%myconfig%gamma)
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
        !fully atomic gas : self%myconfig%mean_mup ~ 1.27
        self%myconfig%pressure = self%myconfig%density*&
                              kB* self%myconfig%temperature/(self%myconfig%mean_mup*mp)
        if(trim(self%myconfig%profile_density)=='Ulrich1976')then
          self%myconfig%pressure_Ulrich = self%myconfig%profile_rho_0*&
                                kB* self%myconfig%temperature/(self%myconfig%mean_mup*mp)
        end if
       else if(dabs(self%myconfig%pressure)>0.0_dp) then
        self%myconfig%temperature = self%myconfig%mean_mup*mp*self%myconfig%pressure/&
                              (kB*self%myconfig%density)
       end if
       if(trim(self%myconfig%profile_density)=='Ulrich1976')then
         self%myconfig%c_sound = sqrt(self%myconfig%gamma*self%myconfig%pressure_Ulrich/&
              self%myconfig%profile_rho_0)
       else
          self%myconfig%c_sound = sqrt(self%myconfig%gamma*self%myconfig%pressure/self%myconfig%density)
       end if
      end if cond_csound_set
    end if

    !write(*,*) 'pressure in ISM parameters set_complet == ', self%myconfig%pressure

    !set here pressure and temperature parameters according to grackle

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


   self%myconfig%profile_p = 2.0_dp*kB*self%myconfig%temperature/&
   (self%myconfig%mean_mup*mp)

   self%myconfig%profile_Beta = kB*self%myconfig%temperature/&
   (self%myconfig%mean_mup*mp*self%myconfig%profile_rd)

   if(.not.self%myconfig%profile_add_pmgrav_to_hse)then
    self%myconfig%profile_zero_velocities = .true.
   end if

    call self%myboundaries%set_complet



    if(mype==0)then
      write(*,*) '======================ISM settings======================='
      print*, '* at (R,z) = (r_j,0) : '
      print*, '** The ism base temperature is = ',self%myconfig%temperature
      print*, '** The ism base density is =', self%myconfig%density
      print*, '** The ism base meanmum is =', self%myconfig%mean_mup
      print*, '** The ism base mp is =', mp
      print*, '** The ism base kB is =', kB
      print*, '** The ism base pressure is =', self%myconfig%pressure
      print*, '** The ism base number density nH is ', self%myconfig%number_density
      write(*,*) '========================================================='
    end if



write(*,*) 'HI dens 2 =', gr_solv%myconfig%densityHI(i_obj)
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
     self%myconfig%profile_p            = self%myconfig%profile_p/(physunit_inuse%myconfig%pressure/physunit_inuse%myconfig%density)
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
     self%myconfig%profile_Beta  = self%myconfig%profile_Beta  / (physunit_inuse%myconfig%pressure/(physunit_inuse%myconfig%length*&
                                                                  physunit_inuse%myconfig%density))
     self%myconfig%profile_rhomax=self%myconfig%profile_rhomax / physunit_inuse%myconfig%density

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
   subroutine usr_ism_set_w(ixI^L,ixO^L,qt,x,w,isboundary_iB,gr_solv,self)
   ! * According to subroutine initonegrid_usr in which it is called once in mod_obj_usr_yso_jet.t :
   ! in the src code amrini.t, ixI^L = ixG^LL = ixGlo1,ixGlo2,ixGhi1,ixGhi2
   ! = range delimiting the whole domain (including boundary ghost cells), i.e. between integers
   ! 1 and the highest possible indices for the coordinates for the grid for each dimension
   ! and ixO^L = ixM^LL = ixMlo1,ixMlo2,ixMhi1,ixMhi2 =
   ! = range delimiting the mesh (i.e. the grid without boundary layers)
   ! * According to subroutine bc_phys in which it is also called once from
   ! usr_special_bc ==> specialbound_usr by mod_ghostcells_update.t twice:
   ! in this case, the integers are
   ! ixG^L = ixGmin1,ixGmin2,ixGmax1,ixGmax2 (or ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2
   ! when working with a coarsened grid)
   ! = the range delimiting the whole domain (including boundary ghost cells), i.e. between integers
   ! 1 and the highest possible indices for the coordinates for the grid for each dimension
   ! ,ixB^L = ixBmin1,ixBmin2,ixBmax1,ixBmax2 =
   ! the range delimiting the boundary
   ! As for the input integers, they are ixI^L=ixG^L
   ! and ixO^L=ixI^L defined from ixB^L as :
   ! > idims = 1, iside==2 (maximal boundary)
   ! ==> ixImin1=ixBmax1+1-nghostcells;ixImin2=ixBmin2;ixImax1=ixBmax1;ixImax2=ixBmax2;
   ! > idims = 1, iside==1 (minimal boundary)
   ! ==> ixImin1=ixBmin1;ixImin2=ixBmin2;ixImax1=ixBmin1-1+nghostcells;ixImax2=ixBmax2;
   ! > idims = 2, iside==2 (maximal boundary)
   ! ==> ixImin1=ixBmin1;ixImin2=ixBmax2+1-nghostcells;ixImax1=ixBmax1;ixImax2=ixBmax2;
   ! > idims = 2, iside==1 (minimal boundary)
   ! ==> ixImin1=ixBmin1;ixImin2=ixBmin2;ixImax1=ixBmax1;ixImax2=ixBmin2-1+nghostcells;
   ! to sum up, the range delimiting each boundary individually


    implicit none
    integer, intent(in)           :: ixI^L,ixO^L
    real(kind=dp), intent(in)     :: qt
    real(kind=dp), intent(in)     :: x(ixI^S,1:ndim)
    real(kind=dp)                 :: w(ixI^S,1:nw)
    real(kind=dp)                 :: w_copy(ixI^S,1:nw)
    real(kind=dp),allocatable     :: w_tmp(:^D&,:)
    real(kind=dp)                 :: grid_local(1:ndim)
    class(ism)                    :: self
    integer,             optional :: isboundary_iB(2)
    type(gr_solver), optional :: gr_solv

    ! .. local..
    integer                    :: idir,iB,idims,idims_bound,iw,iwfluxbc,idims2,iside2
    integer                    :: iobject
    real(kind=dp)              :: fprofile(ixI^S), mmw_field(ixI^S)

    logical                    :: isboundary,to_fix,bc_to_fix,some_unfixed
    character(len=30)          :: myboundary_cond
    character(len=64)          :: aString
    !----------------------------------

    !write(*,*) 'HI dens 3 =', gr_solv%myconfig%densityHI(1)
    if(.not.self%myconfig%w_already_initialized)self%myconfig%w_already_initialized=.true.


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


      !boundary_cond : if(.not.isboundary.or.&
      !                      trim(myboundary_cond)=='fix') then

           if(trim(self%myconfig%profile_density)=='Ulrich1976')then
             where(self%patch(ixO^S))
                 w(ixO^S,phys_ind%rho_)        =  self%myconfig%profile_rho_0
             end where
           else
             where(self%patch(ixO^S))
                 !rhotot
                 w(ixO^S,phys_ind%rho_)        =  self%myconfig%density
             end where
           end if


          iobject = self%myconfig%myindice + 1
          !write(*,*) 'HI dens = ', gr_solv%myconfig%densityHI(iobject)
          !Grackle uniform fields initialization : 
          if(phys_config%energy)then
            if(phys_config%use_grackle)then
              if(phys_config%primordial_chemistry>0)then
              where(self%patch(ixO^S))
                w(ixO^S,phys_ind%HI_density_)        =  gr_solv%myconfig%densityHI(iobject)
                w(ixO^S,phys_ind%HII_density_)        =  gr_solv%myconfig%densityHII(iobject)
                w(ixO^S,phys_ind%HeI_density_)        =  gr_solv%myconfig%densityHeI(iobject)
                w(ixO^S,phys_ind%HeII_density_)        =  gr_solv%myconfig%densityHeII(iobject)
                w(ixO^S,phys_ind%HeIII_density_)        =  gr_solv%myconfig%densityHeIII(iobject)
                w(ixO^S,phys_ind%e_density_)        =  gr_solv%myconfig%densityElectrons(iobject)
              end where
              end if
              if(phys_config%primordial_chemistry>1)then
              where(self%patch(ixO^S))
                w(ixO^S,phys_ind%HM_density_)        =  gr_solv%myconfig%densityHM(iobject)
                w(ixO^S,phys_ind%H2I_density_)        =  gr_solv%myconfig%densityH2I(iobject)
                w(ixO^S,phys_ind%H2II_density_)        =  gr_solv%myconfig%densityH2II(iobject)
              end where
              end if
              if(phys_config%primordial_chemistry>2)then
              where(self%patch(ixO^S))
                w(ixO^S,phys_ind%DI_density_)  =  gr_solv%myconfig%densityDI(iobject)
                w(ixO^S,phys_ind%DII_density_) =  gr_solv%myconfig%densityDII(iobject)
                w(ixO^S,phys_ind%HDI_density_) =  gr_solv%myconfig%densityHDI(iobject)                
              end where
              end if
              if(phys_config%use_metal_field==1)then
              where(self%patch(ixO^S))
                w(ixO^S,phys_ind%metal_density_) =  gr_solv%myconfig%density_Z(iobject)
              end where
              end if
              if(phys_config%use_dust_density_field==1)then
              where(self%patch(ixO^S))
                w(ixO^S,phys_ind%dust_density_) =  gr_solv%myconfig%density_dust(iobject)
              end where
              end if
            end if
          end if


          !write(*,*) 'HI dens 4= ', w(ixO^S,phys_ind%HI_density_)

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
           !write(*,*) 'pressure in ISM parameters set_w == ', self%myconfig%pressure*&
           !self%myphysunit%myconfig%pressure
           if(phys_config%energy)then
              !write(*,*) 'ISM set_w temperature',self%myconfig%temperature*&
              !self%myphysunit%myconfig%temperature
             where(self%patch(ixO^S))
               w(ixO^S,phys_ind%temperature_) = self%myconfig%temperature
             end where
           end if

          if(phys_config%mean_mup_on) then
            where(self%patch(ixO^S))
                w(ixO^S,phys_ind%mup_) = self%myconfig%mean_mup
            end where
          end if    

where(self%patch(ixO^S))
w(ixO^S,phys_ind%gamma_) = self%myconfig%gamma
end where  



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
                call self%set_profile(ixI^L,ixO^L,x,w,isboundary)
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








         !write(*,*) 'pressure field in ISM = ', w(ixO^S,phys_ind%pressure_)



         !write(*,*) 'HI dens 4= ', w(ixO^S,phys_ind%HI_density_)


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
   subroutine usr_ism_set_profile(ixI^L, ixO^L,x,w,self,isboundary)
    implicit none
    integer, intent(in)             :: ixI^L,ixO^L
    logical, intent(in),optional     :: isboundary
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

          call self%set_ulrich_profile(ixI^L, ixO^L,x,w,p_profile,isboundary)

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
    end if
    if(self%myconfig%profile_density_on)then
      where(self%patch(ixO^S))
        where(d_profile(ixO^S)>self%myconfig%profile_shiftstart)
          w(ixO^S,phys_ind%rho_)      = w(ixO^S,phys_ind%rho_)* p_profile(ixO^S)
        elsewhere
          w(ixO^S,phys_ind%rho_)      = &
             self%myconfig%profile_shift_density(1,1)

        end where
      end where
    end if

    !Grackle stratified fields : 
    if(phys_config%energy)then
      if(phys_config%use_grackle)then
        if(phys_config%primordial_chemistry>0)then
        where(self%patch(ixO^S))
        w(ixO^S,phys_ind%HI_density_)=w(ixO^S,phys_ind%HI_density_)*p_profile(ixO^S)
        w(ixO^S,phys_ind%HII_density_)=w(ixO^S,phys_ind%HII_density_)*p_profile(ixO^S)
        w(ixO^S,phys_ind%HeI_density_)=w(ixO^S,phys_ind%HeI_density_)*p_profile(ixO^S)
        w(ixO^S,phys_ind%HeII_density_)=w(ixO^S,phys_ind%HeII_density_)*p_profile(ixO^S)
        w(ixO^S,phys_ind%HeIII_density_)=w(ixO^S,phys_ind%HeIII_density_)*p_profile(ixO^S)
        w(ixO^S,phys_ind%e_density_)=w(ixO^S,phys_ind%e_density_)*p_profile(ixO^S)
        end where
        end if
        if(phys_config%primordial_chemistry>1)then
        where(self%patch(ixO^S))
        w(ixO^S,phys_ind%HM_density_)=w(ixO^S,phys_ind%HM_density_)*p_profile(ixO^S)
        w(ixO^S,phys_ind%H2I_density_)=w(ixO^S,phys_ind%H2I_density_)*p_profile(ixO^S)
        w(ixO^S,phys_ind%H2II_density_)=w(ixO^S,phys_ind%H2II_density_)*p_profile(ixO^S)                
        end where
        end if
        if(phys_config%primordial_chemistry>2)then
        where(self%patch(ixO^S))
        w(ixO^S,phys_ind%DI_density_)=w(ixO^S,phys_ind%DI_density_)*p_profile(ixO^S)
        w(ixO^S,phys_ind%DII_density_)=w(ixO^S,phys_ind%DII_density_)*p_profile(ixO^S)
        w(ixO^S,phys_ind%HDI_density_)=w(ixO^S,phys_ind%HDI_density_)*p_profile(ixO^S)                
        end where
        end if
        if(phys_config%use_metal_field==1)then
        where(self%patch(ixO^S))
        w(ixO^S,phys_ind%metal_density_)=w(ixO^S,phys_ind%metal_density_)*p_profile(ixO^S)
        end where
        end if
        if(phys_config%use_dust_density_field==1)then
        where(self%patch(ixO^S))
        w(ixO^S,phys_ind%dust_density_)=w(ixO^S,phys_ind%dust_density_)*p_profile(ixO^S)
        end where
        end if
      end if
    end if



   end subroutine usr_ism_set_profile

   !--------------------------------------------------------------------
   subroutine usr_set_ulrich_profile(ixI^L, ixO^L,x,w,p_profile,self,isboundary)
   implicit none
   integer, intent(in)             :: ixI^L,ixO^L
   logical, intent(in),optional    :: isboundary
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
   real(kind=dp)                   :: fce_factor,fce_vphi
   !--------------------------------------------------------------------


   {zero_dim(^D)=0.0_dp\}!-dx(1,1)

   fce_factor = 1.0_dp
   fce_vphi = 1.0_dp

   if(self%myconfig%profile_zero_velocities)then
    fce_factor = 0.0_dp
   end if

   if(.not.self%myconfig%profile_addvphi)then
    fce_vphi = 0.0_dp
   end if

   p_profile(ixO^S) = 1.0_dp
   vr_profile(ixO^S) = 1.0_dp
   vt_profile(ixO^S) = 1.0_dp
   vp_profile(ixO^S) = 1.0_dp

   call self%set_profile_distance(ixI^L,ixO^L,x,d_profile)

   call usr_get_theta(ixI^L,ixO^L,x,theta_profile,self%myboundaries%myconfig%special_origin_theta)

   !w(ixI^S,phys_ind%mythetafield_)=theta_profile(ixI^S)

   r_normalized(ixI^S) = (d_profile(ixI^S)/self%myconfig%profile_rd)

   call usr_ulrich1976_costheta_zero(ixI^L, ixO^L,x,r_normalized,&
   theta_profile,cos_theta_zero)


   ! considering the semi-planes z>=0 and z<0:
   !where(x(ixI^S,z_)>0.0_dp)
    theta_zero(ixI^S) = DACOS(cos_theta_zero(ixI^S))
    !where(x(ixI^S,z_)<0.0_dp)
      !theta_zero(ixI^S) = dpi-theta_zero(ixI^S)
    !end where
   !elsewhere
    !theta_zero(ixI^S) = DACOS(-cos_theta_zero(ixI^S))
   !end where

   ! Do not forget to recompute the cos_theta_zero after correcting theta_zero,
   ! thus considering the semi-planes z>=0 and z<0::
   cos_theta_zero(ixI^S) = DCOS(theta_zero(ixI^S))

   !w(ixI^S,phys_ind%mytheta_zero_)=theta_zero(ixI^S)

   !if(isboundary) then
     !do iw=1,nw
     !print*,'>>>>>>'
     !print*,'theta = ',theta_profile(ixO^S)
     !print*,'cos(theta) = ',DCOS(theta_profile(ixO^S))
     !print*,'sin(theta) = ',DSIN(theta_profile(ixO^S))
     !print*,'theta_zero = ',theta_zero(ixO^S)
     !print*,'cos(theta_zero) = ',cos_theta_zero(ixO^S)
     !print*,'sin(theta_zero) = ',DSIN(theta_zero(ixO^S))
     !end do
   !end if

   p_profile(ixO^S) = (r_normalized(ixO^S)**(-1.5_dp))*&
   ((1.0_dp + (DCOS(theta_profile(ixO^S))/&
   cos_theta_zero(ixO^S)))**(-0.5_dp))/&
   ((1.0_dp + (3.0_dp * cos_theta_zero(ixO^S) * &
   cos_theta_zero(ixO^S) - 1.0_dp)/r_normalized(ixO^S)))

   vr_profile(ixO^S) = - fce_factor * (r_normalized(ixO^S)**(-0.5_dp))*&
   ((1.0_dp + (DCOS(theta_profile(ixO^S))/cos_theta_zero(ixO^S)))**(0.5_dp))

   vt_profile(ixO^S) =  fce_factor * (r_normalized(ixO^S)**(-0.5_dp))*&
   ((1.0_dp + (DCOS(theta_profile(ixO^S))/cos_theta_zero(ixO^S)))**(0.5_dp))*&
   ((cos_theta_zero(ixO^S)-DCOS(theta_profile(ixO^S)))/DSIN(theta_profile(ixO^S)))

   vp_profile(ixO^S) =  fce_vphi * (r_normalized(ixO^S)**(-0.5_dp))*&
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
            where(x(ixO^S,z_)>=0.0_dp.and.x(ixO^S,r_)>=0.0_dp)
              project_speed(ixO^S,r_) = (vr_profile(ixO^S)*&
              DSIN(theta_profile(ixO^S)))+(vt_profile(ixO^S)*&
              DCOS(theta_profile(ixO^S)))
              project_speed(ixO^S,z_) = (vr_profile(ixO^S)*&
              DCOS(theta_profile(ixO^S)))-(vt_profile(ixO^S)*&
              DSIN(theta_profile(ixO^S)))
            elsewhere(x(ixO^S,z_)<0.0_dp.and.x(ixO^S,r_)>=0.0_dp)
              !v_R(z<0)=v_R(z<0)
              project_speed(ixO^S,r_) = (vr_profile(ixO^S)*&
              DSIN(theta_profile(ixO^S)))+(vt_profile(ixO^S)*&
              DCOS(theta_profile(ixO^S)))
              !v_Z(z<0)=-v_Z(z>0)
              project_speed(ixO^S,z_) = -((vr_profile(ixO^S)*&
              DCOS(theta_profile(ixO^S)))-(vt_profile(ixO^S)*&
              DSIN(theta_profile(ixO^S))))
            elsewhere(x(ixO^S,z_)<0.0_dp.and.x(ixO^S,r_)<0.0_dp)
              !v_R(R<0)=v_R(R>0)
              project_speed(ixO^S,r_) = -((vr_profile(ixO^S)*&
              DSIN(theta_profile(ixO^S)))+(vt_profile(ixO^S)*&
              DCOS(theta_profile(ixO^S))))
              !v_Z(R<0)=-v_Z(R>0)
              project_speed(ixO^S,z_) = -((vr_profile(ixO^S)*&
              DCOS(theta_profile(ixO^S)))-(vt_profile(ixO^S)*&
              DSIN(theta_profile(ixO^S))))
            elsewhere(x(ixO^S,z_)>=0.0_dp.and.x(ixO^S,r_)<0.0_dp)
              !v_R(R<0)=v_R(R>0)
              project_speed(ixO^S,r_) = -((vr_profile(ixO^S)*&
              DSIN(theta_profile(ixO^S)))+(vt_profile(ixO^S)*&
              DCOS(theta_profile(ixO^S))))
              !v_Z(R<0)=v_Z(R>0)
              project_speed(ixO^S,z_) = (vr_profile(ixO^S)*&
              DCOS(theta_profile(ixO^S)))-(vt_profile(ixO^S)*&
              DSIN(theta_profile(ixO^S)))
            end where
          end where
        else
          where(self%patch(ixO^S))
            project_speed(ixO^S,r_) = 1.0_dp ! TO DO
            project_speed(ixO^S,z_) = 1.0_dp ! TO DO
          end where
        end if
        if(ndir>2)then
           where(self%patch(ixO^S))
             where(x(ixO^S,r_)>=0.0_dp)
              project_speed(ixO^S,phi_) = vp_profile(ixO^S)
             elsewhere
              project_speed(ixO^S,phi_) = -vp_profile(ixO^S)
             end where
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


       !if(isboundary) then
         !do iw=1,nw
         !print*,'v_phi = ',w(ixO^S,phys_ind%mom(phi_))*self%myphysunit%myconfig%velocity
         !print*,'>>>>'
         !end do
       !end if


       end  subroutine usr_set_ulrich_profile


   !--------------------------------------------------------------------

   subroutine usr_ism_get_pforce_profile(ixI^L,ixO^L,qt,qdt,x,w,f_profile,divpv,gr_solv,grid_dx,self)
    use mod_radiative_cooling
    implicit none
    integer, intent(in)                  :: ixI^L,ixO^L
    real(kind=dp), intent(in)            :: qt,qdt
    real(kind=dp)                        :: maxv,minv
    real(kind=dp), intent(in)            :: x(ixI^S,1:ndim)
    real(kind=dp), intent(in)            :: w(ixI^S,1:nw)
    type(gr_solver)                         :: gr_solv
    real(kind=dp), intent(in)            :: grid_dx(1:ndim)
    class(ism)                           :: self
    real(kind=dp), intent(inout)         :: f_profile(ixI^S,1:ndim)
    real(kind=dp), intent(inout)         :: divpv(ixI^S)
    ! ..local ..
    integer                              :: idims
    real(kind=dp), dimension(ixI^S,1:ndir) :: pv
    real(kind=dp), dimension(ixI^S,1:nw) :: w_init, w_tmp
    real(kind=dp), dimension(ixI^S,1:nw) :: w_init2, w_tmp2
    real(kind=dp), dimension(ixI^S)      :: ptherm,gradp,d_profile,theta_profile
    real(kind=dp), dimension(ixI^S)      :: alpha_analy,alpha_num,fsign
    logical                              :: tmp_active,src_active
    logical, dimension(ixI^S)            :: patch_tosave,all_patch
    real(dp)                 :: mp,kB,Ggrav,alpha
    !----------------------------------------------------

    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
    else
      mp=mp_cgs
      kB=kB_cgs
    end if

    Ggrav = constusr%G

    f_profile(ixO^S,:) = 0.0_dp
    call self%set_profile_distance(ixI^L,ixO^L,x,d_profile)

    patch_tosave(ixI^S)=self%patch(ixI^S)
    self%patch(ixI^S)  =.true.
    call self%set_w(ixI^L, ixI^L,qt,x,w_init,gr_solv=gr_solv)
    self%patch(ixI^S)=patch_tosave(ixI^S)

    call phys_to_conserved(ixI^L,ixI^L,w_init,x)

    w_tmp(ixI^S,1:nw) = w_init(ixI^S,1:nw)

    cond_density_fprofile : if(self%myconfig%profile_density_on) then
     cond_fprofile_numerical : if(self%myconfig%profile_force_gradP_on)then

      if(phys_config%radiative_cooling)then
        tmp_active=.false.
        src_active=.false.
        call radiative_cooling_add_source(qdt,ixI^L,ixI^L,w_init,w_tmp,x,&
                                           src_active,tmp_active)
      end if
      !write(*,*) 'Before gr_solv%grackle_source w(ixO^S,e_)',w(ixO^S,phys_ind%e_)
      if(phys_config%use_grackle)then
        all_patch(ixI^S) = .true.
        !call MPI_BARRIER(icomm, ierrmpi)
        !write(*,*) ' iteration number = ', it
        !call gr_solv%grackle_source(ixI^L,ixO^L,x,qdt,0.0d0,w_init,&
        !qt,w_tmp,grid_dx,&
        !all_patch)
        !call MPI_BARRIER(icomm, ierrmpi)
        if(self%myconfig%myindice+1/=1) then
          call mpistop('beware : self%myconfig%myindice = ',&
          self%myconfig%myindice+1,' /= 1')
        end if
      end if
      !write(*,*) 'After gr_solv%grackle_source w(ixO^S,e_)',w(ixO^S,phys_ind%e_)
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

      call phys_get_pthermal(w_tmp,x,ixI^L,ixI^L,ptherm)

      if(ndim>1)then

        select case(trim(self%myconfig%profile_density))
            case('Lee2001','Lee2001_floored')
             !call usr_get_theta(ixI^L,ixO^L,x,theta_profile)

                  alpha = self%myconfig%profile_p
                  where(self%patch(ixO^S))

                      f_profile(ixO^S,r_) = - (alpha * (x(ixO^S,r_)**2.0_dp-&
                      x(ixO^S,z_)**2.0_dp)/(x(ixO^S,r_)*(x(ixO^S,r_)**2.0_dp+&
                      x(ixO^S,z_)**2.0_dp))) * w(ixO^S,phys_ind%rho_)

                      f_profile(ixO^S,z_) = - (2.0_dp * alpha * x(ixO^S,z_) / &
                      (x(ixO^S,r_)**2.0_dp+x(ixO^S,z_)**2.0_dp)) * w(ixO^S,phys_ind%rho_)
                      !call mpistop('the code arrives here !!! ')

                  end where

            case default
              Loop_idim_default : do idims=1,ndim
               where(self%patch(ixO^S))
                 where(d_profile(ixO^S)>self%myconfig%profile_shiftstart)
                  f_profile(ixO^S,idims) = 0.0_dp
                 end where
               end where
              end do Loop_idim_default
        end select

      !if(ndim==1) then
      !else
      !TO DO
      end if
      ! restore the self%patch to its state before the call of this subroutine
      ! self%patch(ixI^S)=patch_tosave(ixI^S)

      divpv(ixI^S) = 0.0_dp
      ! divpv(ixI^S) = rho * (Vec(grav) . Vec(v))
      {^D&
      divpv(ixI^S) =  divpv(ixI^S) + f_profile(ixO^S,^D) * w(ixI^S,phys_ind%mom(^D))
      }

      call phys_get_pthermal(w_tmp,x,ixI^L,ixI^L,ptherm)
      do idims=1,ndim
        call gradient(ptherm,ixI^L,ixO^L,idims,gradp)

        where(self%patch(ixO^S))
          f_profile(ixO^S,idims) = gradp(ixO^S)
        end where
        pv(ixI^S,idims) =  ptherm(ixI^S)*w(ixI^S,phys_ind%mom(idims))/w(ixI^S,phys_ind%rho_)
      end do
      if(ndir>ndim)pv(ixI^S,ndim+1:ndir)=0.0_dp
      call  divvector(pv,ixI^L,ixO^L,divpv)






    end if cond_fprofile_numerical

      if(self%myconfig%profile_force_with_grav)then
        where((DABS(x(ixO^S,r_))>self%myconfig%escapencells(r_)*&
              dx(r_,1)).and.self%patch(ixO^S))
              ^D&f_profile(ixO^S,^D)=0.0_dp\
        end where
      end if

    end if cond_density_fprofile

   end subroutine usr_ism_get_pforce_profile



   !--------------------------------------------------------------------
   subroutine usr_ism_add_source(ixI^L,ixO^L,iw^LIM,x,qdt,qtC,wCT,qt,w,gr_solv,grid_dx,self,&
                                 use_tracer,escape_patch,source_filter)
     implicit none
     integer, intent(in)                     :: ixI^L,ixO^L,iw^LIM
     real(kind=dp), intent(in)               :: qdt,qtC,qt
     real(kind=dp), intent(in)               :: x(ixI^S,1:ndim)
     real(kind=dp), intent(in)               :: wCT(ixI^S,1:nw)
     real(kind=dp), intent(inout)            :: w(ixI^S,1:nw)
     type(gr_solver)                         :: gr_solv
     real(kind=dp), intent(in)            :: grid_dx(1:ndim)
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
     real(kind=dp),dimension(ixI^S) :: mmw_field,gmmeff_field,pth_field
     !---------------------------------------------------------
     ! This subroutine implicitly uses conservative variables


      call self%alloc_set_patch(ixI^L,ixI^L,qt,x,&
                   use_tracer=use_tracer,w=w,escape_patch=escape_patch)

      cond_add_force : if(self%myconfig%profile_force_on) then

        cond_inside_prof: if(any(self%patch(ixO^S)))then
          call self%get_pforce_profile(ixI^L,ixO^L,qtC,qdt,x,wCT,f_profile,divpv,&
          gr_solv,grid_dx)
          Loop_idim_force_m : do idims = 1,ndim
          where(self%patch(ixO^S))
           ! This subroutine implicitly uses conservative variables
            w(ixO^S,phys_ind%mom(idims)) = w(ixO^S,phys_ind%mom(idims))&
                +qdt*f_profile(ixO^S,idims)

          end where
        end do Loop_idim_force_m
          if(phys_config%energy) then
            coef_eos = 1.0_dp!self%myconfig%gamma/(self%myconfig%gamma-1.0_dp)

            !write(*,*) 'Before add_source w(ixO^S,e_)',w(ixO^S,phys_ind%e_)
            where(self%patch(ixO^S))
            ! This subroutine implicitly uses conservative variables
             w(ixO^S,phys_ind%e_)=w(ixO^S,phys_ind%e_)  &
                + qdt*coef_eos*divpv(ixO^S)
            end where
            !write(*,*) 'After add_source w(ixO^S,e_)',w(ixO^S,phys_ind%e_)

          end if
        end if  cond_inside_prof
      end if cond_add_force

      ! Compute field of H2-corrected gamma
      call phys_get_gamma(w,ixI^L, ixO^L, gmmeff_field)
      ! Compute field of mean molecular weight
      call phys_get_mup(w, x, ixI^L, ixO^L, mmw_field)
      ! Compute field of H2-gamma corrected temperature
      call phys_get_temperature(ixI^L, ixO^L,w, x, temperature)
      ! Compute field of H2-gamma corrected pressure
      !call phys_get_pthermal(w, x, ixI^L, ixO^L, pth_field)


      if(self%myconfig%profile_force_on) then

        if(any(self%patch(ixO^S)))then
          if(phys_config%energy) then
            if(phys_config%use_grackle)then
            where(self%patch(ixO^S))
            ! This subroutine implicitly uses conservative variables
             w(ixO^S,phys_ind%gamma_)=gmmeff_field(ixO^S)
             w(ixO^S,phys_ind%mup_)=mmw_field(ixO^S)/w_convert_factor(phys_ind%mup_)
             w(ixO^S,phys_ind%temperature_)=temperature(ixO^S)
            end where
            end if
          end if
        end if
      end if


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
