module mod_obj_star_envelope
use mod_constants
use mod_global_parameters
use mod_obj_global_parameters
use mod_obj_mat
use mod_obj_dust
use mod_physics
use mod_obj_usr_unit
implicit none

type star_envelope_parameters
  character(len=20)    :: unit               !> physical unit at parameter file
  character(len=78)    :: obj_name           !> Obj name that call it
  integer              :: myindice       !> envelope  indices associated with envelope in use
  logical              :: normalize_done      !> envelope is the normalisation is already done
  character(len=20)    :: type_flow          !> envelope moving direction
  real(dp)             :: center(1:3)        !> envelope -star center
  real(dp)             :: r_in(1:3)          !> envelope inner boundary
  real(dp)             :: r_out (1:3)        !>  envelope outer  boundary
  real(dp)             :: lfac               !> envelope  Lorentz factor
  real(dp)             :: velocity(1:3)      !> envelope  speed
  real(dp)             :: magnetic(1:3)      !> envelope  magnetic field

  real(dp)             :: temperature        !> envelope  initial temperature
  real(dp)             :: pressure           !> envelope  initial pressure
  real(dp)             :: mass               !> envelope  initial mass
  real(dp)             :: density            !> envelope  initial density
  real(dp)             :: number_density     !> envelope  initial density
  real(dp)             :: xisigma0           !> envelope  magnetisation sigma

  real(dp)             :: kappa                !> envelope index power in pressure
  real(dp)             :: z_c                  !> envelope typical value of height  (cm)
  character(len=30)    :: profile_pressure     !> envelope profile pressure
  logical              :: profile_pressure_on  !> envelope pressure profile on
  logical              :: profile_force_on     !> envelope force profile one
  logical              :: profile_pressure_keep_density !> envelope density profile but keep density

  character(len=30)    :: profile_density     !> envelope profile pressure
  logical              :: profile_density_on  !> envelope density profile on
  logical              :: profile_density_keep_pressure !> envelope density profile but keep pressure

  logical              :: profile_on     !> envelope profile set
  integer              :: profile_idir   !> envelope profile direction
  character(len=20)    :: shape           !> envelope shape
  real(dp)             :: time_set           !> envelope  setting time

  logical              :: boundary_on         !> envelope logical to check if it will use implimented boundary
  character(len=30)    :: boundary_cond(3,2)  !> envelope boundary condition

  logical              :: reset_on            !> envelope reset
  real(dp)             :: reset_coef          !> envelope relaxation coefficient
  logical              :: gravity_on          !> envelope associate gration force on
  real(dp)             :: star_mass           !> envelope the mass of associated star

  logical              :: dust_on             !> logical to set dust
  real(dp)             :: dust_frac           !> dust fraction
  logical              :: tracer_on           !> tracer for the envelope
  integer              :: itr                 !> envelope indice
  real(dp)             :: tracer_init_density    !> envelope tracer initial density
  real(dp)             :: tracer_small_density   !> envelope tracer small density cut

end type




type star_envelope
  !Ref : Komissarov et Lyubarsky 2004,  Mon. Not. R. Astron. Soc. 349, 779–792 (2004)
  logical, allocatable              :: patch(:^D&)         !> envelope is in cell
  logical, allocatable              :: patch_escape(:^D&)  !> envelope  not in cell
  type(star_envelope_parameters)    :: myconfig            !> envelope paramter to read
  type(usrphysical_unit), pointer   :: myphysunit          !> envelope physicq unit in use

  type (dust)                       :: mydust                !> envelope dust
  type(usrboundary_type)            :: myboundaries          !> envelope boundary condition
  character(len=78)                 :: subname
  contains

   PROCEDURE, PASS(self)        :: set_default     => usr_envelope_set_default
   PROCEDURE, PASS(self)        :: set_complet     => usr_envelope_set_complet
   PROCEDURE, PASS(self)        :: normalize       => usr_envelope_normalize
   PROCEDURE, PASS(self)        :: set_w           => usr_envelope_set_w
   PROCEDURE, PASS(self)        :: read_parameters => usr_envelope_read_p
   PROCEDURE, PASS(self)        :: write_setting   => usr_envelope_write_setting
   PROCEDURE, PASS(self)        :: spd_rad_to_cart => usr_envelope_spd_rad_to_cart
   PROCEDURE, PASS(self)        :: clean_memory    => usr_envelope_clean_memory
   PROCEDURE, PASS(self)        :: get_dt          => usr_envelope_get_dt
   PROCEDURE, PASS(self)        :: alloc_set_patch => usr_envelope_alloc_set_patch

   PROCEDURE, PASS(self)        :: set_profile        => usr_envelope_set_profile
   PROCEDURE, PASS(self)        :: get_pforce_profile => usr_envelope_get_pforce_profile
   PROCEDURE, PASS(self)        :: add_source         => usr_envelope_add_source
end type star_envelope
contains

  !-------------------------------------------------------------------------
   !> Read the envelope parameters  from a parfile
    subroutine usr_envelope_read_p(self,envelope_config,files)
      class(star_envelope)                   :: self
      character(len=*), intent(in)           :: files(:)
      type(star_envelope_parameters)         :: envelope_config


    ! .. local ..
    integer                            :: i_file,i_error_read

    namelist /usr_envelope_list/  envelope_config
    namelist /usr_envelope1_list/ envelope_config
    namelist /usr_envelope2_list/ envelope_config
    namelist /usr_envelope3_list/ envelope_config

    if(mype==0)write(*,*)'Reading usr_envelope_list'
    do i_file = 1, size(files)
       open(unitpar, file=trim(files(i_file)), status="old")
       select case(envelope_config%myindice)
       case(1)
         read(unitpar, usr_envelope1_list, iostat=i_error_read)
       case(2)
         read(unitpar, usr_envelope2_list, iostat=i_error_read)
       case(3)
         read(unitpar, usr_envelope3_list, iostat=i_error_read)
       case default
         read(unitpar, usr_envelope_list, iostat=i_error_read)
       end select
       call usr_mat_read_error_message(i_error_read,envelope_config%myindice,&
                                       self%myconfig%obj_name)
       close(unitpar)
    end do



    if(envelope_config%boundary_on)then
      self%myboundaries%myconfig%myindice =1
      call self%myboundaries%read_parameters(self%myboundaries%myconfig,files)
    end if
    if(envelope_config%dust_on)then
      self%mydust%myconfig%associated_medium = 'envelope'
      call self%mydust%read_parameters(self%mydust%myconfig,files)
    end if


    end subroutine usr_envelope_read_p

  !------------------------------------------------------------------------
  !> write the cloud setting
  subroutine usr_envelope_write_setting(self,unit_config)
    implicit none
    class(star_envelope)                :: self
    integer,intent(in)                  :: unit_config
    ! .. local ..

    !-----------------------------------

    write(unit_config,*)'************************************'
    write(unit_config,*)'********** Envelope setting ********'
    write(unit_config,*)'************************************'
    write(unit_config,*)'      ****** Code Unit *******      '
    write(unit_config,*) 'Density     = ', self%myconfig%density,  'code unit'
    write(unit_config,*) 'Pressure    = ', self%myconfig%pressure,  'code unit'
    write(unit_config,*) 'Temperature = ', self%myconfig%temperature,  'code unit'
    write(unit_config,*) 'Speed       = ', self%myconfig%velocity,  'code unit'
    write(unit_config,*)'      ****** Physical Unit *******   '
    write(unit_config,*) 'Density     = ', self%myconfig%density*self%myphysunit%myconfig%density,&
                                            '  ',self%myphysunit%myunit%density
    write(unit_config,*) 'Pressure    = ', self%myconfig%pressure*self%myphysunit%myconfig%pressure,&
                                            '  ',self%myphysunit%myunit%pressure
    write(unit_config,*) 'Temperature = ', self%myconfig%temperature*self%myphysunit%myconfig%temperature,&
                                            '  ',self%myphysunit%myunit%temperature
    write(unit_config,*) 'Speed       = ', self%myconfig%velocity*self%myphysunit%myconfig%velocity,&
                                            '  ',self%myphysunit%myunit%velocity
    if(self%myconfig%dust_on) then
      call self%mydust%write_setting(unit_config)
     end if
    write(unit_config,*)'************************************'
    write(unit_config,*)'******* END Envelope setting *******'
    write(unit_config,*)'************************************'
  end    subroutine usr_envelope_write_setting
  !--------------------------------------------------------------------

  !> subroutine default setting for cloud
   subroutine usr_envelope_set_default(self)
    implicit none
    class(star_envelope)          :: self
    !----------------------------------
    self%myconfig%unit                   = 'code'
    self%myconfig%obj_name               = 'envelope'
    self%myconfig%myindice               = 0
    self%myconfig%mass                   = 0.0_dp
    self%myconfig%density                = 0.0_dp
    self%myconfig%number_density         = 0.0_dp
    self%myconfig%temperature            = 0.0_dp
    self%myconfig%pressure               = 0.0_dp
    self%myconfig%center(:)              = 0.0_dp!(box_limit(2,:)-box_limit(1,:))/2.0_dp
    self%myconfig%r_in                   = 0.0_dp!(box_limit(2,:)+box_limit(1,:))/2.0_dp
    self%myconfig%r_out                  = 0.0_dp
    self%myconfig%velocity(:)            = 0.0_dp
    self%myconfig%lfac                   = 0.0_dp
    self%myconfig%magnetic(:)            = 0.0_dp
    self%myconfig%xisigma0               = 0.0_dp



    self%myconfig%kappa                 = 0.0_dp
    self%myconfig%z_c                   = 0.0_dp
    self%myconfig%profile_pressure      = 'none'
    self%myconfig%profile_pressure_on   = .false.
    self%myconfig%profile_force_on      = .false.


    self%myconfig%profile_on            = .false.
    self%myconfig%profile_idir          = 0
    self%myconfig%profile_density       = 'none'
    self%myconfig%profile_density_on    = .false.

    self%myconfig%profile_pressure_keep_density =.true.
    self%myconfig%profile_density_keep_pressure =.true.
    self%myconfig%reset_coef            = 0.0_dp
    self%myconfig%reset_on              = .false.


    self%myconfig%shape                  = 'sphere'
    self%myconfig%time_set               = 0.0_dp

    self%myconfig%gravity_on             = .false.
    self%myconfig%star_mass              = 0.0_dp

    self%myconfig%boundary_on            = .false.
    self%myconfig%boundary_cond          = 'fix'

    self%myconfig%tracer_on              = .false.
    self%myconfig%itr                    = 0
    self%myconfig%tracer_init_density    = 0.0_dp
    self%myconfig%tracer_small_density   = 0.0_dp

    self%myconfig%normalize_done         = .false.


    call self%mydust%set_default

    call self%myboundaries%set_default

   end subroutine usr_envelope_set_default
   !--------------------------------------------------------------------
   !> subroutine check the parfile setting for cloud
   subroutine usr_envelope_set_complet(self)
     implicit none
     class(star_envelope)          :: self
     ! .. local ..
     logical                  :: dust_is_frac
     real(dp)                 :: mp,kb
     real(dp)                 :: envelope_volume
     !-----------------------------------
     if(SI_unit) then
       mp=mp_SI
       kB=kB_SI
     else
       mp=mp_cgs
       kB=kB_cgs
     end if


    call usr_envelope_get_volume_fprofile(self,envelope_volume)


    cond_isrelativiste: if(phys_config%isrel)then
      if (self%myconfig%lfac>0.0) then
       self%myconfig%velocity(r_)     = dsqrt(1.0_dp-1.0_dp/self%myconfig%lfac**2.0_dp)
       self%myconfig%velocity(theta_) = 0.0_dp
       self%myconfig%velocity(phi_)   = 0.0_dp
      select case(trim(self%myconfig%unit))
       case('code')
        ! do no thing
       case default
         self%myconfig%velocity=unit_velocity*self%myconfig%velocity
       end select
      else
       select case(trim(self%myconfig%unit))
        case('code')
         self%myconfig%lfac=1.0_dp/dsqrt(1.0_dp-SUM(self%myconfig%velocity**2.0_dp))
        case default
         self%myconfig%lfac=1.0_dp/dsqrt(1.0_dp-&
            SUM((self%myconfig%velocity/unit_velocity)**2.0_dp))
       end select
      end if
    end if cond_isrelativiste


    cond_mass_usr: if (dabs(self%myconfig%mass)>smalldouble)then
      self%myconfig%density        = self%myconfig%mass/envelope_volume
      self%myconfig%number_density = self%myconfig%density/mp
    else if(dabs(self%myconfig%density)<=0) then cond_mass_usr
      self%myconfig%density =self%myconfig%number_density*mp
      self%myconfig%mass    = self%myconfig%density*envelope_volume
    else if(dabs(self%myconfig%number_density)<=0) then cond_mass_usr
      self%myconfig%number_density =self%myconfig%density/mp
      self%myconfig%mass           = self%myconfig%density*envelope_volume
    end if cond_mass_usr

     if(dabs(self%myconfig%pressure)<=0.and.&
        dabs(self%myconfig%temperature)<=0) then
        write(*,*)' the temperature and pressure of envelope are not set'
        call mpistop('the code stop at envelope object in the user file')
     elseif(dabs(self%myconfig%pressure)<smalldouble) then
      self%myconfig%pressure =self%myconfig%number_density&
                                   *kB*self%myconfig%temperature
     end if

    select case(self%myconfig%profile_pressure)
    case('none')
     self%myconfig%profile_pressure_on = .false.
    case default
     if(.not.(dabs(self%myconfig%kappa)>smalldouble.and.self%myconfig%z_c>smalldouble)&
        .or.self%myconfig%profile_idir<1)then
       self%myconfig%profile_pressure_on = .false.
       self%myconfig%profile_force_on    = .false.
     end if
    end select

    select case(self%myconfig%profile_density)
    case('none')
     self%myconfig%profile_density_on = .false.
    case default
      if(self%myconfig%z_c<smalldouble.or.self%myconfig%profile_idir<1)then
       self%myconfig%profile_density_on = .false.
       self%myconfig%profile_force_on   = .false.
      end if
      if(self%myconfig%profile_force_on)self%myconfig%profile_pressure_on =.true.
    end select


    if( self%myconfig%dust_on)then
      dust_is_frac=.false.
      call  self%mydust%set_complet(dust_is_frac, self%myconfig%dust_frac,&
                              self%myconfig%density, self%myconfig%velocity)
    end if
    if(self%myconfig%profile_density_on.or.self%myconfig%profile_pressure_on)then
     self%myconfig%profile_on   = .true.
    end if

    if(.not.self%myconfig%boundary_on)                                       &
        self%myboundaries%myconfig%boundary_type=self%myconfig%boundary_cond

    call self%myboundaries%set_complet

    if(.not.self%myconfig%tracer_on)self%myconfig%itr=0
    if(self%myconfig%tracer_on)then
       prim_wnames(self%myconfig%itr+(phys_ind%tracer(1)-1)) = 'tracer_enveloppe'
       cons_wnames(self%myconfig%itr+(phys_ind%tracer(1)-1)) = 'tracer_enveloppe'
    end if


   end subroutine usr_envelope_set_complet
  !--------------------------------------------------------------------
   subroutine usr_envelope_normalize(self,physunit_inuse)
    use mod_obj_usr_unit
    implicit none
    class(star_envelope)                           :: self
    type(usrphysical_unit),target, intent(in)      :: physunit_inuse
    !----------------------------------
    self%myphysunit =>physunit_inuse
    if(trim(self%myconfig%unit)=='code'.or.self%myconfig%normalize_done)then
     if(self%myconfig%normalize_done)then
      write(*,*) 'WARNING: Second call for supernovae remnant normalisation ', &
                   'no new normalisation will be done'
     end if
     return
    end if
    self%myconfig%density     = self%myconfig%density     /physunit_inuse%myconfig%density
    self%myconfig%temperature = self%myconfig%temperature /physunit_inuse%myconfig%temperature
    self%myconfig%pressure    = self%myconfig%pressure    /physunit_inuse%myconfig%pressure
    self%myconfig%velocity    = self%myconfig%velocity    /physunit_inuse%myconfig%velocity
    self%myconfig%magnetic    = self%myconfig%magnetic    /physunit_inuse%myconfig%magnetic
    self%myconfig%center      = self%myconfig%center      /physunit_inuse%myconfig%length
    self%myconfig%r_in        = self%myconfig%r_in        /physunit_inuse%myconfig%length
    self%myconfig%r_out       = self%myconfig%r_out       /physunit_inuse%myconfig%length
    self%myconfig%time_set    = self%myconfig%time_set    /physunit_inuse%myconfig%time


   end subroutine usr_envelope_normalize
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
   !> subroutine patch for the cloud
   subroutine usr_envelope_patch(ixI^L,ixO^L,qt,x,self)
    implicit none
    integer, intent(in)        :: ixI^L,ixO^L
    real(kind=dp), intent(in)  :: qt
    real(kind=dp), intent(in)  :: x(ixI^S,1:ndim)
    class(star_envelope)       :: self
    real(dp), dimension(ixI^S) :: dist
    !----------------------------------
    if(allocated(self%patch))deallocate(self%patch)
    allocate(self%patch(ixI^S))

    call usr_distance(ixI^L,ixO^L,typeaxial,self%myconfig%center,x,dist)
    self%patch(ixO^S) = Dist(ixO^S) >self%myconfig%r_in(1) &
                        .and. Dist(ixO^S) <self%myconfig%r_out(1)
    if(allocated(self%patch_escape))then
      self%patch(ixO^S) = self%patch(ixO^S).and.(.not.self%patch_escape(ixO^S))
    end if
   end subroutine usr_envelope_patch
  !--------------------------------------------------------------------
   !> subroutine setting for cloud
   subroutine usr_envelope_set_w(ixI^L,ixO^L,qt,x,w,self,isboundary_iB)
    implicit none
    integer, intent(in)          :: ixI^L,ixO^L
    real(kind=dp), intent(in)    :: qt
    real(kind=dp)                :: x(ixI^S,1:ndim)
    real(kind=dp)                :: w(ixI^S,1:nw)
    class(star_envelope)         :: self
    integer,            optional :: isboundary_iB(2)
    ! .. local..
    integer                      :: idir,IB,idims,idims_bound,iw
    real(kind=dp)                :: fprofile(ixI^S)
    logical                      :: isboundary
    character(len=30)            :: myboundary_cond

    !----------------------------------

    call usr_envelope_patch(ixI^L,ixO^L,qt,x,self)
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

       w(ixO^S,phys_ind%mom(r_))        = self%myconfig%velocity(r_)
       w(ixO^S,phys_ind%mom(theta_))    = self%myconfig%velocity(theta_)
       w(ixO^S,phys_ind%mom(phi_))      = self%myconfig%velocity(phi_)
       !w(ixO^S,lfac_)          = self%myconfig%lfac
      end where

      where(self%patch(ixO^S))
       w(ixO^S,phys_ind%mag(r_))        = self%myconfig%magnetic(1)
       w(ixO^S,phys_ind%mag(theta_))    = self%myconfig%magnetic(2)
       w(ixO^S,phys_ind%mag(phi_))      = self%myconfig%magnetic(3)
      end where

      where(self%patch(ixO^S))
        w(ixO^S,phys_ind%rho_)           = self%myconfig%density
        w(ixO^S,phys_ind%pressure_)      = self%myconfig%pressure

      end where

    ! if typeaxial is not spherical project the spherical speed and magnetic field
    ! to geometry in use
    if(trim(typeaxial)/='spherical')call self%spd_rad_to_cart(ixI^L,ixO^L,x,w)

    ! add profile to envelope density and pressure
    if(self%myconfig%profile_on) then
       call self%set_profile(ixI^L,ixO^L,x,w)
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

    else boundary_cond
     if(any(self%patch(ixO^S))) then
       call self%myboundaries%set_w(ixI^L,ixO^L,iB,isboundary_iB(1),isboundary_iB(2),&
                                self%patch,x,w)
     end if
    end if boundary_cond
   end subroutine usr_envelope_set_w
 !--------------------------------------------------------------------

   subroutine usr_envelope_set_profile(ixI^L, ixO^L,x,w,self)
    implicit none
    integer, intent(in)           :: ixI^L,ixO^L
    real(kind=dp), intent(in)     :: x(ixI^S,1:ndim)
    real(kind=dp), intent(inout)  :: w(ixI^S,1:nw)
    class(star_envelope)                    :: self
    ! ..local ..

    real(kind=dp), dimension(ixI^S) :: distance,gravity_field, p_profile

    !----------------------------------------------------
    cond_pressure_profile : if(self%myconfig%profile_pressure_on) then
      if(index(trim(self%myconfig%profile_pressure),'radial')>0.and.&
          trim(typeaxial)/='spherical')then
        call usr_distance(ixI^L,ixO^L,typeaxial,self%myconfig%center,x,distance)
      else
        distance(ixO^S)  = x(ixO^S,self%myconfig%profile_idir)
      end if
     select case(trim(self%myconfig%profile_pressure))
      case('king')
        if(dabs(self%myconfig%kappa)>smalldouble.and.self%myconfig%z_c>smalldouble)then
          p_profile(ixO^S) = 1.0_dp/(1.0_dp+(distance(ixO^S)&
                            /self%myconfig%z_c)**(-self%myconfig%kappa) )
        else
          p_profile(ixO^S) = 1.0_dp
        end if
      case('komissarov')
        if(dabs(self%myconfig%kappa)>smalldouble.and.self%myconfig%z_c>smalldouble)then
          p_profile(ixO^S) = (distance(ixO^S)&
                              /self%myconfig%z_c)**(-self%myconfig%kappa)
        else
          p_profile(ixO^S) = 1.0_dp
        end if
      case default
        p_profile(ixO^S) = 1.0_dp
      end select
    end if   cond_pressure_profile


    cond_density_profile : if(self%myconfig%profile_density_on) then

      if(index(trim(self%myconfig%profile_density),'radial')>0.and.&
          trim(typeaxial)/='spherical')then
        call usr_distance(ixI^L,ixO^L,typeaxial,self%myconfig%center,x,distance)
      else
        distance(ixO^S)  = x(ixO^S,self%myconfig%profile_idir)
      end if
     select case(trim(self%myconfig%profile_density))
     case('cabrit1997')
        if(self%myconfig%z_c>smalldouble)then
          p_profile(ixO^S) = (1.0_dp+distance(ixO^S)/self%myconfig%z_c)**(-self%myconfig%kappa)
        else
          p_profile(ixO^S) = 1.0_dp
        end if
      case('radial_powerlaw')
        if(self%myconfig%z_c>smalldouble)then
          p_profile(ixO^S) = (distance(ixO^S)/self%myconfig%z_c)**(-self%myconfig%kappa)
        else
          p_profile(ixO^S) = 1.0_dp
        end if
      case default
        p_profile(ixO^S) = 1.0_dp
      end select
    end if   cond_density_profile

    if(self%myconfig%profile_pressure_on)then
      w(ixO^S,phys_ind%pressure_) = w(ixO^S,phys_ind%pressure_) * p_profile(ixO^S)
      if(.not.self%myconfig%profile_pressure_keep_density) then
        w(ixO^S,phys_ind%rho_)      = w(ixO^S,phys_ind%rho_) * p_profile(ixO^S)**(1.0_dp/phys_config%gamma)
      end if
    elseif(self%myconfig%profile_density_on)then
      w(ixO^S,phys_ind%rho_)      =  w(ixO^S,phys_ind%rho_)* p_profile(ixO^S)
      if(.not.self%myconfig%profile_density_keep_pressure)then
        if(self%myconfig%gravity_on) then
          gravity_field(ixO^S) = constusr%G*self%myconfig%star_mass/distance(ixO^S)
           w(ixO^S,phys_ind%pressure_) = w(ixO^S,phys_ind%pressure_)+&
                gravity_field(ixO^S)*distance(ixO^S) /&
                (1.0_dp+self%myconfig%kappa)* &
                w(ixO^S,phys_ind%rho_)*(1.0_dp - &
                (self%myconfig%r_out(1)/self%myconfig%r_in(1))**(-1-self%myconfig%kappa))
        else

        w(ixO^S,phys_ind%pressure_) =  w(ixO^S,phys_ind%pressure_)&
                                       *p_profile(ixO^S)**(phys_config%gamma)
        end if
      end if
    end if



   end subroutine usr_envelope_set_profile


   !--------------------------------------------------------------------

   subroutine usr_envelope_get_pforce_profile(ixI^L,ixO^L,qt,x,w,f_profile,self)
    implicit none
    integer, intent(in)           :: ixI^L,ixO^L
    real(kind=dp), intent(in)     :: qt
    real(kind=dp), intent(in)     :: x(ixI^S,1:ndim)
    real(kind=dp), intent(in)     :: w(ixI^S,1:nw)
    class(star_envelope)          :: self
    real(kind=dp), intent(inout)  :: f_profile(ixI^S,1:ndim)
    ! ..local ..
    real(kind=dp)                 :: distance(ixI^S)
    !----------------------------------------------------
    cond_pressure_fprofile : if(self%myconfig%profile_pressure_on) then
      if(index(trim(self%myconfig%profile_pressure),'radial')>0.and.&
         trim(typeaxial)/='spherical')then
        call usr_distance(ixI^L,ixO^L,typeaxial,self%myconfig%center,x,distance)
      else
        distance(ixO^S)  = x(ixO^S,self%myconfig%profile_idir)
      end if
      select case(trim(self%myconfig%profile_pressure))
      case('komissarov')
       cond_force1: if(dabs(self%myconfig%kappa)>smalldouble.and.self%myconfig%z_c>smalldouble)then
        where(self%patch(ixO^S))
         f_profile(ixO^S,self%myconfig%profile_idir) =- self%myconfig%kappa/self%myconfig%z_c*&
           (distance(ixO^S)/self%myconfig%z_c)**(-(self%myconfig%kappa+1))
        end  where
       else cond_force1
        where(self%patch(ixO^S))
         f_profile(ixO^S,self%myconfig%profile_idir) = 0.0_dp
        end where
       end if cond_force1
      case default
       where(self%patch(ixO^S))
        f_profile(ixO^S,self%myconfig%profile_idir) = 0.0_dp
       end where
      end select
    end if cond_pressure_fprofile

    cond_density_fprofile : if(self%myconfig%profile_density_on) then
      if(index(trim(self%myconfig%profile_density),'radial')>0.and.&
         trim(typeaxial)/='spherical')then
        call usr_distance(ixI^L,ixO^L,typeaxial,self%myconfig%center,x,distance)
      else
        distance(ixO^S)  = x(ixO^S,self%myconfig%profile_idir)
      end if
      select case(trim(self%myconfig%profile_density))
      case('cabrit1997')
       cond_force_rho1: if(self%myconfig%z_c>smalldouble)then
        where(self%patch(ixO^S))
         f_profile(ixO^S,self%myconfig%profile_idir) = -phys_config%gamma/self%myconfig%z_c*self%myconfig%kappa*&
           (1.0_dp+distance(ixO^S)/self%myconfig%z_c)**(-phys_config%gamma*self%myconfig%kappa-1)
        end  where
       else cond_force_rho1
        where(self%patch(ixO^S))
         f_profile(ixO^S,self%myconfig%profile_idir) = 0.0_dp
        end where
      end if cond_force_rho1
      case('radial_powerlaw')
       cond_force_rho2: if(self%myconfig%z_c>smalldouble)then
        where(self%patch(ixO^S))
         f_profile(ixO^S,self%myconfig%profile_idir) = &
           -phys_config%gamma/self%myconfig%z_c*self%myconfig%kappa*&
           (1.0_dp+distance(ixO^S)/self%myconfig%z_c)**(-phys_config%gamma*self%myconfig%kappa-1.0_dp)
        end  where
      else cond_force_rho2
        where(self%patch(ixO^S))
         f_profile(ixO^S,self%myconfig%profile_idir) = 0.0_dp
        end where
      end if cond_force_rho2
      case default
       where(self%patch(ixO^S))
        f_profile(ixO^S,self%myconfig%profile_idir) = 0.0_dp
       end where
      end select
    end if cond_density_fprofile
   end subroutine usr_envelope_get_pforce_profile

   !--------------------------------------------------------------------
   subroutine usr_envelope_add_source(ixI^L,ixO^L,iw^LIM,x,qdt,qtC,wCT,qt,w,self,&
                                 use_tracer,patch_escape,source_filter)
     implicit none
     integer, intent(in)                     :: ixI^L,ixO^L,iw^LIM
     real(kind=dp), intent(in)               :: qdt,qtC,qt
     real(kind=dp), intent(in)               :: x(ixI^S,1:ndim)
     real(kind=dp), intent(in)               :: wCT(ixI^S,1:nw)
     real(kind=dp), intent(inout)            :: w(ixI^S,1:nw)
     logical, intent(in), optional           :: use_tracer
     logical, intent(in),optional            :: patch_escape(ixI^S)
     real(kind=dp), intent(in),optional      :: source_filter(ixI^S)

     class(star_envelope)                    :: self
     ! .. local ..
     real(kind=dp)                           :: source_filter_loc(ixI^S)
     real(kind=dp)                           :: f_profile(ixI^S,1:ndim)
     real(kind=dp)                           :: w_init(ixI^S,1:nw)
     integer                                 :: idir,i_idir_prof_,imom_profile_
     !---------------------------------------------------------



      call self%alloc_set_patch(ixI^L,ixO^L,qt,x,&
                                use_tracer=use_tracer,w=w,&
                                patch_escape=patch_escape)

      cond_add_force : if(self%myconfig%profile_force_on) then

        cond_inside_prof: if(any(self%patch(ixO^S)))then
          call self%get_pforce_profile(ixI^L,ixO^L,qt,x,wCT,f_profile)
          i_idir_prof_  = self%myconfig%profile_idir
          imom_profile_ =phys_ind%mom(i_idir_prof_)

          where(self%patch(ixO^S))
            w(ixO^S,imom_profile_) = w(ixO^S,imom_profile_)+&
                    qdt*wCT(ixO^S,phys_ind%rho_)*f_profile(ixO^S,i_idir_prof_)
          end where

         !if(energy .and. .not.block%e_is_internal) then
          where(self%patch(ixO^S))
          w(ixO^S,phys_ind%e_)=w(ixO^S,phys_ind%e_) &
              + qdt * f_profile(ixO^S,i_idir_prof_) * wCT(ixO^S,imom_profile_)!/wCT(ixO^S,phys_ind%rho_)
          end where
         !end if
        end if  cond_inside_prof
      end if cond_add_force



     cond_reset : if(self%myconfig%reset_coef>0.0_dp)then
      cond_inside : if(any(self%patch(ixO^S)))then

       cond_filter : if(present(source_filter))then
         source_filter_loc(ixO^S) = source_filter(ixO^S)
       else cond_filter
         source_filter_loc(ixO^S) = 1.0_dp
       end if cond_filter


       call self%set_w(ixI^L,ixO^L,qt,x,w_init)
       call phys_to_primitive(ixI^L,ixO^L,w,x)

       where(self%patch(ixO^S))
         source_filter_loc(ixO^S) = max(dabs(source_filter_loc(ixO^S)*self%myconfig%reset_coef),1.0_dp)

         w(ixO^S,phys_ind%rho_) =   w(ixO^S,phys_ind%rho_)*(1.0_dp-source_filter_loc(ixO^S)) +&
                          w_init(ixO^S,phys_ind%rho_)*source_filter_loc(ixO^S)
         w(ixO^S,phys_ind%pressure_) =  w(ixO^S,phys_ind%pressure_)*(1.0_dp-source_filter_loc(ixO^S)) +&
                          w_init(ixO^S,phys_ind%pressure_)*source_filter_loc(ixO^S)
       end where
       loop_idir :  do idir = 1,ndir
        where(self%patch(ixO^S))
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

       call phys_to_conserved(ixI^L,ixO^L,w,x)

      end if cond_inside
    end if cond_reset
    !if(any(self%patch(ixO^S)))print*,' test envelope force',maxval(dabs(w(ixO^S,phys_ind%mom(z_))),mask=self%patch(ixO^S))
   end subroutine usr_envelope_add_source

   !--------------------------------------------------------------------
   !===========================================================
   !> Subroutine to set time set before envelope  starts
   subroutine usr_envelope_get_dt(self,ixI^L,ixO^L,dx^D,x,w,qt,dtnew)
     use mod_global_parameters
     class(star_envelope)            :: self
     integer, intent(in)             :: ixI^L, ixO^L
     double precision, intent(in)    :: dx^D,qt, x(ixI^S,1:ndim)
     double precision, intent(in)    :: w(ixI^S,1:nw)
     double precision, intent(inout) :: dtnew
     !--------------------------------------------------------------
     if(qt<self%myconfig%time_set)then
      dtnew=min(self%myconfig%time_set-qt,dtnew)
     end if
   end subroutine usr_envelope_get_dt
!===========================================================


  subroutine  usr_envelope_get_volume_fprofile(self,envelope_volume)
    use mod_global_parameters
    implicit none

    class(star_envelope)            :: self
    real(kind=dp)                   :: envelope_volume
    ! .. local ..
    integer                         :: geo_power
    real(kind=dp)                   :: r_inner,r_out,geo_factor
    !-------------------------------------------------------------
    cond_density_profile : if(self%myconfig%profile_density_on) then
      select case(trim(self%myconfig%shape))
      case('cube')
        geo_power  = 0
        geo_factor = 1.0
      case('cylinder')
        geo_power  = 1
        geo_factor = 2.0_dp*dpi
      case('sphere')
        geo_power  = 2
        geo_factor = 4.0_dp*dpi
      end select


      r_inner  = self%myconfig%r_in(1)
      r_out    = self%myconfig%r_out(1)

     select case(trim(self%myconfig%profile_density))

      case('radial_powerlaw')
        if(self%myconfig%z_c/=geo_power)then
          envelope_volume = geo_factor/&
               ((geo_power-self%myconfig%kappa)*self%myconfig%z_c**geo_power)*&
            ((self%myconfig%r_out(1)/self%myconfig%z_c)**(geo_power-self%myconfig%kappa) &
            -(self%myconfig%r_in(1)/self%myconfig%z_c)**(geo_power-self%myconfig%kappa))
        else
          envelope_volume = geo_factor/&
             ((geo_power-self%myconfig%kappa)&
               *self%myconfig%z_c**(-self%myconfig%kappa))  &
               * dlog(self%myconfig%r_out(1)/self%myconfig%r_in(1))
        end if
      case default
        envelope_volume  = 1.0_dp
      end select
    else    cond_density_profile
      call usr_get_volume((/r_inner,r_out/),'spherical', envelope_volume)




    end if   cond_density_profile
  end subroutine  usr_envelope_get_volume_fprofile
  !===========================================================
  !> Subroutine to convert radial wind from spherical to cartesian coordinates
  subroutine usr_envelope_spd_rad_to_cart(ixI^L,ixO^L,x,w,self)
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: ixI^L,ixO^L
    real(kind=dp), intent(in)       :: x(ixI^S,1:ndim)
    real(kind=dp), intent(inout)    :: w(ixI^S,1:nw)
    class(star_envelope)            :: self
    ! .. local ..
    real(kind=dp), dimension(ixI^S,1:ndir) :: v_spherical
    real(kind=dp), dimension(ixI^S)        :: sin_theta,cos_theta
    real(kind=dp), dimension(ixI^S)        :: sin_phi,cos_phi
    real(kind=dp), dimension(ixI^S)        :: Dist
    !--------------------------------------------------

    call usr_distance(ixI^L,ixO^L,typeaxial,self%myconfig%center,x,dist)

   select case(typeaxial)

    case('slab')
     v_spherical(ixO^S,1:ndir) = w(ixO^S,phys_ind%mom(1):phys_ind%mom(ndir))
     if(z_in) then
       sin_theta(ixO^S)=dsqrt(x(ixO^S,x_)**2.0_dp+x(ixO^S,y_)**2.0_dp)/Dist(ixO^S)
       cos_theta(ixO^S)=x(ixO^S,z_)/Dist(ixO^S)
     else
      sin_theta(ixO^S)=1.0_DP
      cos_theta(ixO^S)=0.0_dp
     end if
     if(y_in) then
       sin_phi(ixO^S)=dsin(x(ixO^S,phi_))
       cos_phi(ixO^S)=dcos(x(ixO^S,phi_))
     else
      sin_phi(ixO^S)=0.0_DP
      cos_phi(ixO^S)=1.0_dp
     end if
     where(self%patch(ixO^S))
      w(ixO^S,phys_ind%mom(x_)) = v_spherical(ixO^S,r_)*cos_phi(ixO^S)*sin_theta(ixO^S)

      where(sin_phi(ixO^S)>0.0_dp)
       w(ixO^S,phys_ind%mom(y_)) = v_spherical(ixO^S,r_)*sin_phi(ixO^S)*sin_theta(ixO^S)
      else where
       w(ixO^S,phys_ind%mom(y_)) = zero
      end where
      where(cos_phi(ixO^S)>0.0_dp)
       w(ixO^S,phys_ind%mom(z_)) = v_spherical(ixO^S,r_)*cos_theta(ixO^S)
      else where
       w(ixO^S,phys_ind%mom(z_)) = zero
      end where

     end where
    case('cylindrical')
     v_spherical(ixO^S,1:ndir) = w(ixO^S,phys_ind%mom(1):phys_ind%mom(ndir))
     if(z_in) then
       where(Dist(ixO^S)>0.0_dp)
        sin_theta(ixO^S)=x(ixO^S,r_)/Dist(ixO^S)
        cos_theta(ixO^S)=x(ixO^S,z_)/Dist(ixO^S)
       elsewhere
        sin_theta(ixO^S)=1.0_DP
        cos_theta(ixO^S)=0.0_dp
       endwhere
     else
      sin_theta(ixO^S)=1.0_DP
      cos_theta(ixO^S)=0.0_dp
     end if
     if(phi_in) then
       sin_phi(ixO^S)=dsin(x(ixO^S,phi_))
       cos_phi(ixO^S)=dcos(x(ixO^S,phi_))
     else
      sin_phi(ixO^S)=0.0_DP
      cos_phi(ixO^S)=1.0_dp
     end if
     where(self%patch(ixO^S))
      w(ixO^S,phys_ind%mom(x_)) = v_spherical(ixO^S,r_)*cos_phi(ixO^S)*sin_theta(ixO^S)
     end where
     !  where(sin_phi(ixO^S)>0.0_dp.and.self%patch(ixO^S))
     !   w(ixO^S,phys_ind%mom(y_)) = w(ixO^S,phys_ind%mom(r_))*sin_phi(ixO^S)*sin_theta(ixO^S)
     ! else where(sin_phi(ixO^S)<=0.0_dp.and.self%patch(ixO^S))
     !   w(ixO^S,phys_ind%mom(y_)) = zero
     !  end where
      where(cos_phi(ixO^S)>0.0_dp.and.self%patch(ixO^S))
       w(ixO^S,phys_ind%mom(z_)) = v_spherical(ixO^S,r_)*cos_theta(ixO^S)
      else where(cos_phi(ixO^S)<=0.0_dp.and.self%patch(ixO^S))
       w(ixO^S,phys_ind%mom(z_)) = zero
      end where
    case('spherical')
    ! Dummy
    case default
    write(*,*) ' is not implimented '
    call mpistop(' stop at mod_usr.t  et usr_envelope_spd_rad_to_cart')
   end select
  end subroutine usr_envelope_spd_rad_to_cart

  !--------------------------------------------------------------------
   !> Subroutine to clean array memory of associated with cloud object
   subroutine usr_envelope_alloc_set_patch(ixI^L,ixO^L,qt,x,self,use_tracer,w,patch_escape)
     implicit none
     integer, intent(in)                     :: ixI^L,ixO^L
     real(kind=dp), intent(in)               :: qt
     real(kind=dp), intent(in)               :: x(ixI^S,1:ndim)
     real(kind=dp), intent(in), optional     :: w(ixI^S,1:nw)
     logical, intent(in), optional           :: use_tracer
     logical, intent(in),optional            :: patch_escape(ixI^S)
     class(star_envelope)                              :: self
     !---------------------------------------------------------
     cond_tracer : if(.not.present(use_tracer).and. .not.present(w)) then
       if(allocated(self%patch))deallocate(self%patch)
       allocate(self%patch(ixI^S))
       self%patch(ixO^S) =.true.

     else cond_tracer

       cond_envelopetracer_on : if(self%myconfig%tracer_on)then
        if(allocated(self%patch))deallocate(self%patch)
        allocate(self%patch(ixI^S))
        where(w(ixO^S,phys_ind%tracer(self%myconfig%itr))>small_density)
          self%patch(ixO^S)=.true.
        else where
          self%patch(ixO^S)=.false.
        end where
       end if cond_envelopetracer_on

     end if cond_tracer

     if(allocated(self%patch_escape))deallocate(self%patch_escape)
     allocate(self%patch_escape(ixI^S))
     if(present(patch_escape))then
      self%patch_escape(ixO^S) = patch_escape(ixO^S)
     else
      self%patch_escape(ixO^S) =.false.
     end if

     where(self%patch(ixO^S))self%patch(ixO^S)         = .not.self%patch_escape(ixO^S)
   end subroutine usr_envelope_alloc_set_patch
  !---------------------------------------------------------------------
   subroutine  usr_envelope_get_patch_escape(ixI^L,ixO^L,need_dealloc,patch_escape,self)
     implicit none
     integer, intent(in)           :: ixI^L,ixO^L
     logical, intent(in)           :: need_dealloc
     logical, intent(in)           :: patch_escape(ixI^S)
     class(star_envelope)                    :: self
     !----------------------------------------------------------
     if(allocated(self%patch_escape))deallocate(self%patch_escape)
      allocate(self%patch_escape(ixI^S))
      self%patch_escape(ixO^S) = .false.

     self%patch_escape(ixO^S)=self%patch_escape(ixO^S).or.patch_escape(ixO^S)
   end subroutine  usr_envelope_get_patch_escape
   !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !> Subroutine to clean array memory of associated with cloud object
  subroutine usr_envelope_clean_memory(self)
    class(star_envelope)    :: self
    if(allocated(self%patch_escape))deallocate(self%patch_escape)
    if(allocated(self%patch))deallocate(self%patch)
  end subroutine usr_envelope_clean_memory
  !--------------------------------------------------------------------
end module mod_obj_star_envelope
