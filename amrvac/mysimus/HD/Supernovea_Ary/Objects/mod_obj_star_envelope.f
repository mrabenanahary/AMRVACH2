module mod_obj_star_envelope
use mod_constants
use mod_global_parameters
use mod_obj_global_parameters
use mod_obj_mat
use mod_physics
use mod_srmhd_parameters
use mod_obj_usr_unit
implicit none

type star_envelope_parameters
  character(len=20)    :: unit !> physical unit at parameter file
  logical              :: normalize_done !> envelope is the normalisation is already done
  character(len=20)    :: type_flow          !> envelope moving direction
  real(dp)             :: center(1:3)        !> envelope -star center
  real(dp)             :: r_in(1:3)          !> envelope inner boundary
  real(dp)             :: r_out (1:3)        !>  envelope outer  boundary
  real(dp)             :: lfac               !> envelope  Lorentz factor
  real(dp)             :: velocity(1:3)      !> envelope  speed
  real(dp)             :: magnetic(1:3)      !> envelope  magnetic field
  real(dp)             :: temperature_init   !> envelope  initial temperature
  real(dp)             :: pressure_init      !> envelope  initial pressure
  real(dp)             :: density_init       !> envelope  initial density
  real(dp)             :: number_density_init!> envelope  initial density
  real(dp)             :: xisigma0           !> envelope  magnetisation sigma
  real(dp)             :: time_set           !> envelope  setting time
  logical              :: tracer_on          !> tracer for the envelope
  integer              :: itr                !> envelope indice
end type




type star_envelope
  !Ref : Komissarov et Lyubarsky 2004,  Mon. Not. R. Astron. Soc. 349, 779–792 (2004)
  logical, allocatable              :: patch(:)         !> envelope is in cell
  logical, allocatable              :: patch_escape(:) !> envelope  not in cell
  type(star_envelope_parameters)    :: myconfig !> envelope paramter to read
  type(usrphysical_unit), pointer   :: myphysunit !> envelope physicq unit in use
  character(len=78)                 :: subname
  contains

   PROCEDURE, PASS(self)        :: set_default     => usr_envelope_set_default
   PROCEDURE, PASS(self)        :: set_complet     => usr_envelope_set_complet
   PROCEDURE, PASS(self)        :: normalize       => usr_envelope_normalize
   PROCEDURE, PASS(self)        :: set_w           => usr_envelope_set_w
   PROCEDURE, PASS(self)        :: read_parameters => usr_envelope_read_p
   PROCEDURE, PASS(self)        :: write_setting   => &
      usr_envelope_write_setting
   PROCEDURE, PASS(self)        :: spd_rad_to_cart => &
      usr_envelope_spd_rad_to_cart
   PROCEDURE, PASS(self)        :: clean_memory    => &
      usr_envelope_clean_memory
   PROCEDURE, PASS(self)        :: get_dt          => usr_envelope_get_dt
end type star_envelope
contains

  !-------------------------------------------------------------------------
   !> Read the ism parameters  from a parfile
    subroutine usr_envelope_read_p(self,envelope_config,files)
      class(star_envelope)                   :: self
      character(len=*), intent(in)           :: files(:)
      type(star_envelope_parameters)         :: envelope_config

      integer  :: i_file
      namelist /usr_envelope_list/ envelope_config

      if(mype==0)write(*,*)'Reading usr_envelope_list'
      do i_file = 1, size(files)
         open(unitpar, file=trim(files(i_file)), status="old")
         read(unitpar, usr_envelope_list)
         close(unitpar)
      end do


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
    write(unit_config,*) 'Density     = ', self%myconfig%density_init,&
         'code unit'
    write(unit_config,*) 'Pressure    = ', self%myconfig%pressure_init,&
         'code unit'
    write(unit_config,*) 'Temperature = ', self%myconfig%temperature_init,&
         'code unit'
    write(unit_config,*) 'Speed       = ', self%myconfig%velocity,&
         'code unit'
    write(unit_config,*)'      ****** Physical Unit *******   '
    write(unit_config,*) 'Density     = ',&
        self%myconfig%density_init*self%myphysunit%myconfig%density,'  ',&
       self%myphysunit%myunit%density
    write(unit_config,*) 'Pressure    = ',&
        self%myconfig%pressure_init*self%myphysunit%myconfig%pressure,'  ',&
       self%myphysunit%myunit%pressure
    write(unit_config,*) 'Temperature = ',&
        self%myconfig%temperature_init*self%myphysunit%myconfig%temperature,&
       '  ',self%myphysunit%myunit%temperature
    write(unit_config,*) 'Speed       = ',&
        self%myconfig%velocity*self%myphysunit%myconfig%velocity,'  ',&
       self%myphysunit%myunit%velocity
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
    self%myconfig%density_init           = 0.0_dp
    self%myconfig%number_density_init    = 0.0_dp
    self%myconfig%temperature_init       = 0.0_dp
    self%myconfig%pressure_init          = 0.0_dp
    self%myconfig%center(:)              = 0.0_dp !(box_limit(2,:)-box_limit(1,:))/2.0_dp
    self%myconfig%r_in                   = 0.0_dp !(box_limit(2,:)+box_limit(1,:))/2.0_dp
    self%myconfig%r_out                  = 0.0_dp
    self%myconfig%velocity(:)            = 0.0_dp
    self%myconfig%lfac                   = 0.0_dp
    self%myconfig%magnetic(:)            = 0.0_dp
    self%myconfig%xisigma0               = 0.0_dp
    self%myconfig%time_set               = 0.0_dp

    self%myconfig%tracer_on              = .false.
    self%myconfig%itr                    = 0
    self%myconfig%normalize_done         = .false.
   end subroutine usr_envelope_set_default
   !--------------------------------------------------------------------
   !> subroutine check the parfile setting for cloud
   subroutine usr_envelope_set_complet(self)
     implicit none
     class(star_envelope)          :: self
     ! .. local ..
     real(dp)                 :: mp,kb
     !-----------------------------------
     if(SI_unit) then
       mp=mp_SI
       kB=kB_SI
     else
       mp=mp_cgs
       kB=kB_cgs
     end if



     if (self%myconfig%lfac>0.0) then
       self%myconfig%velocity(r_)     = dsqrt(1.0_dp-&
          1.0_dp/self%myconfig%lfac**2.0_dp)
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
         self%myconfig%lfac=1.0_dp/dsqrt(1.0_dp-&
            SUM(self%myconfig%velocity**2.0_dp))
        case default
         self%myconfig%lfac=1.0_dp/dsqrt(1.0_dp-&
            SUM((self%myconfig%velocity/unit_velocity)**2.0_dp))
       end select
     end if


     if(dabs(self%myconfig%density_init)<=0) then
      self%myconfig%density_init =self%myconfig%number_density_init*mp
     end if

     if(dabs(self%myconfig%number_density_init)<=0) then
      self%myconfig%number_density_init =self%myconfig%density_init/mp
     end if

     if(dabs(self%myconfig%pressure_init)<=&
        0.and.dabs(self%myconfig%temperature_init)<=0) then
        write(*,*)' the temperature and pressure of envelope are not set'
        call mpistop('the code stop at envelope object in the user file')
     elseif(dabs(self%myconfig%pressure_init)<smalldouble) then
      self%myconfig%pressure_init&
          =self%myconfig%number_density_init*kB*self%myconfig%temperature_init
     end if
     if(.not.self%myconfig%tracer_on)self%myconfig%itr=0
    if(self%myconfig%tracer_on)then
       prim_wnames(self%myconfig%itr+(phys_ind%tracer(1)-1)) = &
          'tracer_enveloppe'
       cons_wnames(self%myconfig%itr+(phys_ind%tracer(1)-1)) = &
          'tracer_enveloppe'
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
      write(*,*) 'WARNING: Second call for supernovae remnant normalisation ',&
          'no new normalisation will be done'
     end if
     return
    end if
    self%myconfig%density_init     = self%myconfig%density_init     &
       /physunit_inuse%myconfig%density
    self%myconfig%temperature_init = self%myconfig%temperature_init &
       /physunit_inuse%myconfig%temperature
    self%myconfig%pressure_init    = self%myconfig%pressure_init    &
       /physunit_inuse%myconfig%pressure
    self%myconfig%velocity         = self%myconfig%velocity         &
       /physunit_inuse%myconfig%velocity
    self%myconfig%magnetic         = self%myconfig%magnetic         &
       /physunit_inuse%myconfig%magnetic
    self%myconfig%center           = self%myconfig%center           &
       /physunit_inuse%myconfig%length
    self%myconfig%r_in             = self%myconfig%r_in             &
       /physunit_inuse%myconfig%length
    self%myconfig%r_out            = self%myconfig%r_out            &
       /physunit_inuse%myconfig%length
    self%myconfig%time_set         = self%myconfig%time_set         &
       /physunit_inuse%myconfig%time


   end subroutine usr_envelope_normalize
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
   !> subroutine patch for the cloud
   subroutine usr_envelope_patch(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,x,self)
    implicit none
    integer, intent(in)        :: ixImin1,ixImax1,ixOmin1,ixOmax1
    real(kind=dp), intent(in)  :: qt
    real(kind=dp), intent(in)  :: x(ixImin1:ixImax1,1:ndir)
    class(star_envelope)       :: self
    real(dp), dimension(ixImin1:ixImax1) :: dist
    !----------------------------------

    allocate(self%patch(ixGlo1:ixGhi1))

    call usr_distance(ixImin1,ixImax1,ixOmin1,ixOmax1,typeaxial,&
       self%myconfig%center,x,dist)
    self%patch(ixOmin1:ixOmax1) = Dist(ixOmin1:ixOmax1) >self%myconfig%r_in(1) &
       .and. Dist(ixOmin1:ixOmax1) <self%myconfig%r_out(1)
    if(allocated(self%patch_escape))then
      self%patch(ixOmin1:ixOmax1) = self%patch(ixOmin1:ixOmax1).and.&
         (.not.self%patch_escape(ixOmin1:ixOmax1))
    end if
   end subroutine usr_envelope_patch
  !--------------------------------------------------------------------
   !> subroutine setting for cloud
   subroutine usr_envelope_set_w(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,x,w,self)
    implicit none
    integer, intent(in)          :: ixImin1,ixImax1,ixOmin1,ixOmax1
    real(kind=dp), intent(in)    :: qt
    real(kind=dp)                :: x(ixImin1:ixImax1,1:ndir)
    real(kind=dp)                :: w(ixImin1:ixImax1,1:nw)
    class(star_envelope)         :: self
    ! .. local..
    integer                      :: idir

    !----------------------------------

    call usr_envelope_patch(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,x,self)


    where(self%patch(ixOmin1:ixOmax1))

     w(ixOmin1:ixOmax1,mom(r_))        = self%myconfig%velocity(r_)
     w(ixOmin1:ixOmax1,mom(theta_))    = self%myconfig%velocity(theta_)
     w(ixOmin1:ixOmax1,mom(phi_))      = self%myconfig%velocity(phi_)
     w(ixOmin1:ixOmax1,lfac_)          = self%myconfig%lfac
  end where

    where(self%patch(ixOmin1:ixOmax1))
     w(ixOmin1:ixOmax1,mag(r_))        = self%myconfig%magnetic(1)
     w(ixOmin1:ixOmax1,mag(theta_))    = self%myconfig%magnetic(2)
     w(ixOmin1:ixOmax1,mag(phi_))      = self%myconfig%magnetic(3)
  end where

    where(self%patch(ixOmin1:ixOmax1))
      w(ixOmin1:ixOmax1,rho_)           = self%myconfig%density_init
      w(ixOmin1:ixOmax1,p_)             = self%myconfig%pressure_init

    end where

    if(trim(typeaxial)/='spherical')call self%spd_rad_to_cart(ixImin1,ixImax1,&
       ixOmin1,ixOmax1,x,w)

    cond_tracer_on :if(self%myconfig%tracer_on.and.phys_config%n_tracer>0.and.&
       itr<=phys_config%n_tracer)then
     where(self%patch(ixOmin1:ixOmax1))
      w(ixOmin1:ixOmax1,phys_ind%tracer(itr)) = 1.0d2 !w(ixOmin1:ixOmax1,rho_)
     elsewhere
      w(ixOmin1:ixOmax1,phys_ind%tracer(itr)) = 0.0_dp
     end where
     itr=itr+1
    end if cond_tracer_on



   end subroutine usr_envelope_set_w

   !===========================================================
   !> Subroutine to set time set before envelope  starts
   subroutine usr_envelope_get_dt(self,ixImin1,ixImax1,ixOmin1,ixOmax1,dx1,x,w,&
      qt,dtnew)
     use mod_global_parameters
     class(star_envelope)            :: self
     integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
     double precision, intent(in)    :: dx1,qt, x(ixImin1:ixImax1,1:ndim)
     double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)
     double precision, intent(inout) :: dtnew
     !--------------------------------------------------------------
     if(qt<self%myconfig%time_set)then
      dtnew=min(self%myconfig%time_set-qt,dtnew)
     end if
   end subroutine usr_envelope_get_dt
!===========================================================
  !===========================================================
  !> Subroutine to convert radial wind from spherical to cartesian coordinates
  subroutine usr_envelope_spd_rad_to_cart(ixImin1,ixImax1,ixOmin1,ixOmax1,x,w,&
     self)
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
    real(kind=dp), intent(in)       :: x(ixImin1:ixImax1,1:ndir)
    real(kind=dp), intent(inout)    :: w(ixImin1:ixImax1,1:nw)
    class(star_envelope)            :: self
    ! .. local ..
    real(kind=dp), dimension(ixImin1:ixImax1,1:ndir) :: v_spherical
    real(kind=dp), dimension(ixImin1:ixImax1)        :: sin_theta,cos_theta
    real(kind=dp), dimension(ixImin1:ixImax1)        :: sin_phi,cos_phi
    real(kind=dp), dimension(ixImin1:ixImax1)        :: Dist
    !--------------------------------------------------

    call usr_distance(ixImin1,ixImax1,ixOmin1,ixOmax1,typeaxial,&
       self%myconfig%center,x,dist)

   select case(typeaxial)

    case('slab')
     v_spherical(ixOmin1:ixOmax1,1:ndir) = w(ixOmin1:ixOmax1,mom(1):mom(ndir))
     if(z_in) then
       sin_theta(ixOmin1:ixOmax1)=dsqrt(x(ixOmin1:ixOmax1,&
          x_)**2.0_dp+x(ixOmin1:ixOmax1,y_)**2.0_dp)/Dist(ixOmin1:ixOmax1)
       cos_theta(ixOmin1:ixOmax1)=x(ixOmin1:ixOmax1,z_)/Dist(ixOmin1:ixOmax1)
     else
      sin_theta(ixOmin1:ixOmax1)=1.0_DP
      cos_theta(ixOmin1:ixOmax1)=0.0_dp
     end if
     if(y_in) then
       sin_phi(ixOmin1:ixOmax1)=dsin(x(ixOmin1:ixOmax1,phi_))
       cos_phi(ixOmin1:ixOmax1)=dcos(x(ixOmin1:ixOmax1,phi_))
     else
      sin_phi(ixOmin1:ixOmax1)=0.0_DP
      cos_phi(ixOmin1:ixOmax1)=1.0_dp
     end if
     where(self%patch(ixOmin1:ixOmax1))
      w(ixOmin1:ixOmax1,mom(x_)) = v_spherical(ixOmin1:ixOmax1,&
         r_)*cos_phi(ixOmin1:ixOmax1)*sin_theta(ixOmin1:ixOmax1)

      where(sin_phi(ixOmin1:ixOmax1)>0.0_dp)
       w(ixOmin1:ixOmax1,mom(y_)) = v_spherical(ixOmin1:ixOmax1,&
          r_)*sin_phi(ixOmin1:ixOmax1)*sin_theta(ixOmin1:ixOmax1)
      else where
       w(ixOmin1:ixOmax1,mom(y_)) = zero
      end where
      where(cos_phi(ixOmin1:ixOmax1)>0.0_dp)
       w(ixOmin1:ixOmax1,mom(z_)) = v_spherical(ixOmin1:ixOmax1,&
          r_)*cos_theta(ixOmin1:ixOmax1)
      else where
       w(ixOmin1:ixOmax1,mom(z_)) = zero
      end where

     end where
    case('cylindrical')
     v_spherical(ixOmin1:ixOmax1,1:ndir) = w(ixOmin1:ixOmax1,mom(1):mom(ndir))
     if(z_in) then
       where(Dist(ixOmin1:ixOmax1)>0.0_dp)
        sin_theta(ixOmin1:ixOmax1)=x(ixOmin1:ixOmax1,r_)/Dist(ixOmin1:ixOmax1)
        cos_theta(ixOmin1:ixOmax1)=x(ixOmin1:ixOmax1,z_)/Dist(ixOmin1:ixOmax1)
       elsewhere
        sin_theta(ixOmin1:ixOmax1)=1.0_DP
        cos_theta(ixOmin1:ixOmax1)=0.0_dp
       endwhere
     else
      sin_theta(ixOmin1:ixOmax1)=1.0_DP
      cos_theta(ixOmin1:ixOmax1)=0.0_dp
     end if
     if(phi_in) then
       sin_phi(ixOmin1:ixOmax1)=dsin(x(ixOmin1:ixOmax1,phi_))
       cos_phi(ixOmin1:ixOmax1)=dcos(x(ixOmin1:ixOmax1,phi_))
     else
      sin_phi(ixOmin1:ixOmax1)=0.0_DP
      cos_phi(ixOmin1:ixOmax1)=1.0_dp
     end if
     where(self%patch(ixOmin1:ixOmax1))
      w(ixOmin1:ixOmax1,mom(x_)) = v_spherical(ixOmin1:ixOmax1,&
         r_)*cos_phi(ixOmin1:ixOmax1)*sin_theta(ixOmin1:ixOmax1)
     end where
     !  where(sin_phi(ixO^S)>0.0_dp.and.self%patch(ixO^S))
     !   w(ixO^S,mom(y_)) = w(ixO^S,mom(r_))*sin_phi(ixO^S)*sin_theta(ixO^S)
     ! else where(sin_phi(ixO^S)<=0.0_dp.and.self%patch(ixO^S))
     !   w(ixO^S,mom(y_)) = zero
     !  end where
      where(cos_phi(ixOmin1:ixOmax1)>0.0_dp.and.self%patch(ixOmin1:ixOmax1))
       w(ixOmin1:ixOmax1,mom(z_)) = v_spherical(ixOmin1:ixOmax1,&
          r_)*cos_theta(ixOmin1:ixOmax1)
      else where(cos_phi(ixOmin1:ixOmax1)<=&
         0.0_dp.and.self%patch(ixOmin1:ixOmax1))
       w(ixOmin1:ixOmax1,mom(z_)) = zero
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
  subroutine usr_envelope_clean_memory(self)
    class(star_envelope)    :: self
    if(allocated(self%patch_escape))deallocate(self%patch_escape)
    if(allocated(self%patch))deallocate(self%patch)
  end subroutine usr_envelope_clean_memory
  !--------------------------------------------------------------------
end module mod_obj_star_envelope
