module mod_obj_cloud
  use mod_constants
  use mod_global_parameters
!  use mod_srmhd_parameters!, only: mag,lfac_,psi_,xi_
  use mod_physics
  use mod_obj_dust, only : dust
  use mod_obj_global_parameters
  use mod_obj_mat
  use mod_obj_usr_unit
  implicit none
  type cloud_parameters
    character(len=20)    :: unit            !> physical unit at parameter file
    character(len=78)    :: obj_name        !> cloud obj name
    logical              :: normalize_done  !> cloud is the normalisation is already done
    integer              :: myindice        !> cloud associated indice
    real(dp)             :: time_cloud_on   !> initial time the cloud is set in the simulation box
    real(dp)             :: density         !> cloud density  (g/cm^3)
    real(dp)             :: number_density  !> cloud number density (1/cm^3)
    real(dp)             :: mass            !> cloud mass (g)
    real(dp)             :: temperature     !> cloud temperature  (K)
    real(dp)             :: pressure        !> cloud pressure
    real(dp)             :: c_sound         !> cloud sound speed
    real(dp)             :: center(3)       !> cloud center position (cm)
    real(dp)             :: extend(3)       !> cloud region in space (cm)
    real(dp)             :: velocity(3)     !> cloud velocity (cm/s)
    real(dp)             :: eject_angle     !> ejection angle max cloud (degre)
    logical              :: tracer_on       !> cloud logical to set tracer
    integer              :: itr             !> cloud tracer indice
    real(dp)             :: tracer_init_density    !> cloud tracer initial density
    real(dp)             :: tracer_small_density   !> cloudtracer small density cut
    character(len=20)    :: shape           !> cloud shape
    character(len=20)    :: profile         !> could profile
    logical              :: reset_on        !> cloud : reset
    real(dp)             :: reset_time      !> cloud : reset characteristic time
    real(dp)             :: reset_speed     !> cloud : reset characteristic speed
    real(dp)             :: reset_Mach      !> cloud : reset characteristic Mach number
    real(dp)             :: reset_coef      !> cloud : reset coefficient
    logical              :: dust_on         !> cloud with dust in is true
    real(dp)             :: dust_frac       !> cloud: dust fraction
    character(len=20)    :: dust_profile    !> could: dust inside profile
  end type cloud_parameters
! cloud features
  type cloud
    type(cloud_parameters)            :: myconfig           !> cloud configuration parameters
    type (dust)                       :: mydust             !> cloud dust
    type(usrphysical_unit), pointer   :: myphysunit         !> cloud physics unit in use
    logical, allocatable              :: patch(:^D&)        !> spatial patch
    logical, allocatable              :: escape_patch(:^D&) !> spatial patch
    character(len=78)                 :: subname            !> subroutine name that call it
    contains
     !PRIVATE
     PROCEDURE, PASS(self) :: set_default     => usr_cloud_set_default
     PROCEDURE, PASS(self) :: set_complet     => usr_cloud_set_complet
     PROCEDURE, PASS(self) :: normalize       => usr_cloud_normalize
     PROCEDURE, PASS(self) :: set_w           => usr_cloud_set_w
     PROCEDURE, PASS(self) :: add_source      => usr_cloud_add_source
     PROCEDURE, PASS(self) :: process_grid    => usr_cloud_process_grid
     PROCEDURE, PASS(self) :: read_parameters => usr_cloud_read_p
     PROCEDURE, PASS(self) :: write_setting   => usr_cloud_write_setting
     PROCEDURE, PASS(self) :: alloc_set_patch => usr_cloud_alloc_set_patch
     PROCEDURE, PASS(self) :: the_patch       => usr_cloud_patch
     PROCEDURE, PASS(self) :: clean_memory    => usr_cloud_clean_memory
  end type

contains

  !-------------------------------------------------------------------------
   !> Read the ism parameters  from a parfile
    subroutine usr_cloud_read_p(files,cloud_config,self)
      use mod_obj_mat
      implicit none
      class(cloud)                          :: self
      character(len=*), intent(in)          :: files(:)
      type(cloud_parameters), intent(inout) :: cloud_config
      ! .. local ..
      integer                               :: i_file,i_error_read
      !----------------------------------------------------------------

      namelist /usr_cloud_list/ cloud_config
      namelist /usr_cloud1_list/ cloud_config
      namelist /usr_cloud2_list/ cloud_config
      namelist /usr_cloud3_list/ cloud_config

      if(mype==0)write(*,*)'Reading usr_cloud_list'
      Loop_ifile : do i_file = 1, size(files)
         open(unitpar, file=trim(files(i_file)), status="old")
         select case(self%myconfig%myindice)
         case(1)
           read(unitpar, usr_cloud1_list, iostat=i_error_read)
         case(2)
           read(unitpar, usr_cloud2_list, iostat=i_error_read)
         case(3)
           read(unitpar, usr_cloud3_list, iostat=i_error_read)
         case default
           read(unitpar, usr_cloud_list, iostat=i_error_read)
         end select
         call usr_mat_read_error_message(i_error_read,self%myconfig%myindice,self%myconfig%obj_name)
         close(unitpar)
      end do Loop_ifile

      if(cloud_config%dust_on)then
        self%mydust%myconfig%associated_medium = 'cloud'
        call self%mydust%read_parameters(self%mydust%myconfig,files)
      end if

    end subroutine usr_cloud_read_p

  !------------------------------------------------------------------------
  !> write the cloud setting
  subroutine usr_cloud_write_setting(self,unit_config)
    implicit none
    class(cloud)                        :: self
    integer,intent(in)                  :: unit_config
    ! .. local ..

    !-----------------------------------

    write(unit_config,*)'************************************'
    write(unit_config,*)'************CLOUD setting ************'
    write(unit_config,*)'************************************'
     write(unit_config,*) 'Density     = ',  self%myconfig%density, '  code unit'
     write(unit_config,*) 'Pressure    = ',  self%myconfig%pressure, '  code unit'
     write(unit_config,*) 'Temperature = ',  self%myconfig%temperature, '  code unit'
     write(unit_config,*) 'Speed       = ',  self%myconfig%velocity, '  code unit'
     write(unit_config,*) 'extend      = ',  self%myconfig%extend, '  code unit'

     write(unit_config,*)'      ****** Physical Unit *******   '
     write(unit_config,*) 'Density     = ',  self%myconfig%density*self%myphysunit%myconfig%density,&
                                             '  ',self%myphysunit%myunit%density
     write(unit_config,*) 'Pressure    = ',  self%myconfig%pressure*self%myphysunit%myconfig%pressure,&
                                             '  ',self%myphysunit%myunit%pressure
     write(unit_config,*) 'Temperature = ',  self%myconfig%temperature*self%myphysunit%myconfig%temperature,&
                                             '  ',self%myphysunit%myunit%temperature
     write(unit_config,*) 'Speed       = ',  self%myconfig%velocity*self%myphysunit%myconfig%velocity,&
                                             '  ',self%myphysunit%myunit%velocity
     write(unit_config,*) 'extend      = ',  self%myconfig%extend*self%myphysunit%myconfig%length

    call self%mydust%write_setting(unit_config)
    write(unit_config,*)'************************************'
    write(unit_config,*)'******** END CLOUD setting **********'
    write(unit_config,*)'************************************'
  end    subroutine usr_cloud_write_setting
  !--------------------------------------------------------------------
!> subroutine default setting for cloud
 subroutine usr_cloud_set_default(self)
  implicit none
  class(cloud)          :: self
  !----------------------------------
  self%myconfig%obj_name         = 'cloud'
  self%myconfig%unit             = 'code'
  self%myconfig%myindice         = 0
  self%myconfig%time_cloud_on    = 0.0_dp
  self%myconfig%density          = 0.0_dp
  self%myconfig%number_density   = 0.0_dp
  self%myconfig%mass             = 0.0_dp
  self%myconfig%temperature      = 0.0_dp
  self%myconfig%pressure         = 0.0_dp
  self%myconfig%c_sound          = 0.0_dp
  self%myconfig%center(:)        = 0.0_dp!(box_limit(2,:)-box_limit(1,:))/2.0_dp
  self%myconfig%extend(:)        = 0.0_dp!(box_limit(2,:)+box_limit(1,:))/2.0_dp
  self%myconfig%shape            = 'sphere'
  self%myconfig%profile          = 'uniform'
  self%myconfig%eject_angle      =  0.0_dp
  self%myconfig%tracer_on        = .false.
  self%myconfig%itr              = 0
  self%myconfig%tracer_init_density    = 0.0_dp
  self%myconfig%tracer_small_density   = 0.0_dp
  self%myconfig%reset_on         = .false.
  self%myconfig%reset_time       = bigdouble
  self%myconfig%reset_speed      = 0.0_dp
  self%myconfig%reset_Mach       = 0.0_dp
  self%myconfig%dust_on          = .false.
  self%myconfig%dust_frac        = 0.0_dp
  self%myconfig%dust_profile     = 'uniform'
  self%myconfig%normalize_done   = .false.
  call self%mydust%set_default()
 end subroutine usr_cloud_set_default
 !--------------------------------------------------------------------
 !> subroutine check the parfile setting for cloud
 subroutine usr_cloud_set_complet(self)
   implicit none
   class(cloud)             :: self
   ! .. local ..
   logical                  :: dust_is_frac
   real(dp)                 :: mp,kb,cloud_volume
   !-----------------------------------
   if(SI_unit) then
     mp=mp_SI
     kB=kB_SI
   else
     mp=mp_cgs
     kB=kB_cgs
   end if

  if(.not.self%myconfig%dust_on)then
    self%myconfig%dust_frac = 0.0_dp
  end if

   call usr_get_volume(self%myconfig%extend,self%myconfig%shape, cloud_volume)

   cond_mass_usr: if (dabs(self%myconfig%mass)>smalldouble)then
    self%myconfig%density        = self%myconfig%mass*(1.0_dp-self%myconfig%dust_frac)/cloud_volume
    self%myconfig%number_density = self%myconfig%density/mp
   else cond_mass_usr
    if (dabs(self%myconfig%density)<smalldouble*mp)then
     self%myconfig%density        = self%myconfig%number_density*mp
    else if (dabs(self%myconfig%number_density)<smalldouble)then
     self%myconfig%number_density        = self%myconfig%density/mp
    end if
    if(self%myconfig%dust_frac<1.0_dp.and.self%myconfig%dust_frac>0.0_dp)then
       self%myconfig%mass        = self%myconfig%density*cloud_volume/(1.0_dp-self%myconfig%dust_frac)
    else
       self%myconfig%mass        = self%myconfig%density*cloud_volume
    end if
   end if cond_mass_usr

   if(dabs(self%myconfig%pressure)<smalldouble) then
    self%myconfig%pressure =(2.d0+3.d0*phys_config%He_abundance)*self%myconfig%number_density*kB*self%myconfig%temperature
   end if
   self%myconfig%c_sound     = sqrt(phys_config%gamma*self%myconfig%pressure/self%myconfig%density)



   cond_duston : if(self%myconfig%dust_on)then
    if(self%myconfig%mass>smalldouble) then
      dust_is_frac=.true.
    else
      dust_is_frac=.false.
    end if
     call self%mydust%set_complet(dust_is_frac,self%myconfig%dust_frac,self%myconfig%density,self%myconfig%velocity)
   end if cond_duston

    if(self%myconfig%tracer_on)then
      prim_wnames(self%myconfig%itr+(phys_ind%tracer(1)-1)) = 'tracer_cloud'
      cons_wnames(self%myconfig%itr+(phys_ind%tracer(1)-1)) = 'tracer_cloud'
    else
      self%myconfig%itr=0
    end if

   cond_t_reset : if(self%myconfig%reset_time<bigdouble/2.0) then
     self%myconfig%reset_speed = (cloud_volume)**(1.0_dp/ndim)/self%myconfig%reset_time
     self%myconfig%reset_mach  = self%myconfig%reset_speed/self%myconfig%c_sound
     self%myconfig%reset_on    = .true.
   else if(dabs(self%myconfig%reset_speed)>0.0_dp)then
     self%myconfig%reset_time  = (cloud_volume)**(1.0_dp/ndim)/self%myconfig%reset_speed
     self%myconfig%reset_mach  = self%myconfig%reset_speed/self%myconfig%c_sound
     self%myconfig%reset_on    = .true.
   else if(dabs(self%myconfig%reset_Mach)>0.0_dp)then
     self%myconfig%reset_speed = self%myconfig%reset_mach*self%myconfig%c_sound
     self%myconfig%reset_time  = (cloud_volume)**(1.0_dp/ndim)/self%myconfig%reset_speed
     self%myconfig%reset_on    = .true.
   else if(dabs(self%myconfig%reset_coef)>0.0_dp)then
     self%myconfig%reset_mach  = self%myconfig%reset_coef
     self%myconfig%reset_speed = self%myconfig%reset_mach*self%myconfig%c_sound
     self%myconfig%reset_time  = (cloud_volume)**(1.0_dp/ndim)/self%myconfig%reset_speed
     self%myconfig%reset_on    = .true.
   end if cond_t_reset


 end subroutine usr_cloud_set_complet
!--------------------------------------------------------------------
 subroutine usr_cloud_normalize(self,physunit_inuse)
  use mod_obj_usr_unit
  implicit none
  class(cloud)                                   :: self
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
  self%myconfig%temperature      = self%myconfig%temperature   /physunit_inuse%myconfig%temperature
  self%myconfig%pressure         = self%myconfig%pressure      /physunit_inuse%myconfig%pressure
  self%myconfig%velocity         = self%myconfig%velocity      /physunit_inuse%myconfig%velocity
  self%myconfig%eject_angle      = self%myconfig%eject_angle   *(dpi/180._dp)
  self%myconfig%center           = self%myconfig%center        /physunit_inuse%myconfig%length
  self%myconfig%extend           = self%myconfig%extend        /physunit_inuse%myconfig%length
  self%myconfig%time_cloud_on    = self%myconfig%time_cloud_on /physunit_inuse%myconfig%time
  if(self%myconfig%dust_on)then
    call self%mydust%normalize(physunit_inuse)
    call self%mydust%to_phys()
  end if
  self%myconfig%normalize_done=.true.
 end subroutine usr_cloud_normalize
!--------------------------------------------------------------------

 !--------------------------------------------------------------------
 !> Subroutine to clean array memory of associated with cloud object
 subroutine usr_cloud_alloc_set_patch(ixI^L,ixO^L,qt,x,self,use_tracer,w,escape_patch)
   implicit none
   integer, intent(in)                     :: ixI^L,ixO^L
   real(kind=dp), intent(in)               :: qt
   real(kind=dp), intent(in)               :: x(ixI^S,1:ndim)
   real(kind=dp), intent(in), optional     :: w(ixI^S,1:nw)
   logical, intent(in), optional           :: use_tracer
   logical, intent(in),optional            :: escape_patch(ixI^S)
   class(cloud)                            :: self
   !---------------------------------------------------------

   cond_tracer : if(.not.present(use_tracer).and. .not.present(w)) then
    if(.not.allocated(self%patch))call self%the_patch(ixI^L,ixO^L,x)
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


   if(.not.allocated(self%patch))call self%the_patch(ixI^L,ixO^L,x)

   self%patch(ixO^S)        = self%patch(ixO^S).and. (.not.self%escape_patch(ixO^S))
 end subroutine usr_cloud_alloc_set_patch
!---------------------------------------------------------------------
!--------------------------------------------------------------------
 !> subroutine patch for the cloud
 subroutine usr_cloud_patch(ixI^L,ixO^L,x,self,usr_patch)
  implicit none
  integer, intent(in)              :: ixI^L,ixO^L
  real(kind=dp), intent(in)        :: x(ixI^S,1:ndim)
  class(cloud)                     :: self
  logical, optional                :: usr_patch(ixO^S)
  !----------------------------------

  if(.not.allocated(self%patch))allocate(self%patch(ixG^T))

  self%patch              = .false.
  select case(trim(self%myconfig%shape))
  case('sphere')
    call usr_set_patch_sphere(ixI^L,ixO^L,typeaxial,self%myconfig%center,    &
                              self%myconfig%extend,x,self%patch)

  case('cylinder')
    call usr_set_patch_cylinder(ixI^L,ixO^L,typeaxial,self%myconfig%center,  &
                                self%myconfig%extend,x,self%patch)
  case('cube')
    call usr_set_patch_cube(ixI^L,ixO^L,typeaxial,self%myconfig%center,      &
                            self%myconfig%extend,x,self%patch)
  case('usr')
    if(present(usr_patch))then
     self%patch(ixO^S) = usr_patch(ixO^S)
    else
     write(*,*) 'the subroutine usr_cloud_patch is called from ',self%subname
     write(*,*) 'this cloud shape ','usr',' is not implimented in mod_usr.t '
     write(*,*) ' at subroutine  usr_cloud_patch'
     call mpistop('This cloud shape is not implimented in  mod_obj_cloud.t')
    end if
  case default
     write(*,*)'this cloud shape : ',trim(self%myconfig%shape),' is not implimented'
     call mpistop('This cloud shape is not implimented in usr_cloud_patch at mod_obj_cloud.t')
  end select
 end subroutine usr_cloud_patch
!--------------------------------------------------------------------
 !> subroutine setting for cloud
 subroutine usr_cloud_set_w(ixI^L,ixO^L,qt,x,w,self,usr_density_profile,&
                            usr_pressure_profile,usr_velocity_profile)
  implicit none
  integer, intent(in)        :: ixI^L,ixO^L
  real(kind=dp), intent(in)  :: qt
  real(kind=dp)              :: x(ixI^S,1:ndim)
  real(kind=dp)              :: w(ixI^S,1:nw)
  class(cloud)               :: self
  real(kind=dp), optional    :: usr_density_profile(ixI^S)
  real(kind=dp), optional    :: usr_pressure_profile(ixI^S)
  real(kind=dp), optional    :: usr_velocity_profile(ixI^S,1:ndir)
  ! .. local..
  integer                    :: idir
  logical                    :: dust_is_frac
  real(kind=dp)              :: fprofile(ixI^S)
  !----------------------------------
!  call usr_cloud_patch(ixI^L,ixO^L,x,self)

  if(.not.allocated(self%patch))then
    call usr_cloud_patch(ixI^L,ixO^L,x,self)
  end if


  cond_inside_cloud: if(any(self%patch(ixO^S))) then
   where(self%patch(ixO^S))
    w(ixO^S,phys_ind%rho_)        = self%myconfig%density
    w(ixO^S,phys_ind%pressure_)   = self%myconfig%pressure
   end where
   Loop_idir : do idir=1,ndir
       where(self%patch(ixO^S))
         w(ixO^S,phys_ind%mom(idir)) = self%myconfig%velocity(idir)
       end where
   end do Loop_idir






   cond_cloud_profile : if(self%myconfig%profile/='none') then
    call usr_mat_profile(ixI^L,ixO^L,typeaxial,self%myconfig%profile,&
                         self%myconfig%center,self%myconfig%extend,&
                         x,fprofile)
    where(self%patch(ixO^S))
      w(ixO^S,phys_ind%rho_)=w(ixO^S,phys_ind%rho_)*fprofile(ixO^S)
      !w(ixO^S,phys_ind%pressure_)=w(ixO^S,phys_ind%pressure_)*fprofile(ixO^S)
    end where
   else cond_cloud_profile
      cond_usr_rhoprofile : if(present(usr_density_profile))then
       where(self%patch(ixO^S))
        w(ixO^S,phys_ind%rho_)=w(ixO^S,phys_ind%rho_)*usr_density_profile(ixO^S)
       end where
      end if cond_usr_rhoprofile
      cond_usr_pprofile : if(present(usr_pressure_profile))then
       where(self%patch(ixO^S))
        w(ixO^S,phys_ind%pressure_)=w(ixO^S,phys_ind%pressure_)*usr_pressure_profile(ixO^S)
       end where
      end if cond_usr_pprofile
      cond_usr_vprofile : if(present(usr_velocity_profile))then
        Loop_idir1 : do idir=1,ndir
         where(self%patch(ixO^S))
          w(ixO^S,phys_ind%mom(idir)) =w(ixO^S,phys_ind%mom(idir))*usr_velocity_profile(ixO^S,idir)
         end where
       end do Loop_idir1
     end if cond_usr_vprofile

   end if cond_cloud_profile



   cond_tracer_on :if(self%myconfig%tracer_on.and.phys_config%n_tracer>0&
                      .and.self%myconfig%itr<=phys_config%n_tracer)then
     cond_cloud_on : if(qt< self%myconfig%time_cloud_on)then
      where(.not.self%patch(ixO^S))
       w(ixO^S,phys_ind%tracer(self%myconfig%itr)) = 0.0_dp
      end where
     else cond_cloud_on
      if(self%myconfig%tracer_init_density>0.0_dp) then
       where(self%patch(ixO^S))
        w(ixO^S,phys_ind%tracer(self%myconfig%itr)) = self%myconfig%tracer_init_density
       end where
      else
       where(self%patch(ixO^S))
        w(ixO^S,phys_ind%tracer(self%myconfig%itr)) =  w(ixO^S,phys_ind%rho_)
       end where
      end if
     end if cond_cloud_on
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
     if(self%myconfig%mass> smalldouble) then
       dust_is_frac=.true.
     else
       dust_is_frac=.false.
     end if
     call self%mydust%set_w(ixI^L,ixO^L,qt,dust_is_frac,self%myconfig%dust_frac,fprofile,x,w)
   end if cond_dust_on
  end if cond_inside_cloud



 end subroutine usr_cloud_set_w

  !--------------------------------------------------------------------
   subroutine usr_cloud_add_source(ixI^L,ixO^L,iw^LIM,x,qdt,qtC,wCT,qt,w,use_tracer,&
                                 escape_patch,self)
     implicit none
     integer, intent(in)                     :: ixI^L,ixO^L,iw^LIM
     real(kind=dp), intent(in)               :: qdt,qtC,qt
     real(kind=dp), intent(in)               :: x(ixI^S,1:ndim)
     real(kind=dp), intent(in)               :: wCT(ixI^S,1:nw)
     real(kind=dp), intent(inout)            :: w(ixI^S,1:nw)
     logical, intent(in), optional           :: use_tracer
     logical, intent(in),optional            :: escape_patch(ixI^S)

     class(cloud)                            :: self
     ! .. local ..
     real(kind=dp)                           :: f_profile(ixI^S,1:ndim)
     real(kind=dp)                           :: w_init(ixI^S,1:nw)
     integer                                 :: idir
     !---------------------------------------------------------



     call self%alloc_set_patch(ixI^L,ixO^L,qt,x,use_tracer,wCT,escape_patch)



     if(self%myconfig%reset_coef>0.0_dp)then
       call self%set_w(ixI^L,ixO^L,qt,x,w_init)
       call phys_to_primitive(ixI^L,ixO^L,w,x)
       where(self%patch(ixO^S))
         w(ixO^S,phys_ind%rho_) = w(ixO^S,phys_ind%rho_)*(1.0_dp-self%myconfig%reset_coef) +&
                          w_init(ixO^S,phys_ind%rho_)*self%myconfig%reset_coef
         w(ixO^S,phys_ind%pressure_) = w(ixO^S,phys_ind%pressure_)*(1.0_dp-self%myconfig%reset_coef) +&
                          w_init(ixO^S,phys_ind%pressure_)*self%myconfig%reset_coef
       end where
       loop_idir :  do idir = 1,ndir
        where(self%patch(ixO^S))
          w(ixO^S,phys_ind%mom(idir)) = w(ixO^S,phys_ind%mom(idir))*(1.0_dp-self%myconfig%reset_coef) +&
                          w_init(ixO^S,phys_ind%mom(idir))*self%myconfig%reset_coef
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
     end if
    !if(any(self%patch(ixO^S)))print*,' test ism force',maxval(dabs(w(ixO^S,phys_ind%mom(z_))),mask=self%patch(ixO^S))
  end subroutine usr_cloud_add_source
!--------------------------------------------------------------------
!> Subroutine to process variables in cloud object
 subroutine usr_cloud_process_grid(ixI^L,ixO^L,qt,x,w,self)
  implicit none
  integer, intent(in)        :: ixI^L,ixO^L
  real(kind=dp), intent(in)  :: qt
  real(kind=dp)              :: x(ixI^S,1:ndim)
  real(kind=dp)              :: w(ixI^S,1:nw)
  class(cloud)               :: self
  ! .. local..
  !----------------------------------------------------------
  cond_dust_on : if(self%myconfig%dust_on)then
    call self%mydust%handel_small_val(ixI^L,ixO^L,qt,x,w)
  end if cond_dust_on
end subroutine usr_cloud_process_grid
!--------------------------------------------------------------------
!> Subroutine to clean array memory of associated with cloud object
subroutine usr_cloud_clean_memory(self)
  class(cloud)    :: self
  if(allocated(self%patch))deallocate(self%patch)
  if(self%myconfig%dust_on)call self%mydust%clean_memory()
end subroutine usr_cloud_clean_memory
end module mod_obj_cloud
