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
    logical              :: normalize_done     !> ism is the normalisation is already done
    integer              :: myindice        !> cloud associated indice
    real(dp)             :: time_cloud_on   !> initial time the cloud is set in the simulation box
    real(dp)             :: density         !> cloud density  (g/cm^3)
    real(dp)             :: number_density  !> cloud number density (1/cm^3)
    real(dp)             :: mass            !> cloud mass (g)
    real(dp)             :: temperature     !> cloud temperature  (K)
    real(dp)             :: pressure        !> cloud pressure
    real(dp)             :: center(3)       !> cloud center position (cm)
    real(dp)             :: extend(3)       !> cloud region in space (cm)
    real(dp)             :: velocity(3)     !> cloud velocity (cm/s)
    real(dp)             :: eject_angle     !> ejection angle max cloud (degre)
    logical              :: tracer_on       !> logical to set tracer
    integer              :: itr            !> ISM tracer indice
    character(len=20)    :: shape           !> cloud shape
    character(len=20)    :: profile         !> could profile
    logical              :: dust_on         !> cloud with dust in is true
    real(dp)             :: dust_frac       !> dust fraction
    character(len=20)    :: dust_profile    !> could dust inside profile
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
      class(cloud)                          :: self
      character(len=*), intent(in)          :: files(:)
      type(cloud_parameters), intent(inout) :: cloud_config
      ! .. local ..
      integer                               :: i_file
      !----------------------------------------------------------------

      namelist /usr_cloud_list/ cloud_config

      if(mype==0)write(*,*)'Reading usr_cloud_list'
      do i_file = 1, size(files)
         open(unitpar, file=trim(files(i_file)), status="old")
         read(unitpar, usr_cloud_list, end=112)
  112    close(unitpar)
      end do
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
    write(unit_config,*) 'Density     = ', self%myconfig%density
    write(unit_config,*) 'Pressure    = ', self%myconfig%pressure
    write(unit_config,*) 'Temperature = ', self%myconfig%temperature
    write(unit_config,*) 'Speed       = ', self%myconfig%velocity
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
  self%myconfig%unit             = 'code'
  self%myconfig%myindice         = 0
  self%myconfig%density          = 0.0_dp
  self%myconfig%number_density   = 0.0_dp
  self%myconfig%mass             = 0.0_dp
  self%myconfig%temperature      = 0.0_dp
  self%myconfig%pressure         = 0.0_dp
  self%myconfig%center(:)        = 0.0_dp!(box_limit(2,:)-box_limit(1,:))/2.0_dp
  self%myconfig%extend(:)        = 0.0_dp!(box_limit(2,:)+box_limit(1,:))/2.0_dp
  self%myconfig%shape            = 'sphere'
  self%myconfig%profile          = 'uniform'
  self%myconfig%eject_angle      =  0.0_dp
  self%myconfig%tracer_on        = .false.
  self%myconfig%itr              = 0
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
   call usr_get_volume(self%myconfig%extend,self%myconfig%shape, cloud_volume)
   PRINT*,' is test volume',cloud_volume,self%myconfig%dust_frac,self%myconfig%mass/cloud_volume,&
    self%myconfig%mass/(cloud_volume*unit_density),cloud_volume/unit_length**3.0_dp
   if (dabs(self%myconfig%mass)>smalldouble)then
    self%myconfig%density        = self%myconfig%mass*(1.0_dp-self%myconfig%dust_frac)/cloud_volume
    self%myconfig%number_density = self%myconfig%density/mp
  else if (dabs(self%myconfig%density)<smalldouble*mp)then
    self%myconfig%density        = self%myconfig%number_density*mp
    if(self%myconfig%dust_frac<1.0_dp)then
       self%myconfig%mass        = self%myconfig%density*cloud_volume/(1.0_dp-self%myconfig%dust_frac)
     else
       self%myconfig%mass        = 0.0_dp
     end if
   else
    if (dabs(self%myconfig%number_density)<smalldouble*mp)then
      self%myconfig%number_density        = self%myconfig%density/mp
    end if
    if(self%myconfig%dust_frac<1.0_dp)then
      self%myconfig%mass        = self%myconfig%density*cloud_volume/(1.0_dp-self%myconfig%dust_frac)
    else
      self%myconfig%mass        = 0.0_dp
    end if
   end if

   if(dabs(self%myconfig%pressure)<smalldouble) then
    self%myconfig%pressure =(2.d0+3.d0*phys_config%He_abundance)*self%myconfig%number_density*kB*self%myconfig%temperature
   end if


   if(.not.self%myconfig%tracer_on)self%myconfig%itr=0


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
   end if
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
  if(self%myconfig%dust_on)then
    call self%mydust%normalize(physunit_inuse)
    call self%mydust%to_phys()
  end if
  self%myconfig%normalize_done=.true.
 end subroutine usr_cloud_normalize
!--------------------------------------------------------------------

 !--------------------------------------------------------------------
 !> Subroutine to clean array memory of associated with cloud object
 subroutine usr_cloud_alloc_set_patch(ixI^L,ixO^L,qt,x,escape_patch,self)
   implicit none
   integer, intent(in)           :: ixI^L,ixO^L
   real(kind=dp), intent(in)     :: qt
   real(kind=dp), intent(in)     :: x(ixI^S,1:ndir)
   logical, intent(in),optional  :: escape_patch(ixI^S)
   class(cloud)                  :: self
   !---------------------------------------------------------
   if(allocated(self%patch))deallocate(self%patch)
   allocate(self%patch(ixI^S))
   if(allocated(self%escape_patch))deallocate(self%escape_patch)
   allocate(self%escape_patch(ixI^S))
   if(present(escape_patch))then
     self%escape_patch(ixO^S) = escape_patch(ixO^S)
   else
     self%escape_patch(ixO^S) =.false.
   end if


   call self%the_patch(ixI^L,ixO^L,x)
   self%patch(ixO^S)        = self%patch(ixO^S).and. (.not.self%escape_patch(ixO^S))
 end subroutine usr_cloud_alloc_set_patch
!---------------------------------------------------------------------
!--------------------------------------------------------------------
 !> subroutine patch for the cloud
 subroutine usr_cloud_patch(ixI^L,ixO^L,x,self)
  implicit none
  integer, intent(in)  :: ixI^L,ixO^L
  real(kind=dp)        :: x(ixI^S,1:ndir)
  class(cloud)         :: self
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
  case default
     write(*,*)'this cloud shape ',trim(self%myconfig%shape),' is not implimented'
     call mpistop('This cloud shape is not implimented in usr_cloud_patch at mod_usr.t')
  end select
 end subroutine usr_cloud_patch
!--------------------------------------------------------------------
 !> subroutine setting for cloud
 subroutine usr_cloud_set_w(ixI^L,ixO^L,qt,x,w,self)
  implicit none
  integer, intent(in)        :: ixI^L,ixO^L
  real(kind=dp), intent(in)  :: qt
  real(kind=dp)              :: x(ixI^S,1:ndir)
  real(kind=dp)              :: w(ixI^S,1:nw)
  class(cloud)               :: self
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
    w(ixO^S,phys_ind%rho_) = self%myconfig%density
    w(ixO^S,phys_ind%pressure_)   = self%myconfig%pressure
   end where
   Loop_idir : do idir=1,ndir
       where(self%patch(ixO^S))
         w(ixO^S,phys_ind%mom(idir)) = self%myconfig%velocity(idir)
       end where
   end do Loop_idir
   if(self%myconfig%profile/='none') then
    call usr_mat_profile(ixI^L,ixO^L,typeaxial,self%myconfig%profile,self%myconfig%center,self%myconfig%extend,&
                         x,fprofile)
    where(self%patch(ixO^S))
      w(ixO^S,phys_ind%rho_)=w(ixO^S,phys_ind%rho_)*fprofile(ixO^S)
      !w(ixO^S,phys_ind%pressure_)=w(ixO^S,phys_ind%pressure_)*fprofile(ixO^S)
    end where
   end if



   cond_tracer_on :if(self%myconfig%tracer_on.and.phys_config%n_tracer>0&
                     .and.self%myconfig%itr<=phys_config%n_tracer)then
    where(self%patch(ixO^S))
     w(ixO^S,phys_ind%tracer(self%myconfig%itr)) = w(ixO^S,phys_ind%rho_)
    elsewhere
     w(ixO^S,phys_ind%tracer(self%myconfig%itr)) = 0.0_dp
    end where
    itr=itr+1
   end if cond_tracer_on

   cond_dust_on : if(self%myconfig%dust_on)then
      if(self%myconfig%dust_profile/='none') then
        call usr_mat_profile(ixI^L,ixO^L,typeaxial,self%myconfig%dust_profile,self%myconfig%center,self%myconfig%extend,&
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
!> Subroutine to clean array memory of associated with cloud object
subroutine usr_cloud_clean_memory(self)
  class(cloud)    :: self
  if(allocated(self%patch))deallocate(self%patch)
  if(self%myconfig%dust_on)call self%mydust%clean_memory()
end subroutine usr_cloud_clean_memory
end module mod_obj_cloud
