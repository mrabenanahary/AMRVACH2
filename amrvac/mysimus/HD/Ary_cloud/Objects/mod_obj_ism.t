module mod_obj_ism
  use mod_constants
  use mod_global_parameters
  use mod_hd
  use mod_srmhd_parameters, only: mag,lfac_,psi_,xi_
  use mod_obj_dust, only : dust
  use mod_obj_global_parameters
  use mod_obj_mat
  implicit none

    ! ISM features
    type ism_parameters
      character(len=20)    :: unit           !> physical unit at parameter file
      real(dp)             :: density        !> ISM density  (g/cm^3)
      real(dp)             :: number_density !> ISM number density (1/cm^3)
      real(dp)             :: temperature    !> ISM temperature  (K)
      real(dp)             :: pressure       !> ISM pressure
      real(dp)             :: extend(2,3)    !> region in space (cm)
      real(dp)             :: velocity(3)    !> ISM velocity (cm/s)
      logical              :: tracer_on      !> logical to set tracer
      logical              :: dust_on        !> logical to set dust
      real(dp)             :: dust_frac      !> dust fraction
    end type ism_parameters

    type ISM
      logical, allocatable :: patch(:^D&)           !> spatial patch
      logical, allocatable :: escape_patch(:^D&)    !> spatial patch
      character(len=78)    :: subname               !> subroutine name that call it
      type(ism_parameters) :: myconfig              !> ISM configuation parameters
      type (dust)          :: mydust                !> ISM dust
     contains
     !PRIVATE
     PROCEDURE, PASS(self) :: set_default     => usr_ism_set_default
     PROCEDURE, PASS(self) :: set_complet     => usr_ism_set_complet
     PROCEDURE, PASS(self) :: normalize       => usr_ism_normalize
     PROCEDURE, PASS(self) :: set_w           => usr_ism_set_w
     PROCEDURE, PASS(self) :: read_parameters => usr_ism_read_p
     PROCEDURE, PASS(self) :: write_setting   => usr_ism_write_setting
     PROCEDURE, PASS(self) :: alloc_set_patch => usr_ism_alloc_set_patch
     PROCEDURE, PASS(self) :: clean_memory    => usr_ism_clean_memory
    end type



contains


  !--------------------------------------------------------------------
   !> Read the ism parameters  from a parfile
    subroutine usr_ism_read_p(self,ism_config,files)
      class(ism)                         :: self
      character(len=*),intent(in)        :: files(:)
      type(ism_parameters), intent(out)  :: ism_config
      integer                            :: i_file

      namelist /usr_ism_list/ ism_config

      if(mype==0)write(*,*)'Reading usr_ism_list'
      do i_file = 1, size(files)
         open(unitpar, file=trim(files(i_file)), status="old")
         read(unitpar, usr_ism_list)
         close(unitpar)
      end do
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
     write(unit_config,*) 'Density     = ',  self%myconfig%density
     write(unit_config,*) 'Pressure    = ',  self%myconfig%pressure
     write(unit_config,*) 'Temperature = ',  self%myconfig%temperature
     write(unit_config,*) 'Speed       = ',  self%myconfig%velocity
     call self%mydust%write_setting(unit_config)
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
     self%myconfig%unit               = 'code'
     self%myconfig%density            = 0.0_dp
     self%myconfig%number_density     = 0.0_dp
     self%myconfig%temperature        = 0.0_dp
     self%myconfig%pressure           = 0.0_dp
     self%myconfig%extend(1:2,1:ndim) = 0.0_dp!box_limit(1:2,1:ndim)
     self%myconfig%velocity           = 0.0_dp


     self%myconfig%tracer_on          = .false.
     self%myconfig%dust_on            = .false.
     self%myconfig%dust_frac          = 0.0_dp

    call self%mydust%set_default

   end subroutine usr_ism_set_default
  !--------------------------------------------------------------------
  !> subroutine check the parfile setting for ism
  subroutine usr_ism_set_complet(self)
    implicit none
    class(ism)               :: self
    ! .. local ..
    logical                  :: dust_is_frac
    real(dp)                 :: mp,kb
    !-----------------------------------
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
    else
      mp=mp_cgs
      kB=kB_cgs
    end if

    if (dabs( self%myconfig%density)<smalldouble*mp)then
         self%myconfig%density= self%myconfig%number_density*mp
    end if

    if(dabs( self%myconfig%pressure)<smalldouble) then
       self%myconfig%pressure =(2.d0+3.d0*He_abundance)* self%myconfig%number_density*kB* self%myconfig%temperature
    end if

    if( self%myconfig%dust_on)then
      dust_is_frac=.false.
      call  self%mydust%set_complet(dust_is_frac, self%myconfig%dust_frac,&
                              self%myconfig%density, self%myconfig%velocity)
    end if
  end subroutine usr_ism_set_complet
  !--------------------------------------------------------------------
  !> subroutine normalize setting for ISM
   subroutine usr_ism_normalize(self)
    use mod_obj_usr_unit
    implicit none
    class(ism)                         :: self
    type(usr_physical_unit_values)     :: physunit_inuse
    !----------------------------------
    if(trim( self%myconfig%unit)=='code')return
     self%myconfig%density          =  self%myconfig%density/unit_density
     self%myconfig%temperature      =  self%myconfig%temperature/unit_temperature
     self%myconfig%pressure         =  self%myconfig%pressure/unit_pressure
     self%myconfig%velocity         =  self%myconfig%velocity/unit_velocity
     self%myconfig%extend           =  self%myconfig%extend/unit_length
    if( self%myconfig%dust_on)then
      call self%mydust%normalize
      call self%mydust%to_phys
    end if
    if(mype==0)then
     WRITE(*,*)'********************************************************'
     WRITE(*,*)' is in ISM',  self%myconfig%pressure ,  self%myconfig%density, self%myconfig%temperature,sqrt( self%myconfig%pressure /  self%myconfig%density)
     WRITE(*,*)'********************************************************'
    end if
   end subroutine usr_ism_normalize
  !--------------------------------------------------------------------
   !> subroutine setting for ISM
   subroutine usr_ism_set_w(ixI^L,ixO^L,qt,x,w,self)
    implicit none
    integer, intent(in)        :: ixI^L,ixO^L
    real(kind=dp), intent(in)  :: qt
    real(kind=dp)              :: x(ixI^S,1:ndir)
    real(kind=dp)              :: w(ixI^S,1:nw)
    class(ism)                 :: self
    ! .. local..
    integer                    :: idir
    real(kind=dp)              :: fprofile(ixI^S)
    !----------------------------------
    if(.not.allocated(self%patch)) then
     allocate(self%patch(ixG^T))
     self%patch              = .true.
    end if
    where(self%patch(ixO^S))
     w(ixO^S,rho_) =  self%myconfig%density
     w(ixO^S,p_)   =  self%myconfig%pressure
    end where
    Loop_idir : do idir=1,ndir
      where(  self%patch(ixO^S))
       w(ixO^S,mom(idir)) =  self%myconfig%velocity(idir)
      end where
    end do Loop_idir
    cond_tracer_on :if( self%myconfig%tracer_on.and.phys_n_tracer>0 &
                       .and.itr<=phys_n_tracer)then
     where(  self%patch(ixO^S))
      w(ixO^S,tracer(itr)) = w(ixO^S,rho_)
     elsewhere
      w(ixO^S,tracer(itr)) = 0.0_dp
     end where
     itr=itr+1
    end if cond_tracer_on

    cond_dust_on : if( self%myconfig%dust_on)then
       !allocate(  self%mydust%patch(ixG^T))
       !self%mydust%patch=  self%patch
       call self%mydust%set_patch(ixI^L,ixO^L,self%patch)
       self%mydust%myconfig%velocity= self%myconfig%velocity
     fprofile = 1.0_dp
     call   self%mydust%set_w(ixI^L,ixO^L,qt,.false., self%myconfig%dust_frac,fprofile,x,w)
     !deallocate(self%mydust%patch)
    end if cond_dust_on


   end subroutine usr_ism_set_w
   !--------------------------------------------------------------------
   !> Subroutine to clean array memory of associated with cloud object
   subroutine usr_ism_alloc_set_patch(ixI^L,ixO^L,qt,x,escape_patch,self)
     implicit none
     integer, intent(in)           :: ixI^L,ixO^L
     real(kind=dp), intent(in)     :: qt
     real(kind=dp), intent(in)     :: x(ixI^S,1:ndir)
     logical, intent(in),optional  :: escape_patch(ixI^S)
     class(ism)           :: self
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
     self%patch(ixO^S)        = .not.self%escape_patch(ixO^S)
   end subroutine usr_ism_alloc_set_patch
  !---------------------------------------------------------------------

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
