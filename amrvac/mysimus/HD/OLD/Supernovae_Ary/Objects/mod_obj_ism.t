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
      character(len=78)    :: obj_name               !> Obj name that call it
      logical              :: normalize_done !> ism is the normalisation is already done
      real(dp)             :: density        !> ISM density  (g/cm^3)
      real(dp)             :: number_density !> ISM number density (1/cm^3)
      real(dp)             :: temperature    !> ISM temperature  (K)
      real(dp)             :: pressure       !> ISM pressure
      real(dp)             :: extend(2,3)    !> region in space (cm)
      integer              :: myindice       !> ism indices associated with ism in use
      real(dp)             :: kappa                !> ISM index power in pressure
      real(dp)             :: z_c                  !> ISM typical value of height  (cm)
      character(len=30)    :: profile_pressure     !> ism profile pressure
      logical              :: profile_pressure_on  !> ISM pressure profile on
      logical              :: profile_force_on     !> ISM force profile one
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
      character(len=30)    :: boundary_cond(3,2)!> ism boundary condition
      logical              :: dust_on        !> logical to set dust
      real(dp)             :: dust_frac      !> dust fraction
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
     PROCEDURE, PASS(self) :: set_default        => usr_ism_set_default
     PROCEDURE, PASS(self) :: set_complet        => usr_ism_set_complet
     PROCEDURE, PASS(self) :: normalize          => usr_ism_normalize
     PROCEDURE, PASS(self) :: set_w              => usr_ism_set_w
     PROCEDURE, PASS(self) :: process_grid       => usr_ism_process_grid
     PROCEDURE, PASS(self) :: read_parameters    => usr_ism_read_p
     PROCEDURE, PASS(self) :: write_setting      => usr_ism_write_setting
     PROCEDURE, PASS(self) :: alloc_set_patch    => usr_ism_alloc_set_patch
     PROCEDURE, PASS(self) :: clean_memory       => usr_ism_clean_memory
     PROCEDURE, PASS(self) :: get_patch_escape   => usr_ism_get_patch_escape

     PROCEDURE, PASS(self) :: get_p_profile      => usr_ism_get_p_profile
     PROCEDURE, PASS(self) :: get_pforce_profile => usr_ism_get_pforce_profile
     PROCEDURE, PASS(self) :: add_source         => usr_ism_add_source

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
     self%myconfig%obj_name           = 'ism'
     self%myconfig%unit               = 'code'
     self%myconfig%myindice           = 0
     self%myconfig%density            = 0.0_dp
     self%myconfig%number_density     = 0.0_dp
     self%myconfig%temperature        = 0.0_dp
     self%myconfig%pressure           = 0.0_dp
     self%myconfig%extend(1:2,1:3)    = 0.0_dp!box_limit(1:2,1:ndim)
     self%myconfig%velocity           = 0.0_dp
     self%myconfig%magnetic           = 0.0_dp
     self%myconfig%xisigma            = 0.0_dp


     self%myconfig%kappa                 = 0.0_dp
     self%myconfig%z_c                   = 0.0_dp
     self%myconfig%profile_pressure      = 'none'
     self%myconfig%profile_pressure_on   = .false.
     self%myconfig%profile_force_on      = .false.
     self%myconfig%reset_coef            = 0.0_dp
     self%myconfig%reset_on              = .false.
     self%myconfig%boundary_cond         = 'fix'
     self%myconfig%c_sound               = 0.0_dp


     self%myconfig%tracer_on             = .false.
     self%myconfig%itr                   = 0
     self%myconfig%tracer_init_density   = 0.0_dp
     self%myconfig%tracer_small_density  = 0.0_dp

     self%myconfig%dust_on               = .false.
     self%myconfig%dust_frac             = 0.0_dp

     self%myconfig%normalize_done        = .false.
     !if(phys_config%dust_on)then
     call self%mydust%set_default

     call self%myboundaries%set_default
     !end if
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
       self%myconfig%pressure =(2.0_dp+3.0_dp*phys_config%He_abundance)* &
                             self%myconfig%number_density*kB* self%myconfig%temperature
    end if

    self%myconfig%c_sound = sqrt(phys_config%gamma*self%myconfig%pressure/self%myconfig%density)


    select case(self%myconfig%profile_pressure)
    case('none')
     self%myconfig%profile_pressure_on = .false.
    case default
     if(.not.(dabs(self%myconfig%kappa)>smalldouble.and.self%myconfig%z_c>smalldouble))then
       self%myconfig%profile_pressure_on = .false.
       self%myconfig%profile_force_on    = .false.
     end if
    end select




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


    self%myboundaries%boundary_type=self%myconfig%boundary_cond
    call self%myboundaries%set_complet
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
    if(trim(self%myconfig%unit)=='code'.or.self%myconfig%normalize_done)then
       if(self%myconfig%normalize_done)then
        write(*,*) 'WARNING: Second call for ISM normalisation', &
                     'no new normalisation will be done'
       end if
       return
    end if

     self%myconfig%density          =  self%myconfig%density       /physunit_inuse%myconfig%density
     self%myconfig%number_density   =  self%myconfig%number_density/physunit_inuse%myconfig%number_density
     self%myconfig%temperature      =  self%myconfig%temperature   /physunit_inuse%myconfig%temperature
     self%myconfig%pressure         =  self%myconfig%pressure      /physunit_inuse%myconfig%pressure
     self%myconfig%velocity         =  self%myconfig%velocity      /physunit_inuse%myconfig%velocity
     self%myconfig%extend           =  self%myconfig%extend        /physunit_inuse%myconfig%length
     self%myconfig%c_sound          =  self%myconfig%c_sound       /physunit_inuse%myconfig%velocity


     self%myconfig%z_c              =  self%myconfig%z_c           /physunit_inuse%myconfig%length
    if( self%myconfig%dust_on)then
      call self%mydust%normalize(physunit_inuse)
      call self%mydust%to_phys
    end if
    self%myconfig%normalize_done =.true.
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
      w(ixO^S,phys_ind%pressure_)   =  self%myconfig%pressure
      end where
      Loop_idir : do idir=1,ndir
       where(  self%patch(ixO^S))
        w(ixO^S,phys_ind%mom(idir)) =  self%myconfig%velocity(idir)
       end where
      end do Loop_idir




      if(self%myconfig%profile_pressure_on) then
         call self%get_p_profile(ixI^L,ixO^L,x,w)
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

   subroutine usr_ism_get_p_profile(ixI^L, ixO^L,x,w,self)
    implicit none
    integer, intent(in)           :: ixI^L,ixO^L
    real(kind=dp), intent(in)     :: x(ixI^S,1:ndim)
    real(kind=dp), intent(inout)  :: w(ixI^S,1:nw)
    class(ism)                    :: self
    ! ..local ..
    real(kind=dp)                  :: p_profile(ixI^S)
    !----------------------------------------------------

    select case(trim(self%myconfig%profile_pressure))
    case('king')
     if(dabs(self%myconfig%kappa)>smalldouble.and.self%myconfig%z_c>smalldouble)then
      p_profile(ixO^S) = 1.0_dp/(1.0_dp+(x(ixO^S,z_)/self%myconfig%z_c)**(-self%myconfig%kappa) )
     else
      p_profile(ixO^S) = 1.0_dp
     end if
    case('komissarov')
     if(dabs(self%myconfig%kappa)>smalldouble.and.self%myconfig%z_c>smalldouble)then
      p_profile(ixO^S) = (x(ixO^S,z_)/self%myconfig%z_c)**(-self%myconfig%kappa)
     else
      p_profile(ixO^S) = 1.0_dp
     end if
    case default
     p_profile(ixO^S) = 1.0_dp
    end select

    w(ixO^S,phys_ind%pressure_) = w(ixO^S,phys_ind%pressure_) * p_profile(ixO^S)
    w(ixO^S,phys_ind%rho_) = w(ixO^S,phys_ind%rho_) * p_profile(ixO^S)**(1.0_dp/phys_config%gamma)
   end subroutine usr_ism_get_p_profile


   !--------------------------------------------------------------------

   subroutine usr_ism_get_pforce_profile(ixI^L,ixO^L,qt,x,w,f_profile,self)
    implicit none
    integer, intent(in)           :: ixI^L,ixO^L
    real(kind=dp), intent(in)     :: qt
    real(kind=dp), intent(in)     :: x(ixI^S,1:ndim)
    real(kind=dp), intent(in)     :: w(ixI^S,1:nw)
    class(ism)                    :: self
    real(kind=dp), intent(inout)  :: f_profile(ixI^S,1:ndim)
    ! ..local ..
    !----------------------------------------------------

    select case(trim(self%myconfig%profile_pressure))
    case('komissarov')
     cond_force1: if(dabs(self%myconfig%kappa)>smalldouble.and.self%myconfig%z_c>smalldouble)then
      where(self%patch(ixO^S))
       f_profile(ixO^S,z_) =- self%myconfig%kappa/self%myconfig%z_c*&
         (x(ixO^S,z_)/self%myconfig%z_c)**(-(self%myconfig%kappa+1))
      end  where
     else cond_force1
      where(self%patch(ixO^S))
       f_profile(ixO^S,z_) = 0.0_dp
      end where
     end if cond_force1
    case default
     where(self%patch(ixO^S))
      f_profile(ixO^S,z_) = 0.0_dp
     end where
    end select



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
     real(kind=dp)                           :: source_filter_loc(ixI^S)
     real(kind=dp)                           :: f_profile(ixI^S,1:ndim)
     real(kind=dp)                           :: w_init(ixI^S,1:nw)
     integer                                 :: idir
     !---------------------------------------------------------



     call self%alloc_set_patch(ixI^L,ixO^L,qt,x,&
                               use_tracer=use_tracer,w=w,escape_patch=escape_patch)




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