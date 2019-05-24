module mod_obj_ism
  use mod_constants
  use mod_global_parameters
  use mod_hd, only: hd_dust
  use mod_srmhd_parameters!, only: mag,lfac_,psi_,xi_
  use mod_obj_dust, only : dust
  use mod_obj_global_parameters
  use mod_obj_mat
  use mod_obj_usr_unit
  implicit none

    ! ISM features
    type ism_parameters
      character(len=20)    :: unit !> physical unit at parameter file
      logical              :: normalize_done !> ism is the normalisation is already done
      real(dp)             :: density        !> ISM density  (g/cm3)
      real(dp)             :: number_density !> ISM number density (1/cm3)
      real(dp)             :: temperature    !> ISM temperature  (K)
      real(dp)             :: pressure       !> ISM pressure
      real(dp)             :: extend(2,3)    !> region in space (cm)
      integer              :: myindice !> ism indices associated with ism in use
      real(dp)             :: kappa !> ISM index power in pressure
      real(dp)             :: z_c !> ISM typical value of height  (cm)
      character(len=30)    :: profile_pressure     !> ism profile pressure
      logical              :: profile_pressure_on  !> ISM pressure profile on
      logical              :: profile_force_on     !> ISM force profile one
      real(dp)             :: velocity(3)    !> ISM velocity (cm/s)
      real(dp)             :: magnetic(3)    !> ISM magnetic field (gauss)
      real(dp)             :: c_sound        !> ISM sound speed
      real(dp)             :: xisigma        !> ISM magentisation
      logical              :: tracer_on      !> logical to set tracer
      integer              :: itr            !> ISM tracer indice
      real(dp)             :: coef_relax     !> ISM relaxation coefficient
      character(len=30)    :: boundary_cond(3,2)!> ism boundary condition
      logical              :: dust_on        !> logical to set dust
      real(dp)             :: dust_frac      !> dust fraction
    end type ism_parameters

    type ISM
      logical, allocatable            :: patch(:)           !> spatial patch
      logical, allocatable            :: escape_patch(:)    !> spatial patch
      character(len=78)               :: subname !> subroutine name that call it
      type(ism_parameters)            :: myconfig !> ISM configuation parameters
      type(usrphysical_unit), pointer :: myphysunit !> ISM physic unity in use
      type (dust)                     :: mydust                !> ISM dust
     contains
     !PRIVATE
     PROCEDURE, PASS(self) :: set_default        => usr_ism_set_default
     PROCEDURE, PASS(self) :: set_complet        => usr_ism_set_complet
     PROCEDURE, PASS(self) :: normalize          => usr_ism_normalize
     PROCEDURE, PASS(self) :: set_w              => usr_ism_set_w
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
      class(ism)                         :: self
      character(len=*),intent(in)        :: files(:)
      type(ism_parameters), intent(out)  :: ism_config
      ! .. local ..
      integer                            :: i_file,i_error_read

      namelist /usr_ism_list/ ism_config
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
         if(i_error_read>0)then
          write(*,*&
             )'At user side in mod_obj_sn_remnant: error at reading parfile'
          write(*,*)'Check input.  Something was wrong, it will stop'
          call mpistop('It stops at reading the parfile')
         elseif(i_error_read<0)then
          write(*,*&
             )'At user side in mod_obj_sn_remnant: error at reading parfile'
          write(*,*)'Reach the end of the file it will stop'
          call mpistop('It stops at reading the parfile')
         else
           write(*,*) 'En of Reading usr_supernovae_remnant_list'
         end if
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
     write(unit_config,*) 'Density     = ',  self%myconfig%density,&
         '  code unit'
     write(unit_config,*) 'Pressure    = ',  self%myconfig%pressure,&
         '  code unit'
     write(unit_config,*) 'Temperature = ',  self%myconfig%temperature,&
         '  code unit'
     write(unit_config,*) 'Speed       = ',  self%myconfig%velocity,&
         '  code unit'
     write(unit_config,*) 'Sound speed =' ,  self%myconfig%c_sound,&
         '  code unit'
     write(unit_config,*)'      ****** Physical Unit *******   '
     write(unit_config,*) 'Density     = ',&
          self%myconfig%density*self%myphysunit%myconfig%density,'  ',&
        self%myphysunit%myunit%density
     write(unit_config,*) 'Pressure    = ',&
          self%myconfig%pressure*self%myphysunit%myconfig%pressure,'  ',&
        self%myphysunit%myunit%pressure
     write(unit_config,*) 'Temperature = ',&
          self%myconfig%temperature*self%myphysunit%myconfig%temperature,'  ',&
        self%myphysunit%myunit%temperature
     write(unit_config,*) 'Speed       = ',&
          self%myconfig%velocity*self%myphysunit%myconfig%velocity,'  ',&
        self%myphysunit%myunit%velocity
     write(unit_config,*) 'Sound speed =' ,&
          self%myconfig%c_sound*self%myphysunit%myconfig%velocity,'  ',&
        self%myphysunit%myunit%velocity
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
     self%myconfig%unit               = 'code'
     self%myconfig%myindice           = 0
     self%myconfig%density            = 0.0_dp
     self%myconfig%number_density     = 0.0_dp
     self%myconfig%temperature        = 0.0_dp
     self%myconfig%pressure           = 0.0_dp
     self%myconfig%extend(1:2,1:ndim) = 0.0_dp!box_limit(1:2,1:ndim)
     self%myconfig%velocity           = 0.0_dp
     self%myconfig%magnetic           = 0.0_dp
     self%myconfig%xisigma            = 0.0_dp


     self%myconfig%kappa                 = 0.0_dp
     self%myconfig%z_c                   = 0.0_dp
     self%myconfig%profile_pressure      = 'none'
     self%myconfig%profile_pressure_on   = .false.
     self%myconfig%profile_force_on      = .false.
     self%myconfig%coef_relax            = 0.0_dp
     self%myconfig%boundary_cond         = 'fix'
     self%myconfig%c_sound               = 0.0_dp
     self%myconfig%tracer_on             = .false.
     self%myconfig%itr                   = 0


     self%myconfig%dust_on               = .false.
     self%myconfig%dust_frac             = 0.0_dp

     self%myconfig%normalize_done        = .false.
     if(self%myconfig%dust_on)then
      call self%mydust%set_default
     end if
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
       self%myconfig%pressure =(2.d0+3.d0*He_abundance)* &
          self%myconfig%number_density*kB* self%myconfig%temperature
    end if

    self%myconfig%c_sound = sqrt&
       (srmhd_gamma*self%myconfig%pressure/self%myconfig%density)


    select case(self%myconfig%profile_pressure)
    case('none')
     self%myconfig%profile_pressure_on = .false.
    case default
     if(.not.(dabs(self%myconfig%kappa)>smalldouble.and.&
        self%myconfig%z_c>smalldouble))then
       self%myconfig%profile_pressure_on = .false.
       self%myconfig%profile_force_on    = .false.
     end if
    end select


    if(.not.self%myconfig%tracer_on)self%myconfig%itr=0


    if( self%myconfig%dust_on)then
      dust_is_frac=.false.
      call  self%mydust%set_complet(dust_is_frac, self%myconfig%dust_frac,&
         self%myconfig%density, self%myconfig%velocity)
    end if

    if(self%myconfig%tracer_on)then
       prim_wnames(self%myconfig%itr+(phys_ind%tracer(1)-1)) = 'tracer_ism'
       cons_wnames(self%myconfig%itr+(phys_ind%tracer(1)-1)) = 'tracer_ism'
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
    if(trim(self%myconfig%unit)=='code'.or.self%myconfig%normalize_done)then
       if(self%myconfig%normalize_done)then
        write(*,*) 'WARNING: Second call for ISM normalisation',&
            'no new normalisation will be done'
       end if
       return
    end if

     self%myconfig%density          =  self%myconfig%density       &
        /physunit_inuse%myconfig%density
     self%myconfig%number_density   = &
         self%myconfig%number_density/physunit_inuse%myconfig%number_density
     self%myconfig%temperature      =  self%myconfig%temperature   &
        /physunit_inuse%myconfig%temperature
     self%myconfig%pressure         =  self%myconfig%pressure      &
        /physunit_inuse%myconfig%pressure
     self%myconfig%velocity         =  self%myconfig%velocity      &
        /physunit_inuse%myconfig%velocity
     self%myconfig%extend           =  self%myconfig%extend        &
        /physunit_inuse%myconfig%length
     self%myconfig%c_sound          =  self%myconfig%c_sound       &
        /physunit_inuse%myconfig%velocity


     self%myconfig%z_c              =  self%myconfig%z_c           &
        /physunit_inuse%myconfig%length
    if( self%myconfig%dust_on)then
      call self%mydust%normalize(physunit_inuse)
      call self%mydust%to_phys
    end if
    self%myconfig%normalize_done =.true.
   end subroutine usr_ism_normalize
  !--------------------------------------------------------------------
   !> subroutine setting for ISM
   subroutine usr_ism_set_w(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,x,w,&
      isboundary_iB,self)
    implicit none
    integer, intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1
    real(kind=dp), intent(in)     :: qt
    real(kind=dp), intent(in)     :: x(ixImin1:ixImax1,1:ndim)
    real(kind=dp)                 :: w(ixImin1:ixImax1,1:nw)
    integer,             optional :: isboundary_iB
    class(ism)                 :: self
    ! .. local..
    integer                    :: idir,IB,idims
    integer                    :: ixBmin1,ixBmax1
    real(kind=dp)              :: fprofile(ixImin1:ixImax1)
    logical                    :: isboundary
    character(len=30)          :: myboundary_cond
    !----------------------------------

    if(.not.allocated(self%patch)) then
     allocate(self%patch(ixGlo1:ixGhi1))
     self%patch              = .true.
    end if
    if(present(isboundary_iB))then
      
        idims=1
        if(isboundary_iB(1)==idims.and.isboundary_iB(2)==1) then
         if(all(x(ixOmin1:ixOmax1,idims)<=xprobmin1))then
           isboundary=.true.
           IB = 2*(idims-1)+1
         end if
        end if
        if(isboundary_iB(1)==idims.and.isboundary_iB(2)==2) then
         if(all(x(ixOmin1:ixOmax1,idims)>=xprobmax1))then
           isboundary=.true.
           IB = 2*idims
         end if
        end if
      
      myboundary_cond = self%myconfig%boundary_cond(isboundary_iB(1),&
         isboundary_iB(2))
    else
      myboundary_cond = 'fix'
    end if

    boundary_cond : if(.not.isboundary.or.myboundary_cond=='fix') then


       where(self%patch(ixOmin1:ixOmax1))
        w(ixOmin1:ixOmax1,phys_ind%rho_) =  self%myconfig%density
        w(ixOmin1:ixOmax1,phys_ind%pressure_)   =  self%myconfig%pressure
       end where
       Loop_idir : do idir=1,ndir
         where(  self%patch(ixOmin1:ixOmax1))
          w(ixOmin1:ixOmax1,phys_ind%mom(idir)) =  &
             self%myconfig%velocity(idir)
         end where
       end do Loop_idir




       if(self%myconfig%profile_pressure_on) then
           call self%get_p_profile(ixImin1,ixImax1,ixOmin1,ixOmax1,x,w)
       end if





       cond_tracer_on :if( self%myconfig%tracer_on.and.phys_config%n_tracer>0 &
          .and.self%myconfig%itr<=phys_config%n_tracer)then

        where(  self%patch(ixOmin1:ixOmax1))
         w(ixOmin1:ixOmax1,phys_ind%tracer(self%myconfig%itr)) = 1.0d2 !w(ixOmin1:ixOmax1,phys_ind%rho_)
        elsewhere
         w(ixOmin1:ixOmax1,phys_ind%tracer(self%myconfig%itr)) = 0.0_dp
        end where

       end if cond_tracer_on





       cond_dust_on : if( self%myconfig%dust_on)then
          call self%mydust%set_patch(ixImin1,ixImax1,ixOmin1,ixOmax1,&
             self%patch)
          self%mydust%myconfig%velocity= self%myconfig%velocity
        fprofile = 1.0_dp
        call   self%mydust%set_w(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,.false.,&
            self%myconfig%dust_frac,fprofile,x,w)
       end if cond_dust_on
    else boundary_cond

     select case(iB)
      case(3)
      ixBmin1=ixOmin1;ixBmax1=ixOmax1;
    !  ixBmin2 = ixOmax2+1
    !  ixBmax2 = ixOmax2+nghostcells

     end select

     ! w(ix^D,phys_ind%rho_) =
    end if boundary_cond

   end subroutine usr_ism_set_w

   !--------------------------------------------------------------------

   subroutine usr_ism_get_p_profile(ixImin1,ixImax1, ixOmin1,ixOmax1,x,w,self)
    implicit none
    integer, intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1
    real(kind=dp), intent(in)     :: x(ixImin1:ixImax1,1:ndim)
    real(kind=dp), intent(inout)  :: w(ixImin1:ixImax1,1:nw)
    class(ism)                    :: self
    ! ..local ..
    real(kind=dp)                  :: p_profile(ixImin1:ixImax1)
    !----------------------------------------------------

    select case(trim(self%myconfig%profile_pressure))
    case('king')
     if(dabs(self%myconfig%kappa)>smalldouble.and.&
        self%myconfig%z_c>smalldouble)then
      p_profile(ixOmin1:ixOmax1) = 1.0_dp/(1.0_dp+(x(ixOmin1:ixOmax1,&
         z_)/self%myconfig%z_c)**(-self%myconfig%kappa) )
     else
      p_profile(ixOmin1:ixOmax1) = 1.0_dp
     end if
    case('komissarov')
     if(dabs(self%myconfig%kappa)>smalldouble.and.&
        self%myconfig%z_c>smalldouble)then
      p_profile(ixOmin1:ixOmax1) = (x(ixOmin1:ixOmax1,&
         z_)/self%myconfig%z_c)**(-self%myconfig%kappa)
     else
      p_profile(ixOmin1:ixOmax1) = 1.0_dp
     end if
    case default
     p_profile(ixOmin1:ixOmax1) = 1.0_dp
    end select

    w(ixOmin1:ixOmax1,phys_ind%pressure_) = w(ixOmin1:ixOmax1,&
       phys_ind%pressure_) * p_profile(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,phys_ind%rho_) = w(ixOmin1:ixOmax1,&
       phys_ind%rho_) * p_profile(ixOmin1:ixOmax1)**(1.0_dp/srmhd_gamma)
   end subroutine usr_ism_get_p_profile


   !--------------------------------------------------------------------

   subroutine usr_ism_get_pforce_profile(ixImin1,ixImax1,ixOmin1,ixOmax1,x,w,&
      f_profile,self)
    implicit none
    integer, intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1
    real(kind=dp), intent(in)     :: x(ixImin1:ixImax1,1:ndim)
    real(kind=dp), intent(in)     :: w(ixImin1:ixImax1,1:nw)
    class(ism)                    :: self
    real(kind=dp), intent(inout)  :: f_profile(ixImin1:ixImax1,1:ndim)
    ! ..local ..
    !----------------------------------------------------

    select case(trim(self%myconfig%profile_pressure))
    case('komissarov')
     cond_force1: if(dabs(self%myconfig%kappa)>smalldouble.and.&
        self%myconfig%z_c>smalldouble)then
      where(self%patch(ixOmin1:ixOmax1))
       f_profile(ixOmin1:ixOmax1,z_) =- self%myconfig%kappa/self%myconfig%z_c*(&
          x(ixOmin1:ixOmax1,z_)/self%myconfig%z_c)**(-(self%myconfig%kappa+1))
      end  where
     else cond_force1
      where(self%patch(ixOmin1:ixOmax1))
       f_profile(ixOmin1:ixOmax1,z_) = 0.0_dp
      end where
     end if cond_force1
    case default
     where(self%patch(ixOmin1:ixOmax1))
      f_profile(ixOmin1:ixOmax1,z_) = 0.0_dp
     end where
    end select


   end subroutine usr_ism_get_pforce_profile

   !--------------------------------------------------------------------
   subroutine usr_ism_add_source(ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,x,&
      qdt,qtC,wCT,qt,w,use_tracer,escape_patch,self)
     implicit none
     integer, intent(in)                     :: ixImin1,ixImax1,ixOmin1,&
        ixOmax1,iwmin,iwmax
     real(kind=dp), intent(in)               :: qdt,qtC,qt
     real(kind=dp), intent(in)               :: x(ixImin1:ixImax1,1:ndim)
     real(kind=dp), intent(in)               :: wCT(ixImin1:ixImax1,1:nw)
     real(kind=dp), intent(inout)            :: w(ixImin1:ixImax1,1:nw)
     logical, intent(in), optional           :: use_tracer
     logical, intent(in),optional            :: escape_patch(ixImin1:ixImax1)

     class(ism)                              :: self
     ! .. local ..
     real(kind=dp)                           :: f_profile(ixImin1:ixImax1,&
        1:ndim)
     real(kind=dp)                           :: w_init(ixImin1:ixImax1,1:nw)
     integer                                 :: idir
     !---------------------------------------------------------



     call self%alloc_set_patch(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,x,use_tracer,&
        wCT,escape_patch)




     call self%get_pforce_profile(ixImin1,ixImax1,ixOmin1,ixOmax1,x,wCT,&
        f_profile)
     where(self%patch(ixOmin1:ixOmax1))
      w(ixOmin1:ixOmax1,phys_ind%mom(z_)) = qdt* f_profile(ixOmin1:ixOmax1,z_)
     end where
     !if(energy .and. .not.block%e_is_internal) then
      where(self%patch(ixOmin1:ixOmax1))
        w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,&
           e_) + qdt * f_profile(ixOmin1:ixOmax1,z_) * wCT(ixOmin1:ixOmax1,&
           phys_ind%mom(z_))/wCT(ixOmin1:ixOmax1,phys_ind%rho_)
      end where
     !end if


     if(self%myconfig%coef_relax>0.0_dp)then
       call self%set_w(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,x,w_init)
       call phys_to_primitive(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
       where(self%patch(ixOmin1:ixOmax1))
         w(ixOmin1:ixOmax1,phys_ind%rho_) = w(ixOmin1:ixOmax1,&
            phys_ind%rho_)*(1.0_dp-self%myconfig%coef_relax) &
            +w_init(ixOmin1:ixOmax1,phys_ind%rho_)*self%myconfig%coef_relax
         w(ixOmin1:ixOmax1,phys_ind%pressure_) = w(ixOmin1:ixOmax1,&
            phys_ind%pressure_)*(1.0_dp-self%myconfig%coef_relax) &
            +w_init(ixOmin1:ixOmax1,phys_ind%pressure_)&
            *self%myconfig%coef_relax
       end where
       loop_idir :  do idir = 1,ndir
        where(self%patch(ixOmin1:ixOmax1))
          w(ixOmin1:ixOmax1,phys_ind%mom(idir)) = w(ixOmin1:ixOmax1,&
             phys_ind%mom(idir))*(1.0_dp-self%myconfig%coef_relax) &
             +w_init(ixOmin1:ixOmax1,phys_ind%mom(&
             idir))*self%myconfig%coef_relax
        end where
       end do loop_idir

       call phys_to_conserved(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
     end if
    !if(any(self%patch(ixO^S)))print*,' test ism force',maxval(dabs(w(ixO^S,phys_ind%mom(z_))),mask=self%patch(ixO^S))
   end subroutine usr_ism_add_source
   !--------------------------------------------------------------------
   !> Subroutine to clean array memory of associated with cloud object
   subroutine usr_ism_alloc_set_patch(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,x,&
      use_tracer,w,escape_patch,self)
     implicit none
     integer, intent(in)                     :: ixImin1,ixImax1,ixOmin1,&
        ixOmax1
     real(kind=dp), intent(in)               :: qt
     real(kind=dp), intent(in)               :: x(ixImin1:ixImax1,1:ndim)
     real(kind=dp), intent(in), optional     :: w(ixImin1:ixImax1,1:nw)
     logical, intent(in), optional           :: use_tracer
     logical, intent(in),optional            :: escape_patch(ixImin1:ixImax1)
     class(ism)                              :: self
     !---------------------------------------------------------
     cond_tracer : if(.not.present(use_tracer).and. .not.present(w)) then
       if(allocated(self%patch))deallocate(self%patch)
       allocate(self%patch(ixImin1:ixImax1))
       self%patch(ixOmin1:ixOmax1) =.true.

     else cond_tracer

       cond_ismtracer_on : if(self%myconfig%tracer_on)then
        if(allocated(self%patch))deallocate(self%patch)
        allocate(self%patch(ixImin1:ixImax1))
        where(w(ixOmin1:ixOmax1,phys_ind%tracer(&
           self%myconfig%itr))>small_density)
          self%patch(ixOmin1:ixOmax1)=.true.
        else where
          self%patch(ixOmin1:ixOmax1)=.false.
        end where
       end if cond_ismtracer_on

     end if cond_tracer

     if(allocated(self%escape_patch))deallocate(self%escape_patch)
     allocate(self%escape_patch(ixImin1:ixImax1))
     if(present(escape_patch))then
      self%escape_patch(ixOmin1:ixOmax1) = escape_patch(ixOmin1:ixOmax1)
     else
      self%escape_patch(ixOmin1:ixOmax1) =.false.
     end if

     where(self%patch(ixOmin1:ixOmax1))self%patch(ixOmin1:ixOmax1)         = &
        .not.self%escape_patch(ixOmin1:ixOmax1)
   end subroutine usr_ism_alloc_set_patch
  !---------------------------------------------------------------------
   subroutine  usr_ism_get_patch_escape(ixImin1,ixImax1,ixOmin1,ixOmax1,&
      need_dealloc,escape_patch,self)
     implicit none
     integer, intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1
     logical, intent(in)           :: need_dealloc
     logical, intent(in)           :: escape_patch(ixImin1:ixImax1)
     class(ism)                    :: self
     !----------------------------------------------------------
     if(allocated(self%escape_patch))deallocate(self%escape_patch)
      allocate(self%escape_patch(ixImin1:ixImax1))
      self%escape_patch(ixOmin1:ixOmax1) = .false.

     self%escape_patch(ixOmin1:ixOmax1)=self%escape_patch(ixOmin1:ixOmax1).or.&
        escape_patch(ixOmin1:ixOmax1)
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
