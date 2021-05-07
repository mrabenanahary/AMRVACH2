module mod_usr
  use mod_hd
  use mod_physics
  use mod_dust
  use mod_global_parameters
  use mod_obj_global_parameters
  use mod_obj_mat
  use mod_obj_cloud
  use mod_obj_ism
  use mod_obj_star
  implicit none
  save
  real(dp) :: theta, kx, ly, vc
  logical  :: ism_on,cloud_on
  integer  :: cloud_number
  integer  :: itr,phys_n_tracer
  integer, parameter  :: n_dust_max = 20
  real(dp) :: SUM_MASS   = 0.0_dp
  real(dp) :: SUM_VOLUME = 0.0_dp
  ! the dust type
  ! type dust
  !  character(len=20)    :: unit               !> physical unit at parameter file
  !  character(len=20)    :: distrub_func       !> dust distribution function
  !  integer              :: n_species          !> number of dust species
  !  real(dp)             :: sizes(1:n_dust_max)!> dust size
  !  real(dp)             :: min_radius         !> dust smaller radius
  !  real(dp)             :: max_radius         !> dust bigger radius
  !  real(dp)             :: power_a            !> dust power indice  distrbution func
  !  real(dp)             :: grain_density(n_dust_max) !> dust grain density (g/cm)
  !  real(dp)             :: density(n_dust_max)!> dust density
  !  real(dp)             :: number_density(n_dust_max)!> dust density
  !  real(dp)             :: temperature        !> dust temperature  (K)
  !  real(dp)             :: min_limit_rel      !> dust floor relative to gas
  !  real(dp)             :: min_limit_abs      !> dust floor abs (g/cm^3)
  !  real(dp)             :: max_limit_rel      !> dust high relative to gas
  !  real(dp)             :: extend(1:2,1:ndim) !> region in space (cm)
  !  real(dp)             :: velocity(1:3)      !> dust velocity (cm/s)
  !  logical, allocatable :: patch(:^D&)        !> dust is on cell
  !  logical              :: tracer_on          !> logical to set tracer
  !  logical              :: associate_med      !> logical to associate dust to medium
  !  integer              :: idust_first        !> indice of first dust species
  !  integer              :: idust_last         !> indice of first dust species
  !  contains
  !  PRIVATE
  !  PROCEDURE, PASS(self) :: set_default     => usr_dust_set_default
  !  PROCEDURE, PASS(self) :: set_complet     => usr_dust_set_complet
  !  PROCEDURE, PASS(self) :: to_phys         => usr_dust_to_phys
  !  PROCEDURE, PASS(self) :: normalize       => usr_dust_normalize
  !  PROCEDURE, PASS(self) :: set_w           => usr_dust_set_w
  !  PROCEDURE, PASS(self) :: set_w_zero      => usr_dust_set_w_zero
  !  PROCEDURE, PASS(self) :: read_parameters => usr_dust_read_p
  !  PROCEDURE, PASS(self) :: write_setting   => usr_dust_write_setting
  !  PROCEDURE, PASS(self) :: clean_memory    => usr_dust_clean_memory
  ! end type
  !
  ! ! ISM features
  ! type ISM
  !   character(len=20)    :: unit           !> physical unit at parameter file
  !   real(dp)             :: density        !> ISM density  (g/cm^3)
  !   real(dp)             :: number_density !> ISM number density (1/cm^3)
  !   real(dp)             :: temperature    !> ISM temperature  (K)
  !   real(dp)             :: pressure       !> ISM pressure  ()
  !   real(dp)             :: extend(2,ndim) !> region in space (cm)
  !   real(dp)             :: velocity(3)    !> ISM velocity (cm/s)
  !   logical, allocatable :: patch(:^D&)    !> spatial patch
  !   logical              :: tracer_on      !> logical to set tracer
  !   logical              :: dust_on        !> logical to set dust
  !   real(dp)             :: dust_frac      !> dust fraction
  !   type (dust)          :: mydust         !> ISM dust
  !  contains
  !  PRIVATE
  !  PROCEDURE, PASS(self) :: set_default     => usr_ism_set_default
  !  PROCEDURE, PASS(self) :: set_complet     => usr_ism_set_complet
  !  PROCEDURE, PASS(self) :: normalize       => usr_ism_normalize
  !  PROCEDURE, PASS(self) :: set_w           => usr_ism_set_w
  !  PROCEDURE, PASS(self) :: read_parameters => usr_ism_read_p
  !  PROCEDURE, PASS(self) :: write_setting   => usr_ism_write_setting
  !  PROCEDURE, PASS(self) :: clean_memory    => usr_ism_clean_memory
  ! end type
  !
  !
  ! ! cloud features
  ! type cloud
  !   character(len=20)    :: unit            !> physical unit at parameter file
  !   real(dp)             :: density         !> cloud density  (g/cm^3)
  !   real(dp)             :: number_density  !> cloud number density (1/cm^3)
  !   real(dp)             :: mass            !> cloud mass (g)
  !   real(dp)             :: temperature     !> cloud temperature  (K)
  !   real(dp)             :: pressure        !> cloud pressure  ()
  !   real(dp)             :: center(ndim)    !> cloud center position (cm)
  !   real(dp)             :: extend(ndim)    !> cloud region in space (cm)
  !   real(dp)             :: velocity(3)     !> cloud velocity (cm/s)
  !   real(dp)             :: eject_angle     !> ejection angle max cloud (degre)
  !   logical, allocatable :: patch(:^D&)     !> spatial patch
  !   logical              :: tracer_on       !> logical to set tracer
  !   character(len=20)    :: shape           !> cloud shape
  !   character(len=20)    :: profile         !> could profile
  !   logical              :: dust_on         !> cloud with dust in is true
  !   real(dp)             :: dust_frac       !> dust fraction
  !   type (dust)          :: mydust          !> cloud dust
  !   character(len=20)    :: dust_profile   !> could dust inside profile
  !   contains
  !    PRIVATE
  !    PROCEDURE, PASS(self) :: set_default     => usr_cloud_set_default
  !    PROCEDURE, PASS(self) :: set_complet     => usr_cloud_set_complet
  !    PROCEDURE, PASS(self) :: normalize       => usr_cloud_normalize
  !    PROCEDURE, PASS(self) :: set_w           => usr_cloud_set_w
  !    PROCEDURE, PASS(self) :: read_parameters => usr_cloud_read_p
  !    PROCEDURE, PASS(self) :: write_setting   => usr_cloud_write_setting
  !    PROCEDURE, PASS(self) :: clean_memory    => usr_cloud_clean_memory
  ! end type



  type (ISM)           :: ism_orion
  type (cloud),target  :: cloud_ary
  type (dust),target   :: dust_ary

  ! Star features
  ! type  star
  ! character(len=20)  ::unit  !> physical unit at parameter file
  ! real(dp)  :: radius     !> star radius (cm)
  ! real(dp)  :: mass       !> star mass  (g)
  ! real(dp)  :: luminosity !> star luminosity in erg/s
  ! real(dp)  :: Eddington  !> Eddington limit erg/s
  ! real(dp)  :: temperature!> star temperature
  ! real(dp)  :: vrotation  !> rotation parameter cm/s
  ! real(dp)  :: magnetic   !> Magnetic field strength at star surface (gauss)
  ! real(dp)  :: eta
  ! real(dp)  :: frac_critical_rotation
  ! end type
  ! type(star) :: star_ms
  ! type(star) :: sun
  !
  !
  !
  ! ! local normalisation
  !
  ! type unit_expand
  ! real(dp)  :: volum
  ! real(dp)  :: mass
  ! real(dp)  :: luminosity
  ! end type  unit_expand
  ! type(unit_expand) :: unit_user
  !
  ! logical   :: use1Dfile
  ! ! physical constantes:
  ! type const_expand
  ! real(dp)  :: G      = 6.67259D-8      ! cm^3 g^-1 s^-2
  ! real(dp)  :: clight = 2.99792458d10   ! cm s^-1
  ! end type const_expand
  ! type(const_expand) :: constusr
  !
  ! ! solar constantes
  ! real(dp) :: solar_mass       = 1.9892d+33 !> solar mass cgs
  ! real(dp) :: solar_luminosity = 3.826d+33  !> solar luminosity (erg/s)
  ! real(dp) :: solar_radius     = 6.95987d+10!> solar radius (cm)
  ! ! geometry
  ! character(len=30):: coordinate_system




contains
  subroutine usr_init()
    ! configuration of procedures to be used in this project
    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_source          => specialsource_usr
    usr_refine_grid     => specialrefine_usr
    usr_special_global  => usr_global_var
    usr_process_grid    => process_grid_usr

    call ism_orion%set_default
    call cloud_ary%set_default



    call usr_params_read(par_files)
    call set_coordinate_system(trim(coordinate_system))
    call hd_activate()
    call usr_physical_unit()!

    call usr_check_confilt()
    phys_n_tracer=hd_n_tracer

  end subroutine usr_init
  !------------------------------------------------------------------
  !> Read this module s parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ ism_on,cloud_on,cloud_number, &
                        coordinate_system


    if(mype==0)write(*,*)'Reading usr_list'
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=110)
110    close(unitpar)
    end do
    call usr_unit_read(files)
    if(ism_on)call ism_orion%read_parameters(files,ism_orion%unit,           &
                                  ism_orion%density,ism_orion%number_density,&
                                  ism_orion%temperature,ism_orion%pressure,  &
                                  ism_orion%velocity,                        &
                                  ism_orion%tracer_on,ism_orion%dust_on,     &
                                  ism_orion%dust_frac)

   if(cloud_on)call cloud_ary%read_parameters(files,cloud_ary%unit,         &
                             cloud_ary%center, cloud_ary%extend,            &
                             cloud_ary%density,cloud_ary%number_density,    &
                             cloud_ary%mass,cloud_ary%temperature,          &
                             cloud_ary%pressure,cloud_ary%velocity,         &
                             cloud_ary%eject_angle,cloud_ary%shape,         &
                             cloud_ary%profile,         &
                             cloud_ary%tracer_on,cloud_ary%dust_on,         &
                             cloud_ary%dust_frac,cloud_ary%dust_profile)
  end subroutine usr_params_read
!-------------------------------------------------------------------
!> subroutine read unit used in the code
  subroutine usr_unit_read(files)
   implicit none
   character(len=*), intent(in) :: files(:)
   integer                      :: n

   namelist /usr_unit_list/ unit_length , unit_time,unit_velocity,          &
                      unit_density, unit_numberdensity,               &
                      unit_pressure,unit_temperature



  if(mype==0)write(*,*)'Reading usr_unit_list'
  do n = 1, size(files)
         open(unitpar, file=trim(files(n)), status="old")
         read(unitpar, usr_unit_list, end=109)
  109    close(unitpar)
  end do
 end subroutine usr_unit_read
  !-----------------------------------------------------------
  !> subroutine to check configuration conflits
  subroutine usr_check_confilt
    implicit none
    cond_dust_on : if(.not.hd_dust)then
      ism_orion%dust_on =.false.
      cloud_ary%dust_on =.false.
    end if  cond_dust_on
  end   subroutine usr_check_confilt
  !-----------------------------------------------------------
  !> subroutine to normalize parameters in the code
  subroutine usr_normalise_parameters()
   implicit none
   integer            :: idust
   constusr%G         = constusr%G*&
                      (unit_density*(unit_length/unit_velocity)**(2.0_dp))


   constusr%clight       = constusr%clight/unit_velocity

    w_convert_factor(rho_)              = unit_density
    if(hd_energy)w_convert_factor(e_)   = unit_density*unit_velocity**2.0
    if(saveprim)then
     w_convert_factor(mom(:))           = unit_velocity
    else
     w_convert_factor(mom(:))           = unit_density*unit_velocity
    end if
    time_convert_factor                 = unit_time
    length_convert_factor               = unit_length
    if(hd_dust)then
          w_convert_factor(dust_rho(:))              = unit_density
          if(saveprim)then
            do idust = 1, dust_n_species
             w_convert_factor(dust_mom(:,idust))     = unit_velocity
            end do
          else
            do idust = 1, dust_n_species
             w_convert_factor(dust_mom(:,idust))     = unit_density*unit_velocity
           end do
          end if
    end if

  end subroutine usr_normalise_parameters


!-------------------------------------------------------------------------
  subroutine initglobaldata_usr()
   use mod_variables


   call ism_orion%set_complet()
   call ism_orion%normalize()


   call cloud_ary%set_complet()
   call cloud_ary%normalize()


   call usr_normalise_parameters()
   if(mype==0)call usr_write_setting()
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    use mod_hd
    implicit none

    integer, intent(in)     :: ixI^L,ixO^L
    real(dp), intent(in)    :: x(ixI^S,1:ndim)
    real(dp), intent(inout) :: w(ixI^S,1:nw)
    !.. local ..
    real(dp)      :: res
    integer       :: ix^D,na,flag(ixG^T)
    logical, save :: first=.true.
    !-----------------------------------------

    if(first)then
      if(mype==0) then
        write(*,*)'cloud propagation for ary :-)'
      endif
      first=.false.
    endif

    itr=1
    ! set the ism
    call ism_orion%set_w(ixI^L,ixO^L,x,w)

    ! set one cloud
    call cloud_ary%set_w(ixI^L,ixO^L,x,w)


    ! put dust to zero in all others zones
    cond_dust_on : if(hd_dust) then
      allocate(dust_ary%patch(ixI^S))
      dust_ary%patch(ixO^S) =.false.
      if(ism_orion%dust_on)then
       ism_orion%mydust%patch(ixO^S)=.not.cloud_ary%patch(ixO^S)
       dust_ary%patch(ixO^S) = dust_ary%patch(ixO^S).or.ism_orion%mydust%patch(ixO^S)
      end if
      if(cloud_ary%dust_on)then
       dust_ary%patch(ixO^S) = (dust_ary%patch(ixO^S).or. &
                               cloud_ary%mydust%patch(ixO^S))
      end if
      call dust_ary%set_w_zero(ixI^L,ixO^L,x,w)
      call dust_ary%clean_memory()
    end if cond_dust_on

    ! check is if initial setting is correct
    call  phys_check_w(.true., ixI^L, ixO^L, w, flag)
    if(any(flag(ixO^S)/=0)) PRINT*,' is error',maxloc(abs(flag(ixO^S)))


    ! get conserved variables to be used in the code
    call phys_to_conserved(ixI^L,ixO^L,w,x)
    call ism_orion%clean_memory()
    call cloud_ary%clean_memory()

  end subroutine initonegrid_usr
!----------------------------------------------------------------

  subroutine specialsource_usr(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_dust
    implicit none

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    real(dp), intent(in)            :: qdt, qtC, qt
    real(dp), intent(in)            :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(dp), intent(inout)         :: w(ixI^S,1:nw)
    ! .. local ..
    real(dp)                        :: small_dust_rho,coef_correct
    integer                         :: ix^D,idir,idust
    logical, dimension(ixI^S)       :: patch_correct,patch
    !----------------------------------------------------------
    return

    small_dust_rho = ism_orion%mydust%min_limit_rel
    coef_correct   = 0.4_dp

    call phys_to_primitive(ixI^L,ixO^L,w,x)

    Loop_idust : do idust =1, dust_n_species
      where(w(ixO^S, dust_rho(idust))<max(small_dust_rho*w(ixO^S,rho_),&
           ism_orion%mydust%min_limit_abs))
          w(ixO^S, dust_rho(idust))= 0.8* min(small_dust_rho*w(ixO^S,rho_),&
               ism_orion%mydust%min_limit_abs)
          patch_correct(ixO^S) = .true.
      elsewhere
          patch_correct(ixO^S) = .false.
      end where





     Loop_idir2 : do idir = 1,ndir
      where(patch_correct(ixO^S))
              w(ixO^S, dust_mom(idir,idust))=0.0_dp
      elsewhere
       where(w(ixO^S, dust_rho(idust))>ism_orion%mydust%min_limit_abs)


       where(w(ixO^S, dust_mom(idir,idust))*w(ixO^S,mom(idir))<0.0_dp.and.&
           (w(ixO^S, dust_rho(idust))/w(ixO^S,rho_))>1.0d2)

          w(ixO^S, dust_rho(idust))=(1.0_dp-coef_correct)*w(ixO^S, dust_rho(idust))&
               +coef_correct*w(ixO^S,rho_)
          w(ixO^S, dust_mom(idir,idust))=(1.0_dp-coef_correct)&
               *w(ixO^S, dust_mom(idir,idust))&
               +coef_correct*w(ixO^S,mom(idir))
        end where
       end where
      end where
     end do Loop_idir2
    end do Loop_idust
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine specialsource_usr
  !-------------------------------------------------------------------------
  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in)     :: ixO^L, iB, ixI^L
    real(dp), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(dp), intent(inout) :: w(ixI^S,1:nw)


  end subroutine specialbound_usr


     !> Enforce additional refinement or coarsening
     !> One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
     !> you must set consistent values for integers refine/coarsen:
     !> refine = -1 enforce to not refine
     !> refine =  0 doesn't enforce anything
     !> refine =  1 enforce refinement
     !> coarsen = -1 enforce to not coarsen
     !> coarsen =  0 doesn't enforce anything
     !> coarsen =  1 enforce coarsen
     !> e.g. refine for negative first coordinate x < 0 as
     !> if (any(x(ix^S,1) < zero)) refine=1
     subroutine specialrefine_usr(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
       use mod_global_parameters
       integer, intent(in)          :: igrid, level, ixI^L, ixO^L
       real(dp), intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
       integer, intent(inout)       :: refine, coarsen
      !----------------------------------------

      ! resolve at initial time the cloud
       cond_first_step : if(qt==0.0_dp)then
         call usr_cloud_patch(ixI^L,ixO^L,x,cloud_ary)
         if(any(cloud_ary%patch(ixO^S))) then
           refine  = 1
           coarsen = -1
         end if
         call cloud_ary%clean_memory()
       else cond_first_step
       ! coarsen in the back of the moving cloud
       if(all(w(ixO^S,rho_)<min(small_density*100.0_dp,ism_orion%density/100.0_dp)&
          .or.all(w(ixO^S,rho_)<min(small_pressure*100.0_dp,ism_orion%pressure/100.0_dp))))then
         refine  =-1
         coarsen = 1
        end if
        if(any(w(ixO^S,mom(2))/w(ixO^S,rho_)>20))then
         refine  = 1
         coarsen = -1
        end if
       end if cond_first_step
       ! coarsen in ism
       !if(all(dabs(get_v2)))
     end subroutine specialrefine_usr
  !> special output
  subroutine specialvar_output(ixI^L,ixO^L,win,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_physics
    use mod_dust
    implicit none
    integer, intent(in)        :: ixI^L,ixO^L
    real(dp), intent(in)       :: x(ixI^S,1:ndim)
    real(dp)                   :: win(ixI^S,nw+nwauxio)
    real(dp)                   :: normconv(0:nw+nwauxio)
    real(dp)                   :: fdrag(ixI^S, 1:ndim, 1:dust_n_species)
    real(dp),dimension(ixI^S,1:ndim,1:dust_n_species)  :: deltav,new_mean_dv
    real(dp),dimension(ixI^S,dust_n_species)           :: vdust
    real(dp),dimension(ixI^S,ndim)                     :: vgas,dv
    real(dp), dimension(ixI^S)                         :: vt2,  fd, ptherm,cmax,mean_dv
    integer                                            :: idir,idust
    logical, dimension(ixI^S,dust_n_species)           :: new_dvflag
    logical, dimension(ixI^S)                          :: mean_dvflag,patch_drag
    integer, dimension(ixI^S,dust_n_species,2)         :: indices_median
    real(dp)                                           :: w(ixI^S,nw)
    !----------------------------------------------------

   w(ixI^S,1:nw)=win(ixI^S,1:nw)
  ! call phys_to_primitive(ixI^L,ixO^L,w,x)
  call phys_get_pthermal(w, x, ixI^L, ixI^L, ptherm)
  if(any(ptherm(ixO^S)<0))STOP 'is wrong'
  Loop_idir0 : do idir=1,ndim
    vgas(ixI^S,idir)=w(ixI^S,mom(idir))/w(ixI^S,rho_)
  end do Loop_idir0
  vt2(ixI^S) = 3.0d0*ptherm(ixI^S)/w(ixI^S, rho_)

  Loop_idir : do idir = 1, ndim
      Loop_idust : do idust = 1, dust_n_species
patch_drag(ixI^S) = w(ixI^S, dust_rho(idust)) > 1.0d-9&
      .and.dabs(w(ixI^S, dust_mom(idir, idust)))>smalldouble
if(any(patch_drag(ixI^S))) then
        where(patch_drag(ixI^S))
          vdust(ixI^S,idust)  = w(ixI^S, dust_mom(idir, idust)) / w(ixI^S, dust_rho(idust))
          deltav(ixI^S,idir,idust) = (vgas(ixI^S, idir)-vdust(ixI^S,idust))

          ! 0.75 from sticking coefficient
    !      where(dabs(deltav(ixI^S,idir,idust))>1e-8)
           fd(ixI^S)     = 0.75d0*w(ixI^S, dust_rho(idust))*w(ixI^S, rho_)*deltav(ixI^S,idir,idust) &
               / (dust_density(idust) * dust_size(idust))

          ! 0.75 from spherical grainvolume
          fd(ixI^S)     = -fd(ixI^S)*0.75d0*dsqrt(vt2(ixI^S) + deltav(ixI^S,idir,idust)**2)
    !    elsewhere
    !      fd(ixI^S)     = 0.0_dp
    !    end where
        elsewhere
          vdust(ixI^S,idust)=0.0
          deltav(ixI^S,idir,idust) =0.0

          fd(ixI^S) = 0.0d0
        end where
      else
                vdust(ixI^S,idust)       = 0.0_dp
                deltav(ixI^S,idir,idust) = 0.0_dp
                fd(ixI^S)                = 0.0_dp
      end if
        fdrag(ixI^S, idir, idust) = fd(ixI^S)
      end do Loop_idust
    end do Loop_idir

    Loop_idust2: do idust = 1, dust_n_species
     new_dvflag(ixO^S,idust)=.true.
     Loop_idir1 : do idir = 1,ndim

            dv(ixI^S,idir)=-deltav(ixI^S,idir,idust)
            !w(ixI^S, dust_mom(idir,idust)) - &
            !                       w(ixI^S,mom(idir))

            if(maxval(dv(ixI^S,idir))-minval(dv(ixI^S,idir))>smalldouble) then
             call usr_medianvalue_of_array(ixI^L, ixO^L,dv(ixI^S,idir),mean_dv,&
                                          indices_median,mean_dvflag)

             new_dvflag(ixO^S,idust)=mean_dvflag(ixO^S).and.new_dvflag(ixO^S,idust)
             new_mean_dv(ixO^S,idir,idust) = mean_dv(ixO^S)
            else
              new_dvflag(ixO^S,idust)=.true.
              new_mean_dv(ixO^S,idir,idust)=dv(ixO^S,idir)
            end if

     end do Loop_idir1

    end do Loop_idust2

    win(ixO^S,nw+1) = deltav(ixO^S,1,1)/unit_velocity
    win(ixO^S,nw+2) = deltav(ixO^S,2,1)/unit_velocity
    win(ixO^S,nw+3) = dsqrt(vt2(ixO^S)+deltav(ixO^S,2,1)**2.0_dp)*deltav(ixO^S,2,1)/unit_velocity
    win(ixO^S,nw+4) = w(ixO^S, dust_rho(1))*w(ixO^S, rho_)!fdrag(ixO^S, 1, 1)
    win(ixO^S,nw+5) = fdrag(ixO^S, 2, 1)




    do idir = 1, ndim
        call phys_get_cmax(w,x,ixI^L,ixO^L,idir,cmax)
        win(ixO^S,nw+5+idir) = cmax(ixO^S)
    end do

    where(new_dvflag(ixO^S,1))
         win(ixO^S,nw+8) =  1.0_dp
    elsewhere
         win(ixO^S,nw+8) = -1.0_dp
    end where


    win(ixO^S,nw+9) = new_mean_dv(ixO^S,1,1)
    win(ixO^S,nw+10) = new_mean_dv(ixO^S,2,1)

    normconv(nw+1) = unit_velocity
    normconv(nw+2) = unit_velocity


  end subroutine specialvar_output


  !> this subroutine is ONLY to be used for computing auxiliary variables
  !> which happen to be non-local (like div v), and are in no way used for
  !> flux computations. As auxiliaries, they are also not advanced
  subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! .. local ..
    real(dp), dimension(ixI^S)      :: mean_dv,sum_dv,fd,vt2,mean_fd
    real(dp)                        :: small_dust_rho,coef
    real(dp)                        :: dv(ixI^S,1:ndir)
    integer, dimension(ixI^S)       :: indices_dv
    logical, dimension(ixI^S)       :: mean_dvflag,new_dvflag,patch_correct,patch_slow
    logical, dimension(ixI^S)       :: new_dfflag,mean_dfflag
    integer                         :: idust, idir, ix1,ix2, it_diff,it_start_drag
    integer                         :: indices_median(ixI^S,dust_n_species,2)
    !---------------------------------------------------


    !return
    if(hd_dust) then
     small_dust_rho = ism_orion%mydust%min_limit_rel

     call phys_to_primitive(ixI^L,ixI^L,w,x)
     ! handel small density dust
     Loop_idust : do idust =1, dust_n_species
      where(w(ixI^S, dust_rho(idust))<max(small_dust_rho*w(ixI^S,rho_),&
         ism_orion%mydust%min_limit_abs))
        w(ixI^S, dust_rho(idust))= 0.8* min(small_dust_rho*w(ixI^S,rho_),&
             ism_orion%mydust%min_limit_abs)
        patch_correct(ixI^S) = .true.
      elsewhere
        patch_correct(ixI^S) = .false.
      end where
      ! handel large density dust
      where(w(ixI^S,rho_)<0.9*ism_orion%density)
       where(w(ixI^S, dust_rho(idust))>ism_orion%mydust%max_limit_rel*w(ixI^S,rho_))
        w(ixI^S, dust_rho(idust))=0.8*ism_orion%mydust%max_limit_rel*w(ixI^S,rho_)
        patch_slow(ixI^S) = .true.
       elsewhere
        patch_slow(ixI^S) =.false.
       end where
      end where
      ! handel large cmax in rarefied region
      ! do idir = 1, ndim
      !   call phys_get_cmax(w,x,ixI^L,ixO^L,idir,cmax)
      !
      ! end do


      new_dvflag(ixI^S)=.true.
      new_dfflag(ixI^S)=.true.

      vt2(ixI^S) = 3.0d0*w(ixI^S,e_)/w(ixI^S, rho_)
      Loop_idir1 : do idir = 1,ndim
       where(patch_correct(ixI^S))
        w(ixI^S, dust_mom(idir,idust))=0.0_dp
       end where
       where(patch_slow(ixI^S))
               w(ixI^S, dust_mom(idir,idust))=w(ixI^S,mom(idir))
       end where
      !   dv(ixI^S,idir)=w(ixI^S, dust_mom(idir,idust)) - &
      !                                 w(ixI^S,mom(idir))
      !
      !
      ! !   it_start_drag=5000
      ! !   ! where(dabs(dv(ixI^S,idir))<1.0d-11)
      ! !   !   w(ixI^S, dust_mom(idir,idust)) =w(ixI^S,mom(idir))
      ! !   ! endwhere
      ! !   test_todelete_1 : if(it<it_start_drag) then
      ! !     w(ixI^S,dust_mom(idir,idust))  = w(ixI^S,mom(idir))
      ! !   else test_todelete_1
      ! !    it_diff=20000000
      ! !    test_todelete : if(it<=it_diff)then
      ! !     coef=real((it-it_start_drag)-it_diff,kind=dp)/real(it_diff,kind=dp)
      ! !
      ! !     where(w(ixI^S,dust_rho(idust))>small_dust_rho.and.dabs(dv(ixI^S,idir))>smalldouble)
      ! !       w(ixI^S,dust_mom(idir,idust))  =   (w(ixI^S, dust_mom(idir,idust))  +&
      ! !              (dv(ixI^S,idir))*coef)
      ! !     end where
      ! !   end if test_todelete
      ! ! end if test_todelete_1
      !
      !   dv(ixI^S,idir)=w(ixI^S,mom(idir))- w(ixI^S, dust_mom(idir,idust))
      !
      !   call usr_medianvalue_of_array(ixI^L, ixO^L,dv(ixI^S,idir),mean_dv,&
      !                               indices_median,mean_dvflag)
      !
      !   where(.not.mean_dvflag(ixO^S).and.&
      !        .not.patch_correct(ixO^S).and.dabs(w(ixO^S,dust_mom(idir,idust)))>smalldouble)
      !    w(ixO^S,dust_mom(idir,idust)) = -mean_dv(ixO^S)+w(ixO^S,mom(idir))
      !   end where
      !
      !  dv(ixI^S,idir)=w(ixI^S,mom(idir)) - w(ixI^S, dust_mom(idir,idust))
      ! new_dvflag(ixI^S)=mean_dvflag(ixI^S).and.new_dvflag(ixI^S)
          ! 0.75 from sticking coefficient
       !
       !  where(dabs(dv(ixI^S,idir))>0.0_dp)
       !     fd(ixI^S)     = 0.75d0*w(ixI^S, dust_rho(idust))*w(ixI^S, rho_)*dv(ixI^S,idir) &
       !         / (dust_density(idust) * dust_size(idust))
       !
       !    ! 0.75 from spherical grainvolume
       !    fd(ixI^S)     = -fd(ixI^S)*0.75d0*dsqrt(vt2(ixI^S) + dv(ixI^S,idir)**2)
       !  elsewhere
       !    fd(ixI^S)     =0.0
       !  end where
       !  call usr_medianvalue_of_array(ixI^L, ixO^L,fd,mean_fd,&
       !                                      indices_median,mean_dfflag)
       !
       ! new_dfflag(ixI^S)=mean_dfflag(ixI^S).and.new_dfflag(ixI^S)
      end do Loop_idir1
    !  if(any(.not.new_dfflag(ixO^S)))call dust_average_dustdensity(ixI^L,ixO^L,idust,new_dfflag,w)
      !if(any(.not.new_dvflag(ixO^S)))call dust_average_dspeed(ixI^L,ixO^L,idust,new_dfflag,w)
      !if(any(.not.new_dvflag(ixO^S)))call dust_average_dspeed(ixI^L,ixO^L,idust,new_dfflag,w)
     end do  Loop_idust
     call phys_to_conserved(ixI^L,ixI^L,w,x)
   end if

  end subroutine process_grid_usr
  !---------------------------------------------------------------------
  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*), intent(inout) :: varnames(:)
    varnames(1)  = 'deltav1'
    varnames(2)  = 'deltav2'
    varnames(3)  = 'vtherm2'
    varnames(4)  = 'fdrag1'
    varnames(5)  = 'fdrag2'
    varnames(6)  = 'cmax1'
    varnames(7)  = 'cmax2'
    varnames(8)  = 'flag_dv'
    varnames(9)  = 'mean_dv11'
    varnames(10) = 'mean_dv12'
  end subroutine specialvarnames_output


  !---------------------------------------------------------------------
  !> subroutine to write simulation configuration
  subroutine usr_write_setting()
    integer,parameter   :: unit_config =12
    character(len=50)   :: filename_config
    !-------------------------------------
    filename_config=trim(base_filename)//'_config'
    open(unit_config,file=trim(filename_config), status='replace')
    call ism_orion%write_setting(unit_config)
    call cloud_ary%write_setting(unit_config)
    close(unit_config)
  end subroutine usr_write_setting
!---------------------------------------------------------------------

!  !> Read the ism parameters  from a parfile
!   subroutine usr_ism_read_p(self,files,unit,density,number_density,       &
!                             temperature,pressure,                         &
!                             velocity,                                     &
!                             tracer_on,dust_on,dust_frac)
!     class(ism)                     :: self
!     character(len=*),intent(in)    :: files(:)
!     character(len=20),intent(inout):: unit
!     logical, intent(inout)         :: tracer_on,dust_on
!     real(kind=dp), intent(inout)   :: density,number_density,            &
!                                       temperature,pressure,              &
!                                       velocity(1:3),                  &
!                                       dust_frac
!     integer                        :: n
!
!     namelist /usr_ism_list/ unit,density,number_density,pressure,           &
!                             velocity,temperature,            &
!                             tracer_on,dust_on,dust_frac
!
!     if(mype==0)write(*,*)'Reading usr_ism_list'
!     do n = 1, size(files)
!        open(unitpar, file=trim(files(n)), status="old")
!        read(unitpar, usr_ism_list, end=113)
! 113    close(unitpar)
!     end do
!     if(dust_on)call self%mydust%read_parameters(files,                       &
!                       self%mydust%unit,self%mydust%n_species,                &
!                       self%mydust%sizes,                                     &
!                       self%mydust%min_radius,self%mydust%max_radius,         &
!                       self%mydust%power_a,                                   &
!                       self%mydust%distrub_func,                              &
!                       self%mydust%idust_first, self%mydust%idust_last,       &
!                       self%mydust%grain_density,                             &
!                       self%mydust%density,self%mydust%number_density        ,&
!                       self%mydust%temperature,           &
!                       self%mydust%velocity,             &
!                       self%mydust%tracer_on,self%mydust%associate_med,       &
!                       self%mydust%min_limit_rel,self%mydust%min_limit_abs,   &
!                       self%mydust%max_limit_rel,self%mydust%extend)
!
!
!   end subroutine usr_ism_read_p
!  !------------------------------------------------------------------------
!  subroutine usr_ism_write_setting(self,unit_config)
!    implicit none
!    class(ism)                          :: self
!    integer,intent(in)                  :: unit_config
!    ! .. local ..
!
!    !-----------------------------------
!
!    write(unit_config,*)'************************************'
!    write(unit_config,*)'************ISM setting ************'
!    write(unit_config,*)'************************************'
!    write(unit_config,*) 'Density     = ', self%density
!    write(unit_config,*) 'Pressure    = ', self%pressure
!    write(unit_config,*) 'Temperature = ', self%temperature
!    write(unit_config,*) 'Speed       = ', self%velocity
!    call self%mydust%write_setting(unit_config)
!    write(unit_config,*)'************************************'
!    write(unit_config,*)'******** END ISM setting **********'
!    write(unit_config,*)'************************************'
!  end    subroutine usr_ism_write_setting
!  !-------------------------------------------------------------------------
!  !> subroutine default setting for ISM
!  subroutine usr_ism_set_default(self)
!   implicit none
!   class(ism)            :: self
!   !----------------------------------
!   self%unit               = 'code'
!   self%density            = 0.0_dp
!   self%number_density     = 0.0_dp
!   self%temperature        = 0.0_dp
!   self%pressure           = 0.0_dp
!   self%extend(1:2,1:ndim) = 0.0_dp!box_limit(1:2,1:ndim)
!   self%velocity           = 0.0_dp
!
!
!   self%tracer_on          = .false.
!   self%dust_on            = .false.
!   self%dust_frac          = 0.0_dp
!
!   call self%mydust%set_default()
!
!  end subroutine usr_ism_set_default
! !--------------------------------------------------------------------
! !> subroutine check the parfile setting for ism
! subroutine usr_ism_set_complet(self)
!   implicit none
!   class(ism)               :: self
!   ! .. local ..
!   logical                  :: dust_is_frac
!   real(dp)                 :: mp,kb
!   !-----------------------------------
!   if(SI_unit) then
!     mp=mp_SI
!     kB=kB_SI
!   else
!     mp=mp_cgs
!     kB=kB_cgs
!   end if
!
!   if (dabs(self%density)<smalldouble*mp)then
!       self%density=self%number_density*mp
!   end if
!
!   if(dabs(self%pressure)<smalldouble) then
!     self%pressure =(2.d0+3.d0*He_abundance)*self%number_density*kB*self%temperature
!   end if
!
!   if(self%dust_on)then
!     dust_is_frac=.false.
!     call self%mydust%set_complet(dust_is_frac,self%dust_frac,self%density,self%velocity)
!   end if
! end subroutine usr_ism_set_complet
! !--------------------------------------------------------------------
! !> subroutine normalize setting for ISM
!  subroutine usr_ism_normalize(self)
!   implicit none
!   class(ism)            :: self
!   !----------------------------------
!   if(trim(self%unit)=='code')return
!   self%density          = self%density/unit_density
!   self%temperature      = self%temperature/unit_temperature
!   self%pressure         = self%pressure/unit_pressure
!   self%velocity         = self%velocity/unit_velocity
!   self%extend           = self%extend/unit_length
!   if(self%dust_on)then
!     call self%mydust%normalize()
!     call self%mydust%to_phys()
!   end if
!   if(mype==0)then
!    WRITE(*,*)'********************************************************'
!    WRITE(*,*)' is in ISM', self%pressure , self%density,self%temperature,sqrt(self%pressure / self%density)
!    WRITE(*,*)'********************************************************'
!   end if
!  end subroutine usr_ism_normalize
! !--------------------------------------------------------------------
!  !> subroutine setting for ISM
!  subroutine usr_ism_set_w(ixI^L,ixO^L,x,w,self)
!   implicit none
!   integer, intent(in)  :: ixI^L,ixO^L
!   real(kind=dp)        :: x(ixI^S,1:ndir)
!   real(kind=dp)        :: w(ixI^S,1:nw)
!   class(ism)           :: self
!   ! .. local..
!   integer              :: idir
!   real(kind=dp)        :: fprofile(ixI^S)
!   !----------------------------------
!   allocate(self%patch(ixG^T))
!
!   self%patch              = .true.
!   where(self%patch(ixO^S))
!    w(ixO^S,rho_) = self%density
!    w(ixO^S,p_)   = self%pressure
!   end where
!   Loop_idir : do idir=1,ndir
!     where(self%patch(ixO^S))
!      w(ixO^S,mom(idir)) = self%velocity(idir)
!     end where
!   end do Loop_idir
!   cond_tracer_on :if(self%tracer_on.and.phys_n_tracer>0 &
!                      .and.itr<=phys_n_tracer)then
!    where(self%patch(ixO^S))
!     w(ixO^S,tracer(itr)) = w(ixO^S,rho_)
!    elsewhere
!     w(ixO^S,tracer(itr)) = 0.0_dp
!    end where
!    itr=itr+1
!   end if cond_tracer_on
!
!   cond_dust_on : if(self%dust_on)then
!    allocate(self%mydust%patch(ixG^T))
!    self%mydust%patch=self%patch
!    self%mydust%velocity=self%velocity
!    fprofile = 1.0_dp
!    call self%mydust%set_w(ixI^L,ixO^L,.false.,self%dust_frac,fprofile,x,w)
!    !deallocate(self%mydust%patch)
!   end if cond_dust_on
!
!
!  end subroutine usr_ism_set_w
!
!
!  !--------------------------------------------------------------------
!  !> Subroutine to clean array memory of associated with cloud object
!  subroutine usr_ism_clean_memory(self)
!    class(ism)    :: self
!    if(allocated(self%patch))deallocate(self%patch)
!    if(self%dust_on)call self%mydust%clean_memory()
!  end subroutine usr_ism_clean_memory
! !---------------------------------------------------------------------
!
!  !> Read the ism parameters  from a parfile
!   subroutine usr_dust_read_p(self,files,unit,n_species,sizes,              &
!                              min_radius,max_radius,power_a,distrub_func,   &
!                              idust_first,idust_last                    ,   &
!                              grain_density,density,number_Density      ,   &
!                              temperature,velocity                      ,   &
!                              tracer_on,associate_med,                      &
!                              min_limit_rel,min_limit_abs,max_limit_rel ,   &
!                              extend)
!     class(dust)                     :: self
!     character(len=*), intent(in)    :: files(:)
!     character(len=20), intent(inout):: unit,distrub_func
!     integer, intent(inout)          :: n_species,idust_first,idust_last
!     logical, intent(inout)          :: tracer_on,associate_med
!     real(kind=dp), intent(inout)    :: grain_Density(n_dust_max),         &
!                                        density(n_dust_max),               &
!                                        number_Density(n_dust_max),        &
!                                        temperature,                       &
!                                        velocity(3),power_a,               &
!                                        min_limit_rel,min_limit_abs,       &
!                                        max_limit_rel ,                    &
!                                        extend(ndim),sizes(n_dust_max),    &
!                                        min_radius,max_radius
!
!     integer                         :: n
!
!     namelist /usr_dust_list/ unit,grain_Density,density,       &
!                             velocity,temperature,              &
!                             tracer_on,power_a,                 &
!                             min_limit_rel,min_limit_abs,       &
!                             max_limit_rel,                     &
!                             extend,n_species,sizes,            &
!                             min_radius,max_radius,distrub_func,&
!                             idust_first, idust_last,associate_med
!
!     if(mype==0)write(*,*)'Reading usr_list'
!     do n = 1, size(files)
!        open(unitpar, file=trim(files(n)), status="old")
!        read(unitpar, usr_dust_list, end=111)
! 111    close(unitpar)
!     end do
!
!   end subroutine usr_dust_read_p
! !------------------------------------------------------------------------
! !> write the cloud setting
! subroutine usr_dust_write_setting(self,unit_config)
!   implicit none
!   class(dust)                         :: self
!   integer,intent(in)                  :: unit_config
!   ! .. local ..
!
!   !-----------------------------------
!
!   write(unit_config,*)'**************************'
!   write(unit_config,*)'**** Dust  setting *******'
!   write(unit_config,*)'**************************'
!   write(unit_config,*) 'Density       = ', self%density(1:self%idust_last-self%idust_first+1)
!   write(unit_config,*) 'Grain density = ', self%grain_density(1:self%idust_last-self%idust_first+1)
!   write(unit_config,*) 'Temperature   = ', self%temperature
!   write(unit_config,*) 'Speed         = ', self%velocity
!   write(unit_config,*) 'n_species     = ', self%n_species
!   write(unit_config,*)'**************************'
!   write(unit_config,*)'**** End dust setting ***'
!   write(unit_config,*)'**************************'
! end    subroutine usr_dust_write_setting
! !--------------------------------------------------------------------
! !> subroutine default setting for ISM
!  subroutine usr_dust_set_default(self)
!   implicit none
!   class(dust)          :: self
!   !----------------------------------
!
!   self%unit               = 'code'
!   self%distrub_func       = 'uniform'
!   self%n_species          = 0
!   self%sizes              = 0.0_dp
!   self%density            = 0.0_dp
!   self%number_density     = 0.0_dp
!   self%temperature        = 0.0_dp
!   self%velocity           = 0.0_dp
!   self%min_limit_rel      = 0.0_dp
!   self%min_limit_abs      = 0.0_dp
!   self%max_limit_rel      = 1.0d9
!   self%extend(1:2,1:ndim) = 0.0_dp!box_limit(1:2,1:ndim) ! it zero at this time
!   self%associate_med      = .true.
!   self%grain_Density      = 0.0_dp
!   self%min_radius         =-1.0_dp
!   self%max_radius         =-1.0_dp
!   self%idust_first        = 1
!   self%idust_last         = 1
!   self%tracer_on          = .false.
!  end subroutine usr_dust_set_default
! !--------------------------------------------------------------------
! !> subroutine setting for dust
! subroutine usr_dust_set_w_zero(ixI^L,ixO^L,x,w,self)
!   use mod_dust
!  implicit none
!  integer, intent(in)          :: ixI^L,ixO^L
!  real(kind=dp), intent(in)    :: x(ixI^S,1:ndir)
!  real(kind=dp), intent(inout) :: w(ixI^S,1:nw)
!  class(dust)                  :: self
!  ! .. local ..
!  integer                      :: idust,idir
!  !----------------------------------------------------
!  if(all(self%patch(ixO^S)))return
!   Loop_idust01 : do idust = 1,dust_n_species
!     where(.not.self%patch(ixO^S))
!       w(ixO^S,dust_rho(idust)) = 1.1d-2*self%min_limit_abs
!     end where
!     Loop_idir : do idir=1,ndir
!         where(.not.self%patch(ixO^S))
!          w(ixO^S,dust_mom(idir,idust)) = self%velocity(idir)
!         end where
!     end do Loop_idir
!   end do Loop_idust01
!
! end  subroutine usr_dust_set_w_zero
! !--------------------------------------------------------------------
!  !> subroutine setting W for dust
!  subroutine usr_dust_set_w(ixI^L,ixO^L,is_frac,dust_frac,fprofile,x,w,self)
!   implicit none
!   integer, intent(in)          :: ixI^L,ixO^L
!   logical                      :: is_frac
!   real(kind=dp), intent(in)    :: dust_frac
!   real(kind=dp), intent(in)    :: fprofile(ixI^S)
!   real(kind=dp), intent(in)    :: x(ixI^S,1:ndir)
!   real(kind=dp), intent(inout) :: w(ixI^S,1:nw)
!   class(dust)                  :: self
!   ! .. local..
!   real(kind=dp)                :: slop,dr_dust,dust_frac_loc
!   integer                      :: idir,idust,istart_dust,istop_dust,i_dust_loc
!   !----------------------------------
!
!   if(.not.any(self%patch(ixO^S)))return
!   istart_dust=1
!   istop_dust = (istart_dust-1)+self%n_species
!   if(is_frac) then
!     dust_frac_loc=dust_frac/(1.0_dp-dust_frac)
!   else
!     dust_frac_loc=dust_frac
!   end if
!
!
!    select case(self%distrub_func)
!    case('uniform')
!     cond_fracdust : if(dust_frac_loc<0.0_dp) then
!
!      Loop_idust2 : do idust = self%idust_first, self%idust_last
!       i_dust_loc = idust-self%idust_first+1
!       where(self%patch(ixO^S))
!        w(ixO^S,dust_rho(idust)) = self%density(i_dust_loc)/real(self%n_species,kind=dp)
!       end where
!     end do Loop_idust2
!     else cond_fracdust
!
!      Loop_idust3 : do idust = self%idust_first, self%idust_last
!       where(self%patch(ixO^S))
!        w(ixO^S,dust_rho(idust)) = dust_frac_loc*w(ixO^S,rho_)&
!                                   / real(self%n_species,kind=dp)
!       end where
!     end do Loop_idust3
!
!
!     end if cond_fracdust
!
!    case('powerlaw')
!
!
!
!     Loop_idust4 : do idust = self%idust_first, self%idust_last
!      dr_dust=0.5_dp*(dust_size(max(idust-1,self%idust_first))-&
!                      dust_size(min(idust+1,self%idust_last)))
!      where(self%patch(ixO^S))
!       w(ixO^S,dust_rho(idust)) = w(ixO^S,rho_)*dust_frac_loc*&
!          half*(dr_dust/dsqrt(dust_size(idust))&
!          / (dsqrt(dust_size(self%idust_first))-dsqrt(dust_size(self%idust_last))))
!       w(ixO^S,dust_rho(idust)) = w(ixO^S,dust_rho(idust)) *fprofile(ixO^S)
!      end where
!     end do Loop_idust4
!
!    case default
!     write(*,*)'The ', self%distrub_func, ' is not implimented'
!     call mpistop('is not implimented in usr_dust_set_w at mod_usr.t')
!    end select
!
!
!
!   cond_isassociated : if(self%associate_med)then
!    Loop_idust5 : do idust = self%idust_first, self%idust_last
!     Loop_idir0 : do idir=1,ndir
!     where(self%patch(ixO^S))
!      w(ixO^S,dust_mom(idir,idust)) = w(ixO^S,mom(idir))
!     end where
!    end do Loop_idir0
!   end do Loop_idust5
!   else cond_isassociated
!    Loop_idust6 : do idust = self%idust_first, self%idust_last
!     Loop_idir : do idir=1,ndir
!      where(self%patch(ixO^S))
!      w(ixO^S,dust_mom(idir,idust)) = self%velocity(idir)
!      end where
!     end do Loop_idir
!    end do Loop_idust6
!   end if cond_isassociated
!
!  end subroutine usr_dust_set_w
!  !--------------------------------------------------------------------
!  !> subroutine check the parfile setting for cloud
!  subroutine usr_dust_set_complet(self,is_frac,dust_frac,med_density,med_velocity)
!    implicit none
!    class(dust)              :: self
!    logical                  :: is_frac
!    real(dp), intent(in)     :: dust_frac,med_density,med_velocity(1:ndir)
!    ! .. local ..
!    integer                  :: idust,i_dust_loc
!    real(dp)                 :: mp,kb,slop,dr_dust,dust_frac_loc
!    !-----------------------------------
!    if(SI_unit) then
!      mp=mp_SI
!      kB=kB_SI
!    else
!      mp=mp_cgs
!      kB=kB_cgs
!    end if
!    self%idust_last=(self%idust_first-1) +self%n_species
!   if(is_frac) then
!     dust_frac_loc=dust_frac/(1.0_dp-dust_frac)
!   else
!     dust_frac_loc=dust_frac
!   end if
!    cond_size_off : if(any(self%sizes(1:self%n_species)<=smalldouble))then
!
!     select case(self%distrub_func)
!     case('uniform')
!
!      slop=(self%max_radius-self%min_radius)/dble(self%n_species-1)
!      Loop_idust1 : do idust = self%idust_first, self%idust_last
!        i_dust_loc = idust-self%idust_first+1
!        self%sizes(i_dust_loc)   = self%min_radius+slop*real(idust-1,kind=dp)
!        self%density(i_dust_loc) = &
!             dust_frac_loc* med_density/real(self%n_species,kind=dp)
!      end do Loop_idust1
!
!     case('powerlaw')
!      call usr_mat_powerlaw_withX(self%n_species,self%power_a,self%min_radius,&
!                                  self%max_radius,   &
!                                  self%sizes)
!
!
!      Loop_idust2  : do idust = self%idust_first, self%idust_last
!       i_dust_loc = idust-self%idust_first+1
!       dr_dust=0.5_dp*(self%sizes(max(i_dust_loc-1,1))-&
!                       self%sizes(min(i_dust_loc+1,self%idust_last-self%idust_first+1)))
!       self%density(i_dust_loc) =med_density*dust_frac_loc*&
!              half*(dr_dust/dsqrt(self%sizes(i_dust_loc))&
!         / (dsqrt(self%sizes(1))-dsqrt(self%sizes(self%idust_last-self%idust_first+1))))
!
!      end do  Loop_idust2
!     case default
!      write(*,*)'The ', self%distrub_func, ' is not implimented'
!      call mpistop('is not implimented in usr_dust_set_w at mod_usr.t')
!     end select
!    end if cond_size_off
!    if(self%associate_med)then
!      self%velocity = med_velocity
!    end if
!  end subroutine usr_dust_set_complet
!  !--------------------------------------------------------------------
!  !> subroutine normalize setting for ISM
!   subroutine usr_dust_to_phys(self)
!     implicit none
!     class(dust)             :: self
!     ! .. local ..
!     integer                 :: idust,i_dust_loc
!     !----------------------------------
!
!     Loop_idust1 : do idust = self%idust_first, self%idust_last
!       i_dust_loc            = idust-self%idust_first+1
!       dust_size(idust)      = self%sizes(i_dust_loc)
!       dust_density(idust)   = self%grain_density(i_dust_loc)
!     end do Loop_idust1
!
!
!
!   end subroutine usr_dust_to_phys
! !--------------------------------------------------------------------
! !> subroutine normalize setting for ISM
!  subroutine usr_dust_normalize(self)
!   implicit none
!   class(dust)          :: self
!   !----------------------------------
!   if(trim(self%unit)=='code')return
!   self%sizes            = self%sizes         / unit_length
!   self%density          = self%density       / unit_density
!   self%number_density   = self%number_density/ unit_length**3.0_dp
!   self%grain_density    = self%grain_density / unit_density
!   self%temperature      = self%temperature   / unit_temperature
!   self%velocity         = self%velocity      / unit_velocity
!   self%min_limit_abs    = self%min_limit_abs / unit_density
!   self%extend           = self%extend        / unit_length
!  end subroutine usr_dust_normalize
!  !--------------------------------------------------------------------
!  !> Subroutine to clean array memory of associated with cloud object
!  subroutine usr_dust_clean_memory(self)
!    class(dust)    :: self
!    if(allocated(self%patch))deallocate(self%patch)
!  end subroutine usr_dust_clean_memory
! !-------------------------------------------------------------------------
!  !> Read the ism parameters  from a parfile
!   subroutine usr_cloud_read_p(self,files,unit,center, extend,                 &
!                               density,number_density, mass,                   &
!                               temperature,pressure,                           &
!                               velocity,eject_angle, shape, profile,           &
!                               tracer_on,dust_on,dust_frac,dust_profile)
!     class(cloud)                    :: self
!     character(len=*), intent(in)    :: files(:)
!     character(len=20), intent(inout):: unit
!     character(len=20),intent(inout) :: shape,profile,dust_profile
!     real(kind=dp), intent(inout)    :: center(ndim), extend(ndim),            &
!                                        density,number_Density,mass,           &
!                                        temperature,pressure,                  &
!                                        velocity(3),                           &
!                                        dust_frac,eject_angle
!     logical, intent(inout)          :: dust_on,tracer_on
!     integer                         :: n
!
!     namelist /usr_cloud_list/ unit,center, extend,                            &
!                               density,number_density,mass,                    &
!                               pressure,velocity,temperature,                  &
!                               tracer_on,dust_on,dust_frac,shape,profile,      &
!                               eject_angle,dust_profile
!
!     if(mype==0)write(*,*)'Reading usr_cloud_list'
!     do n = 1, size(files)
!        open(unitpar, file=trim(files(n)), status="old")
!        read(unitpar, usr_cloud_list, end=112)
! 112    close(unitpar)
!     end do
!
!     if(dust_on)call self%mydust%read_parameters(files,&
!                       self%mydust%unit,self%mydust%n_species,                &
!                       self%mydust%sizes,                                     &
!                       self%mydust%min_radius,self%mydust%max_radius,         &
!                       self%mydust%power_a,                                   &
!                       self%mydust%distrub_func,                              &
!                       self%mydust%idust_first, self%mydust%idust_last,       &
!                       self%mydust%grain_density,    &
!                       self%mydust%density,self%mydust%number_density        ,&
!                       self%mydust%temperature,           &
!                       self%mydust%velocity,          &
!                       self%mydust%tracer_on,self%mydust%associate_med,       &
!                       self%mydust%min_limit_rel,self%mydust%min_limit_abs,   &
!                       self%mydust%max_limit_rel,self%mydust%extend)
!   end subroutine usr_cloud_read_p
!
! !------------------------------------------------------------------------
! !> write the cloud setting
! subroutine usr_cloud_write_setting(self,unit_config)
!   implicit none
!   class(cloud)                        :: self
!   integer,intent(in)                  :: unit_config
!   ! .. local ..
!
!   !-----------------------------------
!
!   write(unit_config,*)'************************************'
!   write(unit_config,*)'************CLOUD setting ************'
!   write(unit_config,*)'************************************'
!   write(unit_config,*) 'Density     = ', self%density
!   write(unit_config,*) 'Pressure    = ', self%pressure
!   write(unit_config,*) 'Temperature = ', self%temperature
!   write(unit_config,*) 'Speed       = ', self%velocity
!   call self%mydust%write_setting(unit_config)
!   write(unit_config,*)'************************************'
!   write(unit_config,*)'******** END CLOUD setting **********'
!   write(unit_config,*)'************************************'
! end    subroutine usr_cloud_write_setting
! !--------------------------------------------------------------------
!
! !> subroutine default setting for cloud
!  subroutine usr_cloud_set_default(self)
!   implicit none
!   class(cloud)          :: self
!   !----------------------------------
!   self%unit             = 'code'
!   self%density          = 0.0_dp
!   self%number_density   = 0.0_dp
!   self%mass             = 0.0_dp
!   self%temperature      = 0.0_dp
!   self%pressure         = 0.0_dp
!   self%center(:)        = 0.0_dp!(box_limit(2,:)-box_limit(1,:))/2.0_dp
!   self%extend(:)        = 0.0_dp!(box_limit(2,:)+box_limit(1,:))/2.0_dp
!   self%shape            = 'sphere'
!   self%profile          = 'uniform'
!   self%eject_angle      =  0.0_dp
!   self%tracer_on        = .false.
!   self%dust_on          = .false.
!   self%dust_frac        = 0.0_dp
!   self%dust_profile     = 'uniform'
!   call self%mydust%set_default()
!  end subroutine usr_cloud_set_default
!  !--------------------------------------------------------------------
!  !> subroutine check the parfile setting for cloud
!  subroutine usr_cloud_set_complet(self)
!    implicit none
!    class(cloud)             :: self
!    ! .. local ..
!    logical                  :: dust_is_frac
!    real(dp)                 :: mp,kb,cloud_volume
!    !-----------------------------------
!    if(SI_unit) then
!      mp=mp_SI
!      kB=kB_SI
!    else
!      mp=mp_cgs
!      kB=kB_cgs
!    end if
!    call usr_get_volume(self%extend,self%shape, cloud_volume)
!    PRINT*,' is test volume',cloud_volume,self%dust_frac,self%mass/cloud_volume,&
!     self%mass/(cloud_volume*unit_density),cloud_volume/unit_length**3.0_dp
!    if (dabs(self%mass)>smalldouble)then
!     self%density        = self%mass*(1.0_dp-self%dust_frac)/cloud_volume
!     self%number_density = self%density/mp
!   else if (dabs(self%density)<smalldouble*mp)then
!     self%density        = self%number_density*mp
!     if(self%dust_frac<1.0_dp)then
!        self%mass        = self%density*cloud_volume/(1.0_dp-self%dust_frac)
!      else
!        self%mass        = 0.0_dp
!      end if
!    else
!     if (dabs(self%number_density)<smalldouble*mp)then
!       self%number_density        = self%density/mp
!     end if
!     if(self%dust_frac<1.0_dp)then
!       self%mass        = self%density*cloud_volume/(1.0_dp-self%dust_frac)
!     else
!       self%mass        = 0.0_dp
!     end if
!    end if
!
!    if(dabs(self%pressure)<smalldouble) then
!     self%pressure =(2.d0+3.d0*He_abundance)*self%number_density*kB*self%temperature
!    end if
!
!
!
!
!
!    if(self%dust_on)then
!     if(self%mass>smalldouble) then
!       dust_is_frac=.true.
!     else
!       dust_is_frac=.false.
!     end if
!      call self%mydust%set_complet(dust_is_frac,self%dust_frac,self%density,self%velocity)
!    end if
!
!  end subroutine usr_cloud_set_complet
! !--------------------------------------------------------------------
!  subroutine usr_cloud_normalize(self)
!   implicit none
!   class(cloud)          :: self
!   !----------------------------------
!   if(trim(self%unit)=='code')return
!   self%density          = self%density       /unit_density
!   self%temperature      = self%temperature   /unit_temperature
!   self%pressure         = self%pressure      /unit_pressure
!   self%velocity         = self%velocity      /unit_velocity
!   self%eject_angle      = self%eject_angle   *(dpi/180._dp)
!   self%center           = self%center        /unit_length
!   self%extend           = self%extend        /unit_length
!   if(self%dust_on)then
!     call self%mydust%normalize()
!     call self%mydust%to_phys()
!   end if
!
!  end subroutine usr_cloud_normalize
! !--------------------------------------------------------------------
! !--------------------------------------------------------------------
!  !> subroutine patch for the cloud
!  subroutine usr_cloud_patch(ixI^L,ixO^L,x,self)
!   implicit none
!   integer, intent(in)  :: ixI^L,ixO^L
!   real(kind=dp)        :: x(ixI^S,1:ndir)
!   class(cloud)         :: self
!   !----------------------------------
!
!   allocate(self%patch(ixG^T))
!
!   self%patch              = .false.
!   select case(trim(self%shape))
!   case('sphere')
!     call usr_set_patch_sphere(ixI^L,ixO^L,typeaxial,self%center,    &
!                               self%extend,x,self%patch)
!
!   case('cylinder')
!     call usr_set_patch_cylinder(ixI^L,ixO^L,typeaxial,self%center,  &
!                                 self%extend,x,self%patch)
!   case('cube')
!     call usr_set_patch_cube(ixI^L,ixO^L,typeaxial,self%center,      &
!                             self%extend,x,self%patch)
!   case default
!      write(*,*)'this cloud shape ',trim(self%shape),' is not implimented'
!      call mpistop('This cloud shape is not implimented in usr_cloud_patch at mod_usr.t')
!   end select
!  end subroutine usr_cloud_patch
! !--------------------------------------------------------------------
!  !> subroutine setting for cloud
!  subroutine usr_cloud_set_w(ixI^L,ixO^L,x,w,self)
!   implicit none
!   integer, intent(in)  :: ixI^L,ixO^L
!   real(kind=dp)        :: x(ixI^S,1:ndir)
!   real(kind=dp)        :: w(ixI^S,1:nw)
!   class(cloud)         :: self
!   ! .. local..
!   integer              :: idir
!   logical              :: dust_is_frac
!   real(kind=dp)        :: fprofile(ixI^S)
!   !----------------------------------
!   call usr_cloud_patch(ixI^L,ixO^L,x,self)
!   !cond_inside_cloud: if(any(self%patch(ixO^S))) then
!   where(self%patch(ixO^S))
!    w(ixO^S,rho_) = self%density
!    w(ixO^S,p_)   = self%pressure
!   end where
!   Loop_idir : do idir=1,ndir
!       where(self%patch(ixO^S))
!         w(ixO^S,mom(idir)) = self%velocity(idir)
!       end where
!   end do Loop_idir
!   if(self%profile/='none') then
!    call usr_mat_profile(ixI^L,ixO^L,typeaxial,self%profile,self%center,self%extend,&
!                         x,fprofile)
!    where(self%patch(ixO^S))
!      w(ixO^S,rho_)=w(ixO^S,rho_)*fprofile(ixO^S)
!      !w(ixO^S,p_)=w(ixO^S,p_)*fprofile(ixO^S)
!    end where
!   end if
!   if(self%dust_profile/='none') then
!    call usr_mat_profile(ixI^L,ixO^L,typeaxial,self%dust_profile,self%center,self%extend,&
!                         x,fprofile)
!   end if
!
!   cond_tracer_on :if(self%tracer_on.and.phys_n_tracer>0&
!                     .and.itr<=phys_n_tracer)then
!    where(self%patch(ixO^S))
!     w(ixO^S,tracer(itr)) = w(ixO^S,rho_)
!    elsewhere
!     w(ixO^S,tracer(itr)) = 0.0_dp
!    end where
!    itr=itr+1
!   end if cond_tracer_on
!
!   cond_dust_on : if(self%dust_on)then
!    allocate(self%mydust%patch(ixG^T))
!    self%mydust%patch=self%patch
!    self%mydust%velocity=self%velocity
!    if(self%mass>smalldouble) then
!      dust_is_frac=.true.
!    else
!      dust_is_frac=.false.
!    end if
!    call self%mydust%set_w(ixI^L,ixO^L,dust_is_frac,self%dust_frac,fprofile,x,w)
!   end if cond_dust_on
!
!   !deallocate(self%patch)
!
!  end subroutine usr_cloud_set_w
!
! !--------------------------------------------------------------------
! !> Subroutine to clean array memory of associated with cloud object
! subroutine usr_cloud_clean_memory(self)
!   class(cloud)    :: self
!   if(allocated(self%patch))deallocate(self%patch)
!   if(self%dust_on)call self%mydust%clean_memory()
! end subroutine usr_cloud_clean_memory
! !--------------------------------------------------------------------
! !> Subtroutine for power law distribution,here the dust sizes are defined. Ndust bins, with all bins having equal total mass.
!
! subroutine usr_mat_powerlaw_withX(n_point,power_a,min_var,max_var,   &
!                                   var_r)
!  implicit none
!  integer, intent(in)    :: n_point
!  real(dp), intent(in)   :: power_a
!  real(dp), intent(in)   :: min_var,max_var
!  real(dp), intent(inout):: var_r(n_point)
!  ! .. local ..
!  integer                :: i
!  real(dp)               :: r(0:n_point)
!  !------------------------------------
!  if(dabs(power_a-1.0_dp)<smalldouble)then
!   call mpistop('power_a ==1 at usr_mat_powerlaw_withX in mod_usr.t')
!  end if
!  r(0) = min_var
!
!  Loop_point : do i = 1,n_point
!   r(i)      = (dsqrt(r(i-1))+(dsqrt(max_var)-dsqrt(min_var))/n_point)**2.0_dp
!   !dvar_r(i) = r(i)-r(i-1)
!   var_r(i)  = -1.0_dp/(power_a-1.0_dp)*&
!               (r(i)**(power_a+2.0_dp)-r(i-1)**(power_a+2.0_dp))/&
!               (r(i)**(power_a+1)-r(i-1)**(power_a+1))
!  end do Loop_point
! end subroutine usr_mat_powerlaw_withX
!
! !--------------------------------------------------------------------
! !> subroutien to set profile de distance r
! subroutine usr_mat_profile(ixI^L,ixO^L,typeaxial_loc,profile, &
!                           center,extend,x,fprofile)
!  implicit none
!  integer,intent(in)           :: ixI^L,ixO^L
!  character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
!  character(len=*), intent(in) :: profile
!  real(kind=dp), intent(in)    :: center(1:ndim),extend(1:ndim)
!  real(kind=dp), intent(in)    :: x(ixI^S,1:ndim)
!  real(kind=dp), intent(inout) :: fprofile(ixI^S)
!  !.. local ..
!  real(dp)                     :: Dist(ixI^S)
!  !--------------------------------
!  call usr_distance(ixI^L,ixO^L,typeaxial_loc,center,x,dist)
!  select case(trim(profile))
!    case('gaussian')
!      ! where(Dist(ixO^S)<extend(1)*0.2)
!      !   fprofile(ixO^S) = 1.0_dp
!      ! else where
!      !21.04d5 pour extend(1)
!      !21.025d5
!     fprofile(ixO^S) = 21.025d5*dexp(-Dist(ixO^S)**2.0/(2.0*(extend(1))**2.0_dp))&
!                    / (2.0*dpi*(extend(1))**2.0_dp)
!      ! end where
!   case('tanh')
!     fprofile(ixO^S) = tanh(x(ixO^S,r_)/(3*extend(1)))
!    case default
!     fprofile(ixO^S) =1.0_dp
!  end select
! end subroutine usr_mat_profile
! !--------------------------------------------------------------------
!   subroutine usr_setsun()
!    sun%unit       = 'cgs'
!    sun%mass       = 1.9892d+33  !> solar mass cgs
!    sun%luminosity = 3.826d+33   !> solar luminosity (erg/s)
!    sun%radius     = 6.95987d+10 !> solar radius (cm)
!   end subroutine usr_setsun
! !--------------------------------------------------------------------
!   subroutine usr_physical_unit()
!
!    unit_user%volum          = unit_length**3.0_dp
!    unit_user%mass           = unit_density*unit_user%volum
!    unit_user%luminosity     = unit_length**5.0_dp*unit_density &
!                               *unit_time**(-3.0_dp)
!
!
!   end subroutine usr_physical_unit
! !------------------------------------------------------------------------
! !-----------------------------------------------------------------------
!
!
!
! subroutine usr_set_patch_sphere(ixI^L,ixO^L,typeaxial_loc,center,extend,x,patch)
!   implicit none
!   integer, intent(in)          :: ixI^L,ixO^L
!   character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
!   real(dp), intent(in)         :: x(ixI^S,1:ndim)
!   real(dp), intent(in)         :: center(1:ndim),extend(1:ndim)
!   logical, intent(inout)       :: patch(ixI^S)
!   !.. local ..
!   real(dp)                     :: Dist(ixI^S)
!   !-----------------------------------------------
!   call usr_distance(ixI^L,ixO^L,typeaxial_loc,center,x,dist)
!   patch(ixO^S)=(Dist(ixO^S)<=extend(r_))
! end  subroutine usr_set_patch_sphere
! !-------------------------------------------------------------------------
!
! subroutine usr_set_patch_cylinder(ixI^L,ixO^L,typeaxial_loc,center,extend,x,patch)
!   implicit none
!   integer, intent(in)          :: ixI^L,ixO^L
!   character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
!   real(dp), intent(in)         :: x(ixI^S,1:ndim)
!   real(dp), intent(in)         :: center(1:ndim),extend(1:ndim)
!   logical, intent(inout)       :: patch(ixI^S)
!   !.. local ..
!   real(dp)                     :: Dist(ixI^S)
!   !-----------------------------------------------
!   call usr_distance(ixI^L,ixO^L,typeaxial_loc,center,x,dist)
!   patch(ixO^S)=(Dist(ixI^S)<=extend(r_))
! end  subroutine usr_set_patch_cylinder
! !-------------------------------------------------------------------------
!
! subroutine usr_set_patch_cube(ixI^L,ixO^L,typeaxial_loc,center,extend,x,patch)
!   implicit none
!   integer, intent(in)          :: ixI^L,ixO^L
!   character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
!   real(dp), intent(in)         :: x(ixI^S,1:ndim)
!   real(dp), intent(in)         :: center(1:ndim),extend(1:ndim)
!   logical, intent(inout)       :: patch(ixI^S)
!   !.. local ..
!   integer                      :: idims
!   real(dp)                     :: Dist(ixI^S)
!   !-----------------------------------------------
!   select case(trim(typeaxial_loc))
!   case('slab')
!     patch(ixO^S)=dabs(x(ixO^S,1)-center(1))<extend(1)
!     if(ndim>=2)then
!       Loop_idim : do idims=2,ndim
!         patch(ixO^S)=patch(ixO^S).and.dabs(x(ixO^S,idims)-center(idims))<extend(idims)
!       end do Loop_idim
!     end if
!   case('cylindrical')
!     call mpistop('not implimented')
!   case('spherical')
!     call mpistop('not implimented')
!   end select
! end  subroutine usr_set_patch_cube
! !-------------------------------------------------------------------------
!   !> Distance between a cells and point 'center'
!   subroutine usr_distance(ixI^L,ixO^L,typeaxial_loc,center,x,dist)
!     implicit none
!     integer, intent(in)          :: ixI^L,ixO^L
!     character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
!     real(dp), intent(in)         :: x(ixI^S,1:ndim)
!     real(dp), intent(in)         :: center(1:ndim)
!     real(dp), intent(inout)      :: Dist(ixI^S)
!     !.. local ..
!     integer                      :: idims
!     real(dp)                     :: XX(ixI^S,1:ndim)
!     !------------------------------------------------
!   select case(trim(typeaxial_loc))
!   case('slab')
!     FORALL (idims=1:ndim)XX(ixO^S,idims) = x(ixO^S,idims)-center(idims)
!     Dist(ixO^S) = dsqrt(SUM(xx(ixO^S,1:ndim)**2.0_DP,dim=ndim+1))
!   case('cylindrical')
!     if(phi_<=ndim)then
!       Dist(ixO^S) =   x(ixO^S,r_)**2.0_dp+center(r_)**2.0_dp - &
!                       2.0_dp*x(ixO^S,r_)*center(r_)*dcos(x(ixO^S,phi_)*center(phi_))
!     else
!       Dist(ixO^S) =   (x(ixO^S,r_)-center(r_))**2.0_dp
!     end if
!     if(z_<=ndim) then
!       Dist(ixO^S) =Dist(ixO^S) +(x(ixO^S,z_)-center(z_))**2.0_dp
!     end if
!
!     Dist(ixO^S) = dsqrt(Dist(ixO^S))
!   case('spherical')
!     if(ndim == 1)    then
!       Dist(ixO^S) =   (x(ixO^S,r_)-center(r_))**2.0_dp
!     else if(ndim==2) then
!       if(phi_<=ndim)then
!         Dist(ixO^S) =   x(ixO^S,r_)**2.0_dp+center(r_)**2.0_dp - &
!                         2.0_dp*x(ixO^S,r_)*center(r_)&
!                         *dcos(x(ixO^S,phi_)*center(phi_))
!       else if(z_<=ndim) then
!         Dist(ixO^S) =   x(ixO^S,r_)**2.0_dp+center(r_)**2.0_dp - &
!                         2.0_dp*x(ixO^S,r_)*center(r_)&
!                         *dcos(x(ixO^S,z_)*center(z_))
!       end if
!     else if(ndim==3) then
!       Dist(ixO^S) = x(ixO^S,r_)**2.0_DP+center(r_)**2.0_DP&
!                    -2.0_DP*x(ixO^S,r_)*center(r_)*&
!                    ( dsin(x(ixO^S,z_))*dsin(center(z_))&
!                    *dcos(x(ixO^S,phi_)-center(phi_))&
!                     +dcos(x(ixO^S,z_))*dcos(center(z_)) )
!     end if
!   case default
!   end select
! end subroutine usr_distance
! !------------------------------------------------------------------------
! !> subroutine get the volume
! subroutine usr_get_volume(extend,shape, volume)
!  implicit none
!  real(dp), intent(in)         :: extend(:)
!  character(len=*), intent(in) :: shape
!  real(dp), intent(inout)      :: volume
!  !---------------------------------------------------
!   select case(shape)
!   case('sphere')
!    volume = 4.0_dp/3.0_dp*dpi*extend(1)**3.0_dp
!  case ('cube')
!     select case(ndim)
!      case(1)
!       volume =extend(1)
!      case(2)
!       volume =(extend(1)*extend(2))
!      case(3)
!       volume =extend(1)*extend(2)*extend(3)
!    end select
!  case('cylinder')
!    select case(ndim)
!      case(1)
!       volume = dpi*extend(r_)**2.0_dp
!     case(2)
!       if(phi_<=ndim)then
!         volume = dpi*extend(r_)**2.0_dp
!       elseif(z_<=ndim) then
!         volume = dpi*extend(r_)**2.0_dp *extend(z_)
!       end if
!     case(3)
!       volume = dpi*extend(r_)**2.0_dp *extend(z_)
!    end select
!  end select
! end subroutine usr_get_volume
! !----------------------------------------------------------------------
! !> compute the total mass and volume in the cloud
! subroutine usr_global_var
!   use mod_global_parameters
!   integer       :: igrid,iigrid
!   real(kind=dp) :: all_mass
!   !-----------------------
!   if(it==0)then
!     SUM_MASS   = 0.0_dp
!     SUM_Volume = 0.0_dp
!    do iigrid=1,igridstail; igrid=igrids(iigrid);
!      saveigrid=igrid
!      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
!      if(allocated(cloud_ary%patch))deallocate(cloud_ary%patch)
!      call  usr_cloud_patch(ixG^LL,ixM^LL,pw(igrid)%x,cloud_ary)
!      if(any(cloud_ary%patch(ixM^T)))then
!        SUM_MASS = SUM_MASS+sum(pw(igrid)%w(ixM^T,rho_)*2.0*dpi*&
!                   pw(igrid)%x(ixM^T,1)*(dxlevel(1)*dxlevel(2))&
!                   ,mask=(cloud_ary%patch(ixM^T).and.pw(igrid)%x(ixM^T,1)>0.0_dp))&
!                   *unit_density*unit_length**3.0_dp
!        SUM_Volume= SUM_Volume+ sum(2.0*dpi*pw(igrid)%x(ixM^T,1)*(dxlevel(1)*dxlevel(2))&
!                     ,mask=(cloud_ary%patch(ixM^T).and.pw(igrid)%x(ixM^T,1)>0.0_dp))
!
!      end if
!    end do
!    if(npe>1) then
!      call MPI_Reduce(SUM_MASS,all_mass,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
!    else
!      all_mass=SUM_MASS
!    end if
!    if(mype==0)   PRINT*,' code mass/phys mass :',SUM_MASS/cloud_ary%mass, &
!     'the code volume/phys volume : ',SUM_Volume/(4.0_dp*dpi/3.0_dp*cloud_ary%extend(1)**3.0_dp),&
!     cloud_ary%density*SUM_Volume/(cloud_ary%mass/(unit_density*unit_length**3.0_dp))
!   end if
! end subroutine usr_global_var
!
! !------------------------------------------------------------------------------------------
!   !> subroutine to compute large variation in of an array
!   subroutine usr_meanvalue_of_array(ixI^L,ixO^L,array,mean_array,mean_flag)
!     use mod_global_parameters
!     use mod_small_values
!     implicit none
!     integer, intent(in)             :: ixI^L, ixO^L
!     real(kind=dp), intent(in)       :: array(ixI^S)
!     real(kind=dp), intent(out)      :: mean_array(ixI^S)
!     logical, intent(out)            :: mean_flag(ixI^S)
!     ! .. local ..
!     integer                         :: ix^D,kxO^L,n_cell_daverage,n_cell_true
!     real(kind=dp)                   :: sum_array,usr_max_deviation
!     !----------------------------------------------------------------
!     usr_max_deviation = 100.0d0
!     n_cell_daverage=2
!     mean_flag(ixI^S)  = .true.
!
!     {do ix^DB= ixO^LIM^DB\}
!
!       {kxOmin^D= max(ix^D-n_cell_daverage, ixImin^D);
!       kxOmax^D= min(ix^D+n_cell_daverage, ixImax^D);\}
!
!       SUM_array=sum(array(kxO^S),mean_flag(kxO^S))
!
!       n_cell_true = count(mean_flag(kxO^S))
!
!       if(n_cell_true>1) then
!         mean_array(ix^D)=(sum_array-array(ix^D))&
!            /(count(mean_flag(kxO^S))-1.0_dp)
!       else
!         mean_array(ix^D)= array(ix^D)
!       end if
!
!       if(dabs(array(ix^D))>smalldouble&
!          .and.dabs(mean_array(ix^D))>smalldouble)then
!
!           mean_flag(ix^D)=dabs((array(ix^D)-mean_array(ix^D))/mean_array(ix^D))&
!                            <usr_max_deviation
!       else
!           mean_flag(ix^D)=.true.
!       end if
!
!     {end do^D&\}
!
!   end subroutine usr_meanvalue_of_array
!   !------------------------------------------------------------------------------------------
!
!
!   !------------------------------------------------------------------------------------------
!   !> subroutine to compute large variation in of an array
!   subroutine usr_medianvalue_of_array(ixI^L,ixO^L,array,median_array,indices_median,&
!                                        mean_flag)
!     use mod_global_parameters
!     use mod_small_values
!     implicit none
!     integer, intent(in)             :: ixI^L, ixO^L
!     real(kind=dp), intent(in)       :: array(ixI^S)
!     real(kind=dp), intent(out)      :: median_array(ixI^S)
!     logical, intent(out)            :: mean_flag(ixI^S)
!     integer, intent(out)            :: indices_median(ixI^S,1:ndim,2)
!     ! .. local ..
!     integer                         :: ix^D,kxO^L,n_cell_daverage,n_cell_true
!     integer                         :: loc_indices(1:ndim,1:2)
!     real(kind=dp)                   :: median_value,usr_max_deviation
!     real(kind=dp)                   :: usr_min_deviation_abs,usr_max_variation_abs
!     logical                         :: logic_sort
!     !----------------------------------------------------------------
!     usr_max_deviation     = 2.0d0
!     usr_min_deviation_abs = 1.0e-4
!     usr_max_variation_abs = 2.0d-2
!     n_cell_daverage       = 1
!     mean_flag(ixI^S)      = .true.
!
!
!
!
!     {do ix^DB= ixO^LIM^DB\}
!
!       {kxOmin^D= max(ix^D-n_cell_daverage, ixImin^D);
!       kxOmax^D= min(ix^D+n_cell_daverage, ixImax^D);\}
!       if(dabs(array(ix^D))<usr_min_deviation_abs)then
!         median_array(ix^D)=array(ix^D)
!         mean_flag(ix^D)   =.true.
!         cycle
!       end if
!       if(maxval(array(kxO^S))-minval(array(kxO^S))<usr_max_variation_abs)then
!         median_array(ix^D)=array(ix^D)
!         mean_flag(ix^D)   =.true.
!         cycle
!       end if
!
!       call usr_median_from_w(ixI^L,kxO^L,array,median_value,&
!                              loc_indices,logic_sort)
!
!       if(logic_sort)then
!         indices_median(ix^D,1:ndim,1:2) =  loc_indices(1:ndim,1:2)
!          median_array(ix^D) = median_value
!
!        if(dabs(array(ix^D))>dabs(median_array(ix^D))&
!          .and.dabs(median_value)>smalldouble)then
!
!           mean_flag(ix^D)=dabs((array(ix^D)-median_value)/median_value)&
!                            <usr_max_deviation
!
!        else if(dabs(median_value)<smalldouble) then
!           mean_flag(ix^D)=dabs(array(ix^D))<100.0_dp
!        else
!           mean_flag(ix^D)=.true.
!        end if
!       else
!         median_array(ix^D) =array(ix^D)
!       end if
!     {end do^D&\}
!
! end subroutine usr_medianvalue_of_array
!   !------------------------------------------------------------------------------------------
!   subroutine usr_median_from_w(ixI^L,ixO^L,w_array,median_value,indice_median,logic_sort)
!     use mod_global_parameters
!     use mod_small_values
!     implicit none
!     integer, intent(in)             :: ixI^L, ixO^L
!     real(kind=dp), intent(in)       :: w_array(ixI^S)
!     real(kind=dp), intent(out)      :: median_value
!     integer, intent(out)            :: indice_median(1:ndim,2)
!     logical, intent(out)            :: logic_sort
!     ! .. local ..
!     integer                         :: n_point, n_half,n_R,n_colonne
!     integer                         :: n_cells(ndim),ix_half(ndim)
!     real(kind=dp)                   :: array({^D&(ixOmax^D-ixOmin^D+1)*})
!     integer                         :: array_ind({^D&(ixOmax^D-ixOmin^D+1)*})
!     !----------------------------------------------------
!      n_point={^D&(ixOmax^D-ixOmin^D+1)*}
!      array=reshape(w_array(ixO^S),(/n_point/))
!      call usr_median_array(n_point,array,array_ind,median_value,logic_sort)
!      if(.not.logic_sort)return
!
!      {^D&n_cells(^D)=ixOmax^D-ixOmin^D+1;}
!      n_half    = n_point/2
!
!      n_colonne = array_ind(n_half)/n_cells(1)
!      if(mod(array_ind(n_half),n_cells(1))==0) then
!       indice_median(1,1) = n_cells(1)
!       indice_median(2,1) = n_colonne
!      else
!       indice_median(1,1) = array_ind(n_half)-n_colonne*n_cells(1)
!       indice_median(2,1) = n_colonne +1
!      end if
!
!      if(mod(n_point,2)==0) then
!       n_R       = n_half+1
!       n_colonne = array_ind(n_R)/n_cells(1)
!       if(mod(array_ind(n_R),n_cells(1))==0) then
!        indice_median(1,2) = n_cells(1)
!        indice_median(2,2) = n_colonne
!       else
!        indice_median(1,1) = array_ind(n_R)-n_colonne*n_cells(1)
!        indice_median(2,1) = n_colonne +1
!       end if
!      else
!       indice_median(:,2) = indice_median(:,1)
!      end if
!   end subroutine usr_median_from_w
!
!   !------------------------------------------------------------------------------------------
!   subroutine usr_median_array(n_point,array,array_ind,median_value,logic_sort)
!     implicit none
!     integer, intent(in)             :: n_point
!     real(kind=dp), intent(inout)    :: array(n_point)
!     real(kind=dp), intent(out)      :: median_value
!     integer, intent(out)            :: array_ind(n_point)
!     logical, intent(out)            :: logic_sort
!     ! .. local ..
!
!     integer                         :: i
!     !----------------------------------------------------
!     array_ind=(/(i,i=1,n_point)/)
!     call  usr_sort_array(n_point,array,array_ind,logic_sort)
!     if(mod(n_point,2)==0)then
!       median_value  =(array(n_point/2+1)+array(n_point/2-1))/2.0_dp
!     else
!       median_value=array(n_point/2+1)
!     end if
!
!   end subroutine usr_median_array
!   !-------------------------------------------------
!   subroutine usr_sort_array(n_point,array,array_ind,logic_sort)
!
!     implicit none
!     integer, intent(in)             :: n_point
!     real(kind=dp), intent(inout)    :: array(n_point)
!     integer, intent(out)            :: array_ind(n_point)
!     logical, intent(out)            :: logic_sort
!
!     ! .. local ..
!     real(kind=dp)                   :: temp_data
!     integer                         :: i, iter,temp_indice
!     !----------------------------------------------------
!     if(maxval(array(:))-minval(array(:))<smalldouble)then
!      logic_sort=.false.
!     else
!      logic_sort=.true.
!      iter=0
!     Loop_while : do while(iter<n_point-2)
!       Loop_sort : do i=1,n_point-1
!        if(array(i)>array(i+1))then
!         temp_data = array(i)
!         array(i)  = array(i+1)
!         array(i+1)=temp_data
!         temp_indice    = array_ind(i)
!         array_ind(i)   = array_ind(i+1)
!         array_ind(i+1) = temp_indice
!         iter = 0
!        else
!         iter=iter+1
!        end if
!       end do Loop_sort
!      end do Loop_while
!     end if
!   end subroutine usr_sort_array
end module mod_usr
