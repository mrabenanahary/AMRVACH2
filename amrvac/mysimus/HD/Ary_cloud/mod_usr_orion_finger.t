module mod_usr
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
  integer, parameter  :: n_dust_max = 20
  real(dp) :: SUM_MASS   = 0.0_dp
  real(dp) :: SUM_VOLUME = 0.0_dp




  type (ISM)           :: ism_orion
  type (cloud),target  :: cloud_ary
  type (dust),target   :: dust_ary






contains
  subroutine usr_init
    use mod_hd
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
    call hd_activate
    call usr_physical_unit!

    call usr_check_confilt
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
    if(ism_on)call ism_orion%read_parameters(ism_orion%myconfig,files)

   if(cloud_on)call cloud_ary%read_parameters(cloud_ary%myconfig,files)
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
    use mod_hd
    implicit none
    cond_dust_on : if(.not.hd_dust)then
      ism_orion%myconfig%dust_on =.false.
      cloud_ary%myconfig%dust_on =.false.
    end if  cond_dust_on
  end   subroutine usr_check_confilt
  !-----------------------------------------------------------
  !> subroutine to normalize parameters in the code
  subroutine usr_normalise_parameters
    use mod_hd
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
  subroutine initglobaldata_usr
   use mod_variables


   call ism_orion%set_complet
   call ism_orion%normalize


   call cloud_ary%set_complet
   call cloud_ary%normalize


   call usr_normalise_parameters
   if(mype==0)call usr_write_setting
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
    type(dust)    :: dust_dummy

    !-----------------------------------------

    if(first)then
      if(mype==0) then
        write(*,*)'cloud propagation for ary :-)'
      endif
      first=.false.
    endif

    itr=1

    call cloud_ary%alloc_set_patch(ixI^L,ixO^L,0.0_dp,x)
    call ism_orion%alloc_set_patch(ixI^L,ixO^L,0.0_dp,x,&
                                   escape_patch=cloud_ary%patch)


    ! set one cloud
    call cloud_ary%set_w(ixI^L,ixO^L,0.0_dp,x,w)

    ! set the ism
    call ism_orion%set_w(ixI^L,ixO^L,0.0_dp,x,w)




    ! put dust to zero in all other zones
    cond_dust_on : if(hd_dust) then

      call dust_dummy%set_allpatch(ixI^L,ixO^L,&
                  (/cloud_ary%mydust,ism_orion%mydust/))

      call dust_dummy%set_w_zero(ixI^L,ixO^L,x,w)
      call dust_dummy%clean_memory

      ! allocate(dust_ary%patch(ixI^S))
      ! dust_ary%patch(ixO^S) =.false.
      ! if(ism_orion%myconfig%dust_on)then
      !  ism_orion%mydust%patch(ixO^S)=.not.cloud_ary%patch(ixO^S)
      !  dust_ary%patch(ixO^S) = dust_ary%patch(ixO^S).or.ism_orion%mydust%patch(ixO^S)
      ! end if
      ! if(cloud_ary%myconfig%dust_on)then
      !  dust_ary%patch(ixO^S) = (dust_ary%patch(ixO^S).or. &
      !                          cloud_ary%mydust%patch(ixO^S))
      ! end if
      ! call dust_ary%set_w_zero(ixI^L,ixO^L,x,w)
      ! call dust_ary%clean_memory
    end if cond_dust_on

    ! check is if initial setting is correct
    call  phys_check_w(.true., ixI^L, ixO^L, w, flag)
    if(any(flag(ixO^S)/=0)) PRINT*,' is error',maxloc(abs(flag(ixO^S)))


    ! get conserved variables to be used in the code
    call phys_to_conserved(ixI^L,ixO^L,w,x)
    call ism_orion%clean_memory
    call cloud_ary%clean_memory
    
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

    small_dust_rho = ism_orion%mydust%myconfig%min_limit_rel
    coef_correct   = 0.4_dp

    call phys_to_primitive(ixI^L,ixO^L,w,x)

    Loop_idust : do idust =1, dust_n_species
      where(w(ixO^S, dust_rho(idust))<max(small_dust_rho*w(ixO^S,rho_),&
           ism_orion%mydust%myconfig%min_limit_abs))
          w(ixO^S, dust_rho(idust))= 0.8* min(small_dust_rho*w(ixO^S,rho_),&
               ism_orion%mydust%myconfig%min_limit_abs)
          patch_correct(ixO^S) = .true.
      elsewhere
          patch_correct(ixO^S) = .false.
      end where





     Loop_idir2 : do idir = 1,ndir
      where(patch_correct(ixO^S))
              w(ixO^S, dust_mom(idir,idust))=0.0_dp
      elsewhere
       where(w(ixO^S, dust_rho(idust))>ism_orion%mydust%myconfig%min_limit_abs)


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
         call cloud_ary%clean_memory
       else cond_first_step
       ! coarsen in the back of the moving cloud
       if(all(w(ixO^S,rho_)<min(small_density*100.0_dp,ism_orion%myconfig%density/100.0_dp)&
          .or.all(w(ixO^S,rho_)<min(small_pressure*100.0_dp,ism_orion%myconfig%pressure/100.0_dp))))then
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
    use mod_hd
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
     small_dust_rho = ism_orion%mydust%myconfig%min_limit_rel

     call phys_to_primitive(ixI^L,ixI^L,w,x)
     ! handel small density dust
     Loop_idust : do idust =1, dust_n_species
      where(w(ixI^S, dust_rho(idust))<max(small_dust_rho*w(ixI^S,rho_),&
         ism_orion%mydust%myconfig%min_limit_abs))
        w(ixI^S, dust_rho(idust))= 0.8* min(small_dust_rho*w(ixI^S,rho_),&
             ism_orion%mydust%myconfig%min_limit_abs)
        patch_correct(ixI^S) = .true.
      elsewhere
        patch_correct(ixI^S) = .false.
      end where
      ! handel large density dust
      where(w(ixI^S,rho_)<0.9*ism_orion%myconfig%density)
       where(w(ixI^S, dust_rho(idust))>ism_orion%mydust%myconfig%max_limit_rel*w(ixI^S,rho_))
        w(ixI^S, dust_rho(idust))=0.8*ism_orion%mydust%myconfig%max_limit_rel*w(ixI^S,rho_)
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
  subroutine usr_write_setting
    integer,parameter   :: unit_config =12
    character(len=50)   :: filename_config
    !-------------------------------------
    filename_config=trim(base_filename)//'_config'
    open(unit_config,file=trim(filename_config), status='replace')
    call ism_orion%write_setting(unit_config)
    call cloud_ary%write_setting(unit_config)
    close(unit_config)
  end subroutine usr_write_setting
  !----------------------------------------------------------------------
  !> compute the total mass and volume in the cloud
  subroutine usr_global_var
    use mod_global_parameters
    integer       :: igrid,iigrid
    real(kind=dp) :: all_mass
    !-----------------------
    if(it==0)then
      SUM_MASS   = 0.0_dp
      SUM_Volume = 0.0_dp
     do iigrid=1,igridstail; igrid=igrids(iigrid);
       saveigrid=igrid
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       if(allocated(cloud_ary%patch))deallocate(cloud_ary%patch)
       call  usr_cloud_patch(ixG^LL,ixM^LL,pw(igrid)%x,cloud_ary)
       if(any(cloud_ary%patch(ixM^T)))then
         SUM_MASS = SUM_MASS+sum(pw(igrid)%w(ixM^T,rho_)*2.0*dpi*&
                    pw(igrid)%x(ixM^T,1)*(dxlevel(1)*dxlevel(2))&
                    ,mask=(cloud_ary%patch(ixM^T).and.pw(igrid)%x(ixM^T,1)>0.0_dp))&
                    *unit_density*unit_length**3.0_dp
         SUM_Volume= SUM_Volume+ sum(2.0*dpi*pw(igrid)%x(ixM^T,1)*(dxlevel(1)*dxlevel(2))&
                      ,mask=(cloud_ary%patch(ixM^T).and.pw(igrid)%x(ixM^T,1)>0.0_dp))

       end if
     end do
     if(npe>1) then
       call MPI_Reduce(SUM_MASS,all_mass,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
     else
       all_mass=SUM_MASS
     end if
     if(mype==0)   PRINT*,' code mass/phys mass :',SUM_MASS/cloud_ary%myconfig%mass, &
      'the code volume/phys volume : ',SUM_Volume/(4.0_dp*dpi/3.0_dp*cloud_ary%myconfig%extend(1)**3.0_dp),&
      cloud_ary%myconfig%density*SUM_Volume/(cloud_ary%myconfig%mass/(unit_density*unit_length**3.0_dp))
    end if
  end subroutine usr_global_var

!---------------------------------------------------------------------
end module mod_usr
