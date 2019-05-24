!> Module containing the physics routines for scalar nonlinear equation
module mod_nonlinear_phys

  implicit none
  private

  !> index of the single scalar unknown
  integer, protected, public :: rho_       = 1

  !> switch between burgers (i.e. rho**2) 
  !> or nonconvex flux (i.e. rho**3)
  integer, protected, public :: nonlinear_flux_type = 1

  !> whether the KdV source term is added
  logical, protected, public :: kdv_source_term = .false.

  ! Public methods
  public :: nonlinear_phys_init
  public :: nonlinear_get_v

contains

  !> Read this module's parameters from a file
  subroutine nonlinear_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /nonlinear_list/ nonlinear_flux_type, kdv_source_term

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, nonlinear_list, end=111)
111    close(unitpar)
    end do

  end subroutine nonlinear_params_read

  !> Write this module's parameters to a snapshot
  subroutine nonlinear_write_info(fh)
  ! for nonlinear scalar equation, nothing to write
  ! note: this is info only stored at end of dat files, 
  !       is never read/used for restarts, only expects
  !       an integer (number of parameters) and 
  !       corresponding double values and character names
  !       and is meant for use in the python tools
    use mod_global_parameters
    integer, intent(in)                 :: fh



  end subroutine nonlinear_write_info

  subroutine nonlinear_phys_init()
    use mod_global_parameters
    use mod_physics
    use mod_kdv, only: kdv_init


    call nonlinear_params_read(par_files)

    physics_type = "nonlinear"
    phys_energy  = .false.
    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

    rho_ = var_set_rho()

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    phys_get_cmax        => nonlinear_get_cmax
    phys_get_cbounds     => nonlinear_get_cbounds
    phys_get_flux        => nonlinear_get_flux
    phys_add_source_geom => nonlinear_add_source_geom
    phys_add_source      => nonlinear_add_source
    phys_to_conserved    => nonlinear_to_conserved
    phys_to_primitive    => nonlinear_to_primitive
    phys_get_dt          => nonlinear_get_dt
    phys_write_info      => nonlinear_write_info

    if (kdv_source_term) call kdv_init()

  end subroutine nonlinear_phys_init

  subroutine nonlinear_to_conserved(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:1)

    ! Do nothing (primitive and conservative are equal for nonlinear module)
  end subroutine nonlinear_to_conserved

  subroutine nonlinear_to_primitive(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:1)

    ! Do nothing (primitive and conservative are equal for nonlinear module)
  end subroutine nonlinear_to_primitive

  subroutine nonlinear_get_v(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, idim, v)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    double precision, intent(in)  :: w(ixImin1:ixImax1, nw), x(ixImin1:ixImax1,&
        1:1)
    double precision, intent(out) :: v(ixImin1:ixImax1)

   select case(nonlinear_flux_type)
    case(1)
       v(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,rho_)
    case(2)
       v(ixOmin1:ixOmax1)=3.0d0*w(ixOmin1:ixOmax1,rho_)**2
    case default
       call mpistop('Undefined fluxtype: set nonlinear_flux_type to 1 or 2')
    end select

  end subroutine nonlinear_get_v

  subroutine nonlinear_get_cmax(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, idim,&
      cmax)
    use mod_global_parameters
    integer, intent(in)                       :: ixImin1,ixImax1, ixOmin1,&
       ixOmax1, idim
    double precision, intent(in)              :: w(ixImin1:ixImax1, nw),&
        x(ixImin1:ixImax1, 1:1)
    double precision, intent(inout)           :: cmax(ixImin1:ixImax1)

    call nonlinear_get_v(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, idim, cmax)

    cmax(ixOmin1:ixOmax1) = abs(cmax(ixOmin1:ixOmax1))

  end subroutine nonlinear_get_cmax

  subroutine nonlinear_get_cbounds(wLC, wRC, wLp, wRp, x, ixImin1,ixImax1,&
      ixOmin1,ixOmax1, idim, cmax, cmin)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    double precision, intent(in)    :: wLC(ixImin1:ixImax1, nw),&
        wRC(ixImin1:ixImax1,nw)
    double precision, intent(in)    :: wLp(ixImin1:ixImax1, nw),&
        wRp(ixImin1:ixImax1,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:1)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1)

    double precision :: wmean(ixImin1:ixImax1,nw)

    ! since get_v depends on w, the first argument should be some average over the
    ! left and right state
    wmean(ixOmin1:ixOmax1,1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,&
       1:nwflux)+wRC(ixOmin1:ixOmax1,1:nwflux))
    call nonlinear_get_v(wmean, x, ixImin1,ixImax1, ixOmin1,ixOmax1, idim,&
        cmax)

    if (present(cmin)) then
       cmin(ixOmin1:ixOmax1) = min(cmax(ixOmin1:ixOmax1), zero)
       cmax(ixOmin1:ixOmax1) = max(cmax(ixOmin1:ixOmax1), zero)
    else
       cmax(ixOmin1:ixOmax1) = maxval(abs(cmax(ixOmin1:ixOmax1)))
    end if

  end subroutine nonlinear_get_cbounds

  subroutine nonlinear_get_dt(w, ixImin1,ixImax1, ixOmin1,ixOmax1, dtnew, dx1,&
      x)
    use mod_global_parameters
    use mod_kdv, only: kdv_get_dt

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: dx1, x(ixImin1:ixImax1, 1:1)
    double precision, intent(in)    :: w(ixImin1:ixImax1, 1:nw)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble

    if(kdv_source_term) then
      call kdv_get_dt(w, ixImin1,ixImax1, ixOmin1,ixOmax1, dtnew, dx1, x)
    endif
  end subroutine nonlinear_get_dt

  ! here we select the flux according to the nonlinear_flux_type parameter
  subroutine nonlinear_get_flux(wC, w, x, ixImin1,ixImax1, ixOmin1,ixOmax1,&
      idim, f)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    double precision, intent(in)    :: wC(ixImin1:ixImax1, 1:nw)
    double precision, intent(in)    :: w(ixImin1:ixImax1, 1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:1)
    double precision, intent(out)   :: f(ixImin1:ixImax1, nwflux)

    select case(nonlinear_flux_type)
    case(1)
       f(ixOmin1:ixOmax1,rho_)=half*w(ixOmin1:ixOmax1,rho_)**2
    case(2)
       f(ixOmin1:ixOmax1,rho_)=w(ixOmin1:ixOmax1,rho_)**3
    case default
       call mpistop('Undefined fluxtype: set nonlinear_flux_type to 1 or 2')
    end select

  end subroutine nonlinear_get_flux

  subroutine nonlinear_add_source_geom(qdt, ixImin1,ixImax1, ixOmin1,ixOmax1,&
      wCT, w, x)

    ! Add geometrical source terms to w
    ! There are no geometrical source terms 

    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: qdt, x(ixImin1:ixImax1, 1:1)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1, 1:nw),&
        w(ixImin1:ixImax1, 1:nw)

  end subroutine nonlinear_add_source_geom

  subroutine nonlinear_add_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x,&
     qsourcesplit,active)
    use mod_global_parameters
    use mod_kdv, only: kdv_add_source

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1, 1:nw),&
        x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1, 1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

    if(kdv_source_term) then
      call kdv_add_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x,&
         qsourcesplit,active)
    end if

  end subroutine nonlinear_add_source

end module mod_nonlinear_phys
