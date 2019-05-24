!> Module containing the physics routines for scalar advection
module mod_rho_phys
  use mod_global_parameters, only: std_len
  use mod_constants
  use mod_physics
  implicit none
  private

  integer, protected, public          :: rho_       = 1
  double precision, protected, public :: rho_v(1) = 1.0d0

  ! Public methods
  public :: rho_phys_init
  public :: rho_get_v

  logical, allocatable,target                    :: rho_iw_average(:)
  type(physconfig),target,public                 :: rho_config
  type(phys_variables_indices),target,public     :: rho_ind
contains

  subroutine rho_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /rho_list/ rho_v

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, rho_list, end=111)
111    close(unitpar)
    end do

  end subroutine rho_params_read

  !> Write this module's parameters to a snapsoht
  subroutine rho_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er
    integer                             :: idim

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    do idim=1,ndim
      write(names(idim),'(a,i1)') "v",idim
      values(idim) = rho_v(idim)
    end do
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine rho_write_info

  subroutine rho_phys_init()
    use mod_global_parameters
    use mod_physics

    call rho_params_read(par_files)

    physics_type              = "rho"
    phys_energy               = .false.
    phys_req_diagonal         = .false.
    rho_config%dust_on        =.false.
    rho_config%dust_n_species = 0
    rho_config%ismhd          = .false.
    rho_config%isrel          = .false.
    rho_config%He_abundance   = 0.0_dp
    !rho_ = var_set_rho()
    call hd_fill_phys_indices
    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    phys_get_cmax        => rho_get_cmax
    phys_get_cbounds     => rho_get_cbounds
    phys_get_flux        => rho_get_flux
    phys_get_v_idim      => rho_get_v_idim
    phys_add_source_geom => rho_add_source_geom
    phys_to_conserved    => rho_to_conserved
    phys_to_primitive    => rho_to_primitive
    phys_get_dt          => rho_get_dt
    phys_write_info      => rho_write_info

    rho_config%energy=.false.
    phys_config          => rho_config
    phys_ind             => rho_ind
  end subroutine rho_phys_init

  subroutine hd_fill_phys_indices
    use mod_global_parameters
    implicit none
    ! .. local ..
    integer :: itr, idir
    !-------------------------------------

    ! Determine flux variables
    rho_ = var_set_rho()

    rho_ind%rho_  =rho_
    rho_ind%e_    =-1
    rho_ind%pressure_   =-1
    rho_ind%lfac_ =-1
    rho_ind%xi_   =-1
    rho_ind%psi_  =-1
  end subroutine hd_fill_phys_indices
  subroutine rho_to_conserved(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:1)

    ! Do nothing (primitive and conservative are equal for rho module)
  end subroutine rho_to_conserved

  subroutine rho_to_primitive(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:1)

    ! Do nothing (primitive and conservative are equal for rho module)
  end subroutine rho_to_primitive

  subroutine rho_get_v(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, idim, v,&
      centered)
    use mod_global_parameters
    logical, intent(in)           :: centered
    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    double precision, intent(in)  :: w(ixImin1:ixImax1, nw), x(ixImin1:ixImax1,&
        1:1)
    double precision, intent(out) :: v(ixImin1:ixImax1)

    double precision :: dtheta, dphi, halfdtheta, halfdphi, invdtheta, invdphi
    

    select case (typeaxial)
    case ("cylindrical")
       
       call mpistop("advection in 1D cylindrical not available")
      
       
       
    case ("spherical")
       
       call mpistop("advection in 1D spherical not available")
      
       
       
    case default
       v(ixOmin1:ixOmax1) = rho_v(idim)
    end select
  end subroutine rho_get_v

  !> Calculate simple v component
  subroutine rho_get_v_idim(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    double precision, intent(in)  :: w(ixImin1:ixImax1,nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1)

    v(ixOmin1:ixOmax1) = rho_v(idim)

  end subroutine rho_get_v_idim

  subroutine rho_get_cmax(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, idim, cmax)
    use mod_global_parameters
    integer, intent(in)                       :: ixImin1,ixImax1, ixOmin1,&
       ixOmax1, idim
    double precision, intent(in)              :: w(ixImin1:ixImax1, nw),&
        x(ixImin1:ixImax1, 1:1)
    double precision, intent(inout)           :: cmax(ixImin1:ixImax1)

    call rho_get_v(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, idim, cmax, .true.)

    cmax(ixOmin1:ixOmax1) = abs(cmax(ixOmin1:ixOmax1))

  end subroutine rho_get_cmax

  subroutine rho_get_cbounds(wLC, wRC, wLp, wRp, x, ixImin1,ixImax1, ixOmin1,&
     ixOmax1, idim, cmax, cmin)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    double precision, intent(in)    :: wLC(ixImin1:ixImax1, nw),&
        wRC(ixImin1:ixImax1,nw)
    double precision, intent(in)    :: wLp(ixImin1:ixImax1, nw),&
        wRp(ixImin1:ixImax1,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:1)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1)

    ! If get_v depends on w, the first argument should be some average over the
    ! left and right state
    call rho_get_v(wLC, x, ixImin1,ixImax1, ixOmin1,ixOmax1, idim, cmax,&
        .false.)

    if (present(cmin)) then
       cmin(ixOmin1:ixOmax1) = min(cmax(ixOmin1:ixOmax1), zero)
       cmax(ixOmin1:ixOmax1) = max(cmax(ixOmin1:ixOmax1), zero)
    else
       cmax(ixOmin1:ixOmax1) = maxval(abs(cmax(ixOmin1:ixOmax1)))
    end if

  end subroutine rho_get_cbounds

  subroutine rho_get_dt(w, ixImin1,ixImax1, ixOmin1,ixOmax1, dtnew, dx1, x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: dx1, x(ixImin1:ixImax1, 1:1)
    double precision, intent(in)    :: w(ixImin1:ixImax1, 1:nw)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble
  end subroutine rho_get_dt

  ! There is nothing to add to the transport flux in the transport equation
  subroutine rho_get_flux(wC, w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, idim, f)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    double precision, intent(in)    :: wC(ixImin1:ixImax1, 1:nw)
    double precision, intent(in)    :: w(ixImin1:ixImax1, 1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:1)
    double precision, intent(out)   :: f(ixImin1:ixImax1, nwflux)
    double precision                :: v(ixImin1:ixImax1)

    call rho_get_v(wC, x, ixImin1,ixImax1, ixOmin1,ixOmax1, idim, v, .false.)

    f(ixOmin1:ixOmax1, rho_) = w(ixOmin1:ixOmax1, rho_) * v(ixOmin1:ixOmax1)
  end subroutine rho_get_flux

  subroutine rho_add_source_geom(qdt, ixImin1,ixImax1, ixOmin1,ixOmax1, wCT, w,&
      x)

    ! Add geometrical source terms to w
    ! There are no geometrical source terms in the transport equation

    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: qdt, x(ixImin1:ixImax1, 1:1)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1, 1:nw),&
        w(ixImin1:ixImax1, 1:nw)

  end subroutine rho_add_source_geom

end module mod_rho_phys
