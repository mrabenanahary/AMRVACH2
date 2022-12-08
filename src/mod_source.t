!> Module for handling split source terms (split from the fluxes)
module mod_source

  implicit none
  private

  public :: add_split_source
  public :: addsource2

contains

  subroutine add_split_source(prior)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_thermal_conduction, only: phys_thermal_conduction
    use mod_physics, only: phys_req_diagonal

    logical, intent(in) :: prior

    double precision :: qdt, qt
    integer :: iigrid, igrid, i^D
    logical :: src_active

    ! add thermal conduction
    if(associated(phys_thermal_conduction)) call phys_thermal_conduction()

    src_active = .false.

    if ((.not.prior).and.&
         (typesourcesplit=='sf' .or. typesourcesplit=='ssf')) return

    if (prior) then
       qt=global_time
    else
       qt=global_time+dt
    end if
    !$OMP PARALLEL DO PRIVATE(igrid,qdt,i^D)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
       qdt=dt_grid(igrid)
       block=>pw(igrid)
       call addsource1_grid(igrid,qdt,qt,pw(igrid)%w,src_active)
    end do
    !$OMP END PARALLEL DO

    if (src_active) then
       call getbc(qt,0.d0,0,nwflux+nwaux, phys_req_diagonal)
    end if

  end subroutine add_split_source

  subroutine addsource1_grid(igrid,qdt,qt,w,src_active)

    use mod_global_parameters

    integer, intent(in) :: igrid
    double precision, intent(in) :: qdt, qt
    double precision, intent(inout) :: w(ixG^T,nw)
    logical, intent(inout) :: src_active

    double precision :: w1(ixG^T,nw)

    saveigrid=igrid
    typelimiter=type_limiter(node(plevel_,igrid))
    typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

    ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

    w1(ixG^T,1:nw)=w(ixG^T,1:nw)

    select case (typesourcesplit)
    case ('sf')
       call addsource2(qdt  ,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,&
            pw(igrid)%x,.true.,src_active)
    case ('sfs')
       call addsource2(qdt/2,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,&
       pw(igrid)%x,.true.,src_active)
    case ('ssf')
       call addsource2(qdt/2,ixG^LL,ixG^LL,1,nw,qt,w,qt,w1,&
            pw(igrid)%x,.true.,src_active)
       call addsource2(qdt  ,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,&
            pw(igrid)%x,.true.,src_active)
    case ('ssfss')
       call addsource2(qdt/4,ixG^LL,ixG^LL,1,nw,qt,w,qt,w1,&
            pw(igrid)%x,.true.,src_active)
       call addsource2(qdt/2,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,&
            pw(igrid)%x,.true.,src_active)
    case default
       write(unitterm,*)'No such typesourcesplit=',typesourcesplit
       call mpistop("Error: Unknown typesourcesplit!")
    end select

  end subroutine addsource1_grid

  subroutine addsource2(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,&
       w,x,qsourcesplit,src_active)

    ! Add source within ixO for iws: w=w+qdt*S[wCT]

    use mod_global_parameters
    use mod_physics, only: phys_add_source
    use mod_usr_methods, only: usr_source
    ! differences with VAC is in iw^LIM and in declaration of ranges for wCT,w

    integer, intent(in)              :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)     :: qdt, qtC, qt
    double precision, intent(in)     :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout)  :: w(ixI^S,1:nw)
    logical, intent(in)              :: qsourcesplit
    logical, intent(inout), optional :: src_active
    logical                          :: tmp_active

    tmp_active = .false.
    ! user defined sources, typically explicitly added
    if ((qsourcesplit .eqv. source_split_usr) .and. associated(usr_source)) then
       tmp_active = .true.
       call usr_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    end if

    ! physics defined sources, typically explicitly added,
    ! along with geometrical source additions
    call phys_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,tmp_active)
    if (present(src_active)) src_active = src_active .or. tmp_active
  end subroutine addsource2

end module mod_source
