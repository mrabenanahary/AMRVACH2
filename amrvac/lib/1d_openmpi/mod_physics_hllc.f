module mod_physics_hllc

  implicit none
  public

  procedure(sub_diffuse_hllcd), pointer   :: phys_diffuse_hllcd => null()
  procedure(sub_get_lCD), pointer         :: phys_get_lCD => null()
  procedure(sub_get_wCD), pointer         :: phys_get_wCD => null()

  abstract interface
     subroutine sub_diffuse_hllcd(ixImin1,ixImax1,ixOmin1,ixOmax1,idims,wLC,&
        wRC,fLC,fRC,patchf)
       use mod_global_parameters
       integer, intent(in)                                    :: ixImin1,&
          ixImax1,ixOmin1,ixOmax1,idims
       double precision, dimension(ixImin1:ixImax1,1:nw), intent(in)    :: wRC,&
          wLC
       double precision, dimension(ixImin1:ixImax1,1:nwflux),intent(in) :: fLC,&
           fRC
       integer, dimension(ixImin1:ixImax1),&
           intent(inout)               :: patchf
     end subroutine sub_diffuse_hllcd

     subroutine sub_get_lCD(wLC,wRC,fLC,fRC,cmin,cmax,idims,ixImin1,ixImax1,&
        ixOmin1,ixOmax1, whll,Fhll,lambdaCD,patchf)
       use mod_global_parameters
       integer, intent(in)                                      :: ixImin1,&
          ixImax1,ixOmin1,ixOmax1,idims
       double precision, dimension(ixImin1:ixImax1,1:nw),&
           intent(in)      :: wLC,wRC
       double precision, dimension(ixImin1:ixImax1,1:nwflux),&
           intent(in)  :: fLC,fRC
       double precision, dimension(ixImin1:ixImax1),&
           intent(in)           :: cmax,cmin
       integer, dimension(ixImin1:ixImax1),&
           intent(inout)                 :: patchf
       double precision, dimension(ixImin1:ixImax1,1:nwflux),&
           intent(out) :: Fhll,whll
       double precision, dimension(ixImin1:ixImax1),&
           intent(out)          :: lambdaCD
     end subroutine sub_get_lCD

     subroutine sub_get_wCD(x,wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,cmin,&
        cmax,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,f)
       use mod_global_parameters
       integer, intent(in)                                     :: ixImin1,&
          ixImax1,ixOmin1,ixOmax1,idims
       double precision, dimension(ixImin1:ixImax1,1:nw),&
           intent(in)     :: wRC,wLC
       double precision, dimension(ixImin1:ixImax1,1:ndim), intent(in)    :: x
       double precision, dimension(ixImin1:ixImax1,1:nwflux),&
           intent(in) :: whll, Fhll
       double precision, dimension(ixImin1:ixImax1),&
           intent(in)          :: lambdaCD
       double precision, dimension(ixImin1:ixImax1),&
           intent(in)          :: cmax,cmin
       double precision, dimension(ixImin1:ixImax1,1:nwflux),&
           intent(in) :: fRC,fLC
       integer, dimension(ixImin1:ixImax1),&
           intent(inout)                :: patchf
       double precision, dimension(ixImin1:ixImax1,1:nwflux),intent(out) :: f
     end subroutine sub_get_wCD
  end interface

contains

  subroutine phys_hllc_check
    if (.not. associated(phys_diffuse_hllcd)) phys_diffuse_hllcd => &
       dummy_diffuse_hllcd

    if (.not. associated(phys_get_lCD)) phys_get_lCD => dummy_get_lCD

    if (.not. associated(phys_get_wCD)) phys_get_wCD => dummy_get_wCD
  end subroutine phys_hllc_check

  ! When method is hllcd or hllcd1 then: this subroutine is to impose enforce
  ! regions where we AVOID HLLC and use TVDLF instead: this is achieved by setting
  ! patchf to 4 in certain regions. An additional input parameter is nxdiffusehllc
  ! which sets the size of the fallback region. This nul version enforces TVDLF
  ! everywhere!
  subroutine dummy_diffuse_hllcd(ixImin1,ixImax1,ixOmin1,ixOmax1,idims,wLC,wRC,&
     fLC,fRC,patchf)
    use mod_global_parameters
    integer, intent(in)                                    :: ixImin1,ixImax1,&
       ixOmin1,ixOmax1,idims
    double precision, dimension(ixImin1:ixImax1,1:nw), intent(in)    :: wRC,&
       wLC
    double precision, dimension(ixImin1:ixImax1,1:nwflux),intent(in) :: fLC,&
        fRC
    integer, dimension(ixImin1:ixImax1), intent(inout)               :: patchf

    patchf(ixOmin1:ixOmax1) = 4 ! enforce TVDLF everywhere
  end subroutine dummy_diffuse_hllcd

  ! Calculate lambda at CD and set the patchf to know the orientation
  ! of the riemann fan and decide on the flux choice
  ! We also compute here the HLL flux and w value, for fallback strategy
  ! In this nul version, we simply compute nothing and ensure TVDLF fallback
  subroutine dummy_get_lCD(wLC,wRC,fLC,fRC,cmin,cmax,idims,ixImin1,ixImax1,&
     ixOmin1,ixOmax1, whll,Fhll,lambdaCD,patchf)
    use mod_global_parameters
    integer, intent(in)                                      :: ixImin1,&
       ixImax1,ixOmin1,ixOmax1,idims
    double precision, dimension(ixImin1:ixImax1,1:nw), intent(in)      :: wLC,&
       wRC
    double precision, dimension(ixImin1:ixImax1,1:nwflux), intent(in)  :: fLC,&
       fRC
    double precision, dimension(ixImin1:ixImax1), intent(in)           :: cmax,&
       cmin
    integer, dimension(ixImin1:ixImax1), intent(inout)                 :: &
       patchf
    double precision, dimension(ixImin1:ixImax1,1:nwflux), intent(out) :: Fhll,&
       whll
    double precision, dimension(ixImin1:ixImax1),&
        intent(out)          :: lambdaCD

    ! Next must normally be computed
    Fhll(ixOmin1:ixOmax1,1:nwflux) = zero
    whll(ixOmin1:ixOmax1,1:nwflux) = zero
    lambdaCD(ixOmin1:ixOmax1)      = zero

    ! This actually ensures fallback to TVDLF
    patchf(ixOmin1:ixOmax1)=4
  end subroutine dummy_get_lCD

  ! compute the intermediate state U*. Only needed where patchf=-1/1. This nul
  ! version simply nullifies all values
  subroutine dummy_get_wCD(x,wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,cmin,&
     cmax,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,f)
    use mod_global_parameters
    integer, intent(in)                                     :: ixImin1,ixImax1,&
       ixOmin1,ixOmax1,idims
    double precision, dimension(ixImin1:ixImax1,1:nw), intent(in)     :: wRC,&
       wLC
    double precision, dimension(ixImin1:ixImax1,1:ndim), intent(in)    :: x
    double precision, dimension(ixImin1:ixImax1,1:nwflux), intent(in) :: whll,&
        Fhll
    double precision, dimension(ixImin1:ixImax1),&
        intent(in)          :: lambdaCD
    double precision, dimension(ixImin1:ixImax1), intent(in)          :: cmax,&
       cmin
    double precision, dimension(ixImin1:ixImax1,1:nwflux), intent(in) :: fRC,&
       fLC
    integer, dimension(ixImin1:ixImax1), intent(inout)                :: &
       patchf
    double precision, dimension(ixImin1:ixImax1,1:nwflux),intent(out) :: f

    ! Next must normally be computed
    f(ixOmin1:ixOmax1,1:nwflux)  = zero
  end subroutine dummy_get_wCD

end module mod_physics_hllc
