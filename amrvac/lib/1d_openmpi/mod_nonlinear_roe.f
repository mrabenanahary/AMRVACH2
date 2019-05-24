!> Module containing Roe solver for scalar nonlinear equation
module mod_nonlinear_roe
  use mod_physics_roe
  use mod_nonlinear_phys

  implicit none
  private

  public :: nonlinear_roe_init

contains

  subroutine nonlinear_roe_init()
    use mod_physics_roe

    nworkroe = 1

    phys_average         => nonlinear_average
    phys_get_eigenjump   => nonlinear_get_eigenjump
    phys_rtimes          => nonlinear_rtimes
  end subroutine nonlinear_roe_init

  subroutine nonlinear_average(wL, wR, x, ixmin1,ixmax1, idim, wroe, workroe)
    use mod_global_parameters
    integer, intent(in)             :: ixmin1,ixmax1, idim
    double precision, intent(in)    :: wL(ixGlo1:ixGhi1, nw), wR(ixGlo1:ixGhi1,&
        nw)
    double precision, intent(inout) :: wroe(ixGlo1:ixGhi1, nw)
    double precision, intent(inout) :: workroe(ixGlo1:ixGhi1, nworkroe)
    double precision, intent(in)    :: x(ixGlo1:ixGhi1, 1:1)

    wroe(ixmin1:ixmax1, rho_)=half*(wL(ixmin1:ixmax1, rho_)+wR(ixmin1:ixmax1,&
        rho_))
  end subroutine nonlinear_average

  subroutine nonlinear_get_eigenjump(wL, wR, wC, x, ixmin1,ixmax1, il, idim,&
      smalla, a, jump, workroe)

    ! Calculate the characteristic speed a and the jump in the
    ! characteristic variable in the idim direction within ixL.
    ! For a scalar equation the characteristic and conservative variables coincide
    ! The characteristic speed is just the velocity, but it should be averaged
    ! for the cell interfaces

    use mod_global_parameters

    integer, intent(in)                          :: ixmin1,ixmax1, il, idim
    double precision, dimension(ixGlo1:ixGhi1, nw)       :: wL, wR, wC
    double precision, dimension(ixGlo1:ixGhi1)           :: smalla, a, jump, v
    double precision, dimension(ixGlo1:ixGhi1, nworkroe) :: workroe
    double precision, intent(in)                 :: x(ixGlo1:ixGhi1, 1:1)
    integer                                      :: jxmin1,jxmax1, ixCmin1,&
       ixCmax1

    jxmin1=ixmin1+kr(idim,1);jxmax1=ixmax1+kr(idim,1);
    ixCmin1=ixmin1; ixCmax1=jxmax1;

    ! No entropy fix
    smalla(ixmin1:ixmax1)= -one
    ! The velocity is dependent of w in the nonlinear scalar equation,
    ! and thus depends on the location
    !> TODO: check this, for advection added argument to get velocity at cell edge!!!
    call nonlinear_get_v(wL, x, ixGlo1,ixGhi1, ixCmin1,ixCmax1, idim, v)

    a(ixmin1:ixmax1)=(v(jxmin1:jxmax1)+v(ixmin1:ixmax1))/2

    jump(ixmin1:ixmax1)=wR(ixmin1:ixmax1, rho_)-wL(ixmin1:ixmax1, rho_)

  end subroutine nonlinear_get_eigenjump

  subroutine nonlinear_rtimes(q, w, ixmin1,ixmax1, iw, il, idim, rq, workroe)

    ! Multiply q by R(il, iw), where R is the right eigenvalue matrix at wC.
    ! For a scalar equation the R matrix is unity

    use mod_global_parameters
    integer, intent(in)             :: ixmin1,ixmax1, iw, il, idim
    double precision, intent(in)    :: w(ixGlo1:ixGhi1, nw), q(ixGlo1:ixGhi1)
    double precision, intent(inout) :: rq(ixGlo1:ixGhi1)
    double precision, intent(inout) :: workroe(ixGlo1:ixGhi1, nworkroe)

    rq(ixmin1:ixmax1)=q(ixmin1:ixmax1)

  end subroutine nonlinear_rtimes

end module mod_nonlinear_roe
