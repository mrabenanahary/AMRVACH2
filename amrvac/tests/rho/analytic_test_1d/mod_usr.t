module mod_usr
  use mod_rho

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer :: i_err
  integer :: i_sol

contains

  subroutine usr_init()
    use mod_variables

    usr_init_one_grid => initonegrid_usr
    usr_process_grid => set_error
    usr_print_log => print_error

    call rho_activate()

    i_err = var_set_extravar("error", "error")
    i_sol = var_set_extravar("solution", "solution")
  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    w(ix^S,rho_) = solution(x(ix^S, 1), 0.0_dp)

  end subroutine initonegrid_usr

  elemental function solution(x, t) result(val)
    real(dp), intent(in) :: x, t
    real(dp)             :: val, xrel
    real(dp), parameter  :: pi = acos(-1.0_dp)

    select case(iprob)
       case default
          val = sin(2 * pi * (x - rho_v(1) * t))
       case (2)
          val = sin(2 * pi * (x - rho_v(1) * t))**2
       case (3)
          xrel = x - rho_v(1) * t
          xrel = modulo(xrel, xprobmax1)

          if (xrel > 0.1_dp .and. xrel < 0.3_dp) then
             val = 1.0_dp
          else if (xrel > 0.7 .and. xrel < 0.9) then
             val = 1 - 10 * abs(xrel - 0.8_dp)
          else
             val = 0.0_dp
          end if
       case (4)
          val = sin(2 * pi * (x - rho_v(1) * t))**8
       end select
  end function solution

  subroutine print_error()
    use mod_input_output, only: get_volume_average, get_global_maxima
    double precision   :: modes(nw, 2), volume
    double precision   :: maxvals(nw)

    call get_global_maxima(maxvals)
    call get_volume_average(1, modes(:, 1), volume)
    call get_volume_average(2, modes(:, 2), volume)

    if (mype == 0) then
       write(*, "(A,4E14.6)") " CONV_TEST (t,err_1,err_2,err_inf):", global_time, &
            modes(i_err, 1), sqrt(modes(i_err, 2)), maxvals(i_err)
    end if
  end subroutine print_error

  subroutine set_error(igrid,level,ixI^L,ixO^L,qt,w,x)
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S,i_sol) = solution(x(ixO^S, 1), qt)
    w(ixO^S,i_err) = abs(w(ixO^S,rho_) - w(ixO^S,i_sol))
  end subroutine set_error

end module mod_usr

