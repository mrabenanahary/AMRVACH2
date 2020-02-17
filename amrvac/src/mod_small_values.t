!> Module for handling problematic values in simulations, such as negative
!> pressures
module mod_small_values

  implicit none
  private

  !> How to handle small values
  character(len=20), public :: small_values_method      = "error"
  logical          , public :: small_values_force_floor =.false.
  !> Average over this many cells in each direction
  integer, public :: small_values_daverage = 1

  public :: small_values_error
  public :: small_values_average
  public :: small_var_average
contains


 subroutine small_var_average(ixI^L, ixO^L, pth, x, w_flag)
     use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    integer, intent(in)             :: w_flag(ixI^S)
    double precision, intent(inout) :: pth(ixI^S)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: iw, kxO^L, ix^D, i

    {do ix^DB= ixO^LIM^DB\}

    ! point with local failure identified by w_flag
    if (w_flag(ix^D) /= 0) then
      ! verify in cube with border width small_values_daverage the presence of
      ! cells where all went ok
      do i = 1, max(small_values_daverage, 1)
        {kxOmin^D= max(ix^D-i, ixOmin^D);
        kxOmax^D= min(ix^D+i, ixOmax^D);\}

        ! in case cells are fine within smaller cube than
        ! the userset small_values_daverage: use that smaller cube
        if (any(w_flag(kxO^S) == 0)) exit
      end do

      if (any(w_flag(kxO^S) == 0)) then
        ! within surrounding cube, cells without problem were found

        ! faulty cells are corrected by averaging here
        ! only average those which were ok and replace faulty cells

          pth(ix^D) = sum(pth(kxO^S), w_flag(kxO^S) == 0)&
               / count(w_flag(kxO^S) == 0)

      else
        write(*,*) "no cells without error were found in cube of size", &
             small_values_daverage
        write(*,*) "at location:", x(ix^D, 1:ndim)*unit_length
        write(*,*) "at index:", ix^D
        write(*,*) "w_flag(ix^D):", w_flag(ix^D)
        write(*,*)  "pth = ",pth(ix^D)
        write(*,*) 'iteration',it
        write(*,*) "Saving status at the previous time step"
        crash=.true.
      end if
    end if
    {enddo^D&\}

  end subroutine  small_var_average

  subroutine small_values_error(w, x, ixI^L, ixO^L, w_flag, subname)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    integer, intent(in)          :: w_flag(ixI^S)
    integer                      :: ix_bad(ndim)
    character(len=*), intent(in) :: subname

    ix_bad = maxloc(abs(w_flag(ixO^S))) + [ ixOmin^D-1 ]

    if (.not.crash) then
      write(*,*) "Error: small value of ", &
                 trim(prim_wnames(maxval(abs(w_flag(ixO^S))))), &
                 " encountered when call ", subname
      write(*,*) "Iteration: ", it, " Time: ", global_time
      write(*,*) "Location: ", x({ix_bad(^D)}, :)
      write(*,*) "Cell number: ", ix_bad(:)
      write(*,*) "w(1:nw): ", w({ix_bad(^D)}, 1:nw)
      write(*,*)"igrid = ", saveigrid, "cpu = ",mype
      write(*,*) "Saving status at the previous time step"
      crash=.true.
    end if
  end subroutine small_values_error

  subroutine small_values_average(ixI^L, ixO^L, subname, w, x, w_flag,ierror)
    use mod_global_parameters
    use mod_variables
    integer, intent(in)             :: ixI^L, ixO^L
    integer, intent(in)             :: w_flag(ixI^S)
    character(len=*), intent(in)    :: subname
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer, intent(out)            :: ierror
    integer                         :: iw, kxO^L, ix^D, i

    ierror = 0
    {do ix^DB= ixO^LIM^DB\}

    ! point with local failure identified by w_flag
    cond_need : if (w_flag(ix^D) > 0) then
      ! verify in cube with border width small_values_daverage the presence of
      ! cells where all went ok
      do i = 1, max(small_values_daverage, 1)
        {kxOmin^D= max(ix^D-i, ixOmin^D);
        kxOmax^D= min(ix^D+i, ixOmax^D);\}

        ! in case cells are fine within smaller cube than
        ! the userset small_values_daverage: use that smaller cube
        if (any(w_flag(kxO^S) <= 0)) exit
      end do

      cond_exist: if (any(w_flag(kxO^S) <= 0)) then
        ! within surrounding cube, cells without problem were found

        ! faulty cells are corrected by averaging here
        ! only average those which were ok and replace faulty cells
        Loop_averiage: do iw = 1, nw
          w(ix^D, iw)=sum(w(kxO^S, iw), w_flag(kxO^S) <= 0)&
               / count(w_flag(kxO^S) <= 0)
        end do Loop_averiage
      else cond_exist
        ierror = 1
        if(small_values_force_floor)return
        write(*,*) "no cells without error were found in cube of size", &
                    small_values_daverage
        write(*,*) "at location:", x(ix^D, 1:ndim)*unit_length
        write(*,*) "at index:", ix^D
        write(*,*) "w_flag(ix^D):", w_flag(ix^D)
        Loop_iw : do iw = 1,nw
         write(*,*)  prim_wnames(iw)," =", w(ix^D, iw)
        end do Loop_iw
        write(*,*)  "variable = ",prim_wnames(abs(w_flag(ix^D)))
        write(*,*) 'iteration',it
        write(*,*) "Saving status at the previous time step"
        write(*,*) "called from : ",subname
        crash=.true.
        call MPISTOP('is STOP HERE ')
      end if cond_exist
    end if cond_need
    {enddo^D&\}

  end subroutine small_values_average

end module mod_small_values
