!> Build up AMR
subroutine settree
  use mod_global_parameters
  use mod_advance, only: advance

  ! create and initialize grids on all levels > 1. On entry, all
  ! level=1 grids have been formed and initialized. This subroutine
  ! creates and initializes new level grids

  integer :: igrid, iigrid, levnew

  ! when only one level allowed, there is nothing to do anymore
  if (refine_max_level == 1) return

  do levnew=2,refine_max_level

     ! advanced needed for refinement based on comparing w_n, w_n-1
     if(refine_criterion==1) then
       call setdt
       call advance(0)
     end if

     call errest

     if(refine_criterion==1) then
       do iigrid=1,igridstail; igrid=igrids(iigrid);
         pw(igrid)%w1 = pw(igrid)%wold
         pw(igrid)%wold = pw(igrid)%w
         pw(igrid)%w = pw(igrid)%w1
       end do
     end if

     call amr_coarsen_refine

     if (.not.reset_grid) then
       ! if no finer level grids created: exit
       if (levmax/=levnew) exit
     end if

  end do

end subroutine settree

!> reset AMR and (de)allocate boundary flux storage at level changes
subroutine resettree
  use mod_global_parameters
  use mod_fix_conserve

  if (levmax>levmin) call deallocateBflux

  call errest

  call amr_coarsen_refine

  ! set up boundary flux conservation arrays
  if (levmax>levmin) call allocateBflux

end subroutine resettree

!> Force the tree to desired level(s) from level_io(_min/_max)
!> used for conversion to vtk output
subroutine resettree_convert
  use mod_global_parameters
  use mod_ghostcells_update
  integer  :: igrid,iigrid, my_levmin, my_levmax

  if(level_io > 0) then
    my_levmin = level_io
    my_levmax = level_io
  else
    my_levmin = max(1,level_io_min)
    my_levmax = min(refine_max_level,level_io_max)
  end if

  do while(levmin<my_levmin.or.levmax>my_levmax)
   call getbc(global_time,0.d0,0,nwflux+nwaux)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call forcedrefine_grid_io(igrid,pw(igrid)%w)
   end do

   call amr_coarsen_refine
  end do

end subroutine resettree_convert
