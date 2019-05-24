module mod_setdt

contains
!>setdt  - set dt for all levels between levmin and levmax.
!>         dtpar>0  --> use fixed dtpar for all level
!>         dtpar<=0 --> determine CFL limited timestep
subroutine setdt()
use mod_global_parameters
use mod_physics
use mod_usr_methods, only: usr_get_dt
use mod_thermal_conduction

integer :: iigrid, igrid, ncycle, ncycle2, ifile, idim
double precision :: dtnew, qdtnew, dtmin_mype, factor, dx1, dxmin1

double precision :: dtmax, dxmin, cmax_mype, v(ixGlo1:ixGhi1)
!----------------------------------------------------------------------------

if (dtpar<=zero) then
   dtmin_mype=bigdouble
   cmax_mype = zero
!$OMP PARALLEL DO PRIVATE(igrid,qdtnew,dtnew,dx1)
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      dtnew=bigdouble
      dx1=rnode(rpdx1_,igrid);
      dxlevel(1)=rnode(rpdx1_,igrid);
      saveigrid = igrid
      block=>pw(igrid)
      block%iw0=0

      if (nwaux>0) then
         call phys_get_aux(.true.,pw(igrid)%w,pw(igrid)%x,ixGlo1,ixGhi1,ixMlo1,&
            ixMhi1,'setdt')
      end if


      call getdt_courant(pw(igrid)%w,ixGlo1,ixGhi1,ixMlo1,ixMhi1,qdtnew,&
         pw(igrid)%x)
      dtnew=min(dtnew,qdtnew)
      if(dtnew<1e-8)then
        print*,' is at setdt, get_courant ',mype,it,igrid,dtnew
      end if

      call phys_get_dt(pw(igrid)%w,ixGlo1,ixGhi1,ixMlo1,ixMhi1,qdtnew,dx1,&
         pw(igrid)%x)
      dtnew=min(dtnew,qdtnew)

      if(dtnew<1e-8)then
        print*,' is at setdt, phys_get_dt ',mype,it,igrid,dtnew
      end if
      if (associated(usr_get_dt)) then
         call usr_get_dt(pw(igrid)%w,ixGlo1,ixGhi1,ixMlo1,ixMhi1,global_time,&
            qdtnew,dx1,pw(igrid)%x)
      end if

      dtnew          = min(dtnew,qdtnew)
      dtmin_mype     = min(dtmin_mype,dtnew)
      dt_grid(igrid) = dtnew
   end do
!$OMP END PARALLEL DO
else
   dtmin_mype=dtpar
end if

if (dtmin_mype<dtmin) then
   write(unitterm,*)"Warning: Time step too small!", dtmin_mype
   write(unitterm,*)"on processor:", mype
   write(unitterm,*)"at time:", global_time," step:", it
   call mpistop("too small timestep")
end if

if (slowsteps>it-it_init+1) then
   factor=one-(one-dble(it-it_init+1)/dble(slowsteps))**2
   dtmin_mype=dtmin_mype*factor
end if


dtmin_mype=min(dtmin_mype,time_max-global_time)

if (dtpar<=zero) then
   call MPI_ALLREDUCE(dtmin_mype,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN, icomm,&
      ierrmpi)
else
   dt=dtmin_mype
end if

if(any(dtsave(1:nfile)<bigdouble).or.any(tsave(isavet(1:nfile),&
   1:nfile)<bigdouble))then
   dtmax = minval(ceiling(global_time/dtsave(1:nfile))*dtsave(1:nfile))-&
      global_time
   do ifile=1,nfile
      dtmax = min(tsave(isavet(ifile),ifile)-global_time,dtmax)
   end do
   if(dtmax > smalldouble)then
     dt=min(dt,dtmax)
   else
     ! dtmax=0 means dtsave is divisible by global_time
     dt=min(dt,minval(dtsave(1:nfile)))
   end if
end if

if(mype==0) then
  if(any(dtsave(1:nfile)<dt)) then
    write(unitterm,*) 'Warning: timesteps: ',dt,' exceeding output intervals ',&
        dtsave(1:nfile)
  endif
endif

! estimate time step of thermal conduction
if(associated(phys_getdt_heatconduct)) then
   dtmin_mype=bigdouble
!$OMP PARALLEL DO PRIVATE(igrid,qdtnew,&
!$OMP& dx1) REDUCTION(min:dtmin_mype)
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      dx1=rnode(rpdx1_,igrid);
      saveigrid = igrid
      block=>pw(igrid)
      qdtnew=bigdouble
      call phys_getdt_heatconduct(pw(igrid)%w,ixGlo1,ixGhi1,ixMlo1,ixMhi1,&
         qdtnew,dx1,pw(igrid)%x)
      dtmin_mype=min(dtmin_mype,qdtnew)
   end do
!$OMP END PARALLEL DO
   call MPI_ALLREDUCE(dtmin_mype,dtnew,1,MPI_DOUBLE_PRECISION,MPI_MIN, icomm,&
      ierrmpi)
   if(all(flux_scheme=='nul')) dt=min(dt,dtnew)
   ncycle=ceiling(dt/dtnew)
   if (ncycle>tc_ncycles) then
     if(mype==0 .and. .false.) then
       write(*,*)&
           'CLF time step is too many times larger than conduction time step',&
          ncycle
       write(*,*) 'reducing dt to',tc_ncycles,'times of dt_impl!!'
     endif
     dt=tc_ncycles*dtnew
   endif
  ! get number of sub-steps of supertime stepping (Meyer 2012 MNRAS 422,2102)
   if(dt/dtnew< 0.5d0) then
     s=1
   else if(dt/dtnew< 2.d0) then
     s=2
   else
     s=ceiling((dsqrt(9.d0+8.d0*dt/dtnew)-1.d0)/2.d0)
     ! only use odd s number
     s=s/2*2+1
   endif
   dt_tc=dt*0.5d0
   if(mype==0 .and. .false.) write(*,*) 'supertime steps:',s,&
      ' normal subcycles:',ceiling(dt/dtnew/2.d0)
endif

!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
   dt_grid(igrid)=dt
end do
!$OMP END PARALLEL DO


! global Lax-Friedrich finite difference flux splitting needs fastest wave-speed
! so does GLM:
if(need_global_cmax) call MPI_ALLREDUCE(cmax_mype,cmax_global,1,&
   MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)
! some scheme need maximal speed of flow
if(need_global_vmax) then
  cmax_mype=0.d0
  do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
     do idim=1,ndim
       call phys_get_v_idim(pw(igrid)%w,pw(igrid)%x,ixGlo1,ixGhi1,ixMlo1,&
          ixMhi1,idim,v)
       cmax_mype=max(cmax_mype,maxval(abs(v(ixMlo1:ixMhi1))))
     end do
  end do
  call MPI_ALLREDUCE(cmax_mype,vmax_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
     icomm,ierrmpi)
  vmax_global=cmax_global-vmax_global
end if

contains

  !> compute CFL limited dt (for variable time stepping)
  subroutine getdt_courant(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dtnew,x)

  use mod_global_parameters
  use mod_physics, only: phys_get_cmax

  integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
  double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
  double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw), dtnew

  integer          :: idims,iter_correct,patchierror(ixImin1:ixImax1)
  double precision :: courantmax, dxinv(1:ndim), courantmaxtot, courantmaxtots
  double precision :: cmax(ixImin1:ixImax1), cmaxtot(ixImin1:ixImax1),&
      tmp(ixImin1:ixImax1)
  logical          :: log_nocrrect
  !-----------------------------------------------------------------------------
  log_nocrrect=.true.
  Loop_cor : do iter_correct=1,2*ndim
    if(saveigrid==63.and.it==890)Print*,' i will try ',iter_correct
   dtnew=bigdouble

   courantmax=zero
   courantmaxtot=zero
   courantmaxtots=zero

   dxinv(1)=one/dx1;

   cmaxtot(ixOmin1:ixOmax1)=zero

   Loop_idims : do idims=1,ndim
     if(saveigrid==63.and.it==890)print*,' is test ss',idims
      call phys_get_cmax(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,cmax)
      if(need_global_cmax) cmax_mype = max(cmax_mype,&
         maxval(cmax(ixOmin1:ixOmax1)))
      if (.not.slab) then
         tmp(ixOmin1:ixOmax1)=cmax(ixOmin1:ixOmax1)/block%ds(ixOmin1:ixOmax1,&
            idims)
         cmaxtot(ixOmin1:ixOmax1)=cmaxtot(ixOmin1:ixOmax1)+&
            tmp(ixOmin1:ixOmax1)
         courantmax=max(courantmax,maxval(tmp(ixOmin1:ixOmax1)))
      else
         cmaxtot(ixOmin1:ixOmax1)=cmaxtot(ixOmin1:ixOmax1)+&
            cmax(ixOmin1:ixOmax1)*dxinv(idims)
         courantmax=max(courantmax,maxval(cmax(ixOmin1:ixOmax1)*dxinv(idims)))
      end if
      if (small_getdt_average) then
       where (cmaxtot(ixOmin1:ixOmax1)>small_dt_coef*courantpar/(dtmin))
        patchierror(ixOmin1:ixOmax1) = 1
       else where
        patchierror(ixOmin1:ixOmax1) = 0
       end where
       if (any(patchierror(ixOmin1:ixOmax1)/=0)) then
          log_nocrrect=.false.
          call setdt_small_dt_values_average(ixImin1,ixImax1, ixOmin1,ixOmax1,&
              phys_iw_average,"getdt_courant", w, x, patchierror)
          cycle  Loop_cor
       end if
      end if
      courantmaxtot=courantmaxtot+courantmax
      log_nocrrect=.true.
    end do Loop_idims

    if(saveigrid==63.and.it==890)then
      PRINT*,'test time step correction ',(courantpar/courantmax),iter_correct,&
         small_getdt_average
    end if

    if(log_nocrrect)exit Loop_cor
  end do Loop_cor
  select case (typecourant)
  case ('minimum')
     ! courantmax='max(c/dx)'
     if (courantmax>smalldouble)     dtnew=min(dtnew,courantpar/courantmax)
  case ('summax')
     ! courantmaxtot='summed max(c/dx)'
     if (courantmaxtot>smalldouble)  dtnew=min(dtnew,courantpar/courantmaxtot)
  case ('maxsum')
     ! courantmaxtots='max(summed c/dx)'
     courantmaxtots=max(courantmaxtots,maxval(cmaxtot(ixOmin1:ixOmax1)))
     if (courantmaxtots>smalldouble) dtnew=min(dtnew,&
        courantpar/courantmaxtots)
  case default
     write(unitterm,*)'Unknown typecourant=',typecourant
     call mpistop("Error from getdt_courant: no such typecourant!")
  end select

  end subroutine getdt_courant

end subroutine setdt




!=============================================================================
! subroutine setdt_detect_smalldt(igrid,dt_step,subname)
! use mod_error_info
!
!
! integer, intent(in)            :: igrid
! double precision, intent(in)   :: dt_step
! character(len=*), intent(in)   :: subname
! !--------------------------------------------------
! cond_dtmin : if (dt_step<dtmin) then
!   call error_stepdt_stop_info(igrid,dt_step,subname)
! end if cond_dtmin
! end subroutine setdt_detect_smalldt
!=============================================================================
subroutine setdt_small_dt_values_average(ixImin1,ixImax1, ixOmin1,ixOmax1,&
    iw_average,subname, w, x, w_flag)
  use mod_global_parameters
  use mod_physics
  use mod_small_values
  integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
  logical, intent(in)             :: iw_average(1:nw)
  integer, intent(in)             :: w_flag(ixImin1:ixImax1)
  character(len=*), intent(in)    :: subname
  double precision, intent(inout) :: w(ixImin1:ixImax1, 1:nw)
  double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)
  integer                         :: iw, kxOmin1,kxOmax1, ix1, i

  do ix1= ixOmin1,ixOmax1

  ! point with local failure identified by w_flag
  if (w_flag(ix1) /= 0) then
    ! verify in cube with border width small_values_daverage the presence of
    ! cells where all went ok
    do i = 1, max(small_values_daverage, 1)
      kxOmin1= max(ix1-i, ixOmin1);
      kxOmax1= min(ix1+i, ixOmax1);

      ! in case cells are fine within smaller cube than
      ! the userset small_values_daverage: use that smaller cube
      if (any(w_flag(kxOmin1:kxOmax1) == 0)) exit
    end do

    if (any(w_flag(kxOmin1:kxOmax1) == 0)) then
      ! within surrounding cube, cells without problem were found

      ! faulty cells are corrected by averaging here
      ! only average those which were ok and replace faulty cells
      call phys_to_primitive(ixImin1,ixImax1,kxOmin1,kxOmax1,w,x)
      do iw = 1, nw
        if(iw_average(iw))w(ix1, iw) = sum(w(kxOmin1:kxOmax1, iw),&
            w_flag(kxOmin1:kxOmax1) == 0)/ count(w_flag(kxOmin1:kxOmax1) == 0)
      end do
      call phys_to_conserved(ixImin1,ixImax1,kxOmin1,kxOmax1,w,x)
    else
      write(*,*) "no cells without error were found in cube of size",&
          small_values_daverage
      write(*,*) "at location:", x(ix1, 1:ndim)*unit_length
      write(*,*) "at index:", ix1
      write(*,*) "w_flag(ix^D):", w_flag(ix1)
      write(*,*)  "w = ",w(ix1, 1:nw)
      write(*,*) 'iteration',it
      write(*,*) "Saving status at the previous time step"
      write(*,*) "called from : ",subname
      crash=.true.
      call MPISTOP('is STOP HERE ')
    end if
  end if
  enddo

end subroutine setdt_small_dt_values_average
!=============================================================================
end module mod_setdt
