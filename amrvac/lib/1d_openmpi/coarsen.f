!=============================================================================
subroutine coarsen_grid_siblings(igrid,ipe,child_igrid,child_ipe,active)

use mod_global_parameters

integer, intent(in) :: igrid, ipe
integer, dimension(2), intent(in) :: child_igrid, child_ipe
logical, intent(in) :: active

integer :: igridFi, ipeFi, ixComin1,ixComax1, ixCoGmin1,ixCoGmax1, ixCoMmin1,&
   ixCoMmax1, ic1

if (ipe==mype) call alloc_node(igrid)

! New passive cell, coarsen from initial condition:
if (.not. active) then

   if (ipe == mype) then
      call initial_condition(igrid)
      
      do ic1=1,2
      igridFi=child_igrid(ic1)
      ipeFi=child_ipe(ic1)
      !if (ipeFi==mype) then
      !   ! remove solution space of child      
      !   call dealloc_node(igridFi)
      !end if
      end do
      
   end if

   return
end if



do ic1=1,2
   igridFi=child_igrid(ic1)
   ipeFi=child_ipe(ic1)

   if (ipeFi==mype) then
      dxlevel(1)=rnode(rpdx1_,igridFi);
      if (ipe==mype) then
         ixComin1=ixMlo1+(ic1-1)*(ixMhi1-ixMlo1+1)/2;
         ixComax1=ixMhi1+(ic1-2)*(ixMhi1-ixMlo1+1)/2;

         call coarsen_grid(pw(igridFi)%w,pw(igridFi)%x,ixGlo1,ixGhi1,ixMlo1,&
            ixMhi1,pw(igrid)%w,pw(igrid)%x,ixGlo1,ixGhi1, ixComin1,ixComax1,&
            igridFi,igrid)

         ! remove solution space of child
         !call dealloc_node(igridFi)
      else
         ixCoGmin1=1;
         ixCoGmax1=ixGhi1/2+nghostcells;
         ixCoMmin1=ixCoGmin1+nghostcells;ixCoMmax1=ixCoGmax1-nghostcells;
         call coarsen_grid(pw(igridFi)%w,pw(igridFi)%x,ixGlo1,ixGhi1,ixMlo1,&
            ixMhi1,pw(igridFi)%wcoarse,pw(igridFi)%xcoarse, ixCoGmin1,&
            ixCoGmax1,ixCoMmin1,ixCoMmax1,igridFi,igridFi)

         itag=ipeFi*max_blocks+igridFi
         isend=isend+1
         call MPI_ISEND(pw(igridFi)%wcoarse,1,type_coarse_block,ipe,itag,&
             icomm,sendrequest(isend),ierrmpi)
      end if
   else
      if (ipe==mype) then
         itag=ipeFi*max_blocks+igridFi
         irecv=irecv+1
         call MPI_IRECV(pw(igrid)%w,1,type_sub_block(ic1),ipeFi,itag, icomm,&
            recvrequest(irecv),ierrmpi)
      end if
   end if
end do

end subroutine coarsen_grid_siblings
!=============================================================================
subroutine coarsen_grid(wFi,xFi,ixFiGmin1,ixFiGmax1,ixFimin1,ixFimax1,wCo,xCo,&
   ixCoGmin1,ixCoGmax1,ixComin1,ixComax1,igridFi,igridCo)

use mod_global_parameters
use mod_physics, only: phys_to_primitive, phys_to_conserved

integer, intent(in) :: ixFiGmin1,ixFiGmax1, ixFimin1,ixFimax1, ixCoGmin1,&
   ixCoGmax1, ixComin1,ixComax1, igridFi, igridCo
double precision, intent(inout) :: wFi(ixFiGmin1:ixFiGmax1,1:nw),&
    xFi(ixFiGmin1:ixFiGmax1,1:ndim)
double precision,intent(inout) :: wCo(ixCoGmin1:ixCoGmax1,1:nw),&
    xCo(ixCoGmin1:ixCoGmax1,1:ndim)

integer :: ixCo1, ixFi1, iw
double precision :: CoFiratio
!-----------------------------------------------------------------------------
! coarsen by 2 in every direction - conservatively

if(coarsenprimitive) call phys_to_primitive(ixFiGmin1,ixFiGmax1,ixFimin1,&
   ixFimax1,wFi,xFi)

if (slab) then
   CoFiratio=one/dble(2**ndim)
   do iw=1,nw
      do ixCo1 = ixComin1,ixComax1
         ixFi1=2*(ixCo1-ixComin1)+ixFimin1
         wCo(ixCo1,iw)=sum(wFi(ixFi1:ixFi1+1,iw))*CoFiratio
      end do
   end do
else
   if(igridFi==igridCo) then
     do iw=1,nw
        do ixCo1 = ixComin1,ixComax1
           ixFi1=2*(ixCo1-ixComin1)+ixFimin1
           wCo(ixCo1,iw)= sum(pw(igridFi)%dvolume(ixFi1:ixFi1+&
              1)*wFi(ixFi1:ixFi1+1,iw)) /pw(igridCo)%dvolumecoarse(ixCo1)
!        if(dabs(sum(pw(igridFi)%dvolume(ixFi^D:ixFi^D+1)) &
!               -pw(igridCo)%dvolumecoarse(ixCo^D))>smalldouble) then
!    print *,'mismatching volumes: Co',ixCo^D, 'for fine', ixFi^D
!    print *,' fine volumes:' ,sum(pw(igridFi)%dvolume(ixFi^D:ixFi^D+1))
!    print *,' coarse volume:' ,pw(igridCo)%dvolumecoarse(ixCo^D)
!    call mpistop("mismatch volume")
!   endif
        end do
     end do
   else
     do iw=1,nw
        do ixCo1 = ixComin1,ixComax1
           ixFi1=2*(ixCo1-ixComin1)+ixFimin1
           wCo(ixCo1,iw)= sum(pw(igridFi)%dvolume(ixFi1:ixFi1+&
              1)*wFi(ixFi1:ixFi1+1,iw)) /pw(igridCo)%dvolume(ixCo1)
!        if(dabs(sum(pw(igridFi)%dvolume(ixFi^D:ixFi^D+1)) &
!               -pw(igridCo)%dvolume(ixCo^D))>smalldouble) then
!    print *,'mismatching volumes: Co',ixCo^D, 'for fine', ixFi^D
!    print *,' fine volumes:' ,sum(pw(igridFi)%dvolume(ixFi^D:ixFi^D+1))
!    print *,'      volume:' ,pw(igridCo)%dvolume(ixCo^D)
!    call mpistop("mismatch volume")
!   endif
        end do
     end do
   end if
end if

if(coarsenprimitive) then
  call phys_to_conserved(ixFiGmin1,ixFiGmax1,ixFimin1,ixFimax1,wFi,xFi)
  call phys_to_conserved(ixCoGmin1,ixCoGmax1,ixComin1,ixComax1,wCo,xCo)
end if

end subroutine coarsen_grid
!=============================================================================
