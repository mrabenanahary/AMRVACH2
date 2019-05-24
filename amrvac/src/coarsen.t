!=============================================================================
subroutine coarsen_grid_siblings(igrid,ipe,child_igrid,child_ipe,active)

use mod_global_parameters

integer, intent(in) :: igrid, ipe
integer, dimension(2^D&), intent(in) :: child_igrid, child_ipe
logical, intent(in) :: active

integer :: igridFi, ipeFi, ixCo^L, ixCoG^L, ixCoM^L, ic^D

if (ipe==mype) call alloc_node(igrid)

! New passive cell, coarsen from initial condition:
if (.not. active) then

   if (ipe == mype) then
      call initial_condition(igrid)
      
      {do ic^DB=1,2\}
      igridFi=child_igrid(ic^D)
      ipeFi=child_ipe(ic^D)
      !if (ipeFi==mype) then
      !   ! remove solution space of child      
      !   call dealloc_node(igridFi)
      !end if
      {end do\}
      
   end if

   return
end if



{do ic^DB=1,2\}
   igridFi=child_igrid(ic^D)
   ipeFi=child_ipe(ic^D)

   if (ipeFi==mype) then
      ^D&dxlevel(^D)=rnode(rpdx^D_,igridFi);
      if (ipe==mype) then
         ixComin^D=ixMlo^D+(ic^D-1)*(ixMhi^D-ixMlo^D+1)/2;
         ixComax^D=ixMhi^D+(ic^D-2)*(ixMhi^D-ixMlo^D+1)/2;

         call coarsen_grid(pw(igridFi)%w,pw(igridFi)%x,ixG^LL,ixM^LL,pw(igrid)%w,pw(igrid)%x,ixG^LL, &
                     ixCo^L,igridFi,igrid)

         ! remove solution space of child
         !call dealloc_node(igridFi)
      else
         ixCoGmin^D=1;
         ixCoGmax^D=ixGhi^D/2+nghostcells;
         ixCoM^L=ixCoG^L^LSUBnghostcells;
         call coarsen_grid(pw(igridFi)%w,pw(igridFi)%x,ixG^LL,ixM^LL,pw(igridFi)%wcoarse,pw(igridFi)%xcoarse, &
                           ixCoG^L,ixCoM^L,igridFi,igridFi)

         itag=ipeFi*max_blocks+igridFi
         isend=isend+1
         call MPI_ISEND(pw(igridFi)%wcoarse,1,type_coarse_block,ipe,itag, &
                        icomm,sendrequest(isend),ierrmpi)
      end if
   else
      if (ipe==mype) then
         itag=ipeFi*max_blocks+igridFi
         irecv=irecv+1
         call MPI_IRECV(pw(igrid)%w,1,type_sub_block(ic^D),ipeFi,itag, &
                        icomm,recvrequest(irecv),ierrmpi)
      end if
   end if
{end do\}

end subroutine coarsen_grid_siblings
!=============================================================================
subroutine coarsen_grid(wFi,xFi,ixFiG^L,ixFi^L,wCo,xCo,ixCoG^L,ixCo^L,&
                        igridFi,igridCo)

use mod_global_parameters
use mod_physics, only: phys_to_primitive, phys_to_conserved

integer, intent(in) :: ixFiG^L, ixFi^L, ixCoG^L, ixCo^L, igridFi, igridCo
double precision, intent(inout) :: wFi(ixFiG^S,1:nw), xFi(ixFiG^S,1:ndim)
double precision,intent(inout) :: wCo(ixCoG^S,1:nw), xCo(ixCoG^S,1:ndim)

integer :: ixCo^D, ixFi^D, iw
double precision :: CoFiratio
!-----------------------------------------------------------------------------
! coarsen by 2 in every direction - conservatively

if(coarsenprimitive) call phys_to_primitive(ixFiG^L,ixFi^L,wFi,xFi)

if (slab) then
   CoFiratio=one/dble(2**ndim)
   do iw=1,nw
      {do ixCo^DB = ixCo^LIM^DB
         ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixFimin^DB\}
         wCo(ixCo^D,iw)=sum(wFi(ixFi^D:ixFi^D+1,iw))*CoFiratio
      {end do\}
   end do
else
   if(igridFi==igridCo) then
     do iw=1,nw
        {do ixCo^DB = ixCo^LIM^DB
           ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixFimin^DB\}
           wCo(ixCo^D,iw)= &
               sum(pw(igridFi)%dvolume(ixFi^D:ixFi^D+1)*wFi(ixFi^D:ixFi^D+1,iw)) &
              /pw(igridCo)%dvolumecoarse(ixCo^D)
!        if(dabs(sum(pw(igridFi)%dvolume(ixFi^D:ixFi^D+1)) &
!               -pw(igridCo)%dvolumecoarse(ixCo^D))>smalldouble) then
!    print *,'mismatching volumes: Co',ixCo^D, 'for fine', ixFi^D
!    print *,' fine volumes:' ,sum(pw(igridFi)%dvolume(ixFi^D:ixFi^D+1))
!    print *,' coarse volume:' ,pw(igridCo)%dvolumecoarse(ixCo^D)
!    call mpistop("mismatch volume")
!   endif
        {end do\}
     end do
   else
     do iw=1,nw
        {do ixCo^DB = ixCo^LIM^DB
           ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixFimin^DB\}
           wCo(ixCo^D,iw)= &
               sum(pw(igridFi)%dvolume(ixFi^D:ixFi^D+1)*wFi(ixFi^D:ixFi^D+1,iw)) &
              /pw(igridCo)%dvolume(ixCo^D)
!        if(dabs(sum(pw(igridFi)%dvolume(ixFi^D:ixFi^D+1)) &
!               -pw(igridCo)%dvolume(ixCo^D))>smalldouble) then
!    print *,'mismatching volumes: Co',ixCo^D, 'for fine', ixFi^D
!    print *,' fine volumes:' ,sum(pw(igridFi)%dvolume(ixFi^D:ixFi^D+1))
!    print *,'      volume:' ,pw(igridCo)%dvolume(ixCo^D)
!    call mpistop("mismatch volume")
!   endif
        {end do\}
     end do
   end if
end if

if(coarsenprimitive) then
  call phys_to_conserved(ixFiG^L,ixFi^L,wFi,xFi)
  call phys_to_conserved(ixCoG^L,ixCo^L,wCo,xCo)
end if

end subroutine coarsen_grid
!=============================================================================
