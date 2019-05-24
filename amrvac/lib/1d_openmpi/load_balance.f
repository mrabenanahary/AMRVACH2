!=============================================================================
subroutine load_balance
use mod_forest
use mod_global_parameters

integer :: Morton_no, recv_igrid, recv_ipe, send_igrid, send_ipe, igrid, ipe

integer, external :: getnode

! Jannis: for now, not using version for passive/active blocks
call get_Morton_range()

if (npe==1) then
   sfc_to_igrid(:)=sfc(1,Morton_start(mype):Morton_stop(mype))
   return
end if

irecv=0
isend=0

do ipe=0,npe-1; do Morton_no=Morton_start(ipe),Morton_stop(ipe)
   recv_ipe=ipe

   send_igrid=sfc(1,Morton_no)
   send_ipe=sfc(2,Morton_no)

   if (recv_ipe/=send_ipe) then
      ! get an igrid number for the new node in recv_ipe processor
      recv_igrid=getnode(recv_ipe)
      ! update node igrid and ipe on the tree
      call change_ipe_tree_leaf(recv_igrid,recv_ipe,send_igrid,send_ipe)
      ! receive physical data of the new node in recv_ipe processor
      if (recv_ipe==mype) call lb_recv
      ! send physical data of the old node in send_ipe processor
      if (send_ipe==mype) call lb_send
   end if
   if (recv_ipe==mype) then
      if (recv_ipe==send_ipe) then
         sfc_to_igrid(Morton_no)=send_igrid
      else
         sfc_to_igrid(Morton_no)=recv_igrid
      end if
   end if
end do; end do

if (irecv>0) call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
if (isend>0) call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)

! post processing
do ipe=0,npe-1; do Morton_no=Morton_start(ipe),Morton_stop(ipe)
   recv_ipe=ipe

   send_igrid=sfc(1,Morton_no)
   send_ipe=sfc(2,Morton_no)

   if (recv_ipe/=send_ipe) then
      !if (send_ipe==mype) call dealloc_node(send_igrid)
      call putnode(send_igrid,send_ipe)
   end if
end do; end do


! Update sfc array: igrid and ipe info in space filling curve
call amr_Morton_order()

!!$if (nwaux>0) then
!!$   do Morton_no=Morton_start(mype),Morton_stop(mype)
!!$      if (sfc(2,Morton_no)/=mype) then
!!$         igrid=sfc_to_igrid(Morton_no)
!!$         saveigrid=igrid
!!$         call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL,"load_balance")
!!$      end if
!!$   end do
!!$end if

contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine lb_recv
!-----------------------------------------------------------------------------
call alloc_node(recv_igrid)

itag=recv_igrid
irecv=irecv+1

call MPI_IRECV(pw(recv_igrid)%w,1,type_block_io,send_ipe,itag, icomm,&
   recvrequest(irecv),ierrmpi)


end subroutine lb_recv
!=============================================================================
subroutine lb_send
!-----------------------------------------------------------------------------
itag=recv_igrid
isend=isend+1

call MPI_ISEND(pw(send_igrid)%w,1,type_block_io,recv_ipe,itag, icomm,&
   sendrequest(isend),ierrmpi)


end subroutine lb_send
!=============================================================================
! end of internal procedures
!=============================================================================
end subroutine load_balance
!=============================================================================
subroutine amr_Morton_order
! Construct Morton-order as a global recursive lexicographic ordering.

use mod_forest
use mod_global_parameters

integer :: ig1, Morton_no, isfc
!-----------------------------------------------------------------------------
Morton_no=0
nglev1=ng1(1)
do isfc=1,nglev1
   ig1=sfc_iglevel1(1,isfc) 
   call get_Morton_number(tree_root(ig1))
end do

if (Morton_no/=nleafs) then
   call mpistop("error in amr_Morton_order: Morton_no/=nleafs")
end if

contains
!=============================================================================
! internal procedures
!=============================================================================
recursive subroutine get_Morton_number(tree)

type(tree_node_ptr) :: tree

integer :: ic1
!-----------------------------------------------------------------------------
if (tree%node%leaf) then
   Morton_no=Morton_no+1
   sfc(1,Morton_no)=tree%node%igrid
   sfc(2,Morton_no)=tree%node%ipe
   if (tree%node%active) then 
      sfc(3,Morton_no)=1 
   else 
      sfc(3,Morton_no)=0 
   end if
   igrid_to_sfc(tree%node%igrid)=Morton_no
else
   do ic1=1,2
      call get_Morton_number(tree%node%child(ic1))
   end do
end if

end subroutine get_Morton_number
!=============================================================================
! end of internal procedures
!=============================================================================
end subroutine amr_Morton_order

!> Set the Morton range for each processor
subroutine get_Morton_range
use mod_forest
use mod_global_parameters

integer :: ipe, blocks_left, procs_left, num_blocks
!-----------------------------------------------------------------------------
if (allocated(sfc_to_igrid)) deallocate(sfc_to_igrid)


blocks_left = nleafs

do ipe = 0, npe-1
  if (ipe == 0) then
    Morton_start(ipe) = 1
  else
    Morton_start(ipe) = Morton_stop(ipe-1) + 1
  end if

  ! Compute how many blocks this cpu should take
  procs_left = npe - ipe
  num_blocks = ceiling(blocks_left / dble(procs_left))
  Morton_stop(ipe) = Morton_start(ipe) + num_blocks - 1
  blocks_left = blocks_left - num_blocks
end do

allocate(sfc_to_igrid(Morton_start(mype):Morton_stop(mype)))


end subroutine get_Morton_range
!=============================================================================
subroutine get_Morton_range_active
use mod_forest
use mod_global_parameters

! Cut the sfc based on weighted decision.  
! Oliver Porth, 02.02.2012

!!!Here you choose the weithts:!!
integer, parameter :: wa=3, wp=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! wp : Weight for passive block
! wa : Weight for active block
! wp=0 : balance load (active blocks) exactly and 
! don't care about memory imbalance
! wp=wa : balance memory exactly and don't care about load
! Maximum possible memory imbalance is X=wa/wp.

! If you run into memory issues, decrease this ratio.
! Scaling should be better if you allow for higher ratio.  
! Note also that passive cells still do regridding and boundary swap.  
! I have best results with a ratio 2:1, but it is problem dependent.
! It can't do magic though...
! Best to make sure that the sfc is properly aligned with your problem. 

integer :: ipe, Morton_no
integer :: Mtot, Mstop, Mcurr
! For debugging: 
integer :: nactive(0:npe-1),npassive(0:npe-1)
!double precision, save :: ptasum=0
!-----------------------------------------------------------------------------
if (allocated(sfc_to_igrid)) deallocate(sfc_to_igrid)


Mtot  = nleafs_active*wa+(nleafs-nleafs_active)*wp
ipe = 0 
Mcurr = 0

nactive=0
npassive=0

Morton_start(0) = 1
do Morton_no=1,nleafs
   ! This is where we ideally would like to make the cuts:
   Mstop  = (ipe+1)*int(Mtot/npe)+min(ipe+1,mod(Mtot,npe))
   ! Build up mass:
   Mcurr = Mcurr + (wa*sfc(3,Morton_no)+wp*(1-sfc(3,Morton_no)))

   if (sfc(3,Morton_no)==1) then 
      nactive(ipe) = nactive(ipe) +1
      else
         npassive(ipe) = npassive(ipe) +1
      end if

   if (Mcurr >= Mstop) then 
      Morton_stop(ipe) = Morton_no
      ipe = ipe +1
      if (ipe>=npe) exit
      Morton_start(ipe) = Morton_no + 1
   end if
end do

Xmemory=dble(maxval(npassive+nactive))/dble(minval(npassive+nactive))
Xload=dble(maxval(nactive))/dble(minval(nactive))

!ptasum = ptasum +dble(nleafs-nleafs_active)/dble(nleafs_active)

!if (mype == 0) print*, 'nleafs_passive:',nleafs-nleafs_active, 'nleafs_active:',nleafs_active,'ratio:',dble(nleafs-nleafs_active)/dble(nleafs_active),'mean ratio:',ptasum/it



if (Morton_stop(mype)>=Morton_start(mype)) then
   allocate(sfc_to_igrid(Morton_start(mype):Morton_stop(mype)))

end if

end subroutine get_Morton_range_active
!=============================================================================
