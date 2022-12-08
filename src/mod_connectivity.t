!> This module contains variables that describe the connectivity of the mesh and
!> also data structures for connectivity-related communication.
module mod_connectivity
   implicit none
   save

   integer, parameter :: neighbor_boundary = 1
   integer, parameter :: neighbor_coarse = 2
   integer, parameter :: neighbor_sibling = 3
   integer, parameter :: neighbor_fine = 4

   integer, dimension(:,:^D&,:), allocatable :: neighbor
   integer, dimension(:,:^D&,:), allocatable :: neighbor_child
   integer, dimension(:^D&,:), allocatable :: neighbor_type
   logical, dimension(:^D&,:), allocatable :: neighbor_active
   integer, dimension(:^D&,:), allocatable :: neighbor_pole

   ! grid number array per processor
   integer, dimension(:), allocatable :: igrids
   integer, dimension(:), allocatable :: igrids_active
   integer, dimension(:), allocatable :: igrids_passive
   ! number of grids on current processor
   integer :: igridstail
   integer :: igridstail_active
   integer :: igridstail_passive

   integer, dimension(^ND) :: nrecv_fc, nsend_fc

   ! MPI receive-send pairs for boundary conditions
   ! srl: same refinement level; r: restrict; p: prolong
   integer :: nrecv_bc_srl, nsend_bc_srl, &
              nrecv_bc_r, nsend_bc_r, &
              nrecv_bc_p, nsend_bc_p
!$OMP THREADPRIVATE(leveljump)

end module mod_connectivity
