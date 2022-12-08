module mod_obj_star
  use mod_constants
  use mod_global_parameters
  use mod_obj_global_parameters
  use mod_obj_mat
  use mod_physics
  use mod_hd_grackle, only :hd_dust
  use mod_srmhd_parameters!, only: mag,lfac_,psi_,xi_
implicit none
  ! Star features
  type  star
  character(len=20)  ::unit  !> physical unit at parameter file
  real(dp)  :: radius     !> star radius (cm)
  real(dp)  :: mass       !> star mass  (g)
  real(dp)  :: luminosity !> star luminosity in erg/s
  real(dp)  :: Eddington  !> Eddington limit erg/s
  real(dp)  :: temperature!> star temperature
  real(dp)  :: vrotation  !> rotation parameter cm/s
  real(dp)  :: magnetic   !> Magnetic field strength at star surface (gauss)
  real(dp)  :: eta
  real(dp)  :: frac_critical_rotation
  end type
contains

end module mod_obj_star
