module network_indices

  ! this module is for use only within this network -- these
  ! quantities should not be accessed in general MAESTRO routines.
  ! Instead the species indices should be queried via
  ! network_species_index()

  implicit none

  integer, parameter :: ic12_ = 1
  integer, parameter :: io16_ = 2
  integer, parameter :: img24_ = 3
  !$acc declare copyin(ic12_, io16_, img24_)

end module network_indices

module rpar_indices

  implicit none

  integer, parameter :: n_rpar_comps = 3

  integer, parameter :: irp_dens = 1
  integer, parameter :: irp_temp = 2
  integer, parameter :: irp_o16  = 3
  !$acc declare copyin(n_rpar_comps, irp_dens, irp_temp, irp_o16)

end module rpar_indices
