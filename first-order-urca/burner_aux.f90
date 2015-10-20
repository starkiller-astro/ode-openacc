module network_indices

  ! this module is for use only within this network -- these
  ! quantities should not be accessed in general MAESTRO routines.
  ! Instead the species indices should be queried via
  ! network_species_index()

  implicit none

  integer, parameter :: ine23_ = 1
  integer, parameter :: ina23_ = 2
  !$acc declare copyin(ine23_, ina23_)

end module network_indices

module rpar_indices

  implicit none

  integer, parameter :: n_rpar_comps = 2

  integer, parameter :: irp_dens = 1
  integer, parameter :: irp_temp = 2
  !$acc declare copyin(n_rpar_comps, irp_dens, irp_temp)

end module rpar_indices
