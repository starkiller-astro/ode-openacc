module network_indices

  ! this module is for use only within this network -- these
  ! quantities should not be accessed in general MAESTRO routines.
  ! Instead the species indices should be queried via
  ! network_species_index()

  implicit none

  integer, parameter :: ine23_ = 1
  integer, parameter :: ina23_ = 2

end module network_indices

module rpar_indices

  implicit none

  integer, parameter :: n_rpar_comps = 3

  integer, parameter :: irp_dens = 1
  integer, parameter :: irp_temp = 2
  integer, parameter :: irp_mu_elec = 3

end module rpar_indices

module phase_par_indices

  implicit none

  integer, parameter :: n_phase_par_inds = 3

  integer, parameter  ::  iqn   = 1 
  integer, parameter  ::  itemp = 2
  integer, parameter  ::  ieta  = 3

end module phase_par_indices
