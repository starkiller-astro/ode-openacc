module rpar_indices

  use network, only : nrat

  implicit none

  integer, save :: n_rpar_comps = 11 + nrat

  integer, parameter :: irp_dens = 1
  integer, parameter :: irp_temp = 2
  integer, parameter :: irp_dlambCNOdh1 = 3
  integer, parameter :: irp_drs1dhe4 = 4
  integer, parameter :: irp_drr1dh1 = 5
  integer, parameter :: irp_dlambda1dhe4 = 6
  integer, parameter :: irp_dlambda2dhe4 = 7
  integer, parameter :: irp_delta1 = 8
  integer, parameter :: irp_delta2 = 9
  integer, parameter :: irp_r56eff = 10
  integer, parameter :: irp_dr56effdt = 11
  integer, parameter :: irp_rates = 12 ! nrat components
  !$acc declare copyin(n_rpar_comps, irp_dens, irp_temp)
  !$acc declare copyin(irp_dlambCNOdh1, irp_drs1dhe4)
  !$acc declare copyin(irp_drr1dh1, irp_dlambda1dhe4, irp_dlambda2dhe4)
  !$acc declare copyin(irp_delta1, irp_delta2, irp_r56eff)
  !$acc declare copyin(irp_dr56effdt, irp_rates)
end module rpar_indices
