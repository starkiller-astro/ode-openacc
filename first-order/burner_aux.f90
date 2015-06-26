module network_indices

  ! this module is for use only within this network -- these
  ! quantities should not be accessed in general MAESTRO routines.
  ! Instead the species indices should be queried via
  ! network_species_index()

  implicit none

  integer, parameter :: ic12_ = 1
  integer, parameter :: io16_ = 2
  integer, parameter :: img24_ = 3

end module network_indices

module rpar_indices

  implicit none

  integer, save :: n_rpar_comps = 0

  integer, save :: irp_dens, irp_temp, irp_o16

contains

  function get_next_rpar_index(num) result (next)

    ! return the next starting index for a plotfile quantity,
    ! and increment the counter of plotfile quantities by num
    integer :: num, next

    next = n_rpar_comps + 1
    n_rpar_comps = n_rpar_comps + num

    return
  end function get_next_rpar_index


  subroutine init_rpar_indices(nspec)

    integer, intent(in) :: nspec

    irp_dens  = get_next_rpar_index(1)
    irp_temp  = get_next_rpar_index(1)
    irp_o16   = get_next_rpar_index(1)

  end subroutine init_rpar_indices

end module rpar_indices
