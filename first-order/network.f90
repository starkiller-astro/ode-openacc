! the network module provides the information about the species we are
! advecting: 
!
! nspec      -- the number of species
!
! aion       -- atomic number
! zion       -- proton number
! eion       -- nuclear binding energy (in erg/g)
!
! spec_names -- the name of the isotope
! short_spec_names -- the abbreviated name of the isotope
!
! This module contains two routines:
!
!  network_init()        -- initialize the isotope properties
!
!  network_species_index -- return the index of the species given its name
!

module network

  use bl_types
  use network_indices

  implicit none

  character (len=*), parameter :: network_name = "ignition_simple"

  ! nspec = number of species this network carries
  ! nspec_advance = the number of species that are explicitly integrated
  !                 in the ODE solve (the others are solved for 
  !                 algebraically).
  integer, parameter :: nspec = 3
  integer, parameter :: nspec_advance = 1
  integer, parameter :: naux  = 0

  character (len=16), parameter :: spec_names(nspec) = [& 
    "carbon-12       ",&
    "oxygen-16       ",&
    "magnesium-24    "]

  character (len= 5), parameter :: short_spec_names(nspec) = [&
    "C12  ",&
    "O16  ",&
    "Mg24 "]

  real(kind=dp_t), parameter :: aion(nspec) = [&
    12.0_dp_t,&
    16.0_dp_t,&
    24.0_dp_t]
  
  real(kind=dp_t), parameter :: zion(nspec) = [&
    6.0_dp_t,& 
    8.0_dp_t,& 
    12.0_dp_t]

  
  real(kind=dp_t), parameter :: ebin(nspec) = [&
    -7.4103097e18_dp_t,&     !  92.16294 MeV
    -7.6959672e18_dp_t,&     ! 127.62093 MeV
    -7.9704080e18_dp_t]     ! 198.2579  MeV

  !$acc declare copyin(short_spec_names, aion, zion, ebin) 
  !$acc declare copyin(network_name, nspec, nspec_advance, naux)
  !$acc declare copyin(spec_names) 
  
contains
  
  function network_species_index(name) result(r)

    character(len=*) :: name
    integer :: r, n

    r = -1

    do n = 1, nspec
       if (name == spec_names(n) .or. name == short_spec_names(n)) then
          r = n
          exit
       endif
    enddo
    return
  end function network_species_index


  subroutine network_finalize()

  end subroutine network_finalize

end module network
