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
  !$acc declare copyin(nspec, nspec_advance, naux)

  character (len=16), allocatable :: spec_names(:)
  !TODO: Commented out because compilers 
  !      don't like Fortran character arrays on GPUs
  !!$acc declare create(spec_names)

  character (len= 5), allocatable :: short_spec_names(:)
  !!$acc declare create(short_spec_names)

  real(kind=dp_t), allocatable :: aion(:)
  !$acc declare create(aion)
  
  real(kind=dp_t), allocatable :: zion(:)
  !$acc declare create(zion)
  
  real(kind=dp_t), allocatable :: ebin(:)
  !$acc declare create(ebin)

contains
  
  subroutine network_init()
    allocate(spec_names(nspec))
    spec_names = [& 
    "carbon-12       ",&
    "oxygen-16       ",&
    "magnesium-24    "]
    !!$acc update device(spec_names)

    allocate(short_spec_names(nspec))
    short_spec_names = [&
    "C12  ",&
    "O16  ",&
    "Mg24 "]
    !!$acc update device(short_spec_names)

    allocate(aion(nspec))
    aion = [&
    12.0_dp_t,&
    16.0_dp_t,&
    24.0_dp_t]
    !$acc update device(aion)
    !!$acc enter data copyin(aion)

    allocate(zion(nspec))
    zion = [&
    6.0_dp_t,& 
    8.0_dp_t,& 
    12.0_dp_t]
    !$acc update device(zion)
    !!$acc enter data copyin(zion)

    allocate(ebin(nspec))
    ebin = [&
    -7.4103097e18_dp_t,&     !  92.16294 MeV
    -7.6959672e18_dp_t,&     ! 127.62093 MeV
    -7.9704080e18_dp_t]      ! 198.2579  MeV
    !$acc update device(ebin)
    !!$acc enter data copyin(ebin)
  end subroutine network_init

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
