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

  character (len=*), parameter :: network_name = "urca_simple"

  ! nspec = number of species this network carries
  ! nspec_advance = the number of species that are explicitly integrated
  !                 in the ODE solve (the others are solved for 
  !                 algebraically).
  integer, parameter :: nspec = 2
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
    "neon-23         ",&
    "sodium-23       "]
    !!$acc update device(spec_names)

    allocate(short_spec_names(nspec))
    short_spec_names = [&
    "Ne23  ",&
    "Na23  "]
    !!$acc update device(short_spec_names)

    allocate(aion(nspec))
    aion = [&
    23.0_dp_t,&
    23.0_dp_t]
    !$acc update device(aion)
    !!$acc enter data copyin(aion)

    allocate(zion(nspec))
    zion = [&
    10.0_dp_t,& 
    11.0_dp_t]
    !$acc update device(zion)
    !!$acc enter data copyin(zion)

    allocate(ebin(nspec))
    ebin = [&
    -3.33911e17_dp_t,&     ! Ne23
    -3.40538e17_dp_t]      ! Na23
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
