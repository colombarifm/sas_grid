!---------------------------------------------------------------------------------------------------
! SAS_GRID: A code to obtain the solvent accessible surface (SAS) around a given molecular structure                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2020 Felippe M. Colombari
!
!---------------------------------------------------------------------------------------------------
!
!   This program is free software: you can redistribute it and/or modify it under the terms of the 
!   GNU General Public License as published by the Free Software Foundation, either version 3 of the 
!   License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!   without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
!   the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License along with this program. If 
!   not, see <https://www.gnu.org/licenses/>.
!
!---------------------------------------------------------------------------------------------------
!> @file   mod_read_molecule.f90
!> @author Felippe M. Colombari
!> @brief  This module reads coordinates for a given .xyz file and sets the vdW radii for each atom
!> @date - Nov, 2019                                                           
!> - independent module created                                                
!> @date - Nov 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!> @date - Sep, 2020
!> - vdW radii read from mod_constants.f90
!---------------------------------------------------------------------------------------------------

module mod_read_molecule
  use mod_constants

  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type atom 
    real( kind = dp )                       :: xyz(3), vdw
    character( len = 2 )                    :: label
  end type atom

  type molecule
    type( atom ), allocatable,dimension(:) :: atoms
    integer                                :: num_atoms
  contains
    procedure, pass                        :: Read_molecule
    procedure, pass                        :: Read_vdw_radii
  end type molecule

  type( molecule )                         :: mol

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type atom_entry
    character( len = 2 )                   :: symbol 
    real( kind = DP )                      :: radii 
  end type atom_entry

  logical, allocatable, dimension(:) :: radii_found

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

  !---------------------------------------------------------------------------
  !> @brief This routine reads coordinates for xyz_files.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Read_molecule ( this, molecule_filename )
    use mod_inquire, only: Inquire_file
    use mod_error_handling

    implicit none

    class( molecule ), intent(inout) :: this
    character( len = * ), intent(in) :: molecule_filename
    integer                          :: i
    integer                          :: ios         = 0
    integer                          :: file_unit   = 10        
    character( len = 15 )            :: file_format = "formatted"
    character( len = 15 )            :: file_access = "sequential"

    integer                          :: ierr
    type(error)                      :: err

    call Inquire_file( file_unit, molecule_filename, file_format, file_access )

    read(file_unit,*,iostat=ios) this % num_atoms
    read(file_unit,*)

    if ( allocated ( this % atoms ) ) deallocate ( this % atoms )
    allocate( this % atoms( this % num_atoms ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    do i = 1, this % num_atoms

      read(file_unit,*,iostat=ios) this % atoms(i) % label, this % atoms(i) % xyz(:)

    enddo

    close(file_unit)

    return
  end subroutine Read_molecule

  subroutine Read_vdw_radii( this )
    use mod_cmd_line , only : radius

    implicit none

    class( molecule ), intent(inout)     :: this
    integer                              :: i

    allocate( radii_found( mol % num_atoms ) )

    radii_found = .false. 

    do i = 1, mol % num_atoms

      if ( mol % atoms(i) % label == "H" ) then
        
        this % atoms(i) % vdw = radii_H + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "He" ) then
        
        this % atoms(i) % vdw = radii_He + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Li" ) then
        
        this % atoms(i) % vdw = radii_Li + radius

        radii_found(i) = .true.
      
      else if ( mol % atoms(i) % label == "Be" ) then
        
        this % atoms(i) % vdw = radii_Li + radius

        radii_found(i) = .true.
      
      else if ( mol % atoms(i) % label == "B" ) then
        
        this % atoms(i) % vdw = radii_Li + radius

        radii_found(i) = .true.
      
      else if ( mol % atoms(i) % label == "C" ) then
        
        this % atoms(i) % vdw = radii_C + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "N" ) then
        
        this % atoms(i) % vdw = radii_N + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "O" ) then

        this % atoms(i) % vdw = radii_O + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "F" ) then

        this % atoms(i) % vdw = radii_F + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Ne" ) then

        this % atoms(i) % vdw = radii_Ne + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Na" ) then

        this % atoms(i) % vdw = radii_Na + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Mg" ) then

        this % atoms(i) % vdw = radii_Mg + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Al" ) then

        this % atoms(i) % vdw = radii_Al + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Si" ) then

        this % atoms(i) % vdw = radii_Si + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "P" ) then

        this % atoms(i) % vdw = radii_P + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "S" ) then

        this % atoms(i) % vdw = radii_S + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Cl" ) then

        this % atoms(i) % vdw = radii_Cl + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Ar" ) then

        this % atoms(i) % vdw = radii_Ar + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "K" ) then

        this % atoms(i) % vdw = radii_K + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Ca" ) then

        this % atoms(i) % vdw = radii_Ca + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Ga" ) then

        this % atoms(i) % vdw = radii_Ga + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Ge" ) then

        this % atoms(i) % vdw = radii_Ge + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "As" ) then

        this % atoms(i) % vdw = radii_As + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Se" ) then

        this % atoms(i) % vdw = radii_Se + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Br" ) then

        this % atoms(i) % vdw = radii_Br + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Kr" ) then

        this % atoms(i) % vdw = radii_Kr + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Rb" ) then

        this % atoms(i) % vdw = radii_Rb + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Sr" ) then

        this % atoms(i) % vdw = radii_Sr + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "In" ) then

        this % atoms(i) % vdw = radii_In + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Sn" ) then

        this % atoms(i) % vdw = radii_Sn + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Sb" ) then

        this % atoms(i) % vdw = radii_Sb + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Te" ) then

        this % atoms(i) % vdw = radii_Te + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "I" ) then

        this % atoms(i) % vdw = radii_I + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Xe" ) then

        this % atoms(i) % vdw = radii_Xe + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Cs" ) then

        this % atoms(i) % vdw = radii_Cs + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Ba" ) then

        this % atoms(i) % vdw = radii_Ba + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Tl" ) then

        this % atoms(i) % vdw = radii_Tl + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Pb" ) then

        this % atoms(i) % vdw = radii_Pb + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Bi" ) then

        this % atoms(i) % vdw = radii_Bi + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Po" ) then

        this % atoms(i) % vdw = radii_Po + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "At" ) then

        this % atoms(i) % vdw = radii_At + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Rn" ) then

        this % atoms(i) % vdw = radii_Rn + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Fr" ) then

        this % atoms(i) % vdw = radii_Fr + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Ra" ) then

        this % atoms(i) % vdw = radii_Ra + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Pt" ) then

        this % atoms(i) % vdw = radii_Pt + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Pd" ) then

        this % atoms(i) % vdw = radii_Pd + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Ag" ) then

        this % atoms(i) % vdw = radii_Ag + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Au" ) then

        this % atoms(i) % vdw = radii_Au + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Cu" ) then

        this % atoms(i) % vdw = radii_Cu + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Zn" ) then

        this % atoms(i) % vdw = radii_Zn + radius

        radii_found(i) = .true.
      
      else if ( mol % atoms(i) % label == "W" ) then

        this % atoms(i) % vdw = radii_W + radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "X" ) then

        this % atoms(i) % vdw = radii_X 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "XX" ) then

        this % atoms(i) % vdw = radii_XX 

        radii_found(i) = .true.

      endif

      if ( radii_found(i) .eqv. .false. ) then

        write(*,'(T5, "Atom ", i4.4,": ", a2, " not found!")') i, this % atoms(i) % label
        write(*,*)
        write(*,'(T5, "Check the mod_constants.f90 file.")')

        stop

      endif

    enddo
      
  end subroutine Read_vdw_radii

end module mod_read_molecule
