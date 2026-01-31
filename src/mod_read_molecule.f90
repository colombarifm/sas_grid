!---------------------------------------------------------------------------------------------------
! SAS_GRID: A code to obtain the solvent accessible surface (SAS) around a given molecular structure                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2026 Felippe M. Colombari
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
!> @brief  This module reads coordinates for a given .xyz file and sets the rvdw radii for each atom
!> @date - Nov, 2019                                                           
!> - independent module created                                                
!> @date - Nov 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!> @date - Sep, 2020
!> - rvdw radii read from mod_constants.f90
!---------------------------------------------------------------------------------------------------

module mod_read_molecule
  use mod_constants

  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type atom 
    real( kind = dp )                       :: xyz(3), rvdw, rvdw_solv
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

  logical, allocatable, dimension(:)       :: radii_found

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

  !---------------------------------------------------------------------------
  !> @brief This routine reads coordinates for xyz_files.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Read_molecule ( this, molecule_filename )
    use mod_inquire        , only : Inquire_file
    use mod_error_handling , only : error

    implicit none

    class( molecule ), intent(inout) :: this
    character( len = * ), intent(in) :: molecule_filename
    integer                          :: i
    integer                          :: ios         = 0
    integer                          :: file_unit   = 10        
    character( len = 15 )            :: file_format = "formatted"
    character( len = 15 )            :: file_access = "sequential"
    integer                          :: ierr
    type( error )                    :: err

    call Inquire_file( file_unit, molecule_filename, file_format, file_access )

    read(file_unit,*,iostat=ios) this % num_atoms
    read(file_unit,*)

    if ( allocated ( this % atoms ) ) deallocate ( this % atoms )
    allocate( this % atoms( this % num_atoms ), stat=ierr )
    if(ierr/=0) call err % error('e',message="abnormal memory allocation")

    do i = 1, this % num_atoms

      read(file_unit,*,iostat=ios) this % atoms(i) % label, this % atoms(i) % xyz(:)

    enddo

    close(file_unit)

    return
  end subroutine Read_molecule

  subroutine Read_vdw_radii( this )
    use mod_cmd_line , only : solv_radius

    implicit none

    class( molecule ), intent(inout)     :: this
    integer                              :: i

    allocate( radii_found( mol % num_atoms ) )

    radii_found = .false. 

    do i = 1, mol % num_atoms

      if ( mol % atoms(i) % label == "H" ) then
        
        this % atoms(i) % rvdw      = radii_H 
        this % atoms(i) % rvdw_solv = radii_H + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "He" ) then
        
        this % atoms(i) % rvdw      = radii_He 
        this % atoms(i) % rvdw_solv = radii_He + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Li" ) then
        
        this % atoms(i) % rvdw      = radii_Li 
        this % atoms(i) % rvdw_solv = radii_Li + solv_radius

        radii_found(i) = .true.
      
      else if ( mol % atoms(i) % label == "Be" ) then
        
        this % atoms(i) % rvdw      = radii_Be 
        this % atoms(i) % rvdw_solv = radii_Be + solv_radius 

        radii_found(i) = .true.
      
      else if ( mol % atoms(i) % label == "B" ) then
        
        this % atoms(i) % rvdw      = radii_B 
        this % atoms(i) % rvdw_solv = radii_B + solv_radius 

        radii_found(i) = .true.
      
      else if ( mol % atoms(i) % label == "C" ) then
        
        this % atoms(i) % rvdw      = radii_C 
        this % atoms(i) % rvdw_solv = radii_C + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "N" ) then
        
        this % atoms(i) % rvdw      = radii_N 
        this % atoms(i) % rvdw_solv = radii_N + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "O" ) then

        this % atoms(i) % rvdw      = radii_O 
        this % atoms(i) % rvdw_solv = radii_O + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "F" ) then

        this % atoms(i) % rvdw      = radii_F 
        this % atoms(i) % rvdw_solv = radii_F + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Ne" ) then

        this % atoms(i) % rvdw      = radii_Ne 
        this % atoms(i) % rvdw_solv = radii_Ne + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Na" ) then

        this % atoms(i) % rvdw      = radii_Na 
        this % atoms(i) % rvdw_solv = radii_Na + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Mg" ) then

        this % atoms(i) % rvdw      = radii_Mg 
        this % atoms(i) % rvdw_solv = radii_Mg + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Al" ) then

        this % atoms(i) % rvdw      = radii_Al 
        this % atoms(i) % rvdw_solv = radii_Al + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Si" ) then

        this % atoms(i) % rvdw      = radii_Si 
        this % atoms(i) % rvdw_solv = radii_Si + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "P" ) then

        this % atoms(i) % rvdw      = radii_P 
        this % atoms(i) % rvdw_solv = radii_P + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "S" ) then

        this % atoms(i) % rvdw      = radii_S 
        this % atoms(i) % rvdw_solv = radii_S + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Cl" ) then

        this % atoms(i) % rvdw      = radii_Cl 
        this % atoms(i) % rvdw_solv = radii_Cl + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Ar" ) then

        this % atoms(i) % rvdw      = radii_Ar 
        this % atoms(i) % rvdw_solv = radii_Ar + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "K" ) then

        this % atoms(i) % rvdw      = radii_K 
        this % atoms(i) % rvdw_solv = radii_K + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Ca" ) then

        this % atoms(i) % rvdw      = radii_Ca 
        this % atoms(i) % rvdw_solv = radii_Ca + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Ga" ) then

        this % atoms(i) % rvdw      = radii_Ga 
        this % atoms(i) % rvdw_solv = radii_Ga + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Ge" ) then

        this % atoms(i) % rvdw      = radii_Ge 
        this % atoms(i) % rvdw_solv = radii_Ge + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "As" ) then

        this % atoms(i) % rvdw      = radii_As 
        this % atoms(i) % rvdw_solv = radii_As + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Se" ) then

        this % atoms(i) % rvdw      = radii_Se 
        this % atoms(i) % rvdw_solv = radii_Se + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Br" ) then

        this % atoms(i) % rvdw      = radii_Br 
        this % atoms(i) % rvdw_solv = radii_Br + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Kr" ) then

        this % atoms(i) % rvdw      = radii_Kr 
        this % atoms(i) % rvdw_solv = radii_Kr + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Rb" ) then

        this % atoms(i) % rvdw      = radii_Rb 
        this % atoms(i) % rvdw_solv = radii_Rb + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Sr" ) then

        this % atoms(i) % rvdw      = radii_Sr 
        this % atoms(i) % rvdw_solv = radii_Sr + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "In" ) then

        this % atoms(i) % rvdw      = radii_In 
        this % atoms(i) % rvdw_solv = radii_In + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Sn" ) then

        this % atoms(i) % rvdw      = radii_Sn 
        this % atoms(i) % rvdw_solv = radii_Sn + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Sb" ) then

        this % atoms(i) % rvdw      = radii_Sb 
        this % atoms(i) % rvdw_solv = radii_Sb + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Te" ) then

        this % atoms(i) % rvdw      = radii_Te 
        this % atoms(i) % rvdw_solv = radii_Te + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "I" ) then

        this % atoms(i) % rvdw      = radii_I 
        this % atoms(i) % rvdw_solv = radii_I + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Xe" ) then

        this % atoms(i) % rvdw      = radii_Xe 
        this % atoms(i) % rvdw_solv = radii_Xe + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Cs" ) then

        this % atoms(i) % rvdw      = radii_Cs 
        this % atoms(i) % rvdw_solv = radii_Cs + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Ba" ) then

        this % atoms(i) % rvdw      = radii_Ba 
        this % atoms(i) % rvdw_solv = radii_Ba + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Tl" ) then

        this % atoms(i) % rvdw      = radii_Tl 
        this % atoms(i) % rvdw_solv = radii_Tl + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Pb" ) then

        this % atoms(i) % rvdw      = radii_Pb 
        this % atoms(i) % rvdw_solv = radii_Pb + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Bi" ) then

        this % atoms(i) % rvdw      = radii_Bi 
        this % atoms(i) % rvdw_solv = radii_Bi + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Po" ) then

        this % atoms(i) % rvdw      = radii_Po 
        this % atoms(i) % rvdw_solv = radii_Po + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "At" ) then

        this % atoms(i) % rvdw      = radii_At 
        this % atoms(i) % rvdw_solv = radii_At + solv_radius 

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Rn" ) then

        this % atoms(i) % rvdw      = radii_Rn 
        this % atoms(i) % rvdw_solv = radii_Rn + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Fr" ) then

        this % atoms(i) % rvdw      = radii_Fr 
        this % atoms(i) % rvdw_solv = radii_Fr + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Ra" ) then

        this % atoms(i) % rvdw      = radii_Ra 
        this % atoms(i) % rvdw_solv = radii_Ra + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Pt" ) then

        this % atoms(i) % rvdw      = radii_Pt 
        this % atoms(i) % rvdw_solv = radii_Pt + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Pd" ) then

        this % atoms(i) % rvdw      = radii_Pd 
        this % atoms(i) % rvdw_solv = radii_Pd + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Ag" ) then

        this % atoms(i) % rvdw      = radii_Ag 
        this % atoms(i) % rvdw_solv = radii_Ag + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Au" ) then

        this % atoms(i) % rvdw      = radii_Au 
        this % atoms(i) % rvdw_solv = radii_Au + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Cu" ) then

        this % atoms(i) % rvdw      = radii_Cu 
        this % atoms(i) % rvdw_solv = radii_Cu + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Zn" ) then

        this % atoms(i) % rvdw      = radii_Zn 
        this % atoms(i) % rvdw_solv = radii_Zn + solv_radius

        radii_found(i) = .true.
      
      else if ( mol % atoms(i) % label == "W" ) then

        this % atoms(i) % rvdw      = radii_W 
        this % atoms(i) % rvdw_solv = radii_W + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Mo" ) then

        this % atoms(i) % rvdw      = radii_Mo 
        this % atoms(i) % rvdw_solv = radii_Mo + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Cd" ) then

        this % atoms(i) % rvdw      = radii_Cd 
        this % atoms(i) % rvdw_solv = radii_Cd + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "Fe" ) then

        this % atoms(i) % rvdw      = radii_Fe 
        this % atoms(i) % rvdw_solv = radii_Fe + solv_radius

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "X" ) then

        this % atoms(i) % rvdw      = radii_X 
        this % atoms(i) % rvdw_solv = radii_X

        radii_found(i) = .true.

      else if ( mol % atoms(i) % label == "XX" ) then

        this % atoms(i) % rvdw      = radii_XX 
        this % atoms(i) % rvdw_solv = radii_XX 

        radii_found(i) = .true.

      endif

      if ( radii_found(i) .eqv. .false. ) then

        write(*,'(T5, "WARNING:")') 
        write(*,'(T5, "Radius value for atom ", i4.4," (", a, ") not found!")') i, trim(this % atoms(i) % label)
        write(*,*)
        write(*,'(T5, "Check the mod_constants.f90 file.")')
        write(*,*)
        write(*,'(T5, "Default value of 2.0 was used!")')

        this % atoms(i) % rvdw      = 2.0_DP 
        this % atoms(i) % rvdw_solv = 2.0_DP + solv_radius

        !stop

      endif

    enddo
      
  end subroutine Read_vdw_radii

end module mod_read_molecule
