!---------------------------------------------------------------------------------------------------
! SAS_GRID: A code to obtain the solvent accessible surface (SAS) around a given molecular structure                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2024 Felippe M. Colombari
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
!> @file   sas_grid.f90
!> @author Felippe M. Colombari
!> @brief  Main module of sas_grid 
!> @date - ~2018
!> - first version
!> @date - Nov, 2019                                                           
!> - modular version                                                
!---------------------------------------------------------------------------------------------------

program sas_grid
  use mod_error_handling    , only : error
  use mod_info              , only : Display_header, Display_date_time
  use mod_constants         , only : dashline
  use mod_cmd_line          , only : Parse_arguments, filename, factor, grid_type, box_min, box_max
  use mod_read_molecule     , only : mol 
  use mod_grids             , only : grid_sphere, grid_sas, grid_box 
  use mod_deallocate_all    , only : Deallocate_arrays

  implicit none

  type( error ) :: err

  call display_header()

  call Parse_arguments

  call mol % Read_molecule( filename )

  call mol % Read_vdw_radii

  select case( grid_type )
  
    case( 'sas' )

      call grid_sphere % Build_sphere_grid( factor )

      call grid_sas % Build_sas_grid

    case( 'box' )

      call grid_box % Build_box_grid( box_min, box_max )

  end select

  call Deallocate_arrays

  write(*,'(/, T3, A)') dashline

  call err % termination(0,'f')

end program sas_grid
