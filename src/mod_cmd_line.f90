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
!> @file   mod_cmd_line.f90
!> @author Felippe M. Colombari
!> @brief  Get command line arguments         
!> @date - Oct, 2019                                                           
!> - independent module created                                                
!> @date - Nov, 2019
!> - update error condition by error_handling module 
!> @date - Out, 2025
!> - added "--type" option to select "box" or "sas" grid types
!> - added "--min" and "--max" options to set coordinate ranges for the box
!---------------------------------------------------------------------------------------------------

module mod_cmd_line
  use iso_fortran_env    , only : stdout => output_unit
  use mod_constants      , only : DP, int_alphabet, float_alphabet, char_alphabet, dashline
  use mod_error_handling , only : error

  implicit none

  private
  public :: Parse_Arguments, solv_radius, factor, grid_type, box_min, box_max, filename

  integer                                          :: factor = 0
  real( kind = DP )                                :: solv_radius = 0.0_DP
  real( kind = DP ), dimension(3)                  :: box_min = 0.0_DP
  real( kind = DP ), dimension(3)                  :: box_max = 0.0_DP
  character( len = 3 )                             :: grid_type = char(0)
  character( len = 64 )                            :: filename = char(0)
  character( len = 20 ), allocatable, dimension(:) :: arg         
  integer                                          :: ierr
  type(error)                                      :: err

contains
  
  !---------------------------------------------------------------------------
  !> @brief Parses the command line arguments
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Parse_arguments

    implicit none

    integer                :: i, j
    integer                :: ios      = 0
    integer                :: narg     = 0
    integer                :: nochar   = 0
    logical                :: is_box   = .false.
    logical                :: is_sas   = .false.
    character( len = 256 ) :: cmd_line = char(0)
      
    narg = command_argument_count()

    call get_command(cmd_line)

    write(stdout,'(/,T5, A, A)') "COMMAND LINE READ: ", trim(cmd_line)
    write(stdout,'(/,T3, A)') dashline

    if ( narg > 0 ) then
    
      ! to avoid allocation errors if one forget the argument "rad"

      allocate( arg(narg+1), stat = ierr )
      if( ierr /= 0 ) call err % error('e',message="abnormal memory allocation")

      arg = char(0)

      i = 1
      do while ( i <= narg )

        call get_command_argument( i, arg(i) )

        if ( arg(i)(1:2) == '--' ) then

          select case( arg(i) )

            case( '--help' )

              call Display_help

            case( '--license' )

              call Display_license

            case( '--version')  

              call display_version

            case( '--radius' )

              if ( i <= narg ) then

                call Get_command_argument( i+1, arg(i+1) )

                nochar = verify( trim( arg(i+1) ), float_alphabet )

                if ( nochar > 0 ) then

                  call err % error('e', message = "while reading command line.")
                  call err % error('e', check = "solvent probe radius around atoms.") 
                  call err % error('e', tip = "Its value (in Angstrom) should be > 0.1.")

                  stop

                else

                  read( arg(i+1), *, iostat = ios ) solv_radius
                    
                  if ( ios > 0 ) then

                    call err % error('e', message = "while reading command line.")
                    call err % error('e', check = "solvent probe radius around atoms.") 
                    call err % error('e', tip = "Its value (in Angstrom) should be > 0.1.")

                    stop

                  endif

                  if ( solv_radius <= 0.0_DP ) then

                    call err % error('e', message = "while reading command line.")
                    call err % error('e', check = "solvent probe radius around atoms.") 
                    call err % error('e', tip = "Its value (in Angstrom) should be > 0.1.")
                    stop

                  endif

                endif

              else

                write(*,*) "some error"

              endif

              i = i + 2

            CASE( '--input' )

              if ( i <= narg ) then

                call Get_command_argument( i+1, arg(i+1) )

                nochar = verify( trim( arg(i+1) ), char_alphabet )

                if ( nochar > 0 ) then

                  call err % error('e', message = "while reading command line.")
                  call err % error('e', check = "molecule coordinate file.") 
                  call err % error('e', tip = "Should be a valid .xyz file.")

                  stop

                else

                  read(arg(i+1), *, iostat = ios) filename

                  if ( ios > 0 ) then

                    call err % error('e', message = "while reading command line.")
                    call err % error('e', check = "molecule coordinate file.") 
                    call err % error('e', tip = "Should be a valid .xyz file.")

                    stop

                  endif

                endif

              else

                write(*,*) "some error 2"

              endif
              
              i = i + 2

            case( '--type' )

              if ( i <= narg ) then

                call Get_command_argument( i+1, arg(i+1) )

                nochar = verify( trim( arg(i+1) ), char_alphabet )

                if ( nochar > 0 ) then

                  call err % error('e', message = "while reading command line.")
                  call err % error('e', check = "resullting grid type.") 
                  call err % error('e', tip = "Its value (a string) should be either 'box' or 'sas'.")

                  stop

                else

                  read(arg(i+1),*,iostat=ios) grid_type

                  if ( ios > 0 ) then                

                    call err % error('e', message = "while reading command line.")
                    call err % error('e', check = "resulting grid type.") 
                    call err % error('e', tip = "Its value (a string) should be either 'box' or 'sas'.")

                    stop

                  endif

                  if ( grid_type == 'box' ) then

                    is_box = .true.

                  else if ( grid_type == 'sas' ) then

                    is_sas = .true.

                  else

                    call err % error('e', message = "while reading command line.")
                    call err % error('e', check = "resulting grid type.") 
                    call err % error('e', tip = "Its value (a string) should be either 'box' or 'sas'.")

                    stop

                  endif

                endif

              endif
                    
              i = i + 2

            case( '--factor' )

              if ( i <= narg ) then

                call Get_command_argument( i+1, arg(i+1) )

                nochar = verify( trim( arg(i+1) ), int_alphabet )

                if ( nochar > 0 ) then

                  call err % error('e', message = "while reading command line.")
                  call err % error('e', check = "sphere tessellation factor.") 
                  call err % error('e', tip = "Its value (an integer) should be > 1.")

                  stop       

                else

                  read(arg(i+1),*,iostat=ios) factor

                  if ( ios > 0 ) then

                    call err % error('e', message = "while reading command line.")
                    call err % error('e', check = "sphere tessellation factor.") 
                    call err % error('e', tip = "Its value (an integer) should be > 1.")

                    stop

                  endif

                  if ( factor <= 1.0_DP ) then

                    call err % error('e', message = "while reading command line.")
                    call err % error('e', check = "sphere tessellation factor.") 
                    call err % error('e', tip = "Its value (an integer) should be > 1.")

                    stop

                  endif

                endif

              else

                write(*,*) "some error 4"

              endif

            i = i + 2

            case( '--min' )

              if ( i + 3 <= narg ) then

                do j = 1, 3

                  call Get_command_argument( i+j, arg(i+j) )

                  nochar = verify( trim( arg(i+j) ), float_alphabet )

                  if ( nochar > 0 ) then

                    call err % error('e', message = "while reading command line.")
                    call err % error('e', check = "min coordinates of box-type grid.") 
                    call err % error('e', tip = "XYZ values (in Angstrom) are expected.")

                    stop               

                  else

                    read(arg(i+j), *, iostat = ios) box_min(j)

                    if ( ios > 0 ) then

                      call err % error('e', message = "while reading command line.")
                      call err % error('e', check = "min coordinates of box-type grid.") 
                      call err % error('e', tip = "XYZ values (in Angstrom) are expected.")

                      stop

                    endif

                  endif

                enddo

              else

                write(*,*) "some error min"

              endif
            
              i = i + 4
 
            case( '--max' )

              if ( i + 3 <= narg ) then

                do j = 1, 3

                  call Get_command_argument( i+j, arg(i+j) )

                  nochar = verify( trim( arg(i+j) ), float_alphabet )

                  if ( nochar > 0 ) then

                    call err % error('e', message = "while reading command line.")
                    call err % error('e', check = "center of box-type grid.") 
                    call err % error('e', tip = "XYZ values (in Angstrom) are expected.")

                    stop     

                  else

                    read(arg(i+j), *, iostat = ios) box_max(j)

                    if ( ios > 0 ) then

                      call err % error('e', message = "while reading command line.")
                      call err % error('e', check = "min coordinates of box-type grid.") 
                      call err % error('e', tip = "XYZ values (in Angstrom) are expected.")

                      stop

                    endif

                  endif

                enddo

              else

                write(*,*) "some error max"

              endif

              i = i + 4

            case default

              call err % error('e', message = "while reading command line.")
              call err % error('e', check = "invalid command line flag '"//trim(adjustl(arg(i)))//"'.")
               
              stop

          end select

        else 

            call err % error('e', message = "while reading command line.")
            call err % error('e', check = "'"//trim(adjustl(arg(i)))//"' flag.")

            stop

          !if ( arg(1)(1:2) /= '--' ) then

          !  call err % error('e',message="while reading command line.")
          !  call err % error('e',check="'"//trim(adjustl(arg(i+1)))//"' argument of '"//trim(adjustl(arg(i)))//"' flag.")

          !  stop

          !endif

          !if ( ( i > 1 ) .and. ( arg(i-1)(1:2) ) /= '--' ) then

          !  call err % error('e',message="while reading command line.")
          !  call err % error('e',check="'"//trim(adjustl(arg(i+1)))//"' argument of '"//trim(adjustl(arg(i)))//"' flag.")
          !  stop

          !endif

        endif

      enddo

      if ( allocated(arg) ) deallocate(arg)

    else if ( narg == 0 ) then 

      call err % error('e', message = "while reading command line.")
      call err % error('e', tip = "Command line arguments are missing.")

      stop

    endif

    if ( ( is_box .eqv. .false. ) .and. ( is_sas .eqv. .false. ) ) then

      call err % error('e', message = "while reading command line.")
      call err % error('e', check = "resulting grid type.") 
      call err % error('e', tip = "either '--type sas' or '--type box' must be set. Aborting....")

      stop
      
    endif

    if ( ( is_box .eqv. .true. ) .and. ( is_sas .eqv. .true. ) ) then

      call err % error('e', message = "while reading command line.")
      call err % error('e', check = "resulting grid type.") 
      call err % error('e', tip = "either '--type sas' or '--type box' must be set. Aborting....")

      stop
      
    endif

    if ( ( is_box .eqv. .true. ) .and. ( factor > 0 ) ) then

      call err % error('e', message = "while reading command line.")
      call err % error('e', check = "resulting grid type.") 
      call err % error('e', tip = "--factor sould be used with '--type sas' only. Aborting....")

      stop
  
    endif

    if ( ( is_sas .eqv. .true. ) .and. ( any(box_min /= 0.0_DP) ) ) then

      call err % error('e', message = "while reading command line.")
      call err % error('e', check = "resulting grid type.") 
      call err % error('e', tip = "--min sould be used with '--type box' only. Aborting....")

      stop

    endif
    
    if ( ( is_sas .eqv. .true. ) .and. ( any(box_max /= 0.0_DP) ) ) then

      call err % error('e', message = "while reading command line.")
      call err % error('e', check = "resulting grid type.") 
      call err % error('e', tip = "--max sould be used with '--type box' only. Aborting....")

      stop

    endif
    
    if ( ( is_sas .eqv. .true. ) .and. ( solv_radius == 0.0_DP ) ) then

      call err % error('e', message = "while reading command line.")
      call err % error('e', check = "radius option.")
      call err % error('e', tip = "values should be > 0.1.")

      stop

    else if ( ( is_sas .eqv. .true. ) .and. ( factor == 0 ) ) then

      call err % error('e', message = "while reading command line.")
      call err % error('e', check = "factor option.")
      call err % error('e', tip = "values should be > 1.")

      stop

    else if ( filename == char(0) ) then

      call err % error('e', message = "while reading command line.")
      call err % error('e', check = "molecule coordinate file.") 
      call err % error('e', tip = "Should be a valid .xyz file.")

      stop

    endif

  end subroutine Parse_arguments

  !---------------------------------------------------------------------------
  !> @brief Displays command line options
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Display_help
            
    implicit none

    write(stdout,'(/,T20, A)')'Usage:  sas_grid --input [FILE] --radius [RADIUS] --factor [FACTOR]     '
    write(stdout,'(/,T3, A)') dashline
    write(stdout,'(/,T25, A)')'[FILE]   is a .xyz coordinate file.'
    write(stdout,'(/,T25, A)')'[RADIUS] is the solvent probe radius, in Angstrom.'
    write(stdout,'(/,T25, A)')'[FACTOR] is an integer factor for the tessellation sphere.' 
    write(stdout,'(T33, A)')' N_points = 2 + factor^2 * 10' 
    write(stdout,'(/,T3, A)') dashline
    
    if ( allocated(arg) ) deallocate(arg)

    call err%termination(0,'f')

  end subroutine Display_help

  !---------------------------------------------------------------------------
  !> @brief Displays the license
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Display_license

    implicit none

    write(stdout,'(T36, A)')'Copyright 2020 Felippe M. Colombari'
    write(stdout,'(/,T33, A)')'License GPLv3+: GNU GPL version 3 or later' 
    write(stdout,'(/,T6, A)')' This program is free software: you can redistribute it and/or modify it &
                      &under the terms of the'
    write(stdout,'(T5, A)')'GNU General Public License as published by the Free Software Foundation, &
                      &either version 3 of the'
    write(stdout,'(T30, A)')'License, or (at your option) any later version.'
    write(stdout,'(/,T5, A)')'This program is distributed in the hope that it will be useful, but &
                      &WITHOUT ANY WARRANTY; without'
    write(stdout,'(T12, A)')'even the implied warranty of MERCHANTABILITY or FITNESS FOR A &
                      &PARTICULAR PURPOSE.'
    write(stdout,'(T26, A)')'See the GNU General Public License for more details.'
    write(stdout,'(/,T4, A)')'You should have received a copy of the GNU General Public License along & 
                      &with this program. If not,'
    write(stdout,'(T34, A)')'see <https://www.gnu.org/licenses/>.'
    write(stdout,'(/,T36, A)')'E-mail: colombarifm@hotmail.com'
    write(stdout,'(/,T3, A,/)') dashline

    if ( allocated(arg) ) deallocate(arg)

    call err % termination(0,'f')

  end subroutine Display_license

  !---------------------------------------------------------------------------
  !> @brief Displays the version
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Display_version()

    implicit none
    
    character( len = : ), allocatable :: version
    
    version = '1.1.0'

    write(stdout,'("sas_grid ",a)') version
     
    call err % termination(0,'f')

  end subroutine Display_version

end module mod_cmd_line
