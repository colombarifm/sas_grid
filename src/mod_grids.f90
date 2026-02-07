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
!> @file   mod_grids.f90
!> @author Felippe M. Colombari
!> @email  colombarifm@hotmail.com
!> @brief  This module manipulates the coordinates of translation and rotation grids
!> @date - Nov, 2019                                                           
!> - independent module created                                                
!> @date - Nov, 2019                                                           
!> - spherical grids generation subroutine added
!> @date - Nov, 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!> @date - Jan, 2024
!> - add option to build a cube - sas grid (cavity within a cube)                    
!---------------------------------------------------------------------------------------------------

module mod_grids
  use mod_constants , only : DP

  implicit none

  private

  public :: Build_sphere_grid, Build_sas_grid, Build_box_grid, grid_sphere, grid_sas, grid_box, check 

  type point
    real( kind = DP )          :: grid_xyz(3)
    character( len = 3 )       :: grid_symbol
  end type point

  type sphere_grid
    type( point ), allocatable :: points_sphere(:)
    integer                    :: numpoints_sphere
  contains
    procedure, pass            :: Build_sphere_grid
  end type sphere_grid

  type( sphere_grid )          :: grid_sphere

  type sas_grid
    integer                    :: numpoints_sas
  contains
    procedure, pass            :: Build_sas_grid
  end type sas_grid

  type( sas_grid )             :: grid_sas

  type box_grid
    type( point ), allocatable :: points_box(:)
    integer                    :: numpoints_box
    integer                    :: numpoints_cavity
  contains
    procedure, pass            :: Build_box_grid
  end type box_grid

  type( box_grid )             :: grid_box
  
  real(kind=dp)                :: x1_sphere, y1_sphere, z1_sphere
  real(kind=dp)                :: x2_sphere, y2_sphere, z2_sphere
  real(kind=dp)                :: dx, dy, dz, r_sqr
  logical, allocatable         :: check(:,:)

contains

  !---------------------------------------------------------------------------
  !> @brief This routine generates the translation sphere by tesselation.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Build_sphere_grid( this, tessellation_factor )
    use mod_spherical_grids , only : sphere_icos1_points
    use mod_error_handling  , only : error

    implicit none

    class( sphere_grid ), intent(inout) :: this
    integer, intent(in)                 :: tessellation_factor
    integer                             :: node_num
    real( kind = DP ), allocatable      :: node_xyz(:,:)
    integer                             :: ierr
    type( error )                       :: err

    node_num = 12 + 10 * 3 * ( tessellation_factor - 1 ) + 10 * ( tessellation_factor - 2 ) * ( tessellation_factor - 1 )
    !edge_num = 30 * factor * factor
    !face_num = 20 * factor * factor

    this % numpoints_sphere = node_num

    allocate( node_xyz(3,node_num), stat=ierr )
    if(ierr/=0) call err % error('e',message="abnormal memory allocation")

    allocate( this % points_sphere( node_num ), stat=ierr )
    if(ierr/=0) call err % error('e',message="abnormal memory allocation")

    call sphere_icos1_points ( tessellation_factor, node_num, node_xyz )

    this % points_sphere(:) % grid_xyz(1) = node_xyz(1,:)
    this % points_sphere(:) % grid_xyz(2) = node_xyz(2,:)
    this % points_sphere(:) % grid_xyz(3) = node_xyz(3,:)

    if ( allocated(node_xyz) ) deallocate(node_xyz)
      
  end subroutine Build_sphere_grid

  subroutine Build_sas_grid( this )
    use mod_read_molecule  , only : mol
    use mod_error_handling , only : error
    use mod_constants      , only : PI

    implicit none

    class( sas_grid ), intent(inout) :: this
    character( len = 4 )             :: char_numpoints_sphere
    integer                          :: i, j, k, num_surf_points
    real( kind = DP )                :: atomic_sasa, total_sasa, total_vol, ratio_surf_points

    num_surf_points   = 0
    total_sasa        = 0.0_DP
    total_vol         = 0.0_DP
    ratio_surf_points = 0.0_DP

    allocate( check ( mol % num_atoms, grid_sphere % numpoints_sphere ) )

    check = .true.

    do i = 1, mol % num_atoms

      if ( ( mol % atoms(i) % label == 'X' ) .or. ( mol % atoms(i) % label == 'XX' ) ) cycle

      do j = 1, grid_sphere % numpoints_sphere

        x1_sphere = mol % atoms(i) % xyz(1) + &
                  & mol % atoms(i) % rvdw_solv * &
                  & grid_sphere % points_sphere(j) % grid_xyz(1)

        y1_sphere = mol % atoms(i) % xyz(2) + &
                  & mol % atoms(i) % rvdw_solv * &
                  & grid_sphere % points_sphere(j) % grid_xyz(2)

        z1_sphere = mol % atoms(i) % xyz(3) + &
                  & mol % atoms(i) % rvdw_solv * &
                  & grid_sphere % points_sphere(j) % grid_xyz(3)

        do k = 1, mol % num_atoms

          if ( i == k ) cycle

          x2_sphere = mol % atoms(k) % xyz(1)
          y2_sphere = mol % atoms(k) % xyz(2)       
          z2_sphere = mol % atoms(k) % xyz(3)
          
          dx = x2_sphere - x1_sphere     
          dy = y2_sphere - y1_sphere        
          dz = z2_sphere - z1_sphere
          
          r_sqr = dx*dx + dy*dy + dz*dz 
          
          if ( r_sqr < mol % atoms(k) % rvdw_solv ** 2 ) then
    
            check(i,j) = .false.
            
            exit ! atom is buried
              
          endif
            
        enddo
        
      enddo
      
    enddo

    grid_sas % numpoints_sas = 0

    do i = 1, mol % num_atoms
    
      if ( ( mol % atoms(i) % label == 'X' ) .or. ( mol % atoms(i) % label == 'XX' ) ) cycle

      do j = 1, grid_sphere % numpoints_sphere
          
        if ( check(i,j) .eqv. .true. ) then
        
          this % numpoints_sas = this % numpoints_sas + 1
          
        endif
        
      enddo
     
    enddo

    write(char_numpoints_sphere,'(i4.4)') grid_sphere % numpoints_sphere

    open( unit = 40, file = 'sas_'//char_numpoints_sphere//'.xyz', status = 'unknown' )
      
    write (40,*) this % numpoints_sas
    
    write (40,'(a)')
    
    write(*,*) '   ATOMIC SASA'
    write(*,*) ' -------------------------------------------'
    write(*,'(a8,3x,a5,3x,a10,3x,a10)') "atom", "npts", " SASA(A^2)", "SASA(nm^2)"

    do i = 1, mol % num_atoms

      if ( ( mol % atoms(i) % label == 'X' ) .or. ( mol % atoms(i) % label == 'XX' ) ) cycle

      do j = 1, grid_sphere % numpoints_sphere
      
        if ( check(i,j) .eqv. .true. ) then
        
          x1_sphere = mol % atoms(i) % xyz(1) + &
                    & mol % atoms(i) % rvdw_solv * &
                    & grid_sphere % points_sphere(j) % grid_xyz(1)

          y1_sphere = mol % atoms(i) % xyz(2) + &
                    & mol % atoms(i) % rvdw_solv * &
                    & grid_sphere % points_sphere(j) % grid_xyz(2)

          z1_sphere = mol % atoms(i) % xyz(3) + &
                    & mol % atoms(i) % rvdw_solv * &
                    & grid_sphere % points_sphere(j) % grid_xyz(3)

          write (40,'(a2,3(3x,f12.8))') 'XX', x1_sphere, y1_sphere, z1_sphere
          
        endif
        
      enddo

      num_surf_points = count(check(i,:))

      ratio_surf_points = dble(num_surf_points) / dble(grid_sphere % numpoints_sphere)

      atomic_sasa = ratio_surf_points * 4.0_DP * PI * (mol % atoms(i) % rvdw_solv ** 2 )
        
      write(*,'(i8,3x,i5,3x,f10.3,3x,f10.4)') i, num_surf_points, atomic_sasa, atomic_sasa/100

      total_sasa = total_sasa + atomic_sasa

    enddo
      
    close(40)

    write(*,*) ' -------------------------------------------'
    write(*,'(a13,f13.3,f13.4)') '   TOTAL SASA   ', total_sasa, total_sasa/100

    write(*,*) ' -------------------------------------------'
    write(*,'(a11,f10.3,a7)') 'Esurf = ', total_sasa * 0.02267, ' kJ/mol '

  end subroutine Build_sas_grid 

  subroutine Build_box_grid( this, box_min, box_max )
    use mod_read_molecule  , only : mol
    use mod_error_handling , only : error

    implicit none

    class( box_grid ), intent(inout) :: this
    real( kind = DP ), intent(in)    :: box_min(3), box_max(3)
    real( kind = DP ), allocatable   :: coords_box(:,:), coords_cavity(:,:)
    real( kind = DP )                :: rij, delta
    integer                          :: npoints_x, npoints_y, npoints_z
    integer                          :: i, ix, iy, iz

    delta = 0.5_DP

    npoints_x = dabs( box_min(1) - box_max(1) ) / delta
    npoints_y = dabs( box_min(2) - box_max(2) ) / delta
    npoints_z = dabs( box_min(3) - box_max(3) ) / delta

    this % numpoints_box = npoints_x * npoints_y * npoints_z

    allocate( coords_box( 3, this % numpoints_box ) )
    allocate( coords_cavity( 3, this % numpoints_box ) )

    coords_box    = 0.0_DP
    coords_cavity = 0.0_DP

    this % numpoints_cavity = 0

a:  do ix = 1, npoints_x

b:    do iy = 1, npoints_y

c:      do iz = 1, npoints_z

          coords_box( 1, ix ) = box_min(1) + ( ix * delta )
          coords_box( 2, iy ) = box_min(2) + ( iy * delta )
          coords_box( 3, iz ) = box_min(3) + ( iz * delta )
          
          do i = 1, mol % num_atoms

            rij = dsqrt( (coords_box(1,ix) - mol % atoms(i) % xyz(1))**2 + &
                         (coords_box(2,iy) - mol % atoms(i) % xyz(2))**2 + &
                         (coords_box(3,iz) - mol % atoms(i) % xyz(3))**2 ) 

            if ( rij < mol % atoms(i) % rvdw_solv ) then 

              cycle c

            endif

          enddo
    
          this % numpoints_cavity = this % numpoints_cavity + 1

          coords_cavity(1, this % numpoints_cavity) = coords_box(1,ix)
          coords_cavity(2, this % numpoints_cavity) = coords_box(2,iy)
          coords_cavity(3, this % numpoints_cavity) = coords_box(3,iz)

        enddo c

      enddo b

    enddo a

    open( unit = 987, file = 'cavity.xyz', status = 'unknown' )

    write(987,'(i5)') this % numpoints_cavity

    write(987,*)

    do i = 1, this % numpoints_cavity

      write(987,'(a2,3f12.6)') "XX", coords_cavity(:,i)

    enddo

    close(987)

  end subroutine Build_box_grid

end module mod_grids
