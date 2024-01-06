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
!> @file   mod_constants.f90
!> @author Felippe M. Colombari
!> @brief  Defines a set constants.
!> @date - Oct, 2019                                                           
!> - independent module created                                                
!> @date - Sep, 2020
!> - vdW radii values added
!---------------------------------------------------------------------------------------------------

module mod_constants

  implicit none

  integer, public, parameter           :: DP = selected_real_kind(15, 307) !  double precision constant for portability
  
  real( kind = DP ), public, parameter :: PI      = 3.14159265358979_DP    !         PI constant with 14 decimal places
  real( kind = DP ), public, parameter :: DEG2RAD = 180.0_DP / PI          !                         degrees to radians
  real( kind = DP ), public, parameter :: FPZERO  = tiny(1.0_DP)           !              define machine-precision ZERO
  real( kind = DP ), public, parameter :: FPINF   = huge(1.0_DP)           !          define machine-precision INFINITY
  
  character( len = 11 ), public, parameter   :: INT_ALPHABET   = '1234567890'    !       allowed character for integers
  character( len = 12 ), public, parameter   :: FLOAT_ALPHABET = '.-1234567890'  !         allowed character for floats
  character( len = 66 ), public, parameter   :: CHAR_ALPHABET  = &
    'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ._-1234567890 '   !              allowed character for strings
  character( len = 100 ), public, parameter  :: DASHLINE = repeat('-',100)      !                       just a dashline

  !https://pubs.acs.org/doi/10.1021/jp8111556
  real( kind = DP ), parameter          :: radii_X  =   0.0000_DP
  real( kind = DP ), parameter          :: radii_XX =   0.0000_DP
  real( kind = DP ), parameter          :: radii_H  =   1.1000_DP
  real( kind = DP ), parameter          :: radii_He =   1.4000_DP
  real( kind = DP ), parameter          :: radii_Li =   1.8100_DP
  real( kind = DP ), parameter          :: radii_Be =   1.5300_DP
  real( kind = DP ), parameter          :: radii_B  =   1.9200_DP
  real( kind = DP ), parameter          :: radii_C  =   1.7000_DP
  real( kind = DP ), parameter          :: radii_N  =   1.5500_DP
  real( kind = DP ), parameter          :: radii_O  =   1.5200_DP
  real( kind = DP ), parameter          :: radii_F  =   1.4700_DP
  real( kind = DP ), parameter          :: radii_Ne =   1.5400_DP
  real( kind = DP ), parameter          :: radii_Na =   2.2700_DP
  real( kind = DP ), parameter          :: radii_Mg =   1.7300_DP
  real( kind = DP ), parameter          :: radii_Al =   1.8400_DP
  real( kind = DP ), parameter          :: radii_Si =   2.1000_DP
  real( kind = DP ), parameter          :: radii_P  =   1.8000_DP
  real( kind = DP ), parameter          :: radii_S  =   1.8000_DP
  real( kind = DP ), parameter          :: radii_Cl =   1.7500_DP
  real( kind = DP ), parameter          :: radii_Ar =   1.8800_DP
  real( kind = DP ), parameter          :: radii_K  =   2.7500_DP
  real( kind = DP ), parameter          :: radii_Ca =   2.3100_DP
  real( kind = DP ), parameter          :: radii_Ga =   1.8700_DP
  real( kind = DP ), parameter          :: radii_Ge =   2.1100_DP
  real( kind = DP ), parameter          :: radii_As =   1.8500_DP
  real( kind = DP ), parameter          :: radii_Se =   1.9000_DP
  real( kind = DP ), parameter          :: radii_Br =   1.8300_DP
  real( kind = DP ), parameter          :: radii_Kr =   2.0200_DP
  real( kind = DP ), parameter          :: radii_Rb =   3.0300_DP
  real( kind = DP ), parameter          :: radii_Sr =   2.4900_DP
  real( kind = DP ), parameter          :: radii_In =   1.9300_DP
  real( kind = DP ), parameter          :: radii_Sn =   2.1700_DP
  real( kind = DP ), parameter          :: radii_Sb =   2.0600_DP
  real( kind = DP ), parameter          :: radii_Te =   2.0600_DP
  real( kind = DP ), parameter          :: radii_I  =   1.9800_DP
  real( kind = DP ), parameter          :: radii_Xe =   2.1600_DP
  real( kind = DP ), parameter          :: radii_Cs =   3.4300_DP
  real( kind = DP ), parameter          :: radii_Ba =   2.6800_DP
  real( kind = DP ), parameter          :: radii_Tl =   1.9600_DP
  real( kind = DP ), parameter          :: radii_Pb =   2.0200_DP
  real( kind = DP ), parameter          :: radii_Bi =   2.0700_DP
  real( kind = DP ), parameter          :: radii_Po =   1.9700_DP
  real( kind = DP ), parameter          :: radii_At =   2.0200_DP
  real( kind = DP ), parameter          :: radii_Rn =   2.2000_DP
  real( kind = DP ), parameter          :: radii_Fr =   3.4800_DP
  real( kind = DP ), parameter          :: radii_Ra =   2.8300_DP

  real( kind = DP ), parameter          :: radii_W  =   2.1000_DP ! pubchem
  real( kind = DP ), parameter          :: radii_Zn =   1.3900_DP
  real( kind = DP ), parameter          :: radii_Cu =   1.4000_DP
  real( kind = DP ), parameter          :: radii_Pt =   1.7500_DP
  real( kind = DP ), parameter          :: radii_Pd =   1.6300_DP
  real( kind = DP ), parameter          :: radii_Ag =   1.7200_DP
  real( kind = DP ), parameter          :: radii_Au =   1.6600_DP
  real( kind = DP ), parameter          :: radii_Mo =   2.0900_DP
  real( kind = DP ), parameter          :: radii_Cd =   2.0900_DP

  ! todo

end module mod_constants
