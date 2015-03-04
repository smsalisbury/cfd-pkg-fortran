program main
use config
use file_io
use navier_stokes_solvers
implicit none

!	-----------------------------------------------------
!	NAVIER-STOKES SOLVER
!	
!	This program solves the 2D Navier-Stokes equations for
!	specified 2D fluid field.
!
!	Author:	Spencer Salisbury
!	Date:	March 2015
!	-----------------------------------------------------

!	DATA DICTIONARY
real(wp)						::	density		! Constant density of the fluid

character(20)					::	bc_top_type,bc_left_type,bc_right_type,bc_bottom_type
character(20),dimension(4)		::	bc_type_names	=	(/'bc_top_type   ','bc_right_type ','bc_bottom_type','bc_left_type  '/)
character(20),dimension(4)		::	bc_types

character(20)					::	ic_type
real(wp)						::	ic_u,ic_v

integer							::	namelist_file,io_stat
character(20)					::	namelist_file_name

!	NAMELIST
namelist /FLUID_PROPERTIES/ density
namelist /BOUNDARY_CONDITIONS/ bc_top_type,bc_bottom_type,bc_left_type,bc_right_type

!	--	Get namelist from user
write(*,*)'Enter input file:'
read(*,*)namelist_file_name

!	--	Default namelist values

!	-- 	Read in namelists
namelist_file = file_open(adjustl(trim(namelist_file_name)))
read(namelist_file,FLUID_PROPERTIES,IOSTAT=io_stat)
read(namelist_file,BOUNDARY_CONDITIONS,IOSTAT=io_stat)
close(namelist_file)

!	--	Arrange namelist variables
bc_types = (/bc_top_type,bc_right_type,bc_bottom_type,bc_left_type/)
write(*,*)bc_types

!	SETUP AND OPEN OUTPUT FILES

end program