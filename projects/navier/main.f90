program main
use config
use file_io
use navier_stokes_solvers
implicit none

!	-----------------------------------------------------
!	NAVIER-STOKES SOLVER
!	
!	This program solves the 2D Navier-Stokes equations for
!	a specified 2D fluid field.
!
!	Author:	Spencer Salisbury
!	Date:	March 2015
!	-----------------------------------------------------

!	DATA DICTIONARY
!	--	Fluid Properties
real(wp)							::	density		! Constant density of the fluid
	
!	--	Boundary Conditions	
character(20)						::	bc_top_type,bc_left_type,bc_right_type,bc_bottom_type
character(20),dimension(4)			::	bc_type_names	=	(/'bc_top_type   ','bc_right_type ','bc_bottom_type','bc_left_type  '/)
character(20),dimension(4)			::	bc_types
	
!	--	Initial Conditions	
character(20)						::	ic_type
real(wp)							::	ic_u,ic_v,ic_P
	
!	--	Mesh Characteristics	
integer								::	x_steps,y_steps
real(wp)							::	x_size,y_size
real(wp)							::	dx,dy
real(wp)							::	x_coord,y_coord

!	--	Things to solve for
real(wp),dimension(:,:),allocatable	::	u,v,P

!	--	Namelist File Stuff
integer								::	namelist_file,io_stat
character(20)						::	namelist_file_name
	
!	--	Output File Stuff	
character(20)						::	output_dir
integer								::	u_axes_file
integer								::	v_axes_file
integer								::	P_axes_file

!	--	Counters
integer								::	k

!	NAMELIST
namelist /FLUID_PROPERTIES/ density
namelist /BOUNDARY_CONDITIONS/ bc_top_type,bc_bottom_type,bc_left_type,bc_right_type
namelist /INITIAL_CONDITIONS/ ic_type,ic_u,ic_v,ic_P
namelist /MESH_PROPERTIES/ x_steps,y_steps,x_size,y_size

!	--	Get namelist from user
!write(*,*)'Enter input file:'
!read(*,*)namelist_file_name
namelist_file_name = 'inputs'

!	--	Default namelist values

!	-- 	Read in namelists
namelist_file = file_open(adjustl(trim(namelist_file_name)))
read(namelist_file,FLUID_PROPERTIES,IOSTAT=io_stat)
read(namelist_file,BOUNDARY_CONDITIONS,IOSTAT=io_stat)
read(namelist_file,INITIAL_CONDITIONS,IOSTAT=io_stat)
read(namelist_file,MESH_PROPERTIES,IOSTAT=io_stat)
close(namelist_file)

!	--	Arrange namelist variables
bc_types = (/bc_top_type,bc_right_type,bc_bottom_type,bc_left_type/)

!	SETUP AND OPEN OUTPUT FILES
output_dir = 'outputs'
call create_dir(output_dir)
u_axes_file		= file_open(trim(output_dir) // '/u_axes.out')
v_axes_file		= file_open(trim(output_dir) // '/v_axes.out')
P_axes_file		= file_open(trim(output_dir) // '/P_axes.out')

!	SETUP MESH
dx = x_size/real(x_steps,wp)
dy = y_size/real(y_steps,wp)

!	--	Allocate u, v, and P
allocate(u(1:x_steps+1,0:x_steps+1))
allocate(v(0:y_steps+1,1:y_steps+1))
allocate(P(1:x_steps,1:y_steps))

!	--	Initialize u, v, and P
select case (ic_type)
	case ("uniform")
		u = ic_u
		v = ic_v
		P = ic_P
	case default
		call prog_error(100)
endselect

!	--	Write axes to files for reference
!		--	Pressure axes
do k=1,max(x_steps,y_steps)
	x_coord = dx*real(k,wp)-0.5_wp*dx
	y_coord = dy*real(k,wp)-0.5_wp*dy
	if (k > x_steps) then
		write(P_axes_file,*)0.0_wp,y_coord
	else if (k > y_steps) then
		write(P_axes_file,*)x_coord,0.0_wp
	else if (k <= min(x_steps,y_steps)) then
		write(P_axes_file,*)x_coord,y_coord
	endif
enddo

!		--	u axes
do k=1,(max(x_steps,y_steps)+1)
	x_coord = dx*real(k-1,wp)
	y_coord = dy*real(k,wp) - 0.5_wp*dy
	if (k > x_steps+1) then
		write(u_axes_file,*)0.0_wp,y_coord
	else if (k > y_steps) then
		write(u_axes_file,*)x_coord,0.0_wp
	else if (k <= min(x_steps,y_steps)+1) then
		write(u_axes_file,*)x_coord,y_coord
	endif
enddo

!		--	v axes
do k=1,(max(x_steps,y_steps)+1)
	x_coord = dx*real(k,wp) - 0.5_wp*dx
	y_coord = dy*real(k-1,wp)
	if (k > x_steps) then
		write(v_axes_file,*)0.0_wp,y_coord
	else if (k > y_steps+1) then
		write(v_axes_file,*)x_coord,0.0_wp
	else if (k <= min(x_steps,y_steps)+1) then
		write(v_axes_file,*)x_coord,y_coord
	endif
enddo

call write_array(P,'P=')
call write_array(u,'u=')
call write_array(v,'v=')

!	DEALLOCATE ARRAYS
deallocate(u)
deallocate(v)
deallocate(P)

!	CLOSE OUTPUT FILES
close(u_axes_file)
close(v_axes_file)
close(P_axes_file)

contains
	subroutine write_array(A,annotation)
		real(wp),intent(in),dimension(:,:)	::	A
		integer	::	k
		character(*),optional	::	annotation
		
		if (present(annotation)) then
			write(*,*)annotation
		endif
		do k=size(A,2),1,-1
			write(*,*)A(:,k)
		enddo
	end subroutine write_array
	
	subroutine prog_error(code)
		integer,intent(in)		::	code
		
		select case (code)
			case (100)	!	Initial Conditions
				write(*,*)'ERROR: FAILED TO INITIALIZE PROPERLY.'
				write(*,*)'Check your initial conditions in your inputs and try again.'
				write(*,*)'Aborting...'
				stop
			case default
				write(*,*)'ERROR: SOMETHING WENT WRONG.'
				write(*,*)"Unfortunately, we don't know what."
				write(*,*)'Aborting...'
		endselect
	end subroutine

end program