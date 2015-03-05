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
real(wp)							::	density,viscosity
	
!	--	Boundary Conditions	
character(20)						::	bc_top_type,bc_left_type,bc_right_type,bc_bottom_type
real(wp)							::	bc_top_u_velocity,bc_top_v_velocity,bc_right_u_velocity,bc_right_v_velocity
real(wp)							::	bc_bottom_u_velocity,bc_bottom_v_velocity,bc_left_u_velocity,bc_left_v_velocity
	
!	--	Initial Conditions	
character(20)						::	ic_type
real(wp)							::	ic_u,ic_v,ic_P
	
!	--	Mesh Properties	
integer								::	x_steps,y_steps
real(wp)							::	x_size,y_size
real(wp)							::	dx,dy
real(wp)							::	x_coord,y_coord

!	--	Iterative Properties
real(wp)							::	relax_mom,relax_press

!	--	Things to solve for
real(wp),dimension(:,:),allocatable	::	u,v,P
real(wp),dimension(:,:),allocatable	::	mu,mv

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
namelist /FLUID_PROPERTIES/ density,viscosity
namelist /BOUNDARY_CONDITIONS/ bc_top_type,bc_bottom_type,bc_left_type,bc_right_type, &
	bc_top_u_velocity,bc_top_v_velocity,bc_right_u_velocity,bc_right_v_velocity, &
	bc_bottom_u_velocity,bc_bottom_v_velocity,bc_left_u_velocity,bc_left_v_velocity
namelist /INITIAL_CONDITIONS/ ic_type,ic_u,ic_v,ic_P
namelist /MESH_PROPERTIES/ x_steps,y_steps,x_size,y_size
namelist /ITERATIVE_PROPERTIES/ relax_mom,relax_press

!	--	Get namelist from user
!write(*,*)'Enter input file:'
!read(*,*)namelist_file_name
namelist_file_name = 'inputs'

!	--	Default namelist values
ic_type = 'uniform'

bc_top_type = 'wall'
bc_right_type = 'wall'
bc_bottom_type = 'wall'
bc_left_type = 'wall'

bc_top_u_velocity = 0.0_wp
bc_top_v_velocity = 0.0_wp
bc_right_u_velocity = 0.0_wp
bc_right_v_velocity = 0.0_wp
bc_bottom_u_velocity = 0.0_wp
bc_bottom_v_velocity = 0.0_wp
bc_left_u_velocity = 0.0_wp
bc_left_v_velocity = 0.0_wp

relax_mom = 1.0_wp
relax_press = 1.2_wp

!	-- 	Read in namelists
namelist_file = file_open(adjustl(trim(namelist_file_name)))
read(namelist_file,FLUID_PROPERTIES,IOSTAT=io_stat)
read(namelist_file,BOUNDARY_CONDITIONS,IOSTAT=io_stat)
read(namelist_file,INITIAL_CONDITIONS,IOSTAT=io_stat)
read(namelist_file,MESH_PROPERTIES,IOSTAT=io_stat)
read(namelist_file,ITERATIVE_PROPERTIES,IOSTAT=io_stat)
close(namelist_file)

!	SETUP AND OPEN OUTPUT FILES
output_dir = 'outputs'
call create_dir(output_dir)
u_axes_file		= file_open(trim(output_dir) // '/u_axes.out')
v_axes_file		= file_open(trim(output_dir) // '/v_axes.out')
P_axes_file		= file_open(trim(output_dir) // '/P_axes.out')

!	SETUP MESH
dx = x_size/real(x_steps,wp)
dy = y_size/real(y_steps,wp)

!	--	Allocate u, v, P, and mass flow rates
allocate(u(0:x_steps+1,0:y_steps+1))
allocate(v(0:x_steps+1,0:y_steps+1))
allocate(P(1:x_steps,1:y_steps))

allocate(mu(0:x_steps+1,0:y_steps+1))
allocate(mv(0:x_steps+1,0:y_steps+1))

!	--	Initialize u, v, and P and mass flow rates
select case (ic_type)
	case ("uniform")
		u = ic_u
		v = ic_v
		P = ic_P
	case default
		call prog_error(100)
endselect
mu = 0.0_wp
mv = 0.0_wp

!	--	Set Dirichlet Boundaries
!		--	Top boundary
select case (bc_top_type)
	case ('wall')
		u(:,y_steps+1) = 0.0_wp
		v(:,y_steps+1) = 0.0_wp
	case ('velocity')
		u(:,y_steps+1) = bc_top_u_velocity
		v(:,y_steps+1) = bc_top_v_velocity
	case default
		call prog_error(101)
endselect
!		--	Right boundary
select case (bc_right_type)
	case ('wall')
		u(x_steps+1,:) = 0.0_wp
		v(x_steps+1,:) = 0.0_wp
	case ('velocity')
		u(x_steps+1,:) = bc_right_u_velocity
		v(x_steps+1,:) = bc_right_v_velocity
	case default
		call prog_error(101)
endselect
!		--	Bottom boundary
select case (bc_bottom_type)
	case ('wall')
		u(:,0) = 0.0_wp
		v(:,0) = 0.0_wp
		v(:,1) = 0.0_wp
	case ('velocity')
		u(:,0) = bc_bottom_u_velocity
		v(:,0) = bc_bottom_v_velocity
		v(:,1) = bc_bottom_v_velocity
	case default
		call prog_error(101)
endselect
!		--	Left boundary
select case (bc_left_type)
	case ('wall')
		u(0,:) = 0.0_wp
		u(1,:) = 0.0_wp
		v(0,:) = 0.0_wp
	case ('velocity')
		u(0,:) = bc_left_u_velocity
		u(1,:) = bc_left_u_velocity
		v(0,:) = bc_left_v_velocity
	case default
		call prog_error(101)
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

!	BEGIN SOLUTION
!	--	MASS FLOW RATES
call mass_flow_rates(mu,mv,u,v,P,density,viscosity,dx,dy)
!	--	U-MOMENTUM
call umomentum(mu,mv,u,v,P,density,viscosity,dx,dy,relax_mom)
!	--	V-MOMENTUM
call vmomentum(mu,mv,u,v,P,density,viscosity,dx,dy,relax_mom)


call write_array(mu,'mu=')
call write_array(mv,'mv=')

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
		!	---------------------------
		!	A quick subroutine to display arrays
		!	---------------------------
		real(wp),intent(in),dimension(:,:)	::	A
		integer	::	k
		character(*),optional	::	annotation
		
		if (present(annotation)) then
			write(*,*)annotation
		endif
		
		do k=size(A,2),1,-1
			write(*,*)k,A(:,k)
		enddo
	end subroutine write_array
	
	subroutine prog_error(code)
		!	---------------------------
		!	A quick subroutine to handle errors
		!	---------------------------
		integer,intent(in)		::	code
		
		select case (code)
			case (100)	!	Initial Conditions
				write(*,*)'ERROR: FAILED TO INITIALIZE PROPERLY.'
				write(*,*)'Check your initial conditions in your inputs and try again.'
				write(*,*)'Aborting...'
				stop
			case (101)	!	Boundary Conditions
				write(*,*)'ERROR: FAILED TO SET BOUNDARIES PROPERLY.'
				write(*,*)'Check your boundary conditions in your inputs and try again.'
				write(*,*)'Aborting...'
				stop
			case default
				write(*,*)'ERROR: SOMETHING WENT WRONG.'
				write(*,*)"Unfortunately, we don't know what."
				write(*,*)'Aborting...'
		endselect
	end subroutine

end program