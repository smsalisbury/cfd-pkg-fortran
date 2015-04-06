program main
use config
use file_io
use navier_stokes_solvers
use gnuplot
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

character(20),dimension(4)			::	bc_list_def,bc_list
	
!	--	Initial Conditions	
character(20)						::	ic_type
real(wp)							::	ic_u,ic_v,ic_P
	
!	--	Mesh Properties	
integer								::	x_steps,y_steps
real(wp)							::	x_size,y_size
real(wp)							::	dx,dy
real(wp)							::	x_coord,y_coord

!	--	Iterative Properties
real(wp)							::	relax_mom,relax_press,relax_press_cor
real(wp)							::	conv_error
integer								::	max_iter

!	--	Things to solve for
real(wp),dimension(:,:),allocatable	::	u,v,P
real(wp),dimension(:,:),allocatable	::	mu,mv

real(wp),dimension(:,:),allocatable	::	AP_u,AP_v
real(wp),dimension(:,:),allocatable	::	continuity
real(wp)							::	cont_rms,cont_rms_old
integer								::	conv_satisfied = 0

!	--	Graphics arrays
real(wp),dimension(:,:),allocatable	::	u_g,v_g
real(wp),dimension(:),allocatable	::	x_g,y_g
real(wp),dimension(:),allocatable	::	continuity_data,iteration_data

!	--	Namelist File Stuff
integer								::	namelist_file,io_stat
character(20)						::	namelist_file_name
	
!	--	Output File Stuff	
character(20)						::	output_dir
integer								::	u_axes_file
integer								::	v_axes_file
integer								::	P_axes_file
integer								::	u_file
integer								::	v_file
integer								::	continuity_file

!	--	Counters
integer								::	i,j,k

!	--	Timers
integer								::	time_scale
integer								::	start_time,time_clicks

!	NAMELIST
namelist /FLUID_PROPERTIES/ density,viscosity
namelist /BOUNDARY_CONDITIONS/ bc_top_type,bc_bottom_type,bc_left_type,bc_right_type, &
	bc_top_u_velocity,bc_top_v_velocity,bc_right_u_velocity,bc_right_v_velocity, &
	bc_bottom_u_velocity,bc_bottom_v_velocity,bc_left_u_velocity,bc_left_v_velocity
namelist /INITIAL_CONDITIONS/ ic_type,ic_u,ic_v,ic_P
namelist /MESH_PROPERTIES/ x_steps,y_steps,x_size,y_size
namelist /ITERATIVE_PROPERTIES/ relax_mom,relax_press,relax_press_cor,conv_error,max_iter

!	START CLOCK
call system_clock(COUNT_RATE=time_scale,COUNT=start_time)

!	--	Get namelist from user
!write(*,*)'Enter input file:'
!read(*,*)namelist_file_name
namelist_file_name = 'inputs_channel'

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

relax_mom = 0.5_wp
relax_press = 0.3_wp
relax_press_cor = 1.2_wp
conv_error = 1.0E-5_wp
max_iter = 3000

!	-- 	Read in namelists
call system_clock(time_clicks)
write(*,"(F5.2,A,3X,A)")real(time_clicks-start_time)/real(time_scale),'s','Reading in namelist data.'
namelist_file = file_open(adjustl(trim(namelist_file_name)))
read(namelist_file,FLUID_PROPERTIES,IOSTAT=io_stat)
read(namelist_file,BOUNDARY_CONDITIONS,IOSTAT=io_stat)
read(namelist_file,INITIAL_CONDITIONS,IOSTAT=io_stat)
read(namelist_file,MESH_PROPERTIES,IOSTAT=io_stat)
read(namelist_file,ITERATIVE_PROPERTIES,IOSTAT=io_stat)
close(namelist_file)

!	SETUP AND OPEN OUTPUT FILES
call system_clock(time_clicks)
write(*,"(F5.2,A,3X,A)")real(time_clicks-start_time)/real(time_scale),'s','Opening output files.'
output_dir = 'outputs'
call create_dir(output_dir)
u_axes_file		= file_open(trim(output_dir) // '/u_axes.dat')
v_axes_file		= file_open(trim(output_dir) // '/v_axes.dat')
P_axes_file		= file_open(trim(output_dir) // '/P_axes.dat')
u_file			= file_open(trim(output_dir) // '/u.dat')
v_file			= file_open(trim(output_dir) // '/v.dat')
continuity_file	= file_open(trim(output_dir) // '/continuity.dat')

!	SETUP MESH
call system_clock(time_clicks)
write(*,"(F5.2,A,3X,A)")real(time_clicks-start_time)/real(time_scale),'s','Setting up mesh.'
dx = x_size/real(x_steps,wp)
dy = y_size/real(y_steps,wp)

!	--	Allocate u, v, P, and mass flow rates
allocate(u(0:x_steps+1,0:y_steps+1))
allocate(v(0:x_steps+1,0:y_steps+1))
allocate(P(1:x_steps,1:y_steps))
allocate(continuity(1:x_steps,1:y_steps))

allocate(AP_u(0:x_steps+1,0:y_steps+1))
allocate(AP_v(0:x_steps+1,0:y_steps+1))

allocate(mu(0:x_steps+1,0:y_steps+1))
allocate(mv(0:x_steps+1,0:y_steps+1))

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

!	--	Initialize u, v, and P and mass flow rates
call system_clock(time_clicks)
write(*,"(F5.2,A,3X,A)")real(time_clicks-start_time)/real(time_scale),'s','Initializing momentum, pressure, and flow rates.'
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

AP_u = 0.0_wp
AP_v = 0.0_wp

call system_clock(time_clicks)
write(*,"(F5.2,A,3X,A)")real(time_clicks-start_time)/real(time_scale),'s','Setting up Dirichlet boundaries.'
!	BOUNDARY CONDITIONS
!	--	Set Dirichlet Boundaries
!		--	Top boundary
select case (bc_top_type)
	case ('wall')
		u(:,y_steps+1) = 0.0_wp
		v(:,y_steps+1) = 0.0_wp
	case ('velocity')
		u(:,y_steps+1) = bc_top_u_velocity
		v(:,y_steps+1) = bc_top_v_velocity
	case ('outlet')
		u(:,y_steps+1) = 0.0_wp
		v(:,y_steps+1) = 0.0_wp
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
	case ('outlet')
		u(:,y_steps+1) = 0.0_wp
		v(:,y_steps+1) = 0.0_wp
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
	case ('outlet')
		u(:,y_steps+1) = 0.0_wp
		v(:,y_steps+1) = 0.0_wp
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
	case ('outlet')
		u(:,y_steps+1) = 0.0_wp
		v(:,y_steps+1) = 0.0_wp
	case default
		call prog_error(101)
endselect

!	-- Store Boundary Conditions in an array
bc_list_def(1) = 'top'
bc_list_def(2) = 'right'
bc_list_def(3) = 'bottom'
bc_list_def(4) = 'left'
bc_list(1) = bc_top_type
bc_list(2) = bc_right_type
bc_list(3) = bc_bottom_type
bc_list(4) = bc_left_type

!	BEGIN SOLUTION
call system_clock(time_clicks)
write(*,"(F5.2,A,3X,A)")real(time_clicks-start_time)/real(time_scale),'s','Solving navier stokes equations.'
cont_rms_old = 0.0_wp
do k=1,max_iter
	!	--	MASS FLOW RATES
	call mass_flow_rates(mu,mv,u,v,P,density,viscosity,dx,dy)
	!	--	U-MOMENTUM
	call umomentum(mu,mv,u,v,P,density,viscosity,dx,dy,relax_mom,bc_list,AP_u)
	!	--	V-MOMENTUM
	call vmomentum(mu,mv,u,v,P,density,viscosity,dx,dy,relax_mom,bc_list,AP_v)
	!	--	PRESSURES
	call pressure(mu,mv,u,v,P,density,viscosity,dx,dy,relax_press_cor,relax_press,AP_u,AP_v,continuity)
	
	cont_rms = rms(continuity)
	write(continuity_file,*)k,cont_rms,abs(cont_rms - cont_rms_old)
	write(*,*)k,cont_rms,abs(cont_rms - cont_rms_old)
	
	if (abs(cont_rms - cont_rms_old) < conv_error) then
		conv_satisfied = conv_satisfied + 1
	else
		conv_satisfied = 0
	endif
	
	if (conv_satisfied >= 3) exit
	
	cont_rms_old = cont_rms
enddo
call system_clock(time_clicks)
write(*,"(F5.2,A,3X,A,ES10.3)")real(time_clicks-start_time)/real(time_scale),'s','-->Continuity converged to ',cont_rms
if (k >= max_iter) k = max_iter
allocate(continuity_data(k))
allocate(iteration_data(k))

close(continuity_file)

!	PREPARE DATA FOR GRAPHICS
call system_clock(time_clicks)
write(*,"(F5.2,A,3X,A)")real(time_clicks-start_time)/real(time_scale),'s','Preparing data for graphics.'
allocate(u_g(x_steps,y_steps))
allocate(v_g(x_steps,y_steps))
allocate(x_g(x_steps))
allocate(y_g(y_steps))
u_g = 0.0_wp
v_g = 0.0_wp

do i=1,x_steps
	do j=1,y_steps
		u_g(i,j) = 0.5_wp*(u(i,j) + u(i+1,j))
		v_g(i,j) = 0.5_wp*(v(i,j) + v(i,j+1))
		x_g(i) = dx*real(i,wp)-0.5_wp*dx
		y_g(j) = dy*real(j,wp)-0.5_wp*dy
	enddo
enddo

do k=1,y_steps
	write(u_file,*)u_g(:,k)
	write(v_file,*)v_g(:,k)
enddo

continuity_file	= file_open(trim(output_dir) // '/continuity.dat')
do k=1,size(continuity_data)
	read(continuity_file,*)iteration_data(k),continuity_data(k)
enddo

!	PLOTS
call system_clock(time_clicks)
write(*,"(F5.2,A,3X,A)")real(time_clicks-start_time)/real(time_scale),'s','Graphing flow field.'
call vector_plot2D(x_g,y_g,u_g,v_g,3.0_wp)
call contour_plot2D(x_g,y_g,arrays_mag(u_g,v_g))

call system_clock(time_clicks)
write(*,"(F5.2,A,3X,A)")real(time_clicks-start_time)/real(time_scale),'s','Graphing convergence.'
call scatter_plot2D(iteration_data,continuity_data,2.0)

!	DEALLOCATE ARRAYS
deallocate(u)
deallocate(v)
deallocate(P)

deallocate(u_g)
deallocate(v_g)
deallocate(x_g)
deallocate(y_g)

deallocate(mu)
deallocate(mv)

deallocate(AP_u)
deallocate(AP_v)

deallocate(continuity)

!	CLOSE OUTPUT FILES
close(u_axes_file)
close(v_axes_file)
close(P_axes_file)
close(u_file)
close(v_file)
close(continuity_file)

call system_clock(time_clicks)
write(*,"(F5.2,A,3X,A)")real(time_clicks-start_time)/real(time_scale),'s','Finished'

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
	end subroutine prog_error
	
	function rms(A)
		!	---------------------------
		!	Calculates the rms of an array
		!	---------------------------
		real(wp),dimension(:,:)		::	A
		real(wp)					::	rms
		
		integer						::	i,j,i_size,j_size
		
		i_size = size(A,1)
		j_size = size(A,2)
		rms = 0.0_wp
		do i=1,i_size
			do j=1,j_size
				rms = rms + A(i,j)**2
			enddo
		enddo
		rms = sqrt(rms)/real(i_size*j_size,wp)
	end function rms
	
	function arrays_mag(u,v)
		!	---------------------------
		!	Calculates the magnitude of two u and v arrays
		!	---------------------------
		real(wp),dimension(:,:)		::	u,v
		real(wp),dimension(:,:),allocatable		::	arrays_mag
		
		integer						::	i,j,i_size,j_size
		
		i_size = size(u,1)
		j_size = size(u,2)
		allocate(arrays_mag(i_size,j_size))
		arrays_mag = 0.0_wp
		do i=1,i_size
			do j=1,j_size
				arrays_mag(i,j) = sqrt(u(i,j)**2 + v(i,j)**2)
			enddo
		enddo
	end function arrays_mag
	
end program main
