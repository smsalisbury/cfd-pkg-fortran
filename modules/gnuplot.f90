module gnuplot
use config
use file_io
implicit none

character(5),parameter			::	plot_dir = "plots"

contains
	subroutine vector_plot(x,y,u,v,mag_in)
		!	--------------------------
		!	This subroutine uses gnuplot to plot
		!	a 2D vector plot. x and y are vectors
		!	of the x and y coordinates, and u and v
		! 	are matricies of the magnitudes
		!	--------------------------
		
		!	INPUTS
		real(wp),dimension(:),intent(in)		::	x,y
		real(wp),dimension(:,:),intent(in)		::	u,v
		real(wp),intent(in),optional			::	mag_in
		
		!	INTERNAL VARIABLES
		integer									::	u_x_size, u_y_size, x_size
		integer									::	v_x_size, v_y_size, y_size
		integer									::	data_file
		integer									::	script_file
		
		integer									::	i,j
		integer									::	ret
		
		real(wp)								::	max_u_mag, max_v_mag, max_mag
		real(wp)								::	avg_dx,avg_dy
		
		real(wp)								::	mag	
		
		!	VALIDATION
		!	--	Sizes must match
		u_x_size = size(u,1)
		u_y_size = size(u,2)
		v_x_size = size(v,1)
		v_y_size = size(v,2)
		
		if (u_x_size /= size(x)) then
			call gnuplot_error(100)
		elseif (u_y_size /= size(y)) then
			call gnuplot_error(101)
		elseif (v_x_size /= size(x)) then
			call gnuplot_error(102)
		elseif (v_y_size /= size(y)) then
			call gnuplot_error(103)
		endif
		
		x_size = size(x)
		y_size = size(y)
		
		!	--	Check and set magnification factor
		if (present(mag_in)) then
			mag = mag_in
		else
			mag = 1.0_wp
		endif
		
		!	FIND MAX u AND v
		max_u_mag = maxval(u)
		max_v_mag = maxval(v)
		max_mag = max(max_u_mag,max_v_mag)
		
		!	FIND avg_dx AND avg_dy
		avg_dx = 0.0_wp
		do i=1,x_size-1
			avg_dx = avg_dx + abs(x(i+1)-x(i))
		enddo
		avg_dx = avg_dx/real(x_size-1,wp)
		
		avg_dy = 0.0_wp
		do i=1,y_size-1
			avg_dy = avg_dy + abs(y(i+1)-y(i))
		enddo
		avg_dy = avg_dy/real(y_size-1,wp)
		
		!	WRITE DATA TO TEMPORARY FILE
		!	--	The format for gnuplot is x,y,deltax,deltay
		call create_dir(trim(plot_dir))
		data_file = file_open(trim(plot_dir) // '/vector.dat')
		do i=1,x_size
			do j=1,y_size
				write(data_file,*)x(i),y(j),mag*(u(i,j)/max_mag)*avg_dx,mag*(v(i,j)/max_mag)*avg_dx
			enddo
		enddo		
		close(data_file)
		
		!	WRITE GNUPLOT SCRIPT
		script_file = file_open(trim(plot_dir) // '/script.dat')
		write(script_file,*)"set terminal png"
		write(script_file,*)"set output '" // trim(plot_dir) // "/vector.png'"
		write(script_file,*)"plot '" // trim(plot_dir) // "/vector.dat' using 1:2:3:4 with vectors head filled lt 2 notitle"
		close(script_file)
		
		!	PLOT GNUPLOT
		ret = system("gnuplot " // trim(plot_dir) // "/script.dat")
	
		!	CLEANUP
		call rm_file(trim(plot_dir) // '/script.dat')
		call rm_file(trim(plot_dir) // '/vector.dat')
	end subroutine vector_plot

	subroutine gnuplot_error(code)
		!	----------------------------
		!	This subroutine handles the errors
		!	for the gnuplot module.
		!	----------------------------

		!	INPUTS
		integer			::	code
		
		!	HANDLE ERRORS
		select case (code)
			case (100)
				write(*,*)"GNUPLOT MODULE ERROR:"
				write(*,*)"Dimensions of u don't match with dimensions of the x vector."
				write(*,*)"Aborting..."
				stop
			case (101)
				write(*,*)"GNUPLOT MODULE ERROR:"
				write(*,*)"Dimensions of u don't match with dimensions of the y vector."
				write(*,*)"Aborting..."
				stop
			case (102)
				write(*,*)"GNUPLOT MODULE ERROR:"
				write(*,*)"Dimensions of v don't match with dimensions of the x vector."
				write(*,*)"Aborting..."
				stop
			case (103)
				write(*,*)"GNUPLOT MODULE ERROR:"
				write(*,*)"Dimensions of v don't match with dimensions of the y vector."
				write(*,*)"Aborting..."
				stop
			case default
				write(*,*)"GNUPLOT MODULE ERROR:"
				write(*,*)"There was an unknown gnuplot error."
				write(*,*)"Aborting..."
				stop
		endselect
	end subroutine gnuplot_error
end module gnuplot