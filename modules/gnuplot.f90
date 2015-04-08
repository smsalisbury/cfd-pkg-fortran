module gnuplot
use config
use file_io
implicit none

character(5),parameter			::	plot_dir = "plots"

contains
	subroutine vector_plot2D(x,y,u,v,mag_in)
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
		real(wp)								::	x_mag, y_mag
		real(wp)								::	avg_dx,avg_dy
		
		real(wp)								::	mag
		
		integer,parameter						::	img_height = 480
		real(wp)								::	aspect_ratio !	x to y
		
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
		
		!	--	Set image aspect ratio
		aspect_ratio = maxval(x)/maxval(y)
		
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
				x_mag = mag*(u(i,j)/max_mag)*avg_dx
				y_mag = mag*(v(i,j)/max_mag)*avg_dy
				write(data_file,*)x(i),y(j),x_mag,y_mag,sqrt(x_mag**2 + y_mag**2)
			enddo
		enddo		
		close(data_file)
		
		!	WRITE GNUPLOT SCRIPT
		script_file = file_open(trim(plot_dir) // '/script.dat')
		write(script_file,*)"set terminal png size", int(img_height*aspect_ratio), "," , img_height
		write(script_file,*)"set output '" // trim(plot_dir) // "/vector.png'"
		write(script_file,*)"plot '" // trim(plot_dir) // "/vector.dat' using 1:2:3:4 with vectors head filled lt 2 notitle"
		close(script_file)
		
		!	PLOT GNUPLOT
		ret = system("gnuplot " // trim(plot_dir) // "/script.dat")
	
		!	CLEANUP
		! call rm_file(trim(plot_dir) // '/script.dat')
		call system('mv ' // trim(plot_dir) // '/script.dat ' // trim(plot_dir) // '/script2.dat')
		call system('mv ' // trim(plot_dir) // '/vector.dat ' // trim(plot_dir) // '/vector2.dat')
		call rm_file(trim(plot_dir) // '/vector.dat')
	end subroutine vector_plot2D

	subroutine contour_plot2D(x,y,values)
		!	--------------------------
		!	This subroutine uses gnuplot to plot
		!	a 2D contour plot. x and y are vectors
		!	of the x and y coordinates, and values is
		! 	a matrix of the values to plot.
		!	--------------------------
		
		!	INPUTS
		real(wp),dimension(:),intent(in)		::	x,y
		real(wp),dimension(:,:),intent(in)		::	values
		
		!	INTERNAL VARIABLES
		integer									::	values_x_size, values_y_size
		integer									::	x_size, y_size
		integer									::	data_file
		integer									::	script_file
		
		integer									::	i,j
		integer									::	ret
		
		real(wp)								::	max_values
		real(wp)								::	min_values
		
		integer,parameter						::	img_height = 480
		real(wp)								::	aspect_ratio !	x to y
		
		!	VALIDATION
		!	--	Sizes must match
		values_x_size = size(values,1)
		values_y_size = size(values,2)
		
		if (values_x_size /= size(x)) then
			call gnuplot_error(104)
		elseif (values_y_size /= size(y)) then
			call gnuplot_error(105)
		endif
		
		x_size = size(x)
		y_size = size(y)
		
		!	FIND MAX AND MIN values
		max_values = maxval(values)
		min_values = minval(values)
		
		!	--	Set image aspect ratio
		aspect_ratio = maxval(x)/maxval(y)
		
		!	WRITE DATA TO TEMPORARY FILE
		!	--	The format for gnuplot is x,y,deltax,deltay
		call create_dir(trim(plot_dir))
		data_file = file_open(trim(plot_dir) // '/contour.dat')
		do i=1,x_size
			do j=1,y_size
				write(data_file,*)x(i),y(j),values(i,j)
			enddo
		enddo		
		close(data_file)
		
		!	WRITE GNUPLOT SCRIPT
		script_file = file_open(trim(plot_dir) // '/script.dat')
		write(script_file,*)"set terminal png size", int(img_height*aspect_ratio), "," , img_height
		write(script_file,*)"set parametric"
		write(script_file,*)"set contour base"
		write(script_file,*)"set dgrid3d"
		write(script_file,*)"set view 0,0,1"
		write(script_file,*)"unset surface"
		write(script_file,*)"set cntrparam levels 10"
		write(script_file,*)"set output '" // trim(plot_dir) // "/contour.png'"
		write(script_file,*)"splot '" // trim(plot_dir) // "/contour.dat' using 1:2:3 with line notitle"
		close(script_file)
		
		!	PLOT GNUPLOT
		ret = system("gnuplot " // trim(plot_dir) // "/script.dat")
	
		!	CLEANUP
		call rm_file(trim(plot_dir) // '/script.dat')
		call rm_file(trim(plot_dir) // '/contour.dat')
	end subroutine contour_plot2D

	subroutine scatter_plot2D(x_values,y_values,aspect_ratio_in)
		!	--------------------------
		!	This subroutine uses gnuplot to plot
		!	a 2D scatter plot. x_values and y_values
		!	are vectors.
		!	--------------------------
		
		!	INPUTS
		real(wp),dimension(:),intent(in)		::	x_values,y_values
		real(wp),intent(in),optional			::	aspect_ratio_in
		
		!	INTERNAL VARIABLES
		integer									::	x_size, y_size
		integer									::	data_file
		integer									::	script_file
		
		integer									::	i,j
		integer									::	ret
		
		integer,parameter						::	img_height = 480
		real(wp)								::	aspect_ratio !	x to y
		
		!	DEFAULTS
		!	--	Aspect Ratio
		if (present(aspect_ratio_in)) then
			aspect_ratio = aspect_ratio_in
		else
			aspect_ratio = 1.0_wp
		endif
		
		!	VALIDATION
		!	--	Sizes must match		
		x_size = size(x_values)
		y_size = size(y_values)
		
		if (x_size /= y_size) then
			call gnuplot_error(106)
		endif
		
		!	WRITE DATA TO TEMPORARY FILE
		!	--	The format for gnuplot is x,y,deltax,deltay
		call create_dir(trim(plot_dir))
		data_file = file_open(trim(plot_dir) // '/scatter.dat')
		do i=1,x_size
			write(data_file,*)x_values(i),y_values(i)
		enddo		
		close(data_file)
		
		!	WRITE GNUPLOT SCRIPT
		script_file = file_open(trim(plot_dir) // '/script.dat')
		write(script_file,*)"set terminal png size", int(real(img_height,wp)*aspect_ratio), "," , img_height
		write(script_file,*)"set output '" // trim(plot_dir) // "/scatter.png'"
		write(script_file,*)"plot '" // trim(plot_dir) // "/scatter.dat' using 1:2 with line notitle"
		close(script_file)
		
		!	PLOT GNUPLOT
		ret = system("gnuplot " // trim(plot_dir) // "/script.dat")
	
		!	CLEANUP
		call rm_file(trim(plot_dir) // '/script.dat')
		call rm_file(trim(plot_dir) // '/scatter.dat')
	end subroutine scatter_plot2D
	
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
				write(*,*)"2D Vector Plot: Dimensions of u don't match with dimensions of the x vector."
				write(*,*)"Aborting..."
				stop
			case (101)
				write(*,*)"GNUPLOT MODULE ERROR:"
				write(*,*)"2D Vector Plot: Dimensions of u don't match with dimensions of the y vector."
				write(*,*)"Aborting..."
				stop
			case (102)
				write(*,*)"GNUPLOT MODULE ERROR:"
				write(*,*)"2D Vector Plot: Dimensions of v don't match with dimensions of the x vector."
				write(*,*)"Aborting..."
				stop
			case (103)
				write(*,*)"GNUPLOT MODULE ERROR:"
				write(*,*)"2D Vector Plot: Dimensions of v don't match with dimensions of the y vector."
				write(*,*)"Aborting..."
				stop
			case (104)
				write(*,*)"GNUPLOT MODULE ERROR:"
				write(*,*)"2D Contour Plot: Dimensions of values don't match with dimensions of the x vector."
				write(*,*)"Aborting..."
				stop
			case (105)
				write(*,*)"GNUPLOT MODULE ERROR:"
				write(*,*)"2D Contour Plot: Dimensions of values don't match with dimensions of the y vector."
				write(*,*)"Aborting..."
				stop
			case (106)
				write(*,*)"GNUPLOT MODULE ERROR:"
				write(*,*)"2D Scatter Plot: Dimensions of x and y don't match."
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