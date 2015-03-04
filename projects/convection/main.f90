program main
! 	2D CONVECTION UPWINDING
!
! 	Author: Spencer Salisbury
!     Date: February 2015
!
!	This code demonstrates two dimensional upwinding. It was written
!	for an MAE 5440 class at Utah State University.
use config
use file_io
implicit none

! DATA DICTIONARY
! -- Inputs
integer										::	n	! The number of cells in one direction
real(wp)									::	b	! The beta value for deferred correction

! -- Parameters						
real(wp),parameter							::	density = 1.0_wp				! Density
real(wp),parameter							::	u = 2.0_wp, v=2.0_wp			! Velocities in the x and y direction
real(wp),parameter							::	xsize = 1.0_wp, ysize = 1.0_wp	! Size of the mesh in each direction
integer,parameter							::	max_iter = 10000				! The maximum number of iterations
real(wp),parameter							::	a_error = 1.0E-7_wp				! The acceptable RSS error to stop at
character(40),parameter 					::  output_dir = 'output'			! The output directory

! -- Intermediate values
real(wp)									::	dx, dy									! The calculated size of the cells in both directions
real(wp)									::	rss,rss_old,drss,phi_old				! Error storage and comparison variables
real(wp),allocatable,dimension(:,:)			::	a_p1, a_n1, a_s1, a_e1, a_w1			! First order upwinding coefficients
real(wp),allocatable,dimension(:,:)			::	a_p2, a_n2, a_s2, a_e2, a_w2			! Second order upwinding coefficients
real(wp),allocatable,dimension(:,:)			::	F_e, F_w, F_n, F_s						! Flux
integer										::	i,j,k									! Counters
integer                 					::  convergence_file,phi_file,axes_file		! File units

! -- Calculated values
real(wp),allocatable,dimension(:,:)			::	phi		! The fluid scalar to be upwinded


! GET VALUES FROM USER
write(*,*)'Enter number of cells:'
read(*,*)n

write(*,*)'Enter a beta value (between 0 and 1):'
read(*,*)b

! VALIDATION
! -- Cell number
if (n <= 0) then
	write(*,*)'The number of cells must be greater than 0!'
	write(*,*)'Aborting...'
	stop
end if
! -- Beta value
if (b < 0 .OR. b > 1) then
	write(*,*)'The beta value must be between 0 and 1!'
	write(*,*)'Aborting...'
	stop
endif

! SET UP OUTPUT FILES
call create_dir(output_dir)
convergence_file = file_open(trim(output_dir) // '/conv.out')
phi_file         = file_open(trim(output_dir) // '/phi.out')
axes_file        = file_open(trim(output_dir) // '/axes.out')

! MESH SIZE CALCULATIONS
dx = xsize/real(n,wp)
dy = ysize/real(n,wp)

! ALLOCATE AND INITIALIZE ARRAYS
allocate(a_p1(n,n))
allocate(a_e1(n,n))
allocate(a_w1(n,n))
allocate(a_n1(n,n))
allocate(a_s1(n,n))
allocate(a_p2(n,n))
allocate(a_e2(n,n))
allocate(a_w2(n,n))
allocate(a_n2(n,n))
allocate(a_s2(n,n))
allocate(F_e(n,n))
allocate(F_w(n,n))
allocate(F_n(n,n))
allocate(F_s(n,n))
allocate(phi(0:n+1,0:n+1))

a_p1 = 0.0_wp
a_e1 = 0.0_wp
a_w1 = 0.0_wp
a_n1 = 0.0_wp
a_s1 = 0.0_wp
a_p2 = 0.0_wp
a_e2 = 0.0_wp
a_w2 = 0.0_wp
a_n2 = 0.0_wp
a_s2 = 0.0_wp
F_e = density*u
F_w = density*u
F_n = density*v
F_s = density*v
phi = 0.0_wp

! SET CONSTANTS
a_w1 = F_w
a_s1 = F_s
a_p1 = a_w1 + a_e1 + a_s1 + a_n1

a_e2 = -F_e*0.5_wp
a_n2 = -F_n*0.5_wp
a_w2 = F_w*0.5_wp
a_s2 = F_s*0.5_wp

! DIRICHLET BOUNDARIES
do i=0,n+1
	phi(i,0) = 0.0_wp
end do
do j=0,n+1
	phi(0,j) = 100.0_wp
end do
phi(0,0) = 50.0_wp

do i=1,n
	a_w2(1,i) = F_w(1,i)
	a_e2(n,i) = F_e(n,i)
enddo
do j=1,n
	a_s2(j,1) = F_s(j,1)
	a_n2(j,n) = F_n(j,n)
enddo
a_p2 = a_w2 + a_e2 + a_s2 + a_n2

! CALCULATIONS
rss = 0.0_wp
do k=1,max_iter
	rss_old = rss
	rss = 0.0_wp
	do i=1,n
		do j=1,n
			phi_old = phi(i,j)
			phi(i,j) = upwind(b)
			if (i == n) then
				phi(i+1,j) = phi(i,j)
			end if
			if (j == n) then
				phi(i,j+1) = phi(i,j)
			end if
			rss = rss + (phi_old - phi(i,j))**2
		end do
	end do
	rss = sqrt(rss)
	drss = abs(rss-rss_old)
	write(convergence_file,*)k,drss
	write(*,*)k,drss
	if (drss < a_error) then
		exit
	end if
end do
phi(n+1,n+1) = 0.5_wp*(phi(n,n+1)+phi(n+1,n))

do j=0,n+1
	write(phi_file,*)phi(j,:)
	if (j == 0) then
		write(axes_file,*)0,0
	else if (j == n+1) then
		write(axes_file,*)xsize,ysize
	else
		write(axes_file,*)(2.0_wp*real(j)-1)*dx/2.0_wp,(2.0_wp*real(j)-1)*dy/2.0_wp
	end if
end do

! DEALLOCATE ARRAYS
deallocate(a_p1)
deallocate(a_e1)
deallocate(a_w1)
deallocate(a_n1)
deallocate(a_s1)
deallocate(a_p2)
deallocate(a_e2)
deallocate(a_w2)
deallocate(a_n2)
deallocate(a_s2)
deallocate(F_e)
deallocate(F_w)
deallocate(F_n)
deallocate(F_s)
deallocate(phi)

contains	
	function upwind(beta)
		! INPUTS
		real(wp)					::	beta
		
		! OUTPUTS
		real(wp)					::	upwind
		
		! INTERNAL
		real(wp)					::	A,B,C,D,E,F,G,H,L,M

		! CALCULATIONS
		A = 0.0_wp
		B = a_w1(i,j)*phi(i-1,j)
		C = a_s1(i,j)*phi(i,j-1)
		D = 0.0_wp
		E = a_p1(i,j)*phi(i,j)
		F = a_p2(i,j)*phi(i,j)
		G = a_e2(i,j)*phi(i+1,j)
		H = a_w2(i,j)*phi(i-1,j)
		L = a_n2(i,j)*phi(i,j+1)
		M = a_s2(i,j)*phi(i,j-1)
		
		upwind = (1.0_wp/a_p1(i,j))*(A + B + C + D) &
			- (beta/a_p1(i,j))*(F - G - H - L - M - E + A + B + C + D)
	end function upwind

end program main