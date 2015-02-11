program main
use config
use file_io
implicit none

! DATA DICTIONARY
! -- Inputs
integer										::	n
real(wp)									::	b

! -- Parameters						
real(wp),parameter							::	density = 1.0_wp
real(wp),parameter							::	u = 2.0_wp, v=2.0_wp
real(wp),parameter							::	xsize = 1.0_wp, ysize = 1.0_wp
integer,parameter							::	max_iter = 1
real(wp),parameter							::	a_error = 1.0E-5_wp
character(40),parameter 					::  output_dir = 'output'

! -- Intermediate values
real(wp)									::	dx, dy
real(wp),allocatable,dimension(:,:)			::	a_p1, a_n1, a_s1, a_e1, a_w1
real(wp),allocatable,dimension(:,:)			::	a_p2, a_n2, a_s2, a_e2, a_w2
real(wp),allocatable,dimension(:,:)			::	F_e, F_w, F_n, F_s
integer										::	i,j,k
integer                 					::  convergence_file,phi_file,axes_file

! -- Calculated values
real(wp),allocatable,dimension(:,:)			::	phi

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

! DIRICHLET BOUNDARIES
do i=0,n+1
	phi(i,0) = 0.0_wp
end do
do j=0,n+1
	phi(0,j) = 100.0_wp
end do
phi(0,0) = 50.0_wp

! SET CONSTANTS
a_w1 = F_w
a_s1 = F_s
a_p1 = a_w1 + a_e1 + a_s1 + a_n1

a_e2 = -F_e*0.5_wp
a_n2 = -F_n*0.5_wp
a_w2 = F_w*0.5_wp
a_s2 = F_s*0.5_wp
a_p2 = a_w2 + a_e2 + a_s2 + a_n2

! CALCULATIONS
do k=1,max_iter
	do i=1,n
		do j=1,n
			if (i == n) then
				phi(i+1,j) = phi(i,j)
			end if
			if (j == n) then
				phi(i,j+1) = phi(i,j)
			end if
			phi(i,j) = (1.0_wp/a_p1(i,j)) * (a_w1(i,j)*phi(i-1,j) + a_s1(i,j)*phi(i,j-1))
		end do
	end do
end do
phi(n,n) = 0.5_wp*(phi(n-1,n)+phi(n,n-1))

do j=n,1,-1
	write(phi_file,*)phi(1:n,j)
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
	function upwind_order(order)
		! INPUTS
		integer						::	order
		
		! OUTPUTS
		real(wp)					::	upwind_order
		
		! INTERNAL VARIABLES
		
		! CALCULATIONS
		if (order == 1) then
		
		else if (order == 1) then
		
		end if
		
	end function upwind_order
	
	function upwind(beta)
		! INPUTS
		real(wp)					::	beta
		
		! OUTPUTS
		real(wp)					::	upwind
		
		! INTERNAL VARIABLES
		
		! CALCULATIONS
		upwind = upwind_order(1) !-beta*(upwind(2) - upwind(1))
	end function upwind














end program main