program main
implicit none

! DATA DICTIONARY
real::mesh_size, x_dim=1.0, y_dim=1.0
integer::n, m
real,allocatable::a_p(:,:), a_n(:,:), a_s(:,:), a_e(:,:), a_w(:,:), phi(:,:), gamma(:,:)
integer::i,j,k
real::delta,relaxation=0.8

! GET THE MESH SIZE
write(*,*)'Enter the mesh size:'
read(*,*)mesh_size
if (mesh_size > 1.0) then
	write(*,*)'Mesh size must be less than 1!'
	write(*,*)'Aborting...'
	stop
endif
n = int(x_dim/mesh_size)-1
m = int(y_dim/mesh_size)-1

! ALLOCATE THE ARRAYS
allocate(a_p(1:n,1:m))
allocate(a_n(1:n,1:m))
allocate(a_s(1:n,1:m))
allocate(a_e(1:n,1:m))
allocate(a_w(1:n,1:m))
allocate(gamma(1:n,1:m))

allocate(phi(0:n+1,0:m+1))

! INITIALIZE THE ARRAY VALUES
gamma = 1.0
phi = 0.0
a_n = gamma
a_s = gamma
a_w = gamma
a_e = gamma

! -- Boundaries
do i=1,n
	phi(i,m+1) = (real(i)-.5)*mesh_size
	a_n(i,m) = 2*a_n(i,m)
	a_s(i,1) = 2*a_s(i,1)
end do

do j=1,m
	phi(n+1,j) = (real(j)-.5)*mesh_size
	a_e(n,j) = 2*a_e(n,j)
	a_w(1,j) = 2*a_w(1,j)
end do

phi(n+1,m+1) = 1.0

! FINISH OFF COEFFICIENTS
a_p = a_w + a_e + a_n + a_s

! ITERATE
do k=1,100
	do i=1,n-1
		do j=1,m-1
			delta = (relaxation/a_p(i,j))*(a_w(i,j)*phi(i-1,j) + a_e(i,j)*phi(i+1,j) + a_n(i,j)*phi(i,j+1) + a_s(i,j)*phi(i,j-1) &
				+ a_p(i,j)*phi(i,j))
			write(*,*)delta

			phi(i,j) = phi(i,j) + delta
		end do
	end do
end do

! WRITE SOLUTION
do j=m,0,-1
	write(*,*)phi(:,j)
end do

! DEALLOCATE THE ARRAYS
deallocate(a_p)
deallocate(a_n)
deallocate(a_s)
deallocate(a_e)
deallocate(a_w)
deallocate(phi)
deallocate(gamma)


end program
