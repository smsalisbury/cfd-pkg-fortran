program main
implicit none

! DATA DICTIONARY
real::mesh_size, x_dim=1.0, y_dim=1.0
integer::n, m
real,allocatable::a_n(:,:), a_s(:,:), a_e(:,:), a_w(:,:), phi(:,:), gamma(:,:)
integer::i,j,k
real::a_p,delta,relaxation=0.8

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
allocate(a_n(0:n,0:m))
allocate(a_s(0:n,0:m))
allocate(a_e(0:n,0:m))
allocate(a_w(0:n,0:m))
allocate(phi(0:n,0:m))
allocate(gamma(0:n,0:m))

! INITIALIZE THE ARRAY VALUES
gamma = 1.0
phi = 0.0
a_n = gamma
a_s = gamma
a_w = gamma
a_e = gamma

! -- Boundaries
do i=0,n
	phi(i,m) = real(i)*mesh_size ! FIX THIS
	a_n(i,m) = 2*a_n(i,m)
	a_s(i,m) = 2*a_s(i,m)
end do

do j=0,m
	phi(n,j) = real(j)*mesh_size ! FIX THIS
	a_e(n,j) = 2*a_e(n,j)
	a_w(n,j) = 2*a_w(n,j)
end do

! ITERATE
do k=1,100
	do i=1,n-1
		do j=1,m-1
			a_p = a_w(i,j) + a_e(i,j) + a_n(i,j) + a_s(i,j)
			delta = (relaxation/a_p)*(a_w(i,j)*phi(i-1,j) + a_e(i,j)*phi(i+1,j) + a_n(i,j)*phi(i,j+1) + a_s(i,j)*phi(i,j-1) &
					 + a_p*phi(i,j))

			phi(i,j) = phi(i,j) + delta
		end do
	end do
end do

! WRITE SOLUTION
do j=m,0,-1
	write(*,*)phi(:,j)
end do

! DEALLOCATE THE ARRAYS
deallocate(a_n)
deallocate(a_s)
deallocate(a_e)
deallocate(a_w)
deallocate(phi)
deallocate(gamma)


end program
