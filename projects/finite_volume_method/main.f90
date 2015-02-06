program main
use config
use file_io
implicit none

! DATA DICTIONARY
real                    ::  mesh_size, x_dim=1.0_wp, y_dim=1.0_wp
integer                 ::  n, m
real,allocatable        ::  a_p(:,:), a_n(:,:), a_s(:,:), a_e(:,:), a_w(:,:), phi(:,:), gamma(:,:)
real,allocatable        ::  x(:),y(:)
integer                 ::  i,j,k,h
real                    ::  delta,relaxation=0.8_wp
real                    ::  rss,rss_old,e,error=1.0E-5_wp

character(40),parameter ::  output_dir = 'output'
integer                 ::  convergence_file,phi_file,axes_file

! SETUP OUTPUT FILES
call create_dir(output_dir)
convergence_file = file_open(trim(output_dir) // '/conv.out')
phi_file         = file_open(trim(output_dir) // '/phi.out')
axes_file        = file_open(trim(output_dir) // '/axes.out')

! GET THE MESH SIZE
write(*,*)'Enter the mesh size:'
read(*,*)mesh_size
if (mesh_size > x_dim .OR. mesh_size > y_dim) then
  write(*,*)'Mesh size must be less than the dimensions!'
  write(*,*)'Aborting...'
  stop
endif
n = int(x_dim/mesh_size)
m = int(y_dim/mesh_size)

! ALLOCATE THE ARRAYS
allocate(a_p(1:n,1:m))
allocate(a_n(1:n,1:m))
allocate(a_s(1:n,1:m))
allocate(a_e(1:n,1:m))
allocate(a_w(1:n,1:m))
allocate(gamma(1:n,1:m))

allocate(phi(0:n+1,0:m+1))

allocate(x(0:n+1))
allocate(y(0:m+1))

! INITIALIZE THE ARRAY VALUES
gamma = 1.0_wp
phi = 0.0_wp
x = 0.0_wp
y = 0.0_wp
a_n = gamma
a_s = gamma
a_w = gamma
a_e = gamma

! -- Boundaries and other initializations
do i=1,n
  phi(i,m+1) = (real(i)-.5_wp)*mesh_size
  x(i) = (real(i)-.5_wp)*mesh_size
  a_n(i,m) = 2.0_wp*a_n(i,m)
  a_s(i,1) = 2.0_wp*a_s(i,1)
end do

do j=1,m
  phi(n+1,j) = (real(j)-.5_wp)*mesh_size
  y(j) = (real(j)-.5_wp)*mesh_size
  a_e(n,j) = 2.0_wp*a_e(n,j)
  a_w(1,j) = 2.0_wp*a_w(1,j)
end do

phi(n+1,m+1) = 1.0_wp
x(n+1) = 1.0_wp
y(m+1) = 1.0_wp

write(axes_file,"(A,' ',A)")'x','y'
do k=0,max(m+1,n+1)
  write(axes_file,*)x(k),y(k)
end do

! FINISH OFF COEFFICIENTS
a_p = a_w + a_e + a_n + a_s

! ITERATE
! Goes through a maximum of 10000 iterations, but stops once
! RSS change is less than the specified value for three 
! consecutive iterations.
h = 0
write(convergence_file,"(3(A,' '),A)")'ITER','RSS','dRSS','H'
do k=1,10000
  rss = 0.0_wp
  do i=1,n
    do j=1,m
      delta = (relaxation/a_p(i,j))*(a_w(i,j)*phi(i-1,j) + a_e(i,j)*phi(i+1,j) + a_n(i,j)*phi(i,j+1) + a_s(i,j)*phi(i,j-1) &
        - a_p(i,j)*phi(i,j))
      phi(i,j) = phi(i,j) + delta
      e = abs(phi(i,j) - x(i)*y(j))
      rss = rss + e**2
    end do
  end do
  rss = sqrt(rss)
  if (abs(rss-rss_old) < error) then
    h = h + 1
  else
    h = 0
  endif
  !write(convergence_file,"(I4,' ',2(F10.7,' '),I1)")k,rss,abs(rss-rss_old),h
  write(convergence_file,*)k,rss,abs(rss-rss_old),h
  if (h >= 3) then
    exit
  endif
  rss_old = rss
end do

! WRITE SOLUTION
do j=0,m+1
  write(phi_file,*)phi(:,j)
end do

! DEALLOCATE THE ARRAYS
deallocate(a_p)
deallocate(a_n)
deallocate(a_s)
deallocate(a_e)
deallocate(a_w)
deallocate(phi)
deallocate(gamma)

! CLOSE OUTPUT FILES
call file_close(convergence_file)
call file_close(phi_file)

end program