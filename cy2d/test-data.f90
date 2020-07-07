program testdata
  implicit none
  real(8), allocatable, dimension(:,:) :: x, y, z, p, u, v
  integer :: i, j, imax, jmax
  integer, parameter :: seed = 86456
  character(len=42) :: fmt
  character(len=2) :: char1, char2, char3, char4, char5, char6
character(len=
  real(8):: d1, d2, d3
  
  
  imax = 158
  jmax = 96
  
  allocate(x(imax,jmax), y(imax,jmax), z(imax,jmax), p(imax,jmax), u(imax,jmax), &
       & v(imax,jmax))


  !open(2, file="simple-input.csv")
  !rewind(2)
  !read(2,*) char1, char2, char3, char4, char5, char6
  !do j = 1, jmax
  !   do i = 1, imax
 !       read(2,*) x(i,j), y(i,j), u(i,j), v(i,j), p(i,j)
 !    end do
 ! end do
 ! close(2)

  
  
  
  do j = 1, jmax
     do i = 1, imax
        z(i,j) = 0.0
     end do
  end do

  fmt = "(F8.4,A,F8.4,A,F8.4,A,F8.4,A,F8.4,A,F8.4)"
  open(1, file="simple.csv")
  write(1,*)"x,y,z,u,v,p"
  do j = 1, jmax
     do i = 1, imax
        write(1,fmt) x(i,j),",",y(i,j),",",z(i,j),",",u(i,j),",",v(i,j),",",p(i,j)
     end do
  end do
  close(1)
  
  
end program testdata
