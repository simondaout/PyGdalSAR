program create_bridge
  implicit none

integer :: i,k,l,n
integer, dimension(:,:), allocatable :: B
integer, dimension(:), allocatable :: col, row, cycl

allocate(B(100,5))
allocate(col(n),row(n),cycl(n))

cycl(:)=0

open(unit=10,file='bridge.tmp',status='old')
Do k = 1,n
  read(10,*) col(k), row(k), cycl(k)
End do
close(10)

k=0
Do i=1,n,2
  k=k+1
  !write(*,*)i 
  B(k,1)=col(i)
  B(k,2)=row(i)
  B(k,3)=col(i+1)
  B(k,4)=row(i+1)
  B(k,5)= cycl(i+1)
Enddo

l=0
open(unit=11,file='bridge.in',status='new')
Do l= 1,k
  write(11,*) B(k,:)
Enddo

close (11)
deallocate(B,col,row,cycl)
  
end program

