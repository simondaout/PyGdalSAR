program transfo_acp
    implicit none
    integer :: nargs, klen
    integer :: nx,ny,i,j
    real :: a,b,c,d
    integer :: jstart1,jstart2,jend1,jend2
    real :: val,t
    real,dimension(:),allocatable :: acp,acp_corr,zero
    character*255 acpfile,outfile

    write(*,*) 'acp file: '
    read(*,*) acpfile
    klen=index(acpfile," ")-1
    print*, acpfile
    read(*,*) nx,ny
    print*,'WIDTH=',nx,'FILE_LENGTH=',ny
    write(*,*) 'transformation parameters jstart1,jstart2,jend1,jend2: '
    read(*,*) jstart1
    read(*,*) jstart2
    read(*,*) jend1
    read(*,*) jend2
    print*, jstart1,jstart2,jend1,jend2
    write(*,*) 'fit parameters: '
    read(*,*) a
    read(*,*) b
    read(*,*) c
    read(*,*) d
    print*, a,b,c,d

    allocate(acp(nx),acp_corr(nx),zero(nx))
    do i=1,nx
        zero(i)=0.
    enddo
    open(1,file=acpfile,status='old',access='direct',recl=4*nx)
    outfile = acpfile(1:klen)//'_corr'
    open(2,file=outfile,status='unknown',access='direct',recl=4*nx)
    do j=1,ny
        read(1,rec=j)(acp(i),i=1,nx)
        if(j.lt.jstart1)then
           write(2,rec=j)(zero(i),i=1,nx)
        elseif(j.lt.jend2)then
           val=a+(b+j*(c+j*d))*j
           do i=1,nx
              acp_corr(i)=acp(i)-val
           enddo
           if(j.lt.jstart2)then
              t=float(j-jstart1)/float(jstart2-jstart1)
              do i=1,nx
                acp_corr(i)=acp_corr(i)*t
              enddo
           elseif(j.gt.jend1)then
              do i=1,nx
                t=float(j-jend1)/float(jend2-jend1)
                acp_corr(i)=acp_corr(i)*(1.-t)
              enddo
          endif
          write(2,rec=j)(acp_corr(i),i=1,nx)
        else
          write(2,rec=j)(zero(i),i=1,nx)
        endif
    enddo
    close(1)
    close(2)

    stop
end program
