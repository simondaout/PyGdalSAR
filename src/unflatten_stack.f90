program unflatten_stack
! add stack pattern to the interferogram after unwrapping 
    implicit none

    integer :: nargs, io, j, i, jrec, klen
    integer :: nxr,nyr
    !parameter(nx1=12000,nsumymax=1000)
    character*255 buff, infile, outfile, acpfile, namersc, coeffile
    character*12 width,filel,nbuff
    
    real,dimension(:),allocatable :: amp,pha,pham,acp
    !complex unw(nx1,nsumymax)
    !real acp(nx1,nsumymax)
    
    complex :: epha0   
    real :: coef

    nargs=iargc()
    if (nargs.lt.3) then
        print*, 'Syntax: unflatten_stack infile outfile stack_file coeffile'
        call exit
    endif

    call getarg(1,infile)
    klen=index(infile," ")-1
    call getarg(2,outfile)
    call getarg(3,acpfile)
    call getarg(4,coeffile)

   namersc = infile(1:klen)//'.rsc'
   print*,namersc

   io=0
   open(2,file=namersc,status='old')
   do i=1,50
         read(2,*)buff,nbuff
         read(buff,'(a5)')width
         if(width.eq.'WIDTH')then
           io=io+1
           read(buff,*)width
           read(nbuff,*)nxr
           print*,width,nxr
         endif
         read(buff,'(a11)')filel
         if(filel.eq.'FILE_LENGTH')then
           io=io+1
           read(buff,*)filel
           read(nbuff,*)nyr
           print*,filel,nyr
         endif
         if(io.eq.2) goto 10
   enddo
10     continue
   close(2)
!
      print*,'WIDTH=',nxr,'FILE_LENGTH=',nyr
!
    allocate(amp(nxr),pha(nxr),pham(nxr),acp(nxr))

    open(1,file=infile,form='unformatted',access='direct',status='old',recl=4*nxr)
    open(2,file=outfile,form='unformatted',access='direct',status='unknown',recl=4*nxr)
    open(3,file=acpfile,form='unformatted',access='direct',status='old',recl=4*nxr)
    
    ! read coef 
    open(4,file=coeffile,status='old')
        read(4,*) coef
        print*, 'coef: ', coef
    close(4)
    !coef=-coef

    do j=1,nyr
        read(3,rec=j)(acp(i),i=1,nxr)
        jrec=j*2-1
        read(1,rec=jrec)(amp(i),i=1,nxr)
        !print*, 'amp: ', amp(1:6)
        jrec=j*2
        read(1,rec=jrec)(pha(i),i=1,nxr)
        !print*, 'pha: ', pha(1:6)
      
        do i=1,nxr
            if(abs(amp(i)).gt.0.001)then
                pham(i)=pha(i)-coef*acp(i)
            else
                pham(i)=pha(i)
            endif
        end do
        !print*, pham(1:5)
        !print*, pha(1:5)
        !print*
        jrec=j*2-1
        write(2,rec=jrec)(amp(i),i=1,nxr)
        jrec=j*2
        write(2,rec=jrec)(pham(i),i=1,nxr)
    
    end do
    close(1)
    close(2)


end program 
