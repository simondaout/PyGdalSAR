program ACP_unw
  implicit none
  

  integer :: nargs                              ! Argument Counts       
  integer, parameter :: nmax_int=300            ! Nb max int for reading  
  integer, parameter :: num_eigen_threshold=4
  
  character(len=255) :: argu                    ! Argument reading variable
  character(len=255) :: liste_int               ! Interferogram list
  character(len=255) :: prefix,suffix,nomdir,namersc,nomint            
  character(len=8) :: cinter(nmax_int,2) 
  integer :: inter(nmax_int,2)                  
  integer :: nlign,ncol                         ! Size unw
  
  character(len=255) :: buff,width,filel        ! Variables for rsc file
  character(len=12) :: nbuff
  integer :: io 
  
  integer :: klen_int,klen_pre,klen_suf,klen_dir,klen_rsc,l,nmax,knmax,istart,iend
  integer :: kmax                               ! Nb max int after reading
  integer :: k,kk                               ! Loop on int
  integer :: kref                               ! image ref
  integer :: i,j,ii                             ! Loop on ncol nlign
  
  integer :: jstart,jend                        ! ACP window size 
  integer :: jstartold,jendold,jmax,jmin,jrecl  
  integer :: nlignint(nmax_int)               

  real,dimension(:,:),allocatable :: deplac 
 
  double precision,dimension(:,:),allocatable :: var_cov  ! Covariance matrix 
  real,dimension(:,:),allocatable :: nsum     ! Sum matrix 
  real,dimension(:,:),allocatable :: proj     ! proj ACP 
  real,dimension(:),allocatable :: amoy,asumc       ! Average 
  integer,dimension(:),allocatable :: nsumi
  character(len=1) :: JOBZ, UPLO                          ! SVD variables
  integer :: INFO, LDA, LIWORK, LWORK
  integer, dimension(:), allocatable :: IWORK
  double precision, dimension(:), allocatable :: W,WORK
  real :: asum
  
  ! read arguments
  nargs=iargc()
  if(nargs.lt.5)then
    print*,'ACP_unw liste_interfero prefix_unw suffix_unw jstart jend (istart iend) '
    print*,'produces ACP decomposition of unwrapped int between jstart and jend'
    stop
  endif
  call getarg(1,liste_int)
  klen_int=index(liste_int," ")-1
  call getarg(2,prefix)
  klen_pre=index(prefix," ")-1
  call getarg(3,suffix)
  klen_suf=index(suffix," ")-1
  call getarg(4,argu)
  read(argu,fmt=*)jstart
  call getarg(5,argu)
  read(argu,fmt=*)jend
  istart=1
  iend=100000
  if(nargs.gt.5)then
  call getarg(6,argu)
  read(argu,fmt=*)istart
  call getarg(7,argu)
  read(argu,fmt=*)iend
  endif

  ! read list of interferograms
  print*,'interferogram list: ',liste_int(1:klen_int)
  open(1,file=liste_int(1:klen_int),status='old')
  do k=1,nmax_int
    read(1,*,end=10)cinter(k,1),cinter(k,2)
    read(cinter(k,1),*)inter(k,1)
    read(cinter(k,2),*)inter(k,2)
  end do
10  kmax=k-1
  print*,'number of interferogram: ', kmax
  !write(*,*) cinter(kmax,1), cinter(kmax,2)
  !write(*,*) inter(kmax,1),inter(kmax,2)
  close(1)
  
  ! read rsc file for each interferograms
  do k=1,kmax
    nomdir='int_'//cinter(k,1)//'_'//cinter(k,2)
    klen_dir=index(nomdir," ")-1
!    print*,'directory du rsc :',nomdir(1:klen_dir)
    namersc = nomdir(1:klen_dir)//'/'//prefix(1:klen_pre)//cinter(k,1)//'-'//cinter(k,2)//suffix(1:klen_suf)//'.unw.rsc'
    klen_rsc=index(namersc," ")-1
!    print*,'fichier rsc lu :',namersc(1:klen_rsc)
    io=0
    open(1,file=namersc,status='old')
    do kk=1,50
      read(1,*)buff,nbuff
      read(buff,'(a5)')width
      if(width.eq.'WIDTH')then
        io=io+1
        read(nbuff,*)ncol
!        print*,width,ncol
      endif
      read(buff,'(a11)')filel
      if(filel.eq.'FILE_LENGTH')then
        io=io+1
        read(nbuff,*)nlign
!        print*,filel,nlign
      endif
      if(io.eq.2) goto 12
    enddo
 12      continue
    close(1)
    jmax=nlign
    kref=k
    nlignint(k)=nlign
  enddo

  do k=1,kmax
    if(jmax.lt.nlignint(k)) then
      jmax=nlignint(k)
      kref=k
    end if
  end do 

  iend=min(iend,ncol)
 
  ! read unw phase for each interferogram
  allocate(deplac(kmax,ncol),nsum(kmax,kmax),nsumi(kmax))
  allocate(amoy(kmax))
  allocate(var_cov(kmax,kmax))
  amoy(:)=0.
  var_cov(:,:)=0.
  nsum(:,:)=0
  nsumi(:)=0

  do k=1,kmax
        nomdir='int_'//cinter(k,1)//'_'//cinter(k,2)
        klen_dir=index(nomdir," ")-1
        nomint=nomdir(1:klen_dir)//'/'//prefix(1:klen_pre)//cinter(k,1)//'-'//cinter(k,2)//suffix(1:klen_suf)//'.unw'
!        print*,'j,jmax: ', k,nlignint(k)
!        print*,'read: ',nomint
        open(10+k,file=nomint,form='unformatted',access='direct',status='old',recl=4*ncol)
  enddo

  do j=jstart,jend,2
    do k=1,kmax
      deplac(k,:)=0.
      if (j.lt.nlignint(k)) then
        jrecl=2*j
        read(10+k,rec=jrecl)(deplac(k,ii),ii=1,ncol)
      endif
    enddo

    ! calcul of the covariance matrix 
    do i=istart,iend,2
      do k=1,kmax
        do kk=k,kmax
          nsum(k,kk)=nsum(k,kk)+1
        enddo
        if (j.lt.nlignint(k)) then
          amoy(k)=amoy(k)+deplac(k,i)
          if (abs(deplac(k,i)).gt.1.e-5)nsumi(k)=nsumi(k)+1
          do kk=k,kmax
!            if (abs(deplac(k,i)).gt.1.e-5.and.abs(deplac(kk,i)).gt.1.e-5) then
            var_cov(k,kk)=var_cov(k,kk)+deplac(k,i)*deplac(kk,i)
!            nsum(k,kk)=nsum(k,kk)+1
!            end if
          end do
        end if
      end do
    end do
  end do

  nmax=0
  do k=1,kmax
    if(nsumi(k).gt.nmax)then
     nmax=nsumi(k)
     knmax=k
    endif
  enddo
  print*,'number of points',nsumi(:)
  do k=1,kmax
    if(nsumi(k).lt.(nmax*2)/3)then
      print*,'variance image',k,'arbitrarily decreased by 10'
      amoy(k)=amoy(k)/10.
      do kk=k,kmax
        var_cov(k,kk)=var_cov(k,kk)/10.
      enddo
      if(k.gt.1)then
      do kk=1,k-1
        var_cov(kk,k)=var_cov(kk,k)/10.
      enddo
      endif
    endif
  enddo

  do k=1,kmax
    close(10+k)
  enddo

  deallocate(deplac)

  do k=1,kmax
     amoy(k)=amoy(k)/nsum(k,k)
  enddo
  do k=1,kmax
    do kk=k,kmax
      var_cov(k,kk)=var_cov(k,kk)/nsum(k,kk)-amoy(k)*amoy(kk)
    enddo
  enddo
  do k=2,kmax
    do kk=1,k-1
      var_cov(k,kk)=var_cov(kk,k)
    enddo
  enddo

! dsyevd.f
! SVD decomposition 
! methode divide and conquer for SVD

  JOBZ='V'
  UPLO='U'
  LDA=kmax
  INFO=0

        LWORK = -1
        LIWORK = 1
        allocate(IWORK(LIWORK),WORK(1),W(kmax))
        call DSYEVD( JOBZ, UPLO, kmax, var_cov, LDA, W, WORK, LWORK, IWORK,LIWORK, INFO )
        LWORK=WORK(1)
        deallocate(WORK)
        allocate(WORK(LWORK))
        print*, LWORK

        LIWORK = -1
        call DSYEVD( JOBZ, UPLO, kmax, var_cov, LDA, W, WORK, LWORK, IWORK,LIWORK, INFO )
        LIWORK=IWORK(1)
        deallocate(IWORK)
        allocate(IWORK(LIWORK))
        print*, LIWORK

        call DSYEVD( JOBZ, UPLO, kmax, var_cov, LDA, W, WORK, LWORK, IWORK,LIWORK, INFO )

  deallocate(IWORK,WORK)
  write(*,*) 'eigenvalues: ', W(:)

      do k=1,kmax
        asum=0.
        do kk=1,kmax
          asum=asum+var_cov(kk,k)*var_cov(kk,k)
        enddo
        do kk=1,kmax
          var_cov(kk,k)=var_cov(kk,k)/sqrt(asum)
        enddo
      enddo

      open(1,file='acp_eigenvalues',status='unknown')
       write(1,*)W(:)
      close(1)

      open(1,file='acp_vecteurs',status='unknown')
        do k=1,kmax
        write(1,'(a8,a1,a8,5f9.4)')cinter(k,1),' ',cinter(k,2),(var_cov(k,kmax-kk),kk=0,4)
        enddo
      close(1)

  ! projection in the base of the eigenvector
  open(4,file='acp_1',status='unknown',form='unformatted',access='direct',recl=ncol*4)
  open(5,file='acp_2',status='unknown',form='unformatted',access='direct',recl=ncol*4)
  open(6,file='acp_3',status='unknown',form='unformatted',access='direct',recl=ncol*4)
  open(7,file='acp_4',status='unknown',form='unformatted',access='direct',recl=ncol*4)
  print*,'ouverture depl_cumule'
  allocate(proj(ncol,4))
!
  allocate(deplac(kmax,ncol),asumc(kmax))
  do k=1,kmax
        nomdir='int_'//cinter(k,1)//'_'//cinter(k,2)
        klen_dir=index(nomdir," ")-1
        nomint=nomdir(1:klen_dir)//'/'//prefix(1:klen_pre)//cinter(k,1)//'-'//cinter(k,2)//suffix(1:klen_suf)//'.unw'
        open(100+k,file=nomint,form='unformatted',access='direct',status='old',recl=4*ncol)
  enddo

  do j=1,jmax
    do k=1,kmax
      if (j.lt.nlignint(k)) then
        deplac(k,:)=0.
        jrecl=2*j
        read(100+k,rec=jrecl)(deplac(k,ii),ii=1,ncol)
      end if
    end do

    proj(:,:)=0.
    do i=1,ncol
      asumc(:)=0.
      do k=1,kmax
        if (j.lt.nlignint(k)) then
          do l=1,num_eigen_threshold
              if(abs(deplac(k,i)).gt.1.e-5) then
                proj(i,l)=proj(i,l)+(deplac(k,i)-amoy(k))*var_cov(k,kmax-l+1)
                asumc(l)=asumc(l)+var_cov(k,kmax-l+1)*var_cov(k,kmax-l+1)
              endif
          enddo
        endif
      enddo
      do l=1,num_eigen_threshold
        if(asumc(l).gt.0.00001)proj(i,l)=proj(i,l)/sqrt(asumc(l))
      enddo
    enddo
    jrecl=j
    write(4,rec=jrecl)(proj(i,1),i=1,ncol)
    write(5,rec=jrecl)(proj(i,2),i=1,ncol)
    write(6,rec=jrecl)(proj(i,3),i=1,ncol)
    write(7,rec=jrecl)(proj(i,4),i=1,ncol)
  enddo
  do k=1,kmax
    close(100+k)
  end do

  close(4)
  close(5)
  close(6)
  close(7)
  deallocate(deplac,amoy,var_cov)

  stop
end program
