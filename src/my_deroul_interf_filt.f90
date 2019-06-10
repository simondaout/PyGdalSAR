! Unwrapping of interferograms following paths with the highest possible
! coherence
!
      parameter(nbrmax=200)
      integer*8 rsc
      real*4, allocatable :: unw(:,:), coh(:,:), seuilf(:,:)
      real*4, allocatable :: ampr(:), phar(:)
      complex, allocatable :: phamod(:,:), epha(:)
      complex dx, dy
      logical, allocatable :: notdone(:)
      logical fbridge(nbrmax)
      integer, allocatable :: tab_iter(:,:)
      integer ibridgex(nbrmax,2),ibridgey(nbrmax,2),ibridjump(nbrmax)
      real, allocatable :: cut(:,:)
      integer jmin,jminnew,opt_cut
      character*12 width,filel,nbuff
      character*50 cseed
      character*255 buff,namersc,nom2,name_cut,name_unw,name_int_raw
      character*255 name_int_filtroi,nom
      logical*1, allocatable :: flag(:,:)

      nargs=iargc()
      print*,'nargs',nargs
      if(nargs.lt.9)then
      print*,'usage : deroul_interf_filt name_int_filtSW name_cut'
      print*,'name_int_raw name_int_filtroi iseedx iseedy' 
      print*,'threshold_raw  threshold_filt'
      print*,' opt_cut (1=yes, 0=no)'
      print*
      print*,'cut: negative values set to zero'
      print*,'then 0 to max value scaled between 0 and 1'
      print*,'then 1 means to mask pixel, otherwise scaling factor'
      print*,'applied to interferogram coherence above threshold'
      endif
      call getarg(1,nom)
      klen=index(nom," ")-1
      call getarg(2,name_cut)
      call getarg(3,name_int_raw)
      call getarg(4,name_int_filtroi)
 
      call getarg(5,cseed)
      read(cseed,*)iseedx
      call getarg(6,cseed)
      read(cseed,*)iseedy
      call getarg(7,cseed)
      read(cseed,*)seuil_amp
      call getarg(8,cseed)
      read(cseed,*)seuilmin
      call getarg(9,cseed)
      read(cseed,*)opt_cut
 
      namersc = nom(1:klen)//'.rsc'
      print*,namersc
      
      io=0
      open(11,file=namersc,status='old')
       do i=1,50
         read(11,*)buff,nbuff
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
      close(11)

      print*,'WIDTH=',nxr,'FILE_LENGTH=',nyr
         
      pi2=acos(-1.)*2.
      print*,pi2
      !seuil_amp=1.e-2
      
      allocate(unw(nxr,nyr))
      allocate(coh(nxr,nyr))
      allocate(seuilf(nxr,nyr))
      allocate(ampr(nxr))
      allocate(phar(nxr))
      allocate(phamod(nxr,nyr))
      allocate(epha(nxr))
      allocate(notdone(nyr))
      allocate(tab_iter(nxr,nyr))
      allocate(cut(nxr,nyr))
      allocate(flag(nxr,nyr))
 
      open(2,file=name_int_raw,form='unformatted',access='direct', &
         status='old',recl=8*nxr)
      do j=1,nyr
        read(2,rec=j)(phamod(i,j),i=1,nxr)
      enddo
      close(2)
 
      do j=1,nyr
       do i=1,nxr
        coh(i,j)=1.
        phaz=sqrt(aimag(phamod(i,j))**2+real(phamod(i,j))**2)
        if(abs(phaz).lt.seuil_amp)then
             coh(i,j)=0.
        endif
       enddo
      enddo
 
      open(2,file=nom,form='unformatted',access='direct', &
         status='old',recl=8*nxr)
      do j=1,nyr
        read(2,rec=j)(phamod(i,j),i=1,nxr)
        do i=1,nxr
          if(coh(i,j).gt.0)then
            coh(i,j)=sqrt(aimag(phamod(i,j))**2+real(phamod(i,j))**2)
          endif
        enddo
      enddo
      close(2)
 
      if(opt_cut.eq.1)then
      cutmax=0
      open(2,file=name_cut,form='unformatted',access='direct', &
         status='old',recl=4*nxr)
      do j=1,nyr
        jj=2*j
        read(2,rec=jj)(cut(i,j),i=1,nxr)
        do i=1,nxr
          if(cutmax.lt.cut(i,j))cutmax=cut(i,j)
          if(cut(i,j).lt.0)cut(i,j)=0.
        enddo
      enddo
      close(2)
      print*,'valeur du fichier cut maximale =',cutmax
      do j=1,nyr
      do i=1,nxr
        if(coh(i,j).gt.0)then
         cut(i,j)=1.-cut(i,j)/cutmax
         coh(i,j)=seuilmin*0.99+(coh(i,j)-seuilmin*0.99)*cut(i,j)
        endif
      enddo
      enddo
      endif
 
      klen=index(nom," ")-5
      name_unw=nom(1:klen)//'.unw'
      open(1,file=name_unw,form='unformatted',access='direct', &
         status='unknown',recl=4*nxr)

      open(3,file='bridge.in',status='old')
       do k=1,nbrmax
         read(3,*,end=11)ibridgex(k,1),ibridgey(k,1), &
                      ibridgex(k,2),ibridgey(k,2),ibridjump(k)
         fbridge(k)=.true.
       enddo
      close(3)
11    nbr=k-1
      print*,'nb de bridge',nbr
 
      seuilmax=1.00
      do j=1,nyr
       notdone(j)=.true.
       do i=1,nxr
        unw(i,j)=99999.
        seuilf(i,j)=0.
        if(coh(i,j).ge.seuilmin)then
          phaz=atan2(aimag(phamod(i,j)),real(phamod(i,j)))
          if(coh(i,j).gt.0.9999)then
            coh(i,j)=0.
            phamod(i,j)=cmplx(0.,0.)
          elseif(abs(phaz).lt.1.e-8)then
            coh(i,j)=0.
            phamod(i,j)=cmplx(0.,0.)
          endif
        else
             coh(i,j)=0.
             phamod(i,j)=cmplx(0.,0.)
        endif
       enddo
      enddo
 
! initialisation avec la phase filtree du premier point
      i=iseedx
      j=iseedy
      unw(iseedx,iseedy)=atan2(aimag(phamod(i,j)),real(phamod(i,j)))
      seuilf(iseedx,iseedy)=seuilmax
      tab_iter(iseedx,iseedy)=0
 
      ich=10
      ibas=0
      seuil=seuilmax
      passeuil=0.018
! domaine de recherche pour le deroulement
      jminnew=max(1,iseedy-1)
      jmaxnew=min(nyr-1,iseedy+1)
      iminnew=max(1,iseedx-1)
      imaxnew=min(nxr,iseedx+1)
! iteration
      do iter=1,500000
       if(ich.eq.0.and.abs(seuil-seuilmin).lt.1.e-6)then
          goto 30
       elseif(ich.eq.0)then
          seuil=seuil-passeuil
          seuil=max(seuil,seuilmin)
          ibas=0
! nb d iteration avant d essayer de remonter le seuil
       elseif(ibas.gt.5)then
          seuil=seuil+passeuil
          seuil=min(seuil,seuilmax)
          ibas=0
           print*,iter,seuil,ich
       endif
       ich=0
       ibas=ibas+1
       jmin=max(2,jminnew)
       jmax=min(nyr-1,jmaxnew)
       imin=max(2,iminnew)
       imax=min(nxr-1,imaxnew)
 
! Les points deroules dans l iteration i ne peuvent servir a derouler
! des points dans la meme iteration
       do i=1,nxr
         do j=1,nyr
           flag(i,j)=.true.
         enddo
       enddo
 
! iteration pour augmenter la priorite des chemins dont le
!  cumul de coherence est le plus haut ????????
 
! Mise des cotés des lignes deja entierement déroulées
       do j=max(1,jmin-1),min(nyr,jmax+1)
       if(notdone(j))then
         notdone(j)=.false.
         do i=1,nxr
           if(unw(i,j).gt.99990..and.coh(i,j).ge.seuilmin) &
                 notdone(j)=.true.
         enddo
       endif 
       enddo
 
       do j=jmin,jmax
       if(notdone(j).or.notdone(j-1).or.notdone(j+1))then
       do i=imin,imax

        if(unw(i,j).lt.99990.and.flag(i,j))then
         if(coh(i-1,j).gt.seuil.and.unw(i-1,j).gt.99990.)then
             dx=phamod(i,j)*conjg(phamod(i-1,j))
             phadx=atan2(aimag(dx),real(dx))
             unw(i-1,j)=unw(i,j)-phadx
             tab_iter(i-1,j)=iter
             seuilf(i-1,j)=1.-coh(i-1,j)+seuilf(i,j)
             ich=ich+1
             iminnew=min(i-1,iminnew)
             flag(i-1,j)=.false.
         endif
         if(coh(i,j-1).gt.seuil.and.unw(i,j-1).gt.99990.)then
             dy=phamod(i,j)*conjg(phamod(i,j-1))
             phady=atan2(aimag(dy),real(dy))
             unw(i,j-1)=unw(i,j)-phady
             tab_iter(i,j-1)=iter
             seuilf(i,j-1)=1.-coh(i,j-1)+seuilf(i,j)
             ich=ich+1
             jminnew=min(j-1,jminnew)
             flag(i,j-1)=.false.
         endif
         if(coh(i+1,j).gt.seuil.and.unw(i+1,j).gt.99990.)then
             dx=phamod(i,j)*conjg(phamod(i+1,j))
             phadx=atan2(aimag(dx),real(dx))
             unw(i+1,j)=unw(i,j)-phadx
             tab_iter(i+1,j)=iter
             seuilf(i+1,j)=1.-coh(i+1,j)+seuilf(i,j)
             imaxnew=max(i+1,imaxnew)
             flag(i+1,j)=.false.
             ich=ich+1
         endif
         if(coh(i,j+1).gt.seuil.and.unw(i,j+1).gt.99990.)then
             dy=phamod(i,j)*conjg(phamod(i,j+1))
             phady=atan2(aimag(dy),real(dy))
             unw(i,j+1)=unw(i,j)-phady
             tab_iter(i,j+1)=iter
             seuilf(i,j+1)=1.-coh(i,j+1)+seuilf(i,j)
             ich=ich+1
             jmaxnew=max(j+1,jmaxnew)
             flag(i,j+1)=.false.
         endif
 
        endif
       enddo
       endif
       enddo 

! N utilise les bridges que si ce qui est a derouler > seuil est fini
! ou essaie les bridges toutes les 6 iterations
       if(ich.eq.0)then
        do k=1,nbr
         if(fbridge(k))then
           i1=ibridgex(k,1)
           j1=ibridgey(k,1)
           i2=ibridgex(k,2)
           j2=ibridgey(k,2)
           if(unw(i1,j1).lt.99990. &
             .and.coh(i2,j2).gt.seuil.and.unw(i2,j2).gt.99990.)then
             dy=phamod(i1,j1)*conjg(phamod(i2,j2))
             phady=atan2(aimag(dy),real(dy))
             unw(i2,j2)=unw(i1,j1)-phady+ibridjump(k)*pi2
             tab_iter(i2,j2)=iter
             seuilf(i2,j2)=1.-coh(i2,j2)+seuilf(i1,j1)
             iminnew=min(i2,imin)
             jminnew=min(j2,jmin)
             imaxnew=max(i2,imax)
             jmaxnew=max(j2,jmax)
             fbridge(k)=.false.
             print*,'utilisation du bridge ',k
             ich=ich+1
           endif
           if(unw(i2,j2).lt.99990. &
             .and.coh(i1,j1).gt.seuil.and.unw(i1,j1).gt.99990.)then
             dy=phamod(i2,j2)*conjg(phamod(i1,j1))
             phady=atan2(aimag(dy),real(dy))
             unw(i1,j1)=unw(i2,j2)-phady-ibridjump(k)*pi2
             tab_iter(i1,j1)=iter
             seuilf(i1,j1)=1.-coh(i1,j1)+seuilf(i2,j2)
             iminnew=min(i1,imin)
             jminnew=min(j1,jmin)
             imaxnew=max(i1,imax)
             jmaxnew=max(j1,jmax)
             print*,'utilisation du bridge ',k
             fbridge(k)=.false.
             ich=ich+1
           endif
         endif
        enddo
       endif

      enddo 
 
30    continue
 
      do j=1,nyr
        do i=1,nxr
         if(unw(i,j).gt.99990.)coh(i,j)=0.
         if(unw(i,j).gt.99990.)unw(i,j)=0.
        enddo
      irec=2*j-1
      write(1,rec=irec)(coh(i,j),i=1,nxr)
      irec=2*j
      write(1,rec=irec)(unw(i,j),i=1,nxr)
      enddo
      close(1)
 
      open(9,file='seuil_min.r4',form='unformatted',access='direct', &
         status='unknown',recl=4*nxr)
      do j=1,nyr
      write(9,rec=j)(seuilf(i,j),i=1,nxr)
      enddo
      close(9)
 
      open(9,file='deroul_iter.i4',form='unformatted',access='direct', &
         status='unknown',recl=4*nxr)
      do j=1,nyr
      write(9,rec=j)(tab_iter(i,j),i=1,nxr)
      enddo
      close(9)
 
      open(2,file=name_int_filtroi,form='unformatted',access='direct', &
         status='old',recl=8*nxr)
      klen=index(name_int_filtroi," ")-5
      nom2=name_int_filtroi(1:klen)//'.unw'
      open(9,file=nom2,form='unformatted',access='direct', &
         status='unknown',recl=4*nxr)
      do j=1,nyr
       read(2,rec=j)(phamod(i,j),i=1,nxr)
        do i=1,nxr
         if(coh(i,j).gt.0.05)then
         ampr(i)=sqrt(real(phamod(i,j))**2+aimag(phamod(i,j))**2) 
         epha(i)=cmplx(cos(unw(i,j)),sin(unw(i,j)))
         phamod(i,j)=phamod(i,j)*conjg(epha(i))
         phar(i)=unw(i,j)+atan2(aimag(phamod(i,j)),real(phamod(i,j)))
         else
         ampr(i)=0.
         phar(i)=0.
         endif
        enddo
      irec=2*j-1
      write(9,rec=irec)(ampr(i),i=1,nxr)
      irec=2*j
      write(9,rec=irec)(phar(i),i=1,nxr)
      enddo
      close(2)
      close(9)

      open(2,file=name_int_raw,form='unformatted',access='direct', &
         status='old',recl=8*nxr)
      klen=index(name_int_raw," ")-5
      nom2=name_int_raw(1:klen)//'.unw'
      open(9,file=nom2,form='unformatted',access='direct', &
         status='unknown',recl=4*nxr)
      do j=1,nyr
       read(2,rec=j)(phamod(i,j),i=1,nxr)
        do i=1,nxr
         if(coh(i,j).gt.0.1)then
         ampr(i)=sqrt(real(phamod(i,j))**2+aimag(phamod(i,j))**2)
         epha(i)=cmplx(cos(unw(i,j)),sin(unw(i,j)))
         phamod(i,j)=phamod(i,j)*conjg(epha(i))
         phar(i)=unw(i,j)+atan2(aimag(phamod(i,j)),real(phamod(i,j)))
         else
         ampr(i)=0.
         phar(i)=0.
         endif
        enddo
      irec=2*j-1
      write(9,rec=irec)(ampr(i),i=1,nxr)
      irec=2*j
      write(9,rec=irec)(phar(i),i=1,nxr)
      enddo
      close(2)
      close(9)
 
      open(2,file=nom,form='unformatted',access='direct', &
         status='old',recl=8*nxr)
      klen=index(nom," ")-5
      nom2=nom(1:klen)//'_bridge.int'
      open(9,file=nom2,form='unformatted',access='direct', &
         status='unknown',recl=8*nxr)
      nom2=nom(1:klen)//'_bridge.unw'
      open(1,file=nom2,form='unformatted',access='direct', &
         status='unknown',recl=4*nxr)
      do j=1,nyr-1
       read(2,rec=j)(phamod(i,j),i=1,nxr)
      enddo
        do k=1,nbr
           i1=ibridgex(k,1)
           j1=ibridgey(k,1)
           i2=ibridgex(k,2)
           j2=ibridgey(k,2)
           print*,i1,i2,j1,j2
        if(abs(i1-i2).gt.abs(j1-j2))then
          do i=min(i1,i2),max(i1,i2)
            j=nint(1.e-7+j1+(i-i1)*(j2-j1+1.e-7)/(i2-i1+1.e-7))
            phamod(i,j)=cmplx(0.,0.)
            unw(i,j)=0.
            coh(i,j)=0.
          enddo
        else
          do j=min(j1,j2),max(j1,j2)
            i=nint(1.e-7+i1+(j-j1)*(i2-i1+1.e-7)/(j2-j1+1.e-7))
            phamod(i,j)=cmplx(0.,0.)
            unw(i,j)=0.
            coh(i,j)=0.
          enddo
        endif
        enddo
      do j=1,nyr
      write(9,rec=j)(phamod(i,j),i=1,nxr)
      irec=2*j-1
      write(1,rec=irec)(coh(i,j),i=1,nxr)
      irec=2*j
      write(1,rec=irec)(unw(i,j),i=1,nxr)
      enddo
      close(1)
      close(2)
      close(9)
 
      stop
      end
