ccherche le meilleur fit des franges de l interfero avec une fonction : d f/d acp = acycl0 
c donne les valeurs locales de df/d acp sur l ensemble de l intefero
c a nettoyer ensuite pour trouver  acycl0 
c
c Fait des calculs sur nregx*nregy regions en x,y
c
c Critere pour selectionner pente phase/topo : ampmax, sigma_topo, coh
c pente phase/topo semble OK si ampmax >0.05 ou si sigma_topo > 100. (produit ampmax*sigma_topo>6.)
c ou produit coh*sigma_topo> 17.
c
c On peut eliminer les zones ou sigma < seuil_sz et ou ampmax < seuil_cor
c
c nsumy : taille des zones en azimuth
c nregx : nb de régions en range
c Dtopo : Range of acp values for cyclmax phase cycles
c
      parameter(nx1=12000,nsumymax=1000,nregx=7,ny1=80000)
      parameter(nregymax=600)
      parameter(nmedmax=800,nvar=3,nvar2=4,nvar3=5)
      real aslop(nregx*nregymax),qual(nregx*nregymax)
      real topoz(nregx*nregymax),ypos(nregx*nregymax)
      real x(nmedmax),y(nmedmax),w(nmedmax)
      real aslop_med(nmedmax),apos_med(nmedmax)
      real qual_med(nmedmax),alt_ref
      real*8 topom,sigmaz,aslopef
      complex phamodf(nx1,nsumymax)
      complex epha0
      complex*16 sump
      real*8 amp
      real topo(nx1,nsumymax)
      character*12 width,filel,nbuff
      character*50 aaa
      character*255 buff,namersc,infile, outfile,demfile
      character*255 infilef, outfilef
      character*5 fit

      nargs=iargc()
      if (nargs.lt.6) then
       print*,
     $  'Syntax: flatten_acp infile infile_filtered acp_file 
     $      outfile outfile_filtered thresh_amp'
       print*,'infile : interferogram'
       print*,'infile_filtered : filtered interf (with sliding windows)'
       print*,'acp projection file'
       print*,'thresh_amp: threshold on amplitude of infile '
       print*,'            (if coherence, 0.2 might be good choice)'
       call exit
      endif
      call getarg(1,infile)
      klen=index(infile," ")-1
      call getarg(2,infilef)
      klenf=index(infilef," ")-1
      call getarg(3,demfile)
      call getarg(4,outfile)
      call getarg(5,outfilef)
      call getarg(6,aaa)
      read(aaa,*)thresh_amp


cParametres:
c thresh_amp : seuil sur l amplitude : 0.2 assez restrictif
c test avec une valeur plus cool: 0.1
c seuil_cor : seuil sur coherence (0.08), depend du seuil sur l amplitude
c prendre plus bas (0.05) si seuil sur amp a 0.17
c test avec valeur plus cool: 0.04
        seuil_cor=0.1
c seuil variabilité des altitudes
        seuil_sz=3.
c acp minimale
        seuil_acp=-50
c seuil sur le nb min de points valables dans chaque cellule
        seuil_nb=50
c       
      altmin=10000.
      altmax=-10000.

      namersc = infile(1:klen)//'.rsc'
      print*,namersc
c
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
c
      print*,'WIDTH=',nxr,'FILE_LENGTH=',nyr
c
      open(1,file=infilef,form='unformatted',access='direct',
     $   status='old',recl=8*nxr)
c
      open(3,file=demfile,form='unformatted',
     $   access='direct',status='old',recl=4*nxr)
c
      open(2,file='ncycle_acp',status='unknown')
      open(5,file='ncycle_acp_az',status='unknown')
      write(5,'(a)')'@    s0 symbol 1'
      write(5,'(a)')'@    s0 symbol size 0.260000'
      write(5,'(a)')'@    s0 symbol color 1'
      write(5,'(a)')'@    s0 line type 0'
c
      nsumx=((nxr-1)*2)/(nregx+1)
      nregy=2*(((nregx+1)/2)*nyr)/nxr+1
      nsumy=((nyr-1)*2)/(nregy+1)
c     nb max de cycles de phase fonction de topo (sur Dtopo)
      Dtopo=40
      cyclmax=2
c     nb de pas par cycle
      npascycl=10
      ismax=nint(cyclmax*npascycl)
      fac=acos(-1.)
      fac=2.*fac*cyclmax/(Dtopo*ismax)
      icc=0
c
c     boucle en azimuth
      kj=0
      do ky=1,nregy
        iyb=((ky-1)*nsumy)/2+1
        iyb=max(1,iyb)
        iye=iyb+nsumy
        iye=min(nyr,iye)
        do j =iyb,iye,1
          jind=j-iyb+1
          kj=kj+1
          read(1,rec=j,err=20)(phamodf(i,jind),i=1,nxr)
          read(3,rec=j)(topo(i,jind),i=1,nxr)
        enddo
c
c     boucle en range
      do kx=1,nregx
        ixb=((kx-1)*nsumx)/2+1
        ixb=max(10,ixb)
        ixe=ixb+nsumx
        ixe=min(nxr-10,ixe)
c
        topom=0.
c        if(ky.lt.10)then
c        write(buff,'(a29,i1,a1,i1)')
c     $    'PIX_PHASE_TOPO/pix_phase_topo',kx,'_',ky
c        else
c        write(buff,'(a29,i1,a1,i2)')
c     $    'PIX_PHASE_TOPO/pix_phase_topo',kx,'_',ky
c        endif
c        open(4,file=buff,status='unknown')
        ic=0
        zmin=10000.
        zmax=0.
        sigmaz=0.
        weight=0.
        topom=0.
        do j =iyb,iye,1
          jind=j-iyb+1
        do i =ixb,ixe,1
c          write(aaa,*)topo(i,jind)
c          if(aaa(3:5).eq.'NAN') topo(i,jind)=0.
          if(topo(i,jind).lt.seuil_acp) phamodf(i,jind)=cmplx(0.,0.)
          if(topo(i,jind).lt.altmin.and.topo(i,jind).ge.seuil_acp) 
     $             altmin=topo(i,jind)
          if(topo(i,jind)>altmax) altmax=topo(i,jind)
          amp=sqrt(aimag(phamodf(i,jind))**2+real(phamodf(i,jind))**2)
          if(amp.gt.thresh_amp)then
            ic=ic+1
            weight=weight+amp
            topom=topom+topo(i,jind)
            sigmaz=sigmaz+topo(i,jind)**2
            pha=atan2(aimag(phamodf(i,jind)),real(phamodf(i,jind)))
c            write(4,*)topo(i,jind),pha
            if(topo(i,jind).gt.zmax)zmax=topo(i,jind)
            if(topo(i,jind).lt.zmin)zmin=topo(i,jind)
          endif
        enddo
        enddo

        iop=0
        if(ky.eq.1.and.kx.gt.(nregx-2))iop=1
        if(ic.gt.seuil_nb.and.iop.eq.0)then
        topom=topom/ic
        sigmaz=sqrt(sigmaz/ic -topom**2)
        amax=-10.
        do islope=-ismax,ismax,1
          sump=(0.,0.)
          aslope=islope*fac
          do j =iyb,iye,1
          jind=j-iyb+1
          do i =ixb,ixe,1
           amp=sqrt(aimag(phamodf(i,jind))**2+real(phamodf(i,jind))**2)
           if(amp.gt.thresh_amp)then
             val=aslope*topo(i,jind)
             epha0=cmplx(cos(val),sin(val))
             sump=sump+phamodf(i,jind)*epha0
           endif
          enddo
          enddo
          ampsump=dreal(sump)**2+dimag(sump)**2
          ampsump=ampsump/ic
          if(ampsump.gt.amax)then
                islopef=islope
                amax=ampsump
                pha=atan2(aimag(sump),real(sump))
           endif
        enddo
        coh=amax/weight
        amax=amax/ic
        islope=islopef
c        write(4,'(a1)')'>'
c        do kz=1,100
c          z=(kz-1.)*(zmax-zmin)/100.+zmin
c          val=islopef*z*fac-pha
c          epha0=cmplx(cos(val),sin(val))
c          pha2=atan2(aimag(epha0),real(epha0))
c          write(4,*)z,-pha2
c        enddo
c
c        close(4)
c
        write(2,'(2f7.0,f8.2,i4,2f6.3,2i3,3f8.2)')
     $    (kx-0.5)*nxr/nregx+0.5,
     $    (ky-0.5)*nyr/nregy+0.5,
     $    topom,islopef,amax,coh,kx,ky,zmin,zmax,sigmaz
c
         if(coh.gt.seuil_cor.and.sigmaz.gt.seuil_sz)then 
           write(5,*)(ky-0.5)*nyr/nregy+0.5,islopef*fac
           icc=icc+1
           aslopef=islopef*fac
           aslop(icc)=aslopef
           qual(icc)=amax
           topoz(icc)=topom
           ypos(icc)=(ky-0.5)*nyr/nregy+0.5
          endif
        endif
c
      enddo
      enddo
      goto 40
20    nyr=kj-1
40    print*,'file length',nyr
      close(2)
      close(1)
      close(3)
c
      print*,'minimum acp',altmin
      print*,'maximum acp',altmax
c
c     mediane sur fenetre glissante en fonction de l azimuth sur ls cycles
c     ponderee par la coherence
      call DSORT (topoz,aslop,qual, icc, 2)
      write(5,'(a8)')'@type xy'
      lmed=10
      ip=0
      do ll=1,icc-lmed+1
         apos=0.
         wtot=0.
         do l=0,lmed-1
            x(l)=aslop(ll+l)
            w(l)=qual(ll+l)
            if(fit_az.eq.0)then
              y(l)=ypos(ll+l)
            else
c fonction de z
              y(l)=topoz(ll+l)
            endif
            apos=apos+y(l)*w(l)
            wtot=wtot+w(l)
         enddo
         call amediane(x,lmed,w,nmedmax,amed)
         apos=apos/wtot
         write(5,*)apos,amed
c         if(apos.lt.5700.or.apos.gt.7700)then
           ip=ip+1
           aslop_med(ip)=amed
           apos_med(ip)=apos/100.
           qual_med(ip)=wtot/lmed
c         endif
      enddo
c
c Mediane
      nfit=-1
      if(nfit.eq.-1)then
      call amediane(aslop_med,ip,qual_med,nmedmax,baz)
      print*,'mediane',baz
      write(5,'(a8)')'@type xy'
      do k=1,100
c fonction de y
        apos=(k-1)*nyr/100.
        val=baz
        write(5,*)apos,val
      enddo
c
      elseif(nfit.eq.0)then
      amean=0.
      ameanp=0.
      do i=1,ip
       amean=amean+aslop_med(i)*qual_med(i)
       ameanp=ameanp+qual_med(i)
      enddo
      baz=amean/ameanp
      print*,'valeur constante',baz
      write(5,'(a8)')'@type xy'
      do k=1,100
        apos=(k-1)*nyr/100.
        val=baz
        write(5,*)apos,val
      enddo
      endif
c
c     ecriture du fit en range
      buff=infile(1:klen-4)//'.acp'
      open(2,file=buff,form='formatted',status='unknown')
           write(2,'(4e14.6)')baz
      close(2)
c
c     correction de l interfero filtre
      open(1,file=infilef,form='unformatted',access='direct',
     $   status='old',recl=8*nxr)
      open(2,file=outfilef,form='unformatted',access='direct',
     $   status='unknown',recl=8*nxr)
      open(3,file=demfile,form='unformatted',
     $   access='direct',status='old',recl=4*nxr)
      do j=1,nyr-1
        read(1,rec=j,err=30)(phamodf(i,1),i=1,nxr)
        read(3,rec=j)(topo(i,1),i=1,nxr)

        aslopef=baz
        do i=1,nxr
c          write(aaa,*)topo(i,1)
c          if(aaa(3:5).eq.'NAN') topo(i,1)=0.
          epha0=cmplx(cos(aslopef*topo(i,1)),
     $              sin(aslopef*topo(i,1)))
          phamodf(i,1)=phamodf(i,1)*epha0
        enddo

        write(2,rec=j)(phamodf(i,1),i=1,nxr)
      enddo
30    nyr=j-1
      print*,'file length',nyr
      close(1)
      close(2)
c
c     correction de l interfero original
      open(1,file=infile,form='unformatted',access='direct',
     $   status='old',recl=8*nxr)
      open(2,file=outfile,form='unformatted',access='direct',
     $   status='unknown',recl=8*nxr)
      open(3,file=demfile,form='unformatted',
     $   access='direct',status='old',recl=4*nxr)
      do j=1,nyr
        read(1,rec=j)(phamodf(i,1),i=1,nxr)
        read(3,rec=j)(topo(i,1),i=1,nxr)

          aslopef=baz
          do i=1,nxr
c            write(aaa,*)topo(i,1)
c            if(aaa(3:5).eq.'NAN') topo(i,1)=0.
            epha0=cmplx(cos(aslopef*topo(i,1)),
     $              sin(aslopef*topo(i,1)))
            phamodf(i,1)=phamodf(i,1)*epha0
          enddo

        write(2,rec=j)(phamodf(i,1),i=1,nxr)
      enddo
      close(1)
      close(2)
c
      stop
      end
c
c********************************************************************************
      subroutine amediane(x,n,weight,nt,amed)
c
      parameter (nmax=10000)
      real x(nt),weight(nt),xtri(nmax),weighttri(nmax)

      loopsize = n
      do 100 i = 1, loopsize
         xtri(i) = x(i)
         weighttri(i)=weight(i)
 100  end do


      do while (loopsize .gt.1)
         do 200 i = 2, loopsize
         if(xtri(i).lt.xtri(i-1)) then
                 temp = xtri(i)
                 xtri(i) = xtri(i-1)
                 xtri(i-1) = temp
                 temp = weighttri(i)
                 weighttri(i)=weighttri(i-1)
                 weighttri(i-1) = temp
         endif
 200     enddo
         loopsize = loopsize - 1
      end do

      do i=2,n
        weighttri(i)=weighttri(i)+weighttri(i-1)
      enddo

      i=1
      do while(weighttri(i).lt.weighttri(n)/2.)
        i=i+1
      enddo

      t=(weighttri(i)-weighttri(n)/2.)/(weighttri(i)-weighttri(i-1))
      amed=t*xtri(i-1)+(1.-t)*xtri(i)

      return
      end

      SUBROUTINE DSORT (DX, DY, DY1, N, KFLAG)
!***BEGIN PROLOGUE  DSORT
!***PURPOSE  Sort an array and optionally make the same interchanges in
!            an auxiliary array.  The array may be sorted in increasing
!            or decreasing order.  A slightly modified QUICKSORT
!            algorithm is used.
!***LIBRARY   SLATEC
!***CATEGORY  N6A2B
!***TYPE      DOUBLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
!***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
!***AUTHOR  Jones, R. E., (SNLA)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   DSORT sorts array DX and optionally makes the same interchanges in
!   array DY.  The array DX may be sorted in increasing order or
!   decreasing order.  A slightly modified quicksort algorithm is used.
!
!   Description of Parameters
!      DX - array of values to be sorted   (usually abscissas)
!      DY - array to be (optionally) carried along
!      N  - number of values in array DX to be sorted
!      KFLAG - control parameter
!            =  2  means sort DX in increasing order and carry DY along.
!            =  1  means sort DX in increasing order (ignoring DY)
!            = -1  means sort DX in decreasing order (ignoring DY)
!            = -2  means sort DX in decreasing order and carry DY along.
!
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891009  Removed unreferenced statement labels.  (WRB)
!   891024  Changed category.  (WRB)
!   891024  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   901012  Declared all variables; changed X,Y to DX,DY; changed
!           code to parallel SSORT. (M. McClain)
!   920501  Reformatted the REFERENCES section.  (DWL, WRB)
!   920519  Clarified error messages.  (DWL)
!   920801  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
!***END PROLOGUE  DSORT
!     .. Scalar Arguments ..
      INTEGER :: KFLAG, N
!     .. Array Arguments ..
      real, dimension(:) :: DX(N)
      real, dimension(:) ::  DY(N),DY1(N)
!     .. Local Scalars ..
      real :: R, T, TT
      real :: TY, TTY, TY1, TTY1
      INTEGER :: I, IJ, J, K, KK, L, M, NN
!     .. Local Arrays ..
      INTEGER, dimension(:) :: IL(21), IU(21)
!***FIRST EXECUTABLE STATEMENT  DSORT
      NN = N
      KK = ABS(KFLAG)
!
!     Alter array DX to get decreasing order if needed
!
      IF (KFLAG <= -1) THEN
         DO 10 I=1,NN
            DX(I) = -DX(I)
   10    CONTINUE
      ENDIF
!
      IF (KK .EQ. 2) GO TO 100
!
!     Sort DX only
!
      M = 1
      I = 1
      J = NN
      R = 0.375D0
!
   20 IF (I .EQ. J) GO TO 60
      IF (R <= 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF
!
   30 K = I
!
!     Select a central element of the array and save it in location T
!
      IJ = I + INT((J-I)*R)
      T = DX(IJ)
!
!     If first element of array is greater than T, interchange with T
!
      IF (DX(I) > T) THEN
         DX(IJ) = DX(I)
         DX(I) = T
         T = DX(IJ)
      ENDIF
      L = J
!
!     If last element of array is less than than T, interchange with T
!
      IF (DX(J) .LT. T) THEN
         DX(IJ) = DX(J)
         DX(J) = T
         T = DX(IJ)
!
!        If first element of array is greater than T, interchange with T
!
         IF (DX(I) > T) THEN
            DX(IJ) = DX(I)
            DX(I) = T
            T = DX(IJ)
         ENDIF
      ENDIF
!
!     Find an element in the second half of the array which is smaller
!     than T
!
   40 L = L-1
      IF (DX(L) > T) GO TO 40
!
!     Find an element in the first half of the array which is greater
!     than T
!
   50 K = K+1
      IF (DX(K) .LT. T) GO TO 50
!
!     Interchange these elements
!
      IF (K <= L) THEN
         TT = DX(L)
         DX(L) = DX(K)
         DX(K) = TT
         GO TO 40
      ENDIF
!
!     Save upper and lower subscripts of the array yet to be sorted
!
      IF (L-I > J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70
!
!     Begin again on another portion of the unsorted array
!
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
!
   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1
!
   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = DX(I+1)
      IF (DX(I) .LE. T) GO TO 80
      K = I
!
   90 DX(K+1) = DX(K)
      K = K-1
      IF (T .LT. DX(K)) GO TO 90
      DX(K+1) = T
      GO TO 80
!
!     Sort DX and carry DY along
!
  100 M = 1
      I = 1
      J = NN
      R = 0.375D0
!
  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF
!
  120 K = I
!
!     Select a central element of the array and save it in location T
!
      IJ = I + INT((J-I)*R)
      T = DX(IJ)
      TY = DY(IJ)
      TY1 = DY1(IJ)
!
!     If first element of array is greater than T, interchange with T
!
      IF (DX(I) > T) THEN
         DX(IJ) = DX(I)
         DX(I) = T
         T = DX(IJ)
         DY(IJ) = DY(I)
         DY(I) = TY
         TY = DY(IJ)
         DY1(IJ) = DY1(I)
         DY1(I) = TY1
         TY1 = DY1(IJ)
      ENDIF
      L = J
!
!     If last element of array is less than T, interchange with T
!
      IF (DX(J) .LT. T) THEN
         DX(IJ) = DX(J)
         DX(J) = T
         T = DX(IJ)
         DY(IJ) = DY(J)
         DY(J) = TY
         TY = DY(IJ)
         DY1(IJ) = DY1(J)
         DY1(J) = TY1
         TY1 = DY1(IJ)
!
!        If first element of array is greater than T, interchange with T
!
         IF (DX(I) > T) THEN
            DX(IJ) = DX(I)
            DX(I) = T
            T = DX(IJ)
            DY(IJ) = DY(I)
            DY(I) = TY
            TY = DY(IJ)
            DY1(IJ) = DY1(I)
            DY1(I) = TY1
            TY1 = DY1(IJ)
         ENDIF
      ENDIF
!
!     Find an element in the second half of the array which is smaller
!     than T
!
  130 L = L-1
      IF (DX(L) > T) GO TO 130
!
!     Find an element in the first half of the array which is greater
!     than T
!
  140 K = K+1
      IF (DX(K) .LT. T) GO TO 140
!
!     Interchange these elements
!
      IF (K <= L) THEN
         TT = DX(L)
         DX(L) = DX(K)
         DX(K) = TT
         TTY = DY(L)
         DY(L) = DY(K)
         DY(K) = TTY
         TTY1 = DY1(L)
         DY1(L) = DY1(K)
         DY1(K) = TTY1
         GO TO 130
      ENDIF
!
!     Save upper and lower subscripts of the array yet to be sorted
!
      IF (L-I > J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160
!
!     Begin again on another portion of the unsorted array
!
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
!
  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1
!
  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = DX(I+1)
      TY = DY(I+1)
      TY1 = DY1(I+1)
      IF (DX(I) <= T) GO TO 170
      K = I
!
  180 DX(K+1) = DX(K)
      DY(K+1) = DY(K)
      DY1(K+1) = DY1(K)
      K = K-1
      IF (T .LT. DX(K)) GO TO 180
      DX(K+1) = T
      DY(K+1) = TY
      DY1(K+1) = TY1
      GO TO 170
!
!     Clean up
!
  190 IF (KFLAG <= -1) THEN
         DO 200 I=1,NN
            DX(I) = -DX(I)
  200    CONTINUE
      ENDIF
      RETURN
      END
