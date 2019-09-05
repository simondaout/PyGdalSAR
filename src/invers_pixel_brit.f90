!*************************************************
!
!   Inversion for each pixel independently as a function of time
!
!    Weighting and smoothing included
!
!    M.-P. Doin, CNRS,    Ecole d'ete des Houches, May 2011
!
!    This code is an updated mix of time series methods presented in the two following papers:
!
!    Cavalié, M.-P. Doin, C. Lasserre, P. Briole, Ground motion measurement in the Lake 
!    Mead area (Nevada, USA), by DInSAR time series analysis : probing the lithosphere 
!    rheological structure,  J. Geophys. Res., 112, B03403, doi:10.1029/2006JB004344, 2007
!
!    P. Lopez Quiroz, M.-P. Doin, F. Tupin, P. Briole, J.-M. Nicolas, 
!    Time series analysis of Mexico city subsidence constrained by radar Interferometry, 
!    Journal of Applied Geophysics, 69 (1), Sp. Iss. SI, 1-15, doi: 10.1016/j.jappgeo.2009.02.006, 2009
!
!***************************************************
! Attention :  internally coded parameters that may need to be changed  
!            maximum image number  nsmax = 300
!            maximum interferogram number  nnmax = 3000
!            threshold on local phase dispersion to eliminate outliers threshold = 1.0
!            mask on phase (or offsets) for input couples : abs(phase(is,irect,i))<0.00001_idp or phase(is,irect,i)<-9990_idp
!                                                               or NaN
!***********************************************************************
! Input files :
! List of images : YYYYMMDD time-time0 other_sort_variable PERP_BASELINE APS_Amplitude
! List of interferograms: YYYYMMDD YYYYMMDD COH_MOYENNE
!**********************************************************************

 module info_communes

  integer,   parameter :: IDP = kind(0.0E0)
  real(IDP),   parameter :: xstep = 1.
  real(IDP),   parameter :: ystep = 1.
  real(IDP),   parameter :: t_threshold=0.001
  real(IDP),   parameter :: nodata_output_value=9999.
  real(IDP) :: seuil_ferm
!  real(IDP),   parameter :: gam_liss2 = 0.
  real(IDP) :: gam_liss2
  character(len=8),   dimension(:), allocatable :: Im
  character(len=8),   dimension(:,:), allocatable :: Im1_Im2
  real,   dimension(:), allocatable :: Im_s
  integer,   dimension(:,:), allocatable :: Taille,Dec
  integer   :: Simage,Ninterf,irep1,irep2,irepIm,irep_comb,iliss,ipondliss
  integer   :: Simage_inv,Ninterf_inv,itri,iponder,icov,ibaseinv,ipondimage
  integer   :: ifuncinv,istep,icos
  integer   :: imed,ider_zero,ifermout,iqual,itype_ifg,iwcoh,itype_coh
  integer   :: ilin,iadd_var,iadd_lin,iadd_lin_inv,iadd_var_inv
  integer   :: nsamp_mask,corr_unw,ponder_rms,mask_rms,weight_rms,iter_rms_max,iter_unw_max
  real :: tstep
  real :: frac_interf,scale_rms,thres_rms
  integer,   dimension(4) :: coins
  real(IDP),   dimension(4) :: coinsll
  real,   dimension(:,:,:), allocatable :: phase,coh
  logical,  dimension(:,:), allocatable :: masque
  real,   dimension(:), allocatable :: shift,niveau_eau,base_im
  real,  dimension(:,:), allocatable :: mat_base
  real(IDP),   dimension(:,:), allocatable :: matd,mats
  logical,   dimension(:), allocatable :: flag_Im_inv,flag_Im1_Im2_comb
  logical,   dimension(:), allocatable :: flag_Interf_inv
  real(IDP),   dimension(:,:), allocatable :: mat
  real(IDP),   dimension(:,:), allocatable :: sigma_int_norm
  double precision,   dimension(:,:), allocatable :: sigma_int
  real(IDP),   dimension(:), allocatable :: pondim
  real,   dimension(:), allocatable :: varqual
  character(len=30) :: liste_image
  character(len=50) :: liste_interfero

 end module info_communes

module file_units

  integer :: unit_conf_in = 5, unit_conf_out = 6 ! may become 15 & 16

end module file_units

program invers_pixel
    use info_communes
    use file_units

    implicit none

    integer :: numargs
    character(len=128) :: arg

    integer   :: irefer

    numargs = iargc()
    if (numargs .eq. 1) then
      call getarg(1, arg)
      open(15,action='read',file=arg)
      unit_conf_in = 15
      open(16,action='write',file="/dev/null")
      unit_conf_out = 16
    endif

    write(unit_conf_out,*)'temporal smoothing weight, gamma liss **2 (if <0.0001, no smoothing) ?'
    read(unit_conf_in,*) gam_liss2
    write(unit_conf_out,*)'mask pixels that have high RMS misclosure ? (y=0;n=1)'
    read(unit_conf_in,*)ifermout
    write(unit_conf_out,*)'threshold for the mask on RMS misclosure (in same unit as input files) ?'
    read(unit_conf_in,*)seuil_ferm
    write(unit_conf_out,*)'range and azimuth downsampling (every n pixel)?'
    read(unit_conf_in,*)nsamp_mask
    write(unit_conf_out,*)'iterations to correct unwrapping errors ? (y:nb_of_iterations,n:0)'
    read(unit_conf_in,*)corr_unw
    write(unit_conf_out,*)'iterations to weight pixels of interferograms with large residual? (y:nb_of_iterations,n:0)'
    read(unit_conf_in,*) ponder_rms
    write(unit_conf_out,*)'Scaling value for weighting residuals (1/(res**2+value**2)) (in same unit as input files)?'
    write(unit_conf_out,*)'Must be approximately equal to standart deviation on measurement noise'
    read(unit_conf_in,*) scale_rms
    write(unit_conf_out,*)'iterations to mask (tiny weight) pixels of interferograms with large residual? (y:nb_of_iterations,n:0)'
    read(unit_conf_in,*) mask_rms
    write(unit_conf_out,*)'threshold on residual, defining clearly wrong values (in same unit as input files) ?'
    read(unit_conf_in,*) thres_rms
    iter_unw_max=corr_unw
    corr_unw=min(1,corr_unw)
    corr_unw=1-corr_unw
    iter_rms_max=max(ponder_rms,mask_rms)
    ponder_rms=min(1,ponder_rms)
    mask_rms=min(1,mask_rms)
    weight_rms=(1-ponder_rms)*(1-mask_rms)
    call listing()
    call taille_dec()
    call masque_commun()

    write(unit_conf_out,*)'referencement of interferogram by band (1) or corner (2) ?'
    read(unit_conf_in,*)irefer
    if(irefer==2)then
    call refer_interferos()
    else
    call refer_interferos_bande()
    endif
    call ponderation_variance()

    call matrice_base()

    call iter_px()
    stop
end

!**********************************************************************

subroutine iter_px

   use info_communes
   use file_units

   implicit none

   interface
     subroutine write_ts_envihdr(path, xsize, ysize, cnt, dates)
       character(*),intent(in)::path
       integer,intent(in)::xsize,ysize,cnt
       character(8),dimension(:),intent(in)::dates
     end subroutine
     subroutine write_rmspixel_envihdr(path, xsize, ysize)
       character(*),intent(in)::path
       integer,intent(in)::xsize,ysize
     end subroutine
   end interface

   integer      :: irec1,irec2,j,jlines,iwidth,irect,is,i
   integer      :: nread,jj,iss,nlig,ncol,nx,ny,ist,ireco,isss,ic
   integer      :: irugeau,ivarmodm,ivari
   integer      :: nxx,nyy
   integer      :: iter_unw,ich_unw,iter_rms,icol
   integer,  parameter :: nreadsimult=5
   integer, dimension(1) :: locmax
   real(IDP),   dimension(:), allocatable :: vect,vp,WORK,vectemp,vectsauv
   real(IDP),   dimension(:), allocatable :: vectcomp,vectrug
   real(IDP),   dimension(:,:), allocatable :: rugp,mat_sauv
   double precision, dimension(:), allocatable :: rmsi,rmsim
   double precision, dimension(:,:), allocatable :: rmsp
   double precision, dimension(:,:), allocatable :: rmspi,rmspim
   double precision :: rms,rugm,modm,varmodm,vari,rmsout
   integer,   dimension(:), allocatable :: IWORK,irmsi,irmsim
   integer,   dimension(:,:), allocatable :: irmsp,irugp,irmspim
   integer, dimension(:), allocatable :: loc_rms
   integer, dimension(:), allocatable :: anb_im_inv,anb_int_inv,rang_inv
   real(IDP),   dimension(:), allocatable :: deplac,pondim2,niveau_eau_diff,deplac_liss
   real,   dimension(:,:), allocatable :: correcDEM,varmod,rho,slope
   real,   dimension(:,:), allocatable :: season_cos,season_sin,step
   real, dimension(:), allocatable :: varip
   integer      :: nrhs,rank,LWORK,LIWORK,INFO,iadd_var_tmp
   real(IDP)    :: rcond,sumd,pondimm,rugeaum,vareau,deplm,niveau_eaum,covme
   real(IDP)    :: varmodmed,rugomed,slopemed,rhomed
   real,   dimension(:), allocatable :: xxx
   character(len=8) :: date1,date2
   character(len=30) :: nomRMS

   jlines=coins(4)-coins(2)+1
   iwidth=coins(3)-coins(1)+1
   nxx=int((iwidth-1.)/nsamp_mask)+1
   nyy=int((jlines-1.)/nsamp_mask)+1
   print*,'number of columns and lines for the inverted area :',nxx,nyy

   write(unit_conf_out,*)'Smoothing : scheme with 3pts (0) or 5pts (1) ?'
   read(unit_conf_in,*)iliss
   write(unit_conf_out,*)'Smoothing : weighting by the average time step (y :0 ; n : 1, int : 2) ?'
   read(unit_conf_in,*)ipondliss
   write(unit_conf_out,*)'Put the derivative of the first time step to zero ? (y :0 ; n : 1)?'
   read(unit_conf_in,*)ider_zero

   allocate(flag_Interf_inv(Ninterf),flag_Im_inv(Simage))
   allocate(deplac(Simage),deplac_liss(Simage))
   allocate(anb_im_inv(nxx),anb_int_inv(nxx),rang_inv(nxx))

   open(1,file='exemp_evol_temp',status='unknown')
   open(2,file='exemp_depl_cum',status='unknown')
   open(3,file='depl_cumule',status='unknown',form='unformatted',access='direct',recl=Simage*4)
   call write_ts_envihdr('depl_cumule', nxx, nyy-1, Simage, Im)
   
   if(ilin==1)then
     open(4,file='depl_cumule_liss',status='unknown',form='unformatted',access='direct',recl=Simage*4)
     call write_ts_envihdr('depl_cumule_liss', nxx, nyy-1, Simage, Im)
   endif

   open(31,file='nb_image_px',status='unknown',form='unformatted',access='direct',recl=4*nxx)
   open(32,file='nb_interf_px',status='unknown',form='unformatted',access='direct',recl=4*nxx)
   open(33,file='lien_manquants',status='unknown',form='unformatted',access='direct',recl=4*nxx)
   open(34,file='var_ini_inter',status='unknown',form='unformatted',access='direct',recl=4*nxx)
   open(35,file='RMS_interf_pixel',status='unknown',form='unformatted',access='direct',recl=4*nxx)

   if(ifuncinv>=1) then
     if(ibaseinv>0) then
         allocate(correcDEM(nxx,nyy-1))
         correcDEM(:,:)=0.
     endif
     if(istep>0) then
         allocate(step(nxx,nyy-1))
         step(:,:)=0.
     endif
     if(icos>0) then
         allocate(season_cos(nxx,nyy-1))
         season_cos(:,:)=0.
         allocate(season_sin(nxx,nyy-1))
         season_sin(:,:)=0.
     endif
   endif

! array pour calculer le rms
   allocate(rmsi(Ninterf),rmsp(nxx,nyy-1),irmsi(Ninterf),irmsp(nxx,nyy-1))
   allocate(rmspi(nxx,Ninterf),loc_rms(nxx),rmspim(nxx,Simage),irmspim(nxx,Simage))
   allocate(rmsim(Simage),irmsim(Simage))
   rmsi(:)=0.d0
   rmsim(:)=0.d0
   rmsp(:,:)=0.d0
   irmsp(:,:)=0
   irmsi(:)=0
   irmsim(:)=0

!ouverture des fichiers de RMS par interfero
   do is=1,Ninterf
    date1 = Im1_Im2(is,1)
    date2 = Im1_Im2(is,2)
    nomRMS='RMSpixel_'//date1//'_'//date2
    print*,' ',1000+is,' ',nomRMS
    open(unit=1000+is,file=nomRMS,status='unknown',form='unformatted',access='direct',recl=4*nxx)
    call write_rmspixel_envihdr(nomRMS, nxx, nyy-1)
   enddo
   do is=1,Simage
    date1 = Im(is)
    nomRMS='RMSpixel_'//date1
    print*,' ',500+is,' ',nomRMS
    open(unit=500+is,file=nomRMS,status='unknown',form='unformatted',access='direct',recl=4*nxx)
    call write_rmspixel_envihdr(nomRMS, nxx, nyy-1)
   enddo
   do is=1,Simage
    date1 = Im(is)
    nomRMS='nb_couple_'//date1
    print*,' ',100+is,' ',nomRMS
    open(unit=100+is,file=nomRMS,status='unknown',form='unformatted',access='direct',recl=4*nxx)
   enddo

!
   vari=0.d0
   ivari=0
   allocate(varmod(nxx,nyy-1),varip(nxx),rho(nxx,nyy-1),slope(nxx,nyy-1))
   rho(:,:)=0.
   slope(:,:)=0.

! array pour calculer la rugosite
   allocate(rugp(nxx,nyy-1),irugp(nxx,nyy-1))
   rugp(:,:)=0.
   irugp(:,:)=0
   irugeau=0
   varmod(:,:)=0.
   varmodm=0.d0
   modm=0.d0
   ivarmodm=0
   ireco=0

! boucle sur l ensemble des pixels de la region
   ny=0
   do j=1,jlines,nreadsimult*nsamp_mask
     print*,'j',j
     nread=int((jlines-j)/nsamp_mask)
     nread=min(nread,nreadsimult)
     irec1=j+coins(2)-1
     irec2=irec1+nread*nsamp_mask
     if(imed==0)then
       call ouv_lect_int_med(irec1,irec2,nsamp_mask,1,Ninterf)
     else
       call ouv_lect_int(irec1,irec2,nsamp_mask,1,Ninterf)
     endif

     do irect=1,nread
       ny=ny+1
       jj=j+(irect-1)*nsamp_mask
       rmspi(:,:)=0.d0
       rmspim(:,:)=0.d0
       irmspim(:,:)=0
       anb_im_inv(:)=0
       anb_int_inv(:)=0
       rang_inv(:)=15
       varip(:)=0.
       loc_rms(:)=0

       nx=0
       do i=1,iwidth,nsamp_mask
         nx=nx+1
         ireco=ireco+1

         if(masque(i,jj))then
           ! selection des interferos
           flag_Interf_inv(:)=.true.
           do is=1,Ninterf
             if(isnan(phase(is,irect,i))) then
               flag_Interf_inv(is)=.false.
             else
               if(abs(phase(is,irect,i))<0.00001_idp)flag_Interf_inv(is)=.false.
               if(phase(is,irect,i)<-9990_idp)flag_Interf_inv(is)=.false.
             endif
           enddo
           ! selection des images sur lesquelles on va inverser
           call select_image()

           if(Ninterf_inv<=Simage_inv)write(*,*)'pbe: less interferograms than images'

           iadd_lin_inv=0
           iadd_var_inv=0
           if(ilin>0)then
             iadd_lin_inv=Simage_inv
             iadd_var_inv=Simage_inv+ifuncinv
           endif
           if(ifuncinv>0.and.ilin==0)then
             iadd_var_inv=ifuncinv+1
             iadd_lin_inv=Simage_inv
           endif
           ncol=Simage_inv-1+iadd_var_inv
           nlig=Ninterf_inv+iadd_lin_inv+Simage_inv*ilin

          ! vecteur donnees d
          allocate(vect(nlig))
          vect(:)=0._idp
          iss=0
          do is=1,Ninterf
            if(flag_Interf_inv(is))then
              iss=iss+1
              vect(iss)=phase(is,irect,i)-shift(is)
            endif
          enddo
          allocate(vectsauv(Ninterf_inv))
          vectsauv(1:Ninterf_inv)= vect(1:Ninterf_inv)

          if(iponder==0) then
            allocate(sigma_int_norm(Ninterf_inv,Ninterf_inv))
          endif

          ! restriction de la matrice G 
          call matrice_px()

          ! On sauve la matrice G pour pouvoir iterer (en pondérant le résidu)
          if(corr_unw==0.or.weight_rms==0)then
            iter_unw=0
            iter_rms=0
            allocate(mat_sauv(nlig,ncol))
            mat_sauv(:,:)=mat(:,:)
          endif
33        continue

          ! ponderation pour que toutes les images aient un poids identique
          if(ipondimage==0)then
            iss=0
            allocate(pondim2(Ninterf_inv))
            do is=1,Ninterf
              if(flag_Interf_inv(is))then
                iss=iss+1
                pondim2(iss)=pondim(is)
              endif
            enddo
            pondimm=sum(pondim2(:))/Ninterf_inv
            pondim2(:)=pondim2(:)/pondimm
            do is=1,Ninterf_inv
              vect(is)=vect(is)*pondim2(is)
              mat(is,:)=mat(is,:)*pondim2(is)
            enddo
            deallocate(pondim2) 
          endif

          !ponderation par la cohérence
          if(iwcoh==0)then
            iss=0
            allocate(pondim2(Ninterf_inv))
            do is=1,Ninterf
              if(flag_Interf_inv(is))then
                iss=iss+1
                pondim2(iss)=coh(is,irect,i)+0.05
              endif
            enddo
            pondimm=sum(pondim2(:))/Ninterf_inv
            pondim2(:)=pondim2(:)/pondimm
            do is=1,Ninterf_inv
              vect(is)=vect(is)*pondim2(is)
              mat(is,:)=mat(is,:)*pondim2(is)
            enddo
            deallocate(pondim2)
          endif

          ! ponderation du vecteur donnees d par la matrice de ponderation sigma_int_norm
          if(iponder==0)then
            allocate(vectemp(Ninterf_inv))
            vectemp=matmul(sigma_int_norm(:,:),vect(1:Ninterf_inv))
            vect(1:Ninterf_inv)=vectemp(1:Ninterf_inv)
            deallocate(vectemp)
          endif

          ! nombre de vecteurs solution a trouver:
          nrhs=1
          ! mise a zero des valeurs propres < rcond * valeur propre max
          rcond=0.0000001_idp
          info=0
          LIWORK=64*ncol
          allocate(IWORK(LIWORK),WORK(1),vp(nlig))
          ! compute the minimum-norm solution to a real linear least squares problem
          ! methode divide and conquer for SVD
          call SGELSD(nlig,ncol,nrhs,mat,nlig,vect,nlig,vp,rcond, &
                      rank,WORK,-1,IWORK,INFO)
          LWORK=WORK(1)
          deallocate(WORK)
          allocate(WORK(LWORK))
          WORK(:)=0.
          call SGELSD(nlig,ncol,nrhs,mat,nlig,vect,nlig,vp,rcond, &
                      rank,WORK,LWORK,IWORK,INFO)
          if (INFO < 0) then
            write(0, *) "error: SGELSD() ",-INFO,"th argument is invalid"
            stop 255
          else if (INFO > 0) then
            write(0, *) "warning: SGELSD() did not converge"
          endif
          ! methode classique QR iteration?
          deallocate(IWORK,WORK,vp)

          do is=1,ncol
            if(isnan(vect(is)).or.(abs(vect(is)).gt.10000).or.info/=0)then
              deplac(:)=nodata_output_value
              write(3,rec=ireco)(deplac(iss),iss=1,Simage)
              if(ilin==1) write(4,rec=ireco)(deplac(iss),iss=1,Simage)
              deallocate(vect,mat,vectsauv)
              if(ilin>0)deallocate(matd)
              if(iponder==0) deallocate(sigma_int_norm)
              if(corr_unw==0.or.weight_rms==0) deallocate(mat_sauv)
              deallocate(mats)
              goto 100
            endif
          enddo

          ! Comparaison interfero original interfero synthetique : sur les premieres lignes
          ! de la matrice inversee. Verifier que le poids relatif des lignes de lissage est faible
          allocate(vectcomp(Ninterf_inv))
          vectcomp(:)=0._idp
          vectcomp(:)=matmul(mats(:,:),vect(1:Simage_inv-1))
          vectcomp(:)=vectcomp(:)-vectsauv(:)

          ! On change la phase de +/-2pi si l erreur de fermeture est superieure a 4.0 rad
          ! et on recommence
          !if(ich_unw==1)print*,'pixel apres unw',nx,ny,'interfero',is_unw,vectcomp(is_unw),vectsauv(is_unw)
          if(corr_unw==0.and.iter_unw.lt.iter_unw_max)then
            iter_unw=iter_unw+1
            ich_unw=0
            is=maxloc(abs(vectcomp(:)),1)
            !do is=1,Ninterf_inv
            ! if(abs(vectcomp(is)).gt.4)then
            !   print*,'pixel avant unw',nx,ny,'interfero',is,vectcomp(is),vectsauv(is)
            !   is_unw=is
            ! endif
            if(vectcomp(is).gt.4.5)then
              vectsauv(is)=vectsauv(is)+2.*3.1416
              ich_unw=1
            elseif(vectcomp(is).lt.-4.5)then
              vectsauv(is)=vectsauv(is)-2.*3.1416
              ich_unw=1
            endif
!         enddo
            if(ich_unw==1)then
            vect(:)=0.
            do is=1,Ninterf_inv
              vect(is)=vectsauv(is)
            enddo
            mat(:,:)=mat_sauv(:,:)
            deallocate(vectcomp)
            iter_rms=min(iter_rms,iter_rms_max)
            goto 33
          endif
        endif

        ! On calcule un poids inversement proportionnel a l erreur de fermeture et on recommence
        if(weight_rms==0.and.iter_rms.lt.iter_rms_max)then
          iter_rms=iter_rms+1
          vect(:)=0.
          do is=1,Ninterf_inv
            vect(is)=vectsauv(is)
          enddo
          mat(:,:)=mat_sauv(:,:)
          ! weight on each pixel of each interferogram depending on residue
          allocate(pondim2(Ninterf_inv))
          pondim2(:)=1.
          if(ponder_rms==1)then
            pondim2(:)=1./(scale_rms**2+vectcomp(:)**2)
          endif
          if(mask_rms==1)then
            do is=1,Ninterf_inv
              if(abs(vectcomp(is)).gt.thres_rms) pondim2(is)=pondim2(is)/1000.
            enddo
          endif
          pondimm=sum(pondim2(:))/Ninterf_inv
          pondim2(:)=pondim2(:)/pondimm

          do is=1,Ninterf_inv
            vect(is)=vect(is)*pondim2(is)
            mat(is,:)=mat(is,:)*pondim2(is)
          enddo
          deallocate(pondim2)
          deallocate(vectcomp)
          if(corr_unw==0)iter_unw=iter_rms-1
          goto 33
        endif

        if(corr_unw==0.or.weight_rms==0) deallocate(mat_sauv)

        anb_im_inv(nx)=Simage_inv
        anb_int_inv(nx)=Ninterf_inv
        rang_inv(nx)=Simage_inv-1+iadd_var_inv-rank

        if(ifuncinv>0) then
           iadd_var_tmp=0
           icol=Simage_inv-1
           if(ilin>0)icol=icol+Simage_inv
           if(ibaseinv>0) then
              iadd_var_tmp=iadd_var_tmp+1
              correcDEM(nx,ny)=vect(icol+iadd_var_tmp)
           endif
           if(istep>0)then
              iadd_var_tmp=iadd_var_tmp+1
              step(nx,ny)=vect(icol+iadd_var_tmp)
           endif
           if(icos>0)then
              iadd_var_tmp=iadd_var_tmp+1
              season_cos(nx,ny)=vect(icol+iadd_var_tmp)
              iadd_var_tmp=iadd_var_tmp+1
              season_sin(nx,ny)=vect(icol+iadd_var_tmp)
           endif
        endif

        ! rms of misclosure for each image
        do iss=1,Simage
          if(flag_Im_inv(iss))then
            ic=0
            isss=0
            do is=1,Ninterf
              if(flag_Interf_inv(is))then
                isss=isss+1
                date1 = Im1_Im2(is,1)
                date2 = Im1_Im2(is,2)
                if(Im(iss).eq.date1.or.Im(iss).eq.date2)then
                  ic=ic+1
                  rmspim(nx,iss)=rmspim(nx,iss)+vectcomp(isss)**2
                  rmsim(iss)=rmsim(iss)+vectcomp(isss)**2
                  irmsim(iss)=irmsim(iss)+1
                endif 
              endif
            enddo
            if(ic.gt.0)rmspim(nx,iss)=sqrt(rmspim(nx,iss)/ic)
            irmspim(nx,iss)=ic
          endif
        enddo

        ! residue for each interferogram
        iss=0
        do is=1,Ninterf
          if(flag_Interf_inv(is))then
            iss=iss+1
            rmspi(nx,is)=vectcomp(iss)
          endif
        enddo
        ! Rms of residues for each interferogram or for each pixel
        vectcomp(:)=vectcomp(:)**2
        rmsout=sqrt(sum(vectcomp(:))/Ninterf_inv)
        rmsp(nx,ny)=sum(vectcomp(:))
        locmax=maxloc(vectcomp(:))
        irmsp(nx,ny)=Ninterf_inv
        iss=0
        do is=1,Ninterf
          if(flag_Interf_inv(is))then
            iss=iss+1
            if(iss.eq.locmax(1))loc_rms(nx)=is
            rmsi(is)=rmsi(is)+vectcomp(iss)
            irmsi(is)=irmsi(is)+1
!           if(abs(vectcomp(iss)).gt.10.)then
!             print*,'pixel apres unw',nx,ny,'interfero',is,iss
!             print*,sqrt(vectcomp(iss)),vectsauv(iss)
!           endif
          endif
        enddo

        deallocate(mats)

        ! variance initiale des interferos
        vari=vari+sum(vectsauv(:)**2)
        varip(nx)=sum(vectsauv(:)**2)/Ninterf_inv
        varip(nx)=sqrt(varip(nx))
        ivari=ivari+Ninterf_inv

        deallocate(vectcomp,vectsauv)

        ! on masque les points qui ferment mal
        if(ifermout==0.and.rmsout>seuil_ferm)then
          deplac(:)=nodata_output_value
          write(3,rec=ireco)(deplac(is),is=1,Simage)
          if(ilin==1) write(4,rec=ireco)(deplac(is),is=1,Simage)
          deallocate(vect,mat)
          if(ilin>0)deallocate(matd)
          if(iponder==0) deallocate(sigma_int_norm)
          goto 100
        endif

        ! calcul de la rugosite de la solution
        if(ilin==1)then
          allocate(vectrug(Simage_inv))
          ! partie lissee de la solution
          vectrug(:)=matmul(matd(:,:),vect(Simage_inv:2*Simage_inv-1))*1.e12
          ! definition rugosite
          !vectrug(:)=vectrug(:)**2
          vectrug(:)=abs(vectrug(:))
          rugp(nx,ny)=sum(vectrug(:))
          irugp(nx,ny)=Simage_inv

          ! rugosite des variations de niveau d eau
          if(itri==0.and.Simage_inv==Simage.and.irugeau==0)then
          irugeau=1
          vectrug(:)=matmul(matd(:,:),niveau_eau(:))*1.e12
          ! definition rugosite
          !vectrug(:)=vectrug(:)**2
          vectrug(:)=abs(vectrug(:))
          rugeaum=SUM(vectrug(:))/Simage
          ! definition rugosite
          !rugeaum=rugeaum*1.e-4
          rugeaum=rugeaum*1.e-2
          vareau=SUM(niveau_eau(:)**2)/Simage - (SUM(niveau_eau(:))/Simage)**2
          vareau=sqrt(vareau)*1.e-2
        endif

        ! deallocation
        deallocate(vectrug,matd)

        ! Amplitude du modele, correlation du modele avec le niveau d eau
        allocate(niveau_eau_diff(Simage_inv))
        allocate(vectrug(Simage_inv))
        iss=0
        do is=1,Simage
          if(flag_Im_inv(is))then
            iss=iss+1
            niveau_eau_diff(iss)=niveau_eau(is)
          endif
        enddo
        niveau_eaum=sum(niveau_eau_diff(:))/Simage_inv
        vectrug(:)=vect(Simage_inv:2*Simage_inv-1)
        covme=0.
        do is=1,Simage_inv
          covme=covme+vectrug(is)*niveau_eau_diff(is)
        enddo
        niveau_eau_diff(:)=niveau_eau_diff(:)**2
        vareau=SUM(niveau_eau_diff(:))/Simage_inv -niveau_eaum**2
        vareau=sqrt(vareau)
        deplm=SUM(vectrug(:))/Simage_inv
        vectrug(:)=vectrug(:)**2
        varmod(nx,ny)=SUM(vectrug(:))/Simage_inv- deplm**2
        varmod(nx,ny)=sqrt(varmod(nx,ny))
        varmodm=varmodm+SUM(vectrug(:))
        covme=covme/Simage_inv -niveau_eaum*deplm
        !print*,niveau_eaum,vareau,deplm,varmod(nx,ny)
        ivarmodm=ivarmodm+Simage_inv
        modm=modm+deplm
        rho(nx,ny)=covme/(varmod(nx,ny)*vareau)
        slope(nx,ny)=covme/vareau**2
        deallocate(vectrug)
        deallocate(niveau_eau_diff)
      endif !fin if sur ilin=1

      ! ecriture du champ de deplacement inverse cumule
      iss=0
      sumd=0.
      deplac(1)=sumd
      ist=1
      do while(.not.flag_Im_inv(ist))
        deplac(ist)=nodata_output_value
        ist=ist+1
        deplac(ist)=sumd
      enddo
      do is=ist,Simage-1
        if(flag_Im_inv(is+1))then
          iss=iss+1
          sumd=sumd+vect(iss)
          deplac(is+1)=sumd
        else
          deplac(is+1)=nodata_output_value
        endif
      enddo
      write(3,rec=ireco)(deplac(is),is=1,Simage)
      if(ilin==1) then
        iss=0
        do is=1,Simage
          if(flag_Im_inv(is))then
            iss=iss+1
            deplac_liss(is)=vect(Simage_inv-1+iss)
          else
            deplac_liss(is)=nodata_output_value
          endif
        enddo
        write(4,rec=ireco)(deplac_liss(is),is=1,Simage)
      endif

      deallocate(vect,mat)
      if(iponder==0)deallocate(sigma_int_norm)
    else
      deplac(:)=nodata_output_value
      write(3,rec=ireco)(deplac(is),is=1,Simage)
      if(ilin==1) write(4,rec=ireco) (deplac(is),is=1,Simage)
    endif

100  continue

  enddo ! fin de boucle sur indice nx

  !ecriture des fichiers de RMS par interfero
  do is=1,Ninterf
    write(unit=1000+is,rec=ny)(real(rmspi(i,is)),i=1,nxx)
  enddo
  !ecriture des fichiers de RMS par date
  do is=1,Simage
    write(unit=500+is,rec=ny)(real(rmspim(i,is)),i=1,nxx)
  enddo
  !ecriture du nb de couple par date
  do is=1,Simage
    write(unit=100+is,rec=ny)(real(irmspim(i,is)),i=1,nxx)
  enddo

  ! ecriture du nb d image par pixel et du rang de la matrice
  write(31,rec=ny)(real(anb_im_inv(i)),i=1,nxx)
  write(32,rec=ny)(real(anb_int_inv(i)),i=1,nxx)
  write(33,rec=ny)(rang_inv(i),i=1,nxx)

  write(34,rec=ny)(varip(i),i=1,nxx)
  write(35,rec=ny)(real(loc_rms(i)),i=1,nxx)
 enddo ! fin de boucle sur indice irect

 deallocate(phase,coh)

   enddo

      write(*,*)nx,ny,nxx,nyy

   deallocate(flag_Interf_inv,flag_Im_inv,deplac,loc_rms)
   close(1)
   close(2)
   close(3)
   close(4)
   close(31)
   close(32)
   close(33)
   close(34)
   close(35)
!
!fermeture des fichiers de RMS par interfero
   do is=1,Ninterf
      close(unit=1000+is)
   enddo
   do is=1,Simage
      close(unit=500+is)
   enddo
   do is=1,Simage
      close(unit=100+is)
   enddo

   if(ipondimage==0)then
     deallocate(pondim)
   endif

!  ecriture du coefficient proportionnel a la baseline, ie lie a l erreur de DEM
   if(ifuncinv>=1)then
     if(ibaseinv>0)then
       open(4,file='correcDEM',status='unknown',form='unformatted',access='direct',recl=nx*4)
        do j=1,ny
         write(4,rec=j)(correcDEM(i,j),i=1,nx)
        enddo
       close(4)
       deallocate(correcDEM)
     endif
     if(istep>0)then
       open(4,file='step.r4',status='unknown',form='unformatted',access='direct',recl=nx*4)
        do j=1,ny
         write(4,rec=j)(step(i,j),i=1,nx)
        enddo
       close(4)
       deallocate(step)
     endif
     if(icos>0)then
       open(4,file='cos.r4',status='unknown',form='unformatted',access='direct',recl=nx*4)
        do j=1,ny
         write(4,rec=j)(season_cos(i,j),i=1,nx)
        enddo
       close(4)
       open(4,file='sin.r4',status='unknown',form='unformatted',access='direct',recl=nx*4)
        do j=1,ny
         write(4,rec=j)(season_sin(i,j),i=1,nx)
        enddo
       close(4)
       deallocate(season_cos,season_sin)
     endif
   endif

! ecriture des rms par interfero /par pixel
   rms=SUM(rmsp(:,:))/SUM(irmsp(:,:))
   rms=sqrt(rms)
   print*,'RMS global inversion en radian', rms
   rms=SUM(rmsi(:))/SUM(irmsi(:))
   rms=sqrt(rms)
   print*,'RMS global inversion en radian', rms
   open(4,file='RMSpixel',status='unknown',form='unformatted',access='direct',recl=nx*4)
!$OMP PARALLEL DO PRIVATE(i,j)
   do i=1,nx
     do j=1,ny
       if(irmsp(i,j)>10)then
         rmsp(i,j)=sqrt(rmsp(i,j)/irmsp(i,j))
       else
         rmsp(i,j)=0.d0
       endif
     enddo
   enddo
!$OMP END PARALLEL DO
   do j=1,ny
     write(4,rec=j)(real(rmsp(i,j)),i=1,nx)
    enddo
   close(4)

   open(4,file='RMSinterfero',status='unknown',form='formatted')
     do i=1,Ninterf
      rmsi(i)=sqrt(rmsi(i)/irmsi(i))
      write(4,*)i,Im1_Im2(i,1),' ',Im1_Im2(i,2),' ',rmsi(i)
     enddo
   close(4)
   open(4,file='RMSdate',status='unknown',form='formatted')
     do i=1,Simage
      rmsim(i)=sqrt(rmsim(i)/irmsim(i))
      write(4,*)i,Im(i),rmsim(i)
     enddo
   close(4)
   deallocate(rmsi,rmsp,irmsi,irmsp,rmsim,irmsim)

! ecriture du champ d amplitude du modele
    modm=modm/ivarmodm
    varmodm=sqrt(varmodm/ivarmodm -modm**2)
    print*,'amplitude moyenne du modele en rad',varmodm
    open(4,file='amp_modele',status='unknown',form='unformatted',access='direct',recl=nx*4)
    do j=1,ny
     write(4,rec=j)(varmod(i,j),i=1,nx)
    enddo
    close(4)
    is=0
    varmodm=0
    do i=1,nx
       do j=1,ny
         if(varmod(i,j)>0.00001)then
           is=is+1
           varmodm=varmodm+varmod(i,j)
         endif
       enddo
    enddo
    varmodm=varmodm/is
    print*,'Average model (cumulated phase) amplitude',varmodm
    allocate(xxx(nx*ny/4))
    is=0
    do i=nx/3,(2*nx)/3
       do j=ny/3,(2*ny/3)
         if(varmod(i,j)>0.00001)then
           is=is+1
           xxx(is)=varmod(i,j)
         endif
       enddo
    enddo
    call median(xxx,is,varmodmed)
    print*,'mediane des  amplitudes du modele sur la zone centrale',varmodmed

! ecriture du champ de rugosite
    rugm=SUM(rugp(:,:))/SUM(irugp(:,:))
    print*,'average rugosity in rad / (1.e8 s)**2 ',rugm
!    print*,'rugosite moyenne en rad**2 / (1.e8 s)**4 ',rugm
    print*,'average rugosity en 1 / (1.e8 s)**2 ',rugm/varmodm
!    print*,'rugosite moyenne en 1 / (1.e8 s)**4 ',rugm/varmodm**2
   open(4,file='rugosite',status='unknown',form='unformatted',access='direct',recl=nx*4)
   open(5,file='rugosite2',status='unknown',form='unformatted',access='direct',recl=nx*4)
    is=0
    rugm=0.
     do i=1,nx
       do j=1,ny
         if(irugp(i,j).gt.1)then
           rugp(i,j)=rugp(i,j)/irugp(i,j)
         endif
         if(varmod(i,j)>0.00001)then
           is=is+1
           varmod(i,j)=rugp(i,j)/varmod(i,j)
           rugm=rugm+varmod(i,j)
         endif
       enddo
     enddo
    do j=1,ny
     write(4,rec=j)(rugp(i,j),i=1,nx)
     write(5,rec=j)(varmod(i,j),i=1,nx)
    enddo
   close(4)
   close(5)
   rugm=rugm/is
   print*,'rugosite moyenne en 1 / (1.e8 s)**2 (chaque pixel est /par amp) ',rugm
!   print*,'rugosite moyenne en 1 / (1.e8 s)**4 (chaque pixel est /par amp) ',rugm
   is=0
    do i=nx/3,(2*nx)/3
       do j=ny/3,(2*ny/3)
         if(varmod(i,j)>0.00001)then
           is=is+1
           xxx(is)=varmod(i,j)
         endif
       enddo
    enddo
    call median(xxx,is,rugomed)
    print*,'mediane des rugosite en 1 / (1.e8 s)**2 ds partie centrale',rugomed

   open(4,file='coeff_corr_var',status='unknown',form='unformatted',access='direct',recl=nx*4)
    do j=1,ny
     write(4,rec=j)(rho(i,j),i=1,nx)
    enddo
   close(4)
   is=0
    do i=nx/3,(2*nx)/3
       do j=ny/3,(2*ny/3)
         if(abs(rho(i,j))>0.00001)then
           is=is+1
           xxx(is)=rho(i,j)
         endif
       enddo
    enddo
    call median(xxx,is,rhomed)
    print*,'mediane of correlation coefficient in central part',rhomed

   open(4,file='slope_corr_var',status='unknown',form='unformatted',access='direct',recl=nx*4)
    do j=1,ny
     write(4,rec=j)(slope(i,j),i=1,nx)
   enddo
   close(4)
   is=0
    do i=nx/3,(2*nx)/3
       do j=ny/3,(2*ny/3)
         if(abs(slope(i,j))>0.00001)then
           is=is+1
           xxx(is)=slope(i,j)
         endif
       enddo
    enddo
    call median(xxx,is,slopemed)
    print*,'mediane des pentes rad/yr ds partie centrale',slopemed


   deallocate(irugp,rugp,rho,slope)
   deallocate(varmod,varip,xxx)
   deallocate(anb_im_inv,rang_inv,anb_int_inv)

   vari=dsqrt(vari/ivari)
   print*,'ecart type initial des interferos en rad',vari

   write(*,*)gam_liss2,real(rms),real(varmodm),real(rugm),varmodmed,rugomed,rhomed,slopemed

   open(111,file='lect.in',status='unknown')
   write(111,*)nx,ny,nxx,nyy-1
   write(111,*)coins(1),coins(2)
   write(111,*)nsamp_mask
   write(111,*)Simage
   close(111)

   deallocate(Im_s,Im,Im1_Im2,Taille,Dec,masque,shift,mat_base)
   deallocate(niveau_eau)

   end subroutine

!************************************************************************

   subroutine matrice_px

   use info_communes

   implicit none

   integer   :: ncol,is,iss,nlig,info
   real(IDP) :: idt1,idt2,idt3
   real :: dxm2,dxm1,dx1,dx2,adt,bdt,cdt,ddt
   integer   :: ibeg,is_avt,is_avt2,is_apr2
   real,  dimension(:,:), allocatable :: mat_base2
   real,   dimension(:), allocatable :: Im_s_temp
   real(IDP),   dimension(:,:), allocatable :: sigma_int_norm2,mat_temp
   real(IDP),   dimension(:), allocatable :: val_prop
   real(IDP) :: moy_prop,valmin,pondliss,dtm2
   character(len=1) :: UPLO,DIAG

   ncol=Simage_inv-1+iadd_var_inv
! 3pts decale sur les bords
   nlig=Ninterf_inv+iadd_lin_inv+Simage_inv*ilin
!   nlig=Ninterf_inv+iadd_lin_inv+Simage_inv-2
   allocate(mat(nlig,ncol),Im_s_temp(Simage_inv),val_prop(Ninterf_inv))
   if(ilin>0)allocate(matd(Simage_inv,Simage_inv))

! si probleme de valeur propre negative, on reprend ici
! soit on ajoute la valeur propre negative minimale calculee aux termes diagnonaux
! et on retourne en 10
! soit on multiplie les termes diagonaux par valmin (on augmente valmin jusqu a ce que ca passe)
   valmin=1.
   if(iponder==0.and.icov==0)then
   valmin=1.60
   endif
!   valmin=0.
!   if(iponder==0.and.icov==0)then
!   valmin=-4.
!   endif
10 continue
   if(iponder==0.and.icov==0)then
   valmin=valmin*1.1
!   valmin=-4.
   endif

   allocate(mat_base2(Ninterf_inv+iadd_lin_inv,Simage-1+iadd_var))
   mat_base2(:,:)=0._idp

! on recopie dans mat et dans sigma_int_norm les lignes et col
! correspondant aux interferos et deplacements OK
   
   if(iponder==0) then
       allocate(sigma_int_norm2(Ninterf_inv,Ninterf))
   endif

! suppression des lignes correspondants aux interferos manquants
   iss=0
   do is=1,Ninterf
     if(flag_Interf_inv(is))then
       iss=iss+1
       mat_base2(iss,:)=mat_base(is,:)
       if(iponder==0) then
         sigma_int_norm2(iss,:)=sigma_int(is,:)
       endif
     endif
   enddo 

! suppression des lignes additionnelles correspondant aux images manquantes
   if(iadd_lin>0)then
   iss=0
   do is=1,Simage
     if(flag_Im_inv(is))then
       iss=iss+1
       mat_base2(Ninterf_inv+iss,:)=mat_base(Ninterf+is,:)
     endif
   enddo
   endif

! suppression des colonnes correspondants aux images manquantes
   iss=0
   do is=1,Simage-1
     if(flag_Im_inv(is))then
       iss=iss+1
       Im_s_temp(iss)=Im_s(is)
       if(iss<=Simage_inv-1)then
       mat(1:Ninterf_inv+iadd_lin_inv,iss)=mat_base2(1:Ninterf_inv+iadd_lin_inv,is)
       endif
      endif
   enddo
   if(flag_Im_inv(Simage))Im_s_temp(iss+1)=Im_s(Simage)

! recopie les colonnes avec les variables additionnelles
   if(ifuncinv>0.and.ilin==0)then
     mat(1:Ninterf_inv+iadd_lin_inv,Simage_inv:Simage_inv+iadd_var_inv-1)=mat_base2(:,Simage:Simage+iadd_var-1)
   endif

! suppression des colonnes avec les variables additionnelles  correspondants aux images manquantes
   if(ilin>0)then
   iss=0
   do is=1,Simage
     if(flag_Im_inv(is))then
       iss=iss+1
       mat(1:Ninterf_inv+iadd_lin_inv,Simage_inv-1+iss)=mat_base2(1:Ninterf_inv+iadd_lin_inv,Simage-1+is)
     endif
   enddo
   if(ifuncinv>0)mat(1:Ninterf_inv+iadd_lin_inv,2*Simage_inv:2*Simage_inv-1+ifuncinv)=mat_base2(:,2*Simage:2*Simage+ifuncinv) 
   endif

! on sauve la matrice G
   allocate(mats(Ninterf_inv,Simage_inv-1))
   mats(:,:)=mat(1:Ninterf_inv,1:Simage_inv-1)

! suppression des colonnes  correspondant aux interferos manquants ds sigma_int_norm
   if(iponder==0)then
   iss=0
   do is=1,Ninterf
     if(flag_Interf_inv(is))then
       iss=iss+1
       sigma_int_norm(:,iss)=sigma_int_norm2(:,is)
     endif
   enddo
   do is=1,Ninterf_inv
     sigma_int_norm(is,is)=sigma_int_norm(is,is)*valmin
!     sigma_int_norm(is,is)=sigma_int_norm(is,is)-valmin
   enddo
   endif

   deallocate(mat_base2)
   if(iponder==0) deallocate(sigma_int_norm2)

   if(iponder==0.and.icov==0)then
! METHODE 1
! decompostion de la matrice de covariance sigma_int_norm en P D P**T
! du coup sigma_d ** -1 = P D**-1 P**T
! on garde la matrice  D**-1/2 P**T et on normalise par la moyenne de D**-1/2
!
!    JOBZ='V'
!    UPLO='U'
!    lwork= -1
!    liwork=-1
!    allocate(WORK(1),IWORK(1))
!    call SSYEVD(JOBZ,UPLO,Ninterf_inv,sigma_int_norm,Ninterf_inv,val_prop,WORK,LWORK,IWORK,LIWORK,INFO)
!    lwork=WORK(1)
!    liwork=IWORK(1)
!    deallocate(WORK,IWORK)
!    allocate(WORK(lwork),IWORK(liwork))
!    call SSYEVD(JOBZ,UPLO,Ninterf_inv,sigma_int_norm,Ninterf_inv,val_prop,WORK,LWORK,IWORK,LIWORK,INFO)
!    deallocate(WORK,IWORK)
!    if(info/=0)print*,'SSYEVD failed', info
!    if(val_prop(1)<0.) then
!!         print*,'val_prop negative',val_prop(1),val_prop(2),val_prop(3)
!         valmin=val_prop(1)*1.005
!         goto 10
!    endif
!    allocate(sigma_int_norm2(Ninterf_inv,Ninterf_inv))
!    sigma_int_norm2=transpose(sigma_int_norm)
!    sigma_int_norm=sigma_int_norm2
!    deallocate(sigma_int_norm2)
!!
!! normalisation
!!     print*,'val_prop'
!!     print*,val_prop
!    moy_prop=0.
!    do is=1,Ninterf_inv
!      val_prop(is)=1./sqrt(val_prop(is))
!      moy_prop=moy_prop+val_prop(is)
!    enddo
!    moy_prop=moy_prop/Ninterf_inv
!    do is=1,Ninterf_inv
!      sigma_int_norm(is,:)=sigma_int_norm(is,:)*val_prop(is)/moy_prop
!    enddo
!    deallocate(val_prop)

! METHODE 2
! decomposition de Cholevsky de la matrice de covariance sigma_int_norm
! et inversion de la matrice puis transpose
   UPLO='U'
   call spotrf(UPLO,Ninterf_inv,sigma_int_norm,Ninterf_inv,info)
   if(abs(info)>0) then
        write(*,*)'error decompostion de Cholevsky'
        write(*,*)'valmin',valmin
        goto 10
   endif

!  sigma_int_norm est une matrice triangulaire sup
   do is=2,Ninterf_inv
      sigma_int_norm(is,1:is-1)=0.
   enddo

   DIAG='N'
   call strtri(UPLO,DIAG,Ninterf_inv,sigma_int_norm,Ninterf_inv,info)
   if(abs(info)>0) then
       write(*,*)'error inversion matrice triangul sup'
       print*,'valmin',valmin
       goto 10
   endif

! transpose, sigma_int_norm devient triangulaire inf
   do is=1,Ninterf_inv-1
   do iss=is+1,Ninterf_inv
      sigma_int_norm(iss,is)=sigma_int_norm(is,iss)
      sigma_int_norm(is,iss)=0.
   enddo
   enddo
! normalisation
    moy_prop=0.
    do is=1,Ninterf_inv
      moy_prop=moy_prop+sigma_int_norm(is,is)
    enddo
    moy_prop=moy_prop/Ninterf_inv
    do is=1,Ninterf_inv
      sigma_int_norm(is,:)=sigma_int_norm(is,:)/moy_prop
    enddo

   endif

! ponderation de la matrice G par A = (U**-1)**T, avec la matrice de cov
! egale a (U**T U)
   if(iponder==0)then
      allocate(mat_temp(Ninterf_inv,Simage_inv-1))
      mat_temp(1:Ninterf_inv,1:Simage_inv-1)=matmul(sigma_int_norm(1:Ninterf_inv,1:Ninterf_inv),mat(1:Ninterf_inv,1:Simage_inv-1))
      mat(1:Ninterf_inv,1:Simage_inv-1)=mat_temp(1:Ninterf_inv,1:Simage_inv-1)
      deallocate(mat_temp)
   endif

! lissage : minimisation de la derivee seconde de la LOS
! on ajoute autant de lignes que d images
! ici, applique au deplacement cumule et non aux increments de deplacement
   if(ilin>0)then

   ibeg=Ninterf_inv+iadd_lin_inv
   mat(ibeg+1:nlig,:)=0._idp
   pondliss=0._idp
   matd(1:Simage_inv,:)=0._idp

   if(iliss==0)then

!3pts + 3pts decales sur bords : autant de lignes que d images
   is=1
   if(ider_zero==0)then
    idt1=Im_s_temp(is+1)-Im_s_temp(is)
    matd(is,is)=-1./idt1**2
    matd(is,is+1)=1./idt1**2
    if(ipondliss==0)then
      dtm2=idt1**2
      pondliss=pondliss+1.
      matd(is,is)=matd(is,is)*dtm2
      matd(is,is+1)=matd(is,is+1)*dtm2
    elseif(ipondliss==2)then
      dtm2=abs(idt1)
      pondliss=pondliss+1./dtm2
      matd(is,is)=matd(is,is)*dtm2
      matd(is,is+1)=matd(is,is+1)*dtm2
    else
      dtm2=idt1**2
      pondliss=pondliss+1./dtm2
    endif
   else
    idt1=Im_s_temp(is+1)-Im_s_temp(is)
    idt2=Im_s_temp(is+2)-Im_s_temp(is)
    idt3=Im_s_temp(is+2)-Im_s_temp(is+1)
    matd(is,is)=-(1./idt2/idt3-1./idt1/idt3)
    matd(is,is+1)=-1./idt1/idt3
    matd(is,is+2)= 1./idt2/idt3
    if(ipondliss==0)then
        dtm2=((idt1+idt3)/2.)**2
        pondliss=pondliss+1.
        matd(is,is)=matd(is,is)*dtm2
        matd(is,is+1)=matd(is,is+1)*dtm2
        matd(is,is+2)=matd(is,is+2)*dtm2
    elseif(ipondliss==2)then
        dtm2=abs((idt1+idt3)/2.)
        pondliss=pondliss+1./dtm2
        matd(is,is)=matd(is,is)*dtm2
        matd(is,is+1)=matd(is,is+1)*dtm2
        matd(is,is+2)=matd(is,is+2)*dtm2
    else
      dtm2=((idt1+idt3)/2.)**2
      pondliss=pondliss+1./dtm2
    endif
   endif
!
   do is=2,Simage_inv-1
     idt1=Im_s_temp(is)-Im_s_temp(is-1)
     idt2=Im_s_temp(is+1)-Im_s_temp(is-1)
     idt3=Im_s_temp(is+1)-Im_s_temp(is)
     is_avt=is-1
     if(idt1.lt.t_threshold)then
       is_avt=is-2
       idt1=Im_s_temp(is)-Im_s_temp(is-2)
       idt2=Im_s_temp(is+1)-Im_s_temp(is-2)
     endif
     if(idt3.lt.t_threshold)then
       if(ipondliss==0)then
         matd(is,is)=-3.
         matd(is,is+1)=3.
       elseif(ipondliss==2)then
         matd(is,is)=-3./idt3
         matd(is,is+1)=3./idt3
       else
         matd(is,is)=-1/idt3**2
         matd(is,is+1)=1./idt3**2
       endif
     else
     matd(is,is_avt)=1./idt2/idt1
     matd(is,is)=-1./idt2/idt1-1./idt2/idt3
     matd(is,is+1)=1./idt2/idt3
     if(ipondliss==0)then
        dtm2=((idt1+idt3)/2.)**2
        pondliss=pondliss+1.
        matd(is,is_avt)=matd(is,is_avt)*dtm2
        matd(is,is)=matd(is,is)*dtm2
        matd(is,is+1)=matd(is,is+1)*dtm2
     elseif(ipondliss==2)then
        dtm2=abs((idt1+idt3)/2.)
        pondliss=pondliss+1./dtm2
        matd(is,is_avt)=matd(is,is_avt)*dtm2
        matd(is,is)=matd(is,is)*dtm2
        matd(is,is+1)=matd(is,is+1)*dtm2
     else
        dtm2=((idt1+idt3)/2.)**2
        pondliss=pondliss+1./dtm2
     endif
     endif
   enddo
   is=Simage_inv
   idt1=Im_s_temp(is)-Im_s_temp(is-1)
   idt2=Im_s_temp(is)-Im_s_temp(is-2)
   idt3=Im_s_temp(is-1)-Im_s_temp(is-2)
   matd(is,is-1)=-1./idt1/idt3
   matd(is,is)=(-1./idt2/idt3+1./idt1/idt3)
   matd(is,is-2)=1./idt2/idt3
   if(ipondliss==0)then
        dtm2=((idt1+idt3)/2.)**2
        pondliss=pondliss+1.
        matd(is,is-1)=matd(is,is-1)*dtm2
        matd(is,is-2)=matd(is,is-2)*dtm2
        matd(is,is)=matd(is,is)*dtm2
   elseif(ipondliss==2)then
        dtm2=abs((idt1+idt3)/2.)
        pondliss=pondliss+1./dtm2
        matd(is,is-1)=matd(is,is-1)*dtm2
        matd(is,is-2)=matd(is,is-2)*dtm2
        matd(is,is)=matd(is,is)*dtm2
    else
      dtm2=((idt1+idt3)/2.)**2
      pondliss=pondliss+1./dtm2
   endif
   pondliss=pondliss/Simage_inv
   matd(:,:)=matd(:,:)/pondliss

!lissage a 5pts
   else

   is=1
   if(ider_zero==0)then
    idt1=Im_s_temp(is+1)-Im_s_temp(is)
    matd(is,is)=-1./idt1**2
    matd(is,is+1)=1./idt1**2
    if(ipondliss==0)then
      dtm2=idt1**2
      pondliss=pondliss+1.
      matd(is,is)=matd(is,is)*dtm2
      matd(is,is+1)=matd(is,is+1)*dtm2
    elseif(ipondliss==2)then
      dtm2=abs(idt1)
      pondliss=pondliss+1./dtm2
      matd(is,is)=matd(is,is)*dtm2 
      matd(is,is+1)=matd(is,is+1)*dtm2 
    else
      dtm2=idt1**2
      pondliss=pondliss+1./dtm2
    endif
   else
! 3 pts decales
    idt1=Im_s_temp(is+1)-Im_s_temp(is)
    idt2=Im_s_temp(is+2)-Im_s_temp(is)
    idt3=Im_s_temp(is+2)-Im_s_temp(is+1)
    matd(is,is)=-(1./idt2/idt3-1./idt1/idt3)
    matd(is,is+1)=-1./idt1/idt3
    matd(is,is+2)=1./idt2/idt3
    if(ipondliss==0)then
        dtm2=((idt1+idt3)/2.)**2
        pondliss=pondliss+1.
        matd(is,is)=matd(is,is)*dtm2
        matd(is,is+1)=matd(is,is+1)*dtm2
        matd(is,is+2)=matd(is,is+2)*dtm2
    elseif(ipondliss==2)then
        dtm2=abs((idt1+idt3)/2.)
        pondliss=pondliss+1./dtm2
        matd(is,is)=matd(is,is)*dtm2
        matd(is,is+1)=matd(is,is+1)*dtm2
        matd(is,is+2)=matd(is,is+2)*dtm2
    else
      dtm2=((idt1+idt3)/2.)**2
      pondliss=pondliss+1./dtm2
    endif
   endif

!3 pts centres
   is=2
   idt1=Im_s_temp(is)-Im_s_temp(is-1)
   idt2=Im_s_temp(is+1)-Im_s_temp(is-1)
   idt3=Im_s_temp(is+1)-Im_s_temp(is)
   matd(is,is-1)=1./idt2/idt1
   matd(is,is)=-1./idt2/idt1-1./idt2/idt3
   matd(is,is+1)=1./idt2/idt3
   if(ipondliss==0)then
        dtm2=((idt1+idt3)/2.)**2
        pondliss=pondliss+1.
        matd(is,is-1)=matd(is,is-1)*dtm2
        matd(is,is)=matd(is,is)*dtm2
        matd(is,is+1)=matd(is,is+1)*dtm2
   elseif(ipondliss==2)then
        dtm2=abs((idt1+idt3)/2.)
        pondliss=pondliss+1./dtm2
        matd(is,is-1)=matd(is,is-1)*dtm2
        matd(is,is)=matd(is,is)*dtm2
        matd(is,is+1)=matd(is,is+1)*dtm2
    else
      dtm2=((idt1+idt3)/2.)**2
      pondliss=pondliss+1./dtm2
   endif

! 5pts
   do is=3,Simage_inv-2
     idt3=Im_s_temp(is+1)-Im_s_temp(is)
     if(idt3.lt.t_threshold)then
       if(ipondliss==0)then
         matd(is,is)=-3.
         matd(is,is+1)=3.
       elseif(ipondliss==2)then
         matd(is,is)=-3./idt3
         matd(is,is+1)=3./idt3
       else
         matd(is,is)=-1/idt3**2
         matd(is,is+1)=1./idt3**2
       endif
     else
     idt1=Im_s_temp(is)-Im_s_temp(is-1)
     is_avt=is-1
     is_avt2=is-2
     if(idt1.lt.t_threshold)then
       is_avt=is-2
       is_avt2=is-3
     endif
     idt1=Im_s_temp(is_avt)-Im_s_temp(is_avt2)
     if(idt1.lt.t_threshold)then
       is_avt2=is_avt2-1
     endif
     idt1=Im_s_temp(is+2)-Im_s_temp(is+1)
     is_apr2=is+2
     if(idt1.lt.t_threshold)then
       is_apr2=is+3
     endif
     dxm2=Im_s_temp(is_avt2)-Im_s_temp(is)
     dxm1=Im_s_temp(is_avt)-Im_s_temp(is)
     dx1=Im_s_temp(is+1)-Im_s_temp(is)
     dx2=Im_s_temp(is_apr2)-Im_s_temp(is)
     adt=-(dxm1*dx1+dxm1*dx2+dx1*dx2)/(dxm2*(dx2-dxm2)*(dx1-dxm2)*(dxm1-dxm2))
     bdt=-(dxm2*dx1+dxm2*dx2+dx1*dx2)/(dxm1*(dx2-dxm1)*(dx1-dxm1)*(dxm2-dxm1))
     cdt=-(dxm2*dxm1+dxm2*dx2+dxm1*dx2)/(dx1*(dx2-dx1)*(dxm1-dx1)*(dxm2-dx1))
     ddt=-(dxm2*dxm1+dxm2*dx1+dxm1*dx1)/(dx2*(dx1-dx2)*(dxm1-dx2)*(dxm2-dx2))
     matd(is,is_avt2)=adt
     matd(is,is_avt)=bdt
     matd(is,is)=-(adt+bdt)-(cdt+ddt)
     matd(is,is+1)=cdt
     matd(is,is_apr2)=ddt
     if(ipondliss==0)then
        dtm2=((Im_s_temp(is_apr2)-Im_s_temp(is_avt2))/4.)**2
        pondliss=pondliss+1.
        matd(is,is_avt2)=matd(is,is_avt2)*dtm2
        matd(is,is_avt)=matd(is,is_avt)*dtm2
        matd(is,is)=matd(is,is)*dtm2
        matd(is,is+1)=matd(is,is+1)*dtm2
        matd(is,is_apr2)=matd(is,is_apr2)*dtm2
     elseif(ipondliss==2)then
        dtm2=abs((Im_s_temp(is_apr2)-Im_s_temp(is_avt2))/4.)
        pondliss=pondliss+1./dtm2
        matd(is,is_avt2)=matd(is,is_avt2)*dtm2
        matd(is,is_avt)=matd(is,is_avt)*dtm2
        matd(is,is)=matd(is,is)*dtm2
        matd(is,is+1)=matd(is,is+1)*dtm2
        matd(is,is_apr2)=matd(is,is_apr2)*dtm2
     else
        dtm2=((Im_s_temp(is_apr2)-Im_s_temp(is_avt2))/4.)**2
        pondliss=pondliss+1./dtm2
     endif
     endif
   enddo

   is=Simage_inv-1
   idt1=Im_s_temp(is)-Im_s_temp(is-1)
   idt2=Im_s_temp(is+1)-Im_s_temp(is-1)
   idt3=Im_s_temp(is+1)-Im_s_temp(is)
   matd(is,is-1)=1./idt2/idt1
   matd(is,is)=-1./idt2/idt1-1./idt2/idt3
   matd(is,is+1)=1./idt2/idt3
   if(ipondliss==0)then
        dtm2=((idt1+idt3)/2.)**2
        pondliss=pondliss+1.
        matd(is,is-1)=matd(is,is-1)*dtm2
        matd(is,is)=matd(is,is)*dtm2
        matd(is,is+1)=matd(is,is+1)*dtm2
   elseif(ipondliss==2)then
        dtm2=abs((idt1+idt3)/2.)
        pondliss=pondliss+1./dtm2
        matd(is,is-1)=matd(is,is-1)*dtm2
        matd(is,is)=matd(is,is)*dtm2
        matd(is,is+1)=matd(is,is+1)*dtm2
    else
      dtm2=((idt1+idt3)/2.)**2
      pondliss=pondliss+1./dtm2
   endif

   is=Simage_inv
   idt1=Im_s_temp(is)-Im_s_temp(is-1)
   idt2=Im_s_temp(is)-Im_s_temp(is-2)
   idt3=Im_s_temp(is-1)-Im_s_temp(is-2)
   matd(is,is-1)=-1./idt1/idt3
   matd(is,is)=(-1./idt2/idt3+1./idt1/idt3)
   matd(is,is-2)=1./idt2/idt3
   if(ipondliss==0)then
        dtm2=((idt1+idt3)/2.)**2
        pondliss=pondliss+1.
        matd(is,is-1)=matd(is,is-1)*dtm2
        matd(is,is-2)=matd(is,is-2)*dtm2
        matd(is,is)=matd(is,is)*dtm2
   elseif(ipondliss==2)then
        dtm2=abs((idt1+idt3)/2.)
        pondliss=pondliss+1./dtm2
        matd(is,is-1)=matd(is,is-1)*dtm2
        matd(is,is-2)=matd(is,is-2)*dtm2
        matd(is,is)=matd(is,is)*dtm2
    else
      dtm2=((idt1+idt3)/2.)**2
      pondliss=pondliss+1./dtm2
   endif
   pondliss=pondliss/Simage_inv
   matd(:,:)=matd(:,:)/pondliss

   endif

   mat(ibeg+1:ibeg+Simage_inv,Simage_inv:2*Simage_inv-1)=matd(:,:)*gam_liss2

   endif

   deallocate(Im_s_temp,val_prop)

   return
   end subroutine

!************************************************************************

   subroutine matrice_base

   use info_communes
   use file_units

   implicit none

   character(len=8) :: date1,date2
   integer      :: is,i,idate1,idate2,iinv,isens,idatemin,iadd_var_tmp
   integer,  dimension(:), allocatable :: idate
   real :: ieau1,ieau2,pi
   real(kind=8) :: isgn

   pi=acos(-1.)

   ifuncinv=0
   ibaseinv=0
   icos=0
   istep=0
!   inversion avec lissage en temps pour stabiliser le SVD  (y:1;n:0)
   ilin=0
   if(gam_liss2>=0.0001)ilin=1
!  ajout de fonctions en sus du lissage
   write(unit_conf_out,*)'Adjust functions to phase history ? (y:1;n:0) Require to use smoothing option (smoothing coefficient) !'
   read(unit_conf_in,*) ifuncinv
   if(ifuncinv ==1) then
     write(unit_conf_out,*)'compute DEM error proportional to perpendicular baseline ? (y:1;n:0)'
     read(unit_conf_in,*) ibaseinv
     write(unit_conf_out,*)'include a step function ? (y:1;n:0)'
     read(unit_conf_in,*) istep, tstep
     write(unit_conf_out,*)'include a cosinus / sinus function ? (y:1;n:0)'
     read(unit_conf_in,*) icos
     icos=icos*2
     ifuncinv =ibaseinv+istep+icos
   endif
   iadd_var=0
   iadd_lin=0
! on ajoute le coeff prop. a baseline perpendiculaire
   if(ilin>0)then
! on ajoute les lignes des Simage variables lissees 
     iadd_lin=Simage
! on ajoute les Simage variables lissees
     iadd_var=Simage+ifuncinv
   endif
   if(ifuncinv>0.and.ilin==0)then
       iadd_var=ifuncinv+1
       iadd_lin=Simage
   endif
   allocate(idate(Simage),mat_base(Ninterf+iadd_lin,Simage-1+iadd_var))

   mat_base(:,:)=0

! classement par date
   if(itri==0)then

   do is=1,Simage
      date1 = Im(is)
      read(date1,*)idate(is)
!      print*,idate(is)
   enddo

   do is=1,Ninterf
      date1 = Im1_Im2(is,1)
      date2 = Im1_Im2(is,2)
      read(date1,*)idate1
      read(date2,*)idate2
      idatemin=min(idate1,idate2)
      isens=1
      if(idatemin.eq.idate2)isens=-1
      do i=1,Simage-1
       isgn=(idate(i)-idate1)
       isgn=isgn*(idate(i)-idate2)
       if(isgn<0)mat_base(is,i)=isens
       if(idate(i)==idatemin)mat_base(is,i)=isens
      enddo
!      write(*,*)idate1,idate2,(mat_base(is,i),i=1,Simage-1)      
   enddo

! classement par niveau d eau
   else

   do is=1,Simage
      date1 = Im(is)
      read(date1,*)idate(is)
   enddo

   mat_base(:,:)=0

   do is=1,Ninterf
      date1 = Im1_Im2(is,1)
      date2 = Im1_Im2(is,2)
      read(date1,*)idate1
      read(date2,*)idate2
      do i=1,Simage
        if(idate1==idate(i))ieau1=Im_s(i)
        if(idate2==idate(i))ieau2=Im_s(i)
        iinv=1
        if(ieau1>ieau2)iinv=-1
      enddo
      do i=1,Simage-1
       isgn=(Im_s(i)-ieau1)
       isgn=isgn*(Im_s(i)-ieau2)
       if(isgn<0)mat_base(is,i)=iinv
       if(ieau1<ieau2.and.Im_s(i)==ieau1)mat_base(is,i)=1
       if(ieau1>ieau2.and.Im_s(i)==ieau2)mat_base(is,i)=iinv
      enddo
!      write(*,*)ieau1,ieau2,(mat_base(is,i),i=1,Simage-1)
   enddo

   endif

! on ajoute une contrainte delai cumule = somme des increments de delais (- alpha Bperp)
! ou somme des increments de delais = alpha Bperp +cte
! avec un poids faible: ne change les increments que si liens manquants
!
   if(iadd_lin>0)then
     do is=2,Simage
       do i=1,is-1
        mat_base(is+Ninterf,i)=1._idp
       enddo
     enddo
     iadd_var_tmp=0
     if(ilin==1)then
       iadd_var_tmp=Simage
       do is=1,Simage
         mat_base(is+Ninterf,Simage-1+is)=-1._idp
       enddo
     endif
     if(ifuncinv>0)then
       if(ibaseinv>0)then
         iadd_var_tmp=iadd_var_tmp+1
         do is=1,Simage
             mat_base(is+Ninterf,Simage-1+iadd_var_tmp)=-base_im(is)/100.
         enddo
       endif
       if(istep>0)then
         iadd_var_tmp=iadd_var_tmp+1
         do is=1,Simage
             mat_base(is+Ninterf,Simage-1+iadd_var_tmp)=0.
             if(Im_s(is)>tstep)mat_base(is+Ninterf,Simage-1+iadd_var_tmp)=-1
         enddo
       endif
       if(icos>0)then
         iadd_var_tmp=iadd_var_tmp+1
         do is=1,Simage
             mat_base(is+Ninterf,Simage-1+iadd_var_tmp)=-cos(pi*2.*Im_s(is))
         enddo
         iadd_var_tmp=iadd_var_tmp+1
         do is=1,Simage
             mat_base(is+Ninterf,Simage-1+iadd_var_tmp)=-sin(pi*2.*Im_s(is))
         enddo
       endif
       if(ilin==0)then
         iadd_var_tmp=iadd_var_tmp+1
         do is=1,Simage
           mat_base(is+Ninterf,Simage-1+iadd_var_tmp-1)=-1.
         enddo
       endif
     endif
! ponderation : il faut que le poids de ces lignes soit petit / lignes avec inversion des interferos
     mat_base(Ninterf+1:Ninterf+iadd_lin,:)=mat_base(Ninterf+1:Ninterf+iadd_lin,:)*0.0005_idp
! ponderation par la qualite de l image
     if(iqual==0)then
     do is=1,Simage
       mat_base(Ninterf+is,:)=mat_base(Ninterf+is,:)*varqual(is)
     enddo
     endif
   endif

   deallocate(idate)
   if(iqual==0)deallocate(varqual)

   end subroutine

!************************************************************************

   subroutine select_image

   use info_communes

   implicit none

   character(len=8) :: date1,date2
   integer      :: is,i,niter
   integer,   dimension(:), allocatable :: Im_num

   allocate(Im_num(Simage))

   niter=0
3  continue
   do i=1,Simage
     Im_num(i)=0
   do is=1,Ninterf
      if(flag_Interf_inv(is))then
      date1 = Im1_Im2(is,1)
      date2 = Im1_Im2(is,2)
      if(Im(i).eq.date1.or.Im(i).eq.date2)then
        Im_num(i)=Im_num(i)+1
      endif
      endif
   enddo
!   write(*,*)'image ',Im_ini(i),' vue ',Im_num_ini(i),' fois'
   enddo

   if(irepIm.gt.1.and.niter.lt.10)then
   do i=1,Simage
   if(Im_num(i)<irepIm)then
        do is=1,Ninterf
      if(flag_Interf_inv(is))then
      date1 = Im1_Im2(is,1)
      date2 = Im1_Im2(is,2)
      if(Im(i).eq.date1.or.Im(i).eq.date2)then
        flag_Interf_inv(is)=.false.
      endif
      endif
        enddo
    endif
   enddo
   niter=niter+1
   goto 3
   endif

   Simage_inv=COUNT(Im_num>=irepIm)
   do i=1,Simage
     if(Im_num(i)>=irepIm)then
        flag_Im_inv(i)=.true.
      else
        flag_Im_inv(i)=.false.
      endif
   enddo
   Ninterf_inv=COUNT(flag_Interf_inv)
!   write(*,*)'nb d images ',Simage_inv,' nb d interfero ',Ninterf_inv

   deallocate(Im_num)

   end subroutine select_image
!************************************************************************
subroutine refer_interferos_bande

   use info_communes
   use file_units

   implicit none

   integer      :: irec1,irec2,is,irect,iwidth,i,j,k,ndata,jlines,jj,iref
   integer      :: itr,jtr,itl,jtl,ibr,jbr,ibl,jbl,nread
   integer,   parameter :: nsamp_ref = 6
   integer,  parameter :: nreadsimult=20
   integer,   dimension(:), allocatable :: liminf,limsup
   real(IDP),   dimension(:,:), allocatable :: moy_ref,sigma_ref
   real, dimension(:)  :: val(5000)
   real :: xval

   write(unit_conf_out,*) 'interferogram referencement'
   write(unit_conf_out,*) 'band NW -SW(1), band SW- SE (2), band NW-NE (3),  '
   write(unit_conf_out,*) 'or three band average (4) or no referencement (5) ?'
   read(unit_conf_in,*) iref

!  Position des 4 coins
   iwidth=coins(3)-coins(1)+1
   jtl=nsamp_mask+1
   i=1
   do while(.not.masque(i,jtl).and.i.lt.iwidth)
     i=i+nsamp_mask
   enddo
   itl=i

   jbr=1+nsamp_mask*int((coins(4)-coins(2))/nsamp_mask)-nsamp_mask
   i=iwidth
   do while(.not.masque(i,jbr))
     i=i-nsamp_mask
   enddo
   ibr=i 

   jlines= coins(4)-coins(2)+1
   itr=1
   do j=1,jlines,nsamp_mask
   i=iwidth
   do while(.not.masque(i,j).and.i.gt.1)
     i=i-nsamp_mask
   enddo
   if(i.gt.itr)then
     itr=i
     jtr=j
   endif
   enddo
   jtr=jtr+nsamp_mask
   i=iwidth
   do while(.not.masque(i,jtr).and.i.gt.1)
     i=i-nsamp_mask
   enddo
   itr=i

   ibl=iwidth
   do j=1,jlines,nsamp_mask
     i=1
     do while(.not.masque(i,j).and.i.lt.iwidth)
       i=i+nsamp_mask
     enddo
     if(i.lt.ibl)then
       ibl=i
       jbl=j
     endif
   enddo
   jbl=max(jbl-nsamp_mask, 1)
   i=1
   do while(.not.masque(i,jbl).and.i.lt.iwidth)
     i=i+nsamp_mask
   enddo
   ibl=i

   write(*,*)'coins NE,NW,SE,SW'
   write(*,*)itr,itl,ibr,ibl
   write(*,*)jtr,jtl,jbr,jbl

   allocate(moy_ref(Ninterf,3),sigma_ref(Ninterf,3))

   open(1,file='refere_interfer',status='unknown')

   if(iref==1.or.iref==4)then
   write(1,*)'> referencement bande NW -SW'
   write(*,*)' referencement bande NW - SW'

   jlines=jbl-jtl+1

   allocate(liminf(jlines),limsup(jlines))

   do j=1,jlines
   jj=j+jtl-1
   liminf(j)=int((jj*1.-jtl)*(ibl*1.-itl)/(jbl*1.-jtl)+itl)
   limsup(j)=liminf(j)+66
   enddo

   do is=1,Ninterf

    val(1:5000)=0.
    moy_ref(is,1)=0._idp
    sigma_ref(is,1)=0._idp
    ndata=0

    do j=1,jlines,nreadsimult*nsamp_ref
      nread=int((jlines-j)/nsamp_ref)
      nread=min(nread,nreadsimult)
      irec1=j+jtl-1+coins(2)-1
      irec2=irec1+nread*nsamp_ref
      if(imed==0)then
        call ouv_lect_int_med(irec1,irec2,nsamp_ref,is,is)
      else
        call ouv_lect_int(irec1,irec2,nsamp_ref,is,is)
      endif

      do irect=1,nread
        k=(irect-1)*nsamp_ref+j
        do i=liminf(k),limsup(k),nsamp_ref
          if(abs(phase(1,irect,i))>0.00001.and.phase(1,irect,i)>-9990.)then
            moy_ref(is,1)=moy_ref(is,1)+phase(1,irect,i)
            sigma_ref(is,1)=sigma_ref(is,1)+phase(1,irect,i)**2
            ndata=ndata+1
            val(ndata)=phase(1,irect,i)
          endif
        enddo
      enddo
      deallocate(phase,coh)
    enddo

    if(ndata.gt.100)then
      sigma_ref(is,1)=sigma_ref(is,1)-moy_ref(is,1)**2/ndata
      sigma_ref(is,1)=sqrt(sigma_ref(is,1)/ndata)
      moy_ref(is,1)=moy_ref(is,1)/ndata
      call median(val,ndata,xval)
    else
      moy_ref(is,1)=0._idp
      sigma_ref(is,1)=0._idp
      xval=0.
    endif
    write(1,*)Im1_Im2(is,1),' ',Im1_Im2(is,2),ndata,          &
           moy_ref(is,1),sigma_ref(is,1),xval
!    write(*,*)Im1_Im2(is,1),' ',Im1_Im2(is,2),ndata,          &
!           moy_ref(is,1),sigma_ref(is,1),xval
    moy_ref(is,1)=xval

   enddo

   deallocate(liminf,limsup)

   endif

   if(iref==2.or.iref==4)then
   write(1,*)'> referencement bande SW- SE'
   write(*,*)' referencement bande SW- SE'

   jlines=jbr-jbl+1

   allocate(liminf(jlines),limsup(jlines))

   do j=1,jlines
   jj=j+jbl-1
   liminf(j)=int((jj*1.-jbl)*(ibr*1.-ibl)/(jbr*1.-jbl)+ibl)
   limsup(j)=liminf(j)+300
   limsup(j)=min(limsup(j),iwidth)
   enddo

   do is=1,Ninterf

    val(1:5000)=0.
    moy_ref(is,2)=0._idp
    sigma_ref(is,2)=0._idp
    ndata=0

    do j=1,jlines,nreadsimult*nsamp_ref
      nread=int((jlines-j)/nsamp_ref)
      nread=min(nread,nreadsimult)
      irec1=j+jbl-1+coins(2)-1
      irec2=irec1+nread*nsamp_ref
      if(imed==0)then
        call ouv_lect_int_med(irec1,irec2,nsamp_ref,is,is)
      else
        call ouv_lect_int(irec1,irec2,nsamp_ref,is,is)
      endif

      do irect=1,nread
        k=(irect-1)*nsamp_ref+j
        do i=liminf(k),limsup(k),nsamp_ref
          if(abs(phase(1,irect,i))>0.00001.and.phase(1,irect,i)>-9990.)then
            moy_ref(is,2)=moy_ref(is,2)+phase(1,irect,i)
            sigma_ref(is,2)=sigma_ref(is,2)+phase(1,irect,i)**2
            ndata=ndata+1
            val(ndata)=phase(1,irect,i)
          endif
        enddo
      enddo
      deallocate(phase,coh)
    enddo

    if(ndata.gt.100)then
      sigma_ref(is,2)=sigma_ref(is,2)-moy_ref(is,2)**2/ndata
      sigma_ref(is,2)=sqrt(sigma_ref(is,2)/ndata)
      moy_ref(is,2)=moy_ref(is,2)/ndata
      call median(val,ndata,xval)
    else
      moy_ref(is,2)=0._idp
      sigma_ref(is,2)=0._idp
      xval=0.
    endif
    write(1,*)Im1_Im2(is,1),' ',Im1_Im2(is,2),ndata,          &
           moy_ref(is,2),sigma_ref(is,2),xval
!    write(*,*)Im1_Im2(is,1),' ',Im1_Im2(is,2),ndata,          &
!           moy_ref(is,2),sigma_ref(is,2),xval
    moy_ref(is,2)=xval

   enddo

   deallocate(liminf,limsup)

   endif

   if(iref==3.or.iref==4)then
   write(1,*)'> referencement bande NW- NE'
   write(*,*)' referencement bande NW- NE'

   jlines=jtr-jtl+1

   allocate(liminf(jlines),limsup(jlines))

   do j=1,jlines
   jj=j+jtl-1
   limsup(j)=int((jj*1.-jtl)*(itr*1.-itl)/(jtr*1.-jtl)+itl)
   liminf(j)=limsup(j)-300
   liminf(j)=max(1,liminf(j))
   enddo

   do is=1,Ninterf

    val(1:5000)=0.
    moy_ref(is,2)=0._idp
    sigma_ref(is,2)=0._idp
    ndata=0

    do j=1,jlines,nreadsimult*nsamp_ref
      nread=int((jlines-j)/nsamp_ref)
      nread=min(nread,nreadsimult)
      irec1=j+jtl-1+coins(2)-1
      irec2=irec1+nread*nsamp_ref
      if(imed==0)then
        call ouv_lect_int_med(irec1,irec2,nsamp_ref,is,is)
      else
        call ouv_lect_int(irec1,irec2,nsamp_ref,is,is)
      endif

      do irect=1,nread
        k=(irect-1)*nsamp_ref+j
        do i=liminf(k),limsup(k),nsamp_ref
          if(abs(phase(1,irect,i))>0.00001.and.phase(1,irect,i)>-9990.)then
            moy_ref(is,2)=moy_ref(is,2)+phase(1,irect,i)
            sigma_ref(is,2)=sigma_ref(is,2)+phase(1,irect,i)**2
            ndata=ndata+1
            val(ndata)=phase(1,irect,i)
          endif
        enddo
      enddo
      deallocate(phase,coh)
    enddo

    if(ndata.gt.100)then
      sigma_ref(is,2)=sigma_ref(is,2)-moy_ref(is,2)**2/ndata
      sigma_ref(is,2)=sqrt(sigma_ref(is,2)/ndata)
      moy_ref(is,2)=moy_ref(is,2)/ndata
      call median(val,ndata,xval)
    else
      moy_ref(is,2)=0._idp
      sigma_ref(is,2)=0._idp
      xval=0.
    endif
    write(1,*)Im1_Im2(is,1),' ',Im1_Im2(is,2),ndata,          &
           moy_ref(is,2),sigma_ref(is,2),xval
    write(*,*)Im1_Im2(is,1),' ',Im1_Im2(is,2),ndata,          &
           moy_ref(is,2),sigma_ref(is,2),xval
    moy_ref(is,2)=xval

   enddo

   deallocate(liminf,limsup)

   endif

   close(1)

   allocate(shift(Ninterf))

   if(iref==1)then
       shift(1:Ninterf)=moy_ref(1:Ninterf,1)
   elseif(iref==2)then
       shift(1:Ninterf)=moy_ref(1:Ninterf,2)
   elseif(iref==3)then
       shift(1:Ninterf)=moy_ref(1:Ninterf,3)
   elseif(iref==4)then
       shift(1:Ninterf)=(moy_ref(1:Ninterf,1)+      &
                         moy_ref(1:Ninterf,2)+moy_ref(1:Ninterf,3))/3.
   else
      shift(1:Ninterf)=0.
   endif

   deallocate(moy_ref,sigma_ref)

end subroutine refer_interferos_bande

!************************************************************************

subroutine refer_interferos

   use info_communes
   use file_units

   implicit none

   integer      :: irec1,irec2,is,irect,iwidth,i,j,k,ndata,jlines,jj,iref
   integer      :: i0
   integer,   parameter :: nsamp_ref = 1
   integer,   dimension(:), allocatable :: liminf,limsup
   real(IDP),   dimension(:,:), allocatable :: moy_ref,sigma_ref
   real, dimension(:)  :: val(5000)
   real :: xval

   write(unit_conf_out,*) 'Interferogram referencement'
   write(unit_conf_out,*) 'corner NW (1), SE (2), average of the 2 (3), corner SW (4), '
   write(unit_conf_out,*) 'corner NE (5), average of corners NW, SW, NE (6) or no referencement (7) ?'
   read(unit_conf_in,*) iref

   allocate(moy_ref(Ninterf,4),sigma_ref(Ninterf,4))

   open(1,file='refere_interfer',status='unknown')

   if(iref==1.or.iref==3.or.iref==6)then
   write(1,*)'> referencement coin NW'
   write(*,*)' referencement coin NW'

   irec1=coins(2)+int(33/nsamp_mask)*nsamp_mask
   irec2=irec1+66
   irect=int((irec2-irec1)/nsamp_ref)+1
   iwidth=coins(3)-coins(1)+1
   jlines=irec2-irec1+1

   allocate(liminf(jlines+nsamp_mask),limsup(jlines+nsamp_mask))

   do j=1,jlines,nsamp_mask
   jj=j+irec1-coins(2)
   i=1
   do while(.not.masque(i,jj).and.i.lt.iwidth)
     i=i+nsamp_mask
   enddo
   liminf(j)=i-1+33
   limsup(j)=liminf(j)+66
   if(nsamp_mask>1)then
   do i=1,nsamp_mask-1
     liminf(j+i)=liminf(j)
     limsup(j+i)=limsup(j)
   enddo
   endif
   enddo

   do is=1,Ninterf

      if(imed==0)then
        call ouv_lect_int_med(irec1,irec2,nsamp_ref,is,is)
      else
        call ouv_lect_int(irec1,irec2,nsamp_ref,is,is)
      endif

      val(1:5000)=0.
      moy_ref(is,1)=0._idp
      sigma_ref(is,1)=0._idp
      ndata=0
      do j=1,irect
        k=(j-1)*nsamp_ref+1
        do i=liminf(k),limsup(k),nsamp_ref
          if(abs(phase(1,j,i))>0.00001.and.phase(1,j,i)>-9990.)then
            moy_ref(is,1)=moy_ref(is,1)+phase(1,j,i)
            sigma_ref(is,1)=sigma_ref(is,1)+phase(1,j,i)**2
            ndata=ndata+1
            val(ndata)=phase(1,j,i)
          endif
        enddo
      enddo
      if(ndata.gt.100)then
      sigma_ref(is,1)=sigma_ref(is,1)-moy_ref(is,1)**2/ndata
      sigma_ref(is,1)=sqrt(sigma_ref(is,1)/ndata)
      moy_ref(is,1)=moy_ref(is,1)/ndata
      call median(val,ndata,xval)
      else
      moy_ref(is,1)=0._idp
      sigma_ref(is,1)=0._idp
      xval=0.
      endif
      write(1,*)Im1_Im2(is,1),' ',Im1_Im2(is,2),ndata,          &
           moy_ref(is,1),sigma_ref(is,1),xval
!      write(*,*)Im1_Im2(is,1),' ',Im1_Im2(is,2),ndata,          &
!           moy_ref(is,1),sigma_ref(is,1),xval
      moy_ref(is,1)=xval

      deallocate(phase,coh)

   enddo

   deallocate(liminf,limsup)
 
   endif

   if(iref==2.or.iref==3)then
   write(1,*)'> referencement coin SE'
   write(*,*)' referencement coin SE'

   irec2=nsamp_mask*int((coins(4)-coins(2))/nsamp_mask)+coins(2)
   irec2=irec2-4*nsamp_mask
   irec1=irec2-9*nsamp_mask
   irect=int((irec2-irec1)/nsamp_ref)+1
   iwidth=coins(3)-coins(1)+1
   jlines=irec2-irec1+1

   allocate(liminf(jlines+nsamp_mask),limsup(jlines+nsamp_mask))

   do j=1,jlines,nsamp_mask
   jj=j+irec1-coins(2)
   i=iwidth
   do while(.not.masque(i,jj))
     i=i-nsamp_mask
   enddo
   limsup(j)=i-1-33
   liminf(j)=limsup(j)-66
   if(nsamp_mask>1)then
   do i=1,nsamp_mask-1
     liminf(j+i)=liminf(j)
     limsup(j+i)=limsup(j)
   enddo
   endif
   enddo

   do is=1,Ninterf

      if(imed==0)then
        call ouv_lect_int_med(irec1,irec2,nsamp_ref,is,is)
      else
        call ouv_lect_int(irec1,irec2,nsamp_ref,is,is)
      endif

      val(1:5000)=0.
      moy_ref(is,2)=0._idp
      sigma_ref(is,2)=0._idp
      ndata=0
      do j=1,irect
        k=(j-1)*nsamp_ref+1
        do i=liminf(k),limsup(k),nsamp_ref
          if(abs(phase(1,j,i))>0.00001.and.phase(1,j,i)>-9990.)then
            moy_ref(is,2)=moy_ref(is,2)+phase(1,j,i)
            sigma_ref(is,2)=sigma_ref(is,2)+phase(1,j,i)**2
            ndata=ndata+1
            val(ndata)=phase(1,j,i)
          endif
        enddo
      enddo
      if(ndata.gt.100)then
      sigma_ref(is,2)=sigma_ref(is,2)-moy_ref(is,2)**2/ndata
      sigma_ref(is,2)=sqrt(sigma_ref(is,2)/ndata)
      moy_ref(is,2)=moy_ref(is,2)/ndata
      call median(val,ndata,xval)
      else
         moy_ref(is,2)=0._idp
         sigma_ref(is,2)=0._idp
         xval=0.
      endif
      write(1,*)Im1_Im2(is,1),' ',Im1_Im2(is,2),ndata,moy_ref(is,2),sigma_ref(is,2),xval
!      write(*,*)Im1_Im2(is,1),' ',Im1_Im2(is,2),ndata,moy_ref(is,2),sigma_ref(is,2),xval
      moy_ref(is,2)=xval

      deallocate(phase,coh)

   enddo

   deallocate(liminf,limsup)
 
   endif

   if(iref==4.or.iref==6)then
   write(1,*)'> referencement coin SW'
   write(*,*)' referencement coin SW'

   iwidth=coins(3)-coins(1)+1
   jlines= coins(4)-coins(2)+1
   i0=iwidth
   do j=1,jlines,nsamp_mask
   i=1
   do while(.not.masque(i,j).and.i.lt.iwidth)
     i=i+nsamp_mask
   enddo
   if(i.lt.i0)then
     i0=i
     jj=j
   endif
   enddo

   irec1=coins(2)+jj-int(33/nsamp_mask)*nsamp_mask
   irec2=irec1+66
   irect=int((irec2-irec1)/nsamp_ref)+1
   jlines=irec2-irec1+1

   allocate(liminf(jlines+nsamp_mask),limsup(jlines+nsamp_mask))

   do j=1,jlines,nsamp_mask
   jj=j+irec1-1-coins(2)
   i=1
   do while(.not.masque(i,jj).and.i.lt.iwidth)
     i=i+nsamp_mask
   enddo
   liminf(j)=i-1+33
   limsup(j)=liminf(j)+66
   if(nsamp_mask>1)then
   do i=1,nsamp_mask-1
     liminf(j+i)=liminf(j)
     limsup(j+i)=limsup(j)
   enddo
   endif
   enddo

   do is=1,Ninterf

      if(imed==0)then
        call ouv_lect_int_med(irec1,irec2,nsamp_ref,is,is)
      else
        call ouv_lect_int(irec1,irec2,nsamp_ref,is,is)
      endif

      val(1:5000)=0.
      moy_ref(is,3)=0._idp
      sigma_ref(is,3)=0._idp
      ndata=0
      do j=1,irect
        k=(j-1)*nsamp_ref+1
        do i=liminf(k),limsup(k),nsamp_ref
          if(abs(phase(1,j,i))>0.00001.and.phase(1,j,i)>-9990.)then
            moy_ref(is,3)=moy_ref(is,3)+phase(1,j,i)
            sigma_ref(is,3)=sigma_ref(is,3)+phase(1,j,i)**2
            ndata=ndata+1
            val(ndata)=phase(1,j,i)
          endif
        enddo
      enddo
      if(ndata.gt.0)then
      sigma_ref(is,3)=sigma_ref(is,3)-moy_ref(is,3)**2/ndata
      sigma_ref(is,3)=sqrt(sigma_ref(is,3)/ndata)
      moy_ref(is,3)=moy_ref(is,3)/ndata
      call median(val,ndata,xval)
      else
      moy_ref(is,3)=0._idp
      sigma_ref(is,3)=0._idp
      xval=0.
      endif
      write(1,*)Im1_Im2(is,1),' ',Im1_Im2(is,2),ndata,          &
           moy_ref(is,3),sigma_ref(is,3),xval
      moy_ref(is,3)=xval

      deallocate(phase,coh)

   enddo

   deallocate(liminf,limsup)

   endif

   if(iref==5.or.iref==6)then
   write(1,*)'> referencement coin NE'
   write(*,*)' referencement coin NE'

   iwidth=coins(3)-coins(1)+1
   jlines= coins(4)-coins(2)+1
   i0=1
   do j=1,jlines,nsamp_mask
   i=iwidth
   do while(.not.masque(i,j).and.i.gt.1)
     i=i-nsamp_mask
   enddo
   if(i.gt.i0)then
     i0=i
     jj=j
   endif
   enddo

   irec1=coins(2)+jj-int(33/nsamp_mask)*nsamp_mask
   irec2=irec1+66
   irect=int((irec2-irec1)/nsamp_ref)+1
   jlines=irec2-irec1+1

   allocate(liminf(jlines+nsamp_mask),limsup(jlines+nsamp_mask))

   do j=1,jlines,nsamp_mask
   jj=j+irec1-1-coins(2)
   i=iwidth
   do while(.not.masque(i,jj).and.i.gt.1)
     i=i-nsamp_mask
   enddo
   limsup(j)=i-1-33
   liminf(j)=limsup(j)-66
   if(nsamp_mask>1)then
   do i=1,nsamp_mask-1
     liminf(j+i)=liminf(j)
     limsup(j+i)=limsup(j)
   enddo
   endif
   enddo

   do is=1,Ninterf

      if(imed==0)then
        call ouv_lect_int_med(irec1,irec2,nsamp_ref,is,is)
      else
        call ouv_lect_int(irec1,irec2,nsamp_ref,is,is)
      endif

      val(1:5000)=0.
      moy_ref(is,4)=0._idp
      sigma_ref(is,4)=0._idp
      ndata=0
      do j=1,irect
        k=(j-1)*nsamp_ref+1
        do i=liminf(k),limsup(k),nsamp_ref
          if(abs(phase(1,j,i))>0.00001.and.phase(1,j,i)>-9990.)then
            moy_ref(is,4)=moy_ref(is,4)+phase(1,j,i)
            sigma_ref(is,4)=sigma_ref(is,4)+phase(1,j,i)**2
            ndata=ndata+1
            val(ndata)=phase(1,j,i)
          endif
        enddo
      enddo
      if(ndata.gt.0)then
      sigma_ref(is,4)=sigma_ref(is,4)-moy_ref(is,4)**2/ndata
      sigma_ref(is,4)=sqrt(sigma_ref(is,4)/ndata)
      moy_ref(is,4)=moy_ref(is,4)/ndata
      call median(val,ndata,xval)
      else
      moy_ref(is,4)=0._idp
      sigma_ref(is,4)=0._idp
      xval=0.
      endif
      write(1,*)Im1_Im2(is,1),' ',Im1_Im2(is,2),ndata,          &
           moy_ref(is,4),sigma_ref(is,4),xval
      moy_ref(is,4)=xval

      deallocate(phase,coh)

   enddo

   deallocate(liminf,limsup)

   endif

   close(1)

   allocate(shift(Ninterf))

   if(iref==1)then
       shift(1:Ninterf)=moy_ref(1:Ninterf,1)
   elseif(iref==2)then
       shift(1:Ninterf)=moy_ref(1:Ninterf,2)
   elseif(iref==3)then
       shift(1:Ninterf)=(moy_ref(1:Ninterf,1)+moy_ref(1:Ninterf,2))/2.
   elseif(iref==4)then
       shift(1:Ninterf)=moy_ref(1:Ninterf,3)
   elseif(iref==5)then
       shift(1:Ninterf)=moy_ref(1:Ninterf,4)
   elseif(iref==6)then
       shift(1:Ninterf)=(moy_ref(1:Ninterf,1)+      &
                         moy_ref(1:Ninterf,3)+moy_ref(1:Ninterf,4))/3.
   else
      shift(1:Ninterf)=0.
   endif

   deallocate(moy_ref,sigma_ref)

end subroutine refer_interferos

!************************************************************************

subroutine masque_commun

   use info_communes
   use file_units

   implicit none

   integer      :: irec1,irec2,j,jlines,iwidth,irect,is,i
   integer      :: nread,jj,k,nw,i_restrict,ix1,ix2,iy1,iy2
   integer,  parameter :: nreadsimult=5
   logical,  dimension(:), allocatable :: lines_ok,cols_ok
   logical,  dimension(:,:), allocatable :: masque_ini
   integer,  dimension(:,:), allocatable :: icount

   jlines=coins(4)-coins(2)+1
   iwidth=coins(3)-coins(1)+1

   write(unit_conf_out,*)'maximum fraction of discarded interferograms'
   read(unit_conf_in,*)frac_interf

   allocate (masque_ini(iwidth,jlines),lines_ok(jlines),cols_ok(iwidth))
   allocate (icount(iwidth,jlines))

   icount(:,:)=0
   cols_ok=.false.
   lines_ok=.false.

   do j=1,jlines,nreadsimult*nsamp_mask
      nread=int((jlines-j)/nsamp_mask)
      nread=min(nread,nreadsimult)
      irec1=j+coins(2)-1
      irec2=irec1+nread*nsamp_mask

      if(imed==0)then
        call ouv_lect_int_med(irec1,irec2,nsamp_mask,1,Ninterf/2)
      else
        call ouv_lect_int(irec1,irec2,nsamp_mask,1,Ninterf/2)
      endif

!$OMP PARALLEL DO PRIVATE(irect,jj,i,is)
      do irect=1,nread
        jj=j+(irect-1)*nsamp_mask
        do i=1,iwidth,nsamp_mask
          do is=1,Ninterf/2
            if(isnan(phase(is,irect,i)).or.abs(phase(is,irect,i))<0.00001_idp.or.phase(is,irect,i)<-9990.)then
              icount(i,jj)=icount(i,jj)+1
            endif
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      deallocate(phase,coh)

      if(imed==0)then
        call ouv_lect_int_med(irec1,irec2,nsamp_mask,Ninterf/2+1,Ninterf)
      else
        call ouv_lect_int(irec1,irec2,nsamp_mask,Ninterf/2+1,Ninterf)
      endif

!$OMP PARALLEL DO PRIVATE(irect,jj,i,is)
      do irect=1,nread
        jj=j+(irect-1)*nsamp_mask
        do i=1,iwidth,nsamp_mask
          do is=1,Ninterf-Ninterf/2
            if(isnan(phase(is,irect,i)).or.abs(phase(is,irect,i))<0.00001_idp.or.phase(is,irect,i)<-9990.)then
              icount(i,jj)=icount(i,jj)+1
            endif
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

! Here threshold on the minimum number of interferograms
!$OMP PARALLEL DO PRIVATE(irect,jj,i)
      do irect=1,nread
        jj=j+(irect-1)*nsamp_mask
        do i=1,iwidth,nsamp_mask
          if(icount(i,jj)>(int(Ninterf*frac_interf)))then
            masque_ini(i,jj)=.false.
          else
            masque_ini(i,jj)=.true.
            lines_ok(jj)=.true.
            cols_ok(i)=.true.
          endif
        enddo
      enddo
!$OMP END PARALLEL DO

      deallocate(phase,coh)

   enddo
   deallocate(icount)

   i=(int((iwidth-1)/nsamp_mask))*nsamp_mask+1
   do while(.not.cols_ok(i))
     i=i-nsamp_mask
   enddo
   coins(3)=coins(1)+i-1
   coinsll(3)=coinsll(1)+(i-1)*xstep
   i=1
   do while(.not.cols_ok(i))
     i=i+nsamp_mask
   enddo
   coins(1)=coins(1)+i-1
   coinsll(1)=coinsll(1)+(i-1)*xstep
   j=(int((jlines-1)/nsamp_mask))*nsamp_mask+1
   do while(.not.lines_ok(j))
     j=j-nsamp_mask
   enddo
   coins(4)=coins(2)+j-1
   coinsll(4)=coinsll(2)+(j-1)*ystep
   j=1
   do while(.not.lines_ok(j))
     j=j+nsamp_mask
   enddo
   coins(2)=coins(2)+j-1
   coinsll(2)=coinsll(2)+(j-1)*ystep
   write(*,*) 'portion utile des interferos',(coins(k),k=1,4)
   write(*,*) 'portion utile des interferos',(coinsll(k),k=1,4)

! Here we extend the computation to the whole domain to keep the same geometry
! comment if you want to restrain to the rectangle including the non-masked area
   coins(1)=1
   coins(3)=(int((iwidth-1)/nsamp_mask))*nsamp_mask+1
   coins(2)=1
   coins(4)=(int((jlines-1)/nsamp_mask))*nsamp_mask+1

   i_restrict=0
   write(unit_conf_out,*)'Would you like to restrict the area of inversion ?(y=1,n=0)'
   read(unit_conf_in,*)i_restrict
   write(unit_conf_out,*)'Give four corners, lower, left, top, right in file pixel coord'
   read(unit_conf_in,*)ix1,iy1,ix2,iy2
   if(i_restrict.eq.1)then
     coins(1)=max(coins(1),ix1)
     coins(2)=max(coins(2),iy1)
     coins(3)=min(coins(3),ix2)
     coins(4)=min(coins(4),iy2)
   endif

   deallocate(lines_ok,cols_ok)

   jlines=coins(4)-coins(2)+1
   iwidth=coins(3)-coins(1)+1

   allocate (masque(iwidth,jlines))
   masque(:,:)=.false.
   masque(:,:)=masque_ini(coins(1):coins(1)+iwidth-1,coins(2):coins(2)+jlines-1)
!   masque(:,:)=masque_ini(1:iwidth,1:jlines)

   deallocate(masque_ini)

   open(1,file='masque_commun',status='unknown')
   nw=int(iwidth/(126.*nsamp_mask))+1
   do j=1,jlines,nw*nsamp_mask
      write(1,'(126L1)')(masque(i,j),i=1,iwidth,nw*nsamp_mask)
   enddo
   close(1)


end subroutine masque_commun

!************************************************************************

subroutine ouv_lect_int(irec1,irec2,nsamp,nini,nfin) 

   use info_communes

   implicit none

   character(len=8) :: date1,date2
   integer      :: is,ios,filen,irec,irec1,irec2,irect,iwidth
   integer      :: i,iphase_rej,j,nsamp,nini,nfin,Ninterfouv
   integer      :: irec2_cut,irect_cut,njump
   character(len=255) :: nom_interfero,noml
   logical,   dimension(:), allocatable :: flag_Im1_Im2
   real,   dimension(:), allocatable :: phase_rej,coh_rej

   irect=int((irec2-irec1)/nsamp)+1
   iwidth=coins(3)-coins(1)+1
   Ninterfouv=nfin-nini+1

   allocate(flag_Im1_Im2(Ninterf),phase(Ninterfouv,irect,iwidth),coh(Ninterfouv,irect,iwidth))

   do is=nini,nfin

      irec2_cut=min(irec2,Taille(is,2))
      irect_cut=int((irec2_cut-irec1)/nsamp)+1
      phase(is-nini+1,1:irect,1:iwidth)=0.
      filen = Taille(is,1)*4
      date1 = Im1_Im2(is,1)
      date2 = Im1_Im2(is,2)
!      if(irep_comb==0.and.flag_Im1_Im2_comb(is))then
      if(itype_ifg==0)then
      noml='//'//date1//'-'//date2//'_pre_inv.unw'
      njump=2
      else
      noml='//'//date1//'-'//date2//'.r4'
      njump=1
      endif
!      if(irep1==2)then
      nom_interfero='LN_DATA'//noml
!      else
!      nom_interfero='LN_DATA'//noml
!      endif
!      write(*,'(a)')nom_interfero

   open(7,iostat=IOS,file=nom_interfero,status='old',form='unformatted',    &
        access='direct',recl=filen,action='read')
   if (IOS==0) then
      flag_Im1_Im2(is)=.true.
      iphase_rej=coins(1)-Dec(is,1)-1
      if(irect_cut.ge.1.and.irec1.le.Taille(is,2))then
      do i=1,irect_cut
        irec=(irec1-Dec(is,2)+(i-1)*nsamp)*njump
        if(iphase_rej>0)then
          allocate(phase_rej(iphase_rej))
          read(7,rec=irec)(phase_rej(j),j=1,iphase_rej),(phase(is-nini+1,i,j),j=1,iwidth)
          deallocate(phase_rej)
        else
          read(7,rec=irec)(phase(is-nini+1,i,j),j=1,iwidth)
        endif
      enddo
      endif
      close(7)
   else
      flag_Im1_Im2(is)=.false.
   endif
   enddo

!   if(irep2>0)then
!
!   do is=nini,nfin
!
!    if(.not.flag_Im1_Im2(is))then
!
!      irec2_cut=min(irec2,Taille(is,2))
!      irect_cut=int((irec2_cut-irec1)/nsamp)+1
!      phase(is-nini+1,:,:)=0.
!      date1 = Im1_Im2(is,1)
!      date2 = Im1_Im2(is,2)
!      filen = Taille(is,1)*4
!
!      if(irep_comb==0.and.flag_Im1_Im2_comb(is))then
!      noml='//'//date1//'-'//date2//'_pre_inv.unw'
!      else
!      noml='//'//date1//'-'//date2//'_pre_inv.unw'
!      endif
!      if(irep2==2)then
!      nom_interfero='LN_DATA'//noml
!      else
!      nom_interfero='LN_DATA'//noml
!      endif
!      write(*,'(a)')nom_interfero
!
!   open(7,iostat=IOS,file=nom_interfero,status='old',form='unformatted',    &
!        access='direct',recl=filen,action='read')
!   if (IOS==0) then
!      flag_Im1_Im2(is)=.true.
!      if(irect_cut.ge.1.and.irec1.le.Taille(is,2))then
!      do i=1,irect_cut
!        irec=(irec1-Dec(is,2)+(i-1)*nsamp)*2
!        iphase_rej=coins(1)-Dec(is,1)-1
!        if(iphase_rej>0)then
!          allocate(phase_rej(iphase_rej))
!          read(7,rec=irec)(phase_rej(j),j=1,iphase_rej),(phase(is-nini+1,i,j),j=1,iwidth)
!          deallocate(phase_rej)
!        else
!          read(7,rec=irec)(phase(is-nini+1,i,j),j=1,iwidth)
!        endif
!      enddo
!      endif
!      close(7)
!   else
!     write(*,*)nom_interfero,' non trouve' 
!   endif
!
!   endif

!   enddo

!   endif

   if(iwcoh==0)then
   do is=nini,nfin

      irec2_cut=min(irec2,Taille(is,2))
      irect_cut=int((irec2_cut-irec1)/nsamp)+1
      coh(is-nini+1,1:irect,1:iwidth)=0.
      filen = Taille(is,1)*4
      date1 = Im1_Im2(is,1)
      date2 = Im1_Im2(is,2)
      if(itype_coh==0)then
        noml='//'//date1//'-'//date2//'.cor'
        njump=2
      else
        noml='//'//date1//'-'//date2//'-CC.r4'
        njump=1
      endif
      nom_interfero='LN_DATA'//noml
!      write(*,'(a)')nom_interfero

      open(7,iostat=IOS,file=nom_interfero,status='old',form='unformatted',    &
        access='direct',recl=filen,action='read')
      iphase_rej=coins(1)-Dec(is,1)-1
      if(irect_cut.ge.1.and.irec1.le.Taille(is,2))then
      do i=1,irect_cut
        irec=(irec1-Dec(is,2)+(i-1)*nsamp)*njump
        if(iphase_rej>0)then
          allocate(coh_rej(iphase_rej))
          read(7,rec=irec)(coh_rej(j),j=1,iphase_rej),(coh(is-nini+1,i,j),j=1,iwidth)
          deallocate(coh_rej)
        else
          read(7,rec=irec)(coh(is-nini+1,i,j),j=1,iwidth)
        endif
      enddo
      endif
      close(7)
   enddo
   endif

   deallocate(flag_Im1_Im2)

end subroutine ouv_lect_int

!************************************************************************

subroutine taille_dec

   use info_communes

   implicit none

   character(len=8) :: date1,date2
   integer      :: is,ios,io,i,nxr,nyr
   character(len=255) :: nom_interfero_rsc,buff,noml
   character(len=12) ::  width,filel,nbuff
   logical,   dimension(:), allocatable :: flag_Im1_Im2
   real(IDP),   dimension(:,:), allocatable :: orig

   allocate(flag_Im1_Im2(Ninterf),Taille(Ninterf,2),orig(Ninterf,2))
   do is=1,Ninterf

      date1 = Im1_Im2(is,1)
      date2 = Im1_Im2(is,2)

!      if(irep_comb==0.and.flag_Im1_Im2_comb(is))then
      if(itype_ifg==0)then
      noml='//'//date1//'-'//date2//'_pre_inv.unw.rsc'
      else
      noml='//'//date1//'-'//date2//'.r4.rsc'
      endif
!      if(irep1==2)then
      nom_interfero_rsc='LN_DATA'//noml
!      else
!      nom_interfero_rsc='LN_DATA'//noml
!      endif
!      write(*,'(a)')nom_interfero_rsc

   open(7,iostat=IOS,file=nom_interfero_rsc,status='old',form='formatted',action='read')
   if (IOS==0) then
      flag_Im1_Im2(is)=.true.
      io=0
      do i=1,200
         read(7,*)buff,nbuff
         read(buff,'(a5)')width
         if(width.eq.'WIDTH')then
           io=io+1
           read(nbuff,*)nxr
         endif
         read(buff,'(a11)')filel
         if(filel.eq.'FILE_LENGTH')then
           io=io+1
           read(nbuff,*)nyr
         endif
         if(io.eq.2) goto 10
       enddo
       flag_Im1_Im2(is)=.false.
10     continue
!       write(*,*)date1,'_',date2,' width ',nxr,'filelength ',nyr
!       write(*,*)date1,'_',date2,' xfirst ',xfirst,' yfirst ',yfirst
       Taille(is,1)=nxr
       Taille(is,2)=nyr
       orig(is,1)=0.
       orig(is,2)=0.
      close(7)
   else
      flag_Im1_Im2(is)=.false.
      write(*,*)'interfero_rsc ',date1,'_',date2,' non trouve'
   endif
   enddo

!   if(irep2>0)then
!   do is=1,Ninterf
!
!    if(.not.flag_Im1_Im2(is))then
!
!      date1 = Im1_Im2(is,1)
!      date2 = Im1_Im2(is,2)
!
!      if(irep_comb==0.and.flag_Im1_Im2_comb(is))then
!      noml='//'//date1//'-'//date2//'_pre_inv.unw.rsc'
!      else
!      noml='//'//date1//'-'//date2//'_pre_inv.unw.rsc'
!      endif
!      if(irep2==2)then
!      nom_interfero_rsc='LN_DATA'//noml
!      else
!      nom_interfero_rsc='LN_DATA'//noml
!      endif
!      write(*,'(a)')nom_interfero_rsc
!
!   open(7,iostat=IOS,file=nom_interfero_rsc,status='old',form='formatted',action='read')
!   if (IOS==0) then
!      flag_Im1_Im2(is)=.true.
!!      write(*,*)'rsc interfero ',date1,'_',date2,' trouve'
!      io=0
!      do i=1,50
!         read(7,*)buff,nbuff
!         read(buff,'(a5)')width
!         if(width.eq.'WIDTH')then
!           io=io+1
!           read(nbuff,*)nxr
!         endif
!         read(buff,'(a11)')filel
!         if(filel.eq.'FILE_LENGTH')then
!           io=io+1
!           read(nbuff,*)nyr
!         endif
!         if(io.eq.2) goto 11
!       enddo
!       flag_Im1_Im2(is)=.false.
!       write(*,*)date1,'_',date2,' fichier rsc non complet'
!11     continue
!!       write(*,*)date1,'_',date2,' width ',nxr,'filelength ',nyr
!!      write(*,*)date1,'_',date2,' xfirst ',xfirst,' yfirst ',yfirst
!       Taille(is,1)=nxr
!       Taille(is,2)=nyr
!       orig(is,1)=0.
!       orig(is,2)=0.
!      close(7)
!   else
!      write(*,*)'rsc interfero ',date1,'_',date2,' non trouve'
!      write(*,*)'ou fichier rsc non complet'
!   endif

!   endif

!   enddo
!   endif

   allocate(Dec(Ninterf,2))

   coins(1)=1
   coins(2)=1
   coins(3)=Taille(1,1)
   coins(4)=Taille(1,2)
   Dec(1,1)=0
   Dec(1,2)=0
   do is=2,Ninterf
!     write(*,*) nint((orig(is,1)-orig(1,1))/xstep),nint((orig(is,2)-orig(1,2))/ystep)
     Dec(is,1)=nint((orig(is,1)-orig(1,1))/xstep)
     Dec(is,2)=nint((orig(is,2)-orig(1,2))/ystep)
     coins(1)=max(coins(1),Dec(is,1)+1)
     coins(2)=max(coins(2),Dec(is,2)+1)
     coins(3)=min(coins(3),Dec(is,1)+Taille(is,1))
     coins(4)=max(coins(4),Dec(is,2)+Taille(is,2))
   enddo
   coinsll(1)=orig(1,1)+(coins(1)-1)*xstep
   coinsll(2)=orig(1,2)+(coins(2)-1)*ystep
   coinsll(3)=orig(1,1)+(coins(3)-1)*xstep
   coinsll(4)=orig(1,2)+(coins(4)-1)*ystep
   write(*,*) 'portion commune des interferos',(coins(i),i=1,4)
   write(*,*) 'portion commune des interferos',(coinsll(i),i=1,4)

   deallocate(flag_Im1_Im2,orig)
   open(1,file='taille_dec_retenu',status='unknown')
   do is=1,Ninterf
     write(1,*)is,Im1_Im2(is,1),' ',Im1_Im2(is,2),Taille(is,1),Taille(is,2),Dec(is,1),Dec(is,2)
   enddo
   close(1)

end subroutine taille_dec

!*************************************************************************

subroutine listing

   use info_communes
   use file_units

   implicit none

   integer,   parameter :: nsmax = 300
   integer,   parameter :: nnmax = 3000

   character(len=8) :: date1,date2
   character(len=255) :: nom_interfero,noml
   integer      :: is,Simage_ini,Ninterf_ini,ios,i,niter
   integer      :: k1,k2,iss
   real,   dimension(:), allocatable :: Im_s_ini
   character(len=8),   dimension(:), allocatable :: Im_ini
   character(len=8),   dimension(:,:), allocatable :: Im1_Im2_ini
   logical,   dimension(:), allocatable :: flag_Im1_Im2_ini,flag_Im1_Im2_comb_ini
   integer,   dimension(:), allocatable :: Im_num_ini,Im_num
   real,   dimension(:), allocatable :: basel,niveau_eau_ini,base_im_ini,baselr
   real :: basemax,base_moy


   write(unit_conf_out,*)'outliers elimination by the median (only if nsamp>1) ? (y=0,n=1)'
   read(unit_conf_in,*)imed

   allocate(Im_ini(nsmax),Im1_Im2_ini(nnmax,2),Im_s_ini(nsmax),basel(nnmax),base_im_ini(nsmax))

   write(unit_conf_out,*)'name of the image list'
   read(unit_conf_in,'(a)') liste_image

   write(unit_conf_out,*)'sort by date (0) or by another variable (1) ?'
   read(unit_conf_in,*)itri

   allocate(niveau_eau_ini(nsmax))
   open(7,iostat=IOS,file=liste_image,status='old',form='formatted')

   if (IOS==0) then

      do is=1,nsmax
         read(7,*,end=1)Im_ini(is),Im_s_ini(is),niveau_eau_ini(is),base_im_ini(is)
      enddo
1     continue
      Simage_ini=is-1
      close(7)

   else

      write(*,*)'pas de liste d images'
      stop 255

   endif

   if(itri==1)then
      write(*,*)'ATTENTION: changer aussi liste des baselines'
      call dsort(niveau_eau_ini,Im_ini,Simage_ini,2)
      do is=1,Simage_ini
         Im_s_ini(is)=int(niveau_eau_ini(is)*100)
!         write(*,*)Im_ini(is),Im_s_ini(is)
      enddo
   endif

   write(unit_conf_out,*)'enter the name of the list of interferograms'
   read(unit_conf_in,'(a)') liste_interfero

   write(unit_conf_out,*)'interferogram format (RMG : 0; R4 :1) (date1-date2_pre_inv.unw or date1-date2.r4)'
   read(unit_conf_in,*)itype_ifg

   open(7,iostat=IOS,file=liste_interfero,status='old',form='formatted')

   if (IOS==0) then

      do is=1,nnmax
         read(7,*,end=2)Im1_Im2_ini(is,1),Im1_Im2_ini(is,2)
!         read(7,*,end=2)Im1_Im2_ini(is,1),Im1_Im2_ini(is,2),basel(is)
         do iss=1,Simage_ini
           if(Im1_Im2_ini(is,1).eq.Im_ini(iss))k1=iss
           if(Im1_Im2_ini(is,2).eq.Im_ini(iss))k2=iss
         enddo
         basel(is)=base_im_ini(k2)-base_im_ini(k1)
      enddo
2     continue
      Ninterf_ini=is-1
      close(7)
      print*,'Number of interferograms ',Ninterf_ini

   else

      write(*,*)'pas de liste d interferos'
      stop 255

   endif

   write(unit_conf_out,*)'maximum perpendicular baseline ?'
   read(unit_conf_in,*)basemax

!   write(unit_conf_out,*)'utilisation des interferos corriges (oui:0, non:1) ?'
!   read(unit_conf_in,*)irep_comb
   irep_comb=0

   allocate(flag_Im1_Im2_ini(Ninterf_ini),flag_Im1_Im2_comb_ini(Ninterf_ini))

!   write(unit_conf_out,*)'localisation des interferos (data:1, data5 :2)'
!   read(unit_conf_in,*) irep1
   irep1=1

   do is=1,Ninterf_ini

      date1 = Im1_Im2_ini(is,1)
      date2 = Im1_Im2_ini(is,2)

!      if(irep_comb==0)then
      if(itype_ifg==0)then
      noml='//'//date1//'-'//date2//'_pre_inv.unw'
      else
      noml='//'//date1//'-'//date2//'.r4'
      endif
!      if(irep1==2)then
      nom_interfero='LN_DATA'//noml
!      else
!      nom_interfero='LN_DATA'//noml
!      endif
!      write(*,'(a)')nom_interfero

   open(7,iostat=IOS,file=nom_interfero,status='old',form='unformatted',action='read')
   if (IOS==0) then
      flag_Im1_Im2_ini(is)=.true.
      flag_Im1_Im2_comb_ini(is)=.true.
!      write(*,*)'interfero ',date1,'_',date2,' trouve'
      close(7)
   else
      flag_Im1_Im2_ini(is)=.false.
      flag_Im1_Im2_comb_ini(is)=.false.
      write(*,*)'interfero ',date1,'_',date2,' non trouve'
      print*,nom_interfero
!      noml='//'//date1//'-'//date2//'_pre_inv.unw'
!      if(irep1==2.and.irep_comb==0)then
!      nom_interfero='LN_DATA'//noml
!      elseif(irep1==1.and.irep_comb==0)then
!      nom_interfero='LN_DATA'//noml
!      endif
!      if(irep_comb==0)then
!           open(7,iostat=IOS2,file=nom_interfero,status='old',form='unformatted',action='read')
!           if (IOS2==0) then
!              flag_Im1_Im2_ini(is)=.true.
!              write(*,*)'interfero ',date1,'_',date2,' trouve'
!              close(7)
!              flag_Im1_Im2_comb_ini(is)=.false.
!           endif
!      endif
   endif
   enddo

!   write(unit_conf_out,*)'localisation 2 des interferos non trouves(data:1, data5:2, sinon 0)'
!   read(unit_conf_in,*) irep2
   irep2=0

!   if(irep2>0)then

!   do is=1,Ninterf_ini
!
!    if(.not.flag_Im1_Im2_ini(is))then
!
!      date1 = Im1_Im2_ini(is,1)
!      date2 = Im1_Im2_ini(is,2)
!
!      if(irep_comb==0)then
!      noml='//'//date1//'-'//date2//'_pre_inv.unw'
!      else
!      noml='//'//date1//'-'//date2//'_pre_inv.unw'
!      endif
!      if(irep2==2)then
!      nom_interfero='LN_DATA'//noml
!      else
!      nom_interfero='LN_DATA'//noml
!      endif
!      write(*,'(a)')nom_interfero
!
!   open(7,iostat=IOS,file=nom_interfero,status='old',form='unformatted',action='read')
!   if (IOS==0) then
!      flag_Im1_Im2_ini(is)=.true.
!      flag_Im1_Im2_comb_ini(is)=.true.
!      write(*,*)'interfero ',date1,'_',date2,' trouve'
!      close(7)
!   else
!      write(*,*)'interfero ',date1,'_',date2,' non trouve'
!      noml='//'//date1//'-'//date2//'_pre_inv.unw'
!      if(irep2==2.and.irep_comb==0)then
!      nom_interfero='LN_DATA'//noml
!      elseif(irep2==1.and.irep_comb==0)then
!      nom_interfero='LN_DATA'//noml
!      endif
!      if(irep_comb==0)then
!           open(7,iostat=IOS2,file=nom_interfero,status='old',form='unformatted',action='read')
!           if (IOS2==0) then
!              flag_Im1_Im2_ini(is)=.true.
!!              write(*,*)'interfero ',date1,'_',date2,' trouve'
!              close(7)
!              flag_Im1_Im2_comb_ini(is)=.false.
!           else
!             write(*,*)'interfero ',date1,'_',date2,' non trouve'
!           endif
!      else
!         write(*,*)'interfero ',date1,'_',date2,' non trouve'
!      endif
!   endif

!   endif
!
!   enddo
!
!   endif

   if(irep_comb==1)flag_Im1_Im2_comb_ini=.false.
   do is=1,Ninterf_ini
   if(abs(basel(is)).gt.basemax)flag_Im1_Im2_ini(is)=.false.
   enddo

   iwcoh=1
   itype_coh=0
   write(unit_conf_out,*)'Weight input interferograms by coherence or correlation maps ? (y:0,n:1)'
   read(unit_conf_in,*)iwcoh
   write(unit_conf_out,*)'coherence file format (RMG : 0; R4 :1) (date1-date2.cor or date1-date2-CC.r4)'
   read(unit_conf_in,*)itype_coh

   if(iwcoh==0)then
   do is=1,Ninterf_ini

      if(flag_Im1_Im2_ini(i))then
      date1 = Im1_Im2_ini(is,1)
      date2 = Im1_Im2_ini(is,2)

      if(itype_coh==0)then
      noml='//'//date1//'-'//date2//'.cor'
      else
      noml='//'//date1//'-'//date2//'-CC.r4'
      endif
      nom_interfero='LN_DATA'//noml

      open(7,iostat=IOS,file=nom_interfero,status='old',form='unformatted',action='read')
      if (IOS==0) then
         close(7)
      else
         flag_Im1_Im2_ini(is)=.false.
         write(*,*)'interfero ',date1,'_',date2,'non trouve'
      endif
      endif
   enddo
   endif

   write(unit_conf_out,*)'minimum number of interferograms linked with each image'
   read(unit_conf_in,*) irepIm
   if(irepIm<1)irepIm=1

   allocate(Im_num_ini(Simage_ini))

   niter=0
3  continue
   do i=1,Simage_ini 
     Im_num_ini(i)=0
   do is=1,Ninterf_ini
      if(flag_Im1_Im2_ini(is))then
      date1 = Im1_Im2_ini(is,1)
      date2 = Im1_Im2_ini(is,2)
      if(Im_ini(i).eq.date1.or.Im_ini(i).eq.date2)then
        Im_num_ini(i)=Im_num_ini(i)+1
      endif
      endif
   enddo
   write(*,*)'image ',Im_ini(i),' vue ',Im_num_ini(i),' fois'
   enddo
   if(irepIm.gt.1.and.niter.lt.10)then
   do i=1,Simage_ini
   if(Im_num_ini(i)<irepIm)then
        do is=1,Ninterf_ini
      if(flag_Im1_Im2_ini(is))then
      date1 = Im1_Im2_ini(is,1)
      date2 = Im1_Im2_ini(is,2)
      if(Im_ini(i).eq.date1.or.Im_ini(i).eq.date2)then
        flag_Im1_Im2_ini(is)=.false.
      endif
      endif 
        enddo
    endif
   enddo
   niter=niter+1
!   write(*,*)'niter ',niter
   goto 3
   endif

   Simage=COUNT(Im_num_ini>=irepIm)
   Ninterf=COUNT(flag_Im1_Im2_ini)
   write(*,*)'nb d images ',Simage,' nb d interfero ',Ninterf

   allocate(Im(Simage),Im1_Im2(Ninterf,2),Im_s(Simage),flag_Im1_Im2_comb(Ninterf))
   allocate(baselr(Ninterf),Im_num(Simage),niveau_eau(Simage),base_im(Simage))

! images retenues, nb d interferos les incluant, date en 10**2 sec
   open(1,file='images_retenues',status='unknown')
   is=0
   do i=1,Simage_ini
     if(Im_num_ini(i)>=irepIm)then
       is=is+1
       Im(is)=Im_ini(i) 
       Im_s(is)=Im_s_ini(i)
       Im_num(is)=Im_num_ini(i)
       base_im(is)=base_im_ini(i)
       niveau_eau(is)=niveau_eau_ini(i)
       write(1,*)is,' ',Im(is),' ',Im_num(is),Im_s(is),niveau_eau(is),base_im(is)
     endif 
   enddo
   close(1)

!   base_moy=SUM(base_im)/Simage
   base_moy=base_im(1)
   base_im(:)=base_im(:)-base_moy

! ponderation des interferos pour que le poids de chaque image soit id.
   write(unit_conf_out,*)'Weighting of interferograms for getting a uniform weight on each image (y=0;n=1) ?'
   read(unit_conf_in,*)ipondimage
   if(ipondimage==0)then
     call ponder_image(Im_num)
   endif

   open(1,file='interfero_retenus',status='unknown')
   is=0
   do i=1,Ninterf_ini
      if(flag_Im1_Im2_ini(i))then
       is=is+1
       Im1_Im2(is,1)=Im1_Im2_ini(i,1)
       Im1_Im2(is,2)=Im1_Im2_ini(i,2)
       baselr(is)=basel(i)
       flag_Im1_Im2_comb(is)=flag_Im1_Im2_comb_ini(i)
       write(1,*)is,' ',Im1_Im2(is,1),' ',Im1_Im2(is,2),' ',basel(i),flag_Im1_Im2_comb_ini(i)
      else
        write(*,*)'interfero ',Im1_Im2_ini(i,1),'-',Im1_Im2_ini(i,2),' non conserve',basel(i)
      endif
   enddo
   deallocate(baselr)

   deallocate(flag_Im1_Im2_ini,Im1_Im2_ini,Im_num_ini,Im_ini,Im_s_ini)
   deallocate(basel,flag_Im1_Im2_comb_ini,Im_num)
   deallocate(niveau_eau_ini)

return
end subroutine

!*************************************************************************

subroutine ponder_image(Im_num)

   use info_communes

   implicit none

   character(len=8) :: date
   integer,   dimension(Simage) :: Im_num
   real(IDP), dimension(:,:), allocatable :: matp
   real(IDP),   dimension(:), allocatable :: vp,WORK,vectemp
   integer,   dimension(:), allocatable :: IWORK
   integer      :: nrhs,rank,LWORK,LIWORK,INFO
   integer      :: nlig,ncol,is,iss
   real(IDP)    :: rcond,p1,p2

   nlig=Simage+Ninterf
   ncol=Ninterf
   allocate(matp(nlig,ncol),pondim(nlig))

   print*,'attention:  ajuster manuellement la subroutine ponder_image'
   print*,'pour avoir un poids par interfero >0 et un poids par image proche de 1'
   pondim(:)=0._idp
   matp(:,:)=0._idp
   do is=1,Simage
      pondim(is)=1._idp
   enddo
!  valeur a ajuster
   do is=1,Ninterf
      do iss=1,Simage
        date=Im(iss)
        if(date.eq.Im1_Im2(is,1))p1=Im_num(iss)
        if(date.eq.Im1_Im2(is,2))p2=Im_num(iss)
      enddo
      pondim(Simage+is)=2.*(2./(p1+p2))
!      pondim(Simage+is)=0.15_idp
   enddo
   do is=1,Simage
      date=Im(is)
      do iss=1,Ninterf
        if(date.eq.Im1_Im2(iss,1).or.date.eq.Im1_Im2(iss,2)) matp(is,iss)=1._idp
      enddo
   enddo
!  valeur a ajuster
   do is=1,Ninterf
      matp(Simage+is,is)=2.0_idp
!      matp(Simage+is,is)=1.0_idp
   enddo

!       nombre de vecteurs solution a trouver:
        nrhs=1
!       mise a zero des valeurs propres < rcond * valeur propre max
        rcond=0.000001_idp
        LIWORK=32*ncol
        allocate(IWORK(LIWORK),WORK(1),vp(nlig))
        ! compute the minimum-norm solution to a real linear least squares problem
        ! methode divide and conquer for SVD
        call SGELSD(nlig, ncol, nrhs, matp, nlig, pondim, nlig,  vp,  rcond,   &
                          rank, WORK, -1, IWORK, INFO )
        LWORK=WORK(1)
        deallocate(WORK)
        allocate(WORK(LWORK))
        call SGELSD(nlig, ncol, nrhs, matp, nlig, pondim, nlig,  vp,  rcond,   &
                          rank, WORK, LWORK, IWORK, INFO )
       if (INFO < 0) then
         write(0, *) "error: SGELSD() ",-INFO,"th argument is invalid"
         stop 255
       else if (INFO > 0) then
         write(0, *) "warning: SGELSD() did not converge"
       endif
        deallocate(IWORK,WORK,vp)

      do is=1,Ninterf
!        write(*,*) 'poids interfero ',is,pondim(is)
         if(pondim(is).lt.0.03)write(*,*) 'poids interfero ',is,pondim(is)
      enddo

      allocate(vectemp(1:Simage))
   vectemp(:)=0.
   do is=1,Simage
      date=Im(is)
      do iss=1,Ninterf
        if(date.eq.Im1_Im2(iss,1).or.date.eq.Im1_Im2(iss,2)) vectemp(is)=vectemp(is)+pondim(iss)
      enddo
   enddo
      do is=1,Simage
!        write(*,*) 'poids image ',is,vectemp(is)
         if(vectemp(is).lt.0.5.or.vectemp(is).gt.1.5) write(*,*) 'poids image ',is,vectemp(is)
      enddo
      deallocate(vectemp)

   deallocate(matp)

return
end subroutine ponder_image

!************************************************************************
subroutine ponderation_variance

   use info_communes
   use file_units

   implicit none

   integer,   parameter :: nsmax = 300
   integer,   parameter :: nnmax = 3000
   integer      :: irec1,irec2,is,irect,iwidth,i,j,k,ndata,jlines,jj,iss
   integer      :: IOS
   real :: val
   integer,   dimension(:), allocatable :: liminf,limsup
   integer :: nsamp_ref = 2
   character(len=8) :: date1,date2,datep1,datep2
   double precision,   dimension(:), allocatable :: moy_int
   double precision :: moy_sigma
   real,   dimension(:,:), allocatable :: phase1
   real,   dimension(:), allocatable :: weight
   double precision :: moy_int1,moy_int2,var_int1,var_int2
   character(len=8),   dimension(:,:), allocatable :: Im1_Im2_ini
   real :: varqualt

   nsamp_ref = 2*nsamp_mask

   write(unit_conf_out,*) 'Weighting by image noise amplitude (atmosphere) (y:0,n:1) ?'
   read(unit_conf_in,*) iqual
   write(unit_conf_out,*) 'Weighting by interferogram variance (y:0,n:1) or user given weight (2) ?'
   read(unit_conf_in,*) iponder
   write(unit_conf_out,*) 'use of interferogram covariances (y:0,n:1) ?'
   read(unit_conf_in,*) icov

   if(iqual==0)then
     allocate(varqual(Simage))
     open(1,file=liste_image,status='old')
     do iss=1,nsmax
       read(1,*,end=1) date1, val, val, val, varqualt
       do is=1,Simage
         if(date1.eq.Im(is))varqual(is)=varqualt
       enddo
     enddo
1    continue
     close(1)
     print*,(varqual(is),is=1,Simage)
     varqual(:)=1./varqual(:)
     moy_sigma=SUM(varqual(:))/Simage
     varqual(:)=varqual(:)/moy_sigma
   endif

!! A Tester
   if(iponder==2)then
    iponder=1
    ipondimage=0
    print*,'List of weights read in liste_interfero'
    allocate(pondim(Ninterf),weight(nnmax),Im1_Im2_ini(nnmax,2))
    pondim(:)=1.d0
    open(7,iostat=IOS,file=liste_interfero,status='old',form='formatted')
    if (IOS==0) then
      do is=1,nnmax
         read(7,*,end=2)Im1_Im2_ini(is,1),Im1_Im2_ini(is,2),weight(is)
         do iss=1,Ninterf
           if((Im1_Im2_ini(is,1).eq.Im1_Im2(iss,1)).and.(Im1_Im2_ini(is,2).eq.Im1_Im2(iss,2)))then
             pondim(iss)=weight(is)
           endif
         enddo
      enddo
2     continue
      close(7)
      deallocate(Im1_Im2_ini,weight)
    else
     print*,'liste_interfero does not include weights for inversion'
     stop
    endif
   endif

   if(iponder==0)then

   irec1=coins(2)+int(33/nsamp_mask)*nsamp_mask
   irec2=nsamp_mask*int((coins(4)-coins(2)-33)/nsamp_mask)+coins(2)
   irect=int((irec2-irec1)/nsamp_ref)+1
   iwidth=coins(3)-coins(1)+1
   jlines=irec2-irec1+1

   allocate(liminf(jlines+nsamp_mask),limsup(jlines+nsamp_mask))
   allocate(moy_int(Ninterf),sigma_int(Ninterf,Ninterf))

! limites de la zone communes a l ensemble des interferos

   do j=1,jlines,nsamp_mask
   jj=j+irec1-coins(2)
   i=1
   do while(.not.masque(i,jj).and.i.lt.iwidth)
     i=i+nsamp_mask
   enddo
   liminf(j)=i+int(33/nsamp_mask)*nsamp_mask
   if(nsamp_mask>1)then
   do i=1,nsamp_mask-1
     liminf(j+i)=liminf(j)
   enddo
   endif
   enddo

   do j=1,jlines,nsamp_mask
   jj=j+irec1-coins(2)
   i=iwidth
   do while(.not.masque(i,jj))
     i=i-nsamp_mask
   enddo
   limsup(j)=i-int(33/nsamp_mask)*nsamp_mask
   if(nsamp_mask>1)then
   do i=1,nsamp_mask-1
     limsup(j+i)=limsup(j)
   enddo
   endif
   enddo

   moy_int(:)=0.d0
   sigma_int(:,:)=0.d0

   allocate(phase1(irect,iwidth))

!  matrice de covariance : symetrique, elements non nuls si une image commune

   do is=1,Ninterf

! termes diagonaux
      if(imed==0)then
        call ouv_lect_int_med(irec1,irec2,nsamp_ref,is,is)
      else
        call ouv_lect_int(irec1,irec2,nsamp_ref,is,is)
      endif
      
      phase1(:,:)=phase(1,:,:)

      ndata=0
      do j=1,irect
        k=(j-1)*nsamp_ref+1
        jj=k+irec1-coins(2)
        do i=liminf(k),limsup(k),nsamp_ref
          if(masque(i,jj).and.abs(phase(1,j,i))>0.00001.and.phase(1,j,i)>-9990.)then
            moy_int(is)=moy_int(is)+phase(1,j,i)
            sigma_int(is,is)=sigma_int(is,is)+phase(1,j,i)**2
            ndata=ndata+1
          endif
        enddo
      enddo
      sigma_int(is,is)=sigma_int(is,is)-moy_int(is)**2/ndata
      sigma_int(is,is)=sigma_int(is,is)/ndata
      moy_int(is)=moy_int(is)/ndata

      deallocate(phase,coh)

! termes non diagonaux
      if(icov==0.and.is<Ninterf) then
      date1 = Im1_Im2(is,1)
      date2 = Im1_Im2(is,2)

        do iss=is+1,Ninterf
        
          datep1 = Im1_Im2(iss,1) 
          datep2 = Im1_Im2(iss,2)

          if(date1.eq.datep1.or.date1.eq.datep2.or.date2.eq.datep1.or.date2.eq.datep2)then

            if(imed==0)then
               call ouv_lect_int_med(irec1,irec2,nsamp_ref,iss,iss)
            else
               call ouv_lect_int(irec1,irec2,nsamp_ref,iss,iss)
            endif

            ndata=0
            moy_int1=0.d0
            moy_int2=0.d0
            var_int1=0.d0
            var_int2=0.d0
            do j=1,irect
              k=(j-1)*nsamp_ref+1
              jj=k+irec1-coins(2)
              do i=liminf(k),limsup(k),nsamp_ref
                if(masque(i,jj).and.abs(phase(1,j,i))>0.00001.and.abs(phase1(j,i))>0.00001.and.phase(1,j,i)>-9990.and.phase1(j,i)>-9990.)then
                  moy_int1=moy_int1+phase1(j,i)
                  moy_int2=moy_int2+phase(1,j,i)
                  var_int1=var_int1+phase1(j,i)**2
                  var_int2=var_int2+phase(1,j,i)**2
                  sigma_int(is,iss)=sigma_int(is,iss)+phase(1,j,i)*phase1(j,i)
                  ndata=ndata+1
                endif
              enddo
            enddo
            sigma_int(is,iss)=sigma_int(is,iss)-moy_int1*moy_int2/ndata
            sigma_int(is,iss)=sigma_int(is,iss)/ndata
            var_int1=var_int1-moy_int1**2/ndata
            var_int2=var_int2-moy_int2**2/ndata
            var_int1=sqrt(var_int1/ndata)
            var_int2=sqrt(var_int2/ndata)
            sigma_int(is,iss)=sigma_int(is,iss)/(var_int1*var_int2)

            deallocate(phase,coh)

          endif
           
        enddo
        endif

   enddo

   deallocate(liminf,limsup,phase1)

! renormalisation des termes non diagonaux et symetrie
   if(icov==0)then
   do is=1,Ninterf-1
     do iss=is+1,Ninterf
       sigma_int(is,iss)=sigma_int(is,iss)*sqrt(sigma_int(is,is))*sqrt(sigma_int(iss,iss))
       sigma_int(iss,is)=sigma_int(is,iss)
     enddo
   enddo
!   do is=1,Ninterf
!      sigma_int(is,is)=sigma_int(is,is)
!   enddo
   endif

   open(1,file='interf_norm',status='unknown')
   do is=1,Ninterf
      write(1,'(150f5.2)')(sigma_int(is,iss),iss=1,Ninterf)
   enddo
   close(1)

   if(icov==1)then
! inverse et racine carre de la matrice diagonale D**-1/2, normalisee
   moy_sigma=0
   do is=1,Ninterf
       sigma_int(is,is)=1.d0/sqrt(sigma_int(is,is))
       moy_sigma=moy_sigma+ sigma_int(is,is)
   enddo
   moy_sigma=moy_sigma/Ninterf
   sigma_int=sigma_int/moy_sigma
   print*,'moyenne des ecarts types ',1./moy_sigma
   endif

   deallocate(moy_int)

   endif

end subroutine ponderation_variance

!************************************************************************

subroutine ouv_lect_int_med(irec1,irec2,nsamp,nini,nfin)

   use info_communes

   implicit none

! valeur seuil pour eliminer les outliers (en rad LOS)
   real(IDP),   parameter :: threshold = 1.0
   character(len=8) :: date1,date2
   integer      :: is,ios,filen,irec,irec1,irec2,irect,iwidth
   integer      :: i,iphase_rej,j,nsamp,nini,nfin,Ninterfouv
   integer      :: itemp,ireci,irecf,ntemp,j0,j1,jj,ns,nsamp2
   integer      :: nsamp_med,irecl,itempl,j00,j11,njump
   character(len=255) :: nom_interfero,noml
   logical,   dimension(:), allocatable :: flag_Im1_Im2
   real,   dimension(:), allocatable :: phase_rej,coh_rej
   real,   dimension(:), allocatable :: x
   real :: val
   real,   dimension(:,:), allocatable :: phasetemp,cohtemp
   logical       :: flagm

! si nsamp est trop petit pour faire une mediane, on elimine seulement
! les outliers grace a l ecart type autour de la mediane calculee avec une
! zone plus large
   nsamp_med=nsamp
   flagm=.true.
   if(nsamp.lt.4)then
      nsamp_med=4
      flagm=.false.
   endif

   irect=int((irec2-irec1)/nsamp)+1
   iwidth=coins(3)-coins(1)+1
   Ninterfouv=nfin-nini+1
   nsamp2=(nsamp_med+1)**2

   allocate(flag_Im1_Im2(Ninterf),phase(Ninterfouv,irect,iwidth),coh(Ninterfouv,irect,iwidth))
   allocate(phasetemp(nsamp_med+1,iwidth),x(nsamp2),cohtemp(nsamp_med+1,iwidth))

   do is=nini,nfin

      filen = Taille(is,1)*4
      date1 = Im1_Im2(is,1)
      date2 = Im1_Im2(is,2)
!      if(irep_comb==0.and.flag_Im1_Im2_comb(is))then
      if(itype_ifg==0)then
      noml='//'//date1//'-'//date2//'_pre_inv.unw'
      njump=2
      else
      noml='//'//date1//'-'//date2//'.r4'
      njump=1
      endif
!      if(irep1==2)then
      nom_interfero='LN_DATA'//noml
!      else
!      nom_interfero='LN_DATA'//noml
!      endif
!      write(*,'(a)')nom_interfero

   open(7,iostat=IOS,file=nom_interfero,status='old',form='unformatted',    &
        access='direct',recl=filen,action='read')
   if (IOS==0) then
      flag_Im1_Im2(is)=.true.
      do i=1,irect
        ireci=(irec1-Dec(is,2)+(i-1)*nsamp-nsamp_med/2)*njump
        ireci=max(1,ireci)
        irecf=(irec1-Dec(is,2)+(i-1)*nsamp+nsamp_med/2)*njump
        irecf=min(irecf,Taille(is,2)*2)
        irecl=(irec1-Dec(is,2)+(i-1)*nsamp)*njump
        iphase_rej=coins(1)-Dec(is,1)-1
        phasetemp(:,:)=0.
        if(ireci.le.Taille(is,2)*njump)then
        if(iphase_rej>0)then
          allocate(phase_rej(iphase_rej))
          itemp=0
          do irec=ireci,irecf,njump
            itemp=itemp+1
            read(7,rec=irec)(phase_rej(j),j=1,iphase_rej),(phasetemp(itemp,j),j=1,iwidth)
            if(irec==irecl)itempl=itemp
          enddo
          deallocate(phase_rej)
        else
          itemp=0
          do irec=ireci,irecf,njump
            itemp=itemp+1
            read(7,rec=irec)(phasetemp(itemp,j),j=1,iwidth)
            if(irec==irecl)itempl=itemp
          enddo
        endif
!        print*,'lecture effectuee'
        ntemp=itemp
        do j=1,iwidth,nsamp
          ns=0
          j0=max(1,j-nsamp_med/2)
          j1=min(iwidth,j+nsamp_med/2)
          j00=max(1,j-nsamp/2)
          j11=min(iwidth,j+nsamp/2)
          do itemp=1,ntemp
            do jj=j0,j1
              if(abs(phasetemp(itemp,jj)).gt.1.e-5)then
                ns=ns+1
                x(ns)=phasetemp(itemp,jj)
              endif
            enddo
          enddo
          if(ns.gt.(ntemp*(j1-j0+1))/4)then
             call median(x,ns,val)
! mediane
            if(flagm)then
              do jj=j00,j11
               phase(is-nini+1,i,jj)=val
              enddo
            else
!seult elimination outliers
              do jj=j00,j11
               if(abs(val-phasetemp(itempl,jj))<threshold)then
                 phase(is-nini+1,i,jj)=phasetemp(itempl,jj)
               else
                 phase(is-nini+1,i,jj)=0.
               endif
              enddo
            endif
          else
            do jj=j00,j11
             phase(is-nini+1,i,jj)=0.
            enddo
          endif
        enddo
!        print*,'mediane faite'
        endif
      enddo
      close(7)
   else
      flag_Im1_Im2(is)=.false.
   endif
   enddo

!   if(irep2>0)then
!
!   do is=nini,nfin
!
!    if(.not.flag_Im1_Im2(is))then
!
!      date1 = Im1_Im2(is,1)
!      date2 = Im1_Im2(is,2)
!      filen = Taille(is,1)*4
!
!      if(irep_comb==0.and.flag_Im1_Im2_comb(is))then
!      noml='//'//date1//'-'//date2//'_pre_inv.unw'
!      else
!      noml='//'//date1//'-'//date2//'_pre_inv.unw'
!      endif
!      if(irep2==2)then
!      nom_interfero='LN_DATA'//noml
!      else
!      nom_interfero='LN_DATA'//noml
!      endif
!!      write(*,'(a)')nom_interfero
!
!   open(7,iostat=IOS,file=nom_interfero,status='old',form='unformatted',    &
!        access='direct',recl=filen,action='read')
!   if (IOS==0) then
!      flag_Im1_Im2(is)=.true.
!      do i=1,irect
!        ireci=(irec1-Dec(is,2)+(i-1)*nsamp-nsamp_med/2)*2
!        ireci=max(1,ireci)
!        irecf=(irec1-Dec(is,2)+(i-1)*nsamp+nsamp_med/2)*2
!        irecf=min(irecf,Taille(is,2)*2)
!        irecl=(irec1-Dec(is,2)+(i-1)*nsamp)*2
!        iphase_rej=coins(1)-Dec(is,1)-1
!        phasetemp(:,:)=0.
!        if(ireci.le.Taille(is,2)*2)then
!        if(iphase_rej>0)then
!          allocate(phase_rej(iphase_rej))
!          itemp=0
!          do irec=ireci,irecf,2
!            itemp=itemp+1
!            read(7,rec=irec)(phase_rej(j),j=1,iphase_rej),(phasetemp(itemp,j),j=1,iwidth)
!            if(irec==irecl)itempl=itemp
!          enddo
!          deallocate(phase_rej)
!        else
!          itemp=0
!          do irec=ireci,irecf,2
!            itemp=itemp+1
!            read(7,rec=irec)(phasetemp(itemp,j),j=1,iwidth)
!            if(irec==irecl)itempl=itemp
!          enddo
!        endif
!        ntemp=itemp
!        do j=1,iwidth,nsamp
!          ns=0
!          j0=max(1,j-nsamp_med/2)
!          j1=min(iwidth,j+nsamp_med/2)
!          j00=max(1,j-nsamp/2)
!          j11=min(iwidth,j+nsamp/2)
!          do itemp=1,ntemp
!            do jj=j0,j1
!              if(abs(phasetemp(itemp,jj)).gt.1.e-5)then
!                ns=ns+1
!                x(ns)=phasetemp(itemp,jj)
!              endif
!            enddo
!          enddo
!          if(ns.gt.(ntemp*(j1-j0+1))/4)then
!             call median(x,ns,val)
!! mediane
!            if(flagm)then
!              do jj=j00,j11
!               phase(is-nini+1,i,jj)=val
!              enddo
!            else
!!seult elimination outliers
!              do jj=j00,j11
!               if(abs(val-phasetemp(itempl,jj))<threshold)then
!                 phase(is-nini+1,i,jj)=phasetemp(itempl,jj)
!               else
!                 phase(is-nini+1,i,jj)=0.
!               endif
!              enddo
!            endif
!          else
!            do jj=j00,j11
!             phase(is-nini+1,i,jj)=0.
!            enddo
!          endif
!        enddo
!        endif
!      enddo
!      close(7)
!   else
!     write(*,*)nom_interfero,' non trouve'
!   endif
!
!   endif
!
!   enddo
!
!   endif

   if(iwcoh==0)then
   do is=nini,nfin
      filen = Taille(is,1)*4
      date1 = Im1_Im2(is,1)
      date2 = Im1_Im2(is,2)
      if(itype_coh==0)then
        noml='//'//date1//'-'//date2//'.cor'
        njump=2
      else
        noml='//'//date1//'-'//date2//'-CC.r4'
        njump=1
      endif
      nom_interfero='LN_DATA'//noml

      open(7,iostat=IOS,file=nom_interfero,status='old',form='unformatted',    &
           access='direct',recl=filen,action='read')
      do i=1,irect
        ireci=(irec1-Dec(is,2)+(i-1)*nsamp-nsamp_med/2)*njump
        ireci=max(1,ireci)
        irecf=(irec1-Dec(is,2)+(i-1)*nsamp+nsamp_med/2)*njump
        irecf=min(irecf,Taille(is,2)*2)
        irecl=(irec1-Dec(is,2)+(i-1)*nsamp)*njump
        iphase_rej=coins(1)-Dec(is,1)-1
        cohtemp(:,:)=0.
        if(ireci.le.Taille(is,2)*njump)then
           if(iphase_rej>0)then
             allocate(coh_rej(iphase_rej))
             itemp=0
             do irec=ireci,irecf,njump
               itemp=itemp+1
               read(7,rec=irec)(coh_rej(j),j=1,iphase_rej),(cohtemp(itemp,j),j=1,iwidth)
               if(irec==irecl)itempl=itemp
             enddo
             deallocate(coh_rej)
           else
             itemp=0
             do irec=ireci,irecf,njump
               itemp=itemp+1
               read(7,rec=irec)(cohtemp(itemp,j),j=1,iwidth)
               if(irec==irecl)itempl=itemp
             enddo
           endif
!        print*,'lecture effectuee'
           ntemp=itemp
           do j=1,iwidth,nsamp
             ns=0
             j0=max(1,j-nsamp_med/2)
             j1=min(iwidth,j+nsamp_med/2)
             j00=max(1,j-nsamp/2)
             j11=min(iwidth,j+nsamp/2)
             do itemp=1,ntemp
               do jj=j0,j1
                   ns=ns+1
                   x(ns)=cohtemp(itemp,jj)
               enddo
             enddo
             call median(x,ns,val)
             do jj=j00,j11
                coh(is-nini+1,i,jj)=val
             enddo
           enddo
        endif
      enddo
      close(7)
   enddo
   endif

   deallocate(flag_Im1_Im2,phasetemp,cohtemp,x)

end subroutine ouv_lect_int_med

!************************************************************************


subroutine write_ts_envihdr(path, xsize, ysize, cnt, dates)
  use info_communes
  character(*),intent(in)::path
  integer,intent(in)::xsize,ysize,cnt
  character(8),dimension(:),intent(in)::dates

  character(1024)::hdrpath
  integer::i

  hdrpath=trim(path)//'.hdr'
  open(unit=50,file=hdrpath)
  write(50,'(a)')'ENVI'
  write(50,'(a,I6)')'samples = ',xsize
  write(50,'(a,I6)')'lines = ',ysize
  write(50,'(a,I3)')'bands = ',cnt
  write(50,'(a)')'header offset = 0'
  write(50,'(a)')'data type = 4'
  write(50,'(a)')'interleave = bip'
  write(50,'(a)')'byte order = 0'
  write(50,'(a,f9.2)')'data ignore value = ', nodata_output_value
  write(50,'(a)',advance='no')'band names = {'
  do i=1,cnt-1
    write(50,'(a,a)',advance='no'),dates(i),','
  enddo
  write(50,'(a,a)',advance='no')dates(cnt),'}'
  close(50)
end subroutine write_ts_envihdr

subroutine write_rmspixel_envihdr(path, xsize, ysize)
  character(*),intent(in)::path
  integer,intent(in)::xsize,ysize

  character(1024)::hdrpath

  hdrpath=trim(path)//'.hdr'
  open(unit=50,file=hdrpath)
  write(50,'(a)')'ENVI'
  write(50,'(a,I6)')'samples = ',xsize
  write(50,'(a,I6)')'lines = ',ysize
  write(50,'(a)')'bands = 1'
  write(50,'(a)')'header offset = 0'
  write(50,'(a)')'data type = 4'
  write(50,'(a)')'interleave = bip'
  write(50,'(a)')'byte order = 0'
  close(50)
end subroutine write_rmspixel_envihdr
