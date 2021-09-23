c-----------------------------------------------------------------------
      subroutine cp_setup

      ! preprossing on the tensor and call compute cp

      include 'SIZE'
      include 'MOR'

      common /scrtens_norm/ norm_c,norm_c0

      logical ifexist
      real wk(lb+1),cl0(lcglo)
      real norm_c,norm_c0

      jfield=ifield
      if (ifrom(1)) then
         ifield=1
         inquire (file='ops/u0',exist=ifexist)
         if (ifexist) call read_serial(u,nb+1,'ops/u0 ',wk,nid)

         if (ifquad) then
            call force_skew(cul,ic1,ic2,jc1,jc2,kc1,kc2)
         endif
         local_size = (ic2-ic1+1)*(jc2-jc1+1)*(kc2-kc1+1)
         norm_c = vlsc2(cul,cul,local_size)
         if (nid.eq.0) write(6,*)'norm_tens:',norm_c

         if (ifcore) then
            call set_c_slice(cuj0,cu0k,cul,ic1,ic2,jc1,jc2,kc1,kc2,nb)
            call set_c_core(cl0,cul,ic1,ic2,jc1,jc2,kc1,kc2,nb)
            local_size = (ic2-ic1+1)*(jc2-jc1)*(kc2-kc1)
            norm_c0 = vlsc2(cl0,cl0,local_size)
         if (nid.eq.0) write(6,*)'norm_tens:',norm_c0
         endif

         call compute_cp(cua,cub,cuc,cp_uw,cul,cl0,
     $                   cuj0,cu0k,'ops/cu ',u)
      endif
      if (ifrom(2)) then
         ifield=2
         inquire (file='ops/t0',exist=ifexist)
         if (ifexist) call read_serial(ut,nb+1,'ops/t0 ',wk,nid)
         call compute_cp(cta,ctb,ctc,cp_tw,ctl,cl0,
     $                   ctj0,ct0k,'ops/ct ',ut)
      endif
c     call read_cp_weight
c     call read_cp_mode
      ifield=jfield

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_cp(cp_a,cp_b,cp_c,cp_w,cl,cl0,cj0,c0k,fname,uu)

c     Compute the cp factors with given tensor cl,
c     cj0 and c0k represent the first frontal slice
c     and first lateral slice.

c     This subroutine requires cl, cj0, c0k, fname and uu and
c     returns cp_a, cp_b, cp_c ,cp_w where
c     cp_a, cp_b, cp_c, cp_w are the results of CPD (CP Decomposition),
c     cp_a, cp_b and cp_c are the factor matrices and cp_w are
c     the weights.

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrtens_norm/ norm_c,norm_c0

      common /scralstwrk/ tmp1(ltr**2),tmp2(ltr*(lb+1)),tmp3((lb+1)**2),
     $                    tmp4(ltr),tmp5(ltr),tmp6(ltr),tmp7((lb+1)),
     $                    tmp8((lb+1)),tmp9((lb+1))

      common /scrmttkrp/ cmr((lb+1)**2*ltr),crm((lb+1)**2*ltr)
      common /scrsetmode/ wk3(lt)

      real cp_a((nb+1)*ltr),cp_b((nb+1)*ltr),cp_c((nb+1)*ltr),cp_w(ltr)
      real cl(lcglo),cl0(lcglo),cj0((nb+1)**2),c0k((nb+1)**2)

      real fcm(ltr*(lb+1),3),fcmpm(ltr**2,3)
      real lsm(ltr**2,3),lsminv(ltr**2,3)
      real lsr(ltr*(lb+1))
      real uu(0:nb)
      real wk1((nb+1)*ltr),wk2(ltr)
      real cu_err,norm_c,norm_c0

      integer rank_list(2,ltr),mm
      integer glo_i,work

      character*128 fname
      character*128 fnlint

      if (nio.eq.0) write (6,*) 'inside compute_cp'

      call nekgsync
      tcp_time=dnekclock()

      call set_rank(rank_list,mm)
      do kk=1,1
c        ntr = rank_list(2,kk)
         if (ifcore) then
            if (ifquad) then
c              call rand_initial(cp_b,nb,ntr)
c              call rand_initial(cp_c,nb,ntr)
               call als_quad(cp_a,cp_b,cp_c,cp_w,cl0,
     $                       ic1,ic2,jc1+1,jc2,kc1+1,kc2,nb,ntr)
            else
               call rand_initial(cp_b,nb,ntr)
               call rand_initial(cp_c,nb,ntr)
               call als_core(cp_a,cp_b,cp_c,cp_w,cl0,
     $                       ic1,ic2,jc1+1,jc2,kc1+1,kc2,nb,ntr)
            endif
         else
            call rand_initial(cp_b,nb+1,ntr)
            call rand_initial(cp_c,nb+1,ntr)
            call als(cp_a,cp_b,cp_c,cp_w,cl,
     $               ic1,ic2,jc1,jc2,kc1,kc2,nb+1,ntr)

         endif
         call check_conv_err(cu_err,cl,cp_a,cp_b,cp_c,cp_w,uu,
     $                       ic1,ic2,jc1,jc2,kc1,kc2,
     $                       nb+1,ntr,ifcore,cj0,c0k)

         if (cu_err.le.1e-4) then
            glo_i = nid
         else 
            glo_i = -1
         endif
         call igop(glo_i,work,'M  ',1)
         if (glo_i.ge.0) then
            exit
         endif
      enddo

      gloerr = glmin(cu_err,1)
      write(6,*)'gloerr',gloerr,cu_err,nid

      if (abs((gloerr-cu_err)).le.1e-12) then 
         glo_i=nid
      else 
         glo_i=-1
      endif
      call igop(glo_i,work,'M  ',1)

      ! processors which do not have the
      ! lowest error reset factor matrices and weight
      if (nid.ne.glo_i) then 
         write(6,*)glo_i,nid,'check'
         call rzero(cp_a,(nb+1)*ltr)
         call rzero(cp_b,(nb+1)*ltr)
         call rzero(cp_c,(nb+1)*ltr)
         call rzero(cp_w,ltr)
         ntr = 0
      endif
      call nekgsync

      call gop(cp_a,wk1,'+  ',(nb+1)*ltr)
      call gop(cp_b,wk1,'+  ',(nb+1)*ltr)
      call gop(cp_c,wk1,'+  ',(nb+1)*ltr)
      call gop(cp_w,wk2,'+  ',ltr)
      call igop(ntr,work,'M  ',1)

         ! read in the cp decomposition
c        call read_cp_weight
c        call read_cp_mode

         ! debug purpose
         ! forming the tensor
c        do kk=1,ltr
c        do k=1,nb+1
c        do j=1,nb+1
c        do i=1,nb
c           cl(i+(j-1)*(nb)+(k-1)*(nb+1)*(nb)) = 
c    $      cl(i+(j-1)*(nb)+(k-1)*(nb+1)*(nb)) + 
c    $      cp_w(kk)*cua(i+(kk-1)*lub)
c    $      * cub(j+(kk-1)*(lub+1))*cuc(k+(kk-1)*(lub+1))
c        enddo
c        enddo
c        enddo
c        enddo
      if (nid.eq.0) close (unit=100)

      call nekgsync

      cp_time=cp_time+dnekclock()-tcp_time

      if (nio.eq.0) write (6,*) 'cp_time: ',cp_time
      if (nio.eq.0) write (6,*) 'lcglo=',lcglo
      if (nio.eq.0) write (6,*) 'lcloc=',lcloc

      if (nio.eq.0) write (6,*) 'exiting compute_cp'

      return
      end
c-----------------------------------------------------------------------
      subroutine set_rank(rank_list,mm)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      integer rank_list(2,ltr)
      integer mm

      mm = ceiling(real(ltr/np))
      write(6,*)mm,'length for each proc'


      call izero(rank_list,2*ltr)
      if (nid.eq.0) then
         do i=1,10
            ii = i*10
            rank_list(1,i) = mod(i,np) ! destination processor
            rank_list(2,i) = ii         ! return location
         enddo
      endif
      ky = 1
      call nekgsync

      if (nid.eq.0) then
         ni = ltr
         call fgslib_crystal_ituple_transfer(cr_h, rank_list,2,ni,
     $                                    20*ltr,ky)
      else
         ni = 0
         call fgslib_crystal_ituple_transfer(cr_h, rank_list,2,ni,
     $                                    20*ltr,ky)
      endif

c     if (nid.eq.0) then
c     do i=1,mm
c        write(6,*)i,rank_list(1,i),nid
c        write(6,*)i,rank_list(2,i),nid
c     enddo
c     endif
c     call nekgsync
c     if (nid.eq.1) then
c     do i=1,mm
c        write(6,*)i,rank_list(1,i),nid
c        write(6,*)i,rank_list(2,i),nid
c     enddo
c     endif
c     call nekgsync
c     call exitt0
   
      return
      end
c-----------------------------------------------------------------------
      subroutine als(cp_a,cp_b,cp_c,cp_w,cl,
     $               ic1,ic2,jc1,jc2,kc1,kc2,mm,nn)

      real cp_a(mm*nn),cp_b(mm*nn),cp_c(mm*nn),cp_w(nn)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)

      real fcm(mm*nn,3),fcmpm(nn*nn,3)
      real lsm(nn*nn,3),lsminv(nn*nn,3),lsr(mm*nn)
      real cmr(ic1:ic2,jc1:jc2,nn),crm(ic1:ic2,nn,kc1:kc2)
      real cp_res,cp_relres,cp_fit,cp_tol,pre_relres,reldiff
      real tmp1(nn*nn)

      integer tmp2(nn),tmp3(nn)
      integer mode,maxit,mm,nn

      real norm_c,norm_c0
      common /scrtens_norm/ norm_c,norm_c0

      if (nid.eq.0) write(6,*) 'inside als'

      maxit = 500
      cp_tol = 1e-12
      pre_relres = 1

      call copy(fcm(1,2),cp_b,mm*nn)
      call copy(fcm(1,3),cp_c,mm*nn)

      do mode=2,3
         call set_product_matrix(fcmpm(1,mode),fcm(1,mode),mm,nn)
      enddo

      do ii=1,maxit
         do mode=1,3
            call mttkrp(lsr,cl,cmr,crm,ic1,ic2,jc1,jc2,kc1,kc2,
     $                  fcm,mm,nn,mode)
            call set_lsm(lsm,fcmpm,mode,nn)
            call invmat(lsminv(1,mode),tmp1,lsm(1,mode)
     $           ,tmp2,tmp3,nn)
            do jj=1,nn
               call mxm(lsr,mm,lsminv(1+(jj-1)*nn,mode),nn,
     $         fcm(1+(mm)*(jj-1),mode),1)
            enddo
            call compute_cp_weight(cp_w,fcm(1,mode),mm,nn)
            call mode_normalize(fcm(1,mode),mm,nn)
            call set_product_matrix(fcmpm(1,mode),fcm(1,mode),mm,nn)
         enddo

         call compute_residual(cp_res,lsr,fcm(1,3),lsm(1,3),
     $                       fcmpm(1,3),cp_w,norm_c,mm,nn)
         cp_relres = cp_res/sqrt(norm_c)
         cp_fit = 1-cp_relres
         reldiff = abs(pre_relres-cp_relres)/pre_relres

         pre_relres = cp_relres
         if (nid.eq.0) write(6,1)'iter:',ii,'rel difference:',reldiff,
     $   'relative residual:',cp_relres,'residual:',cp_res,'rank:',nn
         if (cp_relres.lt.cp_tol.OR.ii.ge.maxit.OR.
     $       reldiff.le.1e-12) then
            exit
         endif
      enddo

      call copy(cp_a,fcm(1,1),mm*nn)
      call copy(cp_b,fcm(1,2),mm*nn)
      call copy(cp_c,fcm(1,3),mm*nn)

    1 format(a,i4,x,a,1p1e13.5,x,a,1p1e13.5,x,a,1p1e13.5,x,a,i4)

      if (nid.eq.0) write(6,*) 'exitting als'

      return
      end
c-----------------------------------------------------------------------
      subroutine check_conv_err(cu_err,cl,fac_a,fac_b,fac_c,cp_w,uu,
     $                          ic1,ic2,jc1,jc2,kc1,kc2,mm,nn,ifcore,
     $                          cj0,c0k)

      real cu_err
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      real fac_a(mm*nn),fac_b(mm*nn),fac_c(mm*nn),cp_w(nn)
      real uu(0:mm-1),cu(mm),cu_diff(mm),approx_cu(mm)
      real cm(ic1:ic2,jc1:jc2)
      real bcu(nn),cuu(nn),tmp(nn)
      real cj0(ic1:ic2,jc1:jc2),c0k(ic1:ic2,kc1:kc2)
      integer ic1,ic2,jc1,jc2,kc1,kc2,mm,nn,i1
      logical ifcore

      ! approximated tensor-vector multiplication
      call rzero(cu,mm)

      i1 = 0
      mb = mm
      if (ifcore) then
         mb=mm-1
         i1 = 1
      endif

      do kk=1,nn
         bcu(kk) = vlsc2(uu(i1),fac_b(1+(kk-1)*mb),mb)
         cuu(kk) = vlsc2(uu(i1),fac_c(1+(kk-1)*mb),mb)
      enddo
      call col4(tmp,bcu,cuu,cp_w,nn)
      call mxm(fac_a,mb,tmp,nn,approx_cu,1)

      if (ifcore) then
         do k=kc1,kc2
         do i=ic1,ic2
            approx_cu(i)=approx_cu(i)+c0k(i,k)*uu(0)*uu(k)
         enddo
         enddo

         do j=jc1,jc2
         do i=ic1,ic2
            approx_cu(i)=approx_cu(i)+cj0(i,j)*uu(j)*uu(0)
         enddo
         enddo

         do i=ic1,ic2
            approx_cu(i)=approx_cu(i)-cj0(i,0)*uu(0)*uu(0)
         enddo
      endif

      ! exact tensor-vector multiplication
      call rzero(cu,mm)
      do k=kc1,kc2
      do j=jc1,jc2
      do i=ic1,ic2
         cu(i)=cu(i)+cl(i,j,k)*uu(j)*uu(k)
      enddo
      enddo
      enddo

      ! difference
      call sub3(cu_diff,cu,approx_cu,mm)
      ! error
      cu_err = sqrt(vlsc2(cu_diff,cu_diff,mm))

      write(6,2) nid,'cu_err:',cu_err,'cu_relerr:',
     $           cu_err/sqrt(vlsc2(cu,cu,mm))

      do ii=1,mm
         write(6,1)ii,'difference:',cu_diff(ii),
     $             'cu:',cu(ii),'approximated cu:',approx_cu(ii)
      enddo

    2 format(i4,x,a,1p1e16.6,x,a,1p1e16.6)
    1 format(i4,x,a,1p1e16.6,x,a,1p1e16.6,x,a,1p1e16.6)

      return
      end
c-----------------------------------------------------------------------
      subroutine check_conv_err_old(cu_err,cl,fac_a,fac_b,fac_c,cp_w,uu)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real cu_err
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      real fac_a((nb+1)*ntr),fac_b((nb+1)*ntr)
      real fac_c((nb+1)*ntr),cp_w(ntr)
      real uu(0:nb)
      real cu(nb),cu_diff(nb)
      real cm(ic1:ic2,jc1:jc2)
      real bcu(ntr),cuu(ntr),tmp(ntr)
      real tmpcu(0:nb)

      ! approximated tensor-vector multiplication
      call rzero(cu,nb)
      do kk=1,ntr
         bcu(kk) = vlsc2(uu,fac_b(1+(kk-1)*(nb+1)),nb+1)
         cuu(kk) = vlsc2(u,fac_c(1+(kk-1)*(nb+1)),nb+1)
      enddo
      call col4(tmp,bcu,cuu,cp_w,ntr) 
      call mxm(fac_a,nb+1,tmp,ntr,tmpcu,1)

      ! exact tensor-vector multiplication
      call rzero(cu,nb)
      do k=kc1,kc2
      do j=jc1,jc2
      do i=ic1,ic2
         cu(i)=cu(i)+cl(i,j,k)*uu(j)*uu(k)
      enddo
      enddo
      enddo
      
      ! difference
      call sub3(cu_diff,cu,tmpcu(0),nb)
      ! error
      cu_err = sqrt(vlsc2(cu_diff,cu_diff,nb))

      write(6,*) nid,cu_err,ntr,'cu_err',' ntr'
      write(6,*)'compare'

      do ii=1,nb
         write(6,*)ii,cu_diff(ii),cu(ii),tmpcu(ii-1),'tmpcu'
      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine CP_ALS_old(cl,mm,nn)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real fcm(0:mm*nn-1,3)
      real fcmpm(nn*nn,3)
      real lsm(nn*nn,3),lsminv(nn*nn,3)
      real tmp(nn*nn),tmp_wrk(nn)
      real cp_weight(nn)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      real lsr(mm*nn)
      real relerr,norm_c
      integer mode,maxit,local_size
      integer mm,nn

      maxit = 1000

      local_size = (ic2-ic1+1)*(jc2-jc1+1)*(kc2-kc1+1)
      norm_c = vlsc2(cl(ic1,jc1,kc1),cl(ic1,jc1,kc1),local_size)

      call rand_initial(fcm,mm,nn)
c     fcm(0,2) = 5.383424352600571e-01;
c     fcm(2,2) = 7.817552875318368e-02;
c     fcm(1,2) = 9.961347166268855e-01;
c     fcm(3,2) = 4.426782697754463e-01;
c     fcm(0,3) = 1.066527701805844e-01;
c     fcm(2,3) = 4.634224134067444e-03;
c     fcm(1,3) = 9.618980808550537e-01;
c     fcm(3,3) = 7.749104647115024e-01;
c     fcm(0,2) = 5.470088922863450e-01 
c     fcm(1,2) = 2.963208056077732e-01 
c     fcm(0,3) = 3.684845964903365e-01 
c     fcm(1,3) = 6.256185607296904e-01 
c     fcm(2,2) = 7.446928070741562e-01 
c     fcm(3,2) = 1.889550150325445e-01 
c     fcm(2,3) = 7.802274351513768e-01 
c     fcm(3,3) = 8.112576886578526e-02 
c     fcm(4,2) = 6.867754333653150e-01
c     fcm(5,2) = 1.835111557372697e-01
c     fcm(4,3) = 9.293859709687300e-01
c     fcm(5,3) = 7.757126786084023e-01

      do mode=2,3
         call set_product_matrix(fcm(0,mode),fcmpm(1,mode),mm,nn)
         do jj=0,mm*nn-1
         write(6,*)jj,fcm(jj,mode),'check fcm'
         enddo
         do jj=1,nn*nn
         write(6,*)jj,fcmpm(jj,mode),'check fcmpm'
         enddo
      enddo


      do ii=1,10!maxit
         do mode=1,3
            call mttkrp(lsr,cl,fcm,mm,nn,mode)
            call set_lsm(lsm,fcmpm,mode,nn)
            call invmat(lsminv(1,mode),tmp,lsm(1,mode),tmp_wrk,nn)
            if (mode.eq.1) then
               write(6,*)'mode',mode
               do jj=1,nn*nn
               write(6,*)jj,lsminv(jj,mode),'check lsminv results'
               enddo
               do jj=1,mm*nn
               write(6,*)jj,lsr(jj),'check lsr results'
               enddo
               do jj=1,nn
                  call mxm(lsr,mm,lsminv(1+(jj-1)*nn,mode),nn,
     $            fcm(0+(mm)*(jj-1),mode),1)
               enddo
               do jj=0,mm*nn-1
               write(6,*)jj,fcm(jj,mode),'check fcm results'
               enddo
            elseif (mode.ne.1) then
               write(6,*)'mode',mode
               do jj=1,nn*nn
               write(6,*)jj,lsminv(jj,mode),'check lsminv results'
               enddo
               do jj=1,mm*nn
               write(6,*)jj,lsr(jj),'check lsr results'
               enddo
               do jj=1,nn
                  call mxm(lsr,mm,lsminv(1+(jj-1)*nn,mode),nn,
     $            fcm(0+(mm)*(jj-1),mode),1)
               enddo
               do jj=0,mm*nn-1
               write(6,*)jj,fcm(jj,mode),'check fcm results'
               enddo
            endif
            call compute_cp_weight(cp_weight,fcm(0,mode),mm,nn)
            call mode_normalize(fcm(0,mode),mm,nn)
            call set_product_matrix(fcm(0,mode),fcmpm(1,mode),mm,nn)
         enddo
         write(6,*)'check cp_w'
         do jj=1,ltr
            write(6,*)jj,cp_weight(jj),'check'
         enddo
         write(6,*)'check fcm outside'
         do jj=0,mm*nn-1
         write(6,*)jj,fcm(jj,1),'check fcm results'
         enddo
         do jj=0,mm*nn-1
         write(6,*)jj,fcm(jj,2),'check fcm results'
         enddo
         do jj=1,nn*nn
         write(6,*)jj,lsm(jj,mode),'check lsminv results'
         enddo
         call compute_relerr(relerr,lsr,fcm(0,3),lsm(1,3),
     $                       fcmpm(1,3),cp_weight,norm_c,mm,nn)
c        call compute_cp_weight(cp_w,fcm(0,3),lub+1,ltr) 
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mttkrp(lsr,cl,cm,cm2,ic1,ic2,jc1,jc2,kc1,kc2,
     $                  fcm,mm,nn,mode)

      ! matricized tensor times khatri-rhao product

      real lsr(mm*nn)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      real cm(ic1:ic2,jc1:jc2,nn)
      real cm2(ic1:ic2,nn,kc1:kc2)
      real fcm(0:mm*nn-1,3)
      integer ic1,ic2,jc1,jc2,kc1,kc2,mode,tr

      call rzero(lsr,mm*nn)
      if (mode.eq.1) then

         ! construct temporary mttkrp
         do tr=1,nn
            call mxm(cl,(ic2-ic1+1)*(jc2-jc1+1),
     $               fcm(0+(mm)*(tr-1),3),(kc2-kc1+1),
     $               cm(ic1,jc1,tr),1)
         enddo

         ! temporary mttkrp with factor matrix
         do tr=1,nn
            call mxm(cm(ic1,jc1,tr),(ic2-ic1+1),
     $               fcm(0+(mm)*(tr-1),2),(jc2-jc1+1),
     $               lsr(1+(mm)*(tr-1)),1)
         enddo

      elseif (mode.eq.2) then

         ! construct temporary mttkrp
         do tr=1,nn
            call mxm(cl,(ic2-ic1+1)*(jc2-jc1+1),
     $               fcm(0+(mm)*(tr-1),3),(kc2-kc1+1),cm(ic1,jc1,tr),1)
         enddo

         ! temporary mttkrp with factor matrix
         do tr=1,nn
            idx = 1
            do jj=jc1,jc2
               lsr(idx+(mm)*(tr-1)) = vlsc2(cm(1,jj,tr),
     $         fcm(0+(mm)*(tr-1),1),mm)
               idx = idx + 1
            enddo
         enddo

      elseif (mode.eq.3) then

         ! construct temporary mttkrp
         do kk=kc1,kc2
         do tr=1,nn
            call mxm(cl(ic1,jc1,kk),(ic2-ic1+1),
     $               fcm(0+(mm)*(tr-1),2),(jc2-jc1+1),
     $               cm2(1,tr,kk),1)
         enddo
         enddo

         ! temporary mttkrp with factor matrix
         do tr=1,nn
            idx = 1
            do kk=kc1,kc2
               lsr(idx+(mm)*(tr-1)) = vlsc2(cm2(1,tr,kk),
     $         fcm(0+(mm)*(tr-1),1),mm)
               idx = idx + 1
            enddo
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_lsm(lsm,fcmpm,mode,nn)

      real lsm(nn*nn,3)
      real fcmpm(nn*nn,3)
      integer mode,idx,nn

      ! Hadamard product

      if (mode.eq.1) then 
         do ii=1,nn
            idx = 1+(ii-1)*nn
            call col3(lsm(idx,1),fcmpm(idx,2),fcmpm(idx,3),nn)
         enddo
      elseif (mode.eq.2) then
         do ii=1,nn
            idx = 1+(ii-1)*nn
            call col3(lsm(idx,2),fcmpm(idx,1),fcmpm(idx,3),nn)
         enddo
      elseif (mode.eq.3) then
         do ii=1,nn
            idx = 1+(ii-1)*nn
            call col3(lsm(idx,3),fcmpm(idx,1),fcmpm(idx,2),nn)
         enddo
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine compute_relerr(relerr,lsr,fcm,lsm,fcmpm,cp_weight,
     $                          norm_c,norm_c0,mm,nn)

      real relerr,norm_c,norm_c0
      real lsr(mm,nn),lsm(nn,nn)
      real fcm(mm,nn),fcmpm(nn,nn)
      real cp_weight(nn)

      real inner_prod,norm_approx
      real tmp1(mm,nn),tmp2(nn,nn)
      real tmp3(nn),tmp4(nn)
      real tmp5(mm),tmp6(mm)
      integer mm,nn

      inner_prod=0.
c     compute the inner product between approximated tensor and
c     exact tensor
      do ii=1,nn
         call col3(tmp1(1,ii),lsr(1,ii),fcm(1,ii),mm)
      enddo
      call rone(tmp6,mm)
      call mxm(tmp1,mm,cp_weight,nn,tmp5,1)
      inner_prod = vlsc2(tmp6,tmp5,mm)

      norm_approx=0.
      ! compute the frobenius norm of the approximated tensor
      do ii=1,nn
         call col3(tmp2(1,ii),lsm(1,ii),fcmpm(1,ii),nn)
      enddo

      call mxm(tmp2,nn,cp_weight,nn,tmp4,1)
      norm_approx = vlsc2(tmp4,cp_weight,nn)

      relerr = sqrt((norm_c0 - 2*inner_prod + norm_approx)/(norm_c))

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_residual(residual,lsr,fcm,lsm,fcmpm,cp_weight,
     $                          norm_tens,mm,nn)

      real residual
      real lsr(mm,nn),lsm(nn,nn)
      real fcm(mm,nn),fcmpm(nn,nn)
      real cp_weight(nn)
      real norm_tens,inner_prod,norm_approx
      integer mm,nn

      inner_prod=0.

c     compute the inner product between approximated tensor and
c     exact tensor
      do ii=1,nn
         call col2c(lsr(1,ii),fcm(1,ii),cp_weight(ii),mm)
      enddo
      inner_prod = vlsum(lsr,mm*nn)

      norm_approx=0.

c     compute the frobenius norm of the approximated tensor
      do ii=1,nn
         call col2c(lsm(1,ii),fcmpm(1,ii),cp_weight(ii),nn)
      enddo
      do jj=1,nn
      do ii=1,nn
         lsm(ii,jj) = lsm(ii,jj) * cp_weight(ii)
      enddo
      enddo
      norm_approx = vlsum(lsm,nn**2)

      residual = sqrt(abs(norm_tens - 2*inner_prod + norm_approx))
      return
      end
c-----------------------------------------------------------------------
      subroutine compute_cp_weight(cp_weight,aa,m,n)

      real cp_weight(n),aa(m,n)
      integer m,n

      do jj=1,n
         cp_weight(jj) = sqrt(vlsc2(aa(1,jj),aa(1,jj),m))
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mode_normalize(aa,m,n)

      real aa(m,n)
      real vlngth
      integer m,n

      do jj=1,n
         vlngth = vlsc2(aa(1,jj),aa(1,jj),m)
         vlngth = 1./sqrt(vlngth)
         do ii=1,m
            aa(ii,jj) = aa(ii,jj)*vlngth
         enddo
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine rand_initial(fcm,mm,nn)

      real fcm(mm*nn)
      integer,parameter :: seed = 86456
  
      call srand(seed)
      
      do ii=1,mm*nn
         fcm(ii) = rand()
      enddo
      call mode_normalize(fcm,mm,nn)
c     call check_normalize(fcm(1,jj),mm,nn)

      return
      end
c-----------------------------------------------------------------------
      subroutine check_normalize(aa,m,n)

      real aa(m,n)
      real vlngth
      integer m,n

      do jj=1,n
         vlngth = vlsc2(aa(1,jj),aa(1,jj),m)
         write(6,*)jj,vlngth,'jj'
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine set_product_matrix(bb,aa,m,n)

      real bb(n,n),aa(m,n)

      do jj=1,n
         do ii=1,n
            bb(ii,jj) = vlsc2(aa(1,ii),aa(1,jj),m)
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine read_cp_weight

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrsetmode/ wk(ltr+1)

      logical ifexist

      inquire (file='./ops/lambda',exist=ifexist)
      if (ifexist) call read_serial(cp_w,ntr,'./ops/lambda ',wk,nid)

      return
      end
c-----------------------------------------------------------------------
      subroutine read_cp_mode

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrsetmode/ wk1(lt)

      if (nio.eq.0) write (6,*) 'reading A1...'
      call read_mat_serial(cua,nb+1,ntr,'./ops/cua ',mb,ntr,wk1,nid)

      if (nio.eq.0) write (6,*) 'reading A2...'
      call read_mat_serial(cub,nb+1,ntr,'./ops/cub ',mb+1,ntr,wk1,nid)

      if (nio.eq.0) write (6,*) 'reading A3...'
      call read_mat_serial(cuc,nb+1,ntr,'./ops/cuc ',mb+1,ntr,wk1,nid)

      return
      end
c-----------------------------------------------------------------------
      subroutine case1_tensor(tensor_1)

      real tensor_1(2,2,2)

      tensor_1(1,1,1) = 1
      tensor_1(2,2,1) = 1
      tensor_1(1,2,2) = 1
      tensor_1(2,1,2) = -1

c     set this in rom.f for testing
c     call case1_tensor(tensor_1)
c     ic1=1
c     ic2=2
c     jc1=0
c     jc2=1
c     kc1=0
c     kc2=1
      return
      end
c-----------------------------------------------------------------------
      subroutine case2_tensor(tensor_1)

      real tensor_1(3,4,4)
      integer ii,jj,kk

      do kk=1,4
         do jj=1,4
            do ii=1,3
               tensor_1(ii,jj,kk) = ii+(jj-1)*3+(kk-1)*12
            enddo
         enddo
      enddo
c     fcm(0,2) = 3.804458469753567e-01 
c     fcm(1,2) = 5.678216407252211e-01 
c     fcm(2,2) = 7.585428956306361e-02
c     fcm(3,2) = 5.395011866660715e-02

c     fcm(0,3) =1.656487294997809e-01 
c     fcm(1,3) =6.019819414016365e-01 
c     fcm(2,3) =2.629712845401443e-01 
c     fcm(3,3) =6.540790984767823e-01 

c     fcm(4,2) = 5.307975530089727e-01 
c     fcm(5,2) = 7.791672301020112e-01 
c     fcm(6,2) = 9.340106842291830e-01 
c     fcm(7,2) = 1.299062084737301e-01 
c                                      
c     fcm(4,3) = 6.892145031400078e-01 
c     fcm(5,3) = 7.481515928237095e-01 
c     fcm(6,3) = 4.505415985024978e-01 
c     fcm(7,3) = 8.382137799693257e-02 

c     fcm(8,2) = 5.688236608721927e-01 
c     fcm(9,2) = 4.693906410582058e-01
c     fcm(10,2) = 1.190206950124140e-02
c     fcm(11,2) = 3.371226443988815e-01
c                                      
c     fcm(8,3) = 2.289769687168188e-01
c     fcm(9,3) = 9.133373615016696e-01
c     fcm(10,3) = 1.523780189692230e-01
c     fcm(11,3) = 8.258169774895474e-01

c     fcm(12,2) = 1.621823081932428e-01
c     fcm(13,2) = 7.942845406839070e-01
c     fcm(14,2) = 3.112150420448049e-01
c     fcm(15,2) = 5.285331355062127e-01
c                                       
c     fcm(12,3) = 5.383424352600571e-01
c     fcm(13,3) = 9.961347166268855e-01
c     fcm(14,3) = 7.817552875318368e-02
c     fcm(15,3) = 4.426782697754463e-01

      ! sete this in rom.f for testing
c     call case2_tensor(tensor_2)

c     ic1=1
c     ic2=3
c     jc1=0
c     jc2=3
c     kc1=0
c     kc2=3
      return
      end
c-----------------------------------------------------------------------
      subroutine als_core(cp_a,cp_b,cp_c,cp_w,cl,
     $                    ic1,ic2,jc1,jc2,kc1,kc2,mm,nn)

c     This subroutine requires tensor cl, indices ic1, ic2, jc1, jc2,
c     kc1, kc2, mm and nn where mm and nn are the dimensions of the
c     factor matrices. It returns cp_a, cp_b, cp_c ,cp_w.

      real cp_a(mm*nn),cp_b(mm*nn),cp_c(mm*nn),cp_w(nn)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)

      real cmr(ic1:ic2,jc1:jc2,nn),crm(ic1:ic2,nn,kc1:kc2)
      real fcm(mm*nn,3),fcmpm(nn*nn,3)
      real lsm(nn*nn,3),lsminv(nn*nn,3),lsr(mm*nn)
      real tmp1(nn*nn)
      real cp_res,cp_relres,cp_fit,cp_tol,pre_relres,reldiff

      integer tmp2(nn),tmp3(nn)
      integer mode,maxit,mm,nn

      real norm_c,norm_c0
      common /scrtens_norm/ norm_c,norm_c0

      if (nid.eq.0) write(6,*) 'inside als_core'

      maxit = 500
      cp_tol = 1e-12
      pre_relres = 1

      call copy(fcm(1,2),cp_b,mm*nn)
      call copy(fcm(1,3),cp_c,mm*nn)

      ! Compute A^TA
      do mode=2,3
         call set_product_matrix(fcmpm(1,mode),fcm(1,mode),mm,nn)
      enddo

      ! ALS
      do ii=1,maxit
         do mode=1,3
            call mttkrp(lsr,cl,cmr,crm,ic1,ic2,jc1,jc2,kc1,kc2,
     $                  fcm,mm,nn,mode)
            call set_lsm(lsm,fcmpm,mode,nn)
            call invmat(lsminv(1,mode),tmp1,lsm(1,mode)
     $           ,tmp2,tmp3,nn)
            do jj=1,nn
               call mxm(lsr,mm,lsminv(1+(jj-1)*nn,mode),nn,
     $         fcm(1+(mm)*(jj-1),mode),1)
            enddo
            call compute_cp_weight(cp_w,fcm(1,mode),mm,nn)
            call mode_normalize(fcm(1,mode),mm,nn)
            call set_product_matrix(fcmpm(1,mode),fcm(1,mode),mm,nn)
         enddo

         ! compute the relative error
         call compute_residual(cp_res,lsr,fcm(1,3),lsm(1,3),
     $                       fcmpm(1,3),cp_w,norm_c0,mm,nn)
         cp_relres = cp_res/sqrt(norm_c)
         cp_fit = 1-cp_relres
         reldiff = abs(pre_relres-cp_relres)/pre_relres

         pre_relres = cp_relres
         if (nid.eq.0) write(6,1)'iter:',ii,'rel difference:',reldiff,
     $   'relative residual:',cp_relres,'residual:',cp_res,'rank:',nn
         if (cp_relres.lt.cp_tol.OR.ii.ge.maxit.OR.
     $       reldiff.le.1e-12) then
            exit
         endif
      enddo

      call copy(cp_a,fcm(1,1),mm*nn)
      call copy(cp_b,fcm(1,2),mm*nn)
      call copy(cp_c,fcm(1,3),mm*nn)

    1 format(a,i4,x,a,1p1e13.5,x,a,1p1e13.5,x,a,1p1e13.5,x,a,i4)
      if (nid.eq.0) write(6,*) 'exitting als_core'

      return
      end
c-----------------------------------------------------------------------
      subroutine set_c_slice(cj0,c0k,cl,ic1,ic2,jc1,jc2,kc1,kc2,nb)

c     This subroutine extracts the two slices of the tensor C
c     one is c(i,j,0) and the other is c(i,0,k)

      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      real cj0(ic1:ic2,jc1:jc2)
      real c0k(ic1:ic2,kc1:kc2)

      if (nid.eq.0) write(6,*) 'inside set_c_slice'

      do j=0,nb
      do i=1,nb
         cj0(i,j) = cl(i,j,0)
      enddo
      enddo

      do k=0,nb
      do i=1,nb
         c0k(i,k) = cl(i,0,k)
      enddo
      enddo

      if (nid.eq.0) write(6,*) 'exitting set_c_slice'

      return
      end
c-----------------------------------------------------------------------
      subroutine set_c_core(cl0,cl,ic1,ic2,jc1,jc2,kc1,kc2,nb)

      ! TODO: fix the for loop

      real cl0(ic1:ic2,jc1+1:jc2,kc1+1:kc2)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      integer ic1,ic2,jc1,jc2,kc1,kc2,nb

      if (nid.eq.0) write(6,*) 'inside set_c_core'

      do k=1,nb
      do j=1,nb
      do i=1,nb
         cl0(i,j,k) = cl(i,j,k)
      enddo
      enddo
      enddo

      if (nid.eq.0) write(6,*) 'exitting set_c_core'

      return
      end
c-----------------------------------------------------------------------
      subroutine als_quad(cp_a,cp_b,cp_c,cp_w,cl,
     $                    ic1,ic2,jc1,jc2,kc1,kc2,mm,nn)

c     This subroutine requires tensor cl, indices ic1, ic2, jc1, jc2,
c     kc1, kc2, mm and nn where mm and nn are the dimensions of the
c     factor matrices. It returns cp_a, cp_b, cp_c ,cp_w.

      real cp_a(mm*nn),cp_b(mm*nn),cp_c(mm*nn),cp_w(nn)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)

      real cmr(ic1:ic2,jc1:jc2,nn),crm(ic1:ic2,nn,kc1:kc2)
      real fcm(mm*nn,3),fcmpm(nn*nn,3)
      real lsm(nn*nn,3),lsminv(nn*nn,3),lsr(mm*nn)
      real tmp1(nn*nn),tmp2(mm*nn)
      real cp_res,cp_relres,cp_fit,cp_tol,pre_relres,reldiff
      real sub_para

      integer tmp3(nn),tmp4(nn)
      integer mode,maxit,mm,nn

      real norm_c,norm_c0
      common /scrtens_norm/ norm_c,norm_c0

      if (nid.eq.0) write(6,*) 'inside als_quad'

      maxit = 500
      cp_tol = 1e-12
      pre_relres = 1
      sub_para = 0.7

      ! Compute A^TA
c     call rand_initial(fcm(1,1),mm,nn)
c     call copy(fcm(1,2),fcm((mm*nn/2)+1,1),mm*nn/2)
c     call copy(fcm((mm*nn/2)+1,2),fcm(1,1),mm*nn/2)
c     call chsign(fcm((mm*nn/2)+1,2),mm*nn/2)
      call copy(fcm(1,1),cp_a,mm*nn)
      call copy(fcm(1,2),cp_b,mm*nn)

      do mode=1,2
         call set_product_matrix(fcmpm(1,mode),fcm(1,mode),mm,nn)
      enddo

      ! ALS with quad factor
      call copy(tmp2,fcm(1,1),mm*nn)
      do ii=1,maxit
            ! solve for mode 3 first
            call mttkrp(lsr,cl,cmr,crm,ic1,ic2,jc1,jc2,kc1,kc2,
     $                  fcm,mm,nn,3)
            call set_lsm(lsm,fcmpm,3,nn)
            call invmat(lsminv(1,3),tmp1,lsm(1,3)
     $           ,tmp3,tmp4,nn)
            do jj=1,nn
               call mxm(lsr,mm,lsminv(1+(jj-1)*nn,3),nn,
     $         fcm(1+(mm)*(jj-1),3),1)
            enddo
            call mode_normalize(fcm(1,3),mm,nn)
            call set_product_matrix(fcmpm(1,3),fcm(1,3),mm,nn)

         do kk=1,10
               call mttkrp(lsr,cl,cmr,crm,ic1,ic2,jc1,jc2,kc1,kc2,
     $                     fcm,mm,nn,1)
               call set_lsm(lsm,fcmpm,1,nn)
               call invmat(lsminv(1,1),tmp1,lsm(1,1)
     $                     ,tmp3,tmp4,nn)
               do jj=1,nn
                  call mxm(lsr,mm,lsminv(1+(jj-1)*nn,1),nn,
     $            fcm(1+(mm)*(jj-1),1),1)
               enddo
               call compute_cp_weight(cp_w,fcm(1,1),mm,nn)
               call mode_normalize(fcm(1,1),mm,nn)

               call add2sxy(fcm(1,1),sub_para,tmp2,1-sub_para,mm*nn)
               call copy(tmp2,fcm(1,1),mm*nn)

               call copy(fcm(1,2),fcm((mm*nn/2)+1,1),mm*nn/2)
               call copy(fcm((mm*nn/2)+1,2),fcm(1,1),mm*nn/2)
               call chsign(fcm((mm*nn/2)+1,2),mm*nn/2)

               call set_product_matrix(fcmpm(1,2),fcm(1,2),mm,nn)
         enddo
         call set_product_matrix(fcmpm(1,1),fcm(1,1),mm,nn)

         ! compute the relative error
         call compute_residual(cp_res,lsr,fcm(1,1),lsm(1,1),
     $                       fcmpm(1,1),cp_w,norm_c0,mm,nn)
         cp_relres = cp_res/sqrt(norm_c)
         cp_fit = 1-cp_relres
         reldiff = abs(pre_relres-cp_relres)/pre_relres

         pre_relres = cp_relres
         if (nid.eq.0) write(6,1)'iter:',ii,'rel difference:',reldiff,
     $   'relative residual:',cp_relres,'residual:',cp_res,'rank:',nn
         if (cp_relres.lt.cp_tol.OR.ii.ge.maxit.OR.reldiff.le.1e-5) then
            exit
         endif
      enddo

      call copy(cp_a,fcm(1,1),mm*nn)
      call copy(cp_b,fcm(1,2),mm*nn)
      call copy(cp_c,fcm(1,3),mm*nn)

    1 format(a,i4,x,a,1p1e13.5,x,a,1p1e13.5,x,a,1p1e13.5,x,a,i4)
      if (nid.eq.0) write(6,*) 'exitting als_quad'

      return
      end
c-----------------------------------------------------------------------
      subroutine force_skew(cl,ic1,ic2,jc1,jc2,kc1,kc2)

      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      real cl0t(jc1+1:jc2,ic1:ic2,kc1+1:kc2)
      integer ic1,ic2,jc1,jc2,kc1,kc2

      if (nid.eq.0) write(6,*) 'inside force_skew'

      do k=kc1+1,kc2
         call transpose(cl0t(1,1,k),(jc2-jc1),cl(1,1,k),(ic2-ic1+1))
      do j=jc1+1,jc2
      do i=ic1,ic2
         cl(i,j,k) = 0.5*(cl(i,j,k)-cl0t(i,j,k))
      enddo
      enddo
      enddo

      if (nid.eq.0) write(6,*) 'exitting force_skew'

      return
      end
