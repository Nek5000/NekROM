c-----------------------------------------------------------------------
      subroutine set_cp(cl,fname)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real cl(lcglo)

      character*128 fname
      character*128 fnlint

      if (nio.eq.0) write (6,*) 'inside set_cp'

      call nekgsync
      tcp_time=dnekclock()

      if (iffastc) call exitti('fastc not supported in setc_new$',nb)

      ! set global index
      ic1=1
      ic2=nb
      jc1=0
      jc2=nb
      kc1=0
      kc2=nb

      n=lx1*ly1*lz1*nelv

      call lints(fnlint,fname,128)
      if (nid.eq.0) open (unit=100,file=fnlint)
      if (nio.eq.0) write (6,*) 'setc file:',fnlint

      do k=0,nb
      do j=0,mb
      do i=1,mb
         cel=0.
         if (nid.eq.0) read(100,*) cel
         cel=glsum(cel,1)
         call setc_local(cl,cel,ic1,ic2,jc1,jc2,kc1,kc2,i,j,k)
      enddo
      enddo
      enddo

      write(6,*)'check index',ic1,ic2,jc1,jc2,kc1,kc2,nid
      call nekgsync

      call CP_ALS(cl,cua,cub,cuc,cp_uw,u,nb+1,ntr)

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

      if (nio.eq.0) write (6,*) 'exiting set_cp'

      return
      end
c-----------------------------------------------------------------------
      subroutine CP_ALS(cl,aa,bb,cc,cp_weight,uu,mm,nn)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      real fcm(0:mm*nn-1,3),fcmpm(nn*nn,3)
      real aa(1:mm*nn)
      real bb(1:mm*nn)
      real cc(1:mm*nn)
      real uu(1:m)
      real lsm(nn*nn,3),lsminv(nn*nn,3)
      real tmp(nn*nn),tmp_wrk(nn)
      real lsr(mm*nn)
      real cp_weight(nn)
      real relerr,norm_c,fit,cp_tol,pre_err,rel_diff
      real wk(lb+1)
      integer mode,maxit,local_size
      integer mm,nn
      logical ifexist

      if (nid.eq.0) write(6,*) 'inside cp_als'

      maxit = 500
      cp_tol = 1e-5
      pre_err = 1

      local_size = (ic2-ic1+1)*(jc2-jc1+1)*(kc2-kc1+1)
      norm_c = vlsc2(cl(ic1,jc1,kc1),cl(ic1,jc1,kc1),local_size)

      call rand_initial(fcm,mm,nn)
c     if (nid.eq.0) then
c     write(6,*)'processor 0'
c        do ii=0,mm*nn-1
c        write(6,*)ii,fcm(ii,2),fcm(ii,3),nid,'nid'
c        enddo
c     endif
c     call nekgsync
c     if (nid.eq.1) then
c     write(6,*)'processor 1'
c        do ii=0,mm*nn-1
c        write(6,*)ii,fcm(ii,2),fcm(ii,3),nid,'nid'
c        enddo
c     endif
c     call nekgsync
c     call exitt0

      do mode=2,3
         call set_product_matrix(fcm(0,mode),fcmpm(1,mode),mm,nn)
      enddo


      do ii=1,maxit
         do mode=1,3
            call mttkrp(lsr,cl,fcm,mm,nn,mode)
            call set_lsm(lsm,fcmpm,mode,nn)
            call invmat(lsminv(1,mode),tmp,lsm(1,mode),tmp_wrk,nn)
            if (mode.eq.1) then
               do jj=1,nn
                  call mxm(lsr,mm,lsminv(1+(jj-1)*nn,mode),nn,
     $            fcm(0+(mm)*(jj-1),mode),1)
               enddo
            elseif (mode.ne.1) then
               do jj=1,nn
                  call mxm(lsr,mm,lsminv(1+(jj-1)*nn,mode),nn,
     $            fcm(0+(mm)*(jj-1),mode),1)
               enddo
            endif
            call compute_cp_weight(cp_weight,fcm(0,mode),mm,nn)
            call mode_normalize(fcm(0,mode),mm,nn)
            call set_product_matrix(fcm(0,mode),fcmpm(1,mode),mm,nn)
         enddo
         call compute_relerr(relerr,lsr,fcm(0,3),lsm(1,3),
     $                       fcmpm(1,3),cp_weight,norm_c,mm,nn)
         fit = 1-relerr
         rel_diff = abs(pre_err-relerr)/pre_err
         pre_err = relerr
         if (nid.eq.0) write(6,*) ii, rel_diff, relerr, 'relerr'
         if (relerr.lt.cp_tol.OR.ii.ge.maxit.OR.rel_diff.le.1e-4) then
            exit
         endif
      enddo

c     if (ifrom(1)) then
         call copy(aa,fcm(0,1),mm*nn)
         call copy(bb,fcm(0,2),mm*nn)
         call copy(cc,fcm(0,3),mm*nn)
c        call copy(cp_uw,cp_weight,nn)

c        inquire (file='ops/u0',exist=ifexist)
c        if (ifexist) call read_serial(u,nb+1,'ops/u0 ',wk,nid)
         call check_conv_err(cl,aa,bb,cc,cp_weight,uu)
c     endif
c     if (ifrom(2)) then
c        call copy(cta,fcm(0,1),mm*nn)
c        call copy(ctb,fcm(0,2),mm*nn)
c        call copy(ctc,fcm(0,3),mm*nn)
c        call copy(cp_tw,cp_weight,nn)

c        inquire (file='ops/t0',exist=ifexist)
c        if (ifexist) call read_serial(ut,nb+1,'ops/t0 ',wk,nid)
c        call check_conv_err(cl,cta,ctb,ctc,cp_tw,ut)
c     endif

      call nekgsync

      if (nid.eq.0) write(6,*) 'exit cp_als'

      

      return
      end
c-----------------------------------------------------------------------
      subroutine check_conv_err(cl,fac_a,fac_b,fac_c,cp_weight,uu)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real cu(nb)
      real uu(0:nb)
      real cu_diff(nb), cu_err
      real fac_a((nb+1)*ntr),fac_b((nb+1)*ntr)
      real fac_c((nb+1)*ntr),cp_weight(ntr)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      real cm(ic1:ic2,jc1:jc2)
      real bcu(ntr)
      real cuu(ntr)
      real tmp(ntr)
      real tmpcu(0:nb)

      call rzero(cu,nb)
      do kk=1,ntr
         bcu(kk) = vlsc2(uu,fac_b(1+(kk-1)*(nb+1)),nb+1)
         cuu(kk) = vlsc2(u,fac_c(1+(kk-1)*(nb+1)),nb+1)
      enddo
      call col4(tmp,bcu,cuu,cp_weight,ntr) 
      call mxm(fac_a,nb+1,tmp,ntr,tmpcu,1)

      call rzero(cu,nb)
      do k=kc1,kc2
      do j=jc1,jc2
      do i=ic1,ic2
         cu(i)=cu(i)+cl(i,j,k)*uu(j)*u(k)
      enddo
      enddo
      enddo
      
      call sub3(cu_diff,cu,tmpcu(0),nb)
      cu_err = sqrt(vlsc2(cu_diff,cu_diff,nb))
      write(6,*)nid,cu_err,'cu_err'
      
      write(6,*)'compare'
      do ii=1,nb
         write(6,*)ii,cu(ii),tmpcu(ii-1),cu_diff(ii),'cu, tmpcu',nid
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
      subroutine mttkrp(lsr,cl,fcm,mm,nn,mode)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real lsr(mm*nn)
      real fcm(0:mm*nn-1,3)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      real cm(ic1:ic2,jc1:jc2,ltr)
      real cm2(ic1:ic2,ltr,kc1:kc2)
      integer mode,tr

      if (mode.eq.1) then

         call rzero(lsr,mm*nn)
         ! construct temporary mttkrp
         do tr=1,nn
            call mxm(cl,(ic2-ic1+1)*(jc2-jc1+1),
     $               fcm(kc1+(mm)*(tr-1),3),(kc2-kc1+1),cm(1,0,tr),1)
         enddo

         ! temporary mttkrp with factor matrix
         do tr=1,nn
            call mxm(cm(ic1,jc1,tr),(ic2-ic1+1),
     $               fcm(jc1+(mm)*(tr-1),2),(jc2-jc1+1),
     $               lsr(1+(mm)*(tr-1)),1)
         enddo

      elseif (mode.eq.2) then

         call rzero(lsr,mm*nn)
         ! construct temporary mttkrp
         do tr=1,nn
            call mxm(cl,(ic2-ic1+1)*(jc2-jc1+1),
     $               fcm(kc1+(mm)*(tr-1),3),(kc2-kc1+1),cm(1,0,tr),1)
         enddo

         ! temporary mttkrp with factor matrix
         do tr=1,nn
            do jj=jc1,jc2
               lsr((jj+1)+(mm)*(tr-1)) = vlsc2(cm(1,jj,tr),
     $         fcm(0+(mm)*(tr-1),1),mm)
            enddo
         enddo

      elseif (mode.eq.3) then

         call rzero(lsr,mm*nn)

         ! construct temporary mttkrp
         do kk=kc1,kc2
         do tr=1,nn
            call mxm(cl(ic1,jc1,kk),(ic2-ic1+1),
     $               fcm(jc1+(mm)*(tr-1),2),(jc2-jc1+1),
     $               cm2(1,tr,kk),1) 
         enddo
         enddo

         ! temporary mttkrp with factor matrix
         do tr=1,nn
            do kk=kc1,kc2
               lsr((kk+1)+(mm)*(tr-1)) = vlsc2(cm2(1,tr,kk),
     $         fcm(0+(mm)*(tr-1),1),mm)
            enddo
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_lsm(lsm,fcmpm,mode,nn)

      real fcmpm(nn*nn,3)
      real lsm(nn*nn,3)
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
     $                          norm_c,mm,nn)

      real relerr,norm_c
      real inner_prod,norm_approx 
      real lsr(mm,nn),lsm(nn,nn)
      real fcm(mm,nn),fcmpm(nn,nn)
      real tmp1(mm,nn),tmp2(nn,nn)
      real tmp3(nn),tmp4(nn)
      real tmp5(mm),tmp6(mm)
      real cp_weight(nn)
      integer mm,nn

      inner_prod=0.
      norm_approx=0.

      do ii=1,nn
         call col3(tmp1(1,ii),lsr(1,ii),fcm(1,ii),mm)
      enddo

      call rone(tmp6,mm)
      call mxm(tmp1,mm,cp_weight,nn,tmp5,1)

      inner_prod = vlsc2(tmp6,tmp5,mm)

      do ii=1,nn
         call col3(tmp2(1,ii),lsm(1,ii),fcmpm(1,ii),nn)
      enddo

      call mxm(tmp2,nn,cp_weight,nn,tmp4,1)
      norm_approx = vlsc2(tmp4,cp_weight,nn)

      relerr = sqrt((norm_c - 2*inner_prod + norm_approx)/(norm_c))



      return
      end
c-----------------------------------------------------------------------
      subroutine compute_cp_weight(cp_weight,aa,m,n)

      real aa(m,n)
      real cp_weight(n)
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

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real fcm(mm*nn,3)
      integer,parameter :: seed = 86456
  
      call srand(seed)
      
      do jj=2,3
         do ii=1,mm*nn
            fcm(ii,jj) = rand()
         enddo
         call mode_normalize(fcm(1,jj),mm,nn)

         call check_normalize(fcm(1,jj),mm,nn)
      enddo

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
      subroutine set_product_matrix(aa,bb,m,n)

      real aa(m,n)
      real bb(n,n)

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
