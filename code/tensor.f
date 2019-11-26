c-----------------------------------------------------------------------
      subroutine CP_ALS(cl)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real fcm((lub+1)*ltr,3)
      real fcmpm(ltr*ltr,3)
      real lsm(ltr*ltr,3)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      integer mode,maxit

      maxit = 1000

      call rand_initial(fcm)

      do mode=2,3
         call set_product_matrix(fcm(1,mode),fcmpm(1,mode),lub+1,ltr)
      enddo

      do ii=1,2!maxit
         do mode=1,3
            call set_lsm(lsm,fcmpm,mode,ltr)
         enddo
      enddo

      call mttkrp(lsr,cl,fcm,1)
      

      return
      end
c-----------------------------------------------------------------------
      subroutine mttkrp(lsr,cl,fcm,mode)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real lsr((lub+1)*ltr)
      real fcm((lub+1)*ltr,3)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      real cm(ic1:ic2,jc1:jc2,ltr)
      integer mode

      write(6,*)ic1,ic2,jc1,jc2,kc1,kc2,'index'
 
      if (mode.eq.1) then
         do ii=1,ltr
         call mxm(cl,(ic2-ic1+1)*(jc2-jc1+1),
     $            fcm(kc1,3),(kc2-kc1+1),cm(1,1,ii),1)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_lsm(lsm,fcmpm,mode,nn)

      real fcmpm(nn*nn,3)
      real lsm(nn*nn,3)
      integer mode,idx,nn

      if (mode.eq.1) then 
         do ii=1,nn
            idx = 1+(ii-1)*nn
            call col3(lsm(idx,1),fcmpm(idx,2),fcmpm(idx,3),nn)
         enddo
         do ii=1,nn*nn
            write(6,*)ii,lsm(ii,1),fcmpm(ii,2),fcmpm(ii,3),'check'
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
      subroutine rand_initial(fcm)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real fcm((lub+1)*ltr,3)
      integer,parameter :: seed = 86456
  
      call srand(seed)
      
      do jj=2,3
         do ii=1,(lub+1)*ltr
            fcm(ii,jj) = rand()
         enddo
         call mode_normalize(fcm(1,jj),lub+1,ltr)

         call check_normalize(fcm(1,jj),lub+1,ltr)
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
      call read_mat_serial(cua,nb,ntr,'./ops/cua ',mb,ntr,wk1,nid)

      if (nio.eq.0) write (6,*) 'reading A2...'
      call read_mat_serial(cub,nb+1,ntr,'./ops/cub ',mb+1,ntr,wk1,nid)

      if (nio.eq.0) write (6,*) 'reading A3...'
      call read_mat_serial(cuc,nb+1,ntr,'./ops/cuc ',mb+1,ntr,wk1,nid)

      return
      end
