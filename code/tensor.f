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
      subroutine rand_initial()

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real mode_unfold(nb*ltr,3)
      integer,parameter :: seed = 86456
  
      call srand(seed)
      
      do jj=1,3
         do ii=1,nb*ltr
            mode_unfold(ii,jj) = rand()
         enddo
         call mode_normalize(mode_unfold(1,jj),nb,ltr)
      enddo

      do jj=1,3
         call check_normalize(mode_unfold(1,jj),nb,ltr)
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
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine set_cp_weight

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
      subroutine set_cp_mode

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
