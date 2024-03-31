c-----------------------------------------------------------------------
      subroutine breduce(a,n,m)

      include 'SIZE'
      include 'LMOR'

      ! global reduction of vector a of length n (batch size m)

      common /workbr/ v(lbat),w(lbat)

      real a(n)

      if (nid.eq.0) then
         write (6,*) 'calling breduce with n=',n,' m=',m
      endif

      if (m.le.0) return
      if (nid.eq.0) write (6,*) 'wp 0',m,lbat

      m = min(m,lbat)
      if (nid.eq.0) write (6,*) 'wp 1',m

      k = n / m
      if (nid.eq.0) write (6,*) 'wp 2',k,n,m
      if (n.ne.k*m) k=k+1
      if (nid.eq.0) write (6,*) 'wp 3',k
      m = n / k
      if (nid.eq.0) write (6,*) 'wp 4',m

      nrem = n - k*m
      if (nid.eq.0) write (6,*) 'wp 5',nrem

      if (nid.eq.0) then
         write (6,*) 'starting loop with k=',k,' m=',m,' nrem=',nrem
      endif

      ia=1
      do i=1,k-nrem
         call gop(a(ia),w,'+  ',m)
         ia=ia+m
      enddo

      do i=1,nrem
         call gop(a(ia),w,'+  ',m+1)
         ia=ia+m+1
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine brprofile

      include 'SIZE'
      include 'PARALLEL'
      include 'LMOR'

      common /workbr/ v(lbat),w(lbat)

      ntrial = 1024

      do n=1,lbat
         time_min = 1.0e+10
         time_max = 0.0
         time_avg = 0.0
         time_avg2 = 0.0
         do ir=1,n
            v(ir) = rand()
         enddo
         do i=1,ntrial
            start_time=dnekclock()
            call breduce(v,n,n)
            end_time=dnekclock()
            time = end_time - start_time
            time_min = min(time_min,time)
            time_max = max(time_max,time)
            time_avg = time_avg + time
            time_avg2 = time_avg2 + time*time
         enddo
         time_avg = time_avg / ntrial
         time_avg2 = time_avg2 / ntrial
         time_std = sqrt(time_avg2 - time_avg*time_avg)
         if (nio.eq.0) then
            write (6,*) n, time_min, time_avg, time_max, time_std, np
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
