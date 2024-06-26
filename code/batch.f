c-----------------------------------------------------------------------
      subroutine breduce(a,n,m)

      include 'SIZE'
      include 'LMOR'

      ! global reduction of vector a of length n (batch size m)

      common /workbr/ v(lbat),w(lbat)

      real a(n)

      if (m.le.0) return
      if (m.gt.lbat) m = lbat
      i=1

      do while (i.le.n)
         call gop(a(i),w,'+  ',min(m,n-i+1))
         i=i+m
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
