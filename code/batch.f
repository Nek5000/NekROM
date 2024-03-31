c-----------------------------------------------------------------------
      subroutine breduce(a,nn,mm)

      include 'SIZE'
      include 'LMOR'

      ! global reduction of vector a of length n (batch size m)

      common /workbr/ v(lbat),w(lbat)

      real a(nn)

      if (nid.eq.0) then
         write (6,*) 'calling breduce with n=',nn,' m=',mm
      endif

      if (mm.le.0) return
      write (6,*) 'wp 0',mm,lbat

      if (mm.gt.lbat) mm=lbat
      if (nid.eq.0) write (6,*) 'wp 1',mm

      kk = nn / mm
      write (6,*) 'wp 2',kk,nn,mm
      if (nn.ne.kk*mm) kk=kk+1
      write (6,*) 'wp 3',kk
      write (6,*) 'breduce_test_idivide',10000/10000
      mm = nn / kk
      write (6,*) 'wp 4',mm

      nnrem = nn - kk*mm
      if (nid.eq.0) write (6,*) 'wp 5',nnrem

      if (nid.eq.0) then
         write (6,*) 'starting loop with k=',kk,' m=',mm,' nrem=',nnrem
      endif

      ia=1
c     do i=1,kk-nnrem
c        call gop(a(ia),w,'+  ',mm)
c        ia=ia+mm
c     enddo

c     do i=1,nnrem
c        call gop(a(ia),w,'+  ',mm+1)
c        ia=ia+mm+1
c     enddo

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
