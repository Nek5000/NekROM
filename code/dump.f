c-----------------------------------------------------------------------
      subroutine dumpevec(evec,ns,nb)

      include 'SIZE'

      real evec(ns,nb)

      if (nid.eq.0) then
         open (unit=12,file='ops/evec')

         do j=1,nb
         do i=1,ns
            write (12,*) evec(i,j)
         enddo
         enddo

         close (unit=12)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dumpgram(uu,ns)

      include 'SIZE'

      real uu(ns,ns)

      n=lx1*ly1*lz1*nelt

      if (nid.eq.0) then
         open (unit=12,file='ops/gram')

         do j=1,ns
         do i=1,ns
            write (12,*) uu(i,j)
         enddo
         enddo

         close (unit=12)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dumpbases(ub,vb,wb,nb)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ub(lt,0:nb), vb(lt,0:nb), wb(lt,0:nb)

      n=lx1*ly1*lz1*nelt

      if (nid.eq.0) then
         open (unit=12,file='ops/bases')

         do j=0,nb
         do i=1,n
            write (12,*) ub(i,j)
         enddo
         enddo

         do j=0,nb
         do i=1,n
            write (12,*) vb(i,j)
         enddo
         enddo

         if (ldim.eq.3) then
         do j=0,nb
         do i=1,n
            write (12,*) wb(i,j)
         enddo
         enddo
         endif

         close (unit=12)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dumptens(c,nb)

      include 'SIZE'

      real c(0:nb,0:nb,0:nb)

      if (nid.eq.0) then
         open (unit=50,file='ops/cten')

         do i=0,(nb+1)**3-1
            write (50,*) c(i,0,0)
         enddo

         close (unit=50)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dumpmats(a,b,nb)

      include 'SIZE'

      real a(0:nb,0:nb), b(0:nb,0:nb)

      if (nid.eq.0) then
         open (unit=50,file='ops/amat')

         do i=0,(nb+1)**2-1
            write (50,*) a(i,0)
         enddo

         close (unit=50)

         open (unit=50,file='ops/bmat')

         do i=0,(nb+1)**2-1
            write (50,*) b(i,0)
         enddo

         close (unit=50)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dumpic(ic,nb)

      include 'SIZE'

      real ic(0:nb)

      if (nid.eq.0) then
         open (unit=50,file='ops/ic')

         do i=0,nb
            write (50,*) ic(i)
         enddo

         close (unit=50)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dumpcoef(uuu,nb,k)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)
      real uuu(0:nb)
      character*8 fname

      if (nid .eq. 0) then
         write(fname,22) k
   22 format(i4.4,".out")
         open(unit=33,file=fname)

         do i=1,nb
            write(33,33) uuu(i)
   33    format(1p1e16.7)
         enddo
         close(33)
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine dumpeig(lmbda)
      include 'SIZE'
      include 'MOR'

      real eig(ls),lmbda(ls)

      if (nio.eq.0) write (6,*) 'inside dumpeig'

      if (nio.eq.0) write(6,*)'number of mode:',nb
      if (nid.eq.0) then
         open (unit=50,file='./gram_eig')

         do l = 1,ls
            eig(l)=lmbda(ls-l+1) ! reverse order of eigvalues
            write (50,*) eig(l)
         enddo

         close (unit=50)
      endif


      end
