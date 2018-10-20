c-----------------------------------------------------------------------
      subroutine dumpbasis(ub,vb,wb,nb)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ub(lt,0:nb), vb(lt,0:nb), wb(lt,0:nb)

      n=lx1*ly1*lz1*nelt

      if (nid.eq.0) then
         open (unit=12,file='basis')

         do j=0,nb
         do i=1,n
            write (12,*) ub(i,j)
         enddo
         enddo

         do j=0,nb
         do i=1,n
            write (12,*) ub(i,j)
         enddo
         enddo

         if (ldim.eq.3) then
         do j=0,nb
         do i=1,n
            write (12,*) ub(i,j)
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
         open (unit=50,file='cten')

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
         open (unit=50,file='amat')

         do i=0,(nb+1)**2-1
            write (50,*) a(i,1)
         enddo

         close (unit=50)

         open (unit=50,file='bmat')

         do i=0,(nb+1)**2-1
            write (50,*) b(i,1)
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
         open (unit=50,file='ic')

         do i=0,nb
            write (50,*) ic(i)
         enddo

         close (unit=50)
      endif

      return
      end
c-----------------------------------------------------------------------
