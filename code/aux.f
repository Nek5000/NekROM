c-----------------------------------------------------------------------
      subroutine readeig(evec)

      include 'SIZE'
      include 'MOR'

      common /scrk5/ w(lx1*ly1*lz1*lelt)

      real evec(ls,nb)

      call rzero(evec,ls*nb)

      if (nid.eq.0) then
         open (unit=12,file='evectors.dat')
         read (12,*) (evec(k,1),k=1,ls*nb)
         close (unit=12)
      endif

      call gop(evec,w,'+  ',ls*nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine factor3(mq,mp,mr,m)

      integer dmin,d

      n=m
      l=nint(real(n)**(1/3))

      dmin=n
      imin=-1

      do i=1,n
          d=abs(n-i**3)
          if (d.lt.dmin.and.mod(n,i).eq.0) then
              dmin=d
              imin=i
          endif
      enddo

      mp=imin
      n=n/mp

      dmin=n
      imin=-1

      do i=1,n
          d=abs(n-i*i)
          if (d.lt.dmin.and.mod(n,i).eq.0) then
              dmin=d
              imin=i
          endif
      enddo

      mq=imin
      mr=n/mq

      if (nio.eq.0) write (6,*) 'mp,mq,mr,mp*mq*mr',mp,mq,mr,mp*mq*mr

      return
      end
c-----------------------------------------------------------------------
      subroutine setpart(mps,mp,n)

      integer mps(mp)

      do i=0,mp-1
         mps(i+1)=n/mp+max(mod(n,mp)-i,0)/max(mod(n,mp)-i,1)
         mps(i+1)=mps(i+1)+mps(max(i,1))*max(i,0)/max(i,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setpart3(mps,mqs,mrs,mp,mq,mr,nb)

      integer mps(mp),mqs(mq),mrs(mr)

      call setpart(mps,mp,nb)
      call setpart(mqs,mq,nb+1)
      call setpart(mrs,mr,nb+1)

      return
      end
c-----------------------------------------------------------------------
      function i2p(i,mps)

      integer mps(1)

      i2p=0

      do while (mps(i2p+1).ne.0)
         i2p=i2p+1
         if (i.le.mps(i2p)) return
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ijk2pqr(ip,iq,ir,i,j,k,mps,mqs,mrs)

      real mps(1),mqs(1),mrs(1)

      ip=i2p(i,mps)
      iq=i2p(j+1,mqs)
      ir=i2p(k+1,mrs)

      return
      end
c-----------------------------------------------------------------------
      function ijk2pid(i,j,k,mps,mqs,mrs,mp,mq,mr)

      real mps(1),mqs(1),mrs(1)

      call ijk2pqr(ip,iq,ir,i,j,k,mps,mqs,mrs)
      ijk2pid=(ip-1)+mp*(iq-1)+mp*mq*(ir-1)

      return
      end
c-----------------------------------------------------------------------
      subroutine setrange(mps,mqs,mrs,mp,mq,mr)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      integer mps(1),mqs(1),mrs(1)

      ip=mod(nid,mp)
      iq=mod(nid/mp,mq)
      irr=   (nid/mp)/mq

      i0=mps(max(ip,1))*max(ip,0)/max(ip,1)+1
      i1=mps(ip+1)

      j0=mqs(max(iq,1))*max(iq,0)/max(iq,1)
      j1=mqs(iq+1)-1

      k0=mrs(max(irr,1))*max(irr,0)/max(irr,1)
      k1=mrs(irr+1)-1

      return
      end
c-----------------------------------------------------------------------
      subroutine ijk2l(l,i,j,k)
      
      include 'SIZE'
      include 'MOR'

      il=i-i0
      jl=j-j0
      kl=k-k0

      l=il+jl*(i1-i0+1)+kl*(i1-i0+1)*(j1-j0+1)+1

      return
      end
c-----------------------------------------------------------------------
      subroutine opadd3 (a1,a2,a3,b1,b2,b3,c1,c2,c3)
      include 'SIZE'
      real a1(1),a2(1),a3(1),b1(1),b2(1),b3(1)
      real c1(1),c2(1),c3(1)
      ntot1=lx1*ly1*lz1*nelv
      call add2(a1,b1,c1,ntot1)
      call add2(a2,b2,c2,ntot1)
      if (ldim.eq.3) call add2(a3,b3,c3,ntot1)
      return
      end
c-----------------------------------------------------------------------
      subroutine outmatl(uu,v,w,k)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      real uu(lt),v(lt),w(lt)
      character*7 fname

      if (np.gt.1) call exitti('outmatl works for P=1 only!$',np)

      write(fname,22) k
   22 format('out.',i3.3)
      open(unit=33,file=fname)

      n = lx1*ly1*lz1*nelt
      do i=1,n
         write(33,33) xm1(i,1,1,1),ym1(i,1,1,1),uu(i),v(i)
   33    format(1p4e16.7)
      enddo
      close(33)

      return
      end
c-----------------------------------------------------------------------
      subroutine readc0(c0,n)

      include 'SIZE'
      include 'PARALLEL'

      real c0(n)

      if (nid.eq.(np-1)) then
         open (unit=12,file='cten')
         read (12,*) (c0(k),k=1,n)
         close (unit=12)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine readab(a0,b0,n)

      include 'SIZE'

      real a0(n),b0(n)

      common /scrk5/ w(lx1*ly1*lz1*lelt)

      call rzero(a0,n)
      call rzero(b0,n)

      if (nid.eq.0) then
         open (unit=12,file='amat')
         read (12,*) (a0(k),k=1,n)
         close (unit=12)

         open (unit=12,file='bmat')
         read (12,*) (b0(k),k=1,n)
         close (unit=12)
      endif

      call gop(a0,w,'+  ',n)
      call gop(b0,w,'+  ',n)

      return
      end
c-----------------------------------------------------------------------
      subroutine readic(ic,n)

      include 'SIZE'

      real ic(n)

      if (nid.eq.0) then
         open (unit=12,file='ic')
         read (12,*) (ic(k),k=1,n)
         close (unit=12)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine readbasis(ub,vb,wb,nb)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ub(lt,0:nb), vb(lt,0:nb), wb(lt,0:nb)

      n=lx1*ly1*lz1*nelt

      if (nid.eq.0) then
         open (unit=12,file='basis')

         do j=0,nb
            read (12,*) (ub(i,j),i=1,n)
         enddo

         do j=0,nb
            read (12,*) (vb(i,j),i=1,n)
         enddo

         if (ldim.eq.3) then
         do j=0,nb
            read (12,*) (wb(i,j),i=1,n)
         enddo
         endif

         close (unit=12)
      endif

      return
      end
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
      subroutine dumpic(ic)

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
