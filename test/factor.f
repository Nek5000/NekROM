c-----------------------------------------------------------------------
      program main

      parameter (n=2400)

      integer factors(n), parts(3), blocks(n)
      integer mps(n),mqs(n),mrs(n)

      nb = 400

      call izero(blocks,n)

      call izero(mps,n)
      call izero(mqs,n)
      call izero(mrs,n)

      call factor(factors,n)

      m=nint(rand(1)*948)

      call factor3(mp,mq,mr,m)
      call setparts(mps,mqs,mrs,mp,mq,mr,nb)

      end
c-----------------------------------------------------------------------
      subroutine factor3(mp,mq,mr,m)

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

      write (6,*) 'mp,mq,mr,mp*mq*mr',mp,mq,mr,mp*mq*mr

      return
      end
c-----------------------------------------------------------------------
      subroutine setpart(mps,mp,n)

      integer mps(mp)

      do i=0,n-1
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
      subroutine izero(a,n)

      integer a(n)

      do i=1,n
         a(i) = 0
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine clocal(c,cc,u,i0,i1,j0,j1,k0,k1,nb)

      real c(nb),cc(1),u(nb+1),a(1)

      common /scrk1/ work(100)

c     call rzero(c,nb)

      l=0

      do k=k0,k1
      do j=j0,j1
         ajk=a(j)*a(k)
         do i=i0,i1
            l=l+1
            c(i)=c(i)+cc(l)*ajk
         enddo
      enddo
      enddo

c     call gop(c,work,'+  ',nb)

      return
      end
c-----------------------------------------------------------------------
      function i2p(i,mps)

      real mps(1)

      i2p=0

      ic=1

      do while (i.le.mps(ic))
         i2p=mps(ic)
         ic=ic
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine ijk2pqr(ip,iq,ir,i,j,k,mps,mqs,mrs)

      real mps(1),mqs(1),mrs(1)

      ip=i2p(i,mps)
      iq=i2p(j,mqs)
      ir=i2p(k,mrs)

      return
      end
c-----------------------------------------------------------------------
      function ijk2proc(i,j,k,mps,mqs,mrs)

      real mps(1),mqs(1),mrs(1)

      call ijk2pqr(ip,iq,ir,i,j,k,mps,mqs,mrs)

      ipqr2id=(ip-1)+mp*(iq-1)+mp*mq*(ir-1)

      return
      end
c-----------------------------------------------------------------------
