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
      subroutine factor(factors,m)

      integer factors(m)

      n=m

      call izero(factors,n)

      j=1
      i=2

    1 continue

      do while (n.ne.1)
         if (n.eq.(n/i)*i) then
            factors(j)=i
            n = n / i
            write (6,*) 'i,n',i,n
            if (n.eq.1) then
               goto 2
            else
               j=j+1
            endif
         else
            i=i+1
         endif
         goto 1
      enddo

    2 continue

      return
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
      subroutine setparts(mps,mqs,mrs,mp,mq,mr,nb)

      integer mps(mp),mqs(mq),mrs(mr)

      mb=nb+1

      do i=0,mp-1
         mps(i+1)=nb/mp+max(mod(nb,mp)-i,0)/max(mod(nb,mp)-i,1)
         mps(i+1)=mps(i+1)+mps(max(i,1))*max(i,0)/max(i,1)
         write (6,*) 'mp',i+1,mps(i+1)
      enddo

      do i=0,mq-1
         mqs(i+1)=mb/mq+max(mod(mb,mq)-i,0)/max(mod(mb,mq)-i,1)
         mqs(i+1)=mqs(i+1)+mqs(max(i,1))*max(i,0)/max(i,1)
         write (6,*) 'mq',i+1,mqs(i+1)
      enddo

      do i=0,mr-1
         mrs(i+1)=mb/mr+max(mod(mb,mr)-i,0)/max(mod(mb,mr)-i,1)
         mrs(i+1)=mrs(i+1)+mrs(max(i,1))*max(i,0)/max(i,1)
         write (6,*) 'mr',i+1,mrs(i+1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setblocks(blocks,nparts,nb)

      integer blocks(nparts)

      min = nb / nparts
      mb = nb - min * nparts

      do i=1,nparts
         blocks(i) = min
         if (i.lt.mb) blocks(i) = blocks(i) + 1
      enddo
        
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
