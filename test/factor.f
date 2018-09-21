c-----------------------------------------------------------------------
      program main

      parameter (n=2400)

      integer factors(n), parts(3), blocks(n)

      nb = 400

      call izero(blocks,n)

      call factor(factors,n)
c     call factor3(parts,factors)
      m=nint(rand(1)*2048)
      m=2
      write (6,*) 'n=',m
      call factor3(mp,mq,mr,m)
      stop

      do i=1,3
         write (6,*) i,'th partition:'
         write (6,*) parts(i), 'parts'
         call setblocks(blocks,parts(i),nb)
         do j=1,parts(i)
            write (6,*) blocks(j)
         enddo
      enddo

      iflag = -1
      i = 1

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
