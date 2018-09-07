c-----------------------------------------------------------------------
      program main

      parameter (n=2400)

      integer factors(n), parts(3), blocks(n)

      nb = 300

      call izero(blocks,n)

      call factor(factors,n)
      call partition3(parts,factors)

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

c     do while (iflag.ne.0)
c        iflag = factors(i)
c        write (6,*) iflag
c        i=i+1
c     enddo

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
      subroutine partition3(parts,facs)

      integer parts(3), facs(1), n

      parts(1) = 1
      parts(2) = 1
      parts(3) = 1

      n=1
      i=0

      do while (facs(i+1).ne.0)
         m=mod(i,3)+1
         parts(m) = parts(m) * facs(i+1)
         i=i+1
      enddo

      write (6,*) parts(1),parts(2),parts(3)

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

      real c(nb),cc(1),u(nb+1)

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

      call gop(c,work,'+  ',nb)

      return
      end
c-----------------------------------------------------------------------
