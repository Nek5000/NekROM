c-----------------------------------------------------------------------
      program main

      parameter (n=2400)

      integer factors(n), parts(3)

      call factor(factors,n)

      call partition3(parts,factors)

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

      do i=1,n
         factors(i) = 0
      enddo

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

      nb=100

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
