c-----------------------------------------------------------------------
      ! Unit Tests
c-----------------------------------------------------------------------
      subroutine grammian_test

      include 'SIZE'
      include 'MOR'

      real vv(lt,ls)

      call gengram
      call readgram(vv,ls)

      s1=0.
      s2=0.

      do i=1,ls*ls
         s1=s1+(uu(i,1)-vv(i,1))**2
         s2=s2+uu(i,1)**2
      enddo

      efrob=sqrt(s1/s2)
      if (nio.eq.0) write (6,*) 'fnorm',efrob

      iexit=1
      if (efrob.lt.1) iexit=0

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine basis_test

      logical ifl2

      return
      end
c-----------------------------------------------------------------------
      subroutine initial_condition_test(ifl2)

      logical ifl2

      return
      end
c-----------------------------------------------------------------------
      subroutine initial_condition_test(ifl2)

      logical ifl2

      return
      end
c-----------------------------------------------------------------------
      subroutine a_operator_test(ifl2)

      logical ifl2

      return
      end
c-----------------------------------------------------------------------
      subroutine b_operator_test(ifl2)

      logical ifl2

      return
      end
c-----------------------------------------------------------------------
      subroutine c_operator_test(ifl2)

      logical ifl2

      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
