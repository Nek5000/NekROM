c-----------------------------------------------------------------------
      ! Unit Tests
c-----------------------------------------------------------------------
      include 'read.f'
      include 'pod.f'
c-----------------------------------------------------------------------
      subroutine grammian_test

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real vv(lt,ls)

      call gengram
      call readgram(vv,ls)

      s1=0.
      s2=0.

      do j=1,ls
      do i=1,ls
         s1=s1+0.25*(uu(i,j)-uu(j,i))**2
         s2=s2+uu(i,j)**2
      enddo
      enddo

      esym=sqrt(s1/s2)
      if (nio.eq.0) write (6,*) 'esym',esym

      s1=0.
      s2=0.

      do i=1,ls*ls
         s1=s1+(uu(i,1)-vv(i,1))**2
         s2=s2+uu(i,1)**2
      enddo

      edif=sqrt(s1/s2)
      if (nio.eq.0) write (6,*) 'edif',edif

      iexit=1
      if (edif.lt.1e-16.and.esym.lt.2-e15) iexit=0

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
