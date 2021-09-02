c-----------------------------------------------------------------------
      ! Unit Tests
c-----------------------------------------------------------------------
      subroutine gramian_l2_unit
      call gramian_unit(.true.)
      return
      end
c-----------------------------------------------------------------------
      subroutine gramian_h10_unit
      call gramian_unit(.false.)
      return
      end
c-----------------------------------------------------------------------
      subroutine gramian_unit(iflag)

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      common /scrtest/ wk(ls*ls)

      parameter (lt=lx1*ly1*lz1*lelt)

      logical iflag
      real vv(ls,ls)

      param(171) = 1.
      if (iflag) param(171) = 0.
      param(172) = 1.
      param(173) = 0.

      call rom_setup
      call gengram(ug,us0,ns,ldim,ips)

      call read_serial(vv,ls*ls,'tops/gu ',wk,nid)

      iexit=0

      s1=0.
      s2=0.
      s3=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do j=1,ls
      do i=1,ls
         s1=s1+(ug(i,j)-ug(j,i))**2
         s2=s2+(ug(i,j)-vv(i,j))**2
         s3=s3+vv(i,j)**2
         if (nio.eq.0) write (6,*) 'gram',i,j,ug(i,j),vv(i,j)
      enddo
      enddo

      esym=sqrt(s1/s3)
      edif=sqrt(s2/s3)

      if (nio.eq.0) write (6,*) 'esym',esym,s1,s3
      if (nio.eq.0) write (6,*) 'edif',edif,s2,s3

      if (ips.eq.'L2 '.and.esym.gt.1e-16) iexit=iexit+1
      if (ips.eq.'H10'.and.esym.gt.1e-14) iexit=iexit+1
      if (edif.gt.5.e-15) iexit=iexit+2

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine ic_l2_unit
      call ic_unit(.true.)
      return
      end
c-----------------------------------------------------------------------
      subroutine ic_h10_unit
      call ic_unit(.false.)
      return
      end
c-----------------------------------------------------------------------
      subroutine ic_unit(iflag)

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MOR'

      common /scrtest/ wk(lb+1),u0(0:lb)

      logical iflag

      param(171) = 1.
      if (iflag) param(171) = 0.
      param(172) = 1.
      param(173) = 0.

      call rom_setup
      call read_serial(u0,nb+1,'tops/u ',wk,nid)

      s1=0.
      s2=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do i=0,nb
         write (6,*) i,u(i),u0(i),'ic'
      enddo

      do i=1,nb
         s1=s1+(u0(i)-u(i))**2
         s2=s2+u0(i)**2
      enddo

      edif=sqrt(s1/s2)
      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      iexit=1
      if (edif.lt.8.e-13) iexit=0

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine a0_l2_unit
      call a0_unit(.true.)
      return
      end
c-----------------------------------------------------------------------
      subroutine a0_h10_unit
      call a0_unit(.false.)
      return
      end
c-----------------------------------------------------------------------
      subroutine a0_unit(iflag)

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      logical iflag

      param(171) = 1.
      if (iflag) param(171) = 0.
      param(172) = 1.
      param(173) = 0.

      call rom_setup
      call a0_unit_helper(au0)

      return
      end
c-----------------------------------------------------------------------
      subroutine a0_unit_helper(a0)

      include 'SIZE'
      include 'MOR'

      common /scrtest/ wk(lb+1,lb+1),aa(0:lb,0:lb)

      real a0(0:nb,0:nb)
      call read_serial(aa,(nb+1)**2,'tops/au ',wk,nid)

      iexit=0

      s1=0.
      s2=0.
      s3=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do j=0,nb
      do i=0,nb
         if (nio.eq.0) write (6,*) i,j,a0(i,j),aa(i,j),'a'
         s1=s1+(aa(i,j)-a0(i,j))**2
         s2=s2+(a0(i,j)-a0(j,i))**2
         s3=s3+aa(i,j)**2
      enddo
      enddo

      edif=sqrt(s1/s3)
      esym=sqrt(s2/s3)

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s3
      if (nio.eq.0) write (6,*) 'esym',esym,s2,s3

      if (edif.gt.3.e-13) iexit=iexit+1
      if (esym.gt.1.e-15) iexit=iexit+2

      s1=0.
      s2=0.

      do j=1,nb
      do i=1,nb
         if (i.ne.j) s1=s1+a0(i,j)**2
         s2=s2+a0(i,j)**2
      enddo
      enddo

      edia=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'edia',edia,s1,s2
      if (ips.eq.'H10'.and.edia.gt.1.e-11) iexit=iexit+4

      s1=0.

      do i=1,nb
         s1=s1+(a0(i,i)-1.)**2
      enddo

      euni=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'euni',euni,s1,s2
      if (ips.eq.'H10'.and.euni.gt.1.e-13) iexit=iexit+8

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine b0_l2_unit
      call b0_unit(.true.)
      return
      end
c-----------------------------------------------------------------------
      subroutine b0_h10_unit
      call b0_unit(.false.)
      return
      end
c-----------------------------------------------------------------------
      subroutine b0_unit(iflag)

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      logical iflag

      param(171) = 1.
      if (iflag) param(171) = 0.
      param(172) = 1.
      param(173) = 0.

      call rom_setup
      call b0_unit_helper(bu0)

      return
      end
c-----------------------------------------------------------------------
      subroutine b0_unit_helper(b0)

      include 'SIZE'
      include 'SOLN'
      include 'MOR'

      common /scrtest/ wk(lb+1,lb+1),bb(0:lb,0:lb)

      real b0(0:nb,0:nb)

      call read_serial(bb,(nb+1)**2,'tops/bu ',wk,nid)

      iexit=0

      s1=0.
      s2=0.
      s3=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do j=0,nb
      do i=0,nb
         if (nio.eq.0) write (6,*) i,j,b0(i,j),bb(i,j),'b'
         s1=s1+(bb(i,j)-b0(i,j))**2
         s2=s2+(b0(i,j)-b0(j,i))**2
         s3=s3+bb(i,j)**2
      enddo
      enddo

      edif=sqrt(s1/s3)
      esym=sqrt(s2/s3)

      if (edif.gt.7.e-16) iexit=iexit+1
      if (esym.gt.1.e-16) iexit=iexit+2

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s3
      if (nio.eq.0) write (6,*) 'esym',esym,s2,s3

      s1=0.
      s2=0.

      do j=1,nb
      do i=1,nb
         if (i.ne.j) s1=s1+b0(i,j)**2
         s2=s2+bb(i,j)**2
      enddo
      enddo

      edia=sqrt(s1/s2)

      if (ips.eq.'L2 '.and.edia.gt.1.e-9) iexit=iexit+4
      if (nio.eq.0) write (6,*) 'edia',edia,s1,s2

      s1=0.

      do i=1,nb
         s1=s1+(b0(i,i)-1.)**2
      enddo

      euni=sqrt(s1/s2)

      if (ips.eq.'L2 '.and.euni.gt.1.e-14) iexit=iexit+8
      if (nio.eq.0) write (6,*) 'euni',euni,s1,s2

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine c0_l2_unit
      call c0_unit(.true.)
      return
      end
c-----------------------------------------------------------------------
      subroutine c0_h10_unit
      call c0_unit(.false.)
      return
      end
c-----------------------------------------------------------------------
      subroutine c0_unit(iflag)

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      logical iflag

      param(171) = 1.
      if (iflag) param(171) = 0.
      param(172) = 1.
      param(173) = 0.

      call rom_setup
      call c0_unit_helper(cul)

      return
      end
c-----------------------------------------------------------------------
      subroutine c0_unit_helper(c,iflag)

      include 'SIZE'
      include 'SOLN'
      include 'MOR'

      common /scrtest/ wk(lb+1,lb+1,lb+1),cc(lb,0:lb,0:lb)

      real c(nb,0:nb,0:nb)

      call read_serial(cc,nb*(nb+1)**2,'tops/cu ',wk,nid)

      iexit=0

      s1=0.
      s2=0.
      s3=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do k=0,nb
      do j=0,nb
      do i=1,nb
         s1=s1+(cc(i,j,k)-c(i,j,k))**2
         s3=s3+cc(i,j,k)**2
         write (6,*) i,j,k,c(i,j,k),cc(i,j,k),'c'
      enddo
      enddo
      enddo

      do k=1,nb
      do j=1,nb
      do i=1,nb
         s2=s2+(c(i,j,k)+c(j,i,k))**2
      enddo
      enddo
      enddo

      edif=sqrt(s1/s3)
      if (edif.gt.2.e-10) iexit=iexit+1
      if (nio.eq.0) write (6,*) 'edif',edif,s1,s3

      eskew=sqrt(s2/s3)
c     if (eskew.gt.1.e-16) iexit=iexit+2
      if (nio.eq.0) write (6,*) 'eskew',eskew,s2,s3

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
