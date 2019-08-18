c-----------------------------------------------------------------------
      ! Unit Tests
c-----------------------------------------------------------------------
      subroutine grammian_unit(iflag)

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
      call read_serial(vv,ls*ls,'ops/gu ',wk,nid)

      iexit=0

      s1=0.
      s2=0.
      s3=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do j=1,ls
      do i=1,ls
         s1=s1+(ug(i,j,1)-ug(j,i,1))**2
         s2=s2+(ug(i,j,1)-vv(i,j))**2
         s3=s3+vv(i,j)**2
         if (nio.eq.0) write (6,*) 'gram',i,j,ug(i,j,1),vv(i,j)
      enddo
      enddo

      esym=sqrt(s1/s3)
      edif=sqrt(s2/s3)

      if (nio.eq.0) write (6,*) 'esym',esym,s1,s3
      if (nio.eq.0) write (6,*) 'edif',edif,s2,s3

      if (ips.eq.'L2 '.and.esym.gt.1e-16) iexit=iexit+1
      if (ips.eq.'H10'.and.esym.gt.1e-14) iexit=iexit+1
      if (edif.gt.1e-16) iexit=iexit+2

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine initial_condition_unit(iflag)

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
      call read_serial(u0,nb+1,'ops/u ',wk,nid)

      s1=0.
      s2=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do i=0,nb
         write (6,*) i,u(i,1),u0(i),'ic'
      enddo

      do i=1,nb
         s1=s1+(u0(i)-u(i,1))**2
         s2=s2+u0(i)**2
      enddo

      edif=sqrt(s1/s2)
      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      iexit=1
      if (edif.lt.1.e-13) iexit=0

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine a0_unit(iflag)

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      common /scrtest/ wk(lb+1,lb+1)

      logical iflag

      real aa(0:nb,0:nb)

      param(171) = 1.
      if (iflag) param(171) = 0.
      param(172) = 1.
      param(173) = 0.

      call rom_setup
      call read_serial(aa,(nb+1)**2,'ops/au ',wk,nid)

      iexit=0

      s1=0.
      s2=0.
      s3=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do j=0,nb
      do i=0,nb
         if (nio.eq.0) write (6,*) i,j,au0(i,j),aa(i,j),'a'
         s1=s1+(aa(i,j)-au0(i,j))**2
         s2=s2+(au0(i,j)-au0(j,i))**2
         s3=s3+aa(i,j)**2
      enddo
      enddo

      edif=sqrt(s1/s3)
      esym=sqrt(s2/s3)

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s3
      if (nio.eq.0) write (6,*) 'esym',esym,s2,s3

      if (edif.gt.1.e-16) iexit=iexit+1
      if (esym.gt.1.e-15) iexit=iexit+2

      s1=0.
      s2=0.

      do j=1,nb
      do i=1,nb
         if (i.ne.j) s1=s1+au0(i,j)**2
         s2=s2+au0(i,j)**2
      enddo
      enddo

      edia=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'edia',edia,s1,s2
      if (ips.eq.'H10'.and.edia.gt.1.e-11) iexit=iexit+4

      s1=0.

      do i=1,nb
         s1=s1+(au0(i,i)-1.)**2
      enddo

      euni=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'euni',euni,s1,s2
      if (ips.eq.'H10'.and.euni.gt.1.e-13) iexit=iexit+8

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine b0_unit(iflag)

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MOR'

      common /scrtest/ wk(lb+1,lb+1)

      logical iflag

      real bb(0:nb,0:nb)

      param(171) = 1.
      if (iflag) param(171) = 0.
      param(172) = 1.
      param(173) = 0.

      call rom_setup
      call read_serial(bb,(nb+1)**2,'ops/bu ',wk,nid)

      iexit=0

      s1=0.
      s2=0.
      s3=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do j=0,nb
      do i=0,nb
         if (nio.eq.0) write (6,*) i,j,bu0(i,j),bb(i,j),'b'
         s1=s1+(bb(i,j)-bu0(i,j))**2
         s2=s2+(bu0(i,j)-bu0(j,i))**2
         s3=s3+bb(i,j)**2
      enddo
      enddo

      edif=sqrt(s1/s3)
      esym=sqrt(s2/s3)

      if (edif.gt.1.e-16) iexit=iexit+1
      if (esym.gt.1.e-16) iexit=iexit+2

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s3
      if (nio.eq.0) write (6,*) 'esym',esym,s2,s3

      s1=0.
      s2=0.

      do j=1,nb
      do i=1,nb
         if (i.ne.j) s1=s1+bu0(i,j)**2
         s2=s2+bb(i,j)**2
      enddo
      enddo

      edia=sqrt(s1/s2)

      if (ips.eq.'L2 '.and.edia.gt.1.e-9) iexit=iexit+4
      if (nio.eq.0) write (6,*) 'edia',edia,s1,s2

      s1=0.

      do i=1,nb
         s1=s1+(bu0(i,i)-1.)**2
      enddo

      euni=sqrt(s1/s2)

      if (ips.eq.'L2 '.and.euni.gt.1.e-14) iexit=iexit+8
      if (nio.eq.0) write (6,*) 'euni',euni,s1,s2

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine c0_unit(iflag)

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MOR'

      common /scrtest/ wk(lb+1,lb+1,lb+1)

      logical iflag

      real cc(lcloc), cglob(nb,nb+1,nb+1)

      param(171) = 1.
      if (iflag) param(171) = 0.
      param(172) = 1.
      param(173) = 0.

      call rom_setup
      call read_serial(cc,nb*(nb+1)**2,'ops/cu ',wk,nid)

      iexit=0

      s1=0.
      s2=0.
      s3=0.

      call rzero(cglob,nb*(nb+1)**2)

      if (nio.eq.0) write (6,*) 'live | data'

      do jc=1,nb*(nb+1)**2
         i=icul(1,jc)
         j=icul(2,jc)
         k=icul(3,jc)

         cglob(i,j,k)=cglob(i,j,k)+cul(jc)
         cglob(k,j,i)=cglob(k,j,i)-cc(jc)

         s1=s1+(cc(jc)-cul(jc))**2
         s3=s3+cc(jc)**2
         write (6,*) i,j,k,cul(jc),cc(jc),'c'
      enddo

      do k=1,nb
      do j=1,nb
      do i=1,nb
         s2=s2+(cglob(i,j,k)+cglob(j,i,k))**2
      enddo
      enddo
      enddo

      edif=sqrt(s1/s3)
      if (edif.gt.1.e-16) iexit=iexit+1
      if (nio.eq.0) write (6,*) 'edif',edif,s1,s3

      eskew=sqrt(s2/s3)
c     if (eskew.gt.1.e-16) iexit=iexit+2
      if (nio.eq.0) write (6,*) 'eskew',eskew,s2,s3

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
