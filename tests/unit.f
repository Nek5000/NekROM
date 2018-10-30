c-----------------------------------------------------------------------
      ! Unit Tests
c-----------------------------------------------------------------------
      include 'rom.f'
c-----------------------------------------------------------------------
      subroutine grammian_unit(iflag)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      logical iflag
      real vv(ls,ls)

      call rom_init

      ifl2=iflag

      iexit=0

      call gengram
      call readgram(vv,ls)

      s1=0.
      s2=0.
      s3=0.

      do j=1,ls
      do i=1,ls
         s1=s1+(uu(i,j)-uu(j,i))**2
         s2=s2+(uu(i,j)-vv(i,j))**2
         s3=s3+vv(i,j)**2
         if (nio.eq.0) write (6,*) 'gram',i,j,uu(i,j),vv(i,j)
      enddo
      enddo

      esym=sqrt(s1/s3)
      edif=sqrt(s2/s3)

      if (nio.eq.0) write (6,*) 'esym',esym,s1,s3
      if (nio.eq.0) write (6,*) 'edif',edif,s2,s3

      if (esym.gt.1e-14) iexit=iexit+1
      if (edif.gt.1e-14) iexit=iexit+2

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine eigenvector_unit(iflag)

      include 'SIZE'
      include 'MOR'

      logical iflag
      real evec2(ls,nb)

      iexit=0

      call readgram(uu,ls)
      call genevec(evec)
      call readevec(evec2,ls,nb)

      s1=0.
      s2=0.

      do j=1,nb
      do i=1,ls
         s1=s1+(evec(i,j)-evec2(i,j))**2
         s2=s2+evec2(i,j)**2
         if (nio.eq.0) write (6,*) 'evec',i,j,evec(i,j),evec2(i,j)
      enddo
      enddo

      edif=sqrt(s1/s2)
      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      if (edif.gt.1e-16) iexit=1

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine bases_unit(iflag)

      include 'SIZE'
      include 'MOR'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      logical iflag

      real ubb(lt,0:nb), vbb(lt,0:nb), wbb(lt,0:nb)
      real du(lt,0:nb), dv(lt,0:nb), dw(lt,0:nb)

      call rom_init

      ifl2=iflag

      n=lx1*ly1*lz1*nelt

      call readevec(evec,ls,nb)
      call genbases
      call readbases(ubb,vbb,wbb,nb)

      s1=0.
      s2=0.

      ! TODO use H10 norm if(.not.ifl2)

      do i=0,nb
         call opsub3(du(1,i),dv(1,i),dw(1,i),ub(1,i),vb(1,i),wb(1,i),
     $                                       ubb(1,i),vbb(1,i),wbb(1,i))
         s1=s1+op_glsc2_wt(du(1,i),dv(1,i),dw(1,i),
     $                     du(1,i),dv(1,i),dw(1,i),bm1)
         s2=s2+op_glsc2_wt(ubb(1,i),vbb(1,i),wbb(1,i),
     $                     ubb(1,i),vbb(1,i),wbb(1,i),bm1)
      enddo

      edif=sqrt(s1/s2)
      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      iexit=1
      if (edif.lt.1e-16) iexit=0

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine initial_condition_unit(iflag)

      include 'SIZE'
      include 'SOLN'
      include 'MOR'

      logical iflag
      real u0(0:nb)

      call rom_init

      ifl2=iflag

      call readbases(ub,vb,wb,nb)
      call makeic
      call readic(u0,nb+1)

      s1=0.
      s2=0.

      do i=0,nb
         write (6,*) u0(i),u(i,1),'ic'
      enddo

      do i=1,nb
         s1=s1+(u0(i)-u(i,1))**2
         s2=s2+u0(i)**2
      enddo

      edif=sqrt(s1/s2)
      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      iexit=1
      if (edif.lt.1e-16) iexit=0

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine a0_unit(iflag)

      include 'SIZE'
      include 'MOR'

      logical iflag

      real aa(0:nb,0:nb)

      call rom_init
      ifl2=iflag

      iexit=0

      call readbases(ub,vb,wb,nb)

      call makea
      call reada0(aa,(nb+1)**2)

      s1=0.
      s2=0.
      s3=0.

      do j=0,nb
      do i=0,nb
         s1=s1+(aa(i,j)-a0(i,j))**2
         s2=s2+(a0(i,j)-a0(j,i))**2
         s3=s3+aa(i,j)**2
      enddo
      enddo

      edif=sqrt(s1/s3)
      esym=sqrt(s2/s3)

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s3
      if (nio.eq.0) write (6,*) 'esym',esym,s2,s3

      if (edif.gt.1.e-16) iexit=iexit+1
      if (esym.gt.1.3e-15) iexit=iexit+2

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
      if (.not.ifl2.and.edia.gt.3.5e-15) iexit=iexit+4

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine b0_unit(iflag)

      include 'SIZE'
      include 'SOLN'
      include 'MOR'

      logical iflag

      real bb(0:nb,0:nb)

      call rom_init
      ifl2=iflag

      iexit=0

      call readbases(ub,vb,wb,nb)

      call makeb
      call readb0(bb,(nb+1)**2)

      s1=0.
      s2=0.
      s3=0.

      do j=0,nb
      do i=0,nb
         s1=s1+(bb(i,j)-b0(i,j))**2
         s2=s2+(b0(i,j)-b0(j,i))**2
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
         if (i.ne.j) s1=s1+bb(i,j)**2
         s2=s2+bb(i,j)**2
      enddo
      enddo

      edia=sqrt(s1/s2)

      if (ifl2.and.edia.gt.2.e-15) iexit=iexit+4

      if (nio.eq.0) write (6,*) 'edia',edia,s1,s2

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine c0_unit(iflag)

      include 'SIZE'
      include 'SOLN'
      include 'MOR'

      logical iflag

      real cc(0:nb,0:nb,0:nb)

      call rom_init
      ifl2=iflag

      iexit=0

      call readbases(ub,vb,wb,nb)

      call makec
      call readc0(cc,(nb+1)**3)

      s1=0.
      s2=0.
      s3=0.

      do k=0,nb
      do j=0,nb
      do i=0,nb
         s1=s1+(cc(i,j,k)-c0(i,j,k))**2
         s2=s2+(c0(i,j,k)+c0(j,i,k))**2
         s3=s3+cc(i,j,k)**2
      enddo
      enddo
      enddo

      edif=sqrt(s1/s3)
      if (edif.gt.1.e-16) iexit=iexit+1
      if (nio.eq.0) write (6,*) 'edif',edif,s1,s3

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
