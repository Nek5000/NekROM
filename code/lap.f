c-----------------------------------------------------------------------
      subroutine setops_laplace

      include 'SIZE'
      include 'MOR'
      include 'TSTEP'

      jfield=ifield
      ifield=1
      call seta_laplace(at,at0,'ops/at ')
      call setb_laplace(bt,bt0,'ops/bt ')

      return
      end
c-----------------------------------------------------------------------
      subroutine seta_laplace(a,a0,fname)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrseta/ wk1(lt)

      real a0(0:nb,0:nb),a(nb,nb)

      character*128 fname

      if (nio.eq.0) write (6,*) 'inside seta'

      n=lx1*ly1*lz1*nelt

      nio=-1
      do j=0,nb ! form the stiffness matrix of the discrete problem
      do i=0,nb
         a0(i,j)=h10sip_vd(tb(1,i),tb(1,j))
      enddo
      enddo
      nio=nid

      do j=1,nb ! form the stiffness matrix of the discrete problem
      do i=1,nb
         a(i,j)=a0(i,j)
      enddo
      enddo

      if (nio.eq.0) write (6,*) 'exiting seta'

      return
      end
c-----------------------------------------------------------------------
      subroutine setb_laplace(b,b0,fname)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrseta/ wk1(lt)

      real b0(0:nb,0:nb),b(nb,nb)

      character*128 fname

      if (nio.eq.0) write (6,*) 'inside setb'

      n=lx1*ly1*lz1*nelt

      nio=-1
      do j=0,nb ! form the stiffness matrix of the discrete problem
      do i=0,nb
         b0(i,j)=wl2sip(tb(1,i),tb(1,j))
      enddo
      enddo
      nio=nid

      do j=1,nb ! form the stiffness matrix of the discrete problem
      do i=1,nb
         b(i,j)=b0(i,j)
      enddo
      enddo

      if (nio.eq.0) write (6,*) 'exiting setb'

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_solve_laplace

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'
      include 'AVG'

      integer icalld
      save    icalld
      data    icalld /0/

      real rhs(0:nb)

      logical ifmult

      common /romup/ rom_time

      stime=dnekclock()

      call opcopy(uic,vic,wic,vx,vy,vz)

      call rom_init_params
      ifrom(1)=.false.
      ifpod(1)=.false.
      n=lx1*ly1*lz1*nelt
      call rzero(zeros,n)
      call rone(ones,n)

      call setgram
      call setevec

      call setbases

      nbs=nint(sqrt(nb*1.0))
      if (nb.ne.nbs*nbs) then
         call exitti('nb not a square of an integer$',nb)
      endif

      l=1
      do k=1,nbs
      do j=1,nbs
         do i=1,n
            x=xm1(i,1,1,1)
            y=ym1(i,1,1,1)
            tb(i,l)=sin(j*pi*x)*sin(k*pi*y)
         enddo
         l=l+1
      enddo
      enddo

c      do i=1,ns
c         call copy(tb(1,i),ts(1,i),n)
c      enddo

      call setops_laplace
      
      if (ifpod(1)) call pv2k(uk,us,ub,vb,wb)
      if (ifpod(2)) call ps2k(tk,ts,tb)

c      call asnap

      if (ifdumpops) call dump_all

      ad_step = istep
      jfield=ifield
      ifield=1

      rhs(0)=1.
      call setr_laplace(rhs(1),icount)

      do i=1,nb
         write (6,*) rhs(i),'rhs'
      enddo

      do j=1,nb
      do i=1,nb
         write (6,*) at(i,j),'flu'
      enddo
      enddo

      call lu(at,nb,nb,irv,icv)
      call solve(rhs(1),at,1,nb,nb,irv,icv)

      do i=1,nb
         write (6,*) rhs(i),'u'
      enddo

      rhs(0)=0.

      call copy(ut,rhs,nb+1)
      call recont(t,rhs)
      call outpost(vx,vy,vz,pr,t,'sol')

      ifield=jfield

      dtime=dnekclock()-stime
      rom_time=rom_time+dtime

      return
      end
c-----------------------------------------------------------------------
      subroutine setr_laplace(rhs)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrrhs/ tmp1(0:nb),tmp2(0:nb),s(lt)

      real rhs(nb)
      real tmp(nb)

      do i=1,nb
         call savg(s_bar,a_surf,tb(1,i),2,'f  ')
         rhs(i)=s_bar*a_surf
      enddo
      
      call ps2b(tmp,g,tb)
      call mxm(bt,nb,tmp,nb,rhs,1)

      return
      end
c-----------------------------------------------------------------------
