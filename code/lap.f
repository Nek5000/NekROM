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
      ifpod(2)=.true.
      ifrom(2)=.true.
      call rom_init_fields

      n=lx1*ly1*lz1*nelt
      call rzero(zeros,n)
      call rone(ones,n)

      call setgram
c     call setevec

      call setbases

c     nbs=nint(sqrt(nb*1.0))
c     if (nb.ne.nbs*nbs) then
c        call exitti('nb not a square of an integer$',nb)
c     endif

c     l=1
c     do k=1,nbs
c     do j=1,nbs
c        do i=1,n
c           x=xm1(i,1,1,1)
c           y=ym1(i,1,1,1)
c           tb(i,l)=sin(j*pi*x)*sin(k*pi*y)
c        enddo
c        l=l+1
c     enddo
c     enddo

c     call snorm(tb)

      do i=1,ns
         call copy(tb(1,i),ts(1,i),n)
      enddo

      call setops_laplace
      
      if (ifpod(1)) call pv2k(uk,us,ub,vb,wb)
      if (ifpod(2)) call ps2k(tk,ts,tb)

c      call asnap

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ')
     $   call dump_all

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
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrrhs/ tmp1(0:nb),tmp2(0:nb),s(lt)

      real rhs(nb)
      real tmp(0:nb)

      n=lx1*ly1*lz1*nelv

      one=1.
      pi=4.*atan(one)

      do i=1,n
         s(i)=sin(pi*xm1(i,1,1,1))*sin(pi*ym1(i,1,1,1))
      enddo

      do i=1,nb
c        call savg(s_bar,a_surf,tb(1,i),2,'f  ')
c        rhs(i)=s_bar*a_surf
         rhs(i)=glsc3(tb(1,i),s,bm1,n)
      enddo

c     call rzero(rhs,nb)

c     do i=1,n
c        x=xm1(i,1,1,1)
c        y=ym1(i,1,1,1)
c        g(i,1,1,1) = sin(2*pi*x)*sin(2*pi*y)
c     enddo
c     
c     call ps2b(tmp,g,tb)
c     call mxm(bt,nb,tmp(1),nb,rhs,1)

      return
      end
c-----------------------------------------------------------------------
