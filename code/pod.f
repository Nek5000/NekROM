c-----------------------------------------------------------------------
      subroutine setbases_heat

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      n=lx1*ly1*lz1*nelv

      do ib=0,nb
         call rzero(ub(1,ib),n)
         call rzero(vb(1,ib),n)
         call rzero(wb(1,ib),n)
      enddo

      one=1.
      pi=4.*atan(one)

      do i=1,n
         x=xm1(i,1,1,1)
         y=ym1(i,1,1,1)
         ub(i,0)=1.
         k=2
         ub(i,1)=sin(k*pi*x)*sin(k*pi*y)
         if (nb.ge.2) ub(i,2)=cos(k*pi*x)*sin(k*pi*y)
         if (nb.ge.3) ub(i,3)=sin(k*pi*x)*cos(k*pi*y)
         if (nb.ge.4) ub(i,4)=cos(k*pi*x)*cos(k*pi*y)
         if (nb.ge.4) vb(i,0)=0.5
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setbases

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real u0(lt,3)
      common /scrk3/ t4(lt),t5(lt),t6(lt)
      common /scrk4/ h1(lt),h2(lt),bwm1(lt)
      common /scrk5/ au(lt),av(lt),aw(lt)

      if (nio.eq.0) write (6,*) 'inside setbases'

      if (ifread) then
         call loadbases(ub,vb,wb,nb)
      else
         n=lx1*ly1*lz1*nelt

         call rone(h1,n)
         call rzero(h2,n)
         call col3(bwm1,bm1,wm1,n)

         one = 1.
         zero= 0.

         ns = ls ! REQUIRED: get_saved_fields overwrites ns argument
         call opcopy(u0(1,1),u0(1,2),u0(1,3),ub(1,0),vb(1,0),wb(1,0))

         ! ub, vb, wb, are the modes
         call dgemm( 'N','N',n,nb,ls,one,ust,lt,evec,ls,zero,ub(1,1),lt)
         call dgemm( 'N','N',n,nb,ls,one,vst,lt,evec,ls,zero,vb(1,1),lt)
         if (ldim.eq.3)
     $   call dgemm( 'N','N',n,nb,ls,one,wst,lt,evec,ls,zero,wb(1,1),lt)

         call scale_bases
      endif

      if (ifdrago) then
         n=lx1*ly1*lz1*nelt
         call copy(h1,vdiff,n)
         call rzero(h2,n)
         do i=0,nb
            call lap2d(au,ub(1,i))
            call lap2d(av,vb(1,i))
            if (ldim.eq.3) call lap2d(aw,wb(1,i))
            call opcmult(au,av,aw,param(2))
            call comp_pdrag(fd1(1,i),ub(1,i),vb(1,i),wb(1,i))
            call comp_pdrag(fd3(1,i),au,av,aw)
         enddo
      endif

      if (nio.eq.0) write (6,*) 'exiting setbases'

      return
      end
c-----------------------------------------------------------------------
      subroutine proj2bases(coef,ux,uy,uz)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real coef(0:nb),ux(lt),uy(lt),uz(lt)

      if (ifl2) then
         call wl2proj(coef,ux,uy,uz)
      else
         call h10proj(coef,ux,uy,uz)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine h10proj(coef,ux,uy,uz)

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ux(lt),uy(lt),uz(lt)

      common /scrk3/ t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)
      common /scrk4/ h1(lt),h2(lt)

      real coef(0:nb)

      if (nio.eq.0) write (6,*) 'inside h10proj'

      n=lx1*ly1*lz1*nelt

      call rone(h1,n)
      call rzero(h2,n)
      call opsub3(t1,t2,t3,ux,uy,uz,ub,vb,wb)

      coef(0) = 1.
      if (nio.eq.0) write (6,1) coef(0),coef(0),1.

      do i=1,nb
         call axhelm(t4,ub(1,i),h1,h2,1,1)
         call axhelm(t5,vb(1,i),h1,h2,1,1)

         ww = glsc2(t4,ub(1,i),n)+glsc2(t5,vb(1,i),n)
         vv = glsc2(t4,t1,n)+glsc2(t5,t2,n)

         if (ldim.eq.3) then
            call axhelm(t6,wb(1,i),h1,h2,1,1)
            ww = ww + glsc2(t6,wb(1,i),n)
            vv = vv + glsc2(t6,t3,n)
         endif

         coef(i) = vv/ww
         if (nio.eq.0) write (6,1) coef(i),vv,ww
      enddo

      if (nio.eq.0) write (6,*) 'exiting h10proj'

    1 format(' h10coef',1p3e16.8)

      return
      end
c-----------------------------------------------------------------------
      subroutine wl2proj(coef,ux,uy,uz)

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ux(lt),uy(lt),uz(lt)

      common /scrk3/ t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)
      common /scrk4/ h1(lt),h2(lt),bwm1(lt)

      real coef(0:nb)

      if (nio.eq.0) write (6,*) 'inside wl2proj'

      n=lx1*ly1*lz1*nelt

      call col3(bwm1,bm1,wm1,n)
      call opsub3(t1,t2,t3,ux,uy,uz,ub,vb,wb)

      coef(0) = 1.

      if (nio.eq.0) write (6,1) coef(0),coef(0),1.0

      if (ifvort) then
         do i=1,nb
            ww = glsc3(ub(1,i),ub(1,i),bwm1,n)
            vv = glsc3(ub(1,i),t1,bwm1,n)

            coef(i) = vv/ww
            if (nio.eq.0) write (6,1) coef(i),vv,ww
         enddo
      else
         do i=1,nb
            ww = op_glsc2_wt(
     $         ub(1,i),vb(1,i),wb(1,i),ub(1,i),vb(1,i),wb(1,i),bwm1)
            vv = op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),t1,t2,t3,bwm1)

            coef(i) = vv/ww
            if (nio.eq.0) write (6,1) coef(i),vv,ww
         enddo
      endif

      if (nio.eq.0) write (6,*) 'exiting wl2proj'

    1 format(' wl2coef',1p3e16.8)

      return
      end
c-----------------------------------------------------------------------
      function vecprod(t1,t2,t3,t4,t5,t6,h1,h2,space)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)
      real h1(lt),h2(lt)

      character*3 space

      if (space.eq.'L2 ') then
         vecprod=wl2prod(t1,t2,t3,t4,t5,t6,h1,h2)
      else if (space.eq.'H10') then
         vecprod=h10prod(t1,t2,t3,t4,t5,t6,h1,h2)
      else
         call exitti('did not provide supported inner product space$')
      endif

      return
      end
c-----------------------------------------------------------------------
      function h10prod(t1,t2,t3,t4,t5,t6,h1,h2)

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)
      real h1(lt),h2(lt)

      common /scrk3/ t7(lt),t8(lt),t9(lt)

      if (nio.eq.0) write (6,*) 'inside h10prod'

      n=lx1*ly1*lz1*nelt

      call axhelm(t7,t1,h1,h2,1,1)
      h10prod=glsc2(t7,t4,n)

      if (.not.ifvort) then
         call axhelm(t8,t2,h1,h2,1,1)
         h10prod=h10prod+glsc2(t8,t5,n)

         if (ldim.eq.3) then
            call axhelm(t9,t3,h1,h2,1,1)
            h10prod = h10prod + glsc2(t9,t6,n)
         endif
      endif

      if (nio.eq.0) write (6,*) 'exiting h10prod'

      return
      end
c-----------------------------------------------------------------------
      function wl2prod(t1,t2,t3,t4,t5,t6,h1,h2)

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt)
      real h1(lt),h2(lt)

      common /scrk3/ bwm1(lt),t8(lt),t9(lt)

      if (nio.eq.0) write (6,*) 'inside wl2prod'

      n=lx1*ly1*lz1*nelt

      call col3(bwm1,bm1,wm1,n)

      if (ifvort) then
         wl2prod = glsc3(t1,t4,bwm1,n)
      else
         wl2prod = op_glsc2_wt(t1,t2,t3,t4,t5,t6,bwm1)
      endif

      if (nio.eq.0) write (6,*) 'exiting wl2prod'

      return
      end
c-----------------------------------------------------------------------
      subroutine gengram(gram,s,ms,mdim)

      include 'SIZE'
      include 'MOR'

      real s(1)
      real gram(1)

      if (ifl2) then
         call gengraml2(gram,s,ms,mdim)
      else
         call gengramh10(gram,s,ms,mdim)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine gengramh10(gram,s,ms,mdim)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real gram(ms,ms)
      real s(lt,mdim,ms)

      common /scrgram/ uw(lt),vw(lt),ww(lt),h1(lt),h2(lt)

      if (nio.eq.0) write (6,*) 'inside gengramh10'

      n  = lx1*ly1*lz1*nelt

      call rone (h1,n)
      call rzero(h2,n)

      do j=1,ms ! Form the Gramian, U=U_K^T A U_K using H^1_0 Norm
         call axhelm(uw,s(1,1,j),h1,h2,1,1)
         if (mdim.ge.2) call axhelm(vw,s(1,2,j),h1,h2,1,1)
         if (mdim.ge.3) call axhelm(ww,s(1,3,j),h1,h2,1,1)
         do i=1,ms
            gram(i,j) = glsc2(s(1,1,i),uw,n)
            if (mdim.ge.2) gram(i,j)=gram(i,j)+glsc2(s(1,2,i),vw,n)
            if (mdim.eq.3) gram(i,j)=gram(i,j)+glsc2(s(1,3,i),ww,n)
            if (nio.eq.0) write (99,*) gram(i,j)
         enddo
         if (nio.eq.0) write(6,1) j,gram(1,j)
      enddo

      if (nio.eq.0) write (6,*) 'exiting gengramh10'
    1 format (' gram',i5,1p1e16.6)

      return
      end
c-----------------------------------------------------------------------
      subroutine gengraml2(gram,s,ms,mdim)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrgram/ uw(lt),vw(lt),ww(lt),h1(lt),h2(lt)
      real gram(ms,ms)
      real s(lt,mdim,ms)

      if (nio.eq.0) write (6,*) 'inside gengraml2'

      n=lx1*ly1*lz1*nelv

      do j=1,ms ! Form the Gramian, U=U_K^T A U_K using L2 Norm
      do i=1,ms
         gram(i,j)=glsc3(s(1,1,i),s(1,1,j),bm1,n)
         if (mdim.ge.2)
     $      gram(i,j)=gram(i,j)+glsc3(s(1,2,i),s(1,2,j),bm1,n)
         if (mdim.ge.3)
     $      gram(i,j)=gram(i,j)+glsc3(s(1,3,i),s(1,3,j),bm1,n)
      enddo
         if (nio.eq.0) write (6,1) j,gram(1,j)
      enddo

      if (nio.eq.0) write (6,*) 'exiting gengraml2'

    1 format (' gram',i5,' ',1p1e16.6)

      return
      end
c-----------------------------------------------------------------------
      subroutine genevec

      !!! does not work if ns.lt.ls !!!

      include 'SIZE'
      include 'TSTEP'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real identity(ls,ls),eigv(ls,ls),w(ls,ls)
      real vv(ls,ls)

      if (nio.eq.0) write (6,*) 'inside genevec'

      call rzero(identity,ls*ls)

      do j=1,ls
         identity(j,j) = 1
      enddo

      call copy(vv,uu,ls*ls)

      call generalev(vv,identity,eval,ls,w)
      call copy(eigv,vv,ls*ls)

      if (nio.eq.0) write (6,*)'number of modes:',nb

      do l = 1,nb
         call copy(evec(1,l),eigv(1,ls-l+1),ls) ! reverse order of eigv
      enddo

      do i=1,ns
         if (nio.eq.0) write (6,'(i5,1p1e16.6,3x,a)')
     $      i,eval(ns-i+1),'eval'
      enddo

      if (nio.eq.0) write (6,*) 'exiting genevec'

      return
      end
c-----------------------------------------------------------------------
      subroutine scale_bases

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scruz/ h1(lt),h2(lt)

      nio=-1
      if (ifl2) then
         do i=1,nb
            p=wl2prod(ub(1,i),vb(1,i),wb(1,i),ub(1,i),vb(1,i),wb(1,i))
            s=1./sqrt(p)
            call opcmult(ub(1,i),vb(1,i),wb(1,i),s)
         enddo
      else
         n=lx1*ly1*lz1*nelv
         call rone(h1,n)
         call rzero(h2,n)
         do i=1,nb
            p=h10prod(ub(1,i),vb(1,i),wb(1,i),
     $                ub(1,i),vb(1,i),wb(1,i),h1,h2)
            s=1./sqrt(p)
            call opcmult(ub(1,i),vb(1,i),wb(1,i),s)
         enddo
      endif
      nio=nid

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_hyperpar

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ep
      real uw(lt),vw(lt),ww(lt),h1(lt),h2(lt)
      real tmp1(nb),tmp2(nb),delta(nb)
      real work(ls,nb)

      ! eps is the free parameter
      ! 1e-2 is used in the paper
      ep = 1e-2

      n  = lx1*ly1*lz1*nelt

      call rone (h1,n)
      call rzero(h2,n)

      do j=1,nb                    ! compute hyper-parameter
         call axhelm(uw,ub(1,j),h1,h2,1,1)
         call axhelm(vw,vb(1,j),h1,h2,1,1)
         if (ldim.eq.3) call axhelm(ww,wb(1,j),h1,h2,1,1)
         do i=1,ls
            work(i,j) = glsc2(us(1,i),uw,n)+glsc2(vs(1,i),vw,n)
            if (ldim.eq.3) work(i,j) = work(i,j)+glsc2(ws(1,i),ww,n)
         enddo
         tmp1(j) = vlmin(work(:,j),ls)
         tmp2(j) = vlmax(work(:,j),ls)
         delta(j) = tmp2(j)-tmp1(j)
         sample_min(j) = tmp1(j) - ep * delta(j)
         sample_max(j) = tmp2(j) + ep * delta(j)
         write(6,*) j,sample_min(j),sample_max(j)
      enddo

      ! compute distance between sample_max and sample_min
      call sub3(sam_dis,sample_max,sample_min,nb)
      if (nid.eq.0) then
         do i=1,nb
            write(6,*)i,sam_dis(i)
         enddo
      endif

      if (nid.eq.0) then
         open (unit=51,file='sample_min')
         do i=1,nb
            write (51,*) sample_min(i)
         enddo
         close (unit=51)

         open (unit=52,file='sample_max')
         do i=1,nb
            write (52,*) sample_max(i)
         enddo
         close (unit=52)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine opbinv1_nom(out1,out2,out3,inp1,inp2,inp3,SCALE)
C--------------------------------------------------------------------
C
C     Compute OUT = (B)-1 * INP   (explicit)
C
C--------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'

      real out1  (1)
      real out2  (1)
      real out3  (1)
      real inp1  (1)
      real inp2  (1)
      real inp3  (1)

      include 'OPCTR'

c     call opmask  (inp1,inp2,inp3)
      call opdssum (inp1,inp2,inp3)

      ntot=lx1*ly1*lz1*nelv

      if (if3d) then
         do 100 i=1,ntot
            tmp    =binvm1(i,1,1,1)*scale
            out1(i)=inp1(i)*tmp
            out2(i)=inp2(i)*tmp
            out3(i)=inp3(i)*tmp
  100    continue
      else
         do 200 i=1,ntot
            tmp    =binvm1(i,1,1,1)*scale
            out1(i)=inp1(i)*tmp
            out2(i)=inp2(i)*tmp
  200    continue
      endif

      return
      end
c-----------------------------------------------------------------------
