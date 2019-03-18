c-----------------------------------------------------------------------
      subroutine setbases

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real u0(lt,3)
      common /scrk3/ t4(lt),t5(lt),t6(lt)
      common /scrk4/ bwm1(lt)
      common /scrk5/ a1(lt),a2(lt),a3(lt)

      if (nio.eq.0) write (6,*) 'inside setbases'

      n=lx1*ly1*lz1*nelt

      if (ifread) then
         call loadbases(ub,vb,wb,nb)
      else
         n=lx1*ly1*lz1*nelt

         do i=1,nb
            call rzero(ub(1,i),n)
            call rzero(vb(1,i),n)
            if (ldim.eq.3) call rzero(wb(1,i),n)
         enddo

         ! ub, vb, wb, are the modes
         do j=1,ns
         do i=1,nb
            call opadds(ub(1,i),vb(1,i),wb(1,i),
     $         us0(1,1,j),us0(1,2,j),us0(1,ldim,j),evec(j,i,1),n,2)
         enddo
         enddo

         call scale_vbases

         if (ifpod(2)) then
            do i=1,nb
               call rzero(tb(1,1,i),n)
               do j=1,ns
                  call add2s2(tb(1,1,i),ts0(1,j,1),evec(j,i,2),n)
               enddo
            enddo
            call scale_tbases
         endif
      endif

      if (ifdrago) then
         do i=0,nb
            call lap2d(a1,ub(1,i))
            call lap2d(a2,vb(1,i))
            if (ldim.eq.3) call lap2d(a3,wb(1,i))
            call opcmult(a1,a2,a3,param(2))
            call comp_pdrag(fd1(1,i),ub(1,i),vb(1,i),wb(1,i))
            call comp_pdrag(fd3(1,i),a1,a2,a3)
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

      real coef(0:nb)

      if (nio.eq.0) write (6,*) 'inside h10proj'

      n=lx1*ly1*lz1*nelt

      call opsub3(t1,t2,t3,ux,uy,uz,ub,vb,wb)

      coef(0) = 1.
      if (nio.eq.0) write (6,1) coef(0),coef(0),1.

      do i=1,nb
         call axhelm(t4,ub(1,i),ones,zeros,1,1)
         call axhelm(t5,vb(1,i),ones,zeros,1,1)

         ww = glsc2(t4,ub(1,i),n)+glsc2(t5,vb(1,i),n)
         vv = glsc2(t4,t1,n)+glsc2(t5,t2,n)

         if (ldim.eq.3) then
            call axhelm(t6,wb(1,i),ones,zeros,1,1)
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
      common /scrk4/ bwm1(lt)

      real coef(0:nb)

      if (nio.eq.0) write (6,*) 'inside wl2proj'

      n=lx1*ly1*lz1*nelt

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
     $         ub(1,i),vb(1,i),wb(1,i),ub(1,i),vb(1,i),wb(1,i),bm1)
            vv = op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),t1,t2,t3,bm1)

            coef(i) = vv/ww
            if (nio.eq.0) write (6,1) coef(i),vv,ww
         enddo
      endif

      if (nio.eq.0) write (6,*) 'exiting wl2proj'

    1 format(' wl2coef',1p3e16.8)

      return
      end
c-----------------------------------------------------------------------
      function scaprod(t1,t2)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt)

      if (ips.eq.'L2 ') then
         vecprod=wl2sprod(t1,t2)
      else if (ips.eq.'H10') then
         vecprod=h10sprod(t1,t2)
      else
         call exitti('did not provide supported inner product space$')
      endif

      return
      end
c-----------------------------------------------------------------------
      function vecprod(t1,t2,t3,t4,t5,t6)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)

      if (ips.eq.'L2 ') then
         vecprod=wl2vprod(t1,t2,t3,t4,t5,t6)
      else if (ips.eq.'H10') then
         vecprod=h10vprod(t1,t2,t3,t4,t5,t6)
      else
         call exitti('did not provide supported inner product space$')
      endif

      return
      end
c-----------------------------------------------------------------------
      function h10sprod(t1,t2)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt)

      common /scrk3/ t3(lt)

      if (nio.eq.0) write (6,*) 'inside h10sprod'

      call axhelm(t3,t1,ones,zeros,1,1)
      h10sprod=glsc2(t3,t2,lx1*ly1*lz1*nelt)

      if (nio.eq.0) write (6,*) 'exiting h10sprod'

      return
      end
c-----------------------------------------------------------------------
      function h10vprod(t1,t2,t3,t4,t5,t6)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)

      common /scrk3/ t7(lt),t8(lt),t9(lt)

      if (nio.eq.0) write (6,*) 'inside h10vprod'

      n=lx1*ly1*lz1*nelt

      call axhelm(t7,t1,ones,zeros,1,1)
      h10vprod=glsc2(t7,t4,n)

      call axhelm(t8,t2,ones,zeros,1,1)
      h10vprod=h10vprod+glsc2(t8,t5,n)

      if (ldim.eq.3) then
         call axhelm(t9,t3,ones,zeros,1,1)
         h10vprod = h10vprod + glsc2(t9,t6,n)
      endif

      if (nio.eq.0) write (6,*) 'exiting h10vprod'

      return
      end
c-----------------------------------------------------------------------
      function wl2sprod(t1,t2)

      include 'SIZE'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt)

      if (nio.eq.0) write (6,*) 'inside wl2sprod'

      n=lx1*ly1*lz1*nelt

      wl2sprod = glsc3(t1,t2,bm1,n)

      if (nio.eq.0) write (6,*) 'exiting wl2sprod'

      return
      end
c-----------------------------------------------------------------------
      function wl2vprod(t1,t2,t3,t4,t5,t6)

      include 'SIZE'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)

      if (nio.eq.0) write (6,*) 'inside wl2vprod'

      wl2vprod=op_glsc2_wt(t1,t2,t3,t4,t5,t6,bm1)

      if (nio.eq.0) write (6,*) 'exiting wl2vprod'

      return
      end
c-----------------------------------------------------------------------
      subroutine setgram

      include 'SIZE'
      include 'MOR'

      if (.not.ifread) then
         jfield=ifield
         ifield=1
         if (ifpod(0)) call gengram(ug(1,1,0),dps,ns,ldim)
         if (ifpod(1)) call gengram(ug(1,1,1),us0,ns,ldim)
         do i=2,ldimt1
            ifield=i
            if (ifpod(i)) call gengram(ug(1,1,i),ts0(1,1,i-1),ns,1)
         enddo
         ifield=jfield
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine setevec

      include 'SIZE'
      include 'MOR'

      do i=0,ldimt1
         if (ifpod(i)) call genevec(evec(1,1,i),eval(1,i),ug(1,1,i))
      enddo
c     call genevec(evec(1,1,1),eval(1,1),ug(1,1,1))

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

      common /scrgram/ uw(lt),vw(lt),ww(lt)

      if (nio.eq.0) write (6,*) 'inside gengramh10'

      n  = lx1*ly1*lz1*nelt

      do j=1,ms ! Form the Gramian, U=U_K^T A U_K using H^1_0 Norm
         call axhelm(uw,s(1,1,j),ones,zeros,1,1)
         if (mdim.ge.2) call axhelm(vw,s(1,2,j),ones,zeros,1,1)
         if (mdim.ge.3) call axhelm(ww,s(1,3,j),ones,zeros,1,1)
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

      common /scrgram/ uw(lt),vw(lt),ww(lt)
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
      subroutine genevec(vec,val,gram)

      !!! does not work if ns.lt.ls !!!

      include 'SIZE'
      include 'TSTEP'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      integer icalld
      save    icalld
      data    icalld /0/
      
      common /scrgvec/ eye(ls,ls),gc(ls,ls),wk(ls,ls),eigv(ls,ls)

      real gram(ls,ls),vec(ls,nb),val(ls)

      if (nio.eq.0) write (6,*) 'inside genevec'

      if (icalld.eq.0) then
         call rzero(eye,ls*ls)
         do j=1,ls
            eye(j,j) = 1.
         enddo
         icalld=1
      endif

      call copy(gc,gram,ls*ls)

      call generalev(gc,eye,val,ls,wk)
      call copy(eigv,gc,ls*ls)

      do l = 1,nb
         call copy(vec(1,l),eigv(1,ls-l+1),ls) ! reverse order of eigv
      enddo

      do i=1,ns
         if (nio.eq.0) write (6,'(i5,1p1e16.6,3x,a)')
     $      i,val(ns-i+1),'eval'
      enddo

      if (nio.eq.0) write (6,*) 'exiting genevec'

      return
      end
c-----------------------------------------------------------------------
      subroutine scale_vbases

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      jfield=ifield
      ifield=1
      nio=-1
      do i=1,nb
         p=vecprod(ub(1,i),vb(1,i),wb(1,i),ub(1,i),vb(1,i),wb(1,i))
         s=1./sqrt(p)
         call opcmult(ub(1,i),vb(1,i),wb(1,i),s)
      enddo
      nio=nid
      ifield=jfield

      return
      end
c-----------------------------------------------------------------------
      subroutine scale_tbases

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      nio=-1
      do i=1,nb
         p=scaprod(tb(1,1,i),tb(1,1,i))
         s=1./sqrt(p)
         call cmult(tb(1,1,i),s,lx1*ly1*lz1*nelt)
      enddo
      nio=nid

      return
      end
c-----------------------------------------------------------------------
