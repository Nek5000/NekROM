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
            call axhelm(au,ub(1,i),h1,h2,1,1)
            call axhelm(av,vb(1,i),h1,h2,1,1)
            if (ldim.eq.3) call axhelm(aw,wb(1,i),h1,h2,1,1)
            call opbinv1_nom(au,av,aw,au,av,aw,1.)
            call outpost(au,av,aw,pr,t,'aaa')
            call comp_pdrag(fd1(1,i),ub(1,i),vb(1,i),wb(1,i))
            call comp_pdrag(fd3(1,i),au,av,aw)
         enddo
      endif

      do i=0,nb
         write (6,*) 'fd3',fd3(1,i)
      enddo

      if (nio.eq.0) write (6,*) 'exiting setbases'

      return
      end
c-----------------------------------------------------------------------
      subroutine proj2bases(coef,u1,u2,u3)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real coef(0:nb),u1(lt),u2(lt),u3(lt)

      if (ifl2) then
         call wl2proj(coef,u1,u2,u3)
      else
         call h10proj(coef,u1,u2,u3)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine h10proj(coef,u1,u2,u3)

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real u1(lt),u2(lt),u3(lt)

      common /scrk3/ t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)
      common /scrk4/ h1(lt),h2(lt)

      real coef(0:nb)

      if (nio.eq.0) write (6,*) 'inside h10proj'

      n=lx1*ly1*lz1*nelt

      call rone(h1,n)
      call rzero(h2,n)
      call opsub3(t1,t2,t3,u1,u2,u3,ub,vb,wb)

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
      subroutine wl2proj(coef,u1,u2,u3)

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real u1(lt),u2(lt),u3(lt)

      common /scrk3/ t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)
      common /scrk4/ h1(lt),h2(lt),bwm1(lt)

      real coef(0:nb)

      if (nio.eq.0) write (6,*) 'inside wl2proj'

      n=lx1*ly1*lz1*nelt

      call col3(bwm1,bm1,wm1,n)
      call opsub3(t1,t2,t3,u1,u2,u3,ub,vb,wb)

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
      subroutine gengram

      include 'SIZE'
      include 'MOR'

      if (ifl2) then
         call gengraml2
      else
         call gengramh10
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine gengramh10

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real uw(lt),vw(lt),ww(lt),h1(lt),h2(lt)
      real u0(lt,3)

      if (nio.eq.0) write (6,*) 'inside gengramh10'

      n  = lx1*ly1*lz1*nelt
      ns = ls

      ! copy zero mode to u0(1,1:3)
      call opcopy(u0(1,1),u0(1,2),u0(1,3),ub(1,0),vb(1,0),wb(1,0))

      call rone (h1,n)
      call rzero(h2,n)

      do j=1,ns ! Form the Gramian, U=U_K^T A U_K using H^1_0 Norm
         call axhelm(uw,ust(1,j),h1,h2,1,1)
         call axhelm(vw,vst(1,j),h1,h2,1,1)
         if (ldim.eq.3) call axhelm(ww,wst(1,j),h1,h2,1,1)
         do i=1,ns
            uu(i,j) = glsc2(ust(1,i),uw,n)+glsc2(vst(1,i),vw,n)
            if (ldim.eq.3) uu(i,j) = uu(i,j)+glsc2(wst(1,i),ww,n)
            if (nio.eq.0) write (99,*) uu(i,j)
         enddo
         if (nio.eq.0) write(6,1) j,uu(1,j)
      enddo

      if (nio.eq.0) write (6,*) 'exiting gengramh10'
    1 format (' uu',i5,1p1e16.6)

      return
      end
c-----------------------------------------------------------------------
      subroutine gengraml2

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real uw(lt),vw(lt),ww(lt),h1(lt),h2(lt)
      real u0(lt,3)
      real bwm1(lt)

      if (nio.eq.0) write (6,*) 'inside gengraml2'

      n=lx1*ly1*lz1*nelv

      call col3(bwm1,bm1,wm1,n)

      if (nio.eq.0) write (6,*) 'ns',ns

      do j=1,ns ! Form the Gramian, U=U_K^T A U_K using L2 Norm
      do i=1,ns
         if (ifvort) then
            uu(i,j)=glsc3(us(1,i),us(1,j),bwm1,n)
         else
            uu(i,j) = op_glsc2_wt(ust(1,i),vst(1,i),wst(1,i),
     $                            ust(1,j),vst(1,j),wst(1,j),bwm1)
         endif
      enddo
         if (nio.eq.0) write (6,1) j,uu(1,j)
      enddo

      if (nio.eq.0) write (6,*) 'exiting gengraml2'

    1 format (' uu',i5,' ',1p1e16.6)

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
C
      REAL OUT1  (1)
      REAL OUT2  (1)
      REAL OUT3  (1)
      REAL INP1  (1)
      REAL INP2  (1)
      REAL INP3  (1)
C

      include 'OPCTR'
C
#ifdef TIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'opbnv1'
      endif
#endif
C
c     CALL OPMASK  (INP1,INP2,INP3)
      CALL OPDSSUM (INP1,INP2,INP3)
C
      NTOT=lx1*ly1*lz1*NELV
C
#ifdef TIMER
      isbcnt = ntot*(1+ldim)
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
#endif
C
      IF (IF3D) THEN
         DO 100 I=1,NTOT
            TMP    =BINVM1(I,1,1,1)*scale
            OUT1(I)=INP1(I)*TMP
            OUT2(I)=INP2(I)*TMP
            OUT3(I)=INP3(I)*TMP
  100    CONTINUE
      ELSE
         DO 200 I=1,NTOT
            TMP    =BINVM1(I,1,1,1)*scale
            OUT1(I)=INP1(I)*TMP
            OUT2(I)=INP2(I)*TMP
  200    CONTINUE
      ENDIF
C
      return
      END
c-----------------------------------------------------------------------
