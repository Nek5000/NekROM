c-----------------------------------------------------------------------
      subroutine genbases

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real u0(lt,3)
      common /scrk3/ t4(lt),t5(lt),t6(lt)
      common /scrk4/ h1(lt),h2(lt),bwm1(lt)

      if (nio.eq.0) write (6,*) 'inside genbases'

      n  = lx1*ly1*lz1*nelt

      call rone(h1,n)
      call rzero(h2,n)
      call col3(bwm1,bm1,wm1,n)

      ONE = 1.
      ZERO= 0.

      ns = ls ! REQUIRED: get_saved_fields overwrites ns argument
      call opcopy(u0(1,1),u0(1,2),u0(1,3),ub(1,0),vb(1,0),wb(1,0))

      ! ub, vb, wb, are the modes
      call dgemm( 'N','N',n,nb,ls,ONE,us,lt,evec,ls,ZERO,ub(1,1),lt)
      call dgemm( 'N','N',n,nb,ls,ONE,vs,lt,evec,ls,ZERO,vb(1,1),lt)
      if (ldim.eq.3)
     $call dgemm( 'N','N',n,nb,ls,ONE,ws,lt,evec,ls,ZERO,wb(1,1),lt)

      call scale_bases

      itmp = istep
      ttmp = time
      do i=0,nb ! dump the generated modes
         istep = i
         time = real(istep)
         call outpost(ub(1,i),vb(1,i),wb(1,i),pr,t,'bas')
      enddo
      istep = itmp
      time = ttmp

      if (nio.eq.0) write (6,*) 'exiting genbases'

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
      if (nio.eq.0) write (6,*) 'h10coef', coef(0),coef(0),1

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
         if (nio.eq.0) write (6,*) 'h10coef', coef(i),vv,ww
      enddo

      if (nio.eq.0) write (6,*) 'exiting h10proj'

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

      if (nio.eq.0) write (6,*) 'wl2coef', coef(0),coef(0),1

      do i=1,nb
         ww = op_glsc2_wt(
     $      ub(1,i),vb(1,i),wb(1,i),ub(1,i),vb(1,i),wb(1,i),bwm1)
         vv = op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),t1,t2,t3,bwm1)

         coef(i) = vv/ww
      enddo
            if (nio.eq.0) write (6,*) 'wl2coef', coef(i),vv,ww

      if (nio.eq.0) write (6,*) 'exiting wl2proj'

      return
      end
c-----------------------------------------------------------------------
      function h10prod(t1,t2,t3,t4,t5,t6,h1,h2)

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt)
      real h1(lt),h2(lt)

      common /scrk3/ t7(lt),t8(lt),t9(lt)

      if (nio.eq.0) write (6,*) 'inside h10prod'

      n=lx1*ly1*lz1*nelt

      call axhelm(t7,t1,h1,h2,1,1)
      call axhelm(t8,t2,h1,h2,1,1)

      h10prod = glsc2(t7,t4,n)+glsc2(t8,t5,n)

      if (ldim.eq.3) then
         call axhelm(t9,t3,h1,h2,1,1)
         h10prod = h10prod + glsc2(t9,t6,n)
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

      wl2prod = op_glsc2_wt(t1,t2,t3,t4,t5,t6,bwm1)

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

      call opcopy(u0(1,1),u0(1,2),u0(1,3),ub(1,0),vb(1,0),wb(1,0))

      call rone (h1,n)
      call rzero(h2,n)

      do j=1,ns ! Form the Gramian, U=U_K^T A U_K using H^1_0 Norm
         call axhelm(uw,us(1,j),h1,h2,1,1)
         call axhelm(vw,vs(1,j),h1,h2,1,1)
         if (ldim.eq.3) call axhelm(ww,ws(1,j),h1,h2,1,1)
         do i=1,ns
            uu(i,j) = glsc2(us(1,i),uw,n)+glsc2(vs(1,i),vw,n)
            if (ldim.eq.3) uu(i,j) = uu(i,j)+glsc2(ws(1,i),ww,n)
         enddo
         if (nio.eq.0) write(6,*) j,uu(1,j),' uu'
      enddo

      if (nio.eq.0) write (6,*) 'exiting gengramh10'

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

      call outpost(bwm1,wm1,vz,pr,t,'bw1')

      do j=1,ns ! Form the Gramian, U=U_K^T A U_K using L2 Norm
      do i=1,ns
         uu(i,j) = op_glsc2_wt(us(1,i),vs(1,i),ws(1,i),
     $                         us(1,j),vs(1,j),ws(1,j),bwm1)
         write (88,*) uu(i,j)
      enddo
         if (nio.eq.0) write (6,*) 'uu',uu(1,j)
      enddo

      if (nio.eq.0) write (6,*) 'exiting gengraml2'

      return
      end
c-----------------------------------------------------------------------
      subroutine genevec

      !!! does not work if ns.lt.ls !!!

      include 'SIZE'
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

      call generalev(vv,identity,gram_eig,ls,w)
      call copy(eigv,vv,ls*ls)

c     eig = eig(ls:1:-1)

      if (nio.eq.0) write(6,*)'number of mode:',nb

      do l = 1,nb
         call copy(evec(1,l),eigv(1,ls-l+1),ls) ! reverse order of eigv
         if (nio.eq.0) write (6,*) 'eigenvalue',l,gram_eig(l)
      enddo


      if (nio.eq.0) write (6,*) 'exiting genevec'

      return
      end
c-----------------------------------------------------------------------
      subroutine scale_bases

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scruz/ h1(lt),h2(lt)

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
