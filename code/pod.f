c-----------------------------------------------------------------------
      subroutine get_saved_fields(usave,vsave,wsave,nsave,u0)

c     This routine reads files specificed in file.list

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      parameter (lt=lx1*ly1*lz1*lelt)
      real usave(lt,nsave),vsave(lt,nsave),wsave(lt,nsave)
      real u0(lt,3) ! Initial condtion


      ierr = 0
      if (nid.eq.0) open(77,file='file.list',status='old',err=199)
      ierr = iglmax(ierr,1)
      if (ierr.gt.0) goto 199
      n = lx1*ly1*lz1*nelt
      n2= lx2*ly2*lz2*nelt

      icount = 0
      do ipass=1,nsave

         call blank(initc,127)
         initc(1) = 'done '
         if (nid.eq.0) read(77,127,end=998) initc(1)
  998    call bcast(initc,127)
  127    format(a127)

         if (indx1(initc,'done ',5).eq.0) then ! We're not done
            nfiles = 1
            call restart(nfiles)  ! Note -- time is reset.

!           Usave = U_snapshot - U_stokes:

            call opsub3 (usave(1,ipass),vsave(1,ipass),wsave(1,ipass)
     $                  ,vx,vy,vz,u0(1,1),u0(1,2),u0(1,3))

            icount = icount+1
         else
            goto 999
         endif

      enddo

  999 continue  ! clean up averages
      if (nid.eq.0) close(77)

      nsave = icount ! Actual number of files read

      return

  199 continue ! exception handle for file not found
      ierr = 1
      if (nid.eq.0) ierr = iglmax(ierr,1)
      call exitti('Auto averager did not find list file.$',ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine genmodes

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real usave(lt,ls),vsave(lt,ls),wsave(lt,ls)
      real uu(ls,ls),u0(lt,3),evec(ls,nb)

      if (nio.eq.0) write (6,*) 'inside genmodes'

      n  = lx1*ly1*lz1*nelt

      call gengram(uu)
      call genevec(evec,uu)

      ONE = 1.
      ZERO= 0.

      ns = ls ! REQUIRED: get_saved_fields overwrites ns argument
      call opcopy(u0(1,1),u0(1,2),u0(1,3),ub(1,0),vb(1,0),wb(1,0))
      call get_saved_fields(usave,vsave,wsave,ns,u0)

      ! ub, vb, wb, are the modes
      call dgemm( 'N','N',n,nb,ls,ONE,usave,lt,evec,ls,ZERO,ub(1,1),lt)
      call dgemm( 'N','N',n,nb,ls,ONE,vsave,lt,evec,ls,ZERO,vb(1,1),lt)
      if (ldim.eq.3)
     $call dgemm( 'N','N',n,nb,ls,ONE,wsave,lt,evec,ls,ZERO,wb(1,1),lt)

      do i=0,nb ! dump the generated modes
         call outpost(ub(1,i),vb(1,i),wb(1,i),pr,t,'bas')
      enddo

      if (nio.eq.0) write (6,*) 'exiting genmodes'

      return
      end
c-----------------------------------------------------------------------
      subroutine h10proj(coef,t1,t2,t3)

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt)

      common /scrk3/ t4(lt),t5(lt),t6(lt)
      common /scrk4/ h1(lt),h2(lt)

      real coef(nb)

      if (nio.eq.0) write (6,*) 'inside h10proj'

      n=lx1*ly1*lz1*nelt

      call rone(h1,n)
      call rzero(h2,n)

      do i=1,nb
         call axhelm(t4,ub(1,i),h1,h2,1,1)
         call axhelm(t5,vb(1,i),h1,h2,1,1)

         uu = glsc2(t4,ub(1,i),n)+glsc2(t5,vb(1,i),n)
         vv = glsc2(t4,t1,n)+glsc2(t5,t2,n)

         if (ldim.eq.3) then
            call axhelm(t6,wb(1,i),h1,h2,1,1)
            uu = uu + glsc2(t6,wb(1,i),n)
            vv = vv + glsc2(t6,t3,n)
         endif

         coef(i) = vv/uu
      enddo

      if (nio.eq.0) write (6,*) 'exiting h10proj'

      return
      end
c-----------------------------------------------------------------------
      subroutine l2proj(coef,t1,t2,t3)

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt)

      common /scrk3/ t4(lt),t5(lt),t6(lt)
      common /scrk4/ h1(lt),h2(lt)

      real coef(nb)

      if (nio.eq.0) write (6,*) 'inside l2proj'

      n=lx1*ly1*lz1*nelt

      call rone(h1,n)
      call rzero(h2,n)

      do i=1,nb
         uu = op_glsc2_wt(
     $      ub(1,i),vb(1,i),wb(1,i),ub(1,i),vb(1,i),wb(1,i),bm1)
         vv = op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),t1,t2,t3,bm1)

         coef(i) = vv/uu
      enddo

      if (nio.eq.0) write (6,*) 'exiting l2proj'

      return
      end
c-----------------------------------------------------------------------
      subroutine h10prod(prod,t1,t2,t3,t4,t5,t6,h1,h2)

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt)
      real h1(lt),h2(lt)

      common /scrk3/ t7(lt),t8(lt),t9(lt)

      real coef(nb)

      if (nio.eq.0) write (6,*) 'inside h10prod'

      n=lx1*ly1*lz1*nelt

      call axhelm(t7,t1,h1,h2,1,1)
      call axhelm(t8,t2,h1,h2,1,1)
      if (ldim.eq.3) call axhelm(t9,t3,h1,h2,1,1)

      prod = glsc2(t7,t4,n)+glsc2(t8,t5,n)

      if (nio.eq.0) write (6,*) 'exiting h10prod'

      return
      end
c-----------------------------------------------------------------------
      subroutine gengram(uu)

      include 'SIZE'
      include 'MOR'

      if (ifl2) then
         call gengraml2(uu)
      else
         call gengramh10(uu)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine gengramh10(uu)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real usave(lt,ls),vsave(lt,ls),wsave(lt,ls)
      real uw(lt),vw(lt),ww(lt),h1(lt),h2(lt)
      real u0(lt,3)
      real uu(ls,ls)

      if (nio.eq.0) write (6,*) 'inside gengramh10'

      n  = lx1*ly1*lz1*nelt
      ns = ls

      call opcopy(u0(1,1),u0(1,2),u0(1,3),ub(1,0),vb(1,0),wb(1,0))
      call get_saved_fields(usave,vsave,wsave,ns,u0)

      call rone (h1,n)
      call rzero(h2,n)

      do j=1,ns ! Form the Gramian, U=U_K^T A U_K using H^1_0 Norm
         call axhelm(uw,usave(1,j),h1,h2,1,1)
         call axhelm(vw,vsave(1,j),h1,h2,1,1)
         if (ldim.eq.3) call axhelm(ww,wsave(1,j),h1,h2,1,1)
         do i=1,ns
            uu(i,j) = glsc2(usave(1,i),uw,n)+glsc2(vsave(1,i),vw,n)
            if (ldim.eq.3) uu(i,j) = uu(i,j)+glsc2(wsave(1,i),ww,n)
         enddo
         if (nio.eq.0) write(6,*) j,uu(1,j),' uu'
      enddo

      if (nio.eq.0) write (6,*) 'exiting gengramh10'

      return
      end
c-----------------------------------------------------------------------
      subroutine gengraml2(uu)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real usave(lt,ls),vsave(lt,ls),wsave(lt,ls)
      real uw(lt),vw(lt),ww(lt),h1(lt),h2(lt)
      real u0(lt,3)

      real uu(ls,ls),Identity(ls,ls),eig(ls),eigv(ls,ls),w(ls,ls)
      real u0r(ls)

      if (nio.eq.0) write (6,*) 'inside gengraml2'

      ns = ls ! REQUIRED: get_saved_fields overwrites ns argument

      call opcopy(u0(1,1),u0(1,2),u0(1,3),ub(1,0),vb(1,0),wb(1,0))
      call get_saved_fields(usave,vsave,wsave,ns,u0)

      do j=1,ns ! Form the Gramian, U=U_K^T A U_K using L2 Norm
      do i=1,ns
         uu(i,j) = op_glsc2_wt(usave(1,i),vsave(1,i),wsave(1,i),
     $                         usave(1,j),vsave(1,j),wsave(1,j),bm1)
      enddo
      enddo

      if (nio.eq.0) write (6,*) 'exiting gengraml2'

      return
      end
c-----------------------------------------------------------------------
      subroutine genevec(evec,uu)

      !!! does not work if ns.lt.ls !!!

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real usave(lt,ls),vsave(lt,ls),wsave(lt,ls)
      real uu(ls,ls),identity(ls,ls),eig(ls),eigv(ls,ls),w(ls,ls)

      real evec(ls,nb)

      if (nio.eq.0) write (6,*) 'inside genevec'

      call rzero(identity,ls*ls)

      do j=1,ls
         identity(j,j) = 1
      enddo

      call generalev(uu,identity,eig,ls,w)
      call copy(eigv,uu,ls*ls)

c     eig = eig(ls:1:-1)

      if (nio.eq.0) write(6,*)'number of mode:',nb

      do l = 1,nb
         call copy(evec(1,l),eigv(1,ls-l+1),ls) ! reverse order of eigv
      enddo

      if (nio.eq.0) write (6,*) 'exiting genevec'

      return
      end
c-----------------------------------------------------------------------
