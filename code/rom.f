c-----------------------------------------------------------------------
      include 'pod.f'
      include 'read.f'
      include 'aux.f'
      include 'dump.f'
c-----------------------------------------------------------------------
      subroutine readops

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      if (nio.eq.0) write (6,*) 'inside readops'

      call readic(u,nb+1)

      call reada0(a0,(nb+1)**2)
      call readb0(b0,(nb+1)**2)

      if (nid.eq.(np-1)) call readc0(c0,(nb+1)**3)

      if (nio.eq.0) write (6,*) 'exiting readops'

      return
      end
c-----------------------------------------------------------------------
      subroutine genops

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      if (nio.eq.0) write (6,*) 'inside genops'

      call makea0
      call makeb0
      call makec0
      
      call makeic

      if (nio.eq.0) write (6,*) 'exiting genops'

      return
      end
c-----------------------------------------------------------------------
      subroutine setops

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      if (nio.eq.0) write (6,*) 'inside setops'

      do i=1,nb
         call copy(a(1,i),a0(1,i),nb)
         call copy(b(1,i),b0(1,i),nb)
      enddo

      if (nid.eq.(np-1)) then
      do j=1,nb
      do i=1,nb
         call copy(c(1,i,j),c0(1,i,j),nb)
      enddo
      enddo
      endif

      if (np.gt.1) call setcloc

      if (nio.eq.0) write (6,*) 'exiting setops'

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_init_params

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      if (nid.eq.0) write (6,*) 'inside rom_init_params'

      time = 0.

      ad_nsteps=nsteps
      ad_iostep=iostep

      ad_dt = dt
      ad_re = 1/param(2)

      ifl2=.false.
      if (param(50).eq.0) ifl2=.true.

      ifvort=.false. ! default to false for now
      ifdump=.false.

      call compute_BDF_coef(ad_alpha,ad_beta)

      if (nid.eq.0) write (6,*) 'exiting rom_init_params'

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_init_fields

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk0/ t1(lt),t2(lt),t3(lt),u0(lx1*ly1*lz1*lelt,3)

      if (nid.eq.0) write (6,*) 'inside rom_init_fields'

      n=lx1*ly1*lz1*nelv
      call rone(wm1,n)

      if (ifvort) then
         call comp_vort3(u0,t1,t2,ub,vb,wb)
         call comp_vort3(t1,t2,t3,vx,vy,vz)

         call copy(ub,u0,n)
         call rzero(vb,n)
         call rzero(wb,n)

         call copy(vx,t1,n)
         call rzero(vy,n)
         call rzero(vz,n)
      else
         call opcopy(u0,u0(1,2),u0(1,3),ub,vb,wb)
      endif

      ns = ls
      call get_saved_fields(us,vs,ws,ns,u0,ifvort)

      if (nid.eq.0) write (6,*) 'exiting rom_init_fields'

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_setup

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      if (nid.eq.0) write (6,*) 'inside rom_setup'

      ! BDFk/EXTk coefficients ( will change to BD inside Nek)
      call compute_BDF_coef(ad_alpha,ad_beta)

      if (nid.eq.0) write (6,*) 'exiting rom_setup'

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_step

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

c     Matrices and vectors for advance
      real tmp(0:nb),tmat(nb,nb+1)
      real coef(1:nb), e0(0:nb)

      common /scrk3/ work(lt)
      common /scrk1/ t1(lt),t2(lt),t3(lt)

c     Variable for vorticity
      real vort(lt,3)

c     Working arrays for LU

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

c     if (nio.eq.0) write (6,*) 'entering rom_step'

      n  = lx1*ly1*lz1*nelt

      time=time+ad_dt

      count = min0(ad_step,3)

      call rzero(e0,nb+1)
      e0(0) = 1

      if (ad_step.le.3) then
         call cmult2(flu,b,ad_beta(1,count)/ad_dt,nb*nb)
         call add2s2(flu,a,1/ad_re,nb*nb)
      endif

      ONE = 1.
      ZERO= 0.

      call mxm(u,nb+1,ad_beta(2,count),3,tmp,1)
c     call mxm(b0,nb+1,tmp,nb+1,rhs,1)
      call mxm(b,nb,tmp(1),nb,rhs,1)

      call cmult(rhs,-1/ad_dt,nb+1)

      s=-1.0/ad_re

c     call add2s2(rhs,a0,s,nb+1) ! not working...
      do i=0,nb
         rhs(i)=rhs(i)+s*a0(i,0)
      enddo

      call copy(conv(1,3),conv(1,2),nb)
      call copy(conv(1,2),conv(1,1),nb)

      if (np.eq.1) then
         call mxm(c,nb*(nb+1),u,nb+1,tmat,1)
         call mxm(tmat,nb,u,nb+1,conv,1)
      else
         call evalc(conv)
      endif

c     if (nio.eq.0) write (6,*) ''

c     do i=1,nb
c        if (nio.eq.0) write (6,*) i,conv(i,1),'dbg conv'
c     enddo

      call mxm(conv,nb,ad_alpha(1,count),3,tmp,1)

      call sub2(rhs,tmp,nb)

      if (ad_step.le.3) call lu(flu,nb,nb,ir,ic)

      call solve(rhs(1),flu,1,nb,nb,ir,ic)

      call copy(u(1,3),u(1,2),nb)
      call copy(u(1,2),u(1,1),nb)

      call copy(u(1,1),rhs,nb)

      call copy(coef,rhs,nb)

      if (mod(ad_step,ad_iostep).eq.0) then

!        This output is to make sure the ceof matches with matlab code

         if (nio.eq.0) then
            write (6,*)'ad_step:',ad_step,ad_iostep,npp,nid
            if (ad_step.eq.ad_nsteps) then
               do j=1,nb
                  write(6,*) j,u(j,1),'final'
               enddo
            else
               do j=1,nb
                  write(6,*) j,u(j,1)
               enddo
            endif
         endif

         if (ifdump) then
            call dumpcoef(u(:,1),nb,(ad_step/ad_iostep))

            call opzero(vx,vy,vz)
            do j=1,nb
               call opadds(vx,vy,vz,ub(1,j),vb(1,j),wb(1,j),coef(j),n,2)
            enddo
            call opadd2  (vx,vy,vz,ub,vb,wb)

            ! compute the vorticity of the ROM reconstructed field
            call opcopy(t1,t2,t3,vx,vy,vz)
            call comp_vort3(vort,work1,work2,t1,t2,t3)
            ifto = .true. ! turn on temp in fld file
            call copy(t,vort,n)

            call outpost (vx,vy,vz,pr,t,'rom')
         endif
      endif

c     if (nio.eq.0) write (6,*) 'exiting rom_step'

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_BDF_coef(ad_alpha,ad_beta)

      real ad_alpha(3,3), ad_beta(4,3)

      call rzero(ad_alpha,3*3)
      call rzero(ad_beta,3*4)

      ad_beta(1,1) = 1.
      ad_beta(2,1) = -1.

      ad_beta(1,2) = 1.5
      ad_beta(2,2) = -2
      ad_beta(3,2) = 0.5

      ad_beta(1,3) = 11./6
      ad_beta(2,3) = -3
      ad_beta(3,3) = 1.5
      ad_beta(4,3) = -1./3.

      ad_alpha(1,1)=1

      ad_alpha(1,2)=2
      ad_alpha(2,2)=-1

      ad_alpha(1,3)=3
      ad_alpha(2,3)=-3
      ad_alpha(3,3)=1

      return
      end
c-----------------------------------------------------------------------
      subroutine makec0

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real cux(lt), cuy(lt), cuz(lt)

      common /scrk1/ t1(lt), binv(lt)

      call invers2(binv,bm1,lx1*ly1*lz1*nelv)
      call rone(binv,lx1*ly1*lz1*nelv)

      if (nio.eq.0) write (6,*) 'inside makec'

      do k=0,nb
         call setcnv_c(ub(1,k),vb(1,k),wb(1,k))
         do j=0,nb
            call setcnv_u(ub(1,j),vb(1,j),wb(1,j))
            call ccu(cux,cuy,cuz)
            do i=0,nb
               c0(i,j,k) = op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),
     $                                cux,cuy,cuz,binv)
c              if (nid.eq.0) write (6,*) 'c0',i,j,k,c0(i,j,k)
            enddo
         enddo
      enddo

      if (nio.eq.0) write (6,*) 'exiting makec'

      return
      end
c-----------------------------------------------------------------------
      subroutine setcloc ! TODO Fix for 1-proc case

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (m=2)

      real cux(lt), cuy(lt), cuz(lt)
      integer mps(lt),mqs(lt),mrs(lt)

      real vr(lt), pts(lt), ur(lt)
      integer vl(lt),vi(m,lt),ui(m,lt)

      common /scrk1/ t1(lt), binv(lt)
      integer crh, parts(lt)
      common /romi/ crh
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      if (nio.eq.0) write (6,*) 'inside makecloc'

      ntot=nb*(nb+1)*(nb+1)

      nblock=nb*(nb+1)*(nb+1)/npp+1

      call fgslib_crystal_setup(crh,nekcomm,np)

      call izero(parts,lt)
      call setpart(parts,npp,ntot)

      call izero(mps,lt)
      call izero(mqs,lt)
      call izero(mrs,lt)

      call factor3(mp,mq,mr,npp)

      call nekgsync
c     write (6,*) 'factor3',nid,mp,mq,mr,npp

      call nekgsync
      call setpart3(mps,mqs,mrs,mp,mq,mr,nb)

      if (nio.eq.0) then
         write (6,*) ''
         do i=1,mp+1
            write (6,*) i,mps(i),'mps'
         enddo
         write (6,*) ''
         do i=1,mq+1
            write (6,*) i,mqs(i),'mqs'
         enddo
         write (6,*) ''
         do i=1,mr+1
            write (6,*) i,mrs(i),'mrs'
         enddo
      endif

      call nekgsync
      call setrange(mps,mqs,mrs,mp,mq,mr) ! set index range i0,i1,j0,j1,k0,k1

      call nekgsync
      write (6,1) nid,i0,i1,j0,j1,k0,k1

    1 format ('range',7(i4))

      kp=1
      nmax=lt
      ii=1
      ltrack=1

      call rzero(ctmp,nmax)

      if (nio.eq.0) write (6,*) 'first crystal_router'

      ! for testing purposes

      npr=np

      nmax=nb*(nb+1)**2/npr+1

      do ipr=0,npr-1
c        write (6,*) 'ipr',nid,ipr
         n=0
         if (nid.eq.npr-1) then
            call partialc(vr,n,nmax)
            do i=1,n
               vi(1,i)=ipr
               vi(2,i)=ii
               ii=ii+1
            enddo
         else
            n=0
         endif

         call sleep(nid)
c        write (6,*) nid,n,ipr,npr,'ninfo'
         call nekgsync
         call fgslib_crystal_tuple_transfer(
     $           crh,n,nmax,vi,m,vl,0,vr,1,kp)

         if (nid.eq.ipr) then
            call copy(ctmp,vr,n)
            call icopy(ui,vi,nmax*m)
            ncloc=n
c           write (6,*) 'inums',inums
         endif
      enddo

      call nekgsync
      call sleep(nid)

      do i=1,ncloc
         ii=mod(ui(2,i),nb)
         jj=mod((ui(2,i)-1)/nb,nb*(nb+1))
         kk=(ui(2,i)-1)/(nb*(nb+1))
         ui(1,i)=ijk2pid(ii,jj,kk,mps,mqs,mrs,mp,mq,mr)
c        write (6,*) nid,i,ctmp(i),ui(1,i),ui(2,i)
      enddo

      call nekgsync

      ! process vi such to find destination
      ! end process

      if (nio.eq.0) write (6,*) 'second crystal_router'

      call fgslib_crystal_tuple_transfer
     $   (crh,ncloc,lt,ui,m,vl,0,ctmp,1,kp)

      call nekgsync
      call sleep(nid)

c     do i=1,ncloc
c        write (6,*) nid,i,ctmp(i),ui(1,i),ui(2,i)
c     enddo

      call nekgsync

      call fgslib_crystal_free(crh)

      do i=1,ncloc
         clocal(i) = ctmp(i)
      enddo

      if (nio.eq.0) write (6,*) 'exiting makecloc'

      return
      end
c-----------------------------------------------------------------------
      subroutine setcnv_c(cx,cy,cz)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ltd=lxd*lyd*lzd*lelt)

      common /convect/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      real cx(lt), cy(lt), cz(lt)

      call set_convect_new(c1v,c2v,c3v,cx,cy,cz)

      return
      end
c-----------------------------------------------------------------------
      subroutine setcnv_u(ux,uy,uz)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ltd=lxd*lyd*lzd*lelt)

      common /convect/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      real ux(lt), uy(lt), uz(lt)

      call intp_rstd_all(u1v,ux,nelv)
      call intp_rstd_all(u2v,uy,nelv)
      if (ldim.eq.3) call intp_rstd_all(u3v,uz,nelv)

      return
      end
c-----------------------------------------------------------------------
      subroutine ccu(cu1,cu2,cu3) ! compute C(c) * u set by setcnv

      include 'SIZE'

      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(ltd=lxd*lyd*lzd*lelt)

      real cu1(lt), cu2(lt), cu3(lt)

      common /convect/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      call convect_new(cu1,u1v,.true.,c1v,c2v,c3v,.true.)
      call convect_new(cu2,u2v,.true.,c1v,c2v,c3v,.true.)
      if (ldim.eq.3)
     $ call convect_new(cu3,u3v,.true.,c1v,c2v,c3v,.true.)

      return
      end
c-----------------------------------------------------------------------
      subroutine intp_rstd_all(uf,u,nel)

      include 'SIZE'
      include 'INPUT'

      parameter (lxyz1=lx1*ly1*lz1)
      parameter (lxyzd=lxd*lyd*lzd)

      real uf(lxyzd,lelt), u(lxyz1,lelt)

      do i=1,nel
         call intp_rstd(uf(1,i),u(1,i),lx1,lxd,if3d,0) ! 0 --> forward
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine makea0

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrvh/ h1(lt),h2(lt)
      common /scrns/ usave(lt),vsave(lt),wsave(lt)

      if (nio.eq.0) write (6,*) 'inside makea'

      n=lx1*ly1*lz1*nelt
      call rone (h1,n)
      call rzero(h2,n)

      do j=0,nb                  ! Form the A matrix for basis function
         call axhelm(usave,ub(1,j),h1,h2,1,1)
         call axhelm(vsave,vb(1,j),h1,h2,1,1)
         if (ldim.eq.3) call axhelm(wsave,wb(1,j),h1,h2,1,1)
         do i=0,nb
            a0(i,j) = glsc2(ub(1,i),usave,n)+glsc2(vb(1,i),vsave,n)
            if (ldim.eq.3) a0(i,j) = a0(i,j)+glsc2(wb(1,i),wsave,n)
         enddo
      enddo

      if (nio.eq.0) write (6,*) 'exiting makea'

      return
      end
c-----------------------------------------------------------------------
      subroutine makeb0

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrvh/ h1(lt),h2(lt)
      common /scrns/ usave(lt),vsave(lt),wsave(lt)

      if (nio.eq.0) write (6,*) 'inside makeb'

      do j=0,nb
      do i=0,nb
         b0(i,j) = op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),
     $                         ub(1,j),vb(1,j),wb(1,j),bm1)
      enddo
      enddo

      if (nio.eq.0) write (6,*) 'exiting makeb'

      return
      end
c-----------------------------------------------------------------------
      subroutine makeic

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk1/ t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)
      common /scrk2/ h1(lt),h2(lt)

      if (nio.eq.0) write (6,*) 'inside makeic'

      call opcopy(t1,t2,t3,vx,vy,vz)

      n=lx1*ly1*lz1*nelt

      call opzero(vxlag,vylag,vzlag)
      call proj2bases(u,vx,vy,vz)

      do i=0,nb
         call opadds(
     $      vxlag,vylag,vzlag,ub(1,i),vb(1,i),wb(1,i),u(i,1),n,2)
      enddo

      call outpost(vx,vy,vz,pr,t,'ric')
      call outpost(vxlag,vylag,vzlag,pr,t,'ric')
      call opsub3(t1,t2,t3,vx,vy,vz,vxlag,vylag,vzlag)
      call outpost(t1,t2,t3,pr,t,'ric')
      el2 = op_glsc2_wt(t1,t2,t3,t1,t2,t3,bm1)
     $    / op_glsc2_wt(vx,vy,vz,vx,vy,vz,bm1)

      itest=1
      if (abs(el2-5.7993972158270270E-002).lt.1e-4) itest=0
      write (6,*) itest,el2,'ic_test'

    1 format('ic: ',i3,1p3e16.7)

      return
      end
c-----------------------------------------------------------------------
      subroutine evalc(cu)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real cu(nb)

      common /scrk4/ work(lx1*ly1*lz1*lelt)

      l=1

c     if (nio.eq.0) write (6,*) 'calling evalc'

      call rzero(cu,nb)

      do k=k0,k1
         uk=u(k,1)
         do j=j0,j1
            ujk=u(j,1)*uk
            do i=i0,i1
               cu(i)=cu(i)+clocal(l)*ujk
               l=l+1
            enddo
         enddo
      enddo

      call gop(cu,work,'+  ',nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_const
c This subroutine is solving rom with constrains
c The subroutine is based on BFGS method with barrier function
      
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

c     Matrices and vectors for advance
      real tmp(0:nb),tmat(nb,nb+1)
      real coef(1:nb)

      common /scrk3/ work(lt)
      common /scrk1/ t1(lt),t2(lt),t3(lt)

c     Variable for vorticity
      real vort(lt,3)

c     Working arrays for LU

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      n  = lx1*ly1*lz1*nelt

      time=time+ad_dt

      count = min0(ad_step,3)

      if (ad_step.le.3) then
         call cmult2(helm,b,ad_beta(1,count)/ad_dt,nb*nb)
         call add2s2(helm,a,1/ad_re,nb*nb)
      endif

      ONE = 1.
      ZERO= 0.

      call mxm(u,nb+1,ad_beta(2,count),3,tmp,1)

      call dgemv( 'N',nb,nb,ONE,b,nb,tmp(1),1,ZERO,rhs,1)

      call cmult(rhs,-1/ad_dt,nb)
      call add2s2(rhs,a0(1,0),-1/ad_re,nb)

      call copy(conv(1,3),conv(1,2),nb)
      call copy(conv(1,2),conv(1,1),nb)

      if (param(51).eq.0) then
         call mxm(c,nb*(nb+1),u,nb+1,tmat,1)
         call mxm(tmat,nb,u,nb+1,conv,1)
      else
         call evalc(conv)
      endif

      call mxm(conv,nb,ad_alpha(1,count),3,tmp,1)

      call sub2(rhs,tmp,nb)

      call copy(u(1,3),u(1,2),nb)
      call copy(u(1,2),u(1,1),nb)

      call opt_const

      if (mod(ad_step,ad_iostep).eq.0) then

!        This output is to make sure the ceof matches with matlab code

         call sleep(nid)

         write(6,*)'ad_step:',ad_step,ad_iostep,npp,nid

         if (ad_step.eq.ad_nsteps) then
            do j=1,nb
               write(6,*) 'final',j,u(j,1)
            enddo
         else
            do j=1,nb
               write(6,*) j,u(j,1)
            enddo
         endif
         call dumpcoef(u(:,1),nb,(ad_step/ad_iostep))

         call sleep(np-1-nid)

         call opzero(vx,vy,vz)
         do j=1,nb
            call opadds(vx,vy,vz,ub(1,j),vb(1,j),wb(1,j),coef(j),n,2)
         enddo
         call opadd2  (vx,vy,vz,ub,vb,wb)

!        comput the vorticity of the ROM reconstructed field
         call opcopy(t1,t2,t3,vx,vy,vz)
         call comp_vort3(vort,work1,work2,t1,t2,t3)
         ifto = .true. ! turn on temp in fld file
         call copy(t,vort,n)

         call outpost (vx,vy,vz,pr,t,'rom')
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine opt_const

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

c     parameter for barrier function
c     it should starting from value greater than one and decrease
      real B_qn(nb,nb), IBgf(nb), IBy(nb)
      real yIBy,sgf,sy,yBIgf
      real go(nb),fo,qndf
      integer par_step

      par_step = 3
      par = 1 


c     BFGS method with barrier function starts
      do k=1,par_step

c     use helm from BDF3/EXT3 as intial approximation
         do i=1,nb
            call copy(B_qn(1,i),helm(1,i),nb)
         enddo
         call comp_qnf
         call comp_qngradf
         call cmult(qngradf,-1,nb)

c     compute quasi-Newton step
         do j=1,500
            if (j==1) then
               call lu(flu,nb,nb,ir,ic)
               call copy(qns,qngradf,nb)
               call solve(qns,flu,1,nb,nb,ir,ic)
               call add2(u(1,1),qns,nb)
            else
               call copy(qns,qngradf,nb)
               call add2(u(1,1),qns,nb)
               
            endif
            call copy(go,gngraf,nb) ! store old qn-gradf
            call comp_qngradf ! update qn-gradf
            call sub3(qny,qngradf,go,nb)
c     BFGS update
            call copy(IBgf,gngraf,nb) ! comp B^-1 \nabla f
            call copy(IBy,gny,nb) ! compt B^-1 y
            sy = glsc2(qns,qny,nb) 
            

            fo = qnf ! store old qn-f
            call comp_qnf ! update qn-f
            qndf = abs(qnf-fo) 
         
c     update solution

c     update approximate Hessian

         enddo
         par = par*0.1

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_qngradf
      
      include 'SIZE'
      include 'MOR'
      real tmp1(nb),tmp2(nb),tmp3(nb),tmp4(nb)


      call sub3(tmp1,u(1,1),sample_max,nb)  
      call sub3(tmp2,u(1,1),sample_min,nb)  
      call add3(tmp3,tmp1,tmp2,nb)

      call add3s12(qngradf,rhs,tmp3,-1,-par,nb)

      ONE = 1.
      ZERO= 0.
      call dgemv( 'N',nb,nb,ONE,helm,nb,u(1,1),1,ZERO,tmp4,1)
      call add2(qngradf,tmp4,nb)

      return 
      end
c-----------------------------------------------------------------------
      subroutine comp_qnf
      
      include 'SIZE'
      include 'MOR'
      real tmp1(nb),tmp2(nb),tmp3(nb),tmp4(nb),tmp5(nb)
      real term1,term2,term3,term4

c     bar1 and bar2 are the barrier function for two constrains
      real bar1,bar2

      call sub3(tmp1,sample_max,u(1,1),nb)  
      call sub3(tmp2,u(1,1),sample_min,nb)  

c     currently can only come up with this way to compute log for an array
      do i=1,nb
         tmp3(i) = log(tmp1(i))
         tmp4(i) = log(tmp2(i))
      enddo

      bar1 = vlsum(tmp3,nb)
      bar2 = vlsum(tmp4,nb)

c     evaluate quasi-newton f
c     inverse of helm hasn't been implemented yet

      term2 = glsc2(u(1,1),rhs,nb)
      term3 = par*(bar1+bar2)

      ONE = 1.
      ZERO= 0.
      call dgemv( 'N',nb,nb,ONE,helm,nb,u(1,1),1,ZERO,tmp5,1)
      term1 = 0.5 * glsc2(tmp5,u(1,1),nb)

      qnf = term1 - term2 - term3

      return
      end
