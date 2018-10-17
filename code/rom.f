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
      subroutine readops

      include 'SIZE'
      include 'TOTAL'
      include 'POD'

      parameter (lt=lx1*ly1*lz1*lelt)

      real usave(lt,ms),vsave(lt,ms),wsave(lt,ms)
      real u0(lt,3)

      real uu(ms,ms),Identity(ms,ms),eig(ms),eigv(ms,ms),w(ms,ms)
      real uw(lt),vw(lt),ww(lt),h1(lt),h2(lt)
      real u0r(ms)

      real evec(ms,nb)

      character(len=10) fname

      n  = lx1*ly1*lz1*nelt
      ns = ms

      call rzero(vz,n)
      call rzero(wb,n)

      call opcopy   (u0(1,1),u0(1,2),u0(1,3),vx,vy,vz)
      call opcopy   (ub(1,0),vb(1,0),wb(1,0),vx,vy,vz)

      call get_saved_fields(usave,vsave,wsave,ns,u0)

      call rone (h1,n)
      call rzero(h2,n)
      call rzero(Identity,ms*ms)

      do j=1,ns                    ! Form the Gramian, U=U_K^T A U_K
         call axhelm(uw,usave(1,j),h1,h2,1,1)
         call axhelm(vw,vsave(1,j),h1,h2,1,1)
         if (ldim.eq.3) call axhelm(ww,wsave(1,j),h1,h2,1,1)
         Identity(j,j) = 1.
         do i=1,ns
            uu(i,j) = glsc2(usave(1,i),uw,n)+glsc2(vsave(1,i),vw,n)
            if (ldim.eq.3) uu(i,j) = uu(i,j)+glsc2(wsave(1,i),ww,n)
         enddo
      enddo

      call readeig(evec)

      call get_saved_fields(usave,vsave,wsave,ns,u0)

      ONE = 1.
      ZERO= 0.

      ! ub, vb, wb, are the modes
      call dgemm( 'N','N',n,nb,ms,ONE,usave,lt,evec,ms,ZERO,ub(1,1),lt)
      call dgemm( 'N','N',n,nb,ms,ONE,vsave,lt,evec,ms,ZERO,vb(1,1),lt)
      if (ldim.eq.3)
     $call dgemm( 'N','N',n,nb,ms,ONE,wsave,lt,evec,ms,ZERO,wb(1,1),lt)

      call opcopy(u0(1,1),u0(1,2),u0(1,3),vx,vy,vz)
      call opcopy(ub(1,0),vb(1,0),wb(1,0),vx,vy,vz)

      do i=0,nb
         call outpost(ub(1,i),vb(1,i),wb(1,i),pr,t,'bas')
      enddo

      call readc0(c0,(nb+1)**3)
      call readab(a0,b0,(nb+1)**2)
      call readic(ic,nb+1)

      if (np.gt.1) call makecloc

      call rom_setup

      return
      end
c-----------------------------------------------------------------------
      subroutine genops

      include 'SIZE'
      include 'TOTAL'
      include 'POD'

      parameter (lt=lx1*ly1*lz1*lelt)

      real usave(lt,ms),vsave(lt,ms),wsave(lt,ms)
      real uu(ms,ms),Identity(ms,ms),eig(ms),eigv(ms,ms),w(ms,ms)
      real uw(lt),vw(lt),ww(lt),h1(lt),h2(lt)
      real u0(lt,3)
      real u0r(ms)

      real evec(ms,nb)

      character(len=10) fname

      n  = lx1*ly1*lz1*nelt
      ns = ms

      call rzero(vz,n)
      call rzero(wb,n)

      call opcopy   (u0(1,1),u0(1,2),u0(1,3),vx,vy,vz)
      call opcopy   (ub(1,0),vb(1,0),wb(1,0),vx,vy,vz)

      call get_saved_fields(usave,vsave,wsave,ns,u0)

      call rone (h1,n)
      call rzero(h2,n)
      call rzero(Identity,ms*ms)

      do j=1,ns                    ! Form the Gramian, U=U_K^T A U_K
         call axhelm(uw,usave(1,j),h1,h2,1,1)
         call axhelm(vw,vsave(1,j),h1,h2,1,1)
         if (ldim.eq.3) call axhelm(ww,wsave(1,j),h1,h2,1,1)
         do i=1,ns
            uu(i,j) = glsc2(usave(1,i),uw,n)+glsc2(vsave(1,i),vw,n)
            if (ldim.eq.3) uu(i,j) = uu(i,j)+glsc2(wsave(1,i),ww,n)
            if (i == j) Identity(i,j) = 1
         enddo
         if (nio.eq.0) write(6,*) j,uu(1,j),' uu'
      enddo

      call generalev(uu,Identity,eig,ms,w)
      call copy(eigv,uu,ms*ms)
      eig = eig(ms:1:-1)

      nvecs = nb
      if (nio.eq.0) write(6,*)'number of mode:',nb

      do l = 1,nvecs
         call copy(evec(1,l),eigv(1,ms-l+1),ms)
      enddo

      ONE = 1.
      ZERO= 0.

      ! ub, vb, wb, are the modes
      call dgemm( 'N','N',n,nb,ms,ONE,usave,lt,evec,ms,ZERO,ub(1,1),lt)
      call dgemm( 'N','N',n,nb,ms,ONE,vsave,lt,evec,ms,ZERO,vb(1,1),lt)
      if (ldim.eq.3)
     $call dgemm( 'N','N',n,nb,ms,ONE,wsave,lt,evec,ms,ZERO,wb(1,1),lt)

      do i=0,nb
         call outpost(ub(1,i),vb(1,i),wb(1,i),pr,t,'bas')
c        call outmatl(ub(1,i),vb(1,i),wb(1,i),i)
      enddo

      call makec
      call makecloc
      call get_a_b
      call makeic

      call rom_setup

      if (nio.eq.0) write (6,*) 'end in my_pod4$',ns

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_setup

      include 'SIZE'
      include 'TOTAL'
      include 'POD'

      parameter (lt=lx1*ly1*lz1*lelt)

c     Matrices and vectors for advance
      real helm(1:nb,1:nb), rhs(1:nb)
      real tmp(0:nb),tmat(nb,nb+1)
      real coef(1:nb)

c     Working arrays for LU

      time = 0.

      do i=1,nb
         call copy(a(1,i),a0(1,i),nb)
         call copy(b(1,i),b0(1,i),nb)
      enddo

      do i=0,nb
      do j=0,nb
         call copy(c(1,j,i),c0(1,j,i),nb)
      enddo
      enddo

      n  = lx1*ly1*lz1*nelt

      call rzero(u,(nb+1)*3)
      u(0,2) = 1.
      u(0,3) = 1.
      call copy(u,ic,nb+1)

      ad_nsteps=nsteps
      ad_iostep=iostep

      ad_dt = dt
      ad_re = 1/param(2)

      ! BDFk/EXTk coefficients ( will change to BD inside Nek)
      call compute_BDF_coef(ad_alpha,ad_beta)

      if (nio.eq.0) write(6,*)'ad_step:',ad_step,ad_iostep,npp

      do j=1,nb
         if (nio.eq.0) write(6,*) j,u(j,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_step

      include 'SIZE'
      include 'TOTAL'
      include 'POD'

      parameter (lt=lx1*ly1*lz1*lelt)

c     Matrices and vectors for advance
      real helm(1:nb,1:nb), rhs(1:nb)
      real tmp(0:nb),tmat(nb,nb+1)
      real coef(1:nb)

      common /scrk3/ work(lt)

c     Working arrays for LU

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      n  = lx1*ly1*lz1*nelt

      time=time+ad_dt

      count = min0(ad_step,3)

      if (ad_step.le.3) then
         call cmult2(flu,b,ad_beta(1,count)/ad_dt,nb*nb)
         call add2s2(flu,a,1/ad_re,nb*nb)
      endif

      ONE = 1.
      ZERO= 0.

      call mxm(u,nb+1,ad_beta(2,count),3,tmp,1)

      call dgemv( 'N',nb,nb,ONE,b,nb,tmp(1),1,ZERO,rhs,1)

      call cmult(rhs,-1/ad_dt,nb)
      call add2s2(rhs,a0(1,0),-1/ad_re,nb)

      call copy(conv(1,3),conv(1,2),nb)
      call copy(conv(1,2),conv(1,1),nb)

      if (npp.eq.1) then
         call mxm(c,nb*(nb+1),u,nb+1,tmat,1)
         call mxm(tmat,nb,u,nb+1,conv,1)
      else
         call evalc(conv)
      endif

      call mxm(conv,nb,ad_alpha(1,count),3,tmp,1)

      if (ad_step.eq.1.or.ad_step.eq.2.or.ad_step.eq.ad_nsteps) then
      do i=0,nb-1
         if (nio.eq.0) write (6,*) '2 evalc',tmp(i)
      enddo
      endif

      call sub2(rhs,tmp,nb)

      if (ad_step.le.3) call lu(flu,nb,nb,irr,icc)

      call solve(rhs,flu,1,nb,nb,irr,icc)

      call copy(u(1,3),u(1,2),nb)
      call copy(u(1,2),u(1,1),nb)

c     if (npp.ne.1) call gop(rhs,work,'+  ',nb)

      call copy(u(1,1),rhs,nb)
      call copy(coef,rhs,nb)

c     if (mod(ad_step,ad_iostep)==0) then
      if (ad_step.eq.ad_nsteps.or.ad_step.eq.1.or.ad_step.eq.2) then

!        This output is to make sure the ceof matches with matlab code

         call sleep(nid)

         write(6,*)'ad_step:',ad_step,ad_iostep,npp,nid

         do j=1,nb
            write(6,*) j,u(j,1)
         enddo

         call sleep(np-1-nid)

         call opzero(vx,vy,vz)
         do j=1,nb
            call opadds(vx,vy,vz,ub(1,j),vb(1,j),wb(1,j),coef(j),n,2)
         enddo
         call opadd2  (vx,vy,vz,ub,vb,wb)
         call outpost (vx,vy,vz,pr,t,'rom')
      endif

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
      subroutine outmatl(uu,v,w,k)

      include 'SIZE'
      include 'TOTAL'
      include 'POD'

      parameter (lt=lx1*ly1*lz1*lelt)
      real uu(lt),v(lt),w(lt)
      character*7 fname

      if (np.gt.1) call exitti('outmatl works for P=1 only!$',np)

      write(fname,22) k
   22 format('out.',i3.3)
      open(unit=33,file=fname)

      n = lx1*ly1*lz1*nelt
      do i=1,n
         write(33,33) xm1(i,1,1,1),ym1(i,1,1,1),uu(i),v(i)
   33    format(1p4e16.7)
      enddo
      close(33)

      return
      end
c-----------------------------------------------------------------------
      subroutine makec

      include 'SIZE'
      include 'TOTAL'
      include 'POD'

      parameter (lt=lx1*ly1*lz1*lelt)

      real cux(lt), cuy(lt), cuz(lt)

      common /scrk1/ t1(lt), binv(lt)

      call invers2(binv,bm1,lx1*ly1*lz1*nelv)
      call rone(binv,lx1*ly1*lz1*nelv)

      do k=0,nb
         call setcnv_c(ub(1,k),vb(1,k),wb(1,k))
         do j=0,nb
            call setcnv_u(ub(1,j),vb(1,j),wb(1,j))
            call ccu(cux,cuy,cuz)
            do i=0,nb
               c0(i,j,k) = op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),
     $                                cux,cuy,cuz,binv)
               if (nid.eq.0) write (6,*) 'c0',i,j,k,c0(i,j,k)
            enddo
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine makecloc

      include 'SIZE'
      include 'TOTAL'
      include 'POD'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (m=5)

      real cux(lt), cuy(lt), cuz(lt)
      integer mps(lt),mqs(lt),mrs(lt)

      real vr(lt), pts(lt), ur(lt)
      integer vl(lt),vi(m,lt),ui(m,lt)

      common /scrk1/ t1(lt), binv(lt)
      integer crh, parts(lt)
      common /romi/ crh
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      ntot=nb*(nb+1)*(nb+1)

      nblock=nb*(nb+1)*(nb+1)/npp+1

      call fgslib_crystal_setup(crh,nekcomm,np)

      call izero(parts,lt)
      call setpart(parts,npp,ntot)

      call izero(mps,lt)
      call izero(mqs,lt)
      call izero(mrs,lt)

      call factor3(mp,mq,mr,npp)
      if (nio.eq.0) write (6,*) 'factor3',mp,mq,mr,npp

      call setpart3(mps,mqs,mrs,mp,mq,mr,nb)
      if (nio.eq.0) then
         write (6,*) 'setpart3'
         do i=1,mp
            write (6,*) 'mps',mps(i)
         enddo
         do i=1,mq
            write (6,*) 'mqs',mqs(i)
         enddo
         do i=1,mr
            write (6,*) 'mrs',mrs(i)
         enddo
      endif

      call setrange(mps,mqs,mrs,mp,mq,mr) ! set index range i0,i1,j0,j1,k0,k1

      kp=1
      n=1
      ip=0
      l=0

      do k=0,nb
      do j=0,nb
      do i=1,nb
         l=l+1
         vi(1,l)=mod(ip+nid+1,npp)
         vi(2,l)=ijk2pid(i,j,k,mps,mqs,mrs,mp,mq,ms)

         vi(3,l)=i
         vi(4,l)=j
         vi(5,l)=k

         if (nio.eq.0) write (6,*) 'ijk2pid',i,j,k,vi(2,l)
         vr(  l)=c0(i,j,k)
         if (l.eq.nblock.or.(k.eq.nb.and.j.eq.nb.and.i.eq.nb)) then
            call fgslib_crystal_tuple_transfer
     $         (crh,l,lt,vi,m,vl,0,vr,1,kp)
            if (nid.eq.ip) then
               call copy(ctmp,vr,lt)
               call icopy(ui,vi,lt*m)
               inums=l
            endif
            l=0
            ip=ip+1
         endif
      enddo
      enddo
      enddo

      call sleep(nid)

      write (6,*) '1 nid,npp',nid,npp,inums

      do i=1,inums
         write (6,*) nid,ctmp(i),ui(2,i)
      enddo

      call sleep(npp-nid-1)

      kp=2

      call fgslib_crystal_tuple_transfer
     $   (crh,inums,lt,ui,m,vl,0,ctmp,1,kp)

      call sleep(nid)

      write (6,*) '2 nid,npp',nid,npp,inums

      do i=1,inums
         write (6,*) nid,ctmp(i),ui(2,i),ui(3,i),ui(4,i),ui(5,i)
      enddo

      call sleep(npp-nid-1)

      call fgslib_crystal_free(crh)

      do i=1,inums
         call ijk2l(l,ui(3,i),ui(4,i),ui(5,i))
         clocal(l) = ctmp(i)
      enddo

      call sleep(nid)

      do i=1,inums
         write (6,*) nid,clocal(i)
      enddo

      call sleep(npp-nid)

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
      subroutine get_A_B

      include 'SIZE'
      include 'TOTAL'
      include 'POD'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrvh/ h1(lt),h2(lt)
      common /scrns/ usave(lt),vsave(lt),wsave(lt)

      do j=0,nb
      do i=0,nb
         b0(i,j) = op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),
     $                         ub(1,j),vb(1,j),wb(1,j),bm1)
      enddo
      enddo

      n= lx1*ly1*lz1*nelt
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

      return
      end
c-----------------------------------------------------------------------
      subroutine h1prod(prod,u1,v1,w1,u2,v2,w2)

      include 'SIZE'

      parameter(lt=lx1*ly1*lz1*lelt)

      common /scrk1/ au(lt), av(lt), aw(lt)
      common /scrk2/ h1(lt), h2(lt)

      real u1(lt), v1(lt), w1(lt)
      real u2(lt), v2(lt), w2(lt)

      n = lx1*ly1*lz1*nelt

      call rone(h1,n)
      call rzero(h2,n)

      call axhelm(au,u1,h1,h2,1,1)
      call axhelm(av,v1,h1,h2,1,1)

      p1=glsc2(au,u2,n)
      p2=glsc2(av,v2,n)
      prod=p1+p2

      if (ldim.eq.3) then
         call axhelm(aw,w1,h1,h2,1,1)
         prod = prod + glsc2(aw,w2,n)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine makeic

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'POD'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk1/ t1(lt),t2(lt),t3(lt)

      if (nio.eq.0) write (6,*) 'inside makeic'

      call opcopy(t1,t2,t3,vx,vy,vz)

      ic(0) = 1.

      call opsub2(t1,t2,t3,ub(1,0),vb(1,0),wb(1,0))

      do i=1,nb
         t4 = op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),t1,t2,t3,bm1)
         t5 = op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),
     $                    ub(1,i),vb(1,i),wb(1,i),bm1)
c        call h1prod(t4,ub(1,i),vb(1,i),wb(1,i),t1,t2,t3)
c        call h1prod(t5,ub(1,i),vb(1,i),wb(1,i),ub(1,i),vb(1,i),wb(1,i))
         ic(i) = t4 / t5
         if (nio.eq.0) write (6,*) 'find coef: ',i,t4,t5,ic(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine evalc(cu)

      include 'SIZE'
      include 'TOTAL'
      include 'POD'

      real cu(nb)

      common /scrk4/ work(lx1*ly1*lz1*lelt)

      l=1

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

      if (ad_step.eq.1.or.ad_step.eq.2.or.ad_step.eq.ad_nsteps) then
      call sleep(nid)

      do i=1,nb
         write (6,*) '1 evalc',cu(i)
      enddo

      call sleep(np-nid-1)
      endif

      call gop(cu,work,'+  ',nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine readc0(c0,n)

      include 'SIZE'
      include 'PARALLEL'

      real c0(n)

      if (nid.eq.(np-1)) then
         open (unit=12,file='cten')
         read (12,*) (c0(k),k=1,n)
         close (unit=12)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine readab(a0,b0,n)

      include 'SIZE'

      real a0(n),b0(n)

      common /scrk5/ w(lx1*ly1*lz1*lelt)

      call rzero(a0,n)
      call rzero(b0,n)

      if (nid.eq.0) then
         open (unit=12,file='amat')
         read (12,*) (a0(k),k=1,n)
         close (unit=12)

         open (unit=12,file='bmat')
         read (12,*) (b0(k),k=1,n)
         close (unit=12)
      endif

      call gop(a0,w,'+  ',n)
      call gop(b0,w,'+  ',n)

      return
      end
c-----------------------------------------------------------------------
      subroutine readic(ic,n)

      include 'SIZE'

      real ic(n)

      open (unit=12,file='ic')
      read (12,*) (ic(k),k=1,n)
      close (unit=12)

      return
      end
c-----------------------------------------------------------------------
