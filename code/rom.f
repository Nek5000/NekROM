c-----------------------------------------------------------------------
      include 'pod.f'
c-----------------------------------------------------------------------
      subroutine readops

      include 'SIZE'
      include 'TOTAL'
      include 'POD'

      if (nio.eq.0) write (6,*) 'inside readops'

      call readab(a0,b0,(nb+1)**2)
      call readc0(c0,(nb+1)**3)
      call readic(u,nb+1)

      if (np.gt.1) call makecloc

      if (nio.eq.0) write (6,*) 'exiting readops'

      return
      end
c-----------------------------------------------------------------------
      subroutine genops

      include 'SIZE'
      include 'TOTAL'
      include 'POD'

      if (nio.eq.0) write (6,*) 'inside genops'

      call get_a_b
      call makec
      call makecloc
      call makeic

      if (nio.eq.0) write (6,*) 'exiting genops'

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_setup

      include 'SIZE'
      include 'TOTAL'
      include 'POD'

      if (nid.eq.0) write (6,*) 'inside rom_setup'

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

      if (nid.eq.0) write (6,*) 'exiting rom_setup'

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

c     if (ad_step.eq.1.or.ad_step.eq.2.or.ad_step.eq.ad_nsteps) then
c     do i=0,nb-1
c        if (nio.eq.0) write (6,*) '2 evalc',tmp(i)
c     enddo
c     endif

      call sub2(rhs,tmp,nb)

      if (ad_step.le.3) call lu(flu,nb,nb,irr,icc)

      call solve(rhs,flu,1,nb,nb,irr,icc)

      call copy(u(1,3),u(1,2),nb)
      call copy(u(1,2),u(1,1),nb)

c     if (npp.ne.1) call gop(rhs,work,'+  ',nb)

      call copy(u(1,1),rhs,nb)
      call copy(coef,rhs,nb)

      if (mod(ad_step,ad_iostep).eq.0) then

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
c              if (nid.eq.0) write (6,*) 'c0',i,j,k,c0(i,j,k)
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
c     if (nio.eq.0) then
c        write (6,*) 'setpart3'
c        do i=1,mp
c           write (6,*) 'mps',mps(i)
c        enddo
c        do i=1,mq
c           write (6,*) 'mqs',mqs(i)
c        enddo
c        do i=1,mr
c           write (6,*) 'mrs',mrs(i)
c        enddo
c     endif

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

c        if (nio.eq.0) write (6,*) 'ijk2pid',i,j,k,vi(2,l)
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

c     write (6,*) '1 nid,npp',nid,npp,inums

c     do i=1,inums
c        write (6,*) nid,ctmp(i),ui(2,i)
c     enddo

      call sleep(npp-nid-1)

      kp=2

      call fgslib_crystal_tuple_transfer
     $   (crh,inums,lt,ui,m,vl,0,ctmp,1,kp)

      call sleep(nid)

c     write (6,*) '2 nid,npp',nid,npp,inums

c     do i=1,inums
c        write (6,*) nid,ctmp(i),ui(2,i),ui(3,i),ui(4,i),ui(5,i)
c     enddo

      call sleep(npp-nid-1)

      call fgslib_crystal_free(crh)

      do i=1,inums
         call ijk2l(l,ui(3,i),ui(4,i),ui(5,i))
         clocal(l) = ctmp(i)
      enddo

      call sleep(nid)

c     do i=1,inums
c        write (6,*) nid,clocal(i)
c     enddo

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
      subroutine makeic

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'POD'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk1/ t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)
      common /scrk2/ h1(lt),h2(lt)

      if (nio.eq.0) write (6,*) 'inside makeic'

      call opcopy(t1,t2,t3,vx,vy,vz)

      n=lx1*ly1*lz1*nelt

      call rone(h1,n)
      call rzero(h2,n)

      u(0,1) = 1.

      call opsub3(t1,t2,t3,vx,vy,vz,ub(1,0),vb(1,0),wb(1,0))
      call h10proj(u(1,1),t1,t2,t3)
      call opzero(vxlag,vylag,vzlag)

      ii=3

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
      if (abs(el2-3.9525554652664913E-002).lt.1e-4) itest=0
      write (6,*) itest,el2,'ic_test'

    1 format('ic: ',i3,1p3e16.7)

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

c     do i=1,nb
c        write (6,*) '1 evalc',cu(i)
c     enddo

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
