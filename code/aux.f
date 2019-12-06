c-----------------------------------------------------------------------
      subroutine opadd3(a1,a2,a3,b1,b2,b3,c1,c2,c3)

      include 'SIZE'

      real a1(1),a2(1),a3(1),b1(1),b2(1),b3(1)
      real c1(1),c2(1),c3(1)

      ntot1=lx1*ly1*lz1*nelv
      call add3(a1,b1,c1,ntot1)
      call add3(a2,b2,c2,ntot1)
      if (ldim.eq.3) call add3(a3,b3,c3,ntot1)

      return
      end
c-----------------------------------------------------------------------
      subroutine recont(tt,coef)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real tt(lt),coef(0:nb)

      n=lx1*ly1*lz1*nelt

      call rzero(tt,n)

      do i=0,nb
         call add2s2(tt,tb(1,i),coef(i),n)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine recont_rms(tt)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real tt(lt)

      n=lx1*ly1*lz1*nelt

      call rzero(tt,n)

      do j=0,nb
      do i=0,nb
         call col3(tbt,tb(1,i),tb(1,j),n)
         call add2s2(tt,tbt,ut2a(1+i+(nb+1)*j),n)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine reconu_rms(ux,uy,uz)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ux(lt),uy(lt),uz(lt)

      n=lx1*ly1*lz1*nelv

      call opzero(ux,uy,uz)

      do j=0,nb
         if (nio.eq.0) write (6,*) 'reconu_rms:',j,'/',nb
         do i=0,nb
            call admcol3(ux,ub(1,i),ub(1,j),u2a(1+i+(nb+1)*j),n)
            call admcol3(uy,vb(1,i),vb(1,j),u2a(1+i+(nb+1)*j),n)
            if (ldim.eq.3)
     $         call admcol3(uz,wb(1,i),wb(1,j),u2a(1+i+(nb+1)*j),n)
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine reconv(ux,uy,uz,coef)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ux(lt),uy(lt),uz(lt),coef(0:nb)

      n=lx1*ly1*lz1*nelv

      call opzero(ux,uy,uz)

      do i=0,nb
         call add2s2(ux,ub(1,i),coef(i),n)
         call add2s2(uy,vb(1,i),coef(i),n)
         call add2s2(uz,wb(1,i),coef(i),n)
c        call opadds(ux,uy,uz,ub(1,i),vb(1,i),wb(1,i),coef(i),n,2)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ctke_fom(tke,ux,uy,uz,uavg,vavg,wavg)

      include 'SIZE'
      include 'MOR'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrns/ ud(lt),vd(lt),wd(lt)

      real ux(lt),uy(lt),uz(lt),uavg(lt),vavg(lt),wavg(lt)

      call opsub3(ud,vd,wd,ux,uy,uz,uavg,vavg,wavg)
      tke = .5*op_glsc2_wt(ud,vd,wd,ud,vd,wd,bm1)

      return
      end
c-----------------------------------------------------------------------
      subroutine ctke

      include 'SIZE'
      include 'MOR'
      include 'TSTEP'

      common /ctkea/ cdiff(0:lb),ccat(0:lb)

      parameter (lt=lx1*ly1*lz1*lelt)

      tke=0.

      do i=0,nb
         cdiff(i)=u(i)-uas(i)
      enddo

      call mxm(bu0,nb+1,cdiff,nb+1,ccat,1)
      tke=.5*vlsc2(ccat,cdiff,nb+1)

      if (nio.eq.0) write (6,*) time,tke,'tke'

      return
      end
c-----------------------------------------------------------------------
      subroutine lints(s1,s2,l)

      character*1 s1(1)
      character*1 s2(1)

      len=ltruncr(s2,l)
      call blank(s1,l)
      call chcopy(s1,s2,len)

      return
      end
c-----------------------------------------------------------------------
      function ltruncr(string,l)

      character*1 string(l)
      character*1   blnk
      data blnk/' '/

      do i=1,l
         l1=i-1
         if (string(i).eq.blnk) goto 200
      enddo
      l1=0

  200 continue
      ltruncr=l1

      return
      end
c-----------------------------------------------------------------------
      subroutine shift3(u,v,n)

      real u(n,3),v(n)

      call copy(u(1,3),u(1,2),n)
      call copy(u(1,2),u(1,1),n)
      call copy(u(1,1),v,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine lap2d(d2u,u)

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)

      real d2u(1),u(1)

      common /scrl2d/ ux(lt),uy(lt),uxx(lt),uyy(lt),t1(lt),t2(lt)

      call gradm1(ux,uy,t1,u)
      call dsavg(ux)
      call dsavg(uy)

      call gradm1(uxx,t1,t1,ux)
      call gradm1(t1,uyy,t2,uy)

      call add3(d2u,uxx,uyy,lx1*ly1*lz1*nelv)
      call dsavg(d2u)

      return
      end
c-----------------------------------------------------------------------
      subroutine push_op(vx,vy,vz)

      include 'SIZE'

      parameter (lt1=lx1*ly1*lz1*lelt)

      common /pushpop/ ux(lt1),uy(lt1),uz(lt1)

      real vx(lt1),vy(lt1),vz(lt1)

      call opcopy(ux,uy,uz,vx,vy,vz)

      return
      end
c-----------------------------------------------------------------------
      subroutine pop_op(vx,vy,vz)

      include 'SIZE'

      parameter (lt1=lx1*ly1*lz1*lelt)

      common /pushpop/ ux(lt1),uy(lt1),uz(lt1)

      real vx(lt1),vy(lt1),vz(lt1)

      call opcopy(vx,vy,vz,ux,uy,uz)

      return
      end
c-----------------------------------------------------------------------
      subroutine push_sol(vx,vy,vz,pr,t)

      include 'SIZE'

      parameter (lt1=lx1*ly1*lz1*lelt)
      parameter (lt2=lx2*ly2*lz2*lelt)

      common /pushpop/ ux(lt1),uy(lt1),uz(lt1),pp(lt2),tt(lt1,ldimt)

      real vx(lt1),vy(lt1),vz(lt1),pr(lt2),t(lt1,ldimt)

      call opcopy(ux,uy,uz,vx,vy,vz)
      call copy(pp,pr,lx2*ly2*lz2*nelv)

      do idim=1,ldimt
         call copy(tt(1,idim),t(1,idim),lx1*ly1*lz1*nelt)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine pop_sol(vx,vy,vz,pr,t)

      include 'SIZE'

      parameter (lt1=lx1*ly1*lz1*lelt)
      parameter (lt2=lx2*ly2*lz2*lelt)

      common /pushpop/ ux(lt1),uy(lt1),uz(lt1),pp(lt2),tt(lt1,ldimt)

      real vx(lt1),vy(lt1),vz(lt1),pr(lt2),t(lt1,ldimt)

      call opcopy(vx,vy,vz,ux,uy,uz)
      call copy(pr,pp,lx2*ly2*lz2*nelv)

      do idim=1,ldimt
         call copy(t(1,idim),tt(1,idim),lx1*ly1*lz1*nelt)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine zero_sol(vx,vy,vz,pr,t)

      include 'SIZE'

      parameter (lt1=lx1*ly1*lz1*lelt)
      parameter (lt2=lx2*ly2*lz2*lelt)

      real vx(lt1),vy(lt1),vz(lt1),pr(lt2),t(lt1,ldimt)

      call opzero(vx,vy,vz)
      call rzero(pr,lx2*ly2*lz2*nelv)

      do idim=1,ldimt
         call rzero(t(1,idim),lx1*ly1*lz1*nelt)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine add_sol(vx,vy,vz,pr,t,ux,uy,uz,pp,tt)

      include 'SIZE'

      parameter (lt1=lx1*ly1*lz1*lelt)
      parameter (lt2=lx2*ly2*lz2*lelt)

      real vx(lt1),vy(lt1),vz(lt1),pr(lt2),t(lt1,ldimt)
      real ux(lt1),uy(lt1),uz(lt1),pp(lt2),tt(lt1,ldimt)

      call opadd2(vx,vy,vz,ux,uy,uz)
      call add2(pr,pp,lx2*ly2*lz2*nelv)

      do idim=1,1
         call add2(t(1,idim),tt(1,idim),lx1*ly1*lz1*nelt)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine scale_sol(vx,vy,vz,pr,t,s)

      include 'SIZE'

      parameter (lt1=lx1*ly1*lz1*lelt)
      parameter (lt2=lx2*ly2*lz2*lelt)

      real vx(lt1),vy(lt1),vz(lt1),pr(lt2),t(lt1,ldimt)

      call opcmult(vx,vy,vz,s)
      call cmult(pr,s,lx2*ly2*lz2*nelv)

      do idim=1,ldimt
         call cmult(t(1,idim),s,lx1*ly1*lz1*nelt)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_sol(vx,vy,vz,pr,t,ux,uy,uz,pp,tt)

      include 'SIZE'

      parameter (lt1=lx1*ly1*lz1*lelt)
      parameter (lt2=lx2*ly2*lz2*lelt)

      real vx(lt1),vy(lt1),vz(lt1),pr(lt2),t(lt1,ldimt)
      real ux(lt1),uy(lt1),uz(lt1),pp(lt2),tt(lt1,ldimt)

      call opcopy(vx,vy,vz,ux,uy,uz)
      call copy(pr,pp,lx2*ly2*lz2*nelv)
      call copy(t,tt,lx1*ly1*lz1*nelv)

c     do idim=1,ldimt
c        call copy(t(1,idim),tt(1,idim),lx1*ly1*lz1*nelt)
c     enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine conv_sol(ux,uy,uz,vx,vy,vz) ! compute convection

      include 'SIZE'
      include 'INPUT'
      include 'MASS'

      parameter (lt1=lx1*ly1*lz1*lelt)
      parameter (lt2=lx2*ly2*lz2*lelt)

      real vx(lt1),vy(lt1),vz(lt1),pr(lt2),t(lt1,ldimt)
      real ux(lt1),uy(lt1),uz(lt1)
      common /scrctd/ t1(lt1),t2(lt1),t3(lt1),h1(lt1),h2(lt1)

      n1=lx1*ly1*lz1*nelt
      n2=lx2*ly2*lz2*nelt

      call setcnv_c(vx,vy,vz)
      call setcnv_u(vx,vy,vz)
      call ccu(ux,uy,uz)
      call opbinv1(ux,uy,uz,ux,uy,uz,1.)

      return
      end
c-----------------------------------------------------------------------
      subroutine ctd_sol(ux,uy,uz,vx,vy,vz,pr,t) ! compute the time-derivative

      include 'SIZE'
      include 'INPUT'
      include 'MASS'

      parameter (lt1=lx1*ly1*lz1*lelt)
      parameter (lt2=lx2*ly2*lz2*lelt)

      real vx(lt1),vy(lt1),vz(lt1),pr(lt2),t(lt1,ldimt)
      real ux(lt1),uy(lt1),uz(lt1)
      common /scrctd/ t1(lt1),t2(lt1),t3(lt1),h1(lt1),h2(lt1)

      n1=lx1*ly1*lz1*nelt
      n2=lx2*ly2*lz2*nelt

      call gradp(ux,uy,uz,pr)
      call opchsgn(ux,uy,uz)
      call conv_sol(t1,t2,t3,vx,vy,vz)
      call opchsgn(t1,t2,t3)

      call opadd2(ux,uy,uz,t1,t2,t3)

      call rone(h1,n1)
      call rzero(h2,n1)

      call lap2d(t1,vx)
      call lap2d(t2,vy)
      if (ldim.eq.3) call lap2d(t3,vz)
      s=-param(2)
      call opcmult(t1,t2,t3,s)

      call opadd2(ux,uy,uz,t1,t2,t3)

      return
      end
c-----------------------------------------------------------------------
      subroutine gradp(px,py,pz,pr)

      include 'SIZE'

      parameter (lt1=lx1*ly1*lz1*lelt)
      parameter (lt2=lx2*ly2*lz2*lelt)

      real pr(lt2),px(lt1),py(lt1),pz(lt1)

      call opgradt(px,py,pz,pr)
      call opbinv1(px,py,pz,px,py,pz,1.)

      return
      end
c-----------------------------------------------------------------------
      subroutine setdps
      call exitti('called deprecated subroutine setdps$',1)
      return
      end
c-----------------------------------------------------------------------
      subroutine setconvbases ! deprecated
      call exitti('called deprecated subroutine setconvbases$',1)
      return
      end
c-----------------------------------------------------------------------
      subroutine setdtbases
      call exitti('called deprecated subroutine setdtbases$',1)
      return
      end
c-----------------------------------------------------------------------
      subroutine comp_hyperpar

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ep
      real uw(lt),vw(lt),ww(lt)
      real tmp1(nb),tmp2(nb),delta(nb)
      real work(ls,nb)

      ! eps is the free parameter
      ! 1e-2 is used in the paper
      ep = 1e-2

      n  = lx1*ly1*lz1*nelt

      do j=1,nb                    ! compute hyper-parameter
         call axhelm(uw,ub(1,j),ones,zeros,1,1)
         call axhelm(vw,vb(1,j),ones,zeros,1,1)
         if (ldim.eq.3) call axhelm(ww,wb(1,j),ones,zeros,1,1)
         do i=1,ls
            work(i,j) = glsc2(us0(1,1,i),uw,n)+glsc2(us0(1,2,i),vw,n)
            if (ldim.eq.3) work(i,j)=work(i,j)+glsc2(us0(1,ldim,i),ww,n)
         enddo
         tmp1(j) = vlmin(work(1,j),ls)
         tmp2(j) = vlmax(work(1,j),ls)
         delta(j) = tmp2(j)-tmp1(j)
         umin(j) = tmp1(j) - ep * delta(j)
         umax(j) = tmp2(j) + ep * delta(j)
         write(6,*) j,umin(j),umax(j)
      enddo

      ! compute distance between umax and umin
      call sub3(udis,umax,umin,nb)
      if (nid.eq.0) then
         do i=1,nb
            write(6,*)i,udis(i)
         enddo
      endif

      call dump_serial(umin,nb,'ops/umin ',nid)
      call dump_serial(umax,nb,'ops/umax ',nid)

      return
      end
c-----------------------------------------------------------------------
      subroutine hyperpar

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ep
      real wk(nb)

      call nekgsync
      hpar_time=dnekclock()

      ! eps is the free parameter
      ! 1e-2 is used in the paper
      ep = 1.e-2

      n  = lx1*ly1*lz1*nelt
      if (ifpod(1)) then
         if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
            call read_serial(umin,nb,'ops/umin ',wk,nid)
            call read_serial(umax,nb,'ops/umax ',wk,nid)
         else
            call cfill(umin,1.e9,nb)
            call cfill(umax,-1.e9,nb)
            do j=1,ns
            do i=1,nb
               if (uk(i,j).lt.umin(i)) umin(i)=uk(i,j)
               if (uk(i,j).gt.umax(i)) umax(i)=uk(i,j)
            enddo
            enddo
            do j=1,nb                    ! compute hyper-parameter
               d= umax(j)-umin(j)
               umin(j) = umin(j) - ep * d
               umax(j) = umax(j) + ep * d
               if (nio.eq.0) write (6,*) j,umin(j),umax(j)
            enddo
         endif

         ! compute distance between umax and umin
         call sub3(udis,umax,umin,nb)

         if (nio.eq.0) then
            do i=1,nb
               write (6,*) i,udis(i)
            enddo
         endif
         if (rmode.eq.'ALL'.or.rmode.eq.'OFF') then
            call dump_serial(umin,nb,'ops/umin ',nid)
            call dump_serial(umax,nb,'ops/umax ',nid)
         endif
      endif   

      if (ifpod(2)) then
         if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
            call read_serial(tmin,nb,'ops/tmin ',wk,nid)
            call read_serial(tmax,nb,'ops/tmax ',wk,nid)
         else
            call cfill(tmin,1.e9,nb)
            call cfill(tmax,-1.e9,nb)
            do j=1,ns
            do i=1,nb
               if (tk(i,j).lt.tmin(i)) tmin(i)=tk(i,j)
               if (tk(i,j).gt.tmax(i)) tmax(i)=tk(i,j)
            enddo
            enddo
            do j=1,nb                    ! compute hyper-parameter
               d= tmax(j)-tmin(j)
               tmin(j) = tmin(j) - ep * d
               tmax(j) = tmax(j) + ep * d
               if (nio.eq.0) write (6,*) j,tmin(j),tmax(j)
            enddo
         endif

         ! compute distance between tmax and tmin
         call sub3(tdis,tmax,tmin,nb)
         if (nio.eq.0) then
            do i=1,nb
               write (6,*) i,tdis(i)
            enddo
         endif

         if (rmode.eq.'ALL'.or.rmode.eq.'OFF') then
            call dump_serial(tmin,nb,'ops/tmin ',nid)
            call dump_serial(tmax,nb,'ops/tmax ',nid)
         endif
      endif   

      call nekgsync
      if (nio.eq.0) write (6,*) 'hpar_time',dnekclock()-hpar_time

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
      subroutine binv1_nom(out1)
C--------------------------------------------------------------------
C
C     Compute OUT = (B)-1 * INP   (explicit)
C
C--------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'

      parameter (lt=lx1*ly1*lz1*lelt)

      real out1(lt)

      include 'OPCTR'

      n=lx1*ly1*lz1*nelv

      call dssum(out1,lx1,ly1,lz1)
      call col2(out1,binvm1,n)

      return
      end
C--------------------------------------------------------------------
      subroutine binv1(out1)
C--------------------------------------------------------------------
C
C     Compute OUT = (B)-1 * INP   (explicit)
C
C--------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'

      parameter (lt=lx1*ly1*lz1*lelt)

      real out1(lt)

      include 'OPCTR'

      n=lx1*ly1*lz1*nelv
      
      if (ifield.eq.1) then
         call col2(out1,v1mask,n) ! disable mask for now
      else
         call col2(out1,tmask,n) ! disable mask for now
      endif
      call dssum(out1,lx1,ly1,lz1)
      call col2(out1,binvm1,n)

      return
      end
C--------------------------------------------------------------------
      subroutine invmat(a,b,c,iwk,n)

      real a(n,n),b(n,n),c(n,n)
      integer iwk1(n)

      if (n.gt.0) then
         call rzero(a,n*n)

         do i=1,n
            a(i,i)=1.
         enddo

         call copy(b,c,n*n)

         call dgetrf(n,n,b,n,iwk,info)
         call dgetrs('N',n,n,b,n,iwk,a,n,info)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cpsi(uu,vv,ww)

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)

      logical iftemp
      real uu(lt),vv(lt),ww(lt)

      common /cpsiv/ psi(lt),omega(lt,3),rhs(lt),
     $               w1(lt),w2(lt),h1(lt),h2(lt)

      common /cptmp/ t1(lt),t2(lt),t3(lt),psimask(lx1,ly1,lz1,lelt)

      integer icalld
      save    icalld
      data    icalld /0/

      n = lx1*ly1*lz1*nelv

      call opcopy(t1,t2,t3,xm1,ym1,zm1)

      ! assume homogeneous Dirichlet for stream function for now

      call rone(psimask,n)

      do ie=1,nelt
      do if=1,2*ldim
         if (cbc(if,ie,1).ne.'E  ') then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,if)
            do k=kz1,kz2
            do j=ky1,ky2
            do i=kx1,kx2
               psimask(i,j,k,ie)=0.
            enddo
            enddo
            enddo
         endif
      enddo
      enddo

      iftemp=ifaxis
      ifaxis=.false.
      call comp_vort3(omega,w1,w2,uu,vv,ww)
      call col3(rhs,bm1,omega,n)
      call rone(h1,n)
      call rzero(h2,n)
      tol = param(22)
c     call chsign(rhs,n)
      call hmholtz('psi ',psi,rhs,h1,h2,psimask,vmult,1,tol,1000,1)
      call dsavg(psi)
      ifaxis=iftemp

      call outpost(psi,omega,vz,pr,t,'psi')
      ifxyo=iftemp

      call opcopy(xm1,ym1,zm1,t1,t2,t3)

      return
      end
c-----------------------------------------------------------------------
      subroutine csga(ga,sb)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ga(lt,0:nb),sb(lt,0:nb)

      n=lx1*ly1*lz1*nelv

      do ib=0,nb
         call axhelm(ga(1,ib),sb(1,ib),ones,zeros,1,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine emask(g,ie)

      include 'SIZE'
      include 'TOTAL'

      real g(lx1*ly1*lz1,lelt)

      n=lx1*ly1*lz1*nelt
      
      do i=1,nelt
         if (ie.ne.i) call rzero(g(1,i),lx1*ly1*lz1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_surf(s,f) ! set surface term using flux

      include 'SIZE'
      include 'TOTAL'

      real s(lx1,ly1,lz1,lelt)
      real f(lx1,ly1,lz1,lelt)

      n=lx1*ly1*lz1*nelv

      call rzero(s,n)

      do ie=1,nelv
      do ifc=1,2*ldim
         ia=0
         call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,ifc)
         do iz=kz1,kz2
         do iy=ky1,ky2
         do ix=kx1,kx2
            ia=ia+1
            s(ix,iy,iz,ie)=s(ix,iy,iz,ie)
     $                    +f(ix,iy,iz,ie)*area(ia,1,ifc,ie)
         enddo
         enddo
         enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_gradn(gn,s) ! set grad s . n

      include 'SIZE'
      include 'TOTAL'

      common /scrgn/ sx(lx1,ly1,lz1,lelt),
     $               sy(lx1,ly1,lz1,lelt),
     $               sz(lx1,ly1,lz1,lelt)

      real gn(lx1,ly1,lz1,lelt)

      n=lx1*ly1*lz1*nelv

      call rzero(gn,n)
      call gradm1(sx,sy,sz,s)

      do ie=1,nelv
      do ifc=1,2*ldim
         ia=0
         if (cbc(ifc,ie,2).ne.'E  ') then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,ifc)
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
               ia=ia+1
               gn(ix,iy,iz,ie)=sx(ix,iy,iz,ie)*unx(ia,1,ifc,ie)
     $                        +sy(ix,iy,iz,ie)*uny(ia,1,ifc,ie)
     $                        +sz(ix,iy,iz,ie)*unz(ia,1,ifc,ie)
            enddo
            enddo
            enddo
         endif
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_sfld

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ux(lt),uy(lt),uz(lt)

      logical iftmp

      common /scrdump2/ ux1(lt),uy1(lt),uz1(lt),tt(lt),wk(lt)
      common /testb/ ux2(lt),uy2(lt),uz2(lt)

      iftmp=ifxyo
      ifxyo=.true.

      if (ifrom(1)) then
         call reconv(ux1,uy1,uz1,ua)
         if (ifrom(2)) call recont(tt,uta)
         call outpost(ux1,uy1,uz1,pr,tt,'avg')

         call reconu_rms(ux2,uy2,uz2,u2a)
         if (ifrom(2)) call recont_rms(tt)
         call outpost(ux2,uy2,uz2,pr,tt,'rms')

         call opcol2(ux1,uy1,uz1,ux1,uy1,uz1)
         call opsub2(ux2,uy2,uz2,ux1,uy1,uz1)
         do i=1,lx1*ly1*lz1*nelv
            wk(i)=.5*(ux2(i)+uy2(i)+uz2(i))
         enddo

         call outpost(ux2,uy2,uz2,pr,wk,'tke')
      endif

      if (ifrom(1).and.ifrom(2)) then
         n=lx1*ly1*lz1*nelt
         call opzero(ux1,uy1,uz1)
         do j=0,nb
         do i=0,nb
            call admcol3(ux1,ub(1,i),tb(1,j),uuta(1+i+(nb+1)*j),n)
            call admcol3(uy1,vb(1,i),tb(1,j),uuta(1+i+(nb+1)*j),n)
            if (ldim.eq.3)
     $         call admcol3(uz1,wb(1,i),tb(1,j),uuta(1+i+(nb+1)*j),n)
         enddo
         enddo
         call outpost(ux1,uy1,uz1,pr,tt,'tmn')
      endif

      ifxyo=iftmp

      return
      end
c-----------------------------------------------------------------------
      subroutine sol_intp_xline(ux,uy,uz,ut,ytgt,ztgt,nx,pfx)

      include 'SIZE'
      include 'TOTAL'

      parameter (lmax=1024)
      parameter (lt=lx1*ly1*lz1*lelt)

      character*3 pfx
      character*127 fname

      real ux(lt),uy(lt),uz(lt),ut(lt)

      real xi(lmax),yi(lmax),zi(lmax)
      real uxi(lmax),uyi(lmax),uzi(lmax),uti(lmax)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside sol_intp_xline'

      nn=0
      if (nid.eq.0) nn=nx

      s=pi/(2.*nx)
      do i=1,nn
         xi(i)=cos((2*i-1.)*s)*.5
         yi(i)=ytgt
         zi(i)=ztgt
      enddo

      if (nio.eq.0) write (6,*) 'wp1'

      call interp_vec_prop2(uxi,uyi,uzi,ux,uy,uz,xi,yi,zi,nn)

      if (nio.eq.0) write (6,*) 'wp2'

      call interp_sca_prop2(uti,ut,xi,yi,zi,nn)

      if (nio.eq.0) write (6,*) 'wp3'

      if (nid.eq.0) then
c        blank(fname,127)
c        write (fname,'(A3)') i
c        open (unit=10,file='')
         do i=1,nx
            write (6,1) i,xi(i),uxi(i),uyi(i),uzi(i),uti(i),pfx
            write (10,1) i,xi(i),uxi(i),uyi(i),uzi(i),uti(i),pfx
         enddo
      endif

    1 format (i5,1p5e16.8,1x,'intp_result',a3)

      call nekgsync

      return
      end
c-----------------------------------------------------------------------
      subroutine sol_intp_xline_qoi(ux,uy,uz,ut,ytgt,ztgt,nx,j)

      include 'SIZE'
      include 'TOTAL'

      parameter (lmax=1024)
      parameter (lt=lx1*ly1*lz1*lelt)

      character*3 pfx
      character*127 fname

      real ux(lt),uy(lt),uz(lt),ut(lt)

      real xi(lmax),yi(lmax),zi(lmax)
      real uxi(lmax),uyi(lmax),uzi(lmax),uti(lmax)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside sol_intp_xline',j

      nn=0
      if (nid.eq.0) nn=nx

      pfx='ttt'

      s=pi/(2.*nx)
      do i=1,nn
         xi(i)=cos((2*i-1.)*s)*.5
         yi(i)=ytgt
         zi(i)=ztgt
      enddo

      call interp_vec_prop2(uxi,uyi,uzi,ux,uy,uz,xi,yi,zi,nn)
      call interp_sca_prop2(uti,ut,xi,yi,zi,nn)

      if (nid.eq.0) then
         call blank(fname,127)
         write (fname,'(A1,I0,A4)') 'q',j,'.dat'
         open (unit=10,file=fname)
         do i=1,nx
            write (10,1) i,xi(i),uxi(i),uyi(i),uzi(i),uti(i),pfx
         enddo
         close (unit=10)
      endif

    1 format (i5,1p5e16.8,1x,'intp_result',a3)

      call nekgsync

      return
      end
c-----------------------------------------------------------------------
      subroutine sol_intp_xline_fqoi(ux,uy,uz,ut,ytgt,ztgt,nx)

      include 'SIZE'
      include 'TOTAL'

      parameter (lmax=1024)
      parameter (lt=lx1*ly1*lz1*lelt)

      character*3 pfx
      character*127 fname

      real ux(lt),uy(lt),uz(lt),ut(lt)

      real xi(lmax),yi(lmax),zi(lmax)
      real uxi(lmax),uyi(lmax),uzi(lmax),uti(lmax)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside sol_intp_xline'

      nn=0
      if (nid.eq.0) nn=nx

      pfx='fom'

      s=pi/(2.*nx)
      do i=1,nn
         xi(i)=cos((2*i-1.)*s)*.5
         yi(i)=ytgt
         zi(i)=ztgt
      enddo

      call interp_vec_prop2(uxi,uyi,uzi,ux,uy,uz,xi,yi,zi,nn)
      call interp_sca_prop2(uti,ut,xi,yi,zi,nn)

      if (nid.eq.0) then
         call blank(fname,127)
         fname='fq.dat'
         open (unit=10,file=fname)
         do i=1,nx
            write (10,1) i,xi(i),uxi(i),uyi(i),uzi(i),uti(i),pfx
         enddo
         close (unit=10)
      endif

    1 format (i5,1p5e16.8,1x,'intp_result',a3)

      call nekgsync

      return
      end
c-----------------------------------------------------------------------
      subroutine interp_vec_prop2(uu,vv,ww,ux,uy,uz,xx,yy,zz,n)
c
c     evaluate velocity for list of points xyz using separate field arrays
c
      include 'SIZE'
      include 'TOTAL'

      real uvw(ldim,n),xyz(ldim,n)
      real uu(n),vv(n),ww(n),xx(n),yy(n),zz(n)
      real ux(n),uy(n),uz(n)

      do i=1,n
         xyz(1,i)=xx(i)
         xyz(2,i)=yy(i)
         if (ldim.eq.3) xyz(3,i)=zz(i)
      enddo

      call interp_vec2(uvw,ux,uy,uz,xyz,n)

      do i=1,n
         uu(i)=uvw(1,i)
         vv(i)=uvw(2,i)
         if (ldim.eq.3) ww(i)=uvw(3,i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine interp_vec2(uvw,ux,uy,uz,xyz,n)
c
c     evaluate velocity for list of points xyz
c
      include 'SIZE'
      include 'TOTAL'

      parameter (intp_nmax=1024)
      real uvw(ldim,n),ux(n),uy(n),uz(n),xyz(ldim,n)

      real    rwk(INTP_NMAX,ldim+1) ! r, s, t, dist2
      integer iwk(INTP_NMAX,3)      ! code, proc, el
      save    rwk, iwk

      integer intp_h
      save    intp_h

      common /rwk_intp/
     $       fwrk(lx1*ly1*lz1*lelt,ldim),
     $       fpts(ldim*INTP_NMAX),
     $       pts(ldim*INTP_NMAX)

      integer icalld,e
      save    icalld
      data    icalld /0/

      nxyz  = nx1*ny1*nz1
      ntot  = nxyz*nelv

      if (n.gt.INTP_NMAX) call exitti ('n > INTP_NMAX in interp_v!$',n)

      do i=1,n ! ? not moving -> save?
         pts(i)     = xyz(1,i)
         pts(i + n) = xyz(2,i)
         if (ldim.eq.3) pts(i + n*2) = xyz(3,i)
      enddo

      if (icalld.eq.0) then
        icalld = 1
        call interp_setup(intp_h,0,0,nelt) ! prod
c       call interp_setup(0,0,intp_h)      ! v17
        if (nio.eq.0) write (6,*) 'finished interp_setup'
      endif

      ! pack working array
      call opcopy(fwrk(1,1),fwrk(1,2),fwrk(1,3),ux,uy,uz)

      ndim=ldim

      ! interpolate
      call interp_nfld(fpts,fwrk,ndim,pts(1),pts(1+n),pts(2*n+1),
     $                 n,iwk,rwk,INTP_NMAX,.true.,intp_h)

      do i=1,n
         uvw(1,i) = fpts(i)
         uvw(2,i) = fpts(i + n)
         if (ldim.eq.3) uvw(3,i) = fpts(i + n*2)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine interp_sca_prop2(u,s,xx,yy,zz,n)
c
c     evaluate velocity for list of points xyz using separate field arrays
c
      include 'SIZE'
      include 'TOTAL'

      real xyz(ldim,n)
      real u(n),xx(n),yy(n),zz(n)
      real s(1)

      do i=1,n
         xyz(1,i)=xx(i)
         xyz(2,i)=yy(i)
         if (ldim.eq.3) xyz(3,i)=zz(i)
      enddo

      call interp_sca2(u,xyz,s,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine interp_sca2(u,xyz,s,n)
c
c     evaluate scalar field, s, for list of points xyz
c
      include 'SIZE'
      include 'TOTAL'

      parameter (intp_nmax=1024)
      real u(n),s(n),xyz(ldim,n)

      real    rwk(INTP_NMAX,ldim+1) ! r, s, t, dist2
      integer iwk(INTP_NMAX,3)      ! code, proc, el
      save    rwk, iwk

      integer intp_h
      save    intp_h

      common /rwk_intp/
     $       fwrk(lx1*ly1*lz1*lelt),
     $       fpts(INTP_NMAX),
     $       pts(INTP_NMAX)

      integer icalld,e
      save    icalld
      data    icalld /0/

      nxyz  = nx1*ny1*nz1
      ntot  = nxyz*nelv

      if (n.gt.INTP_NMAX) call exitti ('n > INTP_NMAX in interp_s!$',n)

      do i=1,n ! ? not moving -> save?
         pts(i)   = xyz(1,i)
         pts(i+n) = xyz(2,i)
         if (ldim.eq.3) pts(i+2*n) = xyz(3,i)
      enddo

      if (icalld.eq.0) then
        icalld = 1
c       call interp_setup(0,0,intp_h)      ! v17
        call interp_setup(intp_h,0,0,nelt) ! prod
        if (nio.eq.0) write (6,*) 'done with setup'
      endif

      ! pack working array
      call copy(fwrk,s,lx1*ly1*lz1*nelv)

      ! interpolate
      mdim=1
      call interp_nfld(fpts,fwrk,mdim,pts(1),pts(1+n),pts(2*n+1),
     $                 n,iwk,rwk,INTP_NMAX,.true.,intp_h)

      call copy(u,fpts,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine count_gal(cgal,acgal,uu,amax,amin,bctol,n)

      real cgal(n),uu(n)
      real amax(n),amin(n)
      real bctol,acgal
      integer chekbc

      chekbc=0

      do ii=1,n
         if ((uu(ii)-amax(ii)).ge.bctol.OR.(amin(ii)-uu(ii)).ge.bctol) 
     $   then
            chekbc=1
         else
            cgal(ii) = cgal(ii) + 1
         endif
      enddo

      if (chekbc.eq.0) acgal=acgal+1

      return
      end
c-----------------------------------------------------------------------
      subroutine cdump

      include 'SIZE'
      include 'MOR'

      complex ceval,cevec,work

      common /ceigscr/ ceval(lb),cevec(lb*lb),work(lb*lb)

      call gop(ctmp(nb+1),work,'+  ',nb*nb)

      do i=1,nb*nb
         if (nio.eq.0) write (6,*) ad_step,ctmp(nb+i),'ccc'
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_binv2

      include 'SIZE'
      include 'MOR'

      common /invevf/ buinv(lb*lb),btinv(lb*lb)

      call invmat(buinv,rtmp1,bu,itmp1,nb)
      call invmat(btinv,rtmp1,bt,itmp1,nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine trace

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrtr/ ux(lt),uy(lt),uz(lt)

      character*127 fmat

      if (istep.eq.0) then
         call rom_init_params
         call rom_init_fields
         call loadbases
      else
         if (istep.gt.lcs) then
            if (nio.eq.0) write (6,*) 'WARNING: lcs <= nsteps'
         else
            if (ifheat) then
               n=lx1*ly1*lz1*nelt
               call sub3(ux,t,tb,n)
               call ps2b(tk(0,istep),ux,uy,uz,ub,vb,wb)
            endif
            if (ifflow) then
               call opsub3(ux,uy,uz,vx,vy,vz,ub,vb,wb)
               call pv2b(uk(0,istep),ux,uy,uz,ub,vb,wb)
            endif
         endif
      endif

      if (istep.eq.nsteps) then
      if (nio.eq.0) then
         call blank(fmat,127)
         write (fmat,2) nb+1
         if (ifflow) then
            open (unit=10,file='ops/utrace')
            do i=1,min(lcs,nsteps)
               write (10,fmat) (uk(j,i),j=0,nb)
            enddo
            close (unit=10)

            call reconv(ux,uy,uz,uk(0,istep))
            call outpost(vx,vy,vz,pr,t,'err')
            call outpost(ux,uy,uz,pr,t,'err')
            call opsub2(ux,uy,uz,vx,vy,vz)
            call outpost(ux,uy,uz,pr,t,'err')

            err=sqrt(op_glsc2_wt(ux,uy,uz,ux,uy,uz,bm1))
            ul2=sqrt(op_glsc2_wt(vx,vy,vz,vx,vy,vz,bm1))
            if (nio.eq.0) write (6,*) err,ul2,err/ul2,'err'
         endif
         if (ifheat) then
            open (unit=10,file='ops/ttrace')
            do i=1,min(lcs,nsteps)
               write (10,fmat) (tk(j,i),j=0,nb)
            enddo
            close (unit=10)
         endif
      endif
      endif

    2 format ('(1p',i4,'e25.17)')

      return
      end
c-----------------------------------------------------------------------
      subroutine set_trace

      include 'SIZE'
      include 'MOR'

      if (nio.eq.0) write (6,*) 'inside set_trace'

      open (unit=10,file='ops/utrace')

      if (nio.eq.0) write (6,*) ad_nsteps,lcs,nb
      do i=1,min(ad_nsteps,lcs)
         read (10,*) (uk(j,i),j=0,nb)
      enddo

      close (unit=10)

      return
      end
c-----------------------------------------------------------------------
      subroutine diag(h,wt,wk,n)

      real h(n,n),wt(n,n),wk(n)

      call regularev(h,wk,n,wt)

      do j=1,n
      do i=1,n
         wt(i,j)=h(j,i)
      enddo
      enddo

      call rzero(h,n*n)

      do i=1,n
         h(i,i)=1./wk(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine checker(cstr,myind)

      include 'SIZE'
      include 'MOR'

      character*3 cstr

      vmin=vlmin(u,nb+1)
      vmax=vlmax(u,nb+1)
      vmean=vlsum(u,nb+1)/(nb+1.)
      vprod=1.
      do i=0,nb
         vprod=vprod*u(i)
      enddo

      if (nio.eq.0) write (6,1) vmin,vmax,vmean,vprod,cstr,myind
    1 format('check ',1p4e13.5,1x,a,1x,i8)

      return
      end
c-----------------------------------------------------------------------
      subroutine checkera(cstr,a,nn,myind)

      include 'SIZE'
      include 'MOR'

      character*3 cstr

      real a(nn)

      vmin=vlmin(a,nn)
      vmax=vlmax(a,nn)
      vmean=vlsum(a,nn)/(nn*1.)
      vprod=1.
      do i=1,nn
         vprod=vprod*a(i)
      enddo

      if (nio.eq.0) write (6,1) vmin,vmax,vmean,vprod,cstr,myind
    1 format('check ',1p4e13.5,1x,a,1x,i8)

      return
      end
c-----------------------------------------------------------------------
