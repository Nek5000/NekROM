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
         call add2s2(tt,tb(1,i,1),coef(i),n)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine recont_rms(tt)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk1/ tbt(lt)

      real tt(lt)

      n=lx1*ly1*lz1*nelt

      call rzero(tt,n)

      do j=0,nb
      do i=0,nb
         call col3(tbt,tb(1,i,1),tb(1,j,1),n)
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
      subroutine reconu_rm2(ux,uy,uz)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ux(lt),uy(lt),uz(lt)

      n=lx1*ly1*lz1*nelv

      call opzero(ux,uy,uz)

      do j=0,nb
         if (nio.eq.0) write (6,*) 'reconu_rms:',j,'/',nb
         do i=0,nb
            call admcol3(ux,ub(1,i),vb(1,j),u2a(1+i+(nb+1)*j),n)
            call admcol3(uy,vb(1,i),wb(1,j),u2a(1+i+(nb+1)*j),n)
            if (ldim.eq.3)
     $         call admcol3(uz,wb(1,i),ub(1,j),u2a(1+i+(nb+1)*j),n)
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
      subroutine shift(u,v,n,m)

      ! shift v array into u array

      ! u := target array
      ! v := input array
      ! n := length of v
      ! m := number of shifts

      real u(n,m),v(n)

      do i=1,m-1
         call copy(u(1,m+1-i),u(1,m-i),n)
      enddo

      call copy(u,v,n)

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
      subroutine lap3d(d2u,u)

      ! set Laplacian of the input scalar field

      ! d2u := Laplacian field
      ! u   := input scalar field

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)

      real d2u(1),u(1)

      common /scrl2d/ ux(lt),uy(lt),uz(lt),
     $                uxx(lt),uyy(lt),uzz(lt),t1(lt),t2(lt)

      n=lx1*ly1*lz1*nelt

      call gradm1(ux,uy,uz,u)

      call gradm1(uxx,t1,t2,ux)
      call gradm1(t1,uyy,t2,uy)
      call gradm1(t1,t2,uzz,uz)

      call add3(d2u,uxx,uyy,n)
      call add2(d2u,uzz,n)

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

      do idim=1,ldimt
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
      include 'INPUT'

      parameter (lt1=lx1*ly1*lz1*lelt)
      parameter (lt2=lx2*ly2*lz2*lelt)

      real vx(lt1),vy(lt1),vz(lt1),pr(lt2),t(lt1,ldimt)
      real ux(lt1),uy(lt1),uz(lt1),pp(lt2),tt(lt1,ldimt)

      call opcopy(vx,vy,vz,ux,uy,uz)
      call copy(pr,pp,lx2*ly2*lz2*nelv)

      do idim=1,min(1+npscal,ldimt)
         call copy(t(1,idim),tt(1,idim),lx1*ly1*lz1*nelt)
      enddo

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
      subroutine hyperpar(uuk,ttk)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ep
      real wk(nb)
      real uuk(0:nb,ns),ttk(0:nb,ns)

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
               if (uuk(i,j).lt.umin(i)) umin(i)=uuk(i,j)
               if (uuk(i,j).gt.umax(i)) umax(i)=uuk(i,j)
            enddo
            enddo
            do j=1,nb                    ! compute hyper-parameter
               d= umax(j)-umin(j)
               umin(j) = umin(j) - ep * d
               umax(j) = umax(j) + ep * d
               if (nio.eq.0) write (6,*) j,umin(j),umax(j),'uminmax'
            enddo
         endif

         ! compute distance between umax and umin
         call sub3(udis,umax,umin,nb)

         if (nio.eq.0) then
            do i=1,nb
               write (6,*) i,udis(i),'udis'
            enddo
         endif
         if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ') then
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
               if (ttk(i,j).lt.tmin(i)) tmin(i)=ttk(i,j)
               if (ttk(i,j).gt.tmax(i)) tmax(i)=ttk(i,j)
            enddo
            enddo
            do j=1,nb                    ! compute hyper-parameter
               d= tmax(j)-tmin(j)
               tmin(j) = tmin(j) - ep * d
               tmax(j) = tmax(j) + ep * d
               if (nio.eq.0) write (6,*) j,tmin(j),tmax(j),'tminmax'
            enddo
         endif

         ! compute distance between tmax and tmin
         call sub3(tdis,tmax,tmin,nb)
         if (nio.eq.0) then
            do i=1,nb
               write (6,*) i,tdis(i),'tdis'
            enddo
         endif

         if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ') then
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
      subroutine invmat(a,b,c,iwk1,iwk2,n)

      real a(n,n),b(n,n),c(n,n)
      integer iwk1(n),iwk2(n)

      if (n.gt.0) then
         call rzero(a,n*n)

         do i=1,n
            a(i,i)=1.
         enddo

         call copy(b,c,n*n)

         call lu(b,n,n,iwk1,iwk2)
c        call dgetrf(n,n,b,n,iwk1,info)

         call solve(a,b,n,n,n,iwk1,iwk2)
c        call dgetrs('N',n,n,b,n,iwk1,a,n,info)
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
      common /testb/ ux2(lt),uy2(lt),uz2(lt),ux3(lt),uy3(lt),uz3(lt)

      iftmp=ifxyo
      ifxyo=.true.

      if (ifrom(1)) then
         call reconv(ux1,uy1,uz1,ua)
         if (ifrom(2)) call recont(tt,uta)
         call outpost(ux1,uy1,uz1,pr,tt,'avg')

         call reconu_rms(ux2,uy2,uz2,u2a)
         if (ifrom(2)) call recont_rms(tt)
         call outpost(ux2,uy2,uz2,pr,tt,'rms')

         call reconu_rm2(ux3,uy3,uz3,u2a)
         call outpost(ux3,uy3,uz3,pr,tt,'rm2')

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
            call admcol3(ux1,ub(1,i),tb(1,j,1),uuta(1+i+(nb+1)*j),n)
            call admcol3(uy1,vb(1,i),tb(1,j,1),uuta(1+i+(nb+1)*j),n)
            if (ldim.eq.3)
     $         call admcol3(uz1,wb(1,i),tb(1,j,1),uuta(1+i+(nb+1)*j),n)
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

      call invmat(buinv,rtmp1,bu,itmp1,itmp2,nb)
      call invmat(btinv,rtmp1,bt,itmp1,itmp2,nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine diag(h,wt,wk,n)

      real h(n,n),wt(n,n),wk(n)

      call regularev(h,wk,n,wt,max(n*n,2))

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

      return

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

      return

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
      subroutine k_mean(k,nsu,nsp,nst,fn,seed)

      ! K-means Clustering

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      integer seed
      integer k           ! number of clusters
      integer label(k) 
      integer itermax

      real centroid(k)    ! centroid parameter
      real cent_fld(lt,k) ! centroid field
      real sample(ls)
      real dist(k)
      real num_sc(k)
      real rnk(ls,k)      ! binary indicator
      real tmp(lt,k)      ! dummy variable
      real tmpp(ls)       ! dummy variable
      real obj_f          ! distortion measure 

      character*128 fn
      character*128 fnlint

      n=lx1*ly1*lz1*nelt
      itermax = 50

      ! initialize centroid
      ! currently put here for same initialization
      call srand(seed)
      do i=1,k
         centroid(i) = rand()
      enddo

      ! scale up to parameter range
      pmax = 180
      pmin = 0
      pdiff = pmax-pmin
      do i=1,k
         centroid(i) = pmin + pdiff*centroid(i)
         write(6,*)i,centroid(i),nint(centroid(i))
      enddo

      ! read sample.list
      ierr = 0
      call lints(fnlint,fn,128)
      if (nid.eq.0) open(77,file=fnlint,status='old',err=199)
      ierr = iglmax(ierr,1)
      if (ierr.gt.0) goto 199

      nsave=max(max(nsu,nsp),nst)
c     write(6,*)'nsave',nsave
      icount = 0
c     do ipass=1,nsave
      do ipass=1,ls
         call blank(initc,127)
         initc(1) = 'done '
         if (nid.eq.0) read(77,127,end=998)initc(1) 
         read(initc,'(f10.0)') sample(ipass)
c        write(6,*)sample(ipass)
  998    call bcast(initc,127)
  127    format(a127)
      enddo
      close(77)

      do i=1,k
         do j=1,ls
            if (abs(centroid(i)-sample(j))<5) then  
               centroid(i) = sample(j)
                  write(6,*)i,centroid(i),sample(j)
               label(i) = j
            endif
         enddo
         call copy(cent_fld(1,i),ts0(1,label(i),1),n)
      enddo

      ! minimize distortion measure
      do kk=1,itermax

         ! assign each samlpe to cluster
         call rzero(rnk,ls*k)
         do i=1,ls
            do j=1,k
              call sub3(tmp(1,j),ts0(1,i,1),cent_fld(1,j),n)
              dist(j) = glsc2(tmp(1,j),tmp(1,j),n)
            enddo
            write(6,*)ls,minloc(dist),sample(i)
            do j=1,k
               write(6,*)dist(j)
               if (minloc(dist,1).eq.j) rnk(i,j) = 1
            enddo
         enddo

         call c_distortion_measure(obj_f,cent_fld,rnk,k)
         write(6,*)kk,obj_f,'distortion measure E'

         call rone(tmpp,ls)
         do i=1,k
            num_sc(i) = glsc2(rnk(1,i),tmpp,ls)
c           write(6,*)1./num_sc(i),num_sc(i),'num_sc'
         enddo

         ! compute new centroid
c        do i=1,k
c           write(6,*)glmax(cent_fld(1,i),n),'old max'
c           call mxm(ts0,n,rnk(1,i),ls,cent_fld(1,i),1)
c           call cmult(cent_fld(1,i),1./num_sc(i),n)
c           centroid(i) = glsc2(sample,rnk(1,i),ls)/num_sc(i)
c           write(6,*)i,centroid(i),'new centroid'
c           write(6,*)glmax(cent_fld(1,i),n),'new max'
c        enddo
         do i=1,k
            call rzero(cent_fld,n*k)
            do j=1,ls
               call add2s2(cent_fld(1,i),ts0(1,j,1),rnk(j,i),n)
            enddo
            call cmult(cent_fld(1,i),1./num_sc(i),n)
            centroid(i) = glsc2(sample,rnk(1,i),ls)/num_sc(i)
c           write(6,*)i,centroid(i),'new centroid'
         enddo

         ! assign the closest sample to be centroid
         ! currently does not have enough sample
         do i=1,k
            do j=1,ls
               if (abs(centroid(i)-sample(j))<5) then  
                  centroid(i) = sample(j)
c                 write(6,*)i,centroid(i),sample(j)
                  label(i) = j
               endif
            enddo
            call copy(cent_fld(1,i),ts0(1,label(i),1),n)
         enddo
         call c_distortion_measure(obj_f,cent_fld,rnk,k)
         write(6,*)kk,obj_f,'distortion measure M'
      enddo

      ! write out clusters
      do j=1,k
         write(6,*) 'cluster: ',j,' centroid: ', centroid(j)
         do i=1,ls
            if (abs(rnk(i,j)-1).le.1e-8) then
               write(6,*) sample(i) 
            endif
         enddo
      enddo

      write(6,*) obj_f,'distortion for cluster: ',k

       
      return

  199 continue ! exception handle for file not found
      ierr = 1
      if (nid.eq.0) ierr = iglmax(ierr,1)
      write (6,*) fnlint
      call exitti('get_saved_fields did not find list file.$',ierr)


      return
      end
c-----------------------------------------------------------------------
      subroutine c_distortion_measure(obj_f,cent_fld,rnk,k)

      ! Compute objective function (distortion measure)
      ! J = \sum_{n=1}^N sum^K_{k=1} r_{nk} \|x_n - \nu_k\|^2
      ! Sum of the squares of the distances of each data point to its
      ! assigned vector \nu_k

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      integer k           ! number of clusters
      real rnk(ls,k)      ! binary indicator
      real cent_fld(lt,k) ! centroid
      real obj_f          ! objective function value
      real tmp(lt,k)
      real dist(k)

      n=lx1*ly1*lz1*nelt

      ! assign each samlpe to cluster
      obj_f=0
      do i=1,ls
         do j=1,k
            if (abs(rnk(i,j)-1).le.1e-8) then
               call sub3(tmp(1,j),ts0(1,i,1),cent_fld(1,j),n)
               dist(j) = glsc2(tmp(1,j),tmp(1,j),n)
               obj_f = obj_f + dist(j)
            endif
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine projtoprerb(nocp,uuk,ttk)

      ! This subroutine is for p-greedy. It project the RB basis
      ! read in by loadbases onto space that is perpendicular to the
      ! space spanned by bases in pbas.list

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      character*128 fn
      character*128 fnlint

      integer nocp

      real uuk(0:nb,ns),ttk(0:nb,ns)

      n=lx1*ly1*lz1*nelt

      if (nio.eq.0) write (6,*) 'inside projtoprerb'

      ! project onto the previous RB space
      do j=1,ns
      uuk(0,j) = 1.
      do i=1,nocp
         ww=vip(ub(1,i),vb(1,i),wb(1,i),ub(1,i),vb(1,i),wb(1,i))
         vv=vip(ub(1,i),vb(1,i),wb(1,i),
     $          us0(1,1,j),us0(1,2,j),us0(1,3,j))
         uuk(i,j) = vv/ww
      enddo
      enddo

      do j=1,ns
      ttk(0,j) = 1.
      do i=1,nocp
         ww=sip(tb(1,i,1),tb(1,i,1))
         vv=sip(tb(1,i,1),ts0(1,j,1))
         ttk(i,j) = vv/ww
      enddo
      enddo

      ! project onto the space perpendicular to
      ! the previous RB space
      do i=1,ns
         call reconv_wo0(vx,vy,vz,uuk(0,i),nocp)
         call recont_wo0(t,ttk(0,i),nocp)
         if (ifrom(1)) then
            call sub2(us0(1,1,i),vx,n)
            call sub2(us0(1,2,i),vy,n)
            if (ldim.eq.3) call sub2(us0(1,ldim,i),vz,n)
         endif
         if (ifrom(2)) call sub2(ts0(1,i,1),t,n)
      enddo

      if (nio.eq.0) write (6,*) 'exiting projtoprerb'
      
      return 
      end
c-----------------------------------------------------------------------
      subroutine evaldut(ev,ut,vt,wt,u,v,w,t)

      ! compute <div(u't')> where u := (u,v,w)

      ! ev         := result
      ! (ut,vt,wt) := components of mean temperature-velocity product
      ! (u,v,w)    := mean velocity components
      ! t          := mean temperature

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /myblock4/ t1(lt),t2(lt),t3(lt)

      real ev(lt),ut(lt),vt(lt),wt(lt),u(lt),v(lt),w(lt),t(lt)

      n=lx1*ly1*lz1*nelt

      call opcopy(t1,t2,t3,ut,vt,wt)

      call admcol3(t1,u,t,-1.,n)
      call admcol3(t2,v,t,-1.,n)
      if (ldim.eq.3) call admcol3(t3,w,t,-1.,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine evalut(utp,vtp,wtp,ut,vt,wt,u,v,w,t)

      ! compute <u't'> where u := (u,v,w)

      ! (utp,vtp,wtp) := components of mean temperature-velocity fluctuation
      ! (ut,vt,wt)    := components of mean temperature-velocity product
      ! (u,v,w)       := mean velocity components
      ! t             := mean temperature

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)

      real utp(lt),vtp(lt),wtp(lt)
      real ut(lt),vt(lt),wt(lt),u(lt),v(lt),w(lt),t(lt)

      n=lx1*ly1*lz1*nelv

      call opcopy(utp,vtp,wtp,ut,vt,wt)

      call admcol3(utp,u,t,-1.,n)
      call admcol3(vtp,v,t,-1.,n)
      if (ldim.eq.3) call admcol3(wtp,w,t,-1.,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine reconv_wo0(ux,uy,uz,coef,ncop)

      ! reconstruct velocity field without 0th mode

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ux(lt),uy(lt),uz(lt),coef(0:nb)
      integer ncop

      n=lx1*ly1*lz1*nelv

      call opzero(ux,uy,uz)

      do i=1,ncop
         call add2s2(ux,ub(1,i),coef(i),n)
         call add2s2(uy,vb(1,i),coef(i),n)
         call add2s2(uz,wb(1,i),coef(i),n)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine divm1(div,ux,uy,uz)

      ! compute divergence of <ux,uy,uz> -- mesh 1 to mesh 1 (vec. to sca.)

      include 'SIZE'
      include 'DXYZ'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'

      parameter (lxyz=lx1*ly1*lz1)
      real ux(lxyz,1),uy(lxyz,1),uz(lxyz,1),div(lxyz,1)

      common /ctmp1/ ur(lxyz),us(lxyz),ut(lxyz)

      n = lx1-1
      do ie=1,nelt
         if (ldim.eq.3) then
            call local_grad3(ur,us,ut,ux,n,ie,dxm1,dxtm1)
            do i=1,lxyz
               div(i,ie) =         jacmi(i,ie)*(ur(i)*rxm1(i,1,1,ie)
     $                                        + us(i)*sxm1(i,1,1,ie)
     $                                        + ut(i)*txm1(i,1,1,ie))
            enddo
            call local_grad3(ur,us,ut,uy,n,ie,dxm1,dxtm1)
            do i=1,lxyz
               div(i,ie)=div(i,ie)+jacmi(i,ie)*(ur(i)*rym1(i,1,1,ie)
     $                                        + us(i)*sym1(i,1,1,ie)
     $                                        + ut(i)*tym1(i,1,1,ie))
            enddo
            call local_grad3(ur,us,ut,uz,n,ie,dxm1,dxtm1)
            do i=1,lxyz
               div(i,ie)=div(i,ie)+jacmi(i,ie)*(ur(i)*rzm1(i,1,1,ie)
     $                                        + us(i)*szm1(i,1,1,ie)
     $                                        + ut(i)*tzm1(i,1,1,ie))
            enddo
         else
            if (ifaxis) call setaxdy (ifrzer(ie))
            call local_grad2(ur,us,ux,n,ie,dxm1,dytm1)
            do i=1,lxyz
               div(i,ie) =         jacmi(i,ie)*(ur(i)*rxm1(i,1,1,ie)
     $                                        + us(i)*sxm1(i,1,1,ie))
            enddo
            call local_grad2(ur,us,uy,n,ie,dxm1,dytm1)
            do i=1,lxyz
               div(i,ie)=div(i,ie)+jacmi(i,ie)*(ur(i)*rym1(i,1,1,ie)
     $                                        + us(i)*sym1(i,1,1,ie))
            enddo
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine recont_wo0(tt,coef,nocp)

      ! reconstruct temperature field without 0th mode

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real tt(lt),coef(0:nb)
      integer nocp

      n=lx1*ly1*lz1*nelt

      call rzero(tt,n)

      do i=1,nocp
         call add2s2(tt,tb(1,i,1),coef(i),n)
      enddo

      return
      end
c-----------------------------------------------------------------------
      function svint(s1,s2,s3)

      ! compute unit normal integration of the input vector field

      ! <s1,s2,s3> := input vector field

      include 'SIZE'
      include 'TOTAL'

      real s1(lx1,ly1,lz1,lelt),
     $                 s2(lx1,ly1,lz1,lelt),
     $                 s3(lx1,ly1,lz1,lelt)

      integer e,f,eg

      lxyz  = lx1*ly1*lz1
      nface = 2*ldim

      iobj=1

      a=0.
      s=0.

      do ie=1,nelt
      do ifc=1,2*ldim
         if (cbc(ifc,ie,1).ne.'E  '.and.cbc(ifc,ie,1).ne.'P  ') then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,ifc)
            l=1

            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
               a=a+area(l,1,ifc,ie)
               s=s+area(l,1,ifc,ie)*(unx(l,1,ifc,ie)*s1(ix,iy,iz,ie)+
     $                               uny(l,1,ifc,ie)*s2(ix,iy,iz,ie)+
     $                               unz(l,1,ifc,ie)*s3(ix,iy,iz,ie))
               l=l+1
            enddo
            enddo
            enddo
         endif
      enddo
      enddo

      s=glsum(s,1)
      a=glsum(a,1)

      svint=s

      return
      end
c-----------------------------------------------------------------------
      subroutine checkaeq(ux,uy,uz,uxx,uxy,pp,visc)

      ! checker for aeq formulation for velocity

      ! <ux,uy,uz> := mean velocity field
      ! uxx        := mean correlation (<ux ux>,<uy uy>,<uz uz>) field
      ! uxy        := mean correlation (<ux uy>,<uy uz>,<uz ux>) field
      ! pp         := mean pressure field
      ! visc       := viscosity corresponding to input fields

      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'

      real ux(lx1,ly1,lz1,lelt)
      real uy(lx1,ly1,lz1,lelt)
      real uz(lx1,ly1,lz1,lelt)

      real uxx(lx1*ly1*lz1*lelt,ldim)
      real uxy(lx1*ly1*lz1*lelt,ldim)

      real pp(lx2,ly2,lz2,lelt)

      common /myscr/         xm0(lx1,ly1,lz1,lelt)
     $,                      ym0(lx1,ly1,lz1,lelt)
     $,                      zm0(lx1,ly1,lz1,lelt)
     $,                      pm1(lx1,ly1,lz1,lelv)
     $,                      wx(lx1,ly1,lz1,lelv)
     $,                      wy(lx1,ly1,lz1,lelv)
     $,                      wz(lx1,ly1,lz1,lelv)
      
      n=lx1*ly1*lz1*nelv
      time=1./visc

      call divm1(wx,uxx(1,1),uxy(1,1),uxy(1,3))
      call divm1(wy,uxy(1,1),uxx(1,2),uxy(1,2))
      if (ldim.eq.3) call divm1(wz,uxy(1,3),uxy(1,2),uxx(1,3))

      call chsign(wx,n)
      call chsign(wy,n)
      call chsign(wz,n)

      call outpost(wx,wy,wz,pp,t,'cnv')

      call mappr(pm1,pp,xm0,ym0)
      call gradm1(xm0,ym0,zm0,pm1)

      call outpost(xm0,ym0,zm0,pm1,t,'gdp')

      call sub2(wx,xm0,n)
      call sub2(wy,ym0,n)
      if (ldim.eq.3) call sub2(wz,zm0,n)

      if (ldim.eq.2) then
         call lap2d(xm0,ux)
         call add2s2(wx,xm0,visc,n)

         call lap2d(ym0,uy)
         call add2s2(wy,ym0,visc,n)
      else
         call lap3d(xm0,ux)
         call add2s2(wx,xm0,visc,n)

         call lap3d(ym0,uy)
         call add2s2(wy,ym0,visc,n)

         call lap3d(zm0,uz)
         call add2s2(wz,zm0,visc,n)
      endif

      call outpost(xm0,ym0,zm0,pp,t,'lap')

      call outpost(wx,wy,wz,pp,t,'aeq')

      return
      end
c-----------------------------------------------------------------------
      subroutine checkaeqp(ux,uy,uz,uxx,uxy,pp,visc)

      ! checker for aeq formulation for velocity (with fluc. contribution)

      ! <ux,uy,uz> := mean velocity field
      ! uxx        := mean correlation (<ux'ux'>,<uy'uy'>,<uz'uz'>)
      ! uxy        := mean correlation (<ux'uy'>,<uy'uz'>,<uz'ux'>)
      ! pp         := mean pressure field
      ! visc       := viscosity corresponding to input fields

      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'

      real ux(lx1,ly1,lz1,lelt)
      real uy(lx1,ly1,lz1,lelt)
      real uz(lx1,ly1,lz1,lelt)

      real uxx(lx1*ly1*lz1*lelt,ldim)
      real uxy(lx1*ly1*lz1*lelt,ldim)

      real pp(lx2,ly2,lz2,lelt)

      common /myscr/         xm0(lx1,ly1,lz1,lelt)
     $,                      ym0(lx1,ly1,lz1,lelt)
     $,                      zm0(lx1,ly1,lz1,lelt)
     $,                      pm1(lx1,ly1,lz1,lelv)
     $,                      wx(lx1,ly1,lz1,lelv)
     $,                      wy(lx1,ly1,lz1,lelv)
     $,                      wz(lx1,ly1,lz1,lelv)
      
      n=lx1*ly1*lz1*nelv
      time=1./visc

      call divm1(wx,uxx(1,1),uxy(1,1),uxy(1,3))
      call divm1(wy,uxy(1,1),uxx(1,2),uxy(1,2))
      if (ldim.eq.3) call divm1(wz,uxy(1,3),uxy(1,2),uxx(1,3))

      call chsign(wx,n)
      call chsign(wy,n)
      call chsign(wz,n)

      call col3(xm0,ux,ux,n)
      call col3(ym0,ux,uy,n)
      call col3(zm0,ux,uz,n)
      call divm1(pm1,xm0,ym0,zm0)
      call sub2(wx,pm1,n)

      call col3(xm0,uy,ux,n)
      call col3(ym0,uy,uy,n)
      call col3(zm0,uy,uz,n)
      call divm1(pm1,xm0,ym0,zm0)
      call sub2(wy,pm1,n)

      if (ldim.eq.3) then
         call col3(xm0,uz,ux,n)
         call col3(ym0,uz,uy,n)
         call col3(zm0,uz,uz,n)
         call divm1(pm1,xm0,ym0,zm0)
         call sub2(wz,pm1,n)
      endif

      call mappr(pm1,pp,xm0,ym0)
      call gradm1(xm0,ym0,zm0,pm1)

      call sub2(wx,xm0,n)
      call sub2(wy,ym0,n)
      if (ldim.eq.3) call sub2(wz,zm0,n)

      if (ldim.eq.2) then
         call lap2d(xm0,ux)
         call add2s2(wx,xm0,visc,n)

         call lap2d(ym0,uy)
         call add2s2(wy,ym0,visc,n)
      else
         call lap3d(xm0,ux)
         call add2s2(wx,xm0,visc,n)

         call lap3d(ym0,uy)
         call add2s2(wy,ym0,visc,n)

         call lap3d(zm0,uz)
         call add2s2(wz,zm0,visc,n)
      endif

      call outpost(wx,wy,wz,pp,t,'aq2')

      return
      end
c-----------------------------------------------------------------------
      subroutine checkaeqd(ux,uy,uz,uxx,uxy,pp,ud,visc)

      ! checker for aeq formulation for velocity (with eddy model)

      ! <ux,uy,uz> := mean velocity field
      ! uxx        := mean correlation (<ux'ux'>,<uy'uy'>,<uz'uz'>)
      ! uxy        := mean correlation (<ux'uy'>,<uy'uz'>,<uz'ux'>)
      ! pp         := mean pressure field
      ! visc       := viscosity corresponding to input fields

      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'

      real ux(lx1*ly1*lz1*lelt)
      real uy(lx1*ly1*lz1*lelt)
      real uz(lx1*ly1*lz1*lelt)

      real uxx(lx1*ly1*lz1*lelt,ldim)
      real uxy(lx1*ly1*lz1*lelt,ldim)

      real pp(lx2,ly2,lz2,lelt)

      real ud(lx1*ly1*lz1*lelt,ldim)

      common /myscr/         xm0(lx1,ly1,lz1,lelt)
     $,                      ym0(lx1,ly1,lz1,lelt)
     $,                      zm0(lx1,ly1,lz1,lelt)
     $,                      pm1(lx1,ly1,lz1,lelv)
     $,                      wx(lx1,ly1,lz1,lelv)
     $,                      wy(lx1,ly1,lz1,lelv)
     $,                      wz(lx1,ly1,lz1,lelv)
      
      n=lx1*ly1*lz1*nelv
      time=1./visc

      call divm1(wx,uxx(1,1),uxy(1,1),uxy(1,3))
      call divm1(wy,uxy(1,1),uxx(1,2),uxy(1,2))
      if (ldim.eq.3) call divm1(wz,uxy(1,3),uxy(1,2),uxx(1,3))

      call chsign(wx,n)
      call chsign(wy,n)
      call chsign(wz,n)

      call outpost(wx,wy,wz,pp,t,'cv1')

      call gradm1(xm0,ym0,zm0,ux)
      call opcolv(xm0,ym0,zm0,ud(1,1))
      call divm1(wx,xm0,ym0,zm0)

      call gradm1(xm0,ym0,zm0,uy)
      call opcolv(xm0,ym0,zm0,ud(1,2))
      call divm1(wy,xm0,ym0,zm0)

      if (ldim.eq.3) then
         call gradm1(xm0,ym0,zm0,uz)
         call opcolv(xm0,ym0,zm0,ud(1,3))
         call divm1(wz,xm0,ym0,zm0)
      endif

      call outpost(wx,wy,wz,pp,t,'cv2')

      call col3(xm0,ux,ux,n)
      call col3(ym0,ux,uy,n)
      call col3(zm0,ux,uz,n)
      call divm1(pm1,xm0,ym0,zm0)
      call sub2(wx,pm1,n)

      call col3(xm0,uy,ux,n)
      call col3(ym0,uy,uy,n)
      call col3(zm0,uy,uz,n)
      call divm1(pm1,xm0,ym0,zm0)
      call sub2(wy,pm1,n)

      if (ldim.eq.3) then
         call col3(xm0,uz,ux,n)
         call col3(ym0,uz,uy,n)
         call col3(zm0,uz,uz,n)
         call divm1(pm1,xm0,ym0,zm0)
         call sub2(wz,pm1,n)
      endif

      call mappr(pm1,pp,xm0,ym0)
      call gradm1(xm0,ym0,zm0,pm1)

      call sub2(wx,xm0,n)
      call sub2(wy,ym0,n)
      if (ldim.eq.3) call sub2(wz,zm0,n)

      if (ldim.eq.2) then
         call lap2d(xm0,ux)
         call add2s2(wx,xm0,visc,n)

         call lap2d(ym0,uy)
         call add2s2(wy,ym0,visc,n)
      else
         call lap3d(xm0,ux)
         call add2s2(wx,xm0,visc,n)

         call lap3d(ym0,uy)
         call add2s2(wy,ym0,visc,n)

         call lap3d(zm0,uz)
         call add2s2(wz,zm0,visc,n)
      endif

      call outpost(wx,wy,wz,pp,t,'aq2')

      return
      end
c-----------------------------------------------------------------------
      subroutine setdiff(udfld,tdfld,uafld,tafld,upup,upvp,uptp)

      ! checker for aeq formulation for velocity (with eddy model)

      ! <ux,uy,uz> := mean velocity field
      ! uxx        := mean correlation (<ux'ux'>,<uy'uy'>,<uz'uz'>)
      ! uxy        := mean correlation (<ux'uy'>,<uy'uz'>,<uz'ux'>)
      ! pp         := mean pressure field
      ! visc       := viscosity corresponding to input fields

      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'

      real udfld(lx1*ly1*lz1*lelt,ldim)
      real tdfld(lx1*ly1*lz1*lelt)
      real uafld(lx1*ly1*lz1*lelt,ldim),tafld(lx1*ly1*lz1*lelt)
      real upup(lx1*ly1*lz1*lelt,ldim),upvp(lx1*ly1*lz1*lelt,ldim)
      real uptp(lx1*ly1*lz1*lelt,ldim)

      common /myscr/         xm0(lx1,ly1,lz1,lelt)
     $,                      ym0(lx1,ly1,lz1,lelt)
     $,                      zm0(lx1,ly1,lz1,lelt)
     $,                      pm1(lx1,ly1,lz1,lelv)
     $,                      tx(lx1*ly1*lz1*lelv)
     $,                      ty(lx1*ly1*lz1*lelv)
     $,                      tz(lx1*ly1*lz1*lelv)

      ! evaluate eddy viscosity

      call gradm1(tx,ty,tz,uafld(1,1))
      call projvecm(upup(1,1),upvp(1,1),upvp(1,3),tx,ty,tz,udfld(1,1))

      call gradm1(tx,ty,tz,uafld(1,2))
      call projvecm(upvp(1,1),upup(1,2),upvp(1,2),tx,ty,tz,udfld(1,2))

      if (ldim.eq.3) then
         call gradm1(tx,ty,tz,uafld(1,3))
         call projvecm(upvp(1,3),upvp(1,2),upup(1,3),
     $      tx,ty,tz,udfld(1,3))
      endif

      ! evaluate eddy diffusivity

      call gradm1(tx,ty,tz,tafld)
      call projvecm(uptp(1,1),uptp(1,2),uptp(1,ldim),tx,ty,tz,tdfld)

      call gradm1(tx,ty,tz,uafld(1,1))
      call opcolv(tx,ty,tz,udfld(1,1))
      call divm1(xm0,tx,ty,tz)

      call gradm1(tx,ty,tz,uafld(1,2))
      call opcolv(tx,ty,tz,udfld(1,2))
      call divm1(ym0,tx,ty,tz)

      call gradm1(tx,ty,tz,uafld(1,3))
      call opcolv(tx,ty,tz,udfld(1,3))
      call divm1(zm0,tx,ty,tz)

      call gradm1(tx,ty,tz,tafld)
      call opcolv(tx,ty,tz,tdfld)
      call divm1(pm1,tx,ty,tz)

      call outpost(xm0,ym0,zm0,pr,pm1,'dif')

      call divm1(xm0,upup(1,1),upvp(1,1),upvp(1,3))
      call divm1(ym0,upvp(1,1),upup(1,2),upvp(1,2))
      call divm1(zm0,upvp(1,3),upvp(1,2),upup(1,3))
      call divm1(pm1,uptp(1,1),uptp(1,2),uptp(1,3))

      call outpost(xm0,ym0,zm0,pr,pm1,'dif')

      call outpost(upup(1,1),upvp(1,1),upvp(1,3),pr,pm1,'uvw')
      call outpost(upvp(1,1),upup(1,2),upvp(1,2),pr,pm1,'uvw')
      call outpost(upvp(1,3),upvp(1,2),upup(1,3),pr,pm1,'uvw')
      call outpost(uptp(1,1),uptp(1,2),uptp(1,3),pr,pm1,'uvw')

c     call outpost(udfld(1,1),udfld(1,2),udfld(1,ldim),
c    $     pr,tdfld,'dif')

      return
      end
c-----------------------------------------------------------------------
      subroutine projvecm(ax,ay,az,bx,by,bz,pfld)

      ! find the negative of the projection scale

      ! <ax,ay,az> := vector to be projected
      ! <bx,by,bz> := vector to project onto
      ! pfld       := -a.b/b.b

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ax(lt),ay(lt),az(lt),bx(lt),by(lt),bz(lt)
      real pfld(lx1*ly1*lz1*lelt)
      
      n=lx1*ly1*lz1*nelt

      eps=1.e-16

      do i=1,n
         bnorm2=bx(i)**2+by(i)**2+eps
         if (ldim.eq.3) bnorm2=bnorm2+bz(i)**2

         dot=ax(i)*bx(i)+ay(i)*by(i)+eps
         if (ldim.eq.3) dot=dot+az(i)*bz(i)

         pfld(i)=-dot/bnorm2
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setpsi(psit,veca,vecb)

      ! find the transformation tensor such that psi a = b

      ! <ax,ay,az> := starting vector
      ! <bx,by,bz> := ending vector
      ! psi        := transformation tensor (ldim x ldim)

      include 'SIZE'

      common /psitmp/ w1(ldim)

      real veca(lx1*ly1*lz1*lelt,ldim),vecb(lx1*ly1*lz1*lelt,ldim)
      real psit(lx1*ly1*lz1*lelt,ldim*ldim)
      
      n=lx1*ly1*lz1*nelt

      eps=1.e-16
      s1=0.
      s2=0.

      do i=1,n
         anorm2=veca(i,1)**2+veca(i,2)**2+eps
         if (ldim.eq.3) anorm2=anorm2+veca(i,3)**2

         bnorm2=vecb(i,1)**2+vecb(i,2)**2+eps
         if (ldim.eq.3) bnorm2=bnorm2+vecb(i,3)**2

         dot=veca(i,1)*vecb(i,1)+veca(i,2)*vecb(i,2)
         if (ldim.eq.3) dot=dot+veca(i,3)*vecb(i,3)

         theta=acos(dot/sqrt(bnorm2*anorm2))
         ct=cos(theta)
         st=sin(theta)

         if (ldim.eq.2) then
            psit(i,1)=ct
            psit(i,2)=st
            psit(i,3)=ct
            psit(i,4)=-st
         else
            rx=veca(i,2)*vecb(i,3)-veca(i,3)*vecb(i,2)
            ry=veca(i,3)*vecb(i,1)-veca(i,1)*vecb(i,3)
            rz=veca(i,1)*vecb(i,2)-veca(i,2)*vecb(i,1)

            psit(i,1)=ct+rx*rx*(1-ct)
            psit(i,2)=ry*rx*(1-ct)+rz*st
            psit(i,3)=rz*rx*(1-ct)-ry*st
            psit(i,4)=rx*ry*(1-ct)-rz*st
            psit(i,5)=ct+ry*ry*(1-ct)
            psit(i,6)=rz*ry*(1-ct)+rx*st
            psit(i,7)=rx*rz*(1-ct)+ry*st
            psit(i,8)=ry*rz*(1-ct)-rx*st
            psit(i,9)=ct+rz*rz*(1-ct)
         endif

         call mxm(psit,ldim,veca(i,1),ldim,w1,1)
         call sub2(w1,vecb)

         s2=s2+bnorm2
         s1=s1+vlsc2(w1,w1,ldim)
      enddo

      s1=sqrt(glsum(s1,1))
      s2=sqrt(glsum(s2,1))

      if (nid.eq.0) write (6,*) s1,s2,s1/s2,'psierror'


      return
      end
c-----------------------------------------------------------------------
      subroutine splitvec(x,y,z,xyz,ndim,n)

      ! split vector list into component lists

      ! <x,y,z> := vector component lists
      ! xyz     := vector lies

      real x(n),y(n),z(n),xyz(ndim,n)

      do i=1,n
         x(i) = xyz(1,i)
         y(i) = xyz(2,i)
         if (ndim.eq.3) z(i) = xyz(3,i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine combvec(xyz,x,y,z,ndim,n)

      ! combine vector component lists into vector list

      ! xyz     := vector lies
      ! <x,y,z> := vector component lists

      real xyz(ndim,n),x(n),y(n),z(n)

      do i=1,n
         xyz(1,i) = x(i)
         xyz(2,i) = y(i)
         if (ndim.eq.3) xyz(3,i) = z(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine tfunc(u,rad,wt,n)

      ! apply transfer function defined by rad and wt to a vector

      ! u   := input/output coefficient set
      ! rad := radius of the damping (percentage of trailing modes)
      ! wt  := weight of the damping (percent attenuation for mode n)
      ! n   := length of coefficient set

      real u(n)

      m=n*rad

      do i=0,m-1
         u(n-i)=u(n-i)*(1.0-wt*(m*1.0-i)/(m*1.0))
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine evalf(f,uu,tt,g,mdim,ifweak)

      ! evaluate f, the time-derivative based on the current solution

      ! f      := time-derivative of the solution tt
      ! uu     := advection field
      ! tt     := solution field
      ! g      := forcing field
      ! gnn    := inhomogeneous Neumann value field
      ! mdim   := dimension of tt
      ! ifweak := flag to evaluate the terms in a weak sense

      include 'SIZE'  ! dep: SOLN,MOR,INPUT,lt,ltd,ldim

      include 'TSTEP' ! dep: ifield
      include 'MASS'  ! dep: bintm1
      include 'SOLN'  ! dep: prlag

      parameter (lt=lx1*ly1*lz1*lelt,ltd=lxd*lyd*lzd*lelt)

      common /scrns/ wk(lt,ldim)
      common /convect/ c1v(ltd),c2v(ltd),c3v(ltd),
     $                 u1v(ltd),u2v(ltd),u3v(ltd)

      real f(lt,mdim),uu(lt,ldim),tt(lt,mdim),g(lt,mdim)
      logical ifweak

      nv=lx1*ly1*lz1*nelv
      nt=lx1*ly1*lz1*nelt

      jfield=ifield

      ifield=1
      n=nv

      if (mdim.ne.ldim) ifield=2
      if (mdim.ne.ldim) n=nt

      call rzero(f,nt)

      call evalcflds(f,uu,tt,mdim,1,ifweak)
      call evalaflds(wk,tt,mdim,1,ifweak)

      do idim=1,mdim
         call chsign(f(1,idim),n)
         call add2(f(1,idim),wk(1,idim),n)
         call add2(f(1,idim),g(1,idim),n)
      enddo

      if (ifweak) then
         if (ifield.eq.1) then
            call opbinv1(f(1,1),f(1,2),f(1,ldim),
     $         f(1,1),f(1,2),f(1,ldim),1.)
         else
            call col2(f,tmask,n)
            call dssum(f,lx1,ly1,lz1)
            call col2(f,bintm1,n)
         endif
      else
         do idim=1,mdim
            call dsavg(f(1,idim))
         enddo
      endif

      if (ifield.eq.1) call incomprn(f(1,1),f(1,2),f(1,ldim),prlag)

      ifield=jfield

      return
      end
c-----------------------------------------------------------------------
      subroutine evalcflds(cfld,us0,ts0,mdim,ns,ifweak)

      ! evaluate convection field

      ! cfld   := result of convection application
      ! us0    := velocity fields
      ! ts0    := scalar/vector advectee fields
      ! mdim   := dimension of ts0 at each point
      ! ns     := number of convection fields to evaluate
      ! ifweak := weak form of the convection term

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrkk3/ wk(lt,3)

      real us0(lt,ldim,ns),ts0(lt,mdim,ns)
      real cfld(lt,mdim,ns)

      logical ifweak

      nv=lx1*ly1*lz1*nelv

      if (ifweak) then
         do is=1,ns
            call setcnv_c(us0(1,1,is),us0(1,2,is),us0(1,ldim,is))

            if (mdim.eq.ldim) then
               call setcnv_u(ts0(1,1,is),ts0(1,2,is),ts0(1,ldim,is))
            else
               call setcnv_u1(ts0(1,1,is))
            endif

            call cc(cfld(1,1,is),mdim)
         enddo
      else
         do is=1,ns
            do idim=1,mdim
               call gradm1(wk(1,1),wk(1,2),wk(1,3),ts0(1,idim,is))
               call col3(cfld(1,idim,is),us0(1,1,is),wk(1,1),nv)
               call add2col2(cfld(1,idim,is),us0(1,2,is),wk(1,2),nv)
               if (ldim.eq.3)
     $            call add2col2(cfld(1,idim,is),us0(1,3,is),wk(1,3),nv)
            enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine evalaflds(afld,ts0,mdim,ns,ifweak)

      ! evaluate diffusion field

      ! afld   := result of Laplacian application
      ! ts0    := scalar/vector fields
      ! mdim   := dimension of ts0 at each point
      ! ns     := number of diffusion fields to evaluate
      ! ifweak := weak form of the diffusion term

      include 'SIZE'  ! dep: lt, etc.
      include 'LMOR'  ! dep: lelm
      include 'GEOM'  ! dep: area
      include 'INPUT' ! dep: cbc
      include 'SOLN'  ! dep: vdiff

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrkk3/ wk(lt,2),wk3(lx1,ly1,lz1,lelt)
      common /morid/ ones(lx1*ly1*lz1*lelm),zeros(lx1*ly1*lz1*lelm)
      common /morbc/ gn(lx1,ly1,lz1,lelt)

      real ts0(lt,mdim,ns),afld(lt,mdim,ns)

      logical ifweak

      nv=lx1*ly1*lz1*nelv
      nt=lx1*ly1*lz1*nelt

      jfield=ifield
      ifield=1
      n=nv

      if (mdim.ne.ldim) ifield=2
      if (mdim.ne.ldim) n=nt

      if (ifweak) then
         do is=1,ns
         do idim=1,mdim
            call axhelm(afld(1,idim,is),ts0(1,idim,is),
     $         vdiff(1,1,1,1,ifield),zeros,ifield,1)
            call chsign(afld(1,idim,is),n)
         enddo
         enddo

         if (ifield.eq.2) then
            call rzero(wk3,n)
            do ie=1,nelt
            do ifc=1,2*ldim
               if (cbc(ifc,ie,2).eq.'f  ') then
                  call facind(
     $               kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,ifc)
                  l=1
                  do iz=kz1,kz2
                  do iy=ky1,ky2
                  do ix=kx1,kx2
                     wk3(ix,iy,iz,ie)=wk3(ix,iy,iz,ie)
     $                  +gn(ix,iy,iz,ie)*area(l,1,ifc,ie)
                     l=l+1
                  enddo
                  enddo
                  enddo
               endif
            enddo
            enddo

            do is=1,ns
               call add2(afld(1,1,is),wk3,n)
            enddo
         endif
      else
         do is=1,ns
         do idim=1,mdim
            call gradm1(wk(1,1),wk(1,2),wk(1,3),ts0(1,idim,is))
            call col2(wk(1,1),vdiff(1,1,1,1,ifield),n)
            call col2(wk(1,2),vdiff(1,1,1,1,ifield),n)
            if (ldim.eq.3) call col2(
     $         wk(1,3),vdiff(1,1,1,1,ifield),n)

            call divm1(afld(1,idim,is),wk(1,1),wk(1,2),wk(1,3))
         enddo
         enddo
      endif

      ifield=jfield

      return
      end
c-----------------------------------------------------------------------
      subroutine set0flow(uvwb,nb,idirf)

      ! set the flow rate of velocity fields to be 0 in the idirf direction

      ! uvwb  := velocity fields
      ! nb    := total number of velocity fields
      ! idirf := direction of flow (must be consistent with param(54))

      include 'SIZE'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      real uvwb(lt,ldim,nb)

      parameter (kx1=lx1,ky1=ly1,kz1=lz1,kx2=lx2,ky2=ly2,kz2=lz2)

      common /cvflow_a/ vxc(kx1,ky1,kz1,lelv)
     $                , vyc(kx1,ky1,kz1,lelv)
     $                , vzc(kx1,ky1,kz1,lelv)
     $                , prc(kx2,ky2,kz2,lelv)
     $                , vdc(kx1*ky1*kz1*lelv,2)

      common /cvflow_r/ flow_rate,base_flow,domain_length,xsec
     $                , scale_vf(3)

      nv=lx1*ly1*lz1*nelv

      do ib=1,nb
         current_flow=glsc2(uvwb(1,idirf,ib),bm1,nv)/domain_length
         delta_flow = 0.-current_flow
         scale = delta_flow/base_flow

         call add2s2(uvwb(1,1,ib),vxc,scale,nv)
         call add2s2(uvwb(1,2,ib),vyc,scale,nv)
         if (ldim.eq.3) call add2s2(uvwb(1,3,ib),vzc,scale,nv)
         if (nio.eq.0) write (6,1) ib,idirf,current_flow,base_flow
      enddo

    1 format (i8,i8,1p2e13.4,' set0flow')

      return
      end
c-----------------------------------------------------------------------
      subroutine p2b(uvw,uvwb,mdim,nb,iforth,sc)

      include 'SIZE'
      parameter (lt=lx1*ly1*lz1*lelt)

      real uvw(lt,mdim),uvwb(lt,mdim,nb),sc(nb)
      logical iforth

      nv=lx1*ly1*lz1*nelv
      nt=lx1*ly1*lz1*nelt

      n=nt
      if (mdim.eq.ldim) n=nv

      call p2b_helper(sc,uvw,uvwb,mdim,nb)

      if (.not.iforth) then
         do i=1,mdim
            call rzero(uvw(1,i),n)
         enddo
      endif

      do ib=1,nb
         if (iforth) then
            const=-sc(ib)
         else
            const=sc(ib)
         endif
         do idim=1,mdim
            call add2s2(uvw(1,idim),uvwb(1,idim,ib),const,n)
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine p2b_helper(sc,uvw,uvwb,mdim,nb)

      include 'SIZE'
      parameter (lt=lx1*ly1*lz1*lelt)

      real uvw(lt,mdim),uvwb(lt,mdim,nb),sc(nb)

      do i=1,nb
         sc(i)=fip(uvw(1,1),uvwb(1,1,i),mdim)
      enddo

      return
      end
c-----------------------------------------------------------------------
      function fip(uvw,xyz,mdim)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)
      real uvw(lt,mdim),xyz(lt,mdim)

      if (mdim.eq.ldim) then
         fip=vip(uvw(1,1),uvw(1,2),uvw(1,ldim),
     $            xyz(1,1),xyz(1,2),xyz(1,ldim))
      else
         fip=sip(uvw,xyz)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine find_cerr

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrns/ ux(lt),uy(lt),uz(lt),dx(lt),dy(lt),dz(lt)
      common /scruz/ sx(lt),sy(lt),sz(lt),uwx(lt),uwy(lt),uwz(lt),
     $               tx(lt),ty(lt),tz(lt)

      real cex(nb),crom(nb)

      call rzero(cex,nb)
      call rzero(crom,nb)

      if (nio.eq.0) write (6,*) 'starting find_cerr ...'

      nv=lx1*ly1*lz1*nelv

      mint=25

      nnb=(nb/2)/mint

      do is=1,ns
         if (nio.eq.0) write (6,*) 'find_cerr loop ',is
         ! create tilde u
         call copy(ux,us0(1,1,is),nv)
         call add2(ux,ub,nv)

         call copy(uy,us0(1,2,is),nv)
         call add2(uy,vb,nv)
c        if (nio.eq.0) write (6,*) 'wp 2',is

         if (ldim.eq.3) then
            call copy(uz,us0(1,3,is),nv)
            call add2(uz,wb,nv)
         endif

         do iib=1,mint
c           if (nio.eq.0) write (6,*) 'wp 3',iib,mint
            call opzero(dx,dy,dz)

            ! create u
            call opcopy(dx,dy,dz,ub,vb,wb)
            do ib=1,iib*nnb
               cf=op_glsc2_wt(us0(1,1,is),us0(1,2,is),us0(1,ldim,is),
     $                    uvwb(1,1,ib),uvwb(1,2,ib),uvwb(1,ldim,ib),bm1)

               call add2s2(dx,uvwb(1,1,ib),cf,nv)
               call add2s2(dy,uvwb(1,2,ib),cf,nv)
               if (ldim.eq.3) call add2s2(dz,uvwb(1,3,ib),cf,nv)
            enddo
c           if (nio.eq.0) write (6,*) 'wp 4',iib,mint

            ! create u'
            call opsub3(sx,sy,sz,ux,uy,uz,dx,dy,dz)
c           if (nio.eq.0) write (6,*) 'wp 4.1',iib,mint

            call convect_new(uwx,dx,.false.,ux,uy,uz,.false.)
            call convect_new(uwy,dy,.false.,ux,uy,uz,.false.)
            if (ldim.eq.3)
     $         call convect_new(uwz,dz,.false.,ux,uy,uz,.false.)
c           if (nio.eq.0) write (6,*) 'wp 4.2',iib,mint

            ctmp1=glsc2(sx,uwx,nv)+glsc2(sy,uwy,nv)
c           if (nio.eq.0) write (6,*) 'wp 4.3',iib,mint
            if (ldim.eq.3) ctmp1=ctmp1+glsc2(sz,uwz,nv)
c           if (nio.eq.0) write (6,*) 'wp 4.4',iib,mint
            cex(iib)=cex(iib)+ctmp1*ctmp1
c           if (nio.eq.0) write (6,*) 'wp 5',iib,mint

            ! project u' onto basis

            call opzero(tx,ty,tz)
            ctmp2=0.
            do ib=nnb*iib+1,nnb*iib*2
               cf=op_glsc2_wt(sx,sy,sz,ub(1,ib),vb(1,ib),wb(1,ib),bm1)
               call add2s2(tx,ub(1,ib),cf,nv)
               call add2s2(ty,vb(1,ib),cf,nv)
               if (ldim.eq.3) call add2s2(tz,wb(1,ib),cf,nv)

               ctmp2=glsc2(tx,uwx,nv)+glsc2(ty,uwy,nv)
               if (ldim.eq.3) ctmp2=ctmp2+glsc2(sz,uwz,nv)

c              if (nio.eq.0) write (6,*) 'ctmp ',is,iib,ib,ctmp1,ctmp2

            enddo
            crom(iib)=crom(iib)+(ctmp1-ctmp2)*(ctmp1-ctmp2)
         enddo
      enddo

      if (nio.eq.0) then
         open (unit=10,file='c.dat')

         do i=1,mint
            write (10,*) cex(i),crom(i),sqrt(crom(i)/cex(i))
         enddo

         close (unit=10)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine find_cerr2

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrns/ ux(lt),uy(lt),uz(lt),dx(lt),dy(lt),dz(lt)
      common /scruz/ sx(lt),sy(lt),sz(lt),uwx(lt),uwy(lt),uwz(lt),
     $               tx(lt),ty(lt),tz(lt),
     $               ex(lt),ey(lt),ez(lt)

      real cex(nb),crom(nb)

      call rzero(cex,nb)
      call rzero(crom,nb)

      if (nio.eq.0) write (6,*) 'starting find_cerr ...'

      nv=lx1*ly1*lz1*nelv

      mint=25

      nnb=((nb-1)/2)/mint

      ! orthogonalize AB                 
      do ib=nnb*mint+2,2*nnb*mint+1
         call opcopy(ex,ey,ez,ub(1,ib),vb(1,ib),wb(1,ib))
         do jb=nnb*mint+1,ib-1
            cf=-op_glsc2_wt(ub(1,jb),vb(1,jb),wb(1,jb),
     $         ex,ey,ez,bm1)
            call add2s2(ex,ub(1,jb),cf,nv)
            call add2s2(ey,vb(1,jb),cf,nv)
            if (ldim.eq.3) call add2s2(ez,wb(1,jb),cf,nv)
         enddo

         cf=1./sqrt(op_glsc2_wt(ex,ey,ez,
     $       ex,ey,ez,bm1))
         call opcmult(ex,ey,ez,cf)
         call opcopy(ub(1,ib),vb(1,ib),wb(1,ib),ex,ey,ez,cf)
      enddo

      do is=1,ns
         if (nio.eq.0) write (6,*) 'find_cerr loop ',is
         ! create tilde u
         call copy(ux,us0(1,1,is),nv)
         call add2(ux,ub,nv)

         call copy(uy,us0(1,2,is),nv)
         call add2(uy,vb,nv)
c        if (nio.eq.0) write (6,*) 'wp 2',is

         if (ldim.eq.3) then
            call copy(uz,us0(1,3,is),nv)
            call add2(uz,wb,nv)
         endif

         do iib=1,mint
c           if (nio.eq.0) write (6,*) 'wp 3',iib,mint
            call opzero(dx,dy,dz)

            ! create u
            call opcopy(dx,dy,dz,ub,vb,wb)
            do ib=1,iib*nnb
               cf=op_glsc2_wt(us0(1,1,is),us0(1,2,is),us0(1,ldim,is),
     $                    uvwb(1,1,ib),uvwb(1,2,ib),uvwb(1,ldim,ib),bm1)

               call add2s2(dx,uvwb(1,1,ib),cf,nv)
               call add2s2(dy,uvwb(1,2,ib),cf,nv)
               if (ldim.eq.3) call add2s2(dz,uvwb(1,3,ib),cf,nv)
            enddo
c           if (nio.eq.0) write (6,*) 'wp 4',iib,mint

            ! create u'
            call opsub3(sx,sy,sz,ux,uy,uz,dx,dy,dz)
c           if (nio.eq.0) write (6,*) 'wp 4.1',iib,mint

            call convect_new(uwx,dx,.false.,ux,uy,uz,.false.)
            call convect_new(uwy,dy,.false.,ux,uy,uz,.false.)
            if (ldim.eq.3)
     $         call convect_new(uwz,dz,.false.,ux,uy,uz,.false.)
c           if (nio.eq.0) write (6,*) 'wp 4.2',iib,mint

            ctmp1=glsc2(sx,uwx,nv)+glsc2(sy,uwy,nv)
c           if (nio.eq.0) write (6,*) 'wp 4.3',iib,mint
            if (ldim.eq.3) ctmp1=ctmp1+glsc2(sz,uwz,nv)
c           if (nio.eq.0) write (6,*) 'wp 4.4',iib,mint
            cex(iib)=cex(iib)+ctmp1*ctmp1
c           if (nio.eq.0) write (6,*) 'wp 5',iib,mint

            ! project u' onto basis

            call opzero(tx,ty,tz)
            ctmp2=0.
            do ib=nnb*iib+1,nnb*iib*2+1
               ! orthogonalize basis
               call opcopy(ex,ey,ez,ub(1,ib),vb(1,ib),wb(1,ib))
               do jb=1,nnb*iib
                  cf=-op_glsc2_wt(ub(1,jb),vb(1,jb),wb(1,jb),
     $               ex,ey,ez,bm1)
                  call add2s2(ex,ub(1,jb),cf,nv)
                  call add2s2(ey,vb(1,jb),cf,nv)
                  if (ldim.eq.3) call add2s2(ez,wb(1,jb),cf,nv)
               enddo
               cf=1./sqrt(op_glsc2_wt(ex,ey,ez,
     $                ex,ey,ez,bm1))
               call opcmult(ex,ey,ez,cf)

               ! find component
               cf=op_glsc2_wt(sx,sy,sz,ex,ey,ez,bm1)
               call add2s2(tx,ex,cf,nv)
               call add2s2(ty,ey,cf,nv)
               if (ldim.eq.3) call add2s2(tz,ez,cf,nv)

               ctmp2=glsc2(tx,uwx,nv)+glsc2(ty,uwy,nv)
               if (ldim.eq.3) ctmp2=ctmp2+glsc2(sz,uwz,nv)

c              if (nio.eq.0) write (6,*) 'ctmp ',is,iib,ib,ctmp1,ctmp2

            enddo
            crom(iib)=crom(iib)+(ctmp1-ctmp2)*(ctmp1-ctmp2)
         enddo
      enddo

      if (nio.eq.0) then
         open (unit=10,file='c.dat')

         do i=1,mint
            write (10,*) cex(i),crom(i),sqrt(crom(i)/cex(i))
         enddo

         close (unit=10)
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine genrij()

      include 'SIZE'
      include 'TOTAL'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /myblocka/ wk(lt,ldim,ldim)

      character*127 fname 

      call blank(fname,127)
      fname='avg.list'
      fname='rms.list'
      fname='rm2.list'

      nv=lx1*ly1*lz1*nelv

      call col3(wk(1,1,1),uavg,uavg,nv)
      call col3(wk(1,1,2),vavg,uavg,nv)
      call col3(wk(1,2,2),vavg,vavg,nv)

      if (ldim.eq.3) then
         call col3(wk(1,1,3),uavg,wavg,nv)
         call col3(wk(1,2,3),vavg,wavg,nv)
         call col3(wk(1,3,3),wavg,wavg,nv)
      endif

      call sub2(wk(1,1,1),urms,nv)
      call sub2(wk(1,1,2),uvms,nv)
      call sub2(wk(1,2,2),vrms,nv)

      if (ldim.eq.3) then
         call sub2(wk(1,1,3),wums,nv)
         call sub2(wk(1,2,3),vwms,nv)
         call sub2(wk(1,3,3),wrms,nv)
      endif

      call divm1(wk(1,2,1),wk(1,1,1),wk(1,1,2),wk(1,1,3))
      call divm1(wk(1,2,1),wk(1,1,1),wk(1,1,2),wk(1,1,3))

      return
      end
c-----------------------------------------------------------------------
      subroutine orthonormb(b,mdim,nb)

      ! orthonormalize fields with modified Gram-Schmidt

      include 'SIZE'
      include 'TOTAL'

      real b(lx1*ly1*lz1*lelt,mdim,nb)

      nv=lx1*ly1*lz1*nelv
      nt=lx1*ly1*lz1*nelt

      do i=1,nb
         do j=1,i-1
            if (mdim.eq.ldim) then
               sc=-vip(b(1,1,i),b(1,2,i),b(1,ldim,i),
     $                 b(1,1,j),b(1,2,j),b(1,ldim,j))

               call opadds(b(1,1,i),b(1,2,i),b(1,ldim,i),
     $                      b(1,1,j),b(1,2,j),b(1,ldim,j),sc,nv,2)
            else if (mdim.eq.1) then
               sc=-sip(b(1,1,i),b(1,1,j))
               call add2s2(b(1,1,i),b(1,1,j),sc,nt)
            endif
         enddo
         if (mdim.eq.ldim) then
            sc=1./sqrt(vip(b(1,1,i),b(1,2,i),b(1,ldim,i),
     $                     b(1,1,i),b(1,2,i),b(1,ldim,i)))
            call opcmult(b(1,1,i),b(1,2,i),b(1,ldim,i),sc)
         else if (mdim.eq.1) then
            sc=1./sqrt(sip(b(1,1,i),b(1,1,i)))
            call cmult(b(1,1,i),sc,nt)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setab(uxyz,uvwb,uvwb0)

      ! generate augmented basis from ABM

      ! uxyz  := ABM mode
      ! uvwb  := original basis mode
      ! uvwb0 := original basis mode or 0th mode

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrns/ txyz(lt,ldim)

      real uxyz(lt,ldim),uvwb(lt,ldim),uvwb0(lt,ldim)

      call evalcflds(uxyz,uvwb,uvwb0,ldim,1,.true.)
      call evalcflds(txyz,uvwb0,uvwb,ldim,1,.true.)
      call opadd2(uxyz(1,1),uxyz(1,2),uxyz(1,ldim),
     $            txyz(1,1),txyz(1,2),txyz(1,ldim))


      call opbinv1(uxyz(1,1),uxyz(1,2),uxyz(1,ldim),
     $             uxyz(1,1),uxyz(1,2),uxyz(1,ldim),1.)

      call incomprn(uxyz(1,1),uxyz(1,2),uxyz(1,ldim),prlag)

      return
      end
c-----------------------------------------------------------------------
      subroutine rbf_setup(rbfwt,rbf_sigma,edk,anch,ns,nb)

      ! Setup RBF Method, compute sigma and weight based 
      ! on edk and anch arrays

      ! Output : rbfwt, rbf_sigma
      ! Input  : edk, anch, ns, nb

      ! rbfwt := rbf weight
      ! rbf_sigma:= variance used in rbf
      ! (anch, edk) := sample points
      ! ns := number of sample points
      ! nb := number of modes

      real rbfwt(ns,nb)   ! rbf matrix and weight
      real rbf_sigma(nb)  ! rbf parameter
      real edk(0:nb,ns)   ! rbf data value
      real anch(ns)       ! rbf data point
      integer ns,nb

      call set_rbf_sig(rbf_sigma,edk,anch,ns,nb)
      call c_rbfwt(rbfwt,rbf_sigma,edk,anch,ns,nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_rbf_sig(rbf_sigma,edk,anch,ns,nb)

      ! Compute RBF sigma based on edk and anch arrays

      ! Output : rbf_sigma
      ! Input  : edk, anch, ns, nb

      ! rbf_sigma:= variance used in rbf
      ! (anch, edk) := sample points
      ! ns := number of sample points
      ! nb := number of modes

      real rbf_sigma(nb)                 ! rbf parameter
      real edk(0:nb,ns),edkk(ns,0:nb)    ! rbf data value
      real anch(ns)                      ! rbf data point
      real rbftmp1(ns**2)                ! array used in invmat
      integer ns,nb

      real rbf_mean

      open(unit=11,file='anch')
      do i=1,ns
         read(11,*) anch(i)
      enddo
      close(11)

      ! transpose
      call transpose(edkk,ns,edk,nb+1)

      do k=1,nb
         call copy(rbftmp1,edkk(1,k),ns)
         rbf_mean = vlsum(rbftmp1,ns)/ns
         call cadd(rbftmp1,-1.0*rbf_mean,ns)
         call vsq(rbftmp1,ns)
         rbf_sigma(k) = sqrt(vlsum(rbftmp1,ns)/ns)

         write(6,*) 'rbf_sigma k:',k,rbf_sigma(k)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine c_rbfwt(rbfwt,rbf_sigma,edk,anch,ns,nb)

      ! Compute RBF weight based on edk and anch arrays

      ! Output : rbfwt
      ! Input  : rbf_sigma, edk, anch, ns, nb

      ! rbfwt := rbf weight
      ! rbf_sigma:= variance used in rbf
      ! (anch, edk) := sample points
      ! ns := number of sample points
      ! nb := number of modes

      real rbfmat(ns,ns), rbfwt(ns,nb)   ! rbf matrix and weight
      real rbf_sigma(nb)                 ! rbf parameter
      real edk(0:nb,ns),edkk(ns,0:nb)    ! rbf data value
      real anch(ns)                      ! rbf data point
      real rbftmp1(ns**2),rbftmp2(ns**2) ! array used in invmat
      integer itmp3(ns),itmp4(ns)
      integer ns,nb

      ! transpose
      call transpose(edkk,ns,edk,nb+1)

      do k=1,nb
         write(6,*) 'mode k:',k,rbf_sigma(k),'Construct RBF matrix...'

         ! setup rbf matrix
         do j=1,ns
         do i=1,ns
            rbfmat(i,j) = exp(-1.*(anch(i)-anch(j))**2/
     $                           (2*(rbf_sigma(k)**2)))
         enddo
         enddo

         call invmat(rbftmp1,rbftmp2,rbfmat,itmp3,itmp4,ns)
         call mxm(rbftmp1,ns,edkk(1,k),ns,rbfwt(1,k),1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rbf_interp(edv,test_p,rbfwt,rbf_sigma,anch,ns,nb)

      real test_p, dist_p(ns)
      real edv(nb),rbfwt(ns,nb),rbf_sigma(nb)
      real anch(ns)

      do k=1,nb
         do j=1,ns
            dist_p(j)=exp(-1.*(test_p-anch(j))**2/(2*(rbf_sigma(k)**2)))
         enddo
         edv(k) = vlsc2(rbfwt(1,k),dist_p,ns)
      enddo

      return
      end
c-----------------------------------------------------------------------
      function abm_shuffle_helper(indo,nb,nsplit)

      ! return shuffled index for a given original index

      m = nint(1.0*(nb-1)/nsplit)

      if (indo.le.m) then
         abm_shuffle_helper=(indo-1)*nsplit+2
      else if (indo.eq.m+1) then
         abm_shuffle_helper=1
      else if (indo.le.2*m+1) then
         abm_shuffle_helper=(indo-2-m)*nsplit+3
      else
         abm_shuffle_helper=(indo-2-m*2)*nsplit+4
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine abm_shuffle

      include 'SIZE'
      include 'MOR'

      do i=1,nb
         indn=abm_shuffle_helper(i,nb,max(abs(iaug),2))
         call opcopy(ub(1,indn),vb(1,indn),wb(1,indn),
     $               uvwb(1,1,i),uvwb(1,2,i),uvwb(1,ldim,i))
      enddo

      do i=1,nb
         call opcopy(uvwb(1,1,i),uvwb(1,2,i),uvwb(1,ldim,i),
     $               ub(1,i),vb(1,i),wb(1,i))
      enddo

      return
      end
c-----------------------------------------------------------------------
