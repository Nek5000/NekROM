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
      subroutine recon(ux,uy,uz,coef)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ux(lt),uy(lt),uz(lt),coef(0:nb)

      n=lx1*ly1*lz1*nelv

      call opzero(ux,uy,uz)

      do i=0,nb
         call opadds(ux,uy,uz,ub(1,i),vb(1,i),wb(1,i),coef(i),n,2)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ctke_fom(tke,ux,uy,uz)

      include 'SIZE'
      include 'MOR'
      include 'MASS'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrns/ ud(lt),vd(lt),wd(lt)

      real ux(lt),uy(lt),uz(lt)

      call opsub3(ud,vd,wd,ux,uy,uz,uavg,vavg,wavg)
      tke = op_glsc2_wt(ud,vd,wd,ud,vd,wd,bm1)

      return
      end
c-----------------------------------------------------------------------
      subroutine ctke_rom(tke,coef,acoef)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real coef(0:nb), acoef(0:nb), cdiff(0:nb)

      tke=0.

      do i=0,nb
         cdiff(i)=coef(i)-acoef(i)
      enddo

      do j=0,nb
      do i=0,nb
         tke=tke+bv0(i,j)*cdiff(i)*cdiff(j)
      enddo
      enddo

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
      subroutine comp_drag

      include 'SIZE'
      include 'TSTEP'
      include 'MOR'

      real tmp(0:nb)

      if (ifdrago) then
         if (nio.eq.0) then
            vdx=vlsc2(rdgx,u,nb+1)
            vdy=vlsc2(rdgy,u,nb+1)
            pdx1=0.
            pdx2=0.
            pdx3=0.
            call mxm(u,nb+1,ad_beta(2,count),3,tmp,1)
            do j=0,nb
               pdx1=pdx1+fd1(1,j)*tmp(j)
               pdy1=pdy1+fd1(2,j)*tmp(j)
               do i=0,nb
                  pdx2=pdx2+fd2(1,i,j)*u(j,1)*u(i,1)
                  pdy2=pdy2+fd2(2,i,j)*u(j,1)*u(i,1)
               enddo
               pdx3=pdx3+fd3(1,j)*u(j,1)
               pdy3=pdy3+fd3(2,j)*u(j,1)
            enddo
            dx=pdx1+pdx2+pdx3+vdx
            dy=pdy1+pdy2+pdy3+vdy
c           if (nio.eq.0) write (6,*) ad_step*dt,vdx,pdx1,pdx2,pdx3,dx,'dragx'
c           if (nio.eq.0) write (6,*) ad_step*dt,vdy,pdy1,pdy2,pdy3,dy,'dragy'
            if (ldim.eq.3) then
               dz=vlsc2(rdgz,u,nb+1)
               write (6,*) ad_step*dt,vdz,'dragz'
            endif
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine shiftu(v)

      include 'SIZE'
      include 'MOR'

      real v(nb)

      call copy(u(1,3),u(1,2),nb)
      call copy(u(1,2),u(1,1),nb)
      call copy(u(1,1),v,nb)

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

      parameter (lt1=lx1*ly1*lz1*lelt)
      parameter (lt2=lx2*ly2*lz2*lelt)

      real vx(lt1),vy(lt1),vz(lt1),pr(lt2),t(lt1,ldimt)
      real ux(lt1),uy(lt1),uz(lt1),pp(lt2),tt(lt1,ldimt)

      call opcopy(vx,vy,vz,ux,uy,uz)
      call copy(pr,pp,lx2*ly2*lz2*nelv)

      do idim=1,ldimt
         call copy(t(1,idim),tt(1,idim),lx1*ly1*lz1*nelt)
      enddo

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

      call setcnv_c(vx,vy,vz)
      call setcnv_u(vx,vy,vz)
      call ccu(t1,t2,t3)
      call opbinv1(t1,t2,t3,t1,t2,t3,1.)
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
      subroutine setdps ! set grad p of snapshots

      include 'SIZE'
      include 'MOR'

      do i=1,ns
         call gradp(dps(1,1,i),dps(1,2,i),dps(1,ldim,i),ps(1,i))
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setdtbases

      include 'SIZE'
      include 'MOR'
      include 'AVG'

      n=lx1*ly1*lz1*nelt

      do i=1,ns
         call ctd_sol(us0(1,1,i),us0(1,2,i),us0(1,ldim,i),
     $      us(1,1,i),us(1,2,i),us(1,ldim,i),ps(1,i),ts(1,1,i))
      enddo

      call gengram(dug,us0,ns,ldim)
      call genevec(dug)

      do i=1,nb
         call opzero(ub(1,i),vb(1,i),wb(1,i),n)
      enddo

      ! ub, vb, wb, are the modes
      do j=1,ns
      do i=1,nb
         call opadds(ub(1,i),vb(1,i),wb(1,i),
     $      us0(1,1,j),us0(1,2,j),us0(1,ldim,j),evec(j,i,1),n,2)
      enddo
      enddo

      do i=1,nb
         call outpost(ub(1,i),vb(1,i),wb(1,i),pavg,tavg,'bs3')
      enddo

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
            work(i,j) = glsc2(us(1,1,i),uw,n)+glsc2(us(1,2,i),vw,n)
            if (ldim.eq.3) work(i,j)=work(i,j)+glsc2(us(1,ldim,i),ww,n)
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
