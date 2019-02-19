c-----------------------------------------------------------------------
      subroutine opadd3 (a1,a2,a3,b1,b2,b3,c1,c2,c3)

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

      real ux(lt),uy(lt),uz(lt), coef(0:nb)

      n=lx1*ly1*lz1*nelv

      call opzero(ux,uy,uz)

      do i=0,nb
         call opadds(ux,uy,uz,ub(1,i),vb(1,i),wb(1,i),coef(i),n,2)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine fom_analysis

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk5/ t1(lt),t2(lt),t3(lt)
      common /ctrack/ cmax(0:nb), cmin(0:nb)

      integer icalld
      save    icalld
      data    icalld /0/

      character (len=72) fmt1
      character (len=72) fmt2
      character*8 fname

      real err(0:nb)

      if (icalld.eq.0) then
         icalld=1
         call rom_init_params
         call rom_init_fields

         call gengram
         call genevec
c        call genbases

         do i=0,nb
            cmax(i) = -1e10
            cmin(i) =  1e10
         enddo

         tlast=time
      endif

      if (mod(istep,max(iostep,1)).eq.0) then
         nio = -1
         call proj2bases(u,vx,vy,vz)
         nio = nid

         do i=0,nb
            if (u(i,1).lt.cmin(i)) cmin(i)=u(i,1)
            if (u(i,1).gt.cmax(i)) cmax(i)=u(i,1)
         enddo

         write (fmt1,'("(i7,", i0, "(1pe15.7),1x,a4)")') nb+2
         write (fmt2,'("(i7,", i0, "(1pe15.7),1x,a4)")') nb+3

         call opcopy(t1,t2,t3,vx,vy,vz)
         energy=op_glsc2_wt(t1,t2,t3,t1,t2,t3,bm1)

         n=lx1*ly1*lz1*nelv

         ttmp = time
         itmp = istep
         do i=0,nb
            s=-u(i,1)
            call opadds(t1,t2,t3,ub(1,i),vb(1,i),wb(1,i),s,n,2)
            err(i)=op_glsc2_wt(t1,t2,t3,t1,t2,t3,bm1)
            istep = i
            time = err(i)
            call outpost(t1,t2,t3,pr,t,'err')
         enddo
         time = ttmp
         istep = itmp

         if (nio.eq.0) then
            write (6,fmt1) istep,time,(cmax(i),i=0,nb),'cmax'
            write (6,fmt1) istep,time,(u(i,1),i=0,nb),'coef'
            write (6,fmt1) istep,time,(cmin(i),i=0,nb),'cmin'
            write (6,fmt2) istep,time,energy,(err(i),i=0,nb),'eerr'
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine sample_stats(savg,smax,smin,svar)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk5/ t1(lt),t2(lt),t3(lt)
      common /scrk9/ utmp(0:nb)

      real savg(0:nb), smax(0:nb), smin(0:nb), svar(0:nb)

      character (len=72) fmt1
      character (len=72) fmt2
      character*8 fname

      call cfill(smax,-1e10,nb+1)
      call cfill(smin,1e10,nb+1)

      call rzero(svar,nb+1)
      call rzero(savg,nb+1)

      tkes=0

      write (fmt1,'("(", i0, "(1pe15.7),1x,a4)")') nb+1

      do i=1,ns
         if (nio.eq.0) write (6,*) i,'th snapshot:'
         call opadd3(t1,t2,t3,us(1,i),vs(1,i),ws(1,i),ub,vb,wb)
         nio = -1
         call proj2bases(utmp,t1,t2,t3)
         nio = nid
         call add2(savg,utmp,nb+1)

         do j=0,nb
            if (u(j,1).lt.smin(j)) smin(j)=utmp(j)
            if (u(j,1).gt.smax(j)) smax(j)=utmp(j)
         enddo

         ! ctke_fom is used to compute mean TKE
         call ctke_fom(tmp,us(1,i),vs(1,i),ws(1,i))
         tkes=tkes+tmp
      enddo

      tkes=tkes/real(ns)

      s=1/real(ns)
      call cmult(savg,s,nb+1)

      do i=1,ns
         call opadd3(t1,t2,t3,us(1,i),vs(1,i),ws(1,i),ub,vb,wb)
         nio = -1
         call proj2bases(utmp,t1,t2,t3)
         nio = nid
         do j=0,nb
            svar(j)=svar(j)+(savg(j)-utmp(j))**2
         enddo
      enddo

      call cmult(svar,s,nb+1)

      if (nio.eq.0) then
         write (6,fmt1) (smax(i),i=0,nb),'smax'
         write (6,fmt1) (savg(i),i=0,nb),'savg'
         write (6,fmt1) (smin(i),i=0,nb),'smin'
         write (6,fmt1) (svar(i),i=0,nb),'svar'
         write (6,*)                tkes,'tkes'
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_analysis

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk5/ t1(lt),t2(lt),t3(lt)
      common /ctrack/ tlast,tdiff,tke,cmax(0:nb),cmin(0:nb),cavg(0:nb),
     $                cvar(0:nb)
      common /strack/ smax(0:nb),smin(0:nb),savg(0:nb),svar(0:nb)

      common /scrm1/ rt1(0:nb),rt2(0:nb),rt3(0:nb)

      character (len=72) fmt1
      character (len=72) fmt2
      character*8 fname

      integer icalld
      save    icalld
      data    icalld /0/

      if (icalld.eq.0) then
         icalld=1

         call sample_stats(savg,smin,smax,svar)

         call cfill(cmax,-1e10,nb+1)
         call cfill(cmin, 1e10,nb+1)
         call rzero(cavg,nb+1)
         call rzero(cvar,nb+1)

         tlast=time-dt
      endif

      call add2s2(cavg,u,dt,nb+1)

      do i=0,nb
         cvar(i)=cvar(i)+dt*(savg(i)-u(i,1))**2
      enddo

      do i=0,nb
         if (u(i,1).lt.cmin(i)) cmin(i)=u(i,1)
         if (u(i,1).gt.cmax(i)) cmax(i)=u(i,1)
      enddo

      ! ctke_rom is used to compute instantaneous TKE
      call ctke_rom(tke,u,savg)
      if (nio.eq.0) write (6,*) istep,time,tke,'ctke'

      if (mod(ad_step,max(ad_iostep,1)).eq.0) then
         istep=ad_step
         deltat=time-tlast
         tlast=time

         write (fmt1,'("(i7,", i0, "(1pe15.7),1x,a4)")') nb+2
         write (fmt2,'("(i7,", i0, "(1pe15.7),1x,a4)")') nb+3

         s=1./deltat
         call cmult(cavg,s,nb+1)
         call cmult(cvar,s,nb+1)

         if (nio.eq.0) then
            write (6,fmt1) istep,time,(cmax(i),i=0,nb),'cmax'
            write (6,fmt1) istep,time,(u(i,1),i=0,nb),'coef'
            write (6,fmt1) istep,time,(cmin(i),i=0,nb),'cmin'
            write (6,fmt2) istep,time,deltat,(cavg(i),i=0,nb),'cavg'
            write (6,fmt2) istep,time,deltat,(cvar(i),i=0,nb),'cvar'
         endif

         call rzero(cavg,nb+1)
         call rzero(cvar,nb+1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ctke_fom(tke,u1,u2,u3)

      include 'SIZE'
      include 'MOR'
      include 'MASS'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrns/ ud(lt),vd(lt),wd(lt)

      real u1(lt),u2(lt),u3(lt)

      call opsub3(ud,vd,wd,u1,u2,u3,uavg,vavg,wavg)
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
         tke=tke+b0(i,j)*cdiff(i)*cdiff(j)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ce0n(e0n)

      include 'SIZE'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrns/ ud(lt),vd(lt),wd(lt)

      call opcopy(ud,vd,wd,uavg,vavg,wavg)

      n=lx1*ly1*lz1*nelv

      do i=0,nb
         s=-usa(i)
         call opadds(ud,vd,wd,ub(1,i),vb(1,i),wb(1,i),s,n,2)
      enddo

      e0n = op_glsc2_wt(ud,vd,wd,ud,vd,wd,bm1)
     $    / op_glsc2_wt(uavg,vavg,wavg,uavg,vavg,wavg,bm1)

      return
      end
c-----------------------------------------------------------------------
      subroutine ce1n(e1n)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrns/ ud(lt),vd(lt),wd(lt)

      call opcopy(ud,vd,wd,uavg,vavg,wavg)

      n=lx1*ly1*lz1*nelv

      do i=0,nb
         s=-usa(i)
         call opadds(ud,vd,wd,ub(1,i),vb(1,i),wb(1,i),s,n,2)
      enddo

      e1n = h10prod(ud,vd,wd,ud,vd,wd,h1,h2)
     $    / h10prod(uavg,vavg,wavg,uavg,vavg,wavg,h1,h2)

      return
      end
c-----------------------------------------------------------------------
      subroutine bases_analysis

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk5/ t1(lt),t2(lt),t3(lt)
      common /ctrack/ cmax(0:nb), cmin(0:nb)

      character (len=72) fmt1
      character (len=72) fmt2
      character*8 fname

      real err(0:nb)

      if (nio.eq.0) write (6,*) 'inside bases_analysis'

      n=lx1*ly1*lz1*nelv

      nio = -1
      call proj2bases(u,vx,vy,vz)
      nio = nid

      write (fmt1,'("(i7,", i0, "(1pe15.7),1x,a4)")') nb+2
      write (fmt2,'("(i7,", i0, "(1pe15.7),1x,a4)")') nb+3

      call opcopy(t1,t2,t3,vx,vy,vz)

      if (ifvort) then
         energy=glsc3(t1,t1,bm1,n)
      else
         energy=op_glsc2_wt(t1,t2,t3,t1,t2,t3,bm1)
      endif

      ttmp = time
      itmp = istep
      istep = -1
      time = energy
      call outpost(t1,t2,t3,pr,t,'err')
      do i=0,nb
         s=-u(i,1)
         if (ifvort) then
            call add2s2(t1,ub(1,i),s,n)
            err(i)=glsc3(t1,t1,bm1,n)
         else
            call opadds(t1,t2,t3,ub(1,i),vb(1,i),wb(1,i),s,n,2)
            err(i)=op_glsc2_wt(t1,t2,t3,t1,t2,t3,bm1)
         endif
         istep = i
         time = err(i)
         call outpost(t1,t2,t3,pr,t,'err')
         if (nio.eq.0) write (6,*) i,energy,err(i),'l2err'
      enddo
      time = ttmp
      istep = itmp

      if (nio.eq.0) write (6,*) 'exiting bases_analysis'

      return
      end
c-----------------------------------------------------------------------
      subroutine csparsity

      include 'SIZE'
      include 'MOR'
      include 'TOTAL'

      real tmp(1),tmpp(1)

      if (nid.eq.0) open (unit=50,file='ops/cloc')

      tmp(1)=ncloc
      ncmax=glmax(tmp,1)

      nlocmin = lcglo/np
      npmin = np-lcglo+(lcglo/np)*np

      eps=1.e-16

      tmpp(1)=vlamax(clocal,ncloc)

      tol=glamax(tmpp,1)

      ec=eps*tol

      nc=100

      faci=eps**(1./real(nc))

      do ie=1,nc
         nnz=0
         do i=0,np-1
            mcloc = nlocmin + i / npmin
            if (nid.eq.i) then
               call copy(cltmp,clocal,mcloc)
            else
               call rzero(cltmp,mcloc)
            endif

            call gop(cltmp,ctmp,'+  ',mcloc)

            do j=1,mcloc
               if (abs(cltmp(j)).gt.tol) nnz=nnz+1
            enddo
         enddo
         tol=tol*faci
         pnz=real(nnz)/real(nb*(nb+1)**2)
         write (6,*) ie,tol,nnz,pnz,'nonzero'
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
      subroutine comp_rms

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrkk/ ux(lt),uy(lt),uz(lt)

      if (ad_step.eq.0) then
         call rzero(u2,(nb+1)**2)
      else
         do j=0,nb
         do i=0,nb
            ur(i,j)=ur(i,j)+u(i,1)*u(j,1)
         enddo
         enddo
      endif

      if (ad_step.eq.ad_nsteps) then
         n=lx1*ly1*lz1*nelv
         s=1./real(ad_nsteps)
         call cmult(ur,s,(nb+1)**2)
         call opzero(urms,vrms,wrms)
         do j=0,nb
         do i=0,nb
            call col3(ux,ub(1,i),ub(1,j),n)
            call col3(uy,vb(1,i),vb(1,j),n)
            if (ldim.eq.3) call col3(uz,wb(1,i),wb(1,j),n)
            call opadds(urms,vrms,wrms,ux,uy,uz,ur(i,j),n,2)
         enddo
         enddo
         call outpost(urms,vrms,wrms,pr,t,'rrr')
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine drago

      include 'SIZE'
      include 'TSTEP'
      include 'MOR'

      if (ifdrago) then
         if (nio.eq.0) then
            dx=vlsc2(rdgx,u,nb+1)
            dy=vlsc2(rdgy,u,nb+1)
            write (6,*) ad_step*dt,dx,'dragx'
            write (6,*) ad_step*dt,dy,'dragy'
            if (ldim.eq.3) then
               dz=vlsc2(rdgz,u,nb+1)
               write (6,*) ad_step*dt,dz,'dragz'
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
      subroutine rom_sample(coef)
c     This subroutine computes the sample mean and the sample variance

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /ctrack/ cavg(0:nb),cvar(0:nb)

      real coef(0:nb) 

      integer icalld
      save    icalld
      data    icalld /0/

      if (icalld.eq.0) then
         icalld=1
         call rzero(cavg,nb+1)
         call rzero(cvar,nb+1)
      endif

      ! sum up all the coefficients
      do i=0,nb
         cavg(i)=cavg(i)+coef(i)
         cvar(i)=cvar(i)+coef(i)**2
      enddo

      ! cavg stands for sample mean of coefficients

      if (ad_step.eq.ad_nsteps) then

         K=ad_nsteps/ad_iostep
         s=1./K
         do i=0,nb
            cavg(i)=cavg(i)*s
         enddo
         
         if (nid.eq.0) then
            write(6,*)'s',s,'ad_nsteps',ad_nsteps,'K',K
            write(6,*)'cavg'
            do i=0,nb
               write(6,*)i,cavg(i)
            enddo
         endif

         s=1./(K-1)
         do i=0,nb
            cvar(i)=s*(cvar(i)-K*(cavg(i)**2))
         enddo

         if (nid.eq.0) then
            write(6,*)'s',s,'ad_nsteps',ad_nsteps
            write(6,*)'cvar'
            do i=0,nb
               write(6,*)i,cvar(i)
            enddo
         endif

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_avg(coef)
c     This subroutine computes the reduced space coefficient of the long-time average velocity field <u>_g

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real coef(0:nb)

      integer icalld
      save    icalld
      data    icalld /0/

      if (icalld.eq.0) then
         icalld=1
         call rzero(usa,nb+1)
      endif

      ! sum up all the coefficients
      do i=0,nb
         usa(i)=usa(i)+coef(i)
      enddo

      ! compute usa which is cavg*s
      ! s = \Delta t/T-T_0
      ! NOTE: This only correct when initial conidtion is starting
      ! with snapshot
      s=1./ad_nsteps

      if (ad_step.eq.ad_nsteps) then
         if (nid.eq.0) then 
            write(6,*)'checking'
            do i=0,nb
            write(6,*)i,usa(i)
            enddo
            do i=0,nb
               usa(i)=usa(i)*s
            enddo
            write(6,*)'s',s,'ad_nsteps',ad_nsteps
            write(6,*)'usa'
            do i=0,nb
               write(6,*)i,usa(i)
            enddo
c call dumpusa(usa,nb) <- deprecated subroutine
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine snap_sample
c     This subroutine uses the snapshots (us,vs,ws)
c     and computes the corresponding coefficient by 
c     projecting on to the reduce space

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk5/ t1(lt),t2(lt),t3(lt)
      common /scrk9/ utmp(0:nb)

      real savg(0:nb), svar(0:nb)

      character (len=72) fmt1
      character*8 fname

      call rzero(svar,nb+1)
      call rzero(savg,nb+1)

      write (fmt1,'("(", i0, "(1pe15.7),1x,a4)")') nb+1

      do i=1,ns
         if (nio.eq.0) write (6,*) i,'th snapshot:'
         call opadd3(t1,t2,t3,us(1,i),vs(1,i),ws(1,i),ub,vb,wb)
         nio = -1
         call proj2bases(utmp,t1,t2,t3)
c call dumpcoef(utmp,nb,i) <- deprecated subroutine
         nio = nid
         call add2(savg,utmp,nb+1)
      enddo

      s=1/real(ns)
      call cmult(savg,s,nb+1)

      do i=1,ns
         call opadd3(t1,t2,t3,us(1,i),vs(1,i),ws(1,i),ub,vb,wb)
         nio = -1
         call proj2bases(utmp,t1,t2,t3)
         nio = nid
         do j=0,nb
            svar(j)=svar(j)+(savg(j)-utmp(j))**2
         enddo
      enddo

      s=1/real(ns-1)
      call cmult(svar,s,nb+1)

      if (nio.eq.0) then
         write (6,fmt1) (savg(i),i=0,nb),'savg'
         write (6,fmt1) (svar(i),i=0,nb),'svar'
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mtke_rom(coef)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real coef(0:nb), cdiff(0:nb)
      character*27 fname

      integer icalld
      save    icalld
      data    icalld /0/

      if (icalld.eq.0) then
         icalld=1

         ! usa stands for coefficients for <\hat{u}>_g
         if (nid.eq.0) then
            write(fname,37) 
            write(6,*) 'fname',fname
   37 format('./MOR_data/usa')
            open (unit=12,file=fname)
            read (12,*) (usa(i),i=0,nb)
            close (unit=12)
         endif

         mtke=0.
      endif

      do i=0,nb
         cdiff(i)=coef(i)-usa(i)
      enddo

      do j=0,nb
      do i=0,nb
         mtke=mtke+b0(i,j)*cdiff(i)*cdiff(j)
      enddo
      enddo

      if (ad_step.eq.ad_nsteps) then 
         write(6,*)'usa'
         do i=0,nb
            write(6,*)i,usa(i)
         enddo

         mtke = mtke/(2*(ad_nsteps/ad_iostep))
         if (nid.eq.0) write(6,*)'mtke_rom',mtke,(ad_nsteps/ad_iostep)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine testtt(coef)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real coef(0:nb), cdiff(0:nb)
      real tmp
      character*27 fname

      integer icalld
      save    icalld
      data    icalld /0/

      if (icalld.eq.0) then
         icalld=1
         if (nid.eq.0) then
            write(fname,37) 
            write(6,*) 'fname',fname
   37 format('./MOR_data/usa')
            open (unit=12,file=fname)
            read (12,*) (usa(i),i=0,nb)
            close (unit=12)
         endif
         mtke=0.
      endif

      ! Compute \sum\sum a_i a_j b0(i,j)
      do j=0,nb
      do i=0,nb
c         mtke=mtke+b0(i,j)*(coef(i)*coef(j))
c         mtke=mtke+b0(i,j)*(coef(i)*coef(j)-usa(i)*usa(j))
         mtke=mtke+b0(i,j)*(coef(i)*coef(j)-coef(i)*usa(j)
     $               -usa(i)*coef(j)+usa(i)*usa(j))
      enddo
      enddo

      if (ad_step.eq.ad_nsteps) then 
         mtke=mtke/(2*(ad_nsteps/ad_iostep))

         write(6,*)'usa'
         do i=0,nb
c         usa(i)=0
         write(6,*)i,usa(i)
         enddo

c         tmp=0.
c         write(6,*)'tmp',tmp
c         do j=0,nb
c         do i=0,nb
c            tmp=tmp+b0(i,j)*usa(i)*usa(j)
c         enddo
c         enddo

c         mtke=mtke-((ad_nsteps/ad_iostep)*tmp) 
c         write(6,*)'mtke before',mtke
c         write(6,*)'K',(ad_nsteps/ad_iostep)
c         mtke=mtke/(2*(ad_nsteps/ad_iostep)) 
c         mtke=mtke-(tmp/2) 
         if (nid.eq.0) write(6,*)'mtke_rom',mtke
      endif

      return
      end
