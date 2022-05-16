c-----------------------------------------------------------------------
      subroutine rom_analysis

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /ctrack/ tlast,tdiff,tke,cmax(0:lub),cmin(0:lub),
     $                cavg(0:lub),cvar(0:lub)
      common /strack/ smax(0:lub),smin(0:lub),savg(0:lub),svar(0:lub)

      common /scrm1/ rt1(0:lub),rt2(0:lub),rt3(0:lub)

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
         cvar(i)=cvar(i)+dt*(savg(i)-u(i))**2
      enddo

      do i=0,nb
         if (u(i).lt.cmin(i)) cmin(i)=u(i)
         if (u(i).gt.cmax(i)) cmax(i)=u(i)
      enddo

      ! ctke is used to compute instantaneous TKE
      call ctke(tke,u,savg)
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
            write (6,fmt1) istep,time,(u(i),i=0,nb),'coef'
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
      subroutine fom_analysis

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrana/ t1(lt),t2(lt),t3(lt)
      common /ctrack/ cmax(0:lub), cmin(0:lub)

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

c        call gengram
c        call genevec
c        call genbases

         do i=0,nb
            cmax(i) = -1e10
            cmin(i) =  1e10
         enddo

         tlast=time
      endif

      if (mod(istep,max(iostep,1)).eq.0) then
         nio = -1
         call pv2b(u,vx,vy,vz,ub,vb,wb)
         nio = nid

         do i=0,nb
            if (u(i).lt.cmin(i)) cmin(i)=u(i)
            if (u(i).gt.cmax(i)) cmax(i)=u(i)
         enddo

         write (fmt1,'("(i7,", i0, "(1pe15.7),1x,a4)")') nb+2
         write (fmt2,'("(i7,", i0, "(1pe15.7),1x,a4)")') nb+3

         call opcopy(t1,t2,t3,vx,vy,vz)
         energy=op_glsc2_wt(t1,t2,t3,t1,t2,t3,bm1)

         n=lx1*ly1*lz1*nelv

         ttmp = time
         itmp = istep
         do i=0,nb
            s=-u(i)
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
            write (6,fmt1) istep,time,(u(i),i=0,nb),'coef'
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

      common /scrana/ t1(lt),t2(lt),t3(lt)
      common /scrss/ utmp(0:lub)

      real savg(0:lub), smax(0:lub), smin(0:lub), svar(0:lub)

      character (len=72) fmt1
      character (len=72) fmt2
      character*8 fname

      call exitti('called deprecated subroutine$',nb)

      call cfill(smax,-1e10,nb+1)
      call cfill(smin,1e10,nb+1)

      call rzero(svar,nb+1)
      call rzero(savg,nb+1)

      tkes=0

      write (fmt1,'("(", i0, "(1pe15.7),1x,a4)")') nb+1

      do i=1,ns
         if (nio.eq.0) write (6,*) i,'th snapshot:'
c        call opadd3(t1,t2,t3,
c    $      us(1,1,i),us(1,2,i),us(1,ldim,i),ub,vb,wb)
         nio = -1
         call pv2b(utmp,t1,t2,t3,ub,vb,wb)
         nio = nid
         call add2(savg,utmp,nb+1)

         do j=0,nb
            if (u(j).lt.smin(j)) smin(j)=utmp(j)
            if (u(j).gt.smax(j)) smax(j)=utmp(j)
         enddo

         ! ctke_fom is used to compute mean TKE
c        call ctke_fom(tmp,us(1,1,i),us(1,2,i),us(1,ldim,i))
         tkes=tkes+tmp
      enddo

      tkes=tkes/real(ns)

      s=1/real(ns)
      call cmult(savg,s,nb+1)

      do i=1,ns
c        call opadd3(t1,t2,t3,us(1,1,i),us(1,2,i),us(1,ldim,i),ub,vb,wb)
         nio = -1
         call pv2b(utmp,t1,t2,t3,ub,vb,wb)
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
      subroutine bases_analysis

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrana/ t1(lt),t2(lt),t3(lt)
      common /ctrack/ cmax(0:lb), cmin(0:lb)

      character (len=72) fmt1
      character (len=72) fmt2
      character*8 fname

      real err(0:nb)

      if (nio.eq.0) write (6,*) 'inside bases_analysis'

      n=lx1*ly1*lz1*nelv

      nio = -1
      call pv2b(u,vx,vy,vz,ub,vb,wb)
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
         s=-u(i)
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
      subroutine rom_sample(coef)
c     This subroutine computes the sample mean and the sample variance

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /ctrack/ cavg(0:lb),cvar(0:lb)

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
      subroutine snap_analysis

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrana/ t1(lt),t2(lt),t3(lt),t4(lt)
      common /ctrack/ cmax(0:lb), cmin(0:lb)

      integer icalld
      save    icalld
      data    icalld /0/

      character (len=72) fmt1
      character (len=72) fmt2
      character*8 fname

      real err(0:nb),err_t(0:nb)

      call rom_setup

      do j=1,ns
         istep=j
         nio = -1
         n=lx1*ly1*lz1*nelv

         call pv2b(u,us0(1,1,j),us0(1,2,j),us0(1,ldim,j)
     $                    ,ub,vb,wb)

         call add2(us0(1,1,j),ub,n)
         call add2(us0(1,2,j),vb,n)
         if (ldim.eq.3) call add2(us0(1,ldim,j),wb,n)


         if (ifpod(2)) then
            call ps2b(ut,ts0(1,j,1),tb)
            call add2(ts0(1,j,1),tb,n)
         endif

         nio = nid

         write (fmt1,'("(i7,", i0, "(1pe15.7),1x,a4)")') nb+2
         write (fmt2,'("(i7,", i0, "(1pe15.7),1x,a4)")') nb+3

         call opcopy(t1,t2,t3,us0(1,1,j),us0(1,2,j),us0(1,ldim,j))
         energy=op_glsc2_wt(t1,t2,t3,t1,t2,t3,bm1)

         if (ifpod(2)) call copy(t4,ts0(1,j,1),n)

         ttmp = time
         itmp = istep
         do i=0,nb

            s=-u(i)
            call opadds(t1,t2,t3,ub(1,i),vb(1,i),wb(1,i),s,n,2)
            ss = 0
c           err(i)=op_glsc2_wt(t1,t2,t3,t1,t2,t3,bm1)
            if (ldim.eq.3) then
               do ii=1,n
                  ss=ss+bm1(ii,1,1,1)*(t1(ii)*t1(ii)+t2(ii)*t2(ii)
     $            +t3(ii)*t3(ii))
               enddo
            else
               do ii=1,n
                  ss=ss+bm1(ii,1,1,1)*(t1(ii)*t1(ii)+t2(ii)*t2(ii))
               enddo
            endif
            err(i)=sqrt(ss)

            s=-ut(i)
            ss = 0
            if (ifpod(2)) then
               call add2s2(t4,tb(1,i,1),s,n)
               do ii=1,n
                  ss=ss+bm1(ii,1,1,1)*(t4(ii)*t4(ii))
               enddo
               err_t(i)=sqrt(ss)
c              err_t(i)=op_glsc2_wt(t4,zeros,zeros,t4,zeros,zeros,bm1)
            endif

            istep = i
            time = err(i)

            if (j == 1) then
               call outpost(t1,t2,t3,pr,t4,'err')
            endif

         enddo
         time = ttmp
         istep = itmp

         if (nio.eq.0) then
            write (6,fmt1) istep,time,(u(i),i=0,nb),'coef'
            write (6,fmt2) istep,time,energy,(err(i),i=0,nb),'erru'
            write (6,fmt1) istep,time,(ut(i),i=0,nb),'coef'
            write (6,fmt2) istep,time,energy,(err_t(i),i=0,nb),
     $      'errt'
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine tj_analysis

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lint=128)

      common /my_intp/ uint(lint),vint(lint),wint(lint),tint(lint)

      n=lx1*ly1*lz1*nelv

      npts=128

      do i=0,nb
         call sol_intp_xline_qoi(ub(1,i),vb(1,i),wb(1,i),tb(1,i,1),
     $                           2.5,0.,128,i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine tj_analysis2

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lint=128)
      parameter (lmax=128)

      common /my_intp/ uint(lint),vint(lint),wint(lint),tint(lint)

      character *127 fname
      character *3 pfx

      real xxi(lmax)
      real uxi(lmax),uyi(lmax),uzi(lmax),uti(lmax)

      n=lx1*ly1*lz1*nelv

      npts=128
      nx=128

      if (nid.eq.0) then
         do j=0,nb
            uxi(i)=0.
            uyi(i)=0.
            uzi(i)=0.
            uti(i)=0.
            call blank(fname,127)
            write (fname,'(A1,I0,A4)') 'q',j,'.dat'
            open (unit=10,file=fname)
            do i=1,nx
               read (10,1) k,xxi(i),t1,t2,t3,t4
               uxi(i)=uxi(i)+t1*ua(j)
               uyi(i)=uyi(i)+t2*ua(j)
               uzi(i)=uzi(i)+t3*ua(j)
               uti(i)=uti(i)+t4*uta(j)
            enddo
            close (unit=10)
         enddo

         do i=1,nx
            write (6,*) i,xxi(i),uxi(i),uyi(i),uti(i)
         enddo
      endif

    1 format (i5,1p5e16.8,1x,'intp_resultttt')

      return
      end
c-----------------------------------------------------------------------
      subroutine tj_analysis3

      include 'SIZE'
      include 'TOTAL'
      include 'AVG'
      include 'MOR'

      character*127 fname

      parameter (lint=128)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'tj_analysis3 starting'

      call blank(fname,127)
      fname='a.list'
c     call real_averager(fname)
      call opcopy(uavg,vavg,wavg,vx,vy,vz)
      call copy(tavg,t,n)

      call sol_intp_xline_fqoi(uavg,vavg,wavg,tavg,
     $                           2.5,0.,128)

      return
      end
c-----------------------------------------------------------------------
