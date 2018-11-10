c-----------------------------------------------------------------------
      subroutine factor3(mq,mp,mr,m)

      integer dmin,d

      n=m
      l=nint(real(n)**(1/3))

      dmin=n
      imin=-1

      do i=1,n
          d=abs(n-i**3)
          if (d.lt.dmin.and.mod(n,i).eq.0) then
              dmin=d
              imin=i
          endif
      enddo

      mp=imin
      n=n/mp

      dmin=n
      imin=-1

      do i=1,n
          d=abs(n-i*i)
          if (d.lt.dmin.and.mod(n,i).eq.0) then
              dmin=d
              imin=i
          endif
      enddo

      mq=imin
      mr=n/mq

      if (nio.eq.0) write (6,*) 'mp,mq,mr,mp*mq*mr',mp,mq,mr,mp*mq*mr

      return
      end
c-----------------------------------------------------------------------
      subroutine setpart(mps,mp,n)

      integer mps(mp)

      do i=0,mp-1
         mps(i+1)=n/mp+max(mod(n,mp)-i,0)/max(mod(n,mp)-i,1)
         mps(i+1)=mps(i+1)+mps(max(i,1))*max(i,0)/max(i,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setpart3(mps,mqs,mrs,mp,mq,mr,nb)

      integer mps(mp),mqs(mq),mrs(mr)

      call setpart(mps,mp,nb)
      call setpart(mqs,mq,nb+1)
      call setpart(mrs,mr,nb+1)

      return
      end
c-----------------------------------------------------------------------
      function i2p(i,mps)

      integer mps(1)

      i2p=0

      do while (mps(i2p+1).ne.0)
         i2p=i2p+1
         if (i.le.mps(i2p)) return
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ijk2pqr(ip,iq,ir,i,j,k,mps,mqs,mrs)

      real mps(1),mqs(1),mrs(1)

      ip=i2p(i,mps)
      iq=i2p(j+1,mqs)
      ir=i2p(k+1,mrs)

      return
      end
c-----------------------------------------------------------------------
      function ijk2pid(i,j,k,mps,mqs,mrs,mp,mq,mr)

      real mps(1),mqs(1),mrs(1)

      call ijk2pqr(ip,iq,ir,i,j,k,mps,mqs,mrs)
      ijk2pid=(ip-1)+mp*(iq-1)+mp*mq*(ir-1)

      return
      end
c-----------------------------------------------------------------------
      subroutine setrange(mps,mqs,mrs,mp,mq,mr)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      integer mps(1),mqs(1),mrs(1)

      ip=mod(nid,mp)
      iq=mod(nid/mp,mq)
      irr=   (nid/mp)/mq

      i0=mps(max(ip,1))*max(ip,0)/max(ip,1)+1
      i1=mps(ip+1)

      j0=mqs(max(iq,1))*max(iq,0)/max(iq,1)
      j1=mqs(iq+1)-1

      k0=mrs(max(irr,1))*max(irr,0)/max(irr,1)
      k1=mrs(irr+1)-1

      return
      end
c-----------------------------------------------------------------------
      subroutine ijk2l(l,i,j,k)
      
      include 'SIZE'
      include 'MOR'

      il=i-i0
      jl=j-j0
      kl=k-k0

      l=il+jl*(i1-i0+1)+kl*(i1-i0+1)*(j1-j0+1)+1

      return
      end
c-----------------------------------------------------------------------
      subroutine opadd3 (a1,a2,a3,b1,b2,b3,c1,c2,c3)

      include 'SIZE'

      real a1(1),a2(1),a3(1),b1(1),b2(1),b3(1)
      real c1(1),c2(1),c3(1)

      ntot1=lx1*ly1*lz1*nelv
      call add2(a1,b1,c1,ntot1)
      call add2(a2,b2,c2,ntot1)
      if (ldim.eq.3) call add2(a3,b3,c3,ntot1)

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

      common /scrns/ t1(lt),t2(lt),t3(lt)
      common /ctrack/ cmax(0:nb), cmin(0:nb)

      character (len=72) fmt1
      character (len=72) fmt2
      character*8 fname

      real err(0:nb)

      if (istep.eq.0) then
         call rom_init

         call gengram
         call genevec
         call genbases

         do i=0,nb
            cmax(i) = -1e10
            cmin(i) =  1e10
         enddo

         time=0.
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

         do i=0,nb
            s=-u(i,1)
            call opadds(t1,t2,t3,ub(1,i),vb(1,i),wb(1,i),s,n,2)
            err(i)=op_glsc2_wt(t1,t2,t3,t1,t2,t3,bm1)
         enddo

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
      subroutine snap_analysis

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrns/ t1(lt),t2(lt),t3(lt)
      common /ctrack/ cmax(0:nb), cmin(0:nb), cvar(0:nb)

      character (len=72) fmt1
      character (len=72) fmt2
      character*8 fname

      real err(0:nb)

      call rom_init

      call gengram
      call genevec
      call genbases

      call cfill(cmax,-1e10,nb+1)
      call cfill(cmin,1e10,nb+1)

      call load_avg

      if (nio.eq.0) write (6,*) 'generating average coefficients'
      call proj2bases(usa,ua,va,wa)

      write (fmt1,'("(", i0, "(1pe15.7),1x,a4)")') nb+1

      call rzero(cvar,nb+1)
      tkes=0

      do i=1,ns
         if (nio.eq.0) write (6,*) i,'th snapshot:'
         call proj2bases(u,us(1,i),vs(1,i),ws(1,i))

         do j=0,nb
            cvar(j)=cvar(j)+(usa(j)-u(j,1))**2
            if (u(j,1).lt.cmin(j)) cmin(j)=u(j,1)
            if (u(j,1).gt.cmax(j)) cmax(j)=u(j,1)
         enddo

         call ctke_fom(tmp,us(1,i),vs(1,i),ws(1,i))
         tkes=tkes+tmp
      enddo

      tkes=tkes/real(ns)

      if (nio.eq.0) then
         write (6,fmt1) (cmax(i),i=0,nb),'cmax'
         write (6,fmt1) (usa(i) ,i=0,nb),'cavg'
         write (6,fmt1) (cmin(i),i=0,nb),'cmin'
         write (6,fmt1) (cvar(i),i=0,nb),'cvar'
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

      common /scrns/ t1(lt),t2(lt),t3(lt)
      common /ctrack/ tlast,tdiff,
     $                cmax(0:nb),cmin(0:nb),cavg(0:nb),cvar(0:nb)
      common /scrm1/ rt1(0:nb),rt2(0:nb),rt3(0:nb)

      character (len=72) fmt1
      character (len=72) fmt2
      character*8 fname

      if (istep.eq.1) then
         call cfill(cmax,-1e10,nb+1)
         call cfill(cmin, 1e10,nb+1)
         call rzero(cavg,nb+1)
         call rzero(cvar,nb+1)
         tke=0.
         time = 0.

         tlast=time
      endif

      call add2s2(cavg,u,dt,nb+1)

      do i=0,nb
         cvar(i)=cvar(i)+dt*(usa(i)-u(i,1))**2
      enddo

      call ctke_rom(tmp,u)
      tke=tke+dt*tmp

      if (mod(max(istep,1),max(iostep,1)).eq.0) then
         deltat=time-tlast
         tlast=time

         nio = -1
         call proj2bases(u,vx,vy,vz)
         nio = nid

         do i=0,nb
            if (u(i,1).lt.cmin(i)) cmin(i)=u(i,1)
            if (u(i,1).gt.cmax(i)) cmax(i)=u(i,1)
            cvar(i)=cvar(i)
         enddo

         write (fmt1,'("(i7,", i0, "(1pe15.7),1x,a4)")') nb+2
         write (fmt2,'("(i7,", i0, "(1pe15.7),1x,a4)")') nb+3

         s=1./deltat
         call cmult(cavg,s,nb+1)
         tke=tke/deltat

         if (nio.eq.0) then
            write (6,fmt1) istep,time,(cmax(i),i=0,nb),'cmax'
            write (6,fmt1) istep,time,(u(i,1),i=0,nb),'coef'
            write (6,fmt1) istep,time,(cmin(i),i=0,nb),'cmin'
            write (6,fmt2) istep,time,deltat,(cavg(i),i=0,nb),'cavg'
            write (6,fmt2) istep,time,deltat,(cvar(i),i=0,nb),'cvar'
            write (6,'(i7,3(1pe15.7),1x,a3)')
     $                     istep,time,deltat,tke,'tke'
         endif

         tke=0.
         call rzero(cavg,nb+1)
         call rzero(cvar,nb+1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ctke_fom(tke,u1,u2,u3)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrns/ ud(lt),vd(lt),wd(lt)

      real u1(lt),u2(lt),u3(lt)

      call opsub3(ud,vd,wd,u1,u2,u3,ua,va,wa)

      tke = op_glsc2_wt(ud,vd,wd,ud,vd,wd,bm1)

      return
      end
c-----------------------------------------------------------------------
      subroutine ctke_rom(tke,coef)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real coef(0:nb), cdiff(0:nb)

      tke=0.

      do i=0,nb
         cdiff(i)=coef(i)-usa(i)
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

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrns/ ud(lt),vd(lt),wd(lt)

      call opcopy(ud,vd,wd,ua,va,wa)

      n=lx1*ly1*lz1*nelv

      do i=0,nb
         s=-usa(i)
         call opadds(ud,vd,wd,ub(1,i),vb(1,i),wb(1,i),s,n,2)
      enddo

      e0n = op_glsc2_wt(ud,vd,wd,ud,vd,wd,bm1)
     $    / op_glsc2_wt(ua,va,wa,ua,va,wa,bm1)

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
     $    / h10prod(ua,va,wa,ua,va,wa,h1,h2)

      return
      end
c-----------------------------------------------------------------------
      subroutine get_saved_fields(usave,vsave,wsave,nsave,u0)

c     This routine reads files specificed in file.list

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      parameter (lt=lx1*ly1*lz1*lelt)
      real usave(lt,nsave),vsave(lt,nsave),wsave(lt,nsave)
      real uu(lt),vv(lt),ww(lt)
      real u0(lt,3) ! Initial condtion

      ierr = 0
      if (nid.eq.0) open(77,file='file.list',status='old',err=199)
      ierr = iglmax(ierr,1)
      if (ierr.gt.0) goto 199
      n = lx1*ly1*lz1*nelt
      n2= lx2*ly2*lz2*nelt

      call opcopy(uu,vv,ww,vx,vy,vz)

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

      call opcopy(vx,vy,vz,uu,vv,ww)

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
      subroutine load_avg

c     This routine reads average files specificed in avg.list

      include 'SIZE'
      include 'MOR'
      include 'TOTAL'
      include 'ZPER'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrns/ t1(lt),t2(lt),t3(lt)

      ierr = 0

      if (nid.eq.0) open(77,file='avg.list',status='old',err=199)

      ierr = iglmax(ierr,1)

      if (ierr.gt.0) goto 199

      n = lx1*ly1*lz1*nelt

      call opcopy(t1,t2,t3,vx,vy,vz)
      call opzero(us,vs,ws)

      ttime=0.

      icount = 0

      do ipass=1,navg
         call blank(initc,127)
         initc(1) = 'done '
         if (nid.eq.0) read(77,127,end=998) initc(1)
  998    call bcast(initc,127)
  127    format(a127)

         if (indx1(initc,'done ',5).eq.0) then ! We're not done

            nfiles = 1

            call restart(nfiles)  ! Note -- time is reset.
            ttime=ttime+time

            call opadds(us,vs,ws,vx,vy,vz,time,n,2)

            icount = icount+1
         else
            goto 999
         endif
      enddo

      s=1./ttime
      call opcmult(us,vs,ws,s)

      call opcopy(vx,vy,vz,t1,t2,t3)

  999 continue  ! clean up averages
      if (nid.eq.0) close(77)

      nsave = icount ! Actual number of files read

      return

  199 continue ! exception handle for file not found
      ierr = 1
      if (nid.eq.0) ierr = iglmax(ierr,1)
      call exitti('load_avg did not find avg.list$',ierr)

      return
      end
c-----------------------------------------------------------------------
