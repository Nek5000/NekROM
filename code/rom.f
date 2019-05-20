c-----------------------------------------------------------------------
      subroutine rom_update

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      integer icalld
      save    icalld
      data    icalld /0/

      logical ifmult

      parameter (lt=lx1*ly1*lz1*lelt)

      common /romup/ rom_time
      common /poisson/ bqr(lx1*ly1*lz1*lelt)
      common /eires/ xi(lt,lres),theta(lt,lres),sigma(lres,lres)
      common /eiivar/ nres
      common /eivar/ res

      stime=dnekclock()

      if (icalld.eq.0) then
         ttime=time
         rom_time=0.
         icalld=1
         call rom_setup
         time=ttime
      endif

      ad_step = istep
      jfield=ifield
      ifield=1

      ifmult=.not.ifrom(2).and.ifheat

      if (ifrom(2).and..not.ifrom(1)) then
         eqn='HEA'
         ifield=2
         call set_sigma
      endif

      if (ifmult) then
         if (ifflow) call exitti(
     $   'error: running rom_update with ifflow = .true.$',nelv)
         if (istep.gt.0) then
            if (ifrom(2)) call rom_step_t
            call rom_step
            call reconv(vx,vy,vz,u) ! reconstruct velocity to be used in h-t
         endif
      else
         if (nio.eq.0) write (6,*) 'starting rom_step loop',ad_nsteps
         ad_step = 1
         do i=1,ad_nsteps
            time=time+dt
            if (ifrom(2)) call rom_step_t
            call rom_step
            ad_step=ad_step+1
         enddo
         icalld=0
      endif

      if (ifrom(2).and..not.ifrom(1)) then
         ifield=2
         call cres(res,sigma,theta,nres,lres)
      endif

      ifield=jfield

      dtime=dnekclock()-stime
      rom_time=rom_time+dtime

      if (ifmult) then
         if (nio.eq.0) write (6,*) 'romd_time: ',dtime
      endif

      if (.not.ifmult.or.nsteps.eq.istep) then
         call final
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_setup

      include 'SIZE'
      include 'SOLN'
      include 'MOR'
      include 'AVG'

      logical iftmp

      if (nio.eq.0) write (6,*) 'inside rom_setup'

      setup_start=dnekclock()

      n=lx1*ly1*lz1*nelt

      call opcopy(uic,vic,wic,vx,vy,vz)
      call copy(tic,t,n)

      call rom_init_params
      call rom_init_fields

      call setgram
      call setevec

      call setbases

      call setops

      if (nio.eq.0) write (6,*) 'begin setup for qoi'

      if (ifcdrag) call cvdrag_setup
      call cnuss_setup
      call cubar_setup

      if (nio.eq.0) write (6,*) 'end setup for qoi'

      if (nio.eq.0) write (6,*) 'begin range setup'

      if (ifpod(1)) call pv2k(uk,us0,ub,vb,wb)
      if (ifpod(2)) call ps2k(tk,ts0,tb)

      call asnap

      call hyperpar

      if (nio.eq.0) write (6,*) 'end range setup'

      if (ifdumpops) call dump_all

      setup_end=dnekclock()

      if (nio.eq.0) write (6,*) 'exiting rom_setup'
      if (nio.eq.0) write (6,*) 'setup_time:', setup_end-setup_start

      return
      end
c-----------------------------------------------------------------------
      subroutine asnap

      include 'SIZE'
      include 'MOR'
      include 'AVG'

      common /scrasnap/ t1(0:nb)

      call pv2b(uas,uavg,vavg,wavg,ub,vb,wb)
      call rzero(uvs,nb+1)

      do j=1,ns
         do i=0,nb
            uvs(i)=uvs(i)+(uk(i,j)-uas(i))**2
         enddo
      enddo

      s=1./real(ns)
      do i=0,nb
         uvs(i)=uvs(i)*s
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setops

      include 'SIZE'
      include 'MOR'
      include 'TSTEP'

      logical iftmp

      if (nio.eq.0) write (6,*) 'inside setops'

      jfield=ifield
      ifield=1
      call seta(au,au0,'ops/au ')
      call setb(bu,bu0,'ops/bu ')
      call setc(cul,icul,'ops/cu ')
      call setu
      if (ifrom(2)) then
         ifield=2
         call seta(at,at0,'ops/at ')
         call setb(bt,bt0,'ops/bt ')
         call setc(ctl,ictl,'ops/ct ')
      endif
      call setg
      ifield=jfield

      if (nio.eq.0) write (6,*) 'exiting setops'

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_init_params

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      if (nio.eq.0) write (6,*) 'inside rom_init_params'

      ad_nsteps=nsteps
      ad_iostep=iostep

      ad_dt = dt
      ad_re = 1/param(2)
      ad_pe = 1/param(8)

      isolve=nint(param(170))

      if (param(171).eq.0) then
         ips='L2 '
      else if (param(171).eq.1) then
         ips='H10'
      else
         ips='HLM'
      endif

      ifavg0=.false.
      if (param(172).ne.0) ifavg0=.true.

      if (ifavg0.and.(nb.eq.ls))
     $   call exitti('nb == ls results in linear dependent bases$',nb)

      if (nb.gt.ls)
     $   call exitti('nb > ls is undefined configuration$',nb)

      ifdumpops=.false.
      ifread=.false.
      np173=nint(param(173))
      if (np173.eq.1) then
         ifdumpops=.true.
      else if (np173.eq.2) then
         ifread=.true.
      endif

      ad_qstep=nint(param(180))+ad_iostep*max(1-nint(param(180)),0)

      ifctke=.false.
      if (param(181).ne.0) ifctke=.true.

      ifcdrag=.false.
      if (param(182).ne.0) ifcdrag=.true.

      inus=min(max(nint(param(183)),0),2)

      iffastc=.false.
      if (param(191).ne.0) iffastc=.true.

      iffasth=.false.
      if (param(192).ne.0.and.ips.eq.'HLM') iffasth=.true.

      ifcintp=.false.

      bux=param(193)
      buy=param(194)
      buz=param(195)

      ifavisc=.false.
      if (param(196).ne.0.) ifavisc=.true.

      do i=0,ldimt1
         ifpod(i)=.false.
         ifrom(i)=.false.
      enddo
      ifpod(1)=param(174).ge.0.
      ifpod(2)=(ifheat.and..not.ifread.and.param(174).ne.0.)
      ifrom(1)=ifpod(1)
      ifrom(2)=ifpod(2)

      ifpod(1)=ifpod(1).or.ifrom(2)

      ifvort=.false. ! default to false for now
      ifdump=((.not.ifheat).or.ifrom(2))
      ifrecon=(.not.ifread)

      ifpart=.false.
      ifforce=.false.
      bu2=bux*bux+buy*buy+buz*buz
      ifbuoy=bu2.gt.0..and.ifrom(2)
      ifforce=bu2.gt.0..and..not.(ifrom(1).and.ifrom(2))
      ifcintp=.false.

      call compute_BDF_coef(ad_alpha,ad_beta)

      if (nio.eq.0) then
         write (6,*) 'rp_isolve     ',isolve
         write (6,*) 'rp_ips        ',ips
         write (6,*) 'rp_ifavg0     ',ifavg0
         write (6,*) 'rp_ifdumpops  ',ifdumpops
         write (6,*) 'rp_ifread     ',ifread
         write (6,*) 'rp_ad_qstep   ',ad_qstep
         write (6,*) 'rp_ifctke     ',ifctke
         write (6,*) 'rp_ifcdrag    ',ifcdrag
         write (6,*) 'rp_inus       ',inus
         write (6,*) 'rp_iffastc    ',iffastc
         write (6,*) 'rp_iffasth    ',iffasth
         write (6,*) 'rp_ifavisc    ',ifavisc
         write (6,*) ' '
         write (6,*) 'rp_ifforce    ',ifforce
         write (6,*) 'rp_ifpart     ',ifpart
         write (6,*) 'rp_ifrecon    ',ifrecon
         write (6,*) 'rp_ifdump     ',ifdump
         write (6,*) 'rp_ifvort     ',ifvort
         write (6,*) 'rp_ifbuoy     ',ifbuoy
         write (6,*) 'rp_ifcintp    ',ifcintp
         do i=0,ldimt1
            write (6,*) 'rp_ifpod(',i,')   ',ifpod(i)
         enddo
         do i=0,ldimt1
            write (6,*) 'rp_ifrom(',i,')   ',ifrom(i)
         enddo
      endif

      if (nio.eq.0) write (6,*) 'exiting rom_init_params'

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_init_fields

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      logical alist

      character*128 fname1

      if (nio.eq.0) write (6,*) 'inside rom_init_fields'

      n=lx1*ly1*lz1*nelt
      jfield=ifield
      ifield=1

      call rone(wm1,n)
      call rone(ones,n)
      call rzero(zeros,n)

      ns = ls

      if (nio.eq.0) write (6,*) 'call get_saved_fields'

      if (.not.ifread) then
         fname1='file.list '
         call get_saved_fields(us0,ps,ts0,ns,timek,fname1)

         fname1='avg.list'
         inquire (file=fname1,exist=alist)
         if (alist) then
            call push_sol(vx,vy,vz,pr,t)
            call auto_averager(fname1)
            call copy_sol(uavg,vavg,wavg,pavg,tavg,vx,vy,vz,pr,t)
            call pop_sol(vx,vy,vz,pr,t)
         endif

         if (ifavg0) then
            call copy_sol(ub,vb,wb,pb,tb,uavg,vavg,wavg,pavg,tavg)
         endif

         call outpost(uavg,vavg,wavg,pavg,tavg,'avg')
         if (ifforce) then
            call cfill(bgx,bux,n)
            call cfill(bgy,buy,n)
            if (ldim.eq.3) call cfill(bgz,buz,n)
         endif
      endif

      if (ifrecon) then
         do i=1,ns
            call sub2(us0(1,1,i),ub,n)
            call sub2(us0(1,2,i),vb,n)
            if (ldim.eq.3) call sub2(us0(1,ldim,i),wb,n)
            if (ifpod(2)) call sub2(ts0(1,i),tb,n)
         enddo
      endif

      ifield=jfield

      if (nio.eq.0) write (6,*) 'exiting rom_init_fields'

      return
      end
c-----------------------------------------------------------------------
      subroutine setc(cl,icl,fname)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real cux(lt),cuy(lt),cuz(lt)

      common /scrcwk/ wk(lcloc)

      real cl(lcloc),icl(3,lcloc)

      character*128 fname
      character*128 fnlint

      if (nio.eq.0) write (6,*) 'inside setc'

      conv_time=dnekclock()

      call lints(fnlint,fname,128)

      if (iffastc) then
         jj=nb
         ntot=nb*(nb+1)*(nb+2)/2
      else
         jj=0
         ntot=nb*(nb+1)*(nb+1)
      endif

      ntp=ntot/np
      mm=ntot-(ntot/np)*np

      l=0
      mid=0

      n=lx1*ly1*lz1*nelv

      if (ifread.and.nid.eq.0) open (unit=12,file=fnlint)

      if (.not.ifread.and..not.ifaxis) then
         do i=0,nb
            call set_convect_new(c1v(1,i),c2v(1,i),c3v(1,i),
     $                           ub(1,i),vb(1,i),wb(1,i))
            if (ifield.eq.1) then
               call intp_rstd_all(u1v(1,i),ub(1,i),nelv)
               call intp_rstd_all(u2v(1,i),vb(1,i),nelv)
               if (ldim.eq.3) call intp_rstd_all(u3v(1,i),wb(1,i),nelv)
            else
               call intp_rstd_all(u1v(1,i),tb(1,i),nelv)
            endif
         enddo
      endif

      do k=0,nb
         if (nio.eq.0) write (6,*) 'k=',k
         do j=min(k,jj),nb
            if (.not.ifread) then
               if (ifield.eq.1) then
                  if (iffastc) then
                     call ccu_new(cux,cuy,cuz,j,k)
                  else
                     call ccu(cux,cuy,cuz,k,j)
                  endif
               else
                  call cct(cux,k,j)
               endif
            endif
            do i=1,nb
               l=l+1
               if (.not.ifread) then
                  if (ifield.eq.1) then
                     cultmp(l)=op_glsc2_wt(
     $                  ub(1,i),vb(1,i),wb(1,i),cux,cuy,cuz,ones)
                  else
                     cultmp(l)=glsc2(tb(1,i),cux,n)
                  endif
               endif
               icultmp(1,l) = i
               icultmp(2,l) = j
               icultmp(3,l) = k
               mcloc=ntp+max(mm-mid,0)/max(mm-mid,1)
c              if (nio.eq.0) write (6,*) l,mcloc,'mcloc'
               if (l.eq.mcloc) then
                  if (ifread) then
                     if (nid.eq.0) then
                        read (12,*) (cultmp(kk),kk=1,mcloc)
                     else
                        call rzero(cultmp,mcloc)
                     endif
                  endif
                  if (ifread) call gop(cultmp,wk,'+  ',mcloc)
                  if (nid.eq.mid) then
                     ncloc=mcloc
                     call copy(cl,cultmp,ncloc)
                     call icopy(icl,icultmp,ncloc*3)
                  endif
                  mid=mid+1
                  l = 0
               endif
            enddo
         enddo
      enddo

      if (ifread.and.nid.eq.0) close (unit=12)

      if (nio.eq.0) write (6,*) 'conv_time: ',dnekclock()-conv_time

      if (nio.eq.0) write (6,*) 'exiting setc'

      return
      end
c-----------------------------------------------------------------------
      subroutine seta(a,a0,fname)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrseta/ wk1(lt)

      real a0(0:nb,0:nb),a(nb,nb)

      character*128 fname

      if (nio.eq.0) write (6,*) 'inside seta'

      n=lx1*ly1*lz1*nelt

      if (ifread) then
         call read_serial(a0,(nb+1)**2,fname,wk1,nid)
      else
         nio=-1
         do j=0,nb ! Form the A matrix for basis function
         do i=0,nb
            if (ifield.eq.1) then
               a0(i,j)=h10vip(ub(1,i),vb(1,i),wb(1,i),
     $                        ub(1,j),vb(1,j),wb(1,j))
            else
               a0(i,j)=h10sip(tb(1,i),tb(1,j))
            endif
         enddo
         enddo
         nio=nid
      endif

      do j=1,nb
      do i=1,nb
         a(i,j)=a0(i,j)
      enddo
      enddo

      if (nio.eq.0) write (6,*) 'exiting seta'

      return
      end
c-----------------------------------------------------------------------
      subroutine setb(b,b0,fname)

      include 'SIZE'
      include 'TSTEP'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrread/ tab((nb+1)**2)

      real b(nb,nb),b0(0:nb,0:nb)

      character*128 fname

      if (nio.eq.0) write (6,*) 'inside setb'

      if (ifread) then
         call read_serial(b0,(nb+1)**2,fname,tab,nid)
      else
         mio=nio
         nio=-1
         do j=0,nb
         do i=0,nb
            if (ifield.eq.1) then
               b0(i,j)=wl2vip(ub(1,i),vb(1,i),wb(1,i),
     $                        ub(1,j),vb(1,j),wb(1,j))
            else
               b0(i,j)=wl2sip(tb(1,i),tb(1,j))
            endif
         enddo
         enddo
         nio=mio
      endif

      do j=1,nb
      do i=1,nb
         b(i,j)=b0(i,j)
      enddo
      enddo

      if (nio.eq.0) write (6,*) 'exiting setb'

      return
      end
c-----------------------------------------------------------------------
      subroutine setu

      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrsetu/ uu(lt),vv(lt),ww(lt),tt(lt)

      if (nio.eq.0) write (6,*) 'inside setu'

      n=lx1*ly1*lz1*nelv

      jfield=ifield
      ifield=1

      call opsub2(uic,vic,wic,ub,vb,wb)
      if (ifrom(1)) then
         if (ips.eq.'H10') then
            call h10pv2b(u,uic,vic,wic,ub,vb,wb)
         else if (ips.eq.'HLM') then
            call hlmpv2b(u,uic,vic,wic,ub,vb,wb)
         else
            call pv2b(u,uic,vic,wic,ub,vb,wb)
         endif
      else
         call rzero(u,(nb+1)*3)
         u(0,1)=1.
         u(0,2)=1.
         u(0,3)=1.
      endif
      call opadd2(uic,vic,wic,ub,vb,wb)

      call sub2(tic,tb,n)
      if (ifrom(2)) then
         call ps2b(ut,tic,tb)
         do i=0,nb
            if (nio.eq.0) write (6,*) 'ut',ut(i,1)
         enddo
      endif
      call add2(tic,tb,n)

      call reconv(uu,vv,ww,u)
      call recont(tt,ut)
      call outpost(uu,vv,ww,pr,tt,'rom')

      ttime=time
      jstep=istep
      time=1.
      istep=1
      call outpost(uic,vic,wic,pr,tic,'uic')
      time=2.
      istep=2
      call outpost(uu,vv,ww,pr,tt,'uic')
      call opsub2(uu,vv,ww,uic,vic,wic)
      call sub2(tt,tic,n)
      time=3.
      istep=3
      call outpost(uu,vv,ww,pr,tt,'uic')
      time=ttime
      istep=jstep
      ifield=jfield

      if (nio.eq.0) write (6,*) 'exiting setu'

      return
      end
c-----------------------------------------------------------------------
      subroutine setg

      include 'SIZE'
      include 'SOLN'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      if (nio.eq.0) write (6,*) 'inside setg'

      call rzero(bg,nb)
      call rzero(but0,(nb+1)**2)

      if (ifbuoy) then
         do j=0,nb
         do i=0,nb
            but0(i,j)=bux*sip(tb(1,j),ub(1,i))+buy*sip(tb(1,j),vb(1,i))
            if (ldim.eq.3) but0(i,j)=but0(i,j)+buz*sip(tb(1,j),wb(1,i))
         enddo
         enddo
      else if (ifforce) then
         do i=1,nb
            bg(i)=vip(bgx,bgy,bgz,ub(1,i),vb(1,i),wb(1,i))
            if (nio.eq.0) write (6,*) bg(i),i,'bg'
         enddo
         call outpost(bgx,bgy,bz,pavg,tavg,'bgv')
      endif

      if (nio.eq.0) write (6,*) 'exiting setg'

      return
      end
c-----------------------------------------------------------------------
      subroutine final

      include 'SIZE'
      include 'MOR'

      real t1(0:nb),t2(0:nb)

      if (nid.eq.0) then
         write (6,*) 'evalc_time: ',evalc_time
         write (6,*) 'lu_time:    ',lu_time
         write (6,*) 'solve_time: ',solve_time
         write (6,*) 'step_time:  ',step_time
         write (6,*) 'rom_time:   ',rom_time
      endif

      if (ifdumpops) then
         call dump_serial(u,nb+1,'ops/uf ',nid)
         call dump_serial(ua,nb+1,'ops/ua ',nid)
         if (ifrom(2)) call dump_serial(ut,nb+1,'ops/tf ',nid)
         do i=0,nb
            t1(i)=u2a(i,i)-ua(i)*ua(i)
         enddo
         call dump_serial(t1,nb+1,'ops/uv ',nid)
         call dump_serial(uas,nb+1,'ops/uas ',nid)
         call dump_serial(uvs,nb+1,'ops/uvs ',nid)
      endif

      return
      end
c-----------------------------------------------------------------------
