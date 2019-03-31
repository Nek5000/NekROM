c-----------------------------------------------------------------------
      subroutine rom_update_v

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      integer icalld
      save    icalld
      data    icalld /0/

      logical ifmult

      common /rom_update/ rom_time

      stime=dnekclock()

      if (icalld.eq.0) then
         rom_time=0.
         icalld=1
         call rom_setup_v
      endif

      ad_step = istep
      jfield=ifield
      ifield=1

      ifmult=.not.ifrom(2).and.ifheat

      if (ifmult) then
         if (ifflow) call exitti(
     $   'error: running rom_update_v with ifflow = .true.$',nelv)
         if (istep.gt.0) then
            call rom_step_v
            call reconv(vx,vy,vz,u) ! reconstruct velocity to be used in h-t
         endif
      else
         if (nio.eq.0) write (6,*) 'starting rom_step loop',ad_nsteps
         ad_step = 1
         do i=1,ad_nsteps
            call rom_step_v
            time=time+dt
            ad_step=ad_step+1
         enddo
         icalld=0
      endif

      ifield=jfield

      dtime=dnekclock()-stime
      rom_time=rom_time+dtime

      if (ifheat) then
         if (nio.eq.0) write (6,*) 'rom_time: ',dtime
      else
         if (nio.eq.0) write (6,*) 'evalc_time: ',evalc_time
         if (nio.eq.0) write (6,*) 'rom_time: ',rom_time
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_setup_v

      include 'SIZE'
      include 'SOLN'
      include 'MOR'

      if (nio.eq.0) write (6,*) 'inside rom_setup'

      setup_start=dnekclock()

      call opcopy(uic,vic,wic,vx,vy,vz)

      call rom_init_params
      call rom_init_fields

      call setgram
      call setevec

      call setbases
      call setops

      if (ifdumpops) call dump_all

      if (ifcdrag) call cvdrag_setup
      if (isolve.eq.1) then
         if (nio.eq.0) write(6,*) 
     $       'solving with constrained optimization'
         call comp_hyperpar
      endif

      setup_end=dnekclock()

      if (nio.eq.0) write (6,*) 'exiting rom_setup'
      if (nio.eq.0) write (6,*) 'setup_time:', setup_end-setup_start

      return
      end
c-----------------------------------------------------------------------
      subroutine setops

      include 'SIZE'
      include 'MOR'
      include 'TSTEP'

      if (nio.eq.0) write (6,*) 'inside setops'

      jfield=ifield
      ifield=1
      call seta(au,au0,'ops/au ')
      call setb(bu,bu0,'ops/bu ')
      call setc(cul,icul,'ops/cu ')
      call setu
      if (ifpod(2)) then
         ifield=2
         call seta(at,at0,'ops/at ')
         call setb(bt,bt0,'ops/bt ')
         call setc(ctl,ictl,'ops/ct ')
      endif
      if (ifcintp) then
         do j=0,nb
         do i=1,nb
            buc(i,j)=wl2vip(ub(1,i),vb(1,i),wb(1,i),
     $                      cxb(1,j),cyb(1,j),czb(1,j))
         enddo
         enddo
         call setcintp
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

      ifl2=.false.
      ips='H10'
      if (param(171).eq.0) then
         ifl2=.true.
         ips='L2 '
      endif

      ifavg0=.false.
      if (param(172).ne.0) ifavg0=.true.

      if (ifavg0.and.(nb.eq.ls))
     $   call exitti('nb == ls results in linear dependent bases$',nb)

      ifdumpops=.false.
      ifread=.false.
      np172=nint(param(173))
      if (np173.eq.1) then
         ifdumpops=.true.
      else if (np173.eq.2) then
         ifread=.true.
      endif

      ad_qstep=nint(param(180))+ad_iostep*max(1-nint(param(180)),0)

      ifcdrag=.false.
      if (param(182).ne.0) ifcdrag=.true.

      iffastc=.false.
      if (param(191).ne.0) iffastc=.true.

      ifcintp=.false.
      if (param(192).ne.0) ifcintp=.true.

      do i=0,ldimt1
         ifpod(i)=.false.
         ifrom(i)=.false.
      enddo
      ifpod(1)=.true.
      ifpod(2)=(ifheat.and..not.ifread)
      ifrom(1)=.true.

      ifvort=.false. ! default to false for now
      ifdump=((.not.ifheat).or.ifrom(2))
      ifrecon=(.not.ifread)

      ifpart=.false.
      ifforce=.false.
      ifbuoy=.false.
      ifcintp=.false.

      call compute_BDF_coef(ad_alpha,ad_beta)

      if (nio.eq.0) then
         write (6,*) 'rp_isolve     ',isolve
         write (6,*) 'rp_ifl2       ',ifl2
         write (6,*) 'rp_ifavg0     ',ifavg0
         write (6,*) 'rp_ifdumpops  ',ifdumpops
         write (6,*) 'rp_ifread     ',ifread
         write (6,*) 'rp_ad_qstep   ',ad_qstep
         write (6,*) 'rp_ifcdrag    ',ifcdrag
         write (6,*) 'rp_iffastc    ',iffastc
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

      common /scrk0/ t1(lt),t2(lt),t3(lt),u0(lx1*ly1*lz1*lelt,3)

      logical alist

      character*128 fname1

      if (nio.eq.0) write (6,*) 'inside rom_init_fields'

      n=lx1*ly1*lz1*nelt

      call rone(wm1,n)
      call rone(ones,n)
      call rzero(zeros,n)

      ns = ls

      if (nio.eq.0) write (6,*) 'call get_saved_fields'

      if (.not.ifread) then
         fname1='file.list '
         call get_saved_fields(us,ps,ts,ns,fname1)

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

         if (ifcintp) then
            call conv_sol(cxb,cyb,czb,uavg,vavg,wavg)
            do i=1,ns
               call conv_sol(cs(1,1,i),cs(1,2,i),cs(1,ldim,i),
     $                       us(1,1,i),us(1,2,i),us(1,ldim,i))
               call opsub3(cs0(1,1,i),cs0(1,2,i),cs0(1,ldim,i),
     $                     cs(1,1,i),cs(1,2,i),cs(1,ldim,i),cxb,cyb,czb)
               call outpost(cs(1,1,i),cs(1,2,i),cs(1,ldim,i),
     $                      pavg,tavg,'cnv')
            enddo
         endif

         call outpost(uavg,vavg,wavg,pavg,tavg,'avg')
         if (ifforce) call gradp(bgx,bgy,bgz,pavg)
      endif

      if (ifrecon) then
         do i=1,ns
            call sub3(us0(1,1,i),us(1,1,i),ub,n)
            call sub3(us0(1,2,i),us(1,2,i),vb,n)
            if (ldim.eq.3) call sub3(us0(1,ldim,i),us(1,ldim,i),wb,n)
            if (ifpod(2)) call sub3(ts0(1,i,1),ts(1,i,1),tb,n)
         enddo
      endif

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

      common /scrk1/ t1(lt),binv(lt),wk1(lt),wk2(lt),wk3(lt)
      common /scrcwk/ wk(lcloc),wk4(lt),wk5(lt),wk6(lt)

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

      if (.not.ifread) then
         do i=0,nb
            call set_convect_new(c1v(1,i),c2v(1,i),c3v(1,i),
     $                           ub(1,i),vb(1,i),wb(1,i))
            if (ifield.eq.1) then
               call intp_rstd_all(u1v(1,i),ub(1,i),nelv)
               call intp_rstd_all(u2v(1,i),vb(1,i),nelv)
               if (ldim.eq.3) call intp_rstd_all(u3v(1,i),wb(1,i),nelv)
            else
               call intp_rstd_all(u1v(1,i),tb(1,1,i),nelv)
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
                  call cct(cux1)
               endif
            endif
            do i=1,nb
               l=l+1
               if (.not.ifread) then
                  if (ifield.eq.1) then
                     cultmp(l)=op_glsc2_wt(
     $                  ub(1,i),vb(1,i),wb(1,i),cux,cuy,cuz,ones)
                  else
                     cultmp(l)=glsc2(tb(1,1,i),cux,n)
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
               a0(i,j)=h10sip(tb(1,ifield-1,i),tb(1,ifield-1,j))
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
               b0(i,j)=wl2sip(tb(1,ifield-1,i),tb(1,ifield-1,j))
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

      if (nio.eq.0) write (6,*) 'inside setu'

      jfield=ifield
      ifield=1
      call pv2b(u,uic,vic,wic,ub,vb,wb)

      if (ifpod(2)) then
         call ps2b(ut,tic,tb)
         do i=0,nb
            if (nio.eq.0) write (6,*) 'ut',ut(i,1)
         enddo
      endif

      call outpost(uic,vic,wic,pr,t,'uic')
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
            but0(i,j)=tbeta*sip(tb(1,1,j),vb(1,i))
         enddo
         enddo
      else if (ifforce) then
         do i=1,nb
            bg(i)=-vip(bgx,bgy,bgz,ub(1,i),vb(1,i),wb(1,i))
            if (nio.eq.0) write (6,*) bg(i),i,'bg'
         enddo
      endif

      call outpost(bgx,bgy,bz,pavg,tavg,'bgv')

      if (nio.eq.0) write (6,*) 'exiting setg'

      return
      end
c-----------------------------------------------------------------------
