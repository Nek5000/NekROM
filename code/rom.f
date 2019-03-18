c-----------------------------------------------------------------------
      subroutine rom_update_v

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      integer icalld
      save    icalld
      data    icalld /0/

      common /rom_update/ rom_time

      if (icalld.eq.0) then
         call opcopy(uic,vic,wic,vx,vy,vz)
         rom_time=0.
         icalld=1
      endif

      stime=dnekclock()

      ad_step = istep
      jfield=ifield
      ifield=1

      if (ifheat) then
         if (ifflow) call exitti(
     $   'error: running rom_update_v with ifflow = .true.$',nelv)
         if (istep.eq.0) then
            call rom_setup_v
         else
            call rom_step_v
            call reconv(vx,vy,vz,u) ! reconstruct velocity to be used in h-t
         endif
      else
         call rom_setup_v

         if (nio.eq.0) write (6,*) 'starting rom_step loop',ad_nsteps

         ad_step = 1
         do i=1,ad_nsteps
            call rom_step_v
            time=time+dt
            ad_step=ad_step+1
         enddo
      endif

      ifield=jfield

      dtime=dnekclock()-stime
      rom_time=rom_time+dtime

      if (ifheat) then
         if (nio.eq.0) write (6,*) 'rom_time: ',dtime
      else
         if (nio.eq.0) write (6,*) 'rom_time: ',rom_time
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_setup_v

      include 'SIZE'
      include 'MOR'

      if (nio.eq.0) write (6,*) 'inside rom_setup'

      setup_start=dnekclock()

      call rom_init_params
      call rom_init_fields

      call setgram
      call setevec

      call setbases
      call setops

      if (ifdumpops) call dump_all

      call qoisetup

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
      call seta(av,av0,'ops/av ')
      call setb(bv,bv0,'ops/bv ')
      call setc(cvl,icvl,'ops/cv ')
      call setu
      if (ifpod(2)) then
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

      ifl2=.false.
      ips='H10'
      if (param(33).eq.0) then
         ifl2=.true.
         ips='L2 '
      endif

      ifavgic=.false.
      if (param(34).ne.0) ifavgic=.true.

      ifdumpops=.false.
      ifread=.false.
      np35=nint(param(35))
      if (np35.eq.1) then
         ifdumpops=.true.
      else if (np35.eq.2) then
         ifread=.true.
      endif

      ifdrago=.false.
      if (param(36).ne.0) ifdrago=.true.

      ifvort=.false. ! default to false for now
      ifdump=.true.
      ifrecon=.true.
      if (ifread) ifrecon=.false.
      ifpart=.false.
      ifforce=.false.
c     ifforce=.true.
      do i=0,ldimt1
         ifpod(i)=.false.
      enddo
      if (ifflow.and.ifheat) then
         call exitti('ifflow and ifheat are true...$',n)
      else 
c        ifpod(0)=.true.
         ifpod(1)=.true.
c        ifpod(2)=.true.
      endif

      call compute_BDF_coef(ad_alpha,ad_beta)

      if (nio.eq.0) then
         write (6,*) 'rp_ips        ',ips
         write (6,*) 'rp_ifl2       ',ifl2
         write (6,*) 'rp_ifforce    ',ifforce
         write (6,*) 'rp_ifread     ',ifread
         write (6,*) 'rp_ifpart     ',ifpart
         write (6,*) 'rp_ifrecon    ',ifrecon
         write (6,*) 'rp_ifdump     ',ifdump
         write (6,*) 'rp_ifvort     ',ifvort
         write (6,*) 'rp_ifdumpops  ',ifdumpops
         write (6,*) 'rp_ifavgic    ',ifavgic
         write (6,*) 'rp_ifdrago    ',ifdrago
         do i=0,ldimt1
            write (6,*) 'rp_ifpod(',i,')   ',ifpod(i)
         enddo

c        write(6,*) 'rp_if= ',if
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
         if (ifavgic) then
            call copy_sol(ub,vb,wb,pb,tb,uavg,vavg,wavg,pavg,tavg)
         else
            call copy_sol(ub,vb,wb,pb,tb,uic,vic,wic,pic,tic)
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
      subroutine setc(cl,icl)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real cux(lt),cuy(lt),cuz(lt)

      common /scrk1/ t1(lt),binv(lt),wk1(lt),wk2(lt),wk3(lt)
      common /scrcwk/ wk(lcloc)

      real cl(lcloc),icl(3,lcloc)

      character*128 fname

      conv_time=dnekclock()

      call invers2(binv,bm1,lx1*ly1*lz1*nelv)
      call rone(binv,lx1*ly1*lz1*nelv)

      if (nio.eq.0) write (6,*) 'inside setc'

      l=0
      mid = 0
      nlocmin = lcglo/np
      npmin = np-lcglo+(lcglo/np)*np
      n=lx1*ly1*lz1*nelv

      if (ifread.and.nid.eq.0) open (unit=12,file=fname)

      do k=0,nb
         if (nio.eq.0) write (6,*) 'k=',k
         if (.not.ifread) call setcnv_c(ub(1,k),vb(1,k),wb(1,k))
         do j=0,nb
            if (.not.ifread) then
               call setcnv_u(ub(1,j),vb(1,j),wb(1,j))
               if (ifield.eq.1) then
                  call ccu(cux,cuy,cuz)
               else
                  call cct(cux)
               endif
               if (ifdrago.and.ifield.eq.1) then
                  call opbinv1(wk1,wk2,wk3,cux,cuy,cuz,1.)
                  call comp_pdrag(fd2(1,j,k),wk1,wk2,wk3)
               endif
            endif
            do i=1,nb
               l=l+1
               if (.not.ifread) then
                  if (ifield.eq.1) then
                     cvltmp(l)=op_glsc2_wt(
     $                  ub(1,i),vb(1,i),wb(1,i),cux,cuy,cuz,binv)
                  else
                     cvltmp(l)=glsc2(tb(1,1,i),cux,n)
                  endif
               endif
               icvltmp(1,l) = i
               icvltmp(2,l) = j
               icvltmp(3,l) = k
               mcloc = nlocmin + mid / npmin
               if (l.eq.mcloc) then
                  if (ifread) then
                     if (nid.eq.0) then
                        read (12,*) (cvltmp(kk),kk=1,mcloc)
                     else
                        call rzero(cvltmp,mcloc)
                     endif
                  endif
                  if (ifread) call gop(cvltmp,wk,'+  ',mcloc)
                  if (nid.eq.mid) then
                     ncloc = mcloc
                     call copy(cl,cvltmp,ncloc)
                     call icopy(icl,icvltmp,ncloc*3)
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
               a0(i,j)=h10vprod(ub(1,i),vb(1,i),wb(1,i),
     $                          ub(1,j),vb(1,j),wb(1,j))
            else
               a0(i,j)=h10sprod(tb(1,ifield-1,i),tb(1,ifield-1,j))
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
               b0(i,j)=wl2vprod(ub(1,i),vb(1,i),wb(1,i),
     $                          ub(1,j),vb(1,j),wb(1,j))
            else
               b0(i,j)=wl2sprod(tb(1,ifield-1,i),tb(1,ifield-1,j))
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
      call proj2vbases(u,uic,vic,wic,ub,vb,wb)

      if (ifpod(2)) then
         call proj2sbases(ut,tic,tb)
         do i=0,nb
            write (6,*) 'ut',ut(i,1)
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

      do i=1,nb
         bg(i)=-vecprod(bgx,bgy,bgz,ub(1,i),vb(1,i),wb(1,i))
         if (nio.eq.0) write (6,*) bg(i),i,'bg'
      enddo

      call outpost(bgx,bgy,bz,pavg,tavg,'bgv')

      if (nio.eq.0) write (6,*) 'exiting setg'

      return
      end
c-----------------------------------------------------------------------
      subroutine qoisetup

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk1/ ux(lt),uy(lt),uz(lt)

      if (ifdrago) then
         call opcopy(ux,uy,uz,vx,vy,vz)

         do i=0,nb
            nio = -1
            call opcopy(vx,vy,vz,ub(1,i),vb(1,i),wb(1,i))
            call torque_calc(1.,x0,.true.,.false.)
            rdgx(i)=dragvx(1)
            rdgy(i)=dragvy(1)
            rdgz(i)=dragvz(1)
            nio = nid
            if (nio.eq.0) then
               write (6,*) i,rdgx(i),rdgy(i),rdgz(i),'dragi'
            endif
         enddo

         call opcopy(vx,vy,vz,ux,uy,uz)
      endif

      return
      end
c-----------------------------------------------------------------------
