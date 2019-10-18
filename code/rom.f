c-----------------------------------------------------------------------
      subroutine rom_update

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      integer icalld
      save    icalld
      data    icalld /0/

      logical ifmult,iftmp

      parameter (lt=lx1*ly1*lz1*lelt)

      common /romup/ rom_time

      stime=dnekclock()

      if (icalld.eq.0) then
         rom_time=0.
         icalld=1
         call rom_setup
      endif

      ad_step = istep
      jfield=ifield
      ifield=1

      ifmult=.not.ifrom(2).and.ifheat

      if (rmode.ne.'OFF') then
      if (ifmult) then
         if (ifflow) call exitti(
     $   'error: running rom_update with ifflow = .true.$',nelv)
         if (istep.gt.0) then
            if (ifrom(2)) call rom_step_t_legacy
            if (ifrom(1)) call rom_step_legacy
            call postu_legacy
            call postt_legacy
            call reconv(vx,vy,vz,u) ! reconstruct velocity to be used in h-t
         endif
      else
         if (nio.eq.0) write (6,*) 'starting rom_step loop',ad_nsteps
         ad_step = 1

         cts='bdfext'
c        cts='copt  '
c        cts='rk1   '
c        cts='rkmp  '
c        cts='rk4   '
c        cts='rkck  '

         if (cts.ne.'bdfext'.and.cts.ne.'copt  ') then
            tinit=time
            tfinal=ad_nsteps*ad_dt+time
            ndump=nint(1.*ad_nsteps/ad_iostep)
            idump=1
            tnext=(tfinal-time)/ndump+time
         endif

         do i=1,ad_nsteps
            if (cts.eq.'bdfext') then
               call bdfext_step
               call post
            else if (cts.eq.'copt  ') then
               if (ifrom(2)) call rom_step_t_legacy
               if (ifrom(1)) call rom_step_legacy
               call postu_legacy
               call postt_legacy
            else
               call rk_setup(cts)
               call copy(urki(1),u(1),nb)
               if (ifrom(2)) call copy(urki(nb+1),ut(1),nb)

               nrk=nb
               if (ifrom(2)) nrk=nb*2

               call rk_step(urko,rtmp1,urki,
     $            time,ad_dt,grk,rtmp2,nrk)

               call mxm(bu,nb,rtmp1,nb,rtmp2,1)
               call mxm(rtmp2,1,rtmp1,nb,err,1)
               err=sqrt(err)

               if (nio.eq.0) write (6,1)
     $            ad_step,time,ad_dt,err/ad_dt,err

               if (rktol.ne.0.) then
                  ttmp=ad_dt*.95*((rktol*ad_dt)/err)**.2
                  ad_dt=min(2.*ad_dt,max(.5*ad_dt,ttmp))
               endif

               write (6,*) ad_step,time,tnext,ad_dt,'time'

               if (time*(1.+1.e-12).gt.tnext) then
                  idump=idump+1
                  tnext=(tfinal-tinit)*idump/ndump+tinit

                  call reconv(vx,vy,vz,u)
                  call recont(t,ut)

                  ifto = .true. ! turn on temp in fld file
                  if (rmode.ne.'ON ')
     $               call outpost(vx,vy,vz,pr,t,'rom')
               endif

               if (cts.eq.'rkck'.and.rktol.ne.0.) then
               if (time+ad_dt.gt.(1.+1.e-12)*tnext) then
                  ttmp=ad_dt
                  ad_dt=tnext-time
                  write (6,3) ad_step,time,ad_dt,ttmp,tnext
               endif
               endif

               call copy(u(1),urko,nb)
               call copy(ut(1),urko(nb+1),nb)

               if (time*(1.+1.e-14).gt.tfinal) goto 2
            endif

            time=time+ad_dt
            ad_step=ad_step+1
         enddo
         icalld=0
      endif
      endif

    2 continue

      if (ifei) call cres

      ifield=jfield

      dtime=dnekclock()-stime
      rom_time=rom_time+dtime

      if (ifmult.and.nio.eq.0) write (6,*) 'romd_time: ',dtime
      if (.not.ifmult.or.nsteps.eq.istep) call final

    1 format(i8,1p4e13.5,' err')
    3 format(i8,1p4e15.7,' modt')

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

      call nekgsync
      setup_start=dnekclock()

      ttime=time

      n=lx1*ly1*lz1*nelt

      call rom_init_params
      call rom_init_fields

      call setgram
      call setevec

      call setbases
      call setops
      call setu

      call setf

      call setqoi
      call setmisc
      if (ifei) call set_sigma

      if (nio.eq.0) write (6,*) 'end range setup'

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF') call dump_misc

      time=ttime

      call nekgsync
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

      common /scrasnap/ t1(0:lub)

      call nekgsync
      asnap_time=dnekclock()

      if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
         if (ifrom(1)) then
            call read_serial(uas,nb+1,'ops/uas ',t1,nid)
            call read_serial(uvs,nb+1,'ops/uvs ',t1,nid)
         endif
         if (ifrom(2)) then
            call read_serial(tas,nb+1,'ops/tas ',t1,nid)
            call read_serial(tvs,nb+1,'ops/tvs ',t1,nid)
         endif
      else
         s=1./real(ns)

         if (ifrom(1)) then
            call pv2b(uas,uavg,vavg,wavg,ub,vb,wb)
            call rzero(uvs,nb+1)
            do j=1,ns
               do i=0,nb
                  uvs(i)=uvs(i)+(uk(i,j)-uas(i))**2
               enddo
            enddo
            call cmult(uvs,s,nb+1)

            call dump_serial(uas,nb+1,'ops/uas ',nid)
            call dump_serial(uvs,nb+1,'ops/uvs ',nid)
         endif

         if (ifrom(2)) then
            call ps2b(tas,tavg,tb)
            call rzero(tvs,nb+1)
            do j=1,ns
               do i=0,nb
                  tvs(i)=tvs(i)+(tk(i,j)-tas(i))**2
               enddo
            enddo
            call cmult(tvs,s,nb+1)

            call dump_serial(tas,nb+1,'ops/tas ',nid)
            call dump_serial(tvs,nb+1,'ops/tvs ',nid)
         endif
      endif

      call nekgsync
      if (nio.eq.0) write (6,*) 'asnap_time:',dnekclock()-asnap_time

      return
      end
c-----------------------------------------------------------------------
      subroutine setops

      include 'SIZE'
      include 'MOR'
      include 'TSTEP'

      logical iftmp

      if (nio.eq.0) write (6,*) 'inside setops'

      call nekgsync
      ops_time=dnekclock()

      jfield=ifield
      if (ifrom(1)) then
         ifield=1
         call seta(au,au0,'ops/au ')
         call setb(bu,bu0,'ops/bu ')
         call setc(cul,'ops/cu ')
      endif
      if (ifrom(2)) then
         ifield=2
         call seta(at,at0,'ops/at ')
         call setb(bt,bt0,'ops/bt ')
         call setc(ctl,'ops/ct ')
         call sets(st0,tb,'ops/ct ')
      endif

      if (ifbuoy.and.ifrom(1).and.ifrom(2)) call setbut(but0)

      ifield=jfield

      call set_binv2

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF') call dump_ops

      call nekgsync
      if (nio.eq.0) write (6,*) 'ops_time:',dnekclock()-ops_time

      if (nio.eq.0) write (6,*) 'exiting setops'

      return
      end
c-----------------------------------------------------------------------
      subroutine setqoi

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk2/ a1(lt),a2(lt),a3(lt),wk(lub+1)

      if (nio.eq.0) write (6,*) 'begin setup for qoi'

      call cvdrag_setup
      call cnuss_setup
      call cubar_setup

      if (ifcdrag) then
         if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
            call read_serial(fd1,ldim*(nb+1),'qoi/fd1 ',wk,nid)
            call read_serial(fd3,ldim*(nb+1),'qoi/fd3 ',wk,nid)
         else
            do i=0,nb
               call lap2d(a1,ub(1,i))
               call lap2d(a2,vb(1,i))
               if (ldim.eq.3) call lap2d(a3,wb(1,i))
               call cint(fd1(1+ldim*i),ub(1,i),vb(1,i),wb(1,i))
               call cint(fd3(1+ldim*i),a1,a2,a3)
            enddo
         endif
      endif

      if (nio.eq.0) write (6,*) 'end setup for qoi'

      return
      end
c-----------------------------------------------------------------------
      subroutine setmisc

      include 'SIZE'
      include 'MOR'

      if (nio.eq.0) write (6,*) 'begin range setup'

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF') then
         call nekgsync
         proj_time=dnekclock()

         if (ifpod(1)) call pv2k(uk,us0,ub,vb,wb)
         if (ifpod(2)) call ps2k(tk,ts0,tb)

         call nekgsync
         if (nio.eq.0) write (6,*) 'proj_time:',dnekclock()-proj_time
      endif

      call asnap
      call hyperpar

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_init_params

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      character*3 chartmp

      real a(1),b(1)

      if (nio.eq.0) write (6,*) 'inside rom_init_params'

      ad_nsteps=nsteps
      ad_iostep=iostep

      ubarr0=1e-1
      ubarrseq=5
      tbarr0=1e-1
      tbarrseq=5

      ustep_time = 0.
      solve_time=0.
      lu_time=0.
      ucopt_count=0
      copt_time=0.
      quasi_time=0.
      lnsrch_time=0.
      compgf_time=0.
      compf_time=0.

      anum_galu=0.
      anum_galt=0.

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

      np173=nint(param(173))

      if (np173.eq.0) then
         rmode='ALL'
      else if (np173.eq.1) then
         rmode='OFF'
      else if (np173.eq.2) then
         rmode='ON '
      else if (np173.eq.3) then
         rmode='ONB'
      else if (np173.eq.4) then
         rmode='CP '
      else
         call exitti('unsupported param(173), exiting...$',np173)
      endif

      if (rmode.eq.'ON '.or.rmode.eq.'ONB') then
         open (unit=10,file='ops/ips')
         read (10,*) chartmp
         close (unit=10)
         if (chartmp.ne.ips)
     $      call exitti('online ips does not match offline ips$',nb)
      endif

      ifrecon=(rmode.ne.'ON ')

      ifei=nint(param(175)).ne.0

      navg_step=nint(max(1.,param(176)))
      nb=nint(param(177))
      if (nb.eq.0) nb=lb

      rktol=param(179)
      if (rktol.lt.0.) rktol=10.**rktol

      ad_qstep=nint(param(180))+ad_iostep*max(1-nint(param(180)),0)

      iftneu=(param(178).ne.0.and.param(174).ne.0)

      ifctke=.false.
      if (param(181).ne.0) ifctke=.true.

      ifcdrag=.false.
      if (param(182).ne.0) ifcdrag=.true.

      inus=min(max(nint(param(183)),0),4)

      iffastc=.false.
      if (param(191).ne.0) iffastc=.true.

      iffasth=.false.
      if (param(192).ne.0.and.ips.eq.'HLM') iffasth=.true.

      ifcintp=.false.

      ifavisc=.false.
      if (param(196).ne.0.) ifavisc=.true.

      ifsub0=.true.
      if (param(197).ne.0.) ifsub0=.false.

      do i=0,ldimt1
         ifpod(i)=.false.
         ifrom(i)=.false.
      enddo
      ifpod(1)=param(174).ge.0.
      ifpod(2)=(ifheat.and.param(174).ne.0.)
c     ifrom(1)=(ifpod(1).and.eqn.ne.'ADE')
      ifrom(1)=ifpod(1)
      ifrom(2)=ifpod(2)

      ifpod(1)=ifpod(1).or.ifrom(2)

      ifvort=.false. ! default to false for now

      ifforce=param(193).ne.0.
      ifsource=param(194).ne.0.
      ifbuoy=param(195).ne.0.

      ifpart=.false.
      ifcintp=.false.


      ! constrained optimization parameter
      icopt=nint(param(184))
      barr_func=param(185)
      box_tol=param(186)

      ! barrier initial parameter and number of loop
      ubarr0  =param(187)
      ubarrseq=param(188)
      tbarr0  =param(189)
      tbarrseq=param(190)

      ! filter technique
      np198=nint(param(198))
      if (np198.eq.0) then
         rfilter='STD'
      else if (np198.eq.1) then
         rfilter='LER'
      else if (np198.eq.2) then
         rfilter='EF '
      else
         call exitti('unsupported param(198), exiting...$',np198)
      endif

      ! POD spatial filter
      rbf=nint(param(199))

      ! POD radius of the filter for differential filter
      rdft=param(200)

      call compute_BDF_coef(ad_alpha,ad_beta)

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF') then
         a(1)=nb*1.
         call dump_serial(a,1,'ops/nb ',nid)
      else
         call read_serial(a,1,'ops/nb ',b,nid)
         mb=a(1)
         if (mb.lt.nb) call exitti('mb less than nb, exiting...$',mb)
      endif

      if (nio.eq.0) then
         write (6,*) 'rp_nb         ',nb
         write (6,*) 'rp_lub        ',lub
         write (6,*) 'rp_ltb        ',ltb
         write (6,*) ' '
         write (6,*) 'rp_ls         ',ls
         write (6,*) 'rp_lsu        ',lsu
         write (6,*) 'rp_lst        ',lst
         write (6,*) ' '
         write (6,*) 'rp_isolve     ',isolve
         write (6,*) 'rp_ips        ',ips
         write (6,*) 'rp_ifavg0     ',ifavg0
         write (6,*) 'rp_ifsub0     ',ifsub0
         write (6,*) 'rp_rmode      ',rmode
         write (6,*) 'rp_ad_qstep   ',ad_qstep
         write (6,*) 'rp_ifctke     ',ifctke
         write (6,*) 'rp_ifcdrag    ',ifcdrag
         write (6,*) 'rp_inus       ',inus
         write (6,*) 'rp_iffastc    ',iffastc
         write (6,*) 'rp_iffasth    ',iffasth
         write (6,*) 'rp_ifavisc    ',ifavisc
         write (6,*) 'rp_ifei       ',ifei      
         write (6,*) ' '
         write (6,*) 'rp_ifforce    ',ifforce
         write (6,*) 'rp_ifsource   ',ifsource
         write (6,*) 'rp_ifbuoy     ',ifbuoy
         write (6,*) 'rp_ifpart     ',ifpart
         write (6,*) 'rp_ifvort     ',ifvort
         write (6,*) 'rp_ifcintp    ',ifcintp
         write (6,*) ' '
         do i=0,ldimt1
            write (6,*) 'rp_ifpod(',i,')   ',ifpod(i)
         enddo
         do i=0,ldimt1
            write (6,*) 'rp_ifrom(',i,')   ',ifrom(i)
         enddo
         write (6,*) ' '
         write (6,*) 'rp_icopt       ',icopt
         write (6,*) 'rp_barr_func   ',barr_func
         write (6,*) 'rp_box_tol     ',box_tol
         write (6,*) 'ubarr0         ',ubarr0
         write (6,*) 'ubarrseq       ',ubarrseq
         write (6,*) 'tbarr0         ',tbarr0
         write (6,*) 'tbarrseq       ',tbarrseq
         write (6,*) ' '
         write (6,*) 'rp_rfilter     ',rfilter
         write (6,*) 'rp_rbf         ',rbf
         write (6,*) 'rp_rdft        ',rdft
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

      logical alist,iftmp

      character*128 fname1

      if (nio.eq.0) write (6,*) 'inside rom_init_fields'

      n=lx1*ly1*lz1*nelt
      jfield=ifield
      ifield=1

      call rone(wm1,n)
      call rone(ones,n)
      call rzero(zeros,n)

      if (.not.iftneu) call rzero(gn,n)

      call rzero(num_galu,nb)
      call rzero(num_galt,nb)

      if (ifrom(1)) call opcopy(uic,vic,wic,vx,vy,vz)
      if (ifrom(2)) call copy(tic,t,n)

      ns = ls

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF') then
         fname1='file.list '
         nsu=1
         nsp=1
         nst=1
         if (ifrom(0)) nsp=lsp
         if (ifrom(1)) nsu=lsu
         if (ifrom(2)) nst=lst
         call get_saved_fields(us0,ps,ts0,nsu,nsp,nst,timek,fname1)

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

         iftmp=ifxyo
         ifxyo=.true.
         if (ifrecon) then
            call outpost(uavg,vavg,wavg,pavg,tavg,'avg')
            call outpost(urms,vrms,wrms,prms,trms,'rms')
            call outpost(vwms,wums,uvms,prms,trms,'tmn')
         endif
         ifxyo=iftmp

         if (ifsub0) then
            do i=1,ns
               if (ifrom(1)) then
                  call sub2(us0(1,1,i),ub,n)
                  call sub2(us0(1,2,i),vb,n)
                  if (ldim.eq.3) call sub2(us0(1,ldim,i),wb,n)
               endif
               if (ifrom(2)) call sub2(ts0(1,i),tb,n)
            enddo
            call sub2(uavg,ub,n)
            call sub2(vavg,vb,n)
            if (ldim.eq.3) call sub2(wavg,wb,n)
            if (ifpod(2)) call sub2(tavg,tb,n)
         endif
      endif

      if (ifbuoy) call set_ra

      ifield=jfield

      if (nio.eq.0) write (6,*) 'exiting rom_init_fields'

      return
      end
c-----------------------------------------------------------------------
      subroutine setc(cl,fname)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real cux(lt),cuy(lt),cuz(lt)

      common /scrcwk/ wk(lcloc),wk2(0:lub)

      real cl(lcloc),icl(3,lcloc)

      character*128 fname
      character*128 fnlint

      if (nio.eq.0) write (6,*) 'inside setc'

      call nekgsync
      conv_time=dnekclock()

      if (iffastc) call exitti('fastc not supported in setc_new$',nb)

      call cpart(kc1,kc2,jc1,jc2,ic1,ic2,ncloc,nb,np,nid+1) ! old indexing
c     call cpart(ic1,ic2,jc1,jc2,kc1,kc2,ncloc,nb,np,nid+1) ! new indexing

      n=lx1*ly1*lz1*nelv

      call lints(fnlint,fname,128)
      if (nid.eq.0) open (unit=100,file=fnlint)
      if (nio.eq.0) write (6,*) 'setc file:',fnlint

      if (rmode.eq.'ON '.or.rmode.eq.'ONB') then
         do k=0,nb
         do j=0,mb
         do i=1,mb
            cel=0.
            if (nid.eq.0) read(100,*) cel
            cel=glsum(cel,1)
            call setc_local(cl,cel,ic1,ic2,jc1,jc2,kc1,kc2,i,j,k)
         enddo
         enddo
         enddo
      else if (rmode.eq.'CP ') then
         ! read in the cp decomposition
         if (nid.eq.0) open (unit=100,file='./ops/lambda')
         do kk=1,ltr
            cel=0.
            read(100,*) cel 
            cp_w(kk) = cel
         enddo
         call rzero(cua,lub*ltr)
         if (nid.eq.0) open (unit=100,file='./ops/cua')
         do kk=1,ltr
         do i=1,nb
            cel=0.
            read(100,*) cel 
            cua(i+(kk-1)*lub) = cel
         enddo
         enddo
         call rzero(cub,(lub+1)*ltr)
         if (nid.eq.0) open (unit=100,file='./ops/cub')
         do kk=1,ltr
         do i=1,nb+1
            cel=0.
            read(100,*) cel 
            cub(i+(kk-1)*(lub+1)) = cel
         enddo
         enddo
         call rzero(cuc,(lub+1)*ltr)
         if (nid.eq.0) open (unit=100,file='./ops/cuc')
         do kk=1,ltr
         do i=1,nb+1
            cel=0.
            read(100,*) cel 
            cuc(i+(kk-1)*(lub+1)) = cel
         enddo
         enddo

         ! debug purpose
         ! forming the tensor
c        do kk=1,ltr
c        do k=1,nb+1
c        do j=1,nb+1
c        do i=1,nb
c           cl(i+(j-1)*(nb)+(k-1)*(nb+1)*(nb)) = 
c    $      cl(i+(j-1)*(nb)+(k-1)*(nb+1)*(nb)) + 
c    $      cp_w(kk)*cua(i+(kk-1)*lub)
c    $      * cub(j+(kk-1)*(lub+1))*cuc(k+(kk-1)*(lub+1))
c        enddo
c        enddo
c        enddo
c        enddo
      else
         if (.not.ifaxis) then
            do i=0,nb
               call set_convect_new(c1v(1,i),c2v(1,i),c3v(1,i),
     $                              ub(1,i),vb(1,i),wb(1,i))
               if (ifield.eq.1) then
                  call intp_rstd_all(u1v(1,i),ub(1,i),nelv)
                  call intp_rstd_all(u2v(1,i),vb(1,i),nelv)
                  if (ldim.eq.3) 
     $               call intp_rstd_all(u3v(1,i),wb(1,i),nelv)
               else
                  call intp_rstd_all(u1v(1,i),tb(1,i),nelv)
               endif
            enddo
         endif

         do k=0,nb
            if (nio.eq.0) write (6,*) 'setc: ',k,'/',nb
            do j=0,nb
               if (ifield.eq.1) then
                  call ccu(cux,cuy,cuz,k,j)
               else
                  call cct(cux,k,j)
               endif
               do i=1,nb
                  if (ifield.eq.1) then
                     cel=op_glsc2_wt(
     $                  ub(1,i),vb(1,i),wb(1,i),cux,cuy,cuz,ones)
                  else
                     cel=glsc2(tb(1,i),cux,n)
                  endif
                  call setc_local(cl,cel,ic1,ic2,jc1,jc2,kc1,kc2,i,j,k)
                  if (nid.eq.0) write (100,*) cel
               enddo
            enddo
         enddo
      endif

      if (nid.eq.0) close (unit=100)

      call nekgsync
      if (nio.eq.0) write (6,*) 'conv_time: ',dnekclock()-conv_time
      if (nio.eq.0) write (6,*) 'ncloc=',ncloc

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

      n=lx1*ly1*lz1*nelt

      if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
         if (nio.eq.0) write (6,*) 'reading a...'
         call read_mat_serial(a0,nb+1,nb+1,fname,mb+1,nb+1,wk1,nid)
      else
         if (nio.eq.0) write (6,*) 'forming a...'
         do j=0,nb
            if (nio.eq.0) write (6,*) 'seta: ',j,'/',nb
            nio=-1
            do i=0,nb
               if (ifield.eq.1) then
                  a0(i,j)=h10vip(ub(1,i),vb(1,i),wb(1,i),
     $                           ub(1,j),vb(1,j),wb(1,j))
               else
                  a0(i,j)=h10sip(tb(1,i),tb(1,j))
               endif
            enddo
            nio=nid
         enddo
      endif

      do j=1,nb
      do i=1,nb
         a(i,j)=a0(i,j)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine sets(s0,tt,fname)

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'GEOM'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrread/ tab((lub+1)**2)
      common /scrsets/ tx(lx1,ly1,lx1,lelt),
     $                 ty(lx1,ly1,lx1,lelt),
     $                 tz(lx1,ly1,lx1,lelt)

      real s0(0:nb),tt(lx1,ly1,lz1,lelt,0:nb)

      character*128 fname

      if (nio.eq.0) write (6,*) 'inside sets'

      if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
         if (nio.eq.0) write (6,*) 'reading s...'
         call read_mat_serial(s0,nb+1,nb+1,fname,mb+1,nb+1,tab,nid)
      else
         if (nio.eq.0) write (6,*) 'forming s...'
         nio=-1
         do i=0,nb
            if (nio.eq.0) write (6,*) 'sets: ',i,'/',nb
            s=0.
            do ie=1,nelt
            do ifc=1,2*ldim
               if (cbc(ifc,ie,2).ne.'E  '.and.
     $             cbc(ifc,ie,2).ne.'P  ') then
                  call facind(kx1,kx2,ky1,ky2,kz1,kz2,
     $                        lx1,ly1,lz1,ifc)
                  l=0
                  do iz=kz1,kz2
                  do iy=ky1,ky2
                  do ix=kx1,kx2
                     l=l+1
                     s=s+area(l,1,ifc,ie)*tt(ix,iy,iz,ie,i)
     $                                   *gn(ix,iy,iz,ie)
                  enddo
                  enddo
                  enddo
               endif
            enddo
            enddo
            s0(i)=glsum(s,1)
         enddo
         nio=nid
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setb(b,b0,fname)

      include 'SIZE'
      include 'TSTEP'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrread/ tab((lub+1)**2)

      real b(nb,nb),b0(0:nb,0:nb)

      character*128 fname

      if (nio.eq.0) write (6,*) 'inside setb'

      if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
         if (nio.eq.0) write (6,*) 'reading b...'
         call read_mat_serial(b0,nb+1,nb+1,fname,mb+1,nb+1,tab,nid)
      else
         if (nio.eq.0) write (6,*) 'forming b...'
         do j=0,nb
            if (nio.eq.0) write (6,*) 'setb: ',j,'/',nb
            nio=-1
            do i=0,nb
               if (ifield.eq.1) then
                  b0(i,j)=wl2vip(ub(1,i),vb(1,i),wb(1,i),
     $                           ub(1,j),vb(1,j),wb(1,j))
               else
                  b0(i,j)=wl2sip(tb(1,i),tb(1,j))
               endif
            enddo
            nio=nid
         enddo
      endif

      do j=1,nb
      do i=1,nb
         b(i,j)=b0(i,j)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setu

      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
      include 'INPUT'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrsetu/ uu(lt),vv(lt),ww(lt),tt(lt),wk(lb+1)

      logical iftmp,ifexist

      if (nio.eq.0) write (6,*) 'inside setu'

      n=lx1*ly1*lz1*nelv

      if (rmode.eq.'ON ') then
         inquire (file='ops/u0',exist=ifexist)
         if (ifexist) call read_serial(u,nb+1,'ops/u0 ',wk,nid)

         inquire (file='ops/t0',exist=ifexist)
         if (ifexist) call read_serial(ut,nb+1,'ops/t0 ',wk,nid)
      else
         jfield=ifield

         if (ifrom(1)) then
            ifield=1
            call opsub2(uic,vic,wic,ub,vb,wb)
            if (ips.eq.'H10') then
               call h10pv2b(u,uic,vic,wic,ub,vb,wb)
            else if (ips.eq.'HLM') then
               call hlmpv2b(u,uic,vic,wic,ub,vb,wb)
            else
               call pv2b(u,uic,vic,wic,ub,vb,wb)
            endif
            do i=0,nb
               if (nio.eq.0) write (6,*) 'u',u(i)
            enddo
            call opadd2(uic,vic,wic,ub,vb,wb)
         else
            call rzero(u,(nb+1)*3)
            u((nb+1)*0)=1.
            u((nb+1)*1)=1.
            u((nb+1)*2)=1.
         endif

         if (rmode.eq.'ALL'.or.rmode.eq.'OFF')
     $      call dump_serial(u,nb+1,'ops/u0 ',nid)

         if (ifrom(2)) then
            ifield=2
            call sub2(tic,tb,n)
            call ps2b(ut,tic,tb)
            do i=0,nb
               if (nio.eq.0) write (6,*) 'ut',ut(i)
            enddo
            call add2(tic,tb,n)
            if (rmode.eq.'ALL'.or.rmode.eq.'OFF')
     $         call dump_serial(ut,nb+1,'ops/t0 ',nid)
         endif
         ifield=jfield

         call reconv(uu,vv,ww,u)
         call recont(tt,ut)

         iftmp=ifxyo
         ifxyo=.true.
         if (ifrecon) call outpost(uu,vv,ww,pr,tt,'rom')

         ttime=time
         jstep=istep
         time=1.
         istep=1
         call outpost(uic,vic,wic,pr,tic,'uic')
         ifxyo=.false.
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

         ifxyo=iftmp
      endif

      if (nio.eq.0) write (6,*) 'exiting setu'

      return
      end
c-----------------------------------------------------------------------
      subroutine setf

      include 'SIZE'
      include 'SOLN'
      include 'MOR'
      include 'MASS'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrsetf/ wk1(lt),wk2(lt),wk3(lt)

      if (nio.eq.0) write (6,*) 'inside setf'

      call rzero(rf,nb)
      call rzero(rq,nb)
      call rzero(rg,nb)

      n=lx1*ly1*lz1*nelv

      if (ifforce.and.ifrom(1)) then ! assume fx,fy,fz has mass
         do i=1,nb
            rf(i)=glsc2(ub(1,i),fx,n)+glsc2(vb(1,i),fy,n)
            if (ldim.eq.3) rf(i)=rf(i)+glsc2(wb(1,i),fz,n)
            if (nio.eq.0) write (6,*) rf(i),i,'rf'
         enddo
         call opcopy(wk1,wk2,wk3,fx,fy,fz)
         call opbinv1(wk1,wk2,wk3,wk1,wk2,wk3,1.)
         call outpost(wk1,wk2,wk3,pavg,tavg,'fff')
      endif

      if (ifsource.and.ifrom(2)) then ! assume qq has mass
         do i=1,nb
            rq(i)=glsc2(qq,tb(1,i),n)
            if (nio.eq.0) write (6,*) rq(i),i,'rq'
         enddo
         call copy(wk1,qq,n)
         call binv1(wk1)
         call outpost(vx,vy,vz,pavg,wk1,'qqq')
      endif

      if (nio.eq.0) write (6,*) 'exiting setf'

      return
      end
c-----------------------------------------------------------------------
      subroutine final

      include 'SIZE'
      include 'MOR'

      common /romup/ rom_time

      real t1(0:nb),t2(0:nb)

      if (nio.eq.0) write (6,*) 'final...'
      if (nid.eq.0) then
         write (6,*) 'evalc_time:  ',evalc_time
         write (6,*) 'lu_time:     ',lu_time
         write (6,*) 'solve_time:  ',solve_time
         write (6,*) 'ustep_time:  ',ustep_time
         write (6,*) 'copt_time:   ',copt_time
         write (6,*) 'quasi_time:  ',quasi_time
         write (6,*) 'lnsrch_time: ',lnsrch_time
         write (6,*) 'compf_time:  ',compf_time
         write (6,*) 'compgf_time: ',compgf_time
         write (6,*) 'ucopt_active:',ucopt_count,
     $         '/',ad_step-1
         if (ifrom(2)) then
            write (6,*) 'tsolve_time: ',tsolve_time
            write (6,*) 'tstep_time:  ',tstep_time
            write (6,*) 'tcopt_active:',tcopt_count,
     $         '/',ad_step-1
         endif
         write (6,*) 'rom_time:    ',rom_time
         write (6,*) 'postu_time:  ',postu_time
         write (6,*) 'postt_time:  ',postt_time
      endif

      call dump_serial(u,nb+1,'ops/uf ',nid)
      if (ifrom(2)) then
         call dump_serial(ut,nb+1,'ops/tf ',nid)
      endif
      do i=0,nb
         t1(i)=u2a(i+i*(nb+1))-ua(i)*ua(i)
      enddo

      call dump_serial(t1,nb+1,'ops/uv ',nid)

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF') then
         call dump_serial(uas,nb+1,'ops/uas ',nid)
         call dump_serial(uvs,nb+1,'ops/uvs ',nid)
      endif

      call dump_serial(ua,nb+1,'ops/ua ',nid)
      call dump_serial(u2a,(nb+1)**2,'ops/u2a ',nid)
      if (ifrom(2)) call dump_serial(uta,nb+1,'ops/uta ',nid)

      if (ifrecon) call dump_sfld

      return
      end
c-----------------------------------------------------------------------
      subroutine set_ra

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      common /scrra/ binv(lx1,ly1,lz1,lelt)

      n=lx1*ly1*lz1*lelt

      call rone(binv,n)
      call invcol2(binv,bm1,n)

      ad_ra=sqrt(op_glsc2_wt(gx,gy,gz,gx,gy,gz,binv)/glsum(bm1,n))
      if (nio.eq.0) write (6,*) ad_ra,'ad_ra'
      s=1./ad_ra
      call opcmult(gx,gy,gz,s)

      return
      end
c-----------------------------------------------------------------------
      subroutine setbut(b)

      include 'SIZE'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrread/ tab((lb+1)**2)
      common /scruz/ wk1(lt),wk2(lt),wk3(lt)

      real b(0:nb,0:nb)

      character*128 fname

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF') then
         do j=0,nb
         do i=0,nb
            b(i,j)=op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),
     $                            gx,gy,gz,tb(1,j))
            if (nio.eq.0) write (6,*) i,j,b(i,j),'but0'
         enddo
         enddo
         call outpost(gx,gy,gz,pavg,tavg,'ggg')
         call dump_serial(b,(nb+1)**2,'ops/but ',nid)
      else
         fname='ops/but '
         if (nio.eq.0) write (6,*) 'reading but...'
         call read_mat_serial(but0,nb+1,nb+1,fname,mb+1,nb+1,tab,nid)
      endif

      return
      end
c-----------------------------------------------------------------------
