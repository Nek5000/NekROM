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

      call checker('baa',ad_step)

      ad_step = istep
      jfield=ifield
      ifield=1

      ifmult=.not.ifrom(2).and.ifheat


      if (rmode.ne.'OFF') then
      if (ifmult) then
         if (ifflow) call exitti(
     $   'error: running rom_update with ifflow = .true.$',nelv)
         if (istep.gt.0) then
            call bdfext_step
            call checker('bba',ad_step)
            call post
            call checker('bca',ad_step)
            call reconv(vx,vy,vz,u) ! reconstruct velocity to be used in h-t
            call checker('bda',ad_step)
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

         time=0
         do i=1,ad_nsteps
            if (cts.eq.'bdfext') then
               call bdfext_step
               time=time+ad_dt
               call post
            else if (cts.eq.'copt  ') then
               if (ifrom(2)) call rom_step_t_legacy
               if (ifrom(1)) call rom_step_legacy
               call postu_legacy
               call postt_legacy
               time=time+ad_dt
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

               time=time+ad_dt

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

            ad_step=ad_step+1
         enddo
         icalld=0
      endif
      endif

    2 continue

      if (ifei.and.rmode.ne.'OFF'.and.nint(param(175)).eq.1) then
         call cres_uns(sigma_u,sigma_t)
c        call cres
      endif

      if (nint(param(175)).eq.2) then
         call c_rieszrd_uns
      endif

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
      subroutine offline_mode ! offline-wrapper for MOR

      include 'SIZE'
      include 'INPUT'

      param(173)=1.
      call rom_update

      return
      end
c-----------------------------------------------------------------------
      subroutine online_mode ! online-wrapper for MOR

      include 'SIZE'
      include 'INPUT'

      param(173)=2.
      call rom_update

      return
      end
c-----------------------------------------------------------------------
      subroutine recon_mode ! reconstruction online-wrapper for MOR

      include 'SIZE'
      include 'INPUT'

      param(173)=3.
      call rom_update

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_setup

      include 'SIZE'
      include 'SOLN'
      include 'MOR'
      include 'AVG'
      include 'INPUT'

      logical iftmp,ifexist
      integer num_ts

      num_ts=1

      if (nio.eq.0) write (6,*) 'inside rom_setup'

      call nekgsync
      setup_start=dnekclock()

      ttime=time

      n=lx1*ly1*lz1*nelt

      call checker('aaa',ad_step)

      call set_bdf_coef(ad_alpha,ad_beta)
      call mor_init_params

      if (param(170).lt.0) then
         call mor_set_params_par
      else
         call mor_set_params_rea
      endif

      call mor_show_params

      call mor_init_fields
c     do k=2,10
c     call k_mean(k,nsu,nsp,nst,'sample.list ',k)
c     enddo
c     call exitt0

      call checker('aca',ad_step)

      inquire (file='ops/evec',exist=ifexist)
      if (ifexist) then
         open (unit=10,file='ops/evec')
         do i=1,ns
            read (10,*) (evec(i,j,1),j=1,nb)
         enddo
         close (unit=10)
      else
         call setgram
         call setevec
         call checker('ada',ad_step)

         call setbases

         call checker('aea',ad_step)
      endif

c     call average_in_xy
c     call average_in_y
      call setops
      call checker('afa',ad_step)
      call setu
      call checker('aga',ad_step)

      call setf
      call checker('aha',ad_step)

      call setqoi
      call checker('aia',ad_step)
      call setmisc
      call checker('aja',ad_step)

      if (ifei) then
c        call set_sigma
         call set_sigma_new
      endif

      if (ifplay) call set_trace
      call checker('ala',ad_step)

      if (nio.eq.0) write (6,*) 'end range setup'

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF') call dump_misc

      call checker('ama',ad_step)

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

      logical iftmp,ifexist
      real wk(lb+1)

      if (nio.eq.0) write (6,*) 'inside setops'

      call nekgsync
      ops_time=dnekclock()

      jfield=ifield
      if (ifrom(1)) then
         ifield=1
         call seta(au,au0,'ops/au ')
         call setb(bu,bu0,'ops/bu ')
         if (rmode.eq.'CP ') then 
            inquire (file='ops/u0',exist=ifexist)
            if (ifexist) call read_serial(u,nb+1,'ops/u0 ',wk,nid)
            call set_cp(cua,cub,cuc,cp_uw,cuj0,cu0k,cul,'ops/cu ',u)
         else
            call setc(cul,'ops/cu ')
         endif
      endif
      if (ifrom(2)) then
         ifield=2
         call seta(at,at0,'ops/at ')
         call setb(bt,bt0,'ops/bt ')
         if (rmode.eq.'CP ') then 
            inquire (file='ops/t0',exist=ifexist)
            if (ifexist) call read_serial(ut,nb+1,'ops/t0 ',wk,nid)
            call set_cp(cta,ctb,ctc,cp_tw,ctj0,ct0k,ctl,'ops/ct ',ut)
         else
            call setc(ctl,'ops/ct ')
         endif
         call sets(st0,tb,'ops/ct ')
      endif

c     if (ifbuoy.and.ifrom(1).and.ifrom(2)) call setbut(but0)
      if (ifbuoy.and.ifrom(1).and.ifrom(2)) then
         call setbut_xyz(but0_x,but0_y,but0_z)
         call setbut0
      endif


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
      subroutine update_k

      include 'SIZE'
      include 'MOR'

      if (nio.eq.0) write (6,*) 'begin update setup'

      call nekgsync
      proj_time=dnekclock()

      if (ifpod(1)) then
         do i=1,ns
            call mxm(wt,nb,uk(1,i),nb,ukp(1,i),1)
         enddo
      endif

      if (ifpod(2)) then
         do i=1,ns
            call mxm(wt(1,2),nb,tk(1,i),nb,tkp(1,i),1)
         enddo
      endif

      call nekgsync
      if (nio.eq.0) write (6,*) 'proj_time:',dnekclock()-proj_time

c     call hyperpar
      call update_hyper

      return
      end
c-----------------------------------------------------------------------
      subroutine update_hyper

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ep

      call nekgsync
      hpar_time=dnekclock()

      ! eps is the free parameter
      ! 1e-2 is used in the paper
      ep = 1.e-2

      n  = lx1*ly1*lz1*nelt
      if (ifpod(1)) then
         call cfill(upmin,1.e9,nb)
         call cfill(upmax,-1.e9,nb)
         do j=1,ns
         do i=1,nb
            if (ukp(i,j).lt.upmin(i)) upmin(i)=ukp(i,j)
            if (ukp(i,j).gt.upmax(i)) upmax(i)=ukp(i,j)
         enddo
         enddo
         do j=1,nb                    ! compute hyper-parameter
            d= upmax(j)-upmin(j)
            upmin(j) = upmin(j) - ep * d
            upmax(j) = upmax(j) + ep * d
            if (nio.eq.0) write (6,*) j,upmin(j),upmax(j),'upminmax'
         enddo

         ! compute distance between umax and umin
         call sub3(updis,upmax,upmin,nb)

         if (nio.eq.0) then
            do i=1,nb
               write (6,*) i,updis(i),'updis'
            enddo
         endif
      endif   

      if (ifpod(2)) then
         call cfill(tpmin,1.e9,nb)
         call cfill(tpmax,-1.e9,nb)
         do j=1,ns
         do i=1,nb
            if (tkp(i,j).lt.tpmin(i)) tpmin(i)=tkp(i,j)
            if (tkp(i,j).gt.tpmax(i)) tpmax(i)=tkp(i,j)
         enddo
         enddo
         do j=1,nb                    ! compute hyper-parameter
            d= tpmax(j)-tpmin(j)
            tpmin(j) = tpmin(j) - ep * d
            tpmax(j) = tpmax(j) + ep * d
            if (nio.eq.0) write (6,*) j,tpmin(j),tpmax(j),'tpminmax'
         enddo

         ! compute distance between tmax and tmin
         call sub3(tpdis,tpmax,tpmin,nb)
         if (nio.eq.0) then
            do i=1,nb
               write (6,*) i,tpdis(i),'tpdis'
            enddo
         endif

      endif   

      call nekgsync
      if (nio.eq.0) write (6,*) 'hpar_time',dnekclock()-hpar_time

      return
      end
c-----------------------------------------------------------------------
      subroutine setmisc

      include 'SIZE'
      include 'MOR'

      logical ifexist

      if (nio.eq.0) write (6,*) 'begin range setup'

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF') then
         call nekgsync
         proj_time=dnekclock()

         if (ifpod(1)) call pv2k(uk,us0,ub,vb,wb)
         if (ifpod(2)) call ps2k(tk,ts0,tb)

         call nekgsync
         if (nio.eq.0) write (6,*) 'proj_time:',dnekclock()-proj_time
      else if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP') then
         inquire (file='ops/uk',exist=ifexist)
         if (ifexist)
     $      call read_mat_serial(uk,lb+1,ns,'ops/uk ',mb+1,ns,stmp,nid)

         inquire (file='ops/tk',exist=ifexist)
         if (ifexist)
     $      call read_mat_serial(tk,lb+1,ns,'ops/tk ',mb+1,ns,stmp,nid)
      endif

      if (ifpod(1)) then
         do i=1,ns
            call copy(ukp(0,i),uk(0,i),nb+1)
         enddo
      endif

      if (ifpod(2)) then
         do i=1,ns
            call copy(tkp(0,i),tk(0,i),nb+1)
         enddo
      endif

      call asnap
      call hyperpar

      return
      end
c-----------------------------------------------------------------------
      subroutine mor_init_params

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      character*3 chartmp

      if (nio.eq.0) write (6,*) 'inside mor_init_params'

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
      invhm_time=0.

      anum_galu=0.
      anum_galt=0.

      ad_dt = dt
      ad_re = 1/param(2)
      ad_pe = 1/param(8)

      ifavg0=.false.

      ifctke=.false.
      ifcdrag=.false.
      iffastc=.false.
      iffasth=.false.
      ifcintp=.false.
      ifavisc=.false.
      ifdecpl=.false.
      ifsub0=.true.

      do i=0,ldimt1
         ifpod(i)=.false.
         ifrom(i)=.false.
      enddo

      ifvort=.false. ! default to false for now

      ifpart=.false.
      ifcintp=.false.

      ifcore=.true.
      ifquad=.true.

      ips='L2 '
      isolve=0
      rmode='OFF'
      ifei=.false.
      navg_step=1
      nb=lb
      rktol=0.
      ad_qstep=ad_iostep
      iftneu=.false.
      inus=0
      ifsub0=.true.

      ifforce=.false.
      ifsource=.false.
      ifbuoy=.false.
      ifplay=.false.
      nplay=0

      icopt=0

      rfilter='STD'
      rbf=0.5
      rdft=0.5

      if (nio.eq.0) write (6,*) 'exiting mor_init_params'

      return
      end
c-----------------------------------------------------------------------
      subroutine mor_set_params_rea

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      character*3 chartmp

      isolve=nint(param(170))

      if (param(171).eq.0) then
         ips='L2 '
      else if (param(171).eq.1) then
         ips='H10'
      else
         ips='HLM'
      endif

      if (param(172).ne.0) ifavg0=.true.

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
         max_tr = ltr
      else
         call exitti('unsupported param(173), exiting...$',np173)
      endif

      if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
         open (unit=10,file='ops/ips')
         read (10,*) chartmp
         close (unit=10)
         if (chartmp.ne.ips)
     $      call exitti('online ips does not match offline ips$',nb)
      endif

      ifrecon=(rmode.ne.'ON '.and.rmode.ne.'CP ')

      ifei=nint(param(175)).ne.0

      navg_step=nint(max(1.,param(176)))
      nb=min(nint(param(177)),lb)
      if (nb.eq.0) nb=lb

      if (ifavg0.and.(nb.eq.ls))
     $   call exitti('nb == ls results in linear dependent bases$',nb)

      if (nb.gt.ls)
     $   call exitti('nb > ls is undefined configuration$',nb)

      rktol=param(179)
      if (rktol.lt.0.) rktol=10.**rktol

      ad_qstep=nint(param(180))+ad_iostep*max(1-nint(param(180)),0)

      iftneu=(param(178).ne.0.and.param(174).ne.0)

      if (param(181).ne.0) ifctke=.true.

      if (param(182).ne.0) ifcdrag=.true.

      inus=min(max(nint(param(183)),0),5)

      if (param(191).ne.0) iffastc=.true.

      if (param(192).ne.0.and.ips.eq.'HLM') iffasth=.true.

      if (param(197).ne.0.) ifsub0=.false.

      ifpod(1)=param(174).ge.0.
      ifpod(2)=(ifheat.and.param(174).ne.0.)
c     ifrom(1)=(ifpod(1).and.eqn.ne.'ADE')
      ifrom(1)=ifpod(1)
      ifrom(2)=ifpod(2)

      ifpod(1)=ifpod(1).or.ifrom(2)

      ifforce=param(193).ne.0.
      ifsource=param(194).ne.0.
      ifbuoy=param(195).ne.0.

      nplay=nint(param(196))

      if (nplay.ne.0) then
         nplay=max(nplay,0)
         ifplay=.true.
      else
         ifplay=.false.
      endif

      ! constrained optimization parameter
      icopt=nint(param(184))
      barr_func=param(185)
      box_tol=param(186)

      ! barrier initial parameter and number of loop
      ubarr0  =param(187)
      ubarrseq=nint(param(188))
      tbarr0  =param(189)
      tbarrseq=nint(param(190))

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
      if (param(199).ge.0) then
         rbf=min(nint(param(199)),nb-1)
      else
         rbf=nint(nb*min(1.,-param(199)))
      endif

      ! POD radius of the filter for differential filter
      rdft=param(200)

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF') then
         rtmp1(1,1)=nb*1.
         call dump_serial(rtmp1(1,1),1,'ops/nb ',nid)
      else
         call read_serial(rtmp1(1,1),1,'ops/nb ',b,nid)
         mb=rtmp1(1,1)
         if (mb.lt.nb) call exitti('mb less than nb, exiting...$',mb)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mor_show_params

      include 'SIZE'
      include 'MOR'

      if (nio.eq.0) then
         write (6,*) 'mp_nb         ',nb
         write (6,*) 'mp_lub        ',lub
         write (6,*) 'mp_ltb        ',ltb
         write (6,*) ' '
         write (6,*) 'mp_ls         ',ls
         write (6,*) 'mp_lsu        ',lsu
         write (6,*) 'mp_lst        ',lst
         write (6,*) ' '
         write (6,*) 'mp_ad_re      ',ad_re
         write (6,*) 'mp_ad_pe      ',ad_pe
         write (6,*) ' '
         write (6,*) 'mp_isolve     ',isolve
         write (6,*) 'mp_ips        ',ips
         write (6,*) 'mp_ifavg0     ',ifavg0
         write (6,*) 'mp_ifsub0     ',ifsub0
         write (6,*) 'mp_rmode      ',rmode
         write (6,*) 'mp_ad_qstep   ',ad_qstep
         write (6,*) 'mp_ifctke     ',ifctke
         write (6,*) 'mp_ifcdrag    ',ifcdrag
         write (6,*) 'mp_inus       ',inus
         write (6,*) 'mp_iffastc    ',iffastc
         write (6,*) 'mp_iffasth    ',iffasth
         write (6,*) 'mp_nplay      ',nplay
         write (6,*) 'mp_ifplay     ',ifplay
         write (6,*) 'mp_ifei       ',ifei
         write (6,*) ' '
         write (6,*) 'mp_ifforce    ',ifforce
         write (6,*) 'mp_ifsource   ',ifsource
         write (6,*) 'mp_ifbuoy     ',ifbuoy
         write (6,*) 'mp_ifpart     ',ifpart
         write (6,*) 'mp_ifvort     ',ifvort
         write (6,*) 'mp_ifcintp    ',ifcintp
         write (6,*) 'mp_ifdecpl    ',ifdecpl
         write (6,*) ' '
         do i=0,ldimt1
            write (6,*) 'mp_ifpod(',i,')   ',ifpod(i)
         enddo
         do i=0,ldimt1
            write (6,*) 'mp_ifrom(',i,')   ',ifrom(i)
         enddo
         write (6,*) ' '
         write (6,*) 'mp_icopt       ',icopt
         write (6,*) 'mp_barr_func   ',barr_func
         write (6,*) 'mp_box_tol     ',box_tol
         write (6,*) 'mp_ubarr0      ',ubarr0
         write (6,*) 'mp_ubarrseq    ',ubarrseq
         write (6,*) 'mp_tbarr0      ',tbarr0
         write (6,*) 'mp_tbarrseq    ',tbarrseq
         write (6,*) ' '
         write (6,*) 'mp_rfilter     ',rfilter
         write (6,*) 'mp_rbf         ',rbf
         write (6,*) 'mp_rdft        ',rdft
         write (6,*) ' '
         write (6,*) 'mp_max_tr      ',max_tr
         write (6,*) 'mp_navg_step   ',navg_step
         write (6,*) 'mp_rk_tol      ',rk_tol
         write (6,*) 'mp_iftneu      ',iftneu
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mor_init_fields

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
      ! ns should be set in get_saved_fields

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
            ttime=time
            call auto_averager(fname1)
            time=ttime
            call copy_sol(uavg,vavg,wavg,pavg,tavg,vx,vy,vz,pr,t)
            call pop_sol(vx,vy,vz,pr,t)
         endif

         if (ifavg0) then
            idc_u=0
            idc_t=0
            do ie=1,nelt
            do ifc=1,2*ldim
               if (cbc(ifc,ie,1).eq.'W  '.or.
     $             cbc(ifc,ie,1).eq.'V  '.or.
     $             cbc(ifc,ie,1).eq.'v  ') idc_u=idc_u+1
               if (cbc(ifc,ie,2).eq.'T  '.or.
     $             cbc(ifc,ie,2).eq.'t  ') idc_t=idc_t+1
            enddo
            enddo
            idc_u=iglsum(idc_u,1)
            idc_t=iglsum(idc_t,1)
            call copy_sol(ub,vb,wb,pb,tb,uavg,vavg,wavg,pavg,tavg)
c           if (idc_u.gt.0) call opzero(ub,vb,wb)
c           if (idc_t.gt.0) call rzero(tb,n)
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

      real cl(lcloc)

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

      if (rmode.eq.'ON '.or.rmode.eq.'CP ') then
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
         if (ifsrct) then
            do i=1,nb
               rqt(i)=glsc2(qqt,tb(1,i),n)
               if (nio.eq.0) write (6,*) rqt(i),i,'rqt'
            enddo
            call copy(wk1,qqt,n)
            call binv1(wk1)
            call outpost(vx,vy,vz,pavg,wk1,'qqq')
         endif
      endif

      if (nio.eq.0) write (6,*) 'exiting setf'

      return
      end
c-----------------------------------------------------------------------
      subroutine final

      include 'SIZE'
      include 'MOR'

      common /romup/ rom_time,t1(0:lb) ,t2(0:lb)

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
         write (6,*) 'invhm_time:  ',invhm_time
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
         write (6,*) 'cp_time:     ',cp_time
      endif

      call dump_serial(u,nb+1,'ops/uf ',nid)
      if (ifrom(2)) then
         call dump_serial(ut,nb+1,'ops/tf ',nid)
      endif
      do i=0,nb
         t1(i)=u2a(1+i+i*(nb+1))-ua(i)*ua(i)
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

      n=lx1*ly1*lz1*nelt

      call rone(binv,n)
      call invcol2(binv,bm1,n)

c     ad_ra=sqrt(op_glsc2_wt(gx,gy,gz,gx,gy,gz,binv)/glsum(bm1,n))
      ad_ra=sqrt(op_glsc2_wt(gx,gy,gz,gx,gy,gz,binv)
     $          /op_glsc2_wt(ones,ones,ones,ones,ones,ones,bm1))

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
      include 'INPUT'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrread/ tab((lb+1)**2)
      common /scrk/ wk1(lt),wk2(lt),wk3(lt)

      logical iftmp

      real b(0:nb,0:nb)

      character*128 fname

      n=lx1*ly1*lz1*nelv

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF') then
         do j=0,nb
         do i=0,nb
            b(i,j)=op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),
     $                            gx,gy,gz,tb(1,j))
            if (nio.eq.0) write (6,*) i,j,b(i,j),'but0'
         enddo
         enddo

         iftmp=ifxyo
         ifxyo=.true.
         call outpost(gx,gy,gz,pavg,tavg,'ggg')

         call invcol3(wk1,gx,bm1,n)
         call invcol3(wk2,gy,bm1,n)
         call invcol3(wk3,gz,bm1,n)

         call outpost(wk1,wk2,wk3,pavg,tavg,'ggg')
         ifxyo=iftmp
         call dump_serial(b,(nb+1)**2,'ops/but ',nid)
      else
         fname='ops/but '
         if (nio.eq.0) write (6,*) 'reading but...'
         call read_mat_serial(but0,nb+1,nb+1,fname,mb+1,nb+1,tab,nid)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setbut_xyz(b_x,b_y,b_z)

      include 'SIZE'
      include 'MOR'
      include 'AVG'
      include 'INPUT'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrread/ tab((lb+1)**2)
      common /scrk/ wk1(lt),wk2(lt),wk3(lt)

      logical iftmp

      real b_x(0:nb,0:nb)
      real b_y(0:nb,0:nb)
      real b_z(0:nb,0:nb)

      character*128 fname

      n=lx1*ly1*lz1*nelv

c     call cmult(gx,1./sqrt(2.),n)
c     call cmult(gy,1./sqrt(2.),n)

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF') then
         do j=0,nb
         do i=0,nb
            b_x(i,j)=glsc3(ub(1,i),gx,tb(1,j),n)
            b_y(i,j)=glsc3(vb(1,i),gy,tb(1,j),n)
            if (ldim.eq.3) then
               b_z(i,j)=glsc3(wb(1,i),gz,tb(1,j),n)
            endif
            if (nio.eq.0) write (6,*) i,j,b_x(i,j),'but0'
         enddo
         enddo

         iftmp=ifxyo
         ifxyo=.true.
         call outpost(gx,gy,gz,pavg,tavg,'ggg')

         call invcol3(wk1,gx,bm1,n)
         call invcol3(wk2,gy,bm1,n)
         call invcol3(wk3,gz,bm1,n)

         call outpost(wk1,wk2,wk3,pavg,tavg,'ggg')
         ifxyo=iftmp
         call dump_serial(b_x,(nb+1)**2,'ops/but_x ',nid)
         call dump_serial(b_y,(nb+1)**2,'ops/but_y ',nid)
         if (ldim.eq.3) then
            call dump_serial(b_z,(nb+1)**2,'ops/but_z ',nid)
         endif
      else
         fname='ops/but_x '
         if (nio.eq.0) write (6,*) 'reading but_x...'
         call read_mat_serial(but0_x,nb+1,nb+1,fname,mb+1,nb+1,tab,nid)
         fname='ops/but_y '
         if (nio.eq.0) write (6,*) 'reading but_y...'
         call read_mat_serial(but0_y,nb+1,nb+1,fname,mb+1,nb+1,tab,nid)
         if (ldim.eq.3) then
            fname='ops/but_z '
            if (nio.eq.0) write (6,*) 'reading but_z...'
            call read_mat_serial(but0_z,nb+1,nb+1,fname,mb+1,nb+1,
     $                           tab,nid)
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setbut0

      include 'SIZE'
      include 'MOR'
      include 'AVG'
      include 'INPUT'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrread/ tab((lb+1)**2)
      common /scrk/ wk1(lt),wk2(lt),wk3(lt)

      character*128 fname

      n=lx1*ly1*lz1*nelv

      call add3s2(but0,but0_x,but0_y,sin(bu_angle),
     $            cos(bu_angle),(nb+1)**2)

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF') then
         call dump_serial(but0,(nb+1)**2,'ops/but ',nid)
      endif
      return
      end
c-----------------------------------------------------------------------
