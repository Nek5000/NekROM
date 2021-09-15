c-----------------------------------------------------------------------
      subroutine rom_step_legacy

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      character*3 rfilter

      real rhs(0:lb),rhstmp(0:lb)
      logical ifdebug
      integer chekbc

      rfilter='   '
      if (cfloc.eq.'POST'.and.cftype.eq.'DIFF') rfilter='EF '

      chekbc=0

      ifdebug=.true.
      ifdebug=.false.

      ulast_time = dnekclock()

      n=lx1*ly1*lz1*nelt

      icount = min(max(1,ad_step),3)

      rhs(0)=1.
      call setr_v(rhs(1),icount)

      do i=0,nb
         if (ifdebug) write (6,*) i,u(i),'sol'
      enddo

      do i=0,nb
         if (ifdebug) write (6,*) i,rhs(i),'rhs'
      enddo

      ttime=dnekclock()
      call seth(fluv,au,bu,1./ad_re)
      if (ad_step.eq.3) call dump_serial(fluv,nb*nb,'ops/hu ',nid)
      if (ad_step.le.3) then
         call copy(helmu,fluv,nb*nb)
         call dgetrf(nb,nb,fluv,nb,ipiv,info)
         call copy(invhelmu,fluv,nb*nb)
      endif
      lu_time=lu_time+dnekclock()-ttime

      do j=1,nb
      do i=1,nb
         if (ifdebug) write (6,*) i,j,au(i+(j-1)*nb),'au'
      enddo
      enddo

      do j=1,nb
      do i=1,nb
         if (ifdebug) write (6,*) i,j,bu(i+(j-1)*nb),'bu'
      enddo
      enddo

      do j=1,nb
      do i=1,nb
         if (ifdebug) write (6,*) i,j,fluv(i+(j-1)*nb),'LU'
      enddo
      enddo

      ttime=dnekclock()
      if ((isolve.eq.0).or.(icopt.eq.2)) then ! standard matrix inversion
         if (.not.iffasth.or.ad_step.le.3) then
            call dgetrs('N',nb,1,fluv,nb,ipiv,rhs(1),nb,info)
         else
            eps=.20
            damp=1.-eps*ad_dt
            do i=1,nb
            if (rhs(i).gt.umax(i)) rhs(i)=umax(i)+(rhs(i)-umax(i))*damp
            if (rhs(i).lt.umin(i)) rhs(i)=umin(i)+(rhs(i)-umin(i))*damp
            enddo
         endif
      else 
         call constrained_POD(rhs,u(1),helmu,invhelmu,umax,umin,udis,
     $                        ubarr0,ubarrseq,ucopt_count)
      endif
      solve_time=solve_time+dnekclock()-ttime

      do i=0,nb
         if (ifdebug) write (6,*) i,rhs(i),'sol'
      enddo

      if (ifdebug) call exitt0

      if (rfilter.eq.'EF ') then
         if (rbf.gt.0) then
            call pod_proj(rhs(1),rbf)
         else if (rbf.lt.0) then
            call pod_df(rhs(1))
         endif
      endif

      call count_gal(num_galu,anum_galu,rhs(1),umax,umin,1e-16,nb)

      call shift(u,rhs,nb+1,3)

      ustep_time=ustep_time+dnekclock()-ulast_time

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_step_t_legacy

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      real rhs(0:lb),rhstmp(0:lb)
      logical ifdebug
      integer chekbc

      common /scrrstep/ t1(lt),t2(lt),t3(lt),work(lt)

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      chekbc = 0

      ifdebug=.true.
      ifdebug=.false.

      if (ad_step.eq.1) then
         tstep_time = 0.
         copt_time=0.
         quasi_time=0.
         lnsrch_time=0.
         compgf_time=0.
         compf_time=0.
         tcopt_count = 0
      endif

      if (nb.eq.0) then
         rhs(0)=1.
         call shift(ut,rhs,nb+1,3)
         return
      endif

      tlast_time = dnekclock()

      n=lx1*ly1*lz1*nelt

      icount = min0(ad_step,3)

      rhs(0)=1.
      call setr_t(rhs(1),icount)

      do i=0,nb
         if (ifdebug) write (6,*) i,ut(i),'sol'
      enddo

      do i=0,nb
         if (ifdebug) write (6,*) i,rhs(i),'rhs'
      enddo

      call seth(flut,at,bt,1./ad_pe)
      if (ad_step.eq.3) call dump_serial(flut,nb*nb,'ops/ht ',nid)
      if (ad_step.le.3) then
         call copy(helmt,flut,nb*nb)
         call dgetrf(nb,nb,flut,nb,ipiv,info)
         call copy(invhelmt,flut,nb*nb)
      endif

      do j=1,nb
      do i=1,nb
         if (ifdebug) write (6,*) i,j,at(i+(j-1)*nb),'at'
      enddo
      enddo

      do j=1,nb
      do i=1,nb
         if (ifdebug) write (6,*) i,j,bt(i+(j-1)*nb),'bt'
      enddo
      enddo

      do j=1,nb
      do i=1,nb
         if (ifdebug) write (6,*) i,j,flut(i+(j-1)*nb),'LU'
      enddo
      enddo

      ttime=dnekclock()
      if ((isolve.eq.0).or.(icopt.eq.1)) then ! standard matrix inversion
         call dgetrs('N',nb,1,flut,nb,ipiv,rhs(1),nb,info)
      else 
         call constrained_POD(rhs,ut(1),helmt,invhelmt,tmax,tmin,tdis,
     $                        tbarr0,tbarrseq,tcopt_count) 
      endif

      tsolve_time=tsolve_time+dnekclock()-ttime

      do i=0,nb
         if (ifdebug) write (6,*) i,rhs(i),'sol'
      enddo

      call count_gal(num_galt,anum_galt,rhs(1),tmax,tmin,1e-16,nb)

      call shift(ut,rhs,nb+1,3)

      tstep_time=tstep_time+dnekclock()-tlast_time

      return
      end
c-----------------------------------------------------------------------
      subroutine postu_legacy

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrrstep/ t1(lt),t2(lt),t3(lt),work(lt)
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      save icalld
      data icalld /0/

      real vort(lt)

      if (icalld.eq.0) then
         postu_time=0.
         icalld=1
      endif

      call nekgsync
      tttime=dnekclock()

      call setuavg(ua,u2a,u)
      call setuj(uj,u2j,u)

      if (mod(ad_step,ad_qstep).eq.0) then
         if (ifctke) call ctke
         if (ifcdrag) call cdrag
         call cnuss
c        call cubar
      endif

      if (mod(ad_step,ad_iostep).eq.0) then
         if (nio.eq.0) then
            if (ifrom(1)) then
               do j=1,nb
                  write(6,*) j,time,u(j),'romu'
               enddo
            endif
         endif

         if (rmode.ne.'ON ') then
            idump=ad_step/ad_iostep
            call reconv(vx,vy,vz,u)
            call opcopy(t1,t2,t3,vx,vy,vz)

            if (ifrom(2)) then
               call recont(vort,ut)
            else
               call comp_vort3(vort,work1,work2,t1,t2,t3)
            endif

            ifto = .true. ! turn on temp in fld file
            call outpost(vx,vy,vz,pavg,vort,'rom')
         endif
      endif

      if (ad_step.eq.ad_nsteps) then
         if (nio.eq.0) then
            do j=1,nb
               write(6,*)j,num_galu(j)/ad_nsteps,'num_galu'
            enddo
            write(6,*)anum_galu/ad_nsteps,'anum_galu'
         endif
      endif

      call nekgsync
      postu_time=postu_time+dnekclock()-tttime

      return
      end
c-----------------------------------------------------------------------
      subroutine postt_legacy

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrrstep/ t1(lt),t2(lt),t3(lt),work(lt)
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      save icalld
      data icalld /0/

      if (icalld.eq.0) then
         postt_time=0.
         icalld=1
      endif

      call nekgsync
      tttime=dnekclock()

      call settavg(uta,uuta,utua,ut2a,u,ut)
      call settj(utj,uutj,utuj,uj,ut)

      if (mod(ad_step,ad_iostep).eq.0) then
         if (ifrom(2)) then
            if (nio.eq.0) then
               do j=1,nb
                  write(6,*) j,time,ut(j),'romt'
               enddo
            endif
         endif
      endif

      if (ad_step.eq.ad_nsteps) then
         if (nio.eq.0) then
            do j=1,nb
               write(6,*)j,num_galt(j)/ad_nsteps,'num_galt'
            enddo
            write(6,*)anum_galt/ad_nsteps,'anum_galt'
         endif
      endif

      call nekgsync
      postt_time=postt_time+dnekclock()-tttime

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_init_params

      include 'SIZE'
      include 'INPUT'

      if (nid.eq.0) write (6,*) 'rom_init_params has been deprecated...'

      call mor_init_params

      if (param(170).ge.0.) then 
         call mor_set_params_rea
      else
         call mor_set_params_par
      endif

      call mor_show_params

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_init_fields

      include 'SIZE'

      if (nid.eq.0) write (6,*) 'rom_init_fields has been deprecated...'
      call mor_init_fields

      return
      end
c-----------------------------------------------------------------------
      subroutine interp_t(scalar,xyz,n,s)
      real scalar(1),xyz(1),s(1)
      call interp_s(scalar,s,xyz,n)
      return
      end
c-----------------------------------------------------------------------
