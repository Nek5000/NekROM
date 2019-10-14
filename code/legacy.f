c-----------------------------------------------------------------------
      subroutine rom_step_legacy

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real rhs(0:lb),rhstmp(0:lb)
      logical ifdebug
      integer chekbc

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

      call shift3(u,rhs,nb+1)

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
         call shift3(ut,rhs,nb+1)
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

      call shift3(ut,rhs,nb+1)

      tstep_time=tstep_time+dnekclock()-tlast_time

      return
      end
c-----------------------------------------------------------------------
