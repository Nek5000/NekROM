c-----------------------------------------------------------------------
      subroutine rom_step

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real rhs(0:nb),rhstmp(0:nb)
      real bctol
      logical ifdebug
      integer chekbc

      chekbc=0

      ifdebug=.true.
      ifdebug=.false.

      if (ad_step.eq.1) then
         ustep_time = 0.
         solve_time=0.
         lu_time=0.
         ucopt_count=0
         if (.not.ifrom(2)) then
            copt_time=0.
            quasi_time=0.
            lnsrch_time=0.
         endif
      endif

      ulast_time = dnekclock()

      n=lx1*ly1*lz1*nelt

      if (nb.eq.0) then
         rhs(0)=1.
         call shift3(u,rhs,nb+1)
         return
      endif

      icount = min0(max(1,ad_step),3)

      rhs(0)=1.
      call setr_v(rhs(1),icount)

      do i=0,nb
         if (ifdebug) write (6,*) i,u(i,1),'sol'
      enddo

      do i=0,nb
         if (ifdebug) write (6,*) i,rhs(i),'rhs'
      enddo

      ttime=dnekclock()
      call seth(fluv,au,bu,1./ad_re)
      if (ad_step.eq.3) call dump_serial(fluv,nb*nb,'ops/hu ',nid)
      if (ad_step.le.3) then
         call copy(helmu,fluv,nb*nb)
         call dgetrf(nb,nb,fluv,lub,ipiv,info)
         call copy(invhelmu,fluv,nb*nb)
      endif
      lu_time=lu_time+dnekclock()-ttime

      do j=1,nb
      do i=1,nb
         if (ifdebug) write (6,*) i,j,au(i,j),'au'
      enddo
      enddo

      do j=1,nb
      do i=1,nb
         if (ifdebug) write (6,*) i,j,bu(i,j),'bu'
      enddo
      enddo

      do j=1,nb
      do i=1,nb
         if (ifdebug) write (6,*) i,j,fluv(i,j),'LU'
      enddo
      enddo

      ttime=dnekclock()
      if (isolve.eq.0) then ! standard matrix inversion
         if (.not.iffasth.or.ad_step.le.3) then
            call dgetrs('N',nb,1,fluv,lub,ipiv,rhs(1),nb,info)
         else
            eps=.20
            damp=1.-eps*ad_dt
            do i=1,nb
            if (rhs(i).gt.umax(i)) rhs(i)=umax(i)+(rhs(i)-umax(i))*damp
            if (rhs(i).lt.umin(i)) rhs(i)=umin(i)+(rhs(i)-umin(i))*damp
            enddo
         endif
      else if (isolve.eq.1) then ! constrained solve with inverse update
         call BFGS_new(rhs(1),u(1,1),helmu,invhelmu,umax,umin,udis,
     $   ubarr0,ubarrseq)
      else if (isolve.eq.2) then ! constrained solve with inverse update
                                 ! and mix with standard solver
         call copy(rhstmp,rhs,nb+1)
         call dgetrs('N',nb,1,fluv,lub,ipiv,rhstmp(1),nb,info)

         bctol = 1e-12
         do ii=1,nb
            if ((rhstmp(ii)-umax(ii)).ge.bctol) then
               chekbc = 1
            elseif ((umin(ii)-rhstmp(ii)).ge.bctol) then
               chekbc = 1
            endif
         enddo

         if (chekbc.eq.1) then
            ucopt_count = ucopt_count + 1
            call BFGS_new(rhs(1),u(1,1),helmu,invhelmu,umax,umin,udis,
     $      ubarr0,ubarrseq)
         else
            call copy(rhs,rhstmp,nb+1)
         endif

      else if (isolve.eq.3) then ! constrained solve with Hessian update
         call BFGS(rhs(1),u(1,1),helmu,invhelmu,umax,umin,udis,
     $   ubarr0,ubarrseq)
      else if (isolve.eq.4) then ! constrained solve with Hessian update
                                 ! and mix with standard solver

         call hybrid_advance(rhs,u(1,1),helmu,invhelmu,umax,umin,
     $                       udis,ubarr0,ubarrseq,ucopt_count)
      else   
         call exitti('incorrect isolve specified...$',isolve)
      endif
      solve_time=solve_time+dnekclock()-ttime

      do i=0,nb
         if (ifdebug) write (6,*) i,rhs(i),'sol'
      enddo

      if (ifdebug) call exitt0

      call count_gal(num_galu,anum_galu,rhs(1),umax,umin,1e-12,nb)

      call shift3(u,rhs,nb+1)

      ustep_time=ustep_time+dnekclock()-ulast_time

      return
      end
c-----------------------------------------------------------------------
      subroutine postu

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

      call setuavg
      call setuj

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
                  write(6,*) j,time,u(j,1),'romu'
               enddo
            endif
         endif

         if (ifdump) then
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
      subroutine rom_step_t

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      real rhs(0:nb),rhstmp(0:nb)
      real bctol
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
         if (ifdebug) write (6,*) i,ut(i,1),'sol'
      enddo

      do i=0,nb
         if (ifdebug) write (6,*) i,rhs(i),'rhs'
      enddo

      call seth(flut,at,bt,1./ad_pe)
      if (ad_step.eq.3) call dump_serial(flut,nb*nb,'ops/ht ',nid)
      if (ad_step.le.3) then
         call copy(helmt,flut,nb*nb)
         call dgetrf(nb,nb,flut,lub,ipiv,info)
         call copy(invhelmt,flut,nb*nb)
      endif

      do j=1,nb
      do i=1,nb
         if (ifdebug) write (6,*) i,j,at(i,j),'at'
      enddo
      enddo

      do j=1,nb
      do i=1,nb
         if (ifdebug) write (6,*) i,j,bt(i,j),'bt'
      enddo
      enddo

      do j=1,nb
      do i=1,nb
         if (ifdebug) write (6,*) i,j,flut(i,j),'LU'
      enddo
      enddo

      ttime=dnekclock()
      if (isolve.eq.0) then ! standard matrix inversion
         call dgetrs('N',nb,1,flut,lub,ipiv,rhs(1),nb,info)
      else if (isolve.eq.1) then ! constrained solve
         call BFGS_new(rhs(1),ut(1,1),helmt,invhelmt,tmax,tmin,tdis,
     $   tbarr0,tbarrseq) 
      else if (isolve.eq.2) then

         call copy(rhstmp,rhs,nb+1)
         call dgetrs('N',nb,1,flut,lub,ipiv,rhstmp(1),nb,info)

         bctol = 1e-12
         do ii=1,nb
            if ((rhstmp(ii)-tmax(ii)).ge.bctol) then
               chekbc = 1
            elseif ((tmin(ii)-rhstmp(ii)).ge.bctol) then
               chekbc = 1
            endif
         enddo

         if (chekbc.eq.1) then
            tcopt_count = tcopt_count + 1
            call BFGS_new(rhs(1),ut(1,1),helmt,invhelmt,tmax,tmin,tdis,
     $      tbarr0,tbarrseq) 
         else
            call copy(rhs,rhstmp,nb+1)
         endif

      else if (isolve.eq.3) then ! constrained solve
         call BFGS(rhs(1),ut(1,1),helmt,invhelmt,tmax,tmin,tdis,
     $   tbarr0,tbarrseq)
      else if (isolve.eq.4) then

         call hybrid_advance(rhs,ut(1,1),helmt,invhelmt,tmax,tmin,
     $                       tdis,tbarr0,tbarrseq,tcopt_count)
      else
         call exitti('incorrect isolve specified...$',isolve)
      endif
      tsolve_time=tsolve_time+dnekclock()-ttime

      do i=0,nb
         if (ifdebug) write (6,*) i,rhs(i),'sol'
      enddo

      call count_gal(num_galt,anum_galt,rhs(1),tmax,tmin,1e-12,nb)

      call shift3(ut,rhs,nb+1)

      tstep_time=tstep_time+dnekclock()-tlast_time

      return
      end
c-----------------------------------------------------------------------
      subroutine postt

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

      call settavg
      call settj

      if (mod(ad_step,ad_iostep).eq.0) then
         if (ifrom(2)) then
            if (nio.eq.0) then
               do j=1,nb
                  write(6,*) j,time,ut(j,1),'romt'
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
      subroutine compute_BDF_coef(ad_alpha,ad_beta)

      real ad_alpha(3,3), ad_beta(4,3)

      call rzero(ad_alpha,3*3)
      call rzero(ad_beta,3*4)

      ad_beta(1,1) = 1.
      ad_beta(2,1) = -1.

      ad_beta(1,2) = 1.5
      ad_beta(2,2) = -2
      ad_beta(3,2) = 0.5

      ad_beta(1,3) = 11./6
      ad_beta(2,3) = -3
      ad_beta(3,3) = 1.5
      ad_beta(4,3) = -1./3.

      ad_alpha(1,1)=1

      ad_alpha(1,2)=2
      ad_alpha(2,2)=-1

      ad_alpha(1,3)=3
      ad_alpha(2,3)=-3
      ad_alpha(3,3)=1

      return
      end
c-----------------------------------------------------------------------
      subroutine evalc(cu,cl,uu,n)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real cu(n)
      real uu(0:n)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)

      common /scrc/ work(max(lub,ltb))

      integer icalld
      save    icalld
      data    icalld /0/

      if (icalld.eq.0) then
         evalc_time=0.
         icalld=1
      endif

      stime=dnekclock()

      if (ifcintp) then
         call mxm(cintp,n,uu,n+1,cu,1)
      else
         if (ncloc.ne.0) then
            l=1

            call rzero(cu,n)

            do k=kc1,kc2
            do j=jc1,jc2
            do i=ic1,ic2
               cu(i)=cu(i)+cl(i,j,k)*uu(j)*u(k,1)
            enddo
            enddo
            enddo
         endif
         call gop(cu,work,'+  ',n)
      endif

      if (nio.eq.0) then
         do i=1,n
            write (6,*) i,cu(i)
         enddo
      endif

      call nekgsync
      call exitt0

      evalc_time=evalc_time+dnekclock()-stime

      return
      end
c-----------------------------------------------------------------------
      subroutine evalc_legacy(cu,cl,icl,uu)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real cu(nb)

      integer icalld
      save    icalld
      data    icalld /0/

      common /scrc/ work(lx1*ly1*lz1*lelt)

      real cl(lub,0:lub,0:lub)
      real uu(0:nb)
      integer icl(3,lcloc)

      if (icalld.eq.0) then
         evalc_time=0.
         icalld=1
      endif

      stime=dnekclock()

      if (ifcintp) then
         call mxm(cintp,nb,uu,nb+1,cu,1)
      else
         l=1

         call rzero(cu,nb)

         if (np.eq.1) then ! don't use index
            do i=1,nb
            do k=0,nb
            do j=0,nb
               cu(i)=cu(i)+cl(i,j,k)*uu(j)*u(k,1)
            enddo
            enddo
            enddo
         else
            do l=1,ncloc
               i=icl(1,l)
               j=icl(2,l)
               k=icl(3,l)
               cu(i)=cu(i)+cl(l,0,0)*uu(j)*u(k,1)
            enddo
            call gop(cu,work,'+  ',nb)
         endif
      endif

      evalc_time=evalc_time+dnekclock()-stime

      return
      end
c-----------------------------------------------------------------------
      subroutine setcintp
      call exitti('called deprecated subroutine setcintp$',1)
      return
      end
c-----------------------------------------------------------------------
      subroutine setr_t(rhs,icount)

      include 'SIZE'
      include 'MOR'

      common /scrrhs/ tmp(0:nb)

      real rhs(nb)

      call mxm(ut,nb+1,ad_beta(2,icount),3,tmp,1)
c     call mxm(bv0,nb+1,tmp,nb+1,rhs,1)
      call mxm(bt,nb,tmp(1),nb,rhs,1)

      call cmult(rhs,-1.0/ad_dt,nb)

      s=-1.0/ad_pe

c     call add2s2(rhs,av0,s,nb+1) ! not working...
      do i=1,nb
         rhs(i)=rhs(i)+s*at0(i,0)
      enddo

      call evalc(tmp(1),ctl,ut,nb)

      call shift3(ctr,tmp(1),nb)

      call mxm(ctr,nb,ad_alpha(1,icount),3,tmp(1),1)

      call sub2(rhs,tmp(1),nb)

      if (ifsource) then
         call add2(rhs,rq,nb)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setr_v(rhs,icount)

      include 'SIZE'
      include 'MOR'

      common /scrrhs/ tmp1(0:nb),tmp2(0:nb)

      real rhs(nb)

      call mxm(u,nb+1,ad_beta(2,icount),3,tmp1,1)
c     call mxm(bv0,nb+1,tmp,nb+1,rhs,1)
      call mxm(bu,nb,tmp1(1),nb,rhs,1)

      call cmult(rhs,-1.0/ad_dt,nb)

      s=-1.0/ad_re

c     call add2s2(rhs,av0,s,nb+1) ! not working...
      do i=1,nb
         rhs(i)=rhs(i)+s*au0(i,0)
      enddo

      call evalc(tmp1(1),cul,u,nb)
c     call evalc_legacy(tmp1(1),cul,icul,u,nb)
      call chsign(tmp1(1),nb)

      if (ifbuoy) then
         call mxm(but0,nb+1,ut(0,1),nb+1,tmp2(0),1)
         call add2s2(tmp1(1),tmp2(1),ad_ra,nb)
      else if (ifforce) then
         call add2(tmp1(1),rg(1),nb)
      endif

      call shift3(fu,tmp1(1),nb)

      call mxm(fu,nb,ad_alpha(1,icount),3,tmp1(1),1)

      call add2(rhs,tmp1(1),nb)

      ! artificial viscosity

      if (ifavisc) then
c        call mxm(au0,nb+1,u,nb+1,tmp1,1)
         do i=1,nb
            tmp1(i)=au(i,i)*u(i,1)
         enddo

         a=5.
         s=3.
         pad=.05

         s=-s/ad_re

         call cmult(tmp1,s,nb+1)

         call rzero(tmp2,nb+1)

         eps=1.e-2

         do i=1,nb
            um=(umax(i)+umin(i))*.5
            ud=(umax(i)-umin(i))*.5*(1.+pad)
            d=(u(i,1)-um)/ud
c           tmp2(i)=(cosh(d*acosh(2.))-1.)**a
            if (u(i,1).gt.umax(i)) then
               d=(u(i,1)/umax(i)-1.)/(1+pad)
c              tmp2(i)=d*d
c              tmp2(i)=d
               tmp2(i)=exp(d)-1.
c              tmp2(i)=exp(d*d)-1.
c              tmp2(i)=log(d)
            endif
            if (u(i,1).lt.umin(i)) then
               d=(u(i,1)/umin(i)-1.)/(1+pad)
c              tmp2(i)=d*d
c              tmp2(i)=d
               tmp2(i)=exp(d)-1.
c              tmp2(i)=exp(d*d)-1.
c              tmp2(i)=log(d)
            endif
         enddo

         call addcol3(rhs,tmp1(1),tmp2(1),nb)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setuavg

      include 'SIZE'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scravg/ ux(lt),uy(lt),uz(lt)

      if (ad_step.eq.navg_start) then
         call rzero(ua,nb+1)
         call rzero(u2a,(nb+1)**2)
      endif

      call add2(ua,u,nb+1)

      do j=0,nb
      do i=0,nb
         u2a(i,j)=u2a(i,j)+u(i,1)*u(j,1)
      enddo
      enddo

      if (ad_step.eq.ad_nsteps) then
         s=1./real(ad_nsteps-navg_start-1)
         call cmult(ua,s,nb+1)
         call cmult(u2a,s,(nb+1)**2)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine settavg

      include 'SIZE'
      include 'MOR'
      include 'AVG'

      if (ad_step.eq.navg_start) then
         call rzero(uta,nb+1)
         call rzero(uuta,(nb+1)**2)
         call rzero(utua,(nb+1)**2)
         call rzero(ut2a,(nb+1)**2)
      endif

      call add2(uta,ut,nb+1)

      do j=0,nb
      do i=0,nb
         uuta(i,j)=uuta(i,j)+u(i,1)*ut(j,1)
         utua(i,j)=utua(i,j)+u(j,1)*ut(i,1)
         ut2a(i,j)=ut2a(i,j)+ut(j,1)*ut(i,1)
      enddo
      enddo

      if (ad_step.eq.ad_nsteps) then
         s=1./real(ad_nsteps-navg_start+1)
         call cmult(uta,s,nb+1)
         call cmult(uuta,s,(nb+1)**2)
         call cmult(utua,s,(nb+1)**2)
         call cmult(ut2a,s,(nb+1)**2)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setuj

      include 'SIZE'
      include 'MOR'

      if (ad_step.eq.(navg_start+1)) then
         call copy(uj(0,1),u(0,3),nb+1)
         call copy(uj(0,2),u(0,2),nb+1)
         call copy(uj(0,3),u(0,1),nb+1)
      endif
      if (ad_step.eq.ad_nsteps) then
         call copy(uj(0,4),u(0,3),nb+1)
         call copy(uj(0,5),u(0,2),nb+1)
         call copy(uj(0,6),u(0,1),nb+1)
         do k=1,6
         do j=0,nb
         do i=0,nb
            u2j(i,j,k)=uj(i,k)*uj(j,k)
         enddo
         enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine settj

      include 'SIZE'
      include 'MOR'

      if (ad_step.eq.(navg_start+1)) then
         call copy(utj(0,1),ut(0,3),nb+1)
         call copy(utj(0,2),ut(0,2),nb+1)
         call copy(utj(0,3),ut(0,1),nb+1)
      endif

      if (ad_step.eq.ad_nsteps) then
         call copy(utj(0,4),ut(0,3),nb+1)
         call copy(utj(0,5),ut(0,2),nb+1)
         call copy(utj(0,6),ut(0,1),nb+1)

         do k=1,6
         do j=0,nb
         do i=0,nb
            uutj(i,j,k)=uj(i,k)*utj(j,k)
            utuj(i,j,k)=uj(j,k)*utj(i,k)
         enddo
         enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine seth(flu,a,b,ad_diff)

      include 'SIZE'
      include 'MOR'

      real flu(nb,nb),a(nb,nb),b(nb,nb)

      if (ad_step.le.3) then
         call cmult2(flu,b,ad_beta(1,ad_step)/ad_dt,nb*nb)
         call add2s2(flu,a,ad_diff,nb*nb)
      endif
         
      return
      end
c-----------------------------------------------------------------------
      subroutine hybrid_advance(rhs,uu,helm,invhelm,amax,amin,
     $                          adis,bpar,bstep,copt_count) 

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'  

      real helm(nb,nb),invhelm(nb,nb)
      real uu(nb),rhs(0:nb),rhstmp(0:nb)
      real amax(nb),amin(nb),adis(nb)
      real bpar,bctol
      integer bstep,chekbc,copt_count

      chekbc=0

      call copy(rhstmp,rhs,nb+1)
      call dgetrs('N',nb,1,invhelm,lub,ipiv,rhstmp(1),nb,info)

      bctol = 1e-12
      do ii=1,nb
         if ((rhstmp(ii)-amax(ii)).ge.bctol) then
            chekbc = 1
         elseif ((amin(ii)-rhstmp(ii)).ge.bctol) then
            chekbc = 1
         endif
      enddo

      if (chekbc.eq.1) then
         copt_count = copt_count + 1
         call BFGS(rhs(1),uu,helm,invhelm,amax,amin,adis,
     $   bpar,bstep)
      else
         call copy(rhs,rhstmp,nb+1)
      endif

      return
      end
