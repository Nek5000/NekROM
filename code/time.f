c-----------------------------------------------------------------------
      subroutine rom_step

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real rhs(0:nb)
      logical ifdebug

      ifdebug=.true.
      ifdebug=.false.

      if (ad_step.eq.1) then
         step_time = 0.
         solve_time=0.
         lu_time=0.
      endif

      last_time = dnekclock()

      n=lx1*ly1*lz1*nelt

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
         call lu(fluv,nb,nb,irv,icv)
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
            call solve(rhs(1),fluv,1,nb,nb,irv,icv)
         else
            eps=.20
            damp=1.-eps*ad_dt
            do i=1,nb
            if (rhs(i).gt.umax(i)) rhs(i)=umax(i)+(rhs(i)-umax(i))*damp
            if (rhs(i).lt.umin(i)) rhs(i)=umin(i)+(rhs(i)-umin(i))*damp
            enddo
         endif
      else if (isolve.eq.1) then ! constrained solve
         call BFGS_freeze(rhs(1),helmu,invhelmu,umax,umin,udis,1e-3,4) 
c        call BFGS(rhs(1),helmu,invhelmu,umax,umin,udis,1e-3,4) 
      else
         call exitti('incorrect isolve specified...')
      endif
      solve_time=solve_time+dnekclock()-ttime

      do i=0,nb
         if (ifdebug) write (6,*) i,rhs(i),'sol'
      enddo

c     if (ifdebug) call exitt0

      call shift3(u,rhs,nb+1)

      step_time=step_time+dnekclock()-last_time

      call pp

      return
      end
c-----------------------------------------------------------------------
      subroutine pp

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrrstep/ t1(lt),t2(lt),t3(lt),work(lt)
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      real vort(lt)

      call setavg
      call setj

      if (mod(ad_step,ad_qstep).eq.0) then
         if (ifctke) call ctke
         if (ifcdrag) call cdrag
         call cnuss
      endif

      if (mod(ad_step,ad_iostep).eq.0) then
         if (nio.eq.0) then
            write (6,*)'ad_step:',ad_step,ad_iostep,npp,nid,step_time
            if (ad_step.eq.ad_nsteps) then
               do j=1,nb
                  write(6,*) j,u(j,1),'final'
               enddo
            else
               do j=1,nb
                  write(6,*) j,u(j,1)
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

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_step_t

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

c     Matrices and vectors for advance
      real tmp(0:nb),rhs(0:nb)

      common /scrrstep/ t1(lt),t2(lt),t3(lt),work(lt)

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      if (ad_step.eq.1) then
         step_time = 0.
      endif

      last_time = dnekclock()

      n=lx1*ly1*lz1*nelt

      icount = min0(ad_step,3)

      rhs(0)=1.
      call setr_t(rhs(1),icount)

      call seth(flut,at,bt,1./ad_pe)
      if (ad_step.eq.3) call dump_serial(flut,nb*nb,'ops/ht ',nid)
      if (ad_step.le.3) then
         call lu(flut,nb,nb,irt,ict)
         call copy(invhelmt,flut,nb*nb)
      endif


      if (isolve.eq.0) then ! standard matrix inversion
         call solve(rhs(1),flut,1,nb,nb,irt,ict)
      else if (isolve.eq.1) then ! constrained solve
         call BFGS_freeze(rhs(1),helmt,invhelmt,tmax,tmin,tdis,1e-3,4) 
c        call BFGS(rhs(1),helmt,invhelmt,tmax,tmin,tdis,1e-3,4) 
      else
         call exitti('incorrect isolve specified...')
      endif

      call shift3(ut,rhs,nb+1)

      step_time=step_time+dnekclock()-last_time

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
      subroutine evalc(cu,cl,icl,uu)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real cu(nb)

      integer icalld
      save    icalld
      data    icalld /0/

      common /scrc/ work(lx1*ly1*lz1*lelt)

      real cl(lcloc),uu(0:nb)
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

         do l=1,ncloc
            i=icl(1,l)
            j=icl(2,l)
            k=icl(3,l)
            cu(i)=cu(i)+cl(l)*uu(j)*u(k,1)
         enddo

         call gop(cu,work,'+  ',nb)
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

      call evalc(tmp(1),ctl,ictl,ut)
      call shift3(ctr,tmp(1),nb)

      call mxm(ctr,nb,ad_alpha(1,icount),3,tmp(1),1)

      call sub2(rhs,tmp(1),nb)

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

      call evalc(tmp1(1),cul,icul,u)
      call chsign(tmp1(1),nb)

      if (ifbuoy) then
         call mxm(but0,nb+1,ut(0,1),nb+1,tmp2(0),1)
         call add2(tmp1(1),tmp2(1),nb)
      else if (ifforce) then
         call add2(tmp1(1),bg(1),nb)
      endif

      call shift3(fu,tmp1(1),nb)

      call mxm(fu,nb,ad_alpha(1,icount),3,tmp1(1),1)

      call add2(rhs,tmp1(1),nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine setavg

      include 'SIZE'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scravg/ ux(lt),uy(lt),uz(lt)

      if (ad_step.eq.1) then
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
         s=1./real(ad_nsteps)
         call cmult(ua,s,nb+1)
         call cmult(u2a,s,(nb+1)**2)

         call reconv(ux,uy,uz,ua)
         call outpost(ux,uy,uz,pavg,tavg,'avg')

         call reconu_rms(ux,uy,uz,u2a)
         call outpost(ux,uy,uz,pavg,tavg,'rms')
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setj

      include 'SIZE'
      include 'MOR'

      if (ad_step.eq.3) call copy(uj,u,3*(nb+1))
      if (ad_step.eq.ad_nsteps) then
         call copy(uj(0,4),u,3*(nb+1))
         do k=1,6
         do j=0,nb
         do i=0,nb
            u2j(i,j,k)=uj(i,k)*uj(j,k)
         enddo
         enddo
         enddo
         s=1./real(ad_nsteps)
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
