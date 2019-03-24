c-----------------------------------------------------------------------
      subroutine rom_step_v

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

c     Matrices and vectors for advance
      real tmp(0:nb),rhs(nb)

      common /scrrstep/ t1(lt),t2(lt),t3(lt),work(lt)

c     Variable for vorticity
      real vort(lt,3)

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      if (ad_step.eq.1) then
         step_time = 0.
      endif

      last_time = dnekclock()

      n  = lx1*ly1*lz1*nelt

      icount = min0(ad_step,3)

      call setr_v(rhs,icount)

      if (ad_step.le.3) then
         call cmult2(fluv,bv,ad_beta(1,icount)/ad_dt,nb*nb)
         call add2s2(fluv,av,1/ad_re,nb*nb)
         call lu(fluv,nb,nb,irv,icv)
      endif

      if (isolve.eq.0) then ! standard matrix inversion
         call solve(rhs,fluv,1,nb,nb,irv,icv)
      else if (isolve.eq.1) then ! constrained solve
         call BFGS_freeze 
      else
         call exitti('incorrect isolve specified...')
      endif

      call shiftu(rhs)
      call setavg
      call setj

      if (ifcdrag) call cdrag
c     call comp_rms ! old

      step_time=step_time+dnekclock()-last_time

      if (ifdump) then
         time=time+dt
      else
         call reconv(vx,vy,vz,u)
      endif

      if (mod(ad_step,ad_iostep).eq.0) then
         if (nio.eq.0) then
            write (6,*)'ad_step:',ad_step,ad_iostep,npp,nid,step_time
            if (ad_step.eq.ad_nsteps) then
               do j=1,nb
                  write(6,*) j,u(j,1),'final'
               enddo
               write (6,*) 'step_time: ',step_time
            else
               do j=1,nb
                  write(6,*) j,u(j,1)
               enddo
            endif
         endif

         ifdump=ifdump

         if (ifdump) then
            idump=ad_step/ad_iostep
            call reconv(vx,vy,vz,u)
            call opcopy(t1,t2,t3,vx,vy,vz)

            if (ifrom(2)) then
               call recont(vort,ut)
            else if (ifheat) then
               call copy(vort,t,n)
            else
               call comp_vort3(vort,work1,work2,t1,t2,t3)
            endif

            ifto = .true. ! turn on temp in fld file
            call outpost(vx,vy,vz,pavg,vort,'rom')
            if (nio.eq.0) write (6,*) 'inside ifdump'
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
      real tmp(0:nb),rhs(nb)

      common /scrrstep/ t1(lt),t2(lt),t3(lt),work(lt)

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      if (ad_step.eq.1) then
         step_time = 0.
      endif

      last_time = dnekclock()

      n=lx1*ly1*lz1*nelt

      icount = min0(ad_step,3)

      call setr_t(rhs,icount)

      if (ad_step.le.3) then
         call cmult2(flut,bt,ad_beta(1,icount)/ad_dt,nb*nb)
         call add2s2(flut,at,1/ad_pe,nb*nb)
         call lu(flut,nb,nb,irt,ict)
      endif

      if (isolve.eq.0) then ! standard matrix inversion
         call solve(rhs,flut,1,nb,nb,irt,ict)
      else if (isolve.eq.1) then ! constrained solve
         !call csolve(rhs,flu,...
      else
         call exitti('incorrect isolve specified...')
      endif

      call shiftt(rhs)

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

      common /scrk4/ work(lx1*ly1*lz1*lelt)

      real cl(lcloc),uu(0:nb)
      integer icl(3,lcloc)

      if (icalld.eq.0) then
         evalc_time=0.
         icalld=1
      endif

      stime=dnekclock()

      l=1

      call rzero(cu,nb)

      do l=1,ncloc
         i=icl(1,l)
         j=icl(2,l)
         k=icl(3,l)
         cu(i)=cu(i)+cl(l)*uu(j)*u(k,1)
      enddo

      call gop(cu,work,'+  ',nb)

      evalc_time=evalc_time+dnekclock()-stime

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

      call copy(ctr(1,3),ctr(1,2),nb)
      call copy(ctr(1,2),ctr(1,1),nb)

      call evalc(ctr,ctl,ictl,ut)

      call mxm(ctr,nb,ad_alpha(1,icount),3,tmp(1),1)

      call sub2(rhs,tmp(1),nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine setr_v(rhs,icount)

      include 'SIZE'
      include 'MOR'

      common /scrrhs/ tmp(0:nb)

      real rhs(nb)

      call mxm(u,nb+1,ad_beta(2,icount),3,tmp,1)
c     call mxm(bv0,nb+1,tmp,nb+1,rhs,1)
      call mxm(bv,nb,tmp(1),nb,rhs,1)

      call cmult(rhs,-1.0/ad_dt,nb)

      s=-1.0/ad_re

c     call add2s2(rhs,av0,s,nb+1) ! not working...
      do i=1,nb
         rhs(i)=rhs(i)+s*av0(i,0)
      enddo

      call copy(cvr(1,3),cvr(1,2),nb)
      call copy(cvr(1,2),cvr(1,1),nb)

      call evalc(cvr,cvl,icvl,u)

      call mxm(cvr,nb,ad_alpha(1,icount),3,tmp(1),1)

      call sub2(rhs,tmp(1),nb)

      if (ifbuoy) then
         call mxm(bvt0,nb+1,ut(0,1),nb+1,tmp(0),1)
         call add2(rhs,tmp(1),nb)
      else if (ifforce) then
         call add2(rhs,bg(1),nb)
      endif

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

         call opzero(ux,uy,uz)
         n=lx1*ly1*lz1*nelv
         do j=0,nb
         do i=0,nb
            call col3(ubt,ub(1,i),ub(1,j),n)
            call col3(vbt,vb(1,i),vb(1,j),n)
            if (ldim.eq.3) call col3(wbt,wb(1,i),wb(1,j),n)
            call opadds(ux,uy,uz,ubt,vbt,wbt,u2a(i,j),n,1)
         enddo
         enddo
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
