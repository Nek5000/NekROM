c-----------------------------------------------------------------------
      subroutine bdfext_step

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      common /scrbdfext/ rhs(0:lb,2),rhstmp(0:lb),
     $                   utmp1(0:lb),utmp2(0:lb)

      common /eddyma/ cedm((lb**2-1)*ledvis+1),cedvec((lb-1)*ledvis+1)

      logical ifdebug
      integer chekbc

      chekbc=0

      ifdebug=.true.
      ifdebug=.false.
      jfield=ifield

      ulast_time = dnekclock()

      n=lx1*ly1*lz1*nelt

      icount = min(max(1,ad_step),3)

      rhs(0,1)=1.
      rhs(0,2)=1.

c     if (icount.le.2) then
      if (.false.) then
         if (ifrom(1)) call setr_v(rhs(1,1),icount)
         if (ifrom(2)) call setr_t(rhs(1,2),icount)
         call rk4_setup
         call copy(urki(1),u(1),nb)
         if (ifrom(2)) call copy(urki(nb+1),ut(1),nb)
         nrk=nb
         if (ifrom(2)) nrk=nb*2

         call rk_step(urko,rtmp1,urki,time,ad_dt,grk,rtmp2,nrk)

         if (ifrom(1)) then
            call copy(rhs(1,1),urko,nb)
            call shift(u,rhs,nb+1,3)
         endif
         if (ifrom(2)) then
            call copy(rhs(1,2),urko(nb+1),nb)
            call shift(ut,rhs(0,2),nb+1,3)
         endif
         return
      endif

      ifield=2
      if (ifrom(2)) then
         if (ad_step.le.3) then
            ttime=dnekclock()
            call seth(hlm(1,2),at,bt,1./ad_pe)
            if (ad_step.eq.3)
     $         call dump_serial(hlm(1,2),nb*nb,'ops/ht ',nid)
            do j=1,nb-nplay
            do i=1,nb-nplay
               hlm(i+(j-1)*(nb-nplay),2)=hlm(i+nplay+(j+nplay-1)*nb,2)
            enddo
            enddo
            if (.not.ifcomb) then
               if (ifdecpl) then
                  call copy(hinv(1,2),hlm(1,2),(nb-nplay)**2)
                  call diag(hinv(1,2),wt(1,2),rhs(1,2),nb)
               else
                  call invmat(hinv(1,2),hlu(1,2),hlm(1,2),
     $            ihlu(1,2),ihlu2(1,2),nb-nplay)
                  call rzero(wt(1,2),(nb-nplay)**2)
                  do i=1,nb-nplay
                     wt(i+(nb-nplay)*(i-1),2)=1.
                  enddo
               endif
            endif
            lu_time=lu_time+dnekclock()-ttime
            call update_k(uk,ukp,tk,tkp)
         endif

         call rom_userfop

         call setr_t(rhs(1,2),icount)
         call rom_userrhs(rhs(1,2))

         ttime=dnekclock()
         if (.not.ifcomb) then
         if (isolve.eq.0) then
            call mxm(wt(1,2),nb,rhs(1,2),nb,rhstmp(1),1)
            call mxm(hinv(1,2),nb,rhstmp(1),nb,rhs(1,2),1)
            if (ifdecpl) then
            eps=1.e-3
            do i=1,nb
               if (rhs(i,2).gt.tpmax(i)) rhs(i,2)=tpmax(i)-tpdis(i)*eps
               if (rhs(i,2).lt.tpmin(i)) rhs(i,2)=tpmin(i)+tpdis(i)*eps
            enddo
            endif
            call mxm(rhs(1,2),1,wt(1,2),nb,rhstmp(1),nb)
            call copy(rhs(1,2),rhstmp(1),nb)
         else
            scopt='tcopt'
            call mxm(ut,nb+1,ad_alpha(1,icount),icount,utmp1,1)
            call mxm(wt(1,2),nb,utmp1(1),nb,utmp2(1),1)
            eps=1.e-3
            do i=1,nb
               if (utmp2(i).gt.tpmax(i)) utmp2(i)=tpmax(i)-tpdis(i)*eps
               if (utmp2(i).lt.tpmin(i)) utmp2(i)=tpmin(i)+tpdis(i)*eps
            enddo
            call mxm(wt(1,2),nb,rhs(1,2),nb,rhstmp(1),1)
            call constrained_POD(rhstmp,hlm(1,2),hinv(1,2),utmp2(1),
     $                           tpmax,tpmin,tpdis,
     $                           tbarr0,tbarrseq,tcopt_count)
            call mxm(rhstmp(1),1,wt(1,2),nb,rhs(1,2),nb)
         endif

         if (cfloc.eq.'POST') then
         if (cftype.eq.'TFUN') then
            call pod_proj(rhs(1,2),rbf,nb,'step  ')
         else if (cftype.eq.'DIFF') then
            call pod_df(rhs(1,2))
         endif
         endif

         tsolve_time=tsolve_time+dnekclock()-ttime
         endif
      endif

      ifield=1
      if (ifrom(1)) then
         if (ad_step.le.3) then
            ttime=dnekclock()
            call seth(hlm,au,bu,1./ad_re)
            if (ifedvs) then
               edv(0) = 1
               edv(1) = 1
               call addeddy(hlm,cedm,cedvec,cedd,edv)
            endif
            if (ad_step.eq.3) call dump_serial(hlm,nb*nb,'ops/hu ',nid)
            do j=1,nb-nplay
            do i=1,nb-nplay
               hlm(i+(j-1)*(nb-nplay),1)=hlm(i+nplay+(j+nplay-1)*nb,1)
            enddo
            enddo
            if (nio.eq.0) write (6,*) 'check ifdecpl',ifdecpl,'cp3'
            if (.not.ifcomb) then
               if (ifdecpl) then
                  call copy(hinv,hlm,(nb-nplay)**2)
                  call diag(hinv,wt,rhs(1,1),nb)
               else
                  call invmat(hinv,hlu,hlm,ihlu,ihlu2,nb-nplay)
               if (nio.eq.0) write (6,*) 'check nplay',nplay,'cp3'
                  call rzero(wt,(nb-nplay)**2)
                  do i=1,nb-nplay
                     wt(i+(nb-nplay)*(i-1),1)=1.
                  enddo
               endif
            endif
            lu_time=lu_time+dnekclock()-ttime
            call update_k(uk,ukp,tk,tkp)
         endif

         call rom_userfop

         call setr_v(rhs(1,1),icount)
         call rom_userrhs(rhs(1,1))

         ttime=dnekclock()
         if (.not.ifcomb) then
         if ((isolve.eq.0).or.(icopt.eq.2)) then ! standard matrix inversion
            call mxm(wt,nb,rhs(1,1),nb,rhstmp(1),1)
            call mxm(hinv,nb,rhstmp(1),nb,rhs(1,1),1)
            if (ifdecpl) then
            eps=1.e-3
            do i=1,nb
               if (rhs(i,1).gt.upmax(i)) rhs(i,1)=upmax(i)-updis(i)*eps
               if (rhs(i,1).lt.upmin(i)) rhs(i,1)=upmin(i)+updis(i)*eps
            enddo
            endif
            call mxm(rhs(1,1),1,wt,nb,rhstmp(1),nb)
            call copy(rhs(1,1),rhstmp(1),nb)
         else 
            scopt='ucopt'
            call mxm(u,nb+1,ad_alpha(1,icount),icount,utmp1,1)
            call mxm(wt,nb,utmp1(1),nb,utmp2(1),1)
            eps=1.e-3
            do i=1,nb
               if (utmp2(i).gt.upmax(i)) utmp2(i)=upmax(i)-updis(i)*eps
               if (utmp2(i).lt.upmin(i)) utmp2(i)=upmin(i)+updis(i)*eps
            enddo
            call mxm(wt,nb,rhs(1,1),nb,rhstmp(1),1)
            call constrained_POD(rhstmp,hlm,hinv,utmp2(1),upmax,upmin,
     $                           updis,ubarr0,ubarrseq,ucopt_count)
            call mxm(rhstmp(1),1,wt,nb,rhs(1,1),nb)

         endif
         endif

         if (cfloc.eq.'POST') then
         if (cftype.eq.'TFUN') then
            call pod_proj(rhs(1,1),rbf,nb,'step  ')
         else if (cftype.eq.'DIFF') then
            call pod_df(rhs(1,1))
         endif
         endif

         solve_time=solve_time+dnekclock()-ttime
      endif

      if (ifcomb) then
         if (ad_step.le.3) then
            call add2(hlm(1,1),hlm(1,2),(nb+1)**2)
            if (ifdecpl) then
               call copy(hinv(1,1),hlm(1,1),(nb-nplay)**2)
               call diag(hinv(1,1),wt(1,1),utmp1,nb)
            else
               call invmat(hinv(1,1),hlu(1,1),hlm(1,1),
     $         ihlu(1,1),ihlu2(1,1),nb-nplay)
               call rzero(wt(1,1),(nb-nplay)**2)
               do i=1,nb-nplay
                  wt(i+(nb-nplay)*(i-1),1)=1.
               enddo
            endif
         endif
         call add2(rhs(1,1),rhs(1,2),(nb+1)**2)

         call mxm(wt,nb,rhs(1,1),nb,rhstmp(1),1)
         call mxm(hinv,nb,rhstmp(1),nb,rhs(1,1),1)
         if (ifdecpl) then
            do i=1,nb
               if (rhs(i,1).gt.upmax(i)) rhs(i,1)=upmax(i)-updis(i)*eps
               if (rhs(i,1).lt.upmin(i)) rhs(i,1)=upmin(i)+updis(i)*eps
            enddo
         endif
         call mxm(rhs(1,1),1,wt,nb,rhstmp(1),nb)
         call copy(rhs(1,1),rhstmp(1),nb)
         call copy(rhs(1,2),rhs(1,1),nb)
      endif

      if (cfloc.eq.'POST'.and.cftype.eq.'POLY')
     $   call apply_les_imp(rhs(0,1),rhs(0,2),rdft,fles1,fles2,rtmp1)

      if (ifrom(2)) call shift(ut,rhs(0,2),nb+1,5)
      if (ifrom(1)) call shift(u,rhs,nb+1,5)

      ifield=jfield

      ustep_time=ustep_time+dnekclock()-ulast_time

      return
      end
c-----------------------------------------------------------------------
      subroutine post

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrrstep/ t1(lt),t2(lt),t3(lt),work(lt)
      common /scrbdfext/ rhs(0:lb,2),rhstmp(0:lb),
     $                   utmp1(0:lb),utmp2(0:lb)
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      save icalld
      data icalld /0/

      real vort(lt)

      call rom_userchk

      if (icalld.eq.0) then
         post_time=0.
         icalld=1
      endif

      call nekgsync
      tttime=dnekclock()

      if (ifrom(1)) then
         call setuavg(ua,u2a,u)
         call count_gal(num_galu,anum_galu,rhstmp(1),upmax,upmin,
     $   1e-16,nb)
      endif

      if (ifrom(2)) then
         call settavg(uta,uuta,utua,ut2a,u,ut)
         call count_gal(num_galt,anum_galt,ut(1),tmax,tmin,1e-16,nb)
      endif

      if (ifei) call set_time_avg_resid_coef

      if (mod(ad_step,ad_qstep).eq.0) then
         if (ifctke) call ctke
         if (ifcdrag) call cdrag
         call cnuss
c        call cubar
      endif

      if (mod(ad_step,ad_iostep).eq.0) then
         if (ifrom(1)) then
            do j=1,nb
               if (nio.eq.0) write (6,*) j,time,u(j),'romu'
            enddo
         endif

         if (ifrom(2)) then
            do j=1,nb
               if (nio.eq.0) write (6,*) j,time,ut(j),'romt'
            enddo
         endif

         if (ifrecon) then
            jstep=istep
            call reconv(vx,vy,vz,u)
            call opcopy(t1,t2,t3,vx,vy,vz)

            if (ifrom(2)) then
               call recont(vort,ut)
            else
               call comp_vort3(vort,work1,work2,t1,t2,t3)
            endif

            istep=ad_step

            ifto = .true. ! turn on temp in fld file
            call outpost(vx,vy,vz,pavg,vort,'rom')
            istep=jstep
         endif
      endif

      if (ad_step.eq.ad_nsteps) then
         if (nio.eq.0) then
            if (ifrom(1)) then
               do j=1,nb
                  write (6,*) j,num_galu(j)/ad_nsteps,'num_galu'
               enddo
               write (6,*) anum_galu/ad_nsteps,'anum_galu'
            endif
            if (ifrom(2)) then
               do j=1,nb
                  write (6,*) j,num_galt(j)/ad_nsteps,'num_galt'
               enddo
               write (6,*) anum_galt/ad_nsteps,'anum_galt'
            endif
         endif
      endif


      call nekgsync
      postu_time=postu_time+dnekclock()-tttime

    1 format(i8,1p4e13.5,2i5,a)
    2 format(i8,1p1e13.5,i5,1p2e13.5,' play_romu')

      return
      end
c-----------------------------------------------------------------------
      subroutine set_bdf_coef(ad_alpha,ad_beta)

      real ad_alpha(3,3),ad_beta(4,3)

      call rzero(ad_alpha,3*3)
      call rzero(ad_beta,3*4)

      ad_beta(1,1)  =  1.
      ad_beta(2,1)  = -1.

      ad_beta(1,2)  =  1.5
      ad_beta(2,2)  = -2
      ad_beta(3,2)  =  0.5

      ad_beta(1,3)  =  11./6
      ad_beta(2,3)  = -3.
      ad_beta(3,3)  =  1.5
      ad_beta(4,3)  = -1./3

      ad_alpha(1,1) =  1.

      ad_alpha(1,2) =  2.
      ad_alpha(2,2) = -1.

      ad_alpha(1,3) =  3.
      ad_alpha(2,3) = -3.
      ad_alpha(3,3) =  1.

      return
      end
c-----------------------------------------------------------------------
      subroutine evalc(cu,cm,cl,uu,tt)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      common /evalctmp/ ucft(0:lb)

      real cu(nb)
      real uu(0:nb),tt(0:nb)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      real cm(ic1:ic2,jc1:jc2)

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
         call rzero(cu,nb)
         if (ncloc.ne.0) then
            if ((kc2-kc1).lt.64.and.(jc2-jc1).lt.64
     $          .and.cfloc.eq.'NONE') then
               call mxm(cl,(ic2-ic1+1)*(jc2-jc1+1),
     $                  uu(kc1),(kc2-kc1+1),cm,1)
               call mxm(cm,(ic2-ic1+1),tt(jc1),(jc2-jc1+1),cu(ic1),1)
            else
               call copy(ucft,uu,nb+1)

               if (cfloc.eq.'CONV') then
               if (cftype.eq.'TFUN') then
                  call pod_proj(ucft(1),rbf,nb,'step  ')
               else if (cftype.eq.'DIFF') then
                  call pod_df(ucft(1))
               endif
               endif

               do k=kc1,kc2
               do j=jc1,jc2
               do i=ic1,ic2
                  cu(i)=cu(i)+cl(i,j,k)*tt(j)*ucft(k)
               enddo
               enddo
               enddo
            endif
         endif
         call gop(cu,work,'+  ',nb)
      endif

      call nekgsync

      evalc_time=evalc_time+dnekclock()-stime

      return
      end
c-----------------------------------------------------------------------
      subroutine evalc3(cu,fac_a,fac_b,fac_c,cp_weight,uu)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real cu(nb)
      real uu(0:nb)
      real ucft(0:nb)
      real fac_a((nb+1)*ntr),fac_b((nb+1)*ntr)
      real fac_c((nb+1)*ntr),cp_weight(ntr)
      real bcu(ntr)
      real cuu(ntr)
      real tmp(ntr)
      real tmpcu(0:nb)

      common /scrc/ work(max(lub,ltb))

      integer icalld
      save    icalld
      data    icalld /0/

      if (icalld.eq.0) then
         evalc_time=0.
         icalld=1
      endif

      stime=dnekclock()

      call rzero(cu,nb)
      do kk=1,ntr
         bcu(kk) = vlsc2(uu,fac_b(1+(kk-1)*(nb+1)),nb+1)
         cuu(kk) = vlsc2(u,fac_c(1+(kk-1)*(nb+1)),nb+1)
      enddo
      call col4(tmp,bcu,cuu,cp_weight,ntr) 
      call mxm(fac_a,nb+1,tmp,ntr,tmpcu,1)

      call copy(cu,tmpcu,nb)

      call nekgsync

      evalc_time=evalc_time+dnekclock()-stime

      return
      end
c-----------------------------------------------------------------------
      subroutine evalc4(cu,fac_a,fac_b,fac_c,cp_weight,cl,cj0,c0k,uu)

c     This subroutine contracts the approximated tensor from CP
c     decomposition with two vectors

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real cu(nb)
      real fac_a((nb)*ntr),fac_b((nb)*ntr)
      real fac_c((nb)*ntr),cp_weight(ntr)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      real cj0(ic1:ic2,jc1:jc2),c0k(ic1:ic2,kc1:kc2)
      real uu(0:nb)

      common /srcevalc4/ bcu(ltr),cuu(ltr),tmp(ltr),tmpcu(0:lb)

      integer icalld
      save    icalld
      data    icalld /0/

      if (icalld.eq.0) then
         evalc_time=0.
         icalld=1
      endif

      stime=dnekclock()

      call rzero(cu,nb)
      do kk=1,ntr
         bcu(kk) = vlsc2(uu(1),fac_b(1+(kk-1)*(nb)),nb)
         cuu(kk) = vlsc2(u(1),fac_c(1+(kk-1)*(nb)),nb)
      enddo
      call col4(tmp,bcu,cuu,cp_weight,ntr) 
      call mxm(fac_a,nb,tmp,ntr,tmpcu(1),1)

      do k=kc1,kc2
      do i=ic1,ic2
         cu(i)=cu(i)+c0k(i,k)*uu(0)*u(k)
      enddo
      enddo

      do j=jc1,jc2
      do i=ic1,ic2
         cu(i)=cu(i)+cj0(i,j)*uu(j)*u(0)
      enddo
      enddo

      do i=ic1,ic2
         cu(i)=cu(i)-cj0(i,0)*uu(0)*u(0)
      enddo

c     call copy(cu,tmpcu(1),nb)
      call add2(cu,tmpcu(1),nb)

      call nekgsync

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
      include 'INPUT'
      include 'TSTEP'

      common /scrrhs/ tmp(0:lb),tmp2(0:lb)

      real rhs(nb)

      call mxm(ut,nb+1,ad_beta(2,icount),3,tmp,1)
      call mxm(bt,nb,tmp(1),nb,rhs,1)

      call cmult(rhs,-1.0/ad_dt,nb)

      s=-1.0/ad_pe

      do i=1,nb
         rhs(i)=rhs(i)+s*at0(1+i)
      enddo

      if (ifadvc(2)) then
         if (ifcp) then
            if (ifcore) then
               call evalc4(tmp(1),cta,ctb,ctc,cp_tw,ctl,ctj0,ct0k,ut)
            else
               call evalc3(tmp(1),cta,ctb,ctc,cp_tw,ut)
            endif
         else
            call evalc(tmp(1),ctmp,ctl,u,ut)
         endif
c        call add2(tmp(1),st0(1),nb)
         call shift(ctr,tmp(1),nb,3)

         call mxm(ctr,nb,ad_alpha(1,icount),3,tmp(1),1)

         call sub2(rhs,tmp(1),nb)
      endif

      if (ifsource) then
         call add2(rhs,rq,nb)
         if (ifsrct) then
            rqt_time_coef(3) = rqt_time_coef(2)
            rqt_time_coef(2) = rqt_time_coef(1)
            time = (ad_step-1)*ad_dt
            call userq(1,1,1,1)
            rqt_time_coef(1) = ft
            if (ad_step .le. 4) then
               write(6,*)'ad_step', ad_step
               write(6,*)'extrapolation points:', rqt_time_coef(1:3)
            endif
            rqttcext = glsc2(ad_alpha(1,icount),rqt_time_coef,3)
            call add2s2(rhs,rqt,rqttcext,nb)
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setr_v(rhs,icount)

      include 'SIZE'
      include 'MOR'

      common /scrrhs/ tmp1(0:lb),tmp2(0:lb)
      common /eddyma/ cedm((lb**2-1)*ledvis+1),cedvec((lb-1)*ledvis+1)

      real rhs(nb)

      call mxm(u,nb+1,ad_beta(2,icount),3,tmp1,1)
      call mxm(bu,nb,tmp1(1),nb,rhs,1)

      call cmult(rhs,-1.0/ad_dt,nb)

      s=-1.0/ad_re

      do i=1,nb
         rhs(i)=rhs(i)+s*au0(1+i)
      enddo
      if (ifedvs) then
         do i=1,nb
            rhs(i)=rhs(i)-cedvec(i)
         enddo
      endif

      if (ifcp) then
         if (ifcore) then 
            call evalc4(tmp1(1),cua,cub,cuc,cp_uw,cul,cuj0,cu0k,u)
         else
            call evalc3(tmp1(1),cua,cub,cuc,cp_uw,u)
         endif 
      else
         call evalc(tmp1(1),ctmp,cul,u,u)
      endif

      call chsign(tmp1(1),nb)

      if (ifbuoy) then
         call mxm(buxt0,nb+1,ut,nb+1,tmp2(0),1)
         call add2s2(tmp1(1),tmp2(1),-gx,nb)

         call mxm(buyt0,nb+1,ut,nb+1,tmp2(0),1)
         call add2s2(tmp1(1),tmp2(1),-gy,nb)

         if (ldim.eq.3) then
            call mxm(buzt0,nb+1,ut,nb+1,tmp2(0),1)
            call add2s2(tmp1(1),tmp2(1),-gz,nb)
         endif
      else if (ifforce) then
         call add2(tmp1(1),rf(1),nb)
      endif

      call shift(fu,tmp1(1),nb,3)


      call mxm(fu,nb,ad_alpha(1,icount),3,tmp1(1),1)

      call add2(rhs,tmp1(1),nb)

      ! artificial viscosity

      if (ifavisc) then
c        call mxm(au0,nb+1,u,nb+1,tmp1,1)
         do i=1,nb
            tmp1(i)=au(i+nb*(i-1))*u(i)
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
            d=(u(i)-um)/ud
c           tmp2(i)=(cosh(d*acosh(2.))-1.)**a
            if (u(i).gt.umax(i)) then
               d=(u(i)/umax(i)-1.)/(1+pad)
c              tmp2(i)=d*d
c              tmp2(i)=d
               tmp2(i)=exp(d)-1.
c              tmp2(i)=exp(d*d)-1.
c              tmp2(i)=log(d)
            endif
            if (u(i).lt.umin(i)) then
               d=(u(i)/umin(i)-1.)/(1+pad)
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
      subroutine setuavg(s1,s2,t1)

      ! set average quantities involving velocity coefficients

      ! s1 := <uj>
      ! s2 := <ujxuj>
      ! t1 := uj

      include 'SIZE'
      include 'MOR'

      real s1(0:nb),s2(0:nb,0:nb)
      real t1(0:nb)

      if (ad_step.eq.navg_step) then
         call rzero(s1,nb+1)
         call rzero(s2,(nb+1)**2)
      endif

      call add2(s1,t1,nb+1)

      do j=0,nb
      do i=0,nb
         s2(i,j)=s2(i,j)+t1(i)*t1(j)
      enddo
      enddo

      if (ad_step.eq.ad_nsteps) then
         s=1./real(ad_nsteps-(navg_step-1))
         call cmult(s1,s,nb+1)
         call cmult(s2,s,(nb+1)**2)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine settavg(s1,s2,s3,s4,t1,t2)

      ! set average quantities involving temperature coefficients

      ! s1 := <uj>
      ! s2 := <ujxtj>
      ! s3 := <tjxuj>
      ! s4 := <tjxtj>
      ! t1 := uj
      ! t2 := tj

      include 'SIZE'
      include 'MOR'

      real s1(0:nb),s2(0:nb,0:nb),s3(0:nb,0:nb),s4(0:nb,0:nb)
      real t1(0:nb),t2(0:nb)

      if (ad_step.eq.navg_step) then
         call rzero(s1,nb+1)
         call rzero(s2,(nb+1)**2)
         call rzero(s3,(nb+1)**2)
         call rzero(s4,(nb+1)**2)
      endif

      call add2(s1,t2,nb+1)

      do j=0,nb
      do i=0,nb
         s2(i,j)=s2(i,j)+t1(i)*t2(j)
         s3(i,j)=s3(i,j)+t1(j)*t2(i)
         s4(i,j)=s4(i,j)+t2(j)*t2(i)
      enddo
      enddo

      if (ad_step.eq.ad_nsteps) then
         s=1./real(ad_nsteps-navg_step+1)
         call cmult(s1,s,nb+1)
         call cmult(s2,s,(nb+1)**2)
         call cmult(s3,s,(nb+1)**2)
         call cmult(s4,s,(nb+1)**2)
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
         if (ifhelm) then
            call add2s2(flu,b,ad_mu,nb*nb)
         endif
      endif
         
      return
      end
c-----------------------------------------------------------------------
      subroutine hybrid_advance(rhs,helm,invhelm,uu,amax,amin,adis,
     $                        bpar,bstep,copt_count,tol_box,ifdiag,nb)

      include 'SIZE'
      include 'TOTAL'

      real rhs(0:nb)
      real helm(nb,nb),invhelm(nb,nb)
      real uu(nb),amax(nb),amin(nb),adis(nb)
      real bpar,tol_box
      real ipiv(nb)
      integer bstep,chekbc,copt_count,nb
      logical ifdiag

      real rhstmp(0:nb)

      call copy(rhstmp,rhs,nb+1)
      if (ifdiag) then
         call col2(rhstmp(1),invhelm,nb)
      else
         call mxm(invhelm,nb,rhs(1),nb,rhstmp(1),1)
      endif

      call check_box(chekbc,rhstmp(1),amax,amin,tol_box,nb)

      if (chekbc.eq.1) then
         copt_count = copt_count + 1
         call IPM(rhs(1),uu,helm,invhelm,amax,amin,adis,
     $   bpar,bstep,tol_box,ifdiag)
      else
         call copy(rhs,rhstmp,nb+1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine constrained_POD(rhs,hh,invhh,uu,amax,amin,adis,
     $                          bpar,bstep,copt_count) 

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      common /scrcopt/ helm(lb**2,2),invhelm(lb**2,2)

      real rhs(0:nb)
      real hh(nb**2),invhh(nb**2)
      real uu(nb),amax(nb),amin(nb),adis(nb)
      real bpar,invhelm
      real tmp(nb),rhstmp(0:nb)

      integer bstep,chekbc,copt_count
      logical ifdiag
      integer checkdiag

      call check_diag(checkdiag,ifdiag,hh,invhh,nb)

      if (isolve.eq.1) then 
         ! constrained solver with inverse update
         call IPM(rhs(1),uu,helm,invhelm,amax,amin,adis,
     $   bpar,bstep,box_tol,ifdiag)

      else if (isolve.eq.2) then 
         ! mix constrained solver with inverse update
         call hybrid_advance(rhs,helm,invhelm,uu,amax,amin,
     $   adis,bpar,bstep,copt_count,box_tol,ifdiag,nb)
      else   
         call exitti('incorrect isolve specified...$',isolve)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine evalc2(cu,cm,cl,uu,tt)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real cu(nb)
      real uu(0:nb)
      real tt(0:nb)
      real ucft(0:nb)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      real cm(ic1:ic2,jc1:jc2)

      common /scrc/ work(max(lub,ltb))

      integer icalld
      save    icalld
      data    icalld /0/

      if (icalld.eq.0) then
         evalc_time=0.
         icalld=1
      endif

      stime=dnekclock()

      call rzero(cu,nb)
      if (ncloc.ne.0) then
         do k=kc1,kc2
         do j=jc1,jc2
         do i=ic1,ic2
            cu(i)=cu(i)+cl(i,j,k)*tt(j)*uu(k)
         enddo
         enddo
         enddo
         call gop(cu,work,'+  ',nb)
      endif

      call nekgsync

      evalc_time=evalc_time+dnekclock()-stime

      return
      end
c-----------------------------------------------------------------------
      subroutine evf(tt,uu,ff)

      include 'SIZE'
      include 'MOR'

      common /scrrhs/ rhs(0:lb),tmp2(0:lb),tmp3(0:lb)
      common /invevf/ buinv(lb*lb),btinv(lb*lb)
      common /screvf/ t1(0:lb),t2(0:lb),t3(0:lb),t4(0:lb)

      real uu(1),ff(1)

      call copy(t1(1),uu(1),nb)
      t1(0)=1.

      if (ifrom(1)) then
         call mxm(au0,nb+1,t1,nb+1,t2,1)

         s=-1.0/ad_re
         call cmult(t2(1),s,nb)

         call evalc2(t3(1),ctmp,cul,t1,t1)
         call sub2(t2(1),t3(1),nb)

         if (ifbuoy) then
            call copy(t3(1),uu(nb+1),nb)
            t3(0)=1.
            call mxm(but0,nb+1,t3,nb+1,t4,1)
            call add2s2(t2(1),t4(1),ad_ra,nb)
         else if (ifforce) then
            call add2(t2(1),rg(1),nb)
         endif

         call mxm(buinv,nb,t2(1),nb,ff,1)
      endif

      if (ifrom(2)) then
         call copy(t4(1),uu(nb+1),nb)
         t4(0)=1.

         call mxm(at0,nb+1,t4,nb+1,t2,1)

         s=-1.0/ad_pe
         call cmult(t2(1),s,nb)

         call evalc2(t3(1),ctmp,ctl,t1,t4)
         call sub2(t2(1),t3(1),nb)

         call mxm(btinv,nb,t2(1),nb,ff(nb+1),1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine check_diag(checkdiag,ifdiag,aa,bb,n)

      real aa(n,n),bb(n,n)
      real vv(n),tmp(n)
      integer checkdiag
      logical ifdiag

      ifdiag=.false.

      call rone(tmp,n)
      call mxm(bb,n,tmp,n,vv,1)
      checkdiag = 0

      do ii=1,n
         if (abs(vv(ii)-bb(ii,ii)).ge.1e-10) then
            checkdiag=checkdiag+1
         endif
      enddo
      if (checkdiag==0) then
         ifdiag=.true.
      endif
      call update_h(aa,bb,ifdiag)

      return
      end
c-----------------------------------------------------------------------
      subroutine check_box(chekbc,uu,amax,amin,tol_box,n)

      real uu(n)
      real amax(n),amin(n)
      integer chekbc

      chekbc=0

      do ii=1,n
         if ((uu(ii)-amax(ii)).ge.tol_box) then
            chekbc = 1
         elseif ((amin(ii)-uu(ii)).ge.tol_box) then
            chekbc = 1
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine update_h(hh,invhh,ifdiag)
         
      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      common /scrcopt/ helm(lb**2,2),invhelm(lb**2,2)

      real hh(nb**2),invhh(nb**2)
      real invhelm
      logical ifdiag

      if (ifpod(1)) then 
         if (ifdiag) then 
            do jj=1,nb
               helm(jj,1) = 1/invhh(jj+(jj-1)*nb)
            enddo
         else 
            call copy(helm(1,1),hh(1),nb*nb)
         endif
         if (ifdiag) then
            do jj=1,nb
               invhelm(jj,1) = invhh(jj+(jj-1)*nb)
            enddo
         else 
            call copy(invhelm(1,1),invhh(1),nb*nb)
         endif
      elseif (ifpod(2)) then
         if (ifdiag) then 
            do jj=1,nb
               helm(jj,2) = 1/invhh(jj+(jj-1)*nb)
            enddo
         else 
            call copy(helm(1,2),hh(1),nb*nb)
         endif
         if (ifdiag) then
            do jj=1,nb
               invhelm(jj,2) = invhh(jj+(jj-1)*nb)
            enddo
         else 
            call copy(invhelm(1,2),invhh(1),nb*nb)
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_ucoef_in_ext(s1,s2,t1)

      ! Set the coefficient of residual involving summation in
      ! extrapolation term.  Take convection term as an example and
      ! consider only EXT3 is used The coefficient can be written as:
      ! The inner product of [c_0 c_1 c_2 c_3 c_4 c_5] and
      ! [{u**2}^{J_0-2} {u**2}^_{J_0-1} {u**2}^{J_0} {u**2}^{J-2}
      ! {u**2}^{J-1} {u**2}^{J}] plus sum^J_{j=J_0-2} {u**2}^j. This
      ! subroutine computes sum^J_{j=J_0-2} {u**2}^j

      ! s1:=ua_ext
      ! s2:=u2a_ext
      ! t1:=u

      include 'SIZE'
      include 'MOR'

      real s1(0:nb),s2(0:nb,0:nb),t1(0:nb,3)
      real tmp(0:nb)

      character*3 rfilter

      logical ifbdf1,ifbdf2,ifbdf3

      ifbdf1 = navg_step.eq.1.and.ad_step.eq.navg_step+1
      ifbdf2 = navg_step.eq.2.and.ad_step.eq.navg_step
      ifbdf3 = navg_step.ge.3.and.ad_step.eq.navg_step-1

      rfilter='   '
      if (cfloc.eq.'CONV'.and.cftype.eq.'TFUN') rfilter='LER'

      if (ifbdf1.or.ifbdf2.or.ifbdf3) then
         call rzero(s1,(nb+1))
         call rzero(s2,(nb+1)**2)
         call copy(s1,t1(0,3),nb+1)
         call add2(s1,t1(0,2),nb+1)
         call add2(s1,t1(0,1),nb+1)
         if (rfilter.eq.'LER'.and.rbf.gt.0) then
            do k=1,3
            call copy(tmp,t1(0,k),nb+1)
            call pod_proj(tmp(1),rbf,nb,'step  ')
            do j=0,nb
            do i=0,nb
               s2(i,j)=s2(i,j)+t1(i,k)*tmp(j)
            enddo
            enddo
            enddo
         else
            do k=1,3
            do j=0,nb
            do i=0,nb
               s2(i,j)=s2(i,j)+t1(i,k)*t1(j,k)
            enddo
            enddo
            enddo
         endif
      else

         call add2(s1,t1(0,1),nb+1)
         if (rfilter.eq.'LER'.and.rbf.gt.0) then
            call copy(tmp,t1(0,1),nb+1)
            call pod_proj(tmp(1),rbf,nb,'step  ')
            do j=0,nb
            do i=0,nb
               s2(i,j)=s2(i,j)+t1(i,1)*tmp(j)
            enddo
            enddo
         else
            do j=0,nb
            do i=0,nb
               s2(i,j)=s2(i,j)+t1(i,1)*t1(j,1)
            enddo
            enddo
         endif
      endif

      if (ad_step.eq.ad_nsteps) then
         s=1./real(ad_nsteps-navg_step+1)
         call cmult(s1,s,nb+1)
         call cmult(s2,s,(nb+1)**2)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_tcoef_in_ext(s1,s2,t1,t2)

      ! Set the coefficient of residual involving summation in
      ! extrapolation term.  Take buoyance term as an example and
      ! consider only EXT3 is used The coefficient can be written as:
      ! The inner product of [c_0 c_1 c_2 c_3 c_4 c_5] and [T^{J_0-2}
      ! T^_{J_0-1} T^{J_0} T^{J-2} T^{J-1} T^{J}] plus sum^J_{j=J_0-2}
      ! T^j. This subroutine computes sum^J_{j=J_0-2} T^j

      ! s1:=uta_ext
      ! s2:=utua_ext
      ! t1:=u
      ! t2:=ut

      include 'SIZE'
      include 'MOR'

      real s1(0:nb),s2(0:nb,0:nb)
      real t1(0:nb,3),t2(0:nb,3)
      real tmp(0:nb)

      logical ifbdf1,ifbdf2,ifbdf3

      character*3 rfilter

      ifbdf1 = navg_step.eq.1.and.ad_step.eq.navg_step+1
      ifbdf2 = navg_step.eq.2.and.ad_step.eq.navg_step
      ifbdf3 = navg_step.ge.3.and.ad_step.eq.navg_step-1

      rfilter='   '
      if (cfloc.eq.'CONV'.and.cftype.eq.'TFUN') rfilter='LER'

      if (ifbdf1.or.ifbdf2.or.ifbdf3) then
         call rzero(s1,(nb+1))
         call rzero(s2,(nb+1)**2)
         call copy(s1,t2(0,3),nb+1)
         call add2(s1,t2(0,2),nb+1)
         call add2(s1,t2(0,1),nb+1)
         if (rfilter.eq.'LER'.and.rbf.gt.0) then
            do k=1,3
            call copy(tmp,t1(0,k),nb+1)
            call pod_proj(tmp(1),rbf,nb,'step  ')
            do j=0,nb
            do i=0,nb
               s2(i,j)=s2(i,j)+tmp(j)*t2(i,k)
            enddo
            enddo
            enddo
         else
            do k=1,3
            do j=0,nb
            do i=0,nb
               s2(i,j)=s2(i,j)+t1(j,k)*t2(i,k)
            enddo
            enddo
            enddo
         endif
      else

         call add2(s1,t2(0,1),nb+1)
         if (rfilter.eq.'LER'.and.rbf.gt.0) then
            call copy(tmp,t1(0,1),nb+1)
            call pod_proj(tmp(1),rbf,nb,'step  ')
            do j=0,nb
            do i=0,nb
               s2(i,j)=s2(i,j)+tmp(j)*t2(i,1)
            enddo
            enddo
         else
            do j=0,nb
            do i=0,nb
               s2(i,j)=s2(i,j)+t1(j,1)*t2(i,1)
            enddo
            enddo
         endif
      endif


      if (ad_step.eq.ad_nsteps) then
         s=1./real(ad_nsteps-navg_step+1)
         call cmult(s1,s,nb+1)
         call cmult(s2,s,(nb+1)**2)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setuj(s1,s2,s3,t1)

      ! set quantities in velocity residual
      ! store uj and ujfilter at different
      ! time by u and compute u2j by uj or ujfilter

      ! s1 := uj
      ! s2 := u2j
      ! s3 := ujfilter
      ! t1 := u

      include 'SIZE'
      include 'MOR'

      real s1(0:nb,8),s2(0:nb,0:nb,8)
      real s3(0:nb,8),t1(0:nb,5)

      character*3 rfilter

      rfilter='   '
      if (cfloc.eq.'CONV'.and.cftype.eq.'TFUN') rfilter='LER'

      mj=5
      nj=3

      if (ad_step.eq.(navg_step+1)) then
         do i=1,mj
            call copy(s1(0,i),t1(0,mj+1-i),nb+1)
         enddo
      else if (ad_step.eq.ad_nsteps) then
         do i=1,nj
            call copy(s1(0,mj+i),t1(0,nj+1-i),nb+1)
         enddo
         do k=1,mj+nj
            call copy(s3(0,k),s1(0,k),nb+1)
            if (rfilter.eq.'LER'.and.rbf.gt.0)
     $         call pod_proj(s3(1,k),rbf,nb,'step  ')
            call mxm(s1(0,k),nb+1,s3(0,k),1,s2(0,0,k),nb+1)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine settj(s1,s2,s3,t1,t2)

      ! set quantities in temperature residual
      ! store utj at different time by ut and compute
      ! uutj, utuj with utj and ujfilter

      ! s1 := utj
      ! s2 := uutj
      ! s3 := utuj
      ! t1 := ujfilter from setuj
      ! t2 := ut

      include 'SIZE'
      include 'MOR'

      real s1(0:nb,8),s2(0:nb,0:nb,8),s3(0:nb,0:nb,8)
      real t1(0:nb,8),t2(0:nb,5)

      mj=5
      nj=3

      if (ad_step.eq.(navg_step+1)) then
         do i=1,mj
            call copy(s1(0,i),t2(0,mj+1-i),nb+1)
         enddo
      else if (ad_step.eq.ad_nsteps) then
         do i=1,nj
            call copy(s1(0,mj+i),t2(0,nj+1-i),nb+1)
         enddo
         do k=1,mj+nj
         do j=0,nb
         do i=0,nb
            s2(i,j,k)=t1(i,k)*s1(j,k)
            s3(i,j,k)=t1(j,k)*s1(i,k)
         enddo
         enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setfj(s1,s2)

      ! set quantities in source term for residual
      ! store rqj and rqtj at different
      ! time by user specified function ft in userq

      ! s1:=rqj
      ! s2:=rqtj

      include 'SIZE'
      include 'TSTEP'
      include 'MOR'

      real s1(6),s2(6)

      logical ifbdf1,ifbdf2,ifbdf3

      ifbdf1 = navg_step.eq.1.and.ad_step.eq.navg_step+1
      ifbdf2 = navg_step.eq.2.and.ad_step.eq.navg_step
      ifbdf3 = navg_step.ge.3.and.ad_step.eq.navg_step-1

      if (ifbdf1.or.ifbdf2.or.ifbdf3) then
         s1(1) = 1
         s1(2) = 1
         s1(3) = 1
         time = (ad_step-2)*ad_dt
         call userq(1,1,1,1)
         s2(1) = ft
         time = (ad_step-1)*ad_dt
         call userq(1,1,1,1)
         s2(2) = ft
         time = (ad_step)*ad_dt
         call userq(1,1,1,1)
         s2(3) = ft
      endif

      if (ad_step.eq.ad_nsteps) then
         s1(4) = 1
         s1(5) = 1
         s1(6) = 1
         time = (ad_step-2)*ad_dt
         call userq(1,1,1,1)
         s2(4) = ft
         time = (ad_step-1)*ad_dt
         call userq(1,1,1,1)
         s2(5) = ft
         time = (ad_step)*ad_dt
         call userq(1,1,1,1)
         s2(6) = ft
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_favg_in_ext(s1,s2)

      ! Set the coefficient of residual involving summation in
      ! extrapolating source term.

      ! s1:=rqa
      ! s2:=rqta

      include 'SIZE'
      include 'TSTEP'
      include 'MOR'

      real s1,s2

      logical ifbdf1,ifbdf2,ifbdf3

      ifbdf1 = navg_step.eq.1.and.ad_step.eq.navg_step+1
      ifbdf2 = navg_step.eq.2.and.ad_step.eq.navg_step
      ifbdf3 = navg_step.ge.3.and.ad_step.eq.navg_step-1

      if (ifbdf1.or.ifbdf2.or.ifbdf3) then
         s1 = 0.0
         s2 = 0.0
         s1 = 3
         time = (ad_step-2)*ad_dt
         call userq(1,1,1,1)
         s2 = s2 + ft
         time = (ad_step-1)*ad_dt
         call userq(1,1,1,1)
         s2 = s2 + ft
         time = (ad_step)*ad_dt
         call userq(1,1,1,1)
         s2 = s2 + ft
      else
         s1 = s1 + 1.
         time = (ad_step)*ad_dt
         call userq(1,1,1,1)
         s2 = s2 + ft
      endif

      if (ad_step.eq.ad_nsteps) then
         s=1./real(ad_nsteps-navg_step+1)
         s1 = s1*s
         s2 = s2*s
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_time_avg_resid_coef

      ! currently only support for BDF3/EXT3
      ! compute the coefficient in the time-averaged residual

      include 'SIZE'
      include 'MOR'
      if (ifrom(1)) call set_time_avg_resid_coef_u
      if (ifrom(2)) call set_time_avg_resid_coef_t

      return
      end
c-----------------------------------------------------------------------
      subroutine set_time_avg_resid_coef_u

      ! currently only support for BDF3/EXT3
      ! compute the coefficient in the time-averaged velocity residual

      include 'SIZE'
      include 'MOR'

      call setuj(uj,u2j,ujfilter,u)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_time_avg_resid_coef_t

      ! currently only support for BDF3/EXT3
      ! compute the coefficient in the time-averaged temperature
      ! residual

      include 'SIZE'
      include 'MOR'

      call settj(utj,uutj,utuj,ujfilter,ut)

      if (ifsource) then
         call set_favg_in_ext(rqa,rqta)
         call setfj(rqj,rqtj)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine addeddy(flu,cedm,cedvec,cl,uu)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real uu(0:nb), flu(nb,nb)
      real cedm(nb,nb),cedvec(nb)

      real tmp(nb,0:nb)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      real cm(ic1:ic2,jc1:jc2)

      call rzero(tmp,nb*(nb+1))

      call mxm(cl,(ic2-ic1+1)*(jc2-jc1+1),uu(kc1),(kc2-kc1+1),cm,1)

      do j=jc1,jc2
      do i=ic1,ic2
      tmp(i,j) = tmp(i,j)+cm(i,j)
      enddo
      enddo
      call gop(tmp,stmp,'+  ',nb*(nb+1))
      call copy(cedm,tmp(1,1),nb*nb)
      call copy(cedvec,tmp(1,0),nb)

      call add2(flu,cedm,nb*nb)

      return
      end
