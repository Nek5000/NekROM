c-----------------------------------------------------------------------
      subroutine cres

      include 'SIZE'
      include 'MOR'

      common /ccres/ cdiff(0:lb)

      parameter (lt=lx1*ly1*lz1*lelt)

      if (eqn.eq.'POI') then
         call set_theta_poisson
      else if (eqn.eq.'HEA') then
         call set_theta_heat
      else if (eqn.eq.'ADE') then
         call set_theta_ad
      else if (eqn.eq.'NSE') then
         call set_theta_ns
      endif

      res=0.

      do j=1,nres
      do i=1,nres
         res=res+sigma(i,j)*theta(i)*theta(j)
      enddo
      enddo

      res=sqrt(res)

      eierr=0.

      call sub3(cdiff,ua,uas,nb+1)
      call mxm(bu0,nb+1,cdiff,nb+1,ctmp,1)

      do i=0,nb
      do j=0,nb
         eierr=eierr+bu0(1+i+(nb+1)*j)*(ua(i)-uas(i))*(ua(j)-uas(j))
      enddo
      enddo

      eierr=sqrt(eierr)

      if (res.le.0) call exitti('negative semidefinite residual$',n)

      if (nid.eq.0) write (6,*) 'res:',res
      if (nid.eq.0) write (6,*) 'err:',eierr

      return
      end
c-----------------------------------------------------------------------
      subroutine set_xi_poisson

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrxi/ wk(lt)

      n=lx1*ly1*lz1*nelv

      l=1
      if (ifield.eq.1) then
        call exitti('(set_xi_poisson) ifield.eq.1 not supported...$',nb)
      else
         if (ips.eq.'L2 ') then
            do i=1,nb
               call axhelm(xi(1,l),tb(1,i),ones,zeros,1,1)
               call binv1(xi(1,l))
               l=l+1
            enddo
            call copy(xi(1,l),qq,n)
            call binv1(xi(1,l))
            l=l+1
            do i=1,nb
               call set_gradn(wk,tb(1,i))
               call set_surf(xi(1,l),wk)
               call binv1_nom(xi(1,l))
               l=l+1
            enddo
         else
            call exitti('(set_xi_poisson) ips!=L2 not supported...$',nb)
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_xi_heat

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      n=lx1*ly1*lz1*nelv

      l=1
      if (ifield.eq.1) then
         call exitti('(set_xi_heat) ifield.eq.1 not supported...$',nb)
      else
         if (ips.eq.'L2 ') then
            do i=0,nb
               call copy(xi(1,l),tb(1,i),n)
               call col2(xi(1,l),bm1,n)
               call binv1(xi(1,l))
               l=l+1
            enddo
            do i=0,nb
               call axhelm(xi(1,l),tb(1,i),ones,zeros,1,1)
               call binv1(xi(1,l))
               l=l+1
            enddo
            call copy(xi(1,l),qq,n)
            call binv1(xi(1,l))
            l=l+1
         else
            call exitti('(set_xi_heat) ips != L2 not supported...$',ips)
         endif
      endif

      if ((l-1).gt.nres) then
         call exitti('increase nres$',l-1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_xi_ad

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      n=lx1*ly1*lz1*nelv

      l=1
      if (ifield.eq.1) then
         call exitti('(set_xi_ad) ifield.eq.1 not supported...$',nb)
      else
         if (ips.eq.'L2 ') then
            do i=0,nb
               call copy(xi(1,l),tb(1,i),n)
               l=l+1
            enddo
            call push_op(vx,vy,vz)
            if (ifrom(1)) then
               do j=0,nb
                  call opcopy(vx,vy,vz,ub(1,j),vb(1,j),wb(1,j))
                  do i=0,nb
                     if (ifaxis) then
                        call conv1d(xi(1,l),tb(1,i))
                     else
                        call convect_new(xi(1,l),tb(1,i),.false.,
     $                                  ub(1,j),vb(1,j),wb(1,j),.false.)
                        call invcol2(wk1,bm1,n)  ! local mass inverse
                     endif
                     l=l+1
                  enddo
               enddo
            else
               call opcopy(vx,vy,vz,ub,vb,wb)
               do i=0,nb
                  call convop(xi(1,l),tb(1,i))
                  l=l+1
               enddo
            endif
            call pop_op(vx,vy,vz)
            do i=0,nb
               call axhelm(xi(1,l),tb(1,i),ones,zeros,1,1)
               call binv1(xi(1,l))
               l=l+1
            enddo
         else
            call exitti('(set_xi_ad) ips != L2 not supported...$',ips)
         endif
      endif

      if ((l-1).gt.nres) then
         call exitti('increase nres$',l-1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_xi_ns

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /screi/ wk1(lt),wk2(lt),wk3(lt),wk4(lt),wk5(lt)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside set_xi_ns'

      l=1
      if (ifield.eq.1) then
         if (ips.eq.'L2 ') then
            do i=0,nb
c              call opcopy(xi_u(1,1,l),xi_u(1,2,l),xi_u(1,ldim,l),
c    $                     ub(1,i),vb(1,i),wb(1,i))
               call comp_vort3(xi_u(1,1,l),wk1,wk2,
     $                         ub(1,i),vb(1,i),wb(1,i))
c              call outpost(xi_u(1,1,l),wk1,wk2,pr,t,'xib')
               l=l+1
            enddo
            call push_op(vx,vy,vz)
            do j=0,nb
               call opcopy(vx,vy,vz,ub(1,j),vb(1,j),wb(1,j))
               do i=0,nb
                  if (ifaxis) then
                     call conv1d(wk1,ub(1,i))
                     call conv1d(wk2,vb(1,i))
                     call rzero(wk3,n)
                     call comp_vort3(xi_u(1,1,l),wk4,wk5,wk1,wk2,wk3)
c                    call outpost(xi_u(1,1,l),wk1,wk2,pr,t,'xic')
                     l=l+1
                  else
                     call convect_new(wk1,ub(1,i),.false.,
     $                                ub(1,j),vb(1,j),wb(1,j),.false.)
                     call convect_new(wk2,vb(1,i),.false.,
     $                                ub(1,j),vb(1,j),wb(1,j),.false.)
                     call invcol2(wk1,bm1,n)  ! local mass inverse
                     call invcol2(wk2,bm1,n)  ! local mass inverse
                     call comp_vort3(xi_u(1,1,l),wk4,wk5,wk1,wk2,wk3)
c                    call outpost(xi_u(1,1,l),wk1,wk2,pr,t,'xic')
                     l=l+1
                  endif
               enddo
            enddo
            call pop_op(vx,vy,vz)
            do i=0,nb ! todo investigate possible source of error for ifaxis
               call copy(wk1,xi_u(1,1,i+1),n)
               call axhelm(xi_u(1,1,l),wk1,ones,zeros,1,1)
               call binv1(xi_u(1,1,l))
c              call outpost(xi_u(1,1,l),wk1,wk2,pr,t,'xia')
               l=l+1
            enddo
            if (ifbuoy) then
               do i=0,nb
                  call opcopy(wk1,wk2,wk3,gx,gy,gz)
                  call opcolv(wk1,wk2,wk3,tb(1,i))
                  call invcol2(wk1,bm1,n)
                  call invcol2(wk2,bm1,n)
                  if (ldim.eq.3) call invcol2(wk3,bm1,n)
                  call comp_vort3(xi_u(1,1,l),wk4,wk5,wk1,wk2,wk3)
c                 call outpost(xi_u(1,1,l),wk1,wk2,pr,t,'xig')
                  l=l+1
               enddo
            endif
         else
            call exitti('(set_xi_ns) ips != L2 not supported...$',ips)
         endif
      else
         call exitti('(set_xi_ns) ifield.ne.1 not supported...$',nb)
      endif

      if ((l-1).gt.nres) then
         call exitti('increase nres$',l-1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_sigma

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      n=lx1*ly1*lz1*nelv

      if (eqn.eq.'POI') then
         nres=nb+1
         nres=nb+1+nb
      else if (eqn.eq.'HEA') then
         nres=(nb+1)*2+1
      else if (eqn.eq.'ADE') then
         nres=(nb+1)*3
         if (ifrom(1)) nres=(nb+1)*2 + (nb+1)**2
      else if (eqn.eq.'NSE') then
         nres=(nb+1)*2+(nb+1)**2
         if (ifbuoy) nres=nres+nb+1
      endif

      if (nres.gt.lres) call exitti('nres > lres$',nres)
      if (nres.le.0) call exitti('nres <= 0$',nres)

      if (rmode.eq.'ON '.or.rmode.eq.'ONB') then
         call read_serial(sigtmp,nres*nres,'ops/sigma ',sigma,nid)
         l=1
         do j=1,nres
         do i=1,nres
            sigma(i,j)=sigtmp(l,1)
            l=l+1
         enddo
         enddo
      else
         if (eqn.eq.'POI') then
            ifield=2
            call set_xi_poisson
         else if (eqn.eq.'HEA') then
            ifield=2
            call set_xi_heat
         else if (eqn.eq.'ADE') then
            ifield=2
            call set_xi_ad
         else if (eqn.eq.'NSE') then
            ifield=1
            call set_xi_ns
         endif

         if (ifield.eq.2) then
            do i=1,nres
            do j=1,nres
               sigma(i,j)=glsc3(xi(1,i),xi(1,j),bm1,n)
            enddo
            enddo
         else if (ifield.eq.1) then
            if (.true.) then ! if using voritcity residual
               do i=1,nres
                  do j=1,nres
                     sigma(i,j)=glsc3(xi_u(1,1,i),xi_u(1,1,j),bm1,n)
                  enddo
                  write (6,*) i,sigma(1,i),'sigma'
               enddo
            else
               do i=1,nres
               do j=1,nres
                  sigma(i,j)=
     $               op_glsc2_wt(xi_u(1,1,i),xi_u(1,2,i),xi_u(1,ldim,i),
     $                       xi_u(1,1,j),xi_u(1,2,j),xi_u(1,ldim,j),bm1)
               enddo
               enddo
            endif
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_theta_poisson

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      n=lx1*ly1*lz1*nelv

      l=1
      do i=1,nb
         theta(l)=ut(i)
         l=l+1
      enddo

      theta(l)=-1.
      l=l+1

      do i=1,nb
         theta(l)=-ut(i)
         theta(l)=0.
         l=l+1
      enddo

      do i=1,l-1
         write (6,*) i,theta(i),'theta'
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_theta_heat

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      n=lx1*ly1*lz1*nelv

      l=1
      call set_betaj
      call mxm(utj,nb+1,betaj,6,theta(l),1)

      l=l+nb+1
      do i=0,nb
         theta(l)=uta(i)
         l=l+1
      enddo

      theta(l)=-1.

      l=l+1

      return
      end
c-----------------------------------------------------------------------
      subroutine set_theta_ad

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      n=lx1*ly1*lz1*nelv

      l=1
      call set_betaj
      call mxm(utj,nb+1,betaj,6,theta(l),1)

      l=l+nb+1

      call set_alphaj

      if (ifrom(1)) then
c        call mxm(uutj,(nb+1)**2,alphaj,6,theta(l),1)
         call mxm(utuj,(nb+1)**2,alphaj,6,theta(l),1)
      else
         call mxm(utj,nb+1,alphaj,6,theta(l),1)
      endif

      if (ifrom(1)) then
         do j=0,nb
         do i=0,nb
            theta(l)=theta(l)+utua(i+(nb+1)*j)
            l=l+1
         enddo
         enddo
      else
         do i=0,nb
            theta(l)=theta(l)+uta(i)
            l=l+1
         enddo
      endif

      do i=0,nb
         theta(l)=param(8)*uta(i)
         l=l+1
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_theta_ns

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      n=lx1*ly1*lz1*nelv

      l=1
      call set_betaj
      call mxm(utj,nb+1,betaj,6,theta(l),1)

      l=l+nb+1

      call set_alphaj(alphaj)
      call mxm(ut2j,(nb+1)**2,alphaj,6,theta(l),1)
      do j=0,nb
      do i=0,nb
         theta(l)=theta(l)+u2a(1+i+(nb+1)*j)
         l=l+1
      enddo
      enddo

      do i=0,nb
         theta(l)=param(2)*ua(i)
         l=l+1
      enddo

      if (ifbuoy) then
         do i=0,nb
            theta(l)=ad_ra*uta(i)
            l=l+1
         enddo
      endif

      do i=1,nres
         if (nio.eq.0) write (6,*) theta(i),'theta'
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_poisson

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

      eqn='POI'

      ifflow=.false.
      ifheat=.true.

      param(174)=-1.

      call rom_init_params
      call rom_init_fields

      call setgram
      call setevec

      call setbases

      call setops
      call dump_all

      call set_sigma

      rhs(0)=1.
      call setr_poisson(rhs(1),icount)

      call add2sxy(flut,0.,at,1./ad_pe,nb*nb)
      call lu(flut,nb,nb,irt,ict)

      call solve(rhs(1),flut,1,nb,nb,irt,ict)

      call recont(t,rhs)
      call copy(ut,rhs,nb+1)

      call cres

      step_time=step_time+dnekclock()-last_time

      return
      end
c-----------------------------------------------------------------------
      subroutine setr_poisson(rhs)

      include 'SIZE'
      include 'MOR'

      real rhs(nb)

      n=lx1*ly1*lz1*nelv

      do i=1,nb
c        rhs(i)=wl2sip(qq,tb(1,i))
         rhs(i)=glsc2(qq,tb(1,i),n)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_alphaj

      include 'SIZE'
      include 'MOR'

      ! ad_alpha(3,3)

      alphaj(1)=ad_alpha(1,1)+ad_alpha(2,2)+ad_alpha(3,3)
      alphaj(2)=ad_alpha(1,2)-ad_alpha(1,3)
      alphaj(3)=0.
      alphaj(4)=-ad_alpha(3,3)
      alphaj(5)=-ad_alpha(2,3)-ad_alpha(3,3)
      alphaj(6)=0.

      call cmult(alphaj,1./(1.*ad_nsteps),6)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_betaj

      include 'SIZE'
      include 'MOR'

      ! ad_beta(4,3)

      betaj(1)=ad_beta(1+1,1)+ad_beta(2+1,2)+ad_beta(3+1,3)
      betaj(2)=ad_beta(0+1,1)+ad_beta(1+1,2)+ad_beta(2+1,3)
     $        +ad_beta(3+1,3)
      betaj(3)=ad_beta(0+1,2)+ad_beta(1+1,3)+ad_beta(2+1,3)
     $        +ad_beta(3+1,3)

      betaj(4)=ad_beta(0+1,3)+ad_beta(1+1,3)+ad_beta(2+1,3)
      betaj(5)=ad_beta(0+1,3)+ad_beta(1+1,3)
      betaj(6)=ad_beta(0+1,3)

      call cmult(betaj,1./(ad_dt*ad_nsteps),6)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_rhs

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /screi/ wk1(lt),wk2(lt),wk3(lt),wk4(lt),wk5(lt)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside set_rhs'
      
      l1=1
      do i=0,nb
         ! setup rhs for velocity representator
         ifield=1
         call opcopy(wk1,wk2,wk3,ub(1,i),vb(1,i),wb(1,i))
         call axhelm(riesz_ru(1,1,l1),wk1,ones,zeros,1,1)
         call axhelm(riesz_ru(1,2,l1),wk2,ones,zeros,1,1)
         if (ldim.eq.3) then
            call axhelm(riesz_ru(1,ldim,l1),wk1,ones,zeros,1,1)
         endif
         call cmult(riesz_ru(1,1,l1),param(2),n)
         call cmult(riesz_ru(1,2,l1),param(2),n)
         if (ldim.eq.3) call cmult(riesz_ru(1,ldim,l1),param(2),n)
         call opchsgn(riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $                riesz_ru(1,ldim,l1))
         l1=l1+1
      enddo
      write(6,*)l1,'lres_u_1'

      l2=1
      do i=0,nb
         ! setup rhs for temperature representator
         ifield=2
         call copy(wk4,tb(1,i),n)
         call axhelm(riesz_rt(1,l2),wk4,ones,zeros,1,1)
         call cmult(riesz_rt(1,l2),param(8),n)
         call chsign(riesz_rt(1,l2),n) 
         l2=l2+1
      enddo
      write(6,*)l2,'lres_t_1'

      do i=0,nb
         ifield=1
         call opcopy(riesz_ru(1,1,l1),riesz_ru(1,2,l1)
     $               ,riesz_ru(1,ldim,l1)
     $               ,tb(1,i),zeros,zeros)
         l1=l1+1
      enddo
      write(6,*)l1,'lres_u_2'

      do i=0,nb
         ifield=1
         call opcopy(riesz_ru(1,1,l1),riesz_ru(1,2,l1)
     $               ,riesz_ru(1,ldim,l1)
     $               ,zeros,tb(1,i),zeros)
         l1=l1+1
      enddo
      write(6,*)l1,'lres_u_3'

      do j=0,nb
         do i=0,nb
         ifield=1
            call convect_new(riesz_ru(1,1,l1),ub(1,i),.false.,
     $                       ub(1,j),vb(1,j),wb(1,j),.false.)
            call convect_new(riesz_ru(1,2,l1),vb(1,i),.false.,
     $                       ub(1,j),vb(1,j),wb(1,j),.false.)
            if (ldim.eq.3) then
               call convect_new(riesz_ru(1,ldim,l1),wb(1,i),.false.,
     $                          ub(1,j),vb(1,j),wb(1,j),.false.)
            endif
            call opchsgn(riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $                   riesz_ru(1,ldim,l1))

         ifield=2
            call convect_new(riesz_rt(1,l2),tb(1,i),.false.,
     $                       ub(1,j),vb(1,j),wb(1,j),.false.)
            call chsign(riesz_rt(1,l2),n) 
            l1=l1+1
            l2=l2+1
         enddo
      enddo
      write(6,*)l1,'lres_u_4'
      write(6,*)l2,'lres_t_2'

      if (nio.eq.0) write (6,*) 'exit set_rhs'

      return
      end
c-----------------------------------------------------------------------
      subroutine set_rr_ns

c     Compute the riesz representators for steady Boussinesq 
c     incompressible NS (quadratically nonlinear elliptic problem)
c     L2 might not be correct......
c     Currently we consider H1 norm
c     Might be able to use H1 semi-norm depending on the 
c     boundary condition of the problem

      include 'SIZE'
      include 'TOTAL'
      include 'ORTHOT'
      include 'CTIMER'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /screi/ wk1(lt),wk2(lt),wk3(lt),wk4(lt),wk5(lt)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside set_rr_ns'

      call rone(ones,n)
      call rzero(zeros,n)

      l1=1
      l2=1
      do i=0,nb
         ifield=1
         call ophinv(xi_u(1,1,l1),xi_u(1,2,l1),xi_u(1,ldim,l1),
     $               riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $               riesz_ru(1,ldim,l1),
     $               ones,ones,tolhv,nmaxv)               
         write(6,*)'riesz_u',l1,'completed'
         ifield=2
         call hsolve  ('TEMP',xi_t(1,l2),riesz_rt(1,l2),ones,ones
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxt(ifield-1),1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
         write(6,*)'riesz_t',l2,'completed'
         l1=l1+1
         l2=l2+1
      enddo
      write(6,*)l1,'lres_u_1'
      write(6,*)l2,'lres_t_1'

      do i=0,nb
         ifield=1
         call ophinv(xi_u(1,1,l1),xi_u(1,2,l1),xi_u(1,ldim,l1),
     $               riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $               riesz_ru(1,ldim,l1),
     $               ones,ones,tolhv,nmaxv)               
         write(6,*)'riesz_u',l1,'completed'
         l1=l1+1
      enddo
      write(6,*)l1,'lres_u_2'

      do i=0,nb
         ifield=1
         call ophinv(xi_u(1,1,l1),xi_u(1,2,l1),xi_u(1,ldim,l1),
     $               riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $               riesz_ru(1,ldim,l1),
     $               ones,ones,tolhv,nmaxv)               
         write(6,*)'riesz_u',l1,'completed'
         l1=l1+1
      enddo
      write(6,*)l1,'lres_u_3'

      do j=0,nb
         do i=0,nb
            ifield=1
            call ophinv(xi_u(1,1,l1),xi_u(1,2,l1),xi_u(1,ldim,l1),
     $               riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $               riesz_ru(1,ldim,l1),
     $               ones,ones,tolhv,nmaxv)               
            write(6,*)'riesz_u',l1,'completed'
            ifield=2
            call hsolve  ('TEMP',xi_t(1,l2),riesz_rt(1,l2),ones,ones
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxt(ifield-1),1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
            l1=l1+1
            l2=l2+1
            write(6,*)'riesz_t',l2,'completed'
         enddo
      enddo
      write(6,*)l1,'lres_u_4'
      write(6,*)l2,'lres_t_2'

      if ((l1-1).gt.nres_u) then
         call exitti('increase nres_u$',l1-1)
      endif
      if ((l2-1).gt.nres_t) then
         call exitti('increase nres_t$',l2-1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_sigma_new

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      n=lx1*ly1*lz1*nelv
      
      l1=1
      l2=1
      do i=1,nres_u
         do j=1,nres_u
            sigma_u(i,j) = h10vip(riesz_ru(1,1,i),riesz_ru(1,2,i),
     $                          riesz_ru(1,ldim,i),riesz_ru(1,1,j),
     $                          riesz_ru(1,2,j),riesz_ru(1,ldim,j))
            sigma_u(i,j) = sigma_u(i,j) + wl2vip(riesz_ru(1,1,i),
     $                          riesz_ru(1,2,i),
     $                          riesz_ru(1,ldim,i),riesz_ru(1,1,j),
     $                          riesz_ru(1,2,j),riesz_ru(1,ldim,j))
         enddo
            write (6,*) i,sigma_u(1,i),'sigma_u'
      enddo

      do i=1,nres_t
         do j=1,nres_t
            sigma_t(i,j) = h10sip(riesz_rt(1,i),riesz_rt(1,j))
            sigma_t(i,j) = sigma_t(i,j) + 
     $                   wl2sip(riesz_rt(1,i),riesz_rt(1,j))
         enddo
            write (6,*) i,sigma_t(1,i),'sigma_t'
      enddo
      call dump_serial(sigma_u,(nres_u)**2,'ops/sigma_u',nid)
      call dump_serial(sigma_t,(nres_t)**2,'ops/sigma_t',nid)

      return
      end
