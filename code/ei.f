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
      subroutine c_avei

      include 'SIZE'
      include 'MOR'

      common /ccres/ cdiff(0:lb)

      parameter (lt=lx1*ly1*lz1*lelt)

      call set_theta_ns

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
      else if (eqn.eq.'SNSE') then
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
         else if (eqn.eq.'SNSE') then
            call set_rr_ns
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

      call cmult(alphaj,1./(1.*(ad_nsteps-navg_step+1)),6)

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

      call cmult(betaj,1./(ad_dt*(ad_nsteps-navg_step+1)),6)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_residual

c     Compute ROM residual for steady NS with Boussinesq 

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /screi/ wk1(lt),wk2(lt),wk3(lt),wk4(lt),wk5(lt),wk6(lt)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside set_residual'

      call rone(ones,n)
      call rzero(zeros,n)
      
      l1=1
      do i=0,nb
         ! setup rhs for velocity representator
         ifield=1
         call opcopy(wk1,wk2,wk3,ub(1,i),vb(1,i),wb(1,i))
         call axhelm(riesz_ru(1,1,l1),wk1,ones,zeros,1,1)
         call axhelm(riesz_ru(1,2,l1),wk2,ones,zeros,1,2)
         if (ldim.eq.3) then
            call axhelm(riesz_ru(1,ldim,l1),wk3,ones,zeros,1,3)
         endif
         call opcmult(riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $                riesz_ru(1,ldim,l1),param(2))
         call opchsgn(riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $                riesz_ru(1,ldim,l1))
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_1'

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
      if (nid.eq.0) write(6,*)l2,'lres_t_1'

      do i=0,nb
         ifield=1
         call opcopy(riesz_ru(1,1,l1),riesz_ru(1,2,l1)
     $               ,riesz_ru(1,ldim,l1)
     $               ,tb(1,i),zeros,zeros)
         call opcolv(riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $               riesz_ru(1,ldim,l1),bm1)
         call opchsgn(riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $                riesz_ru(1,ldim,l1))
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_2'

      do i=0,nb
         ifield=1
         call opcopy(riesz_ru(1,1,l1),riesz_ru(1,2,l1)
     $               ,riesz_ru(1,ldim,l1)
     $               ,zeros,tb(1,i),zeros)
         call opcolv(riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $               riesz_ru(1,ldim,l1),bm1)
         call opchsgn(riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $                riesz_ru(1,ldim,l1))
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_3'

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
      if (nid.eq.0) write(6,*)l1,'lres_u_4'
      if (nid.eq.0) write(6,*)l2,'lres_t_2'

      ! add pressure residual
      ! Uncomment this part when the riesz is not divergence free
c     l3=1
c     call ifield=1
c     do i=1,nb
c     call opgradt(riesz_rp(1,1,l3),riesz_rp(1,2,l3),
c    $            riesz_rp(1,ldim,l3),pb(1,i))
c     l3=l3+1
c     enddo

      if (nio.eq.0) write (6,*) 'exit set_residual'

      return
      end
c-----------------------------------------------------------------------
      subroutine set_rr_ns

c     Compute the nondivergence free riesz representators for steady 
c     Boussinesq incompressible NS 
c     (quadratically nonlinear elliptic problem)
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
         if (nid.eq.0) write(6,*)'riesz_u',l1,'completed'
         ifield=2
         ifld1 = ifield-1
         napproxt(1,ifld1) = laxtt
c        call hsolve  ('TEMP',xi_t(1,l2),riesz_rt(1,l2),ones,ones
c    $                   ,tmask(1,1,1,1,ifield-1)
c    $                   ,tmult(1,1,1,1,ifield-1)
c    $                   ,imesh,tolht(ifield),nmxt(ifield-1),1
c    $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
         call hsolve  ('TEMP',xi_t(1,l2),riesz_rt(1,l2),ones,ones
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxh,1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
         if (nid.eq.0) write(6,*)'riesz_t',l2,'completed'
         l1=l1+1
         l2=l2+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_1'
      if (nid.eq.0) write(6,*)l2,'lres_t_1'

      do i=0,nb
         ifield=1
         call ophinv(xi_u(1,1,l1),xi_u(1,2,l1),xi_u(1,ldim,l1),
     $               riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $               riesz_ru(1,ldim,l1),
     $               ones,ones,tolhv,nmaxv)               
         if (nid.eq.0) write(6,*)'riesz_u',l1,'completed'
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_2'

      do i=0,nb
         ifield=1
         call ophinv(xi_u(1,1,l1),xi_u(1,2,l1),xi_u(1,ldim,l1),
     $               riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $               riesz_ru(1,ldim,l1),
     $               ones,ones,tolhv,nmaxv)               
         if (nid.eq.0) write(6,*)'riesz_u',l1,'completed'
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_3'

      do j=0,nb
         do i=0,nb
            ifield=1
            call ophinv(xi_u(1,1,l1),xi_u(1,2,l1),xi_u(1,ldim,l1),
     $               riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $               riesz_ru(1,ldim,l1),
     $               ones,ones,tolhv,nmaxv)               
            if (nid.eq.0) write(6,*)'riesz_u',l1,'completed'
            ifield=2
            ifld1 = ifield-1
            napproxt(1,ifld1) = laxtt
c           call hsolve  ('TEMP',xi_t(1,l2),riesz_rt(1,l2),ones,ones
c    $                   ,tmask(1,1,1,1,ifield-1)
c    $                   ,tmult(1,1,1,1,ifield-1)
c    $                   ,imesh,tolht(ifield),nmxt(ifield-1),1
c    $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
            call hsolve  ('TEMP',xi_t(1,l2),riesz_rt(1,l2),ones,ones
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxh,1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
            l1=l1+1
            l2=l2+1
            if (nid.eq.0) write(6,*)'riesz_t',l2,'completed'
         enddo
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_4'
      if (nid.eq.0) write(6,*)l2,'lres_t_2'
      ! need to allocate riesz_rp 
c     ifield=1
c     call ophinv(xi_u(1,1,l1),xi_u(1,2,l1),xi_u(1,ldim,l1),
c    $               riesz_rp(1,1,1),riesz_rp(1,2,1),
c    $               riesz_rp(1,ldim,1),
c    $               ones,ones,tolhv,nmaxv)               
c     l1=l1+1
      if (nid.eq.0) write(6,*)'riesz_u',l1,'completed'

      if ((l1-1).gt.nres_u) then
         call exitti('increase nres_u$',l1-1)
      endif
      if ((l2-1).gt.nres_t) then
         call exitti('increase nres_t$',l2-1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine csigma_u(aa)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real aa(1:nres_u,1:nres_u)

      n=lx1*ly1*lz1*nelv

      call nekgsync
      csigma_u_start=dnekclock()

      do i=1,nres_u
         do j=1,nres_u
            aa(i,j) = h10vip(xi_u(1,1,i),xi_u(1,2,i),
     $                          xi_u(1,ldim,i),xi_u(1,1,j),
     $                          xi_u(1,2,j),xi_u(1,ldim,j))
            aa(i,j) = aa(i,j) + wl2vip(xi_u(1,1,i),
     $                          xi_u(1,2,i),
     $                          xi_u(1,ldim,i),xi_u(1,1,j),
     $                          xi_u(1,2,j),xi_u(1,ldim,j))
         enddo
         if (nid.eq.0) write (6,*) i,aa(1,i),'sigma_u'
      enddo

      call nekgsync
      if (nio.eq.0) write (6,*) 'csigma_u_time:',
     $             dnekclock()-csigma_u_start
      return
      end
c-----------------------------------------------------------------------
      subroutine csigma_t(aa)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real aa(1:nres_t,1:nres_t)

      n=lx1*ly1*lz1*nelv
      
      call nekgsync
      csigma_t_start=dnekclock()

      do i=1,nres_t
         do j=1,nres_t
            aa(i,j) = h10sip(xi_t(1,i),xi_t(1,j))
            aa(i,j) = aa(i,j) + 
     $                wl2sip(xi_t(1,i),xi_t(1,j))
         enddo
         if (nid.eq.0) write (6,*) i,aa(1,i),'sigma_t'
      enddo

      call nekgsync
      if (nio.eq.0) write (6,*) 'csigma_t_time:',
     $             dnekclock()-csigma_t_start
      return
      end
c-----------------------------------------------------------------------
      subroutine set_sigma_new

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      real wk1(lt)
      logical ifsteady, ifdebug

      n=lx1*ly1*lz1*nelv
      
      if (nio.eq.0) write (6,*) 'inside set_sigma_new'

      call nekgsync
      sigma_time=dnekclock()

      ifsteady = .true.
      ifdebug = .false.
      ifsteady = .false.

      if (ifsteady) then
         ! nres_u and nres_t for steady NS + energy equations
         nres_u=(nb+1)*3+(nb+1)**2
         nres_t=(nb+1)+(nb+1)**2
      else
         ! nres_u and nres_t for unsteady NS + energy equations
         nres_u=(nb+1)*4+(nb+1)**2
         nres_t=(nb+1)*2+(nb+1)**2
      endif

      write(6,*)'nres_u: ',nres_u
      write(6,*)'nres_t: ',nres_t

      if (nint(param(175)).ne.2) then
         if (nres_u.gt.lres_u) call exitti('nres_u > lres_u$',nres_u)
         if (nres_u.le.0) call exitti('nres_u <= 0$',nres_u)
         if (nres_t.gt.lres_t) call exitti('nres_t > lres_t$',nres_t)
         if (nres_t.le.0) call exitti('nres_t <= 0$',nres_t)
      endif

      if (rmode.eq.'ON '.or.rmode.eq.'ONB'.and.nint(param(175)).ne.2) 
     $   then

         if (nio.eq.0) write (6,*) 'reading sigma_u...'
         call read_sigma_u_serial(sigma_u,nres_u,nres_u,'ops/sigma_u ',
     $                        lres_u,nres_u,(lb+1),(nb+1),wk1,nid)
         if (nio.eq.0) write (6,*) 'reading sigma_t...'
         call read_sigma_t_serial(sigma_t,nres_t,nres_t,'ops/sigma_t ',
     $                        lres_t,nres_t,(lb+1),(nb+1),wk1,nid)

      else if (nint(param(175)).eq.2) then
         goto 199
      else

         if (ifdebug) then
            num_ts = 10 ! need to specify number of test samples
            call set_residual
            ! crd subroutine computes riesz representators directly
            call crd_divf(num_ts)
            call crd(num_ts)
            call crd_test(num_ts)
         else
            if (ifsteady) then
               ! for steady NS with Bossinesq + energy transport
               call set_residual
               call set_sNS_divfrr
            else
               ! for unsteady NS with Bossinesq + energy transport
               call set_residual_unsteady
               call set_uNS_divfrr
            endif
            call csigma_u(sigma_u)
            call csigma_t(sigma_t)
            call dump_serial(sigma_u,(nres_u)**2,'ops/sigma_u ',nid)
            call dump_serial(sigma_t,(nres_t)**2,'ops/sigma_t ',nid)
         endif

      endif

  199 continue
      call nekgsync

      if (nio.eq.0) write (6,*) 'sigma_time:',dnekclock()-sigma_time

      if (nio.eq.0) write (6,*) 'exiting set_sigma_new'
      return
      end
c-----------------------------------------------------------------------
      subroutine crd(num_ts)

      include 'SIZE'
      include 'TOTAL'
      include 'ORTHOT'
      include 'CTIMER'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real uu((nb+1)*2*num_ts),wk(2*(lb+1)*num_ts)
      real ts_an(num_ts)
      real dual_norm(num_ts),tmp(lt)
      integer num_ts
      logical ifexist

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside crd'

      inquire (file='./rom.dat',exist=ifexist)
      if (ifexist) call read_serial(uu,(nb+1)*2*num_ts,
     $                  './rom.dat ',wk,nid)

      inquire (file='./angle.dat',exist=ifexist)
      if (ifexist) call read_serial(ts_an,num_ts,
     $                  './angle.dat ',wk,nid)

      write(6,*)num_ts,'num_ts'
      call rzero(dual_norm,num_ts)
      do i=1,num_ts

         write(6,*)i,ts_an(i),'angle'
         call cres_new(uu(1+(i-1)*(nb+1)*2),ts_an(i))
c        do ii=1,(nb+1)*2*num_ts
c           write(6,*)ii,uu(ii),ts_an(i)
c        enddo

         call rone(ones,n)
         call rzero(zeros,n)
         ifield=1
         tolhv=1e-8
         tolht(2)=1e-8
c        nmxv=1000
         call ophinv(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
     $               res_u(1,1),res_u(1,2),res_u(1,ldim),
     $               ones,ones,tolhv,nmxv)      
c        call ophinv(eh_p(1,1),eh_p(1,2),eh_p(1,ldim),
c    $               riesz_rp(1,1,1),riesz_rp(1,2,1),riesz_rp(1,ldim,1),
c    $               ones,ones,tolhv,nmxv)      
c        call opadd2(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
c    $               eh_p(1,1),eh_p(1,2),eh_p(1,ldim))

         ifield=2
         ifld1 = ifield-1
         napproxt(1,ifld1) = laxtt
c     call hsolve  ('TEMP',eh_t,res_t,ones,ones
c    $                   ,tmask(1,1,1,1,ifield-1)
c    $                   ,tmult(1,1,1,1,ifield-1)
c    $                   ,imesh,tolht(ifield),nmxt(ifield-1),1
c    $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
         call hsolve  ('TEMP',eh_t,res_t,ones,ones
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxh,1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)

         ifield=1
         t1 = h10vip(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
     $         eh_u(1,1),eh_u(1,2),eh_u(1,ldim)) 
         t2 = wl2vip(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
     $           eh_u(1,1),eh_u(1,2),eh_u(1,ldim)) 
         ifield=2
         t3 = h10sip(eh_t,eh_t) 
         t4 = wl2sip(eh_t,eh_t)  
      
         write(6,*)t1,t2,t3,t4,'t1,t2,t3,t4'
         dual_norm(i) = t1+t2+t3+t4
         dual_norm(i) = sqrt(dual_norm(i))
         write(6,*)i,ts_an(i),'angle'
         write(6,*)i,dual_norm(i),
     $            'dual_norm H1 norm for the coupled system'
         write(6,*)i,sqrt(t1+t2),
     $            'dual_norm H1 norm for the velocity'
         write(6,*)i,sqrt(t3+t4),
     $            'dual_norm H1 norm for the temperature'

      enddo
      call dump_serial(dual_norm,num_ts,'./dual_norm ',nid)

      if (nio.eq.0) write (6,*) 'exiting crd'

      return
      end
c-----------------------------------------------------------------------
      subroutine cres_new(uu,angle)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real wk1(lt)
      real uu((nb+1)*2)
      real angle

      n=lx1*ly1*lz1*nelv
      write(6,*)nelt,nelv,lelm,'nelt,nelv'
      pi  = 4.*atan(1.)

      call rzero3(res_u(1,1),res_u(1,2),res_u(1,ldim),n)
      call rzero(res_t,n)

      l1=1
      do i=1,nb+1
         call cfill(wk1,uu(i),n)
         call add2col2(res_u(1,1),riesz_ru(1,1,l1),wk1,n)
         call add2col2(res_u(1,2),riesz_ru(1,2,l1),wk1,n)
         if (ldim.eq.3) then 
            call add2col2(res_u(1,ldim),riesz_ru(1,ldim,l1),wk1,n)
         endif
         ! add pressure contribution
c        call add2col2(res_u(1,1),riesz_rp(1,1,l1),wk1,n)
c        call add2col2(res_u(1,2),riesz_rp(1,2,l1),wk1,n)
c        if (ldim.eq.3) then 
c           call add2col2(res_u(1,ldim),riesz_rp(1,ldim,l1),wk1,n)
c        endif
         l1=l1+1 
      enddo
      write(6,*)l1,'l0'
      ! step 1
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,1),riesz_rp(1,2,1),riesz_rp(1,ldim,1))      
      ! step 2 
c     if (angle.gt.-45) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,2),riesz_rp(1,2,2),riesz_rp(1,ldim,2))
c     else
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,1),riesz_rp(1,2,1),riesz_rp(1,ldim,1))
c     endif
      ! step 3
c     if (angle.gt.-45.and.angle.lt.45) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,3),riesz_rp(1,2,3),riesz_rp(1,ldim,3))
c     elseif(angle>45) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,1),riesz_rp(1,2,1),riesz_rp(1,ldim,1))
c     elseif(angle<-45) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,2),riesz_rp(1,2,2),riesz_rp(1,ldim,2))
c     endif
      ! step 4
c     if (angle.gt.-20.and.angle.lt.45) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,4),riesz_rp(1,2,4),riesz_rp(1,ldim,4))
c     elseif(angle>45) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,2),riesz_rp(1,2,2),riesz_rp(1,ldim,2))
c     elseif(angle<-65) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,3),riesz_rp(1,2,3),riesz_rp(1,ldim,3))
c     elseif(angle.gt.-65.and.angle<-20) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,1),riesz_rp(1,2,1),riesz_rp(1,ldim,1))
c     endif
      ! step 5
c     if (angle.gt.-20.and.angle.lt.20) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,5),riesz_rp(1,2,5),riesz_rp(1,ldim,5))
c     elseif(angle>65) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,3),riesz_rp(1,2,3),riesz_rp(1,ldim,3))
c     elseif(angle<-65) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,4),riesz_rp(1,2,4),riesz_rp(1,ldim,4))
c     elseif(angle.gt.-65.and.angle<-20) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,2),riesz_rp(1,2,2),riesz_rp(1,ldim,2))
c     elseif(angle.lt.65.and.angle>20) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,1),riesz_rp(1,2,1),riesz_rp(1,ldim,1))
c     endif

      ! step 8
c     if (angle.le.22.5.and.angle.gt.-22.5) then 
c     if (angle.gt.-45) then 
c     if (angle.gt.-45.and.angle.lt.45) then 
c     if (angle>-30.and.angle.lt.-10) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,1),riesz_rp(1,2,1),riesz_rp(1,ldim,1))
c     elseif (angle.le.30.and.angle>10) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,2),riesz_rp(1,2,2),riesz_rp(1,ldim,2))
c     elseif (angle.le.-65) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,6),riesz_rp(1,2,6),riesz_rp(1,ldim,6))
c     elseif (angle>22.5.and.angle.lt.67.5) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,4),riesz_rp(1,2,4),riesz_rp(1,ldim,4))      
c     elseif (angle>67.5) then 
c        write(6,*)'angle here',angle
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,3),riesz_rp(1,2,3),riesz_rp(1,ldim,3))      
c     elseif (angle<-22.5.and.angle.gt.-62.5) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,5),riesz_rp(1,2,5),riesz_rp(1,ldim,5))
c     elseif (angle<-62.5.and.angle.gt.-85) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,6),riesz_rp(1,2,6),riesz_rp(1,ldim,6))
c     elseif (angle<-85) then 
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,2),riesz_rp(1,2,2),riesz_rp(1,ldim,2))      
c     else
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,1),riesz_rp(1,2,1),riesz_rp(1,ldim,1))      
c     elseif (angle<-45) then
c     elseif (angle.ge.30.and.angle<65) then
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,3),riesz_rp(1,2,3),riesz_rp(1,ldim,3))
c     elseif (angle>-10.and.angle.lt.10) then
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,7),riesz_rp(1,2,7),riesz_rp(1,ldim,7))
c     elseif (angle>-65.and.angle.lt.-30) then
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,4),riesz_rp(1,2,4),riesz_rp(1,ldim,4))
c     elseif (angle>65) then
c        call opadd2(res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $        riesz_rp(1,1,5),riesz_rp(1,2,5),riesz_rp(1,ldim,5))
c     endif

      write(6,*)sin(angle*pi/180),'sin(angle*pi/180)'
      do i=1,nb+1
         call cfill(wk1,uu(i+(nb+1)),n)
         call cmult(wk1,-sin(angle*pi/180),n)
         call add2col2(res_u(1,1),riesz_ru(1,1,l1),wk1,n)
         call add2col2(res_u(1,2),riesz_ru(1,2,l1),wk1,n)
         if (ldim.eq.3) then 
            call add2col2(res_u(1,ldim),riesz_ru(1,ldim,l1),wk1,n)
         endif
         l1=l1+1 
      enddo
      write(6,*)l1,'l1'

      write(6,*)cos(angle*pi/180),'cos(angle*pi/180)'
      do i=1,nb+1
         call cfill(wk1,uu(i+(nb+1)),n)
         call cmult(wk1,-cos(angle*pi/180),n)
         call add2col2(res_u(1,1),riesz_ru(1,1,l1),wk1,n)
         call add2col2(res_u(1,2),riesz_ru(1,2,l1),wk1,n)
         if (ldim.eq.3) then 
            call add2col2(res_u(1,ldim),riesz_ru(1,ldim,l1),wk1,n)
         endif
c        call opadd2col (res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $                   riesz_ru(1,1,l1),riesz_ru(1,2,l1),
c    $                   riesz_ru(1,ldim,l1),wk1)
         l1=l1+1 
      enddo
      write(6,*)l1,'l2'

      do j=1,nb+1
         do i=1,nb+1
            call cfill(wk1,uu(j),n)
            call cmult(wk1,uu(i),n)
            call add2col2(res_u(1,1),riesz_ru(1,1,l1),wk1,n)
            call add2col2(res_u(1,2),riesz_ru(1,2,l1),wk1,n)
            if (ldim.eq.3) then 
               call add2col2(res_u(1,ldim),riesz_ru(1,ldim,l1),wk1,n)
            endif
c           call opadd2col (res_u(1,1),res_u(1,2),res_u(1,ldim),
c    $                   riesz_ru(1,1,l1),riesz_ru(1,2,l1),
c    $                   riesz_ru(1,ldim,l1),wk1)
            l1=l1+1 
         enddo
      enddo
      write(6,*)l1,'l3'

      l2=1
      do i=1,nb+1
         call cfill(wk1,uu(i+nb+1),n)
         call add2col2(res_t,riesz_rt(1,l2),wk1,n)
         l2=l2+1
      enddo
      write(6,*)l2,'l4'

      do j=1,nb+1
         do i=1,nb+1
            call cfill(wk1,uu(j),n)
            call cmult(wk1,uu(i+nb+1),n)
            call add2col2(res_t,riesz_rt(1,l2),wk1,n)
            l2=l2+1 
         enddo
      enddo
      write(6,*)l2,'l5'

      return
      end
c-----------------------------------------------------------------------
      subroutine cres_new1(angle)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real wk1(lt)
      real angle

      n=lx1*ly1*lz1*nelv
      write(6,*)nelt,nelv,lelm,'nelt,nelv'
      pi  = 4.*atan(1.)

      call rzero3(res_u(1,1),res_u(1,2),res_u(1,ldim),n)
      call rzero(res_t,n)

      l1=1
      do i=1,nb+1
         call opcopy(res_u(1,1),res_u(1,2),res_u(1,ldim),
     $               riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $               riesz_ru(1,ldim,l1))
         l1=l1+1 
      enddo

      write(6,*)sin(angle*pi/180),'sin(angle*pi/180)'
      do i=1,nb+1
         call cfill(wk1,-sin(angle*pi/180),n)
         call add2col2(res_u(1,1),riesz_ru(1,1,l1),wk1,n)
         call add2col2(res_u(1,2),riesz_ru(1,2,l1),wk1,n)
         if (ldim.eq.3) then 
            call add2col2(res_u(1,ldim),riesz_ru(1,ldim,l1),wk1,n)
         endif
         l1=l1+1 
      enddo
      write(6,*)l1,'l1'

      write(6,*)cos(angle*pi/180),'cos(angle*pi/180)'
      do i=1,nb+1
         call cfill(wk1,-cos(angle*pi/180),n)
         call add2col2(res_u(1,1),riesz_ru(1,1,l1),wk1,n)
         call add2col2(res_u(1,2),riesz_ru(1,2,l1),wk1,n)
         if (ldim.eq.3) then 
            call add2col2(res_u(1,ldim),riesz_ru(1,ldim,l1),wk1,n)
         endif
         l1=l1+1 
      enddo
      write(6,*)l1,'l2'

      do j=1,nb+1
         do i=1,nb+1
            call rone(wk1,n)
            call add2col2(res_u(1,1),riesz_ru(1,1,l1),wk1,n)
            call add2col2(res_u(1,2),riesz_ru(1,2,l1),wk1,n)
            if (ldim.eq.3) then 
               call add2col2(res_u(1,ldim),riesz_ru(1,ldim,l1),wk1,n)
            endif
            l1=l1+1 
         enddo
      enddo
      write(6,*)l1,'l3'

      l2=1
      do i=1,nb+1
         call rone(wk1,n)
         call add2col2(res_t,riesz_rt(1,l2),wk1,n)
         l2=l2+1
      enddo
      write(6,*)l2,'l4'

      do j=1,nb+1
         do i=1,nb+1
            call rone(wk1,n)
            call add2col2(res_t,riesz_rt(1,l2),wk1,n)
            l2=l2+1 
         enddo
      enddo
      write(6,*)l2,'l5'

      return
      end
c-----------------------------------------------------------------------
      subroutine crd_test(num_ts)

      include 'SIZE'
      include 'TOTAL'
      include 'ORTHOT'
      include 'CTIMER'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      integer num_ts
      real uu((nb+1)*2*num_ts),wk(2*(lb+1)*num_ts)
      real ts_an(num_ts)
      real dual_norm(num_ts),tmp(lt)
      real wk1(lt),wk2(lt),wk3(lt)
      real wk4(lt),wk5(lt),wk6(lt)
      real tt
      logical ifexist
      pi  = 4.*atan(1.)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside crd'

      inquire (file='./rom.dat',exist=ifexist)
      if (ifexist) call read_serial(uu,(nb+1)*2*num_ts,
     $                  './rom.dat ',wk,nid)

      inquire (file='./angle.dat',exist=ifexist)
      if (ifexist) call read_serial(ts_an,num_ts,
     $                  './angle.dat ',wk,nid)

      call rzero(dual_norm,num_ts)
      do i=1,num_ts

      write(6,*)i,ts_an(i),'angle'
      call cres_new(uu(1+(i-1)*(nb+1)*2),ts_an(i))
c     call cres_new1(angle)

      call rone(ones,n)
      call rzero(zeros,n)

      ifield=1
      tolhv=1e-8
      tolht(2)=1e-8
c     nmxv=1000
      call ophinv(wk1,wk2,wk3,
     $            riesz_ru(1,1,2),riesz_ru(1,2,2),riesz_ru(1,ldim,2),
     $            ones,ones,tolhv,nmxv)      
      ifield=1
      t1 = h10vip(wk1,wk2,wk3,wk1,wk2,wk3)
      t2 = wl2vip(wk1,wk2,wk3,wk1,wk2,wk3)
      write(6,*)t1,t2,t1+t2,'first term'
      call opcopy(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),wk1,wk2,wk3)

      call cmult(riesz_ru(1,1,4),-sin(ts_an(i)*pi/180),n)
      call cmult(riesz_ru(1,2,4),-sin(ts_an(i)*pi/180),n)
      call ophinv(wk1,wk2,wk3,
     $            riesz_ru(1,1,4),riesz_ru(1,2,4),riesz_ru(1,ldim,4),
     $            ones,ones,tolhv,nmxv)      
      ifield=1
      t1 = h10vip(wk1,wk2,wk3,wk1,wk2,wk3)
      t2 = wl2vip(wk1,wk2,wk3,wk1,wk2,wk3)
      write(6,*)t1,t2,t1+t2,'second term'
      call opadd2(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),wk1,wk2,wk3)

      call cmult(riesz_ru(1,1,6),-cos(ts_an(i)*pi/180),n)
      call cmult(riesz_ru(1,2,6),-cos(ts_an(i)*pi/180),n)
      call ophinv(wk1,wk2,wk3,
     $            riesz_ru(1,1,6),riesz_ru(1,2,6),riesz_ru(1,ldim,6),
     $            ones,ones,tolhv,nmxv)      
      ifield=1
      t1 = h10vip(wk1,wk2,wk3,wk1,wk2,wk3)
      t2 = wl2vip(wk1,wk2,wk3,wk1,wk2,wk3)
      write(6,*)t1,t2,t1+t2,'third term'
      call opadd2(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),wk1,wk2,wk3)

      call ophinv(wk1,wk2,wk3,
     $            riesz_ru(1,1,10),riesz_ru(1,2,10),riesz_ru(1,ldim,10),
     $            ones,ones,tolhv,nmxv)      
      ifield=1
      t1 = h10vip(wk1,wk2,wk3,wk1,wk2,wk3)
      t2 = wl2vip(wk1,wk2,wk3,wk1,wk2,wk3)
      write(6,*)t1,t2,t1+t2,'fourth term'
      call opadd2(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),wk1,wk2,wk3)

c     call opgrad(wk4,wk5,wk6,pb(1,0))
c     call ophinv(wk1,wk2,wk3,
c    $            wk4,wk5,wk6,
c    $            ones,ones,tolhv,nmxv)      
      ifield=1
c     t1 = h10vip(wk1,wk2,wk3,wk1,wk2,wk3)
c     t2 = wl2vip(wk1,wk2,wk3,wk1,wk2,wk3)
c     write(6,*)t1,t2,t1+t2,'pressure term'
c     call opsub2(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),wk1,wk2,wk3)

      ifield=2
      ifld1 = ifield-1
      napproxt(1,ifld1) = laxtt
c     call hsolve  ('TEMP',eh_t,res_t,ones,ones
c    $                   ,tmask(1,1,1,1,ifield-1)
c    $                   ,tmult(1,1,1,1,ifield-1)
c    $                   ,imesh,tolht(ifield),nmxt(ifield-1),1
c    $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
      call hsolve  ('TEMP',eh_t,res_t,ones,ones
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxh,1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)

      ifield=1
      t1 = h10vip(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
     $         eh_u(1,1),eh_u(1,2),eh_u(1,ldim)) 
      t2 = wl2vip(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
     $           eh_u(1,1),eh_u(1,2),eh_u(1,ldim)) 
      ifield=2
      t3 = h10sip(eh_t,eh_t) 
      t4 = wl2sip(eh_t,eh_t)  
      
      dual_norm(i) = t1+t2+t3+t4
      write(6,*)t1,t2,t3,t4,'t1,t2,t3,t4'
c     dual_norm(i) = h10vip(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
c    $                   eh_u(1,1),eh_u(1,2),eh_u(1,ldim)) +
c    $               h10sip(eh_t,eh_t) 
      dual_norm(i) = sqrt(dual_norm(i))
      write(6,*)i,dual_norm(i),'dual_norm in for v and t'
      write(6,*)i,sqrt(t1+t2),'dual_norm in for v'
      write(6,*)i,sqrt(t3+t4),'dual_norm in for t'

      enddo
      call dump_serial(dual_norm,num_ts,'./dual_norm ',nid)

      if (nio.eq.0) write (6,*) 'exiting crd'

      return
      end
c-----------------------------------------------------------------------
      subroutine set_sNS_divfrr

c     Compute the divergence free
c     riesz representators for steady Boussinesq 
c     incompressible NS (quadratically nonlinear elliptic problem)

      include 'SIZE'
      include 'TOTAL'
      include 'ORTHOT'
      include 'CTIMER'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      real work(lx2*ly2*lz2*lelt)

      common /screi/ wk1(lt),wk2(lt),wk3(lt),wk4(lt),wk5(lt)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside set_sNS_divfrr'

      call rone(ones,n)
      call rzero(zeros,n)

      l1=1
      l2=1
      do i=0,nb
         ifield=1
         tolhv=1e-8
         tolht(2)=1e-8
         call steady_stoke_solve(xi_u(1,1,l1),xi_u(1,2,l1),
     $       xi_u(1,ldim,l1),work,riesz_ru(1,1,l1),    
     $       riesz_ru(1,2,l1),riesz_ru(1,ldim,l1))
         if (nid.eq.0) write(6,*)'riesz_u',l1,'completed'

         ifield=2
         ifld1 = ifield-1
         napproxt(1,ifld1) = laxtt
         ! nmxt does not supported in nek5000 master branch
         call hsolve  ('TEMP',xi_t(1,l2),riesz_rt(1,l2),ones,ones
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxh,1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
         if (nid.eq.0) write(6,*)'riesz_t',l2,'completed'
         l1=l1+1
         l2=l2+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_1'
      if (nid.eq.0) write(6,*)l2,'lres_t_1'

      do i=0,nb
         ifield=1
         tolhv=1e-8
         tolht(2)=1e-8
         call steady_stoke_solve(xi_u(1,1,l1),xi_u(1,2,l1),
     $       xi_u(1,ldim,l1),work,riesz_ru(1,1,l1),    
     $       riesz_ru(1,2,l1),riesz_ru(1,ldim,l1))
         if (nid.eq.0) write(6,*)'riesz_u',l1,'completed'
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_2'

      do i=0,nb
         ifield=1
         tolhv=1e-8
         tolht(2)=1e-8
         call steady_stoke_solve(xi_u(1,1,l1),xi_u(1,2,l1),
     $       xi_u(1,ldim,l1),work,riesz_ru(1,1,l1),    
     $       riesz_ru(1,2,l1),riesz_ru(1,ldim,l1))
         if (nid.eq.0) write(6,*)'riesz_u',l1,'completed'
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_3'

      do j=0,nb
         do i=0,nb
            ifield=1
            tolhv=1e-8
            tolht(2)=1e-8
            call steady_stoke_solve(xi_u(1,1,l1),xi_u(1,2,l1),
     $          xi_u(1,ldim,l1),work,riesz_ru(1,1,l1), 
     $          riesz_ru(1,2,l1),riesz_ru(1,ldim,l1))
            if (nid.eq.0) write(6,*)'riesz_u',l1,'completed'

            ifield=2
            ifld1 = ifield-1
            napproxt(1,ifld1) = laxtt
            call hsolve  ('TEMP',xi_t(1,l2),riesz_rt(1,l2),ones,ones
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxh,1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
            l1=l1+1
            l2=l2+1
            if (nid.eq.0) write(6,*)'riesz_t',l2,'completed'
         enddo
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_4'
      if (nid.eq.0) write(6,*)l2,'lres_t_2'

      ! umcomment this part if the riesz is not divergence free
c     ifield=1
c     call ophinv(xi_u(1,1,l1),xi_u(1,2,l1),xi_u(1,ldim,l1),
c    $               riesz_rp(1,1,1),riesz_rp(1,2,1),
c    $               riesz_rp(1,ldim,1),
c    $               ones,ones,tolhv,nmaxv)               
c     l1=l1+1
c     if (nid.eq.0) write(6,*)'riesz_u',l1,'completed'

      if ((l1-1).gt.nres_u) then
         call exitti('increase nres_u$',l1-1)
      endif
      if ((l2-1).gt.nres_t) then
         call exitti('increase nres_t$',l2-1)
      endif

      if (nio.eq.0) write (6,*) 'exiting set_sNS_divfrr'

      return
      end
c-----------------------------------------------------------------------
      subroutine steady_stoke_solve(ehu,ehv,ehw,ehp,rhs1,rhs2,rhs3)

      include 'SIZE'
      include 'TOTAL'
      include 'ORTHOT'
      include 'CTIMER'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ltt=lx2*ly2*lz2*lelt)

      real ehu(lt),ehv(lt),ehw(lt),ehp(ltt)
      real dv1(lt),dv2(lt),dv3(lt)
      real h1(lt),h2(lt),h2inv(lt)
      real rhs1(lt),rhs2(lt),rhs3(lt)
      real tmp1(lt),tmp2(lt),tmp3(lt)
      real g1(lt),g2(lt),g3(lt)
      real wp(ltt)
      real div_check(lt)
      real tmp(lt)
      logical IFSTUZ

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside steady_stoke solver'

      IFTRAN = .FALSE.
      IFADVC(IFIELD) = .FALSE.

      ! set vdiff to be one so that the constant for stokes operator 
      ! is correct
      call copy(tmp,vdiff(1,1,1,1,ifield),n)
      call rone(vdiff(1,1,1,1,ifield),n)

      call rzero(vx,n)
      call rzero(vy,n)
      call rzero(vz,n)
c     call bcdirvc (vx,vy,vz,v1mask,v2mask,v3mask)

C     Check if steady state
C
      IFSTUZ = .FALSE.
      CALL CONVUZ (IFSTUZ)
C... no steady state
      IFSTUZ = .FALSE.
      IF (IFSTUZ) THEN
         IF (NIO.EQ.0) WRITE (6,*) 
     $      'Steady state reached in the fluid solver'
            return
         ENDIF
C
C     Uzawa decoupling: First, compute pressure.....
C
      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv
      call rzero(pr,ntot2)

      intype = 0
      if (iftran) intype = -1
      call sethlm  (h1,h2,intype)
      call rone(h1,ntot1)
      if (iftran)  call invers2 (h2inv,h2,ntot1)
      call opcopy(bfx,bfy,bfz,rhs1,rhs2,rhs3)
      call makeg   (   g1,g2,g3,h1,h2,intype)
      call crespuz (wp,g1,g2,g3,h1,h2,h2inv,intype)
      call uzawa   (wp,h1,h2,h2inv,intype,icg)
      if (icg.gt.0) call add2 (pr,wp,ntot2)

C     .... then, compute velocity:
      call cresvuz (tmp1,tmp2,tmp3)
      call ophinv  (dv1,dv2,dv3,tmp1,tmp2,tmp3,h1,h2,tolhv,nmxv)

      call opcopy(ehu,ehv,ehw,dv1,dv2,dv3)
      call copy(ehp,wp,ntot2)

c     call opdiv(div_check,ehu,ehv,ehw)
c     write(6,*)'check divergence of riesz'
c     write(6,*)glmax(div_check,ntot1),glmin(div_check,ntot1)

      ! set vdiff back to original value
      call copy(vdiff(1,1,1,1,ifield),tmp,ntot1)

      if (nio.eq.0) write (6,*) 'exiting steady_stoke solver'

      return
      end
c-----------------------------------------------------------------------
      subroutine crd_divf(num_ts)

      ! make representation to be divergence free

      include 'SIZE'
      include 'TOTAL'
      include 'ORTHOT'
      include 'CTIMER'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real uu((nb+1)*2*num_ts),wk(2*(lb+1)*num_ts)
      real ts_an(num_ts)
      real dual_norm(num_ts),tmp(lt)
      integer num_ts
      logical ifexist

      n=lx1*ly1*lz1*nelv
      nn=lx2*ly2*lz2*nelv

      if (nio.eq.0) write (6,*) 'inside crd'

      inquire (file='./rom.dat',exist=ifexist)
      if (ifexist) call read_serial(uu,(nb+1)*2*num_ts,
     $                  './rom.dat ',wk,nid)

      inquire (file='./angle.dat',exist=ifexist)
      if (ifexist) call read_serial(ts_an,num_ts,
     $                  './angle.dat ',wk,nid)

      write(6,*)num_ts,'num_ts'
      call rzero(dual_norm,num_ts)
      do i=1,num_ts

         write(6,*)i,ts_an(i),'angle'
         call cres_new(uu(1+(i-1)*(nb+1)*2),ts_an(i))
c        do ii=1,(nb+1)*2*num_ts
c           write(6,*)ii,uu(ii),ts_an(i)
c        enddo

         call rone(ones,n)
         call rzero(zeros,n)
         ifield=1
         tolhv=1e-8
         tolht(2)=1e-8
         call steady_stoke_solve(eh_u(1,1),eh_u(1,2),
     $       eh_u(1,ldim),work,res_u(1,1),
     $       res_u(1,2),res_u(1,ldim))

c        t5 = glsc3(xi_p(1,1),xi_p(1,1),bm2,nn)

         ifield=2
         ifld1 = ifield-1
         napproxt(1,ifld1) = laxtt
         call hsolve  ('TEMP',eh_t,res_t,ones,ones
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxh,1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)

         ifield=1
         t1 = h10vip(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
     $         eh_u(1,1),eh_u(1,2),eh_u(1,ldim)) 
         t2 = wl2vip(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
     $           eh_u(1,1),eh_u(1,2),eh_u(1,ldim)) 
         ifield=2
         t3 = h10sip(eh_t,eh_t) 
         t4 = wl2sip(eh_t,eh_t)  
      
         write(6,*)t1,t2,t3,t4,t5,'t1,t2,t3,t4,t5'
         dual_norm(i) = t1+t2+t3+t4+t5
         dual_norm(i) = sqrt(dual_norm(i))
         write(6,*)i,ts_an(i),'angle'
         write(6,*)i,dual_norm(i),
     $            'dual_norm H1 norm for the coupled system'
         write(6,*)i,sqrt(t1+t2),
     $            'dual_norm H1 norm for the velocity'
         write(6,*)i,sqrt(t3+t4),
     $            'dual_norm H1 norm for the temperature'

      enddo
      call dump_serial(dual_norm,num_ts,'./dual_norm ',nid)

      if (nio.eq.0) write (6,*) 'exiting crd'

      return
      end
c-----------------------------------------------------------------------
      subroutine set_residual_unsteady

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /screi/ wk1(lt),wk2(lt),wk3(lt),wk4(lt),wk5(lt),wk6(lt)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside set_residual_unsteady'

      call nekgsync
      sru_start=dnekclock()

      call rone(ones,n)
      call rzero(zeros,n)
      
      ! setup rhs for velocity representator
      l1=1
      do i=0,nb
         ifield=1
         call opcopy(riesz_ru(1,1,l1),riesz_ru(1,2,l1)
     $               ,riesz_ru(1,ldim,l1)
     $               ,ub(1,i),vb(1,i),wb(1,i))
         call opcolv(riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $               riesz_ru(1,ldim,l1),bm1)
         call opchsgn(riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $                riesz_ru(1,ldim,l1))
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_1'

      do i=0,nb
         ifield=1
         call opcopy(wk1,wk2,wk3,ub(1,i),vb(1,i),wb(1,i))
         call axhelm(riesz_ru(1,1,l1),wk1,ones,zeros,1,1)
         call axhelm(riesz_ru(1,2,l1),wk2,ones,zeros,1,2)
         if (ldim.eq.3) then
            call axhelm(riesz_ru(1,ldim,l1),wk3,ones,zeros,1,3)
         endif
         call opcmult(riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $                riesz_ru(1,ldim,l1),param(2))
         call opchsgn(riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $                riesz_ru(1,ldim,l1))
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_2'

      do i=0,nb
         ifield=1
         call opcopy(riesz_ru(1,1,l1),riesz_ru(1,2,l1)
     $               ,riesz_ru(1,ldim,l1)
     $               ,tb(1,i),zeros,zeros)
         call opcolv(riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $               riesz_ru(1,ldim,l1),bm1)
         call opchsgn(riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $                riesz_ru(1,ldim,l1))
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_3'

      do i=0,nb
         ifield=1
         call opcopy(riesz_ru(1,1,l1),riesz_ru(1,2,l1)
     $               ,riesz_ru(1,ldim,l1)
     $               ,zeros,tb(1,i),zeros)
         call opcolv(riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $               riesz_ru(1,ldim,l1),bm1)
         call opchsgn(riesz_ru(1,1,l1),riesz_ru(1,2,l1),
     $                riesz_ru(1,ldim,l1))
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_4'

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
            l1=l1+1
         enddo
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_5'

      ! setup rhs for temperature representator
      l2=1
      do i=0,nb
         ifield=2
         call copy(riesz_rt(1,l2),tb(1,i),n)
         call col2(riesz_rt(1,l2),bm1,n)
         call chsign(riesz_rt(1,l2),n)
         l2=l2+1
      enddo
      if (nid.eq.0) write(6,*)l2,'lres_t_1'
      do i=0,nb
         ifield=2
         call copy(wk4,tb(1,i),n)
         call axhelm(riesz_rt(1,l2),wk4,ones,zeros,1,1)
         call cmult(riesz_rt(1,l2),param(8),n)
         call chsign(riesz_rt(1,l2),n) 
         l2=l2+1
      enddo
      if (nid.eq.0) write(6,*)l2,'lres_t_2'
      do j=0,nb
         do i=0,nb
         ifield=2
            call convect_new(riesz_rt(1,l2),tb(1,i),.false.,
     $                       ub(1,j),vb(1,j),wb(1,j),.false.)
            call chsign(riesz_rt(1,l2),n) 
            l2=l2+1
         enddo
      enddo
      if (nid.eq.0) write(6,*)l2,'lres_t_3'

      ! add pressure residual
      ! Uncomment this part when the riesz is not divergence free
      ! if the riesz is divergence free, we do not need pressure
      ! residual
c     l3=1
c     call ifield=1
c     do i=1,nb
c     call opgradt(riesz_rp(1,1,l3),riesz_rp(1,2,l3),
c    $            riesz_rp(1,ldim,l3),pb(1,i))
c     l3=l3+1
c     enddo

      call nekgsync
      if (nio.eq.0) write (6,*) 'set_residual_unsteady_time:',
     $             dnekclock()-sru_start

      if (nio.eq.0) write (6,*) 'exiting set_residual_unsteady'

      return
      end
c-----------------------------------------------------------------------
      subroutine set_uNS_divfrr

c     Compute the divergence free
c     riesz representators for unsteady Boussinesq 
c     incompressible NS (quadratically nonlinear elliptic problem)

      include 'SIZE'
      include 'TOTAL'
      include 'ORTHOT'
      include 'CTIMER'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      real work(lx2*ly2*lz2*lelt)

      common /screi/ wk1(lt),wk2(lt),wk3(lt),wk4(lt),wk5(lt)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside set_uNS_divfrr'

      call nekgsync
      srru_start=dnekclock()

      call rone(ones,n)
      call rzero(zeros,n)

      l1=1
      l2=1
      do i=0,nb
         ifield=1
         tolhv=1e-8
         tolht(2)=1e-8
c        nmxh=1000
         call steady_stoke_solve(xi_u(1,1,l1),xi_u(1,2,l1),
     $       xi_u(1,ldim,l1),work,riesz_ru(1,1,l1),    
     $       riesz_ru(1,2,l1),riesz_ru(1,ldim,l1))
         if (nid.eq.0) write(6,*)'riesz_u',l1,'completed'

         ifield=2
         ifld1 = ifield-1
         napproxt(1,ifld1) = laxtt
         call hsolve  ('TEMP',xi_t(1,l2),riesz_rt(1,l2),ones,ones
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxh,1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
         if (nid.eq.0) write(6,*)'riesz_t',l2,'completed'
         l1=l1+1
         l2=l2+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_1'
      if (nid.eq.0) write(6,*)l2,'lres_t_1'

      do i=0,nb
         ifield=1
         tolhv=1e-8
         tolht(2)=1e-8
c        nmxh=1000
         call steady_stoke_solve(xi_u(1,1,l1),xi_u(1,2,l1),
     $       xi_u(1,ldim,l1),work,riesz_ru(1,1,l1),    
     $       riesz_ru(1,2,l1),riesz_ru(1,ldim,l1))
         if (nid.eq.0) write(6,*)'riesz_u',l1,'completed'

         ifield=2
         ifld1 = ifield-1
         napproxt(1,ifld1) = laxtt
         call hsolve  ('TEMP',xi_t(1,l2),riesz_rt(1,l2),ones,ones
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxh,1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
         if (nid.eq.0) write(6,*)'riesz_t',l2,'completed'
         l1=l1+1
         l2=l2+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_2'
      if (nid.eq.0) write(6,*)l2,'lres_t_2'

      do i=0,nb
         ifield=1
         tolhv=1e-8
         tolht(2)=1e-8
         call steady_stoke_solve(xi_u(1,1,l1),xi_u(1,2,l1),
     $       xi_u(1,ldim,l1),work,riesz_ru(1,1,l1),    
     $       riesz_ru(1,2,l1),riesz_ru(1,ldim,l1))
         if (nid.eq.0) write(6,*)'riesz_u',l1,'completed'
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_3'

      do i=0,nb
         ifield=1
         tolhv=1e-8
         tolht(2)=1e-8
         call steady_stoke_solve(xi_u(1,1,l1),xi_u(1,2,l1),
     $       xi_u(1,ldim,l1),work,riesz_ru(1,1,l1),    
     $       riesz_ru(1,2,l1),riesz_ru(1,ldim,l1))
         if (nid.eq.0) write(6,*)'riesz_u',l1,'completed'
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_4'

      do j=0,nb
         do i=0,nb
            ifield=1
            tolhv=1e-8
            tolht(2)=1e-8
            call steady_stoke_solve(xi_u(1,1,l1),xi_u(1,2,l1),
     $          xi_u(1,ldim,l1),work,riesz_ru(1,1,l1), 
     $          riesz_ru(1,2,l1),riesz_ru(1,ldim,l1))
            if (nid.eq.0) write(6,*)'riesz_u',l1,'completed'

            ifield=2
            ifld1 = ifield-1
            napproxt(1,ifld1) = laxtt
            call hsolve  ('TEMP',xi_t(1,l2),riesz_rt(1,l2),ones,ones
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxh,1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
            l1=l1+1
            l2=l2+1
            if (nid.eq.0) write(6,*)'riesz_t',l2,'completed'
         enddo
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_5'
      if (nid.eq.0) write(6,*)l2,'lres_t_3'

      if (nid.eq.0) write(6,*)'riesz_u',l1,'completed'

      if ((l1-1).gt.nres_u) then
         call exitti('increase nres_u$',l1-1)
      endif
      if ((l2-1).gt.nres_t) then
         call exitti('increase nres_t$',l2-1)
      endif

      call nekgsync
      if (nio.eq.0) write (6,*) 'set_uNS_divfrr_time:',
     $             dnekclock()-srru_start

      return
      end
c-----------------------------------------------------------------------
      subroutine set_theta_uns

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside set_theta_uns'

      call rzero(theta_u,lres_u)
      call rzero(theta_t,lres_t)

      call set_betaj
      call set_alphaj

      ! time derivative (done)
      l1=1
      call mxm(uj,nb+1,betaj,6,theta_u(l1),1)

      l1=l1+nb+1

      ! diffusion term (done)
      do i=0,nb
         theta_u(l1)=ua(i)
         l1=l1+1
      enddo

      ! buoyancy term (done)
      if (ifbuoy) then
         call mxm(utj,(nb+1),alphaj,6,theta_u(l1),1)
         do i=0,nb
            theta_u(l1)=-sin(bu_angle)*ad_ra*(theta_u(l1)+uta_wol(i))
            l1=l1+1
         enddo
         call mxm(utj,(nb+1),alphaj,6,theta_u(l1),1)
         do i=0,nb
            theta_u(l1)=-cos(bu_angle)*ad_ra*(theta_u(l1)+uta_wol(i))
            l1=l1+1
         enddo
      endif

      ! convection (done)
      ! Order of u2j does not matter since it is symmetric
      call mxm(u2j,(nb+1)**2,alphaj,6,theta_u(l1),1)
      do j=0,nb
      do i=0,nb
         theta_u(l1)=theta_u(l1)+u2a_wol(1+i+(nb+1)*j)
         l1=l1+1
      enddo
      enddo

      do i=1,nres_u
         if (nio.eq.0) write (6,*) theta_u(i),'theta_u'
      enddo

      l2=1
      call mxm(utj,nb+1,betaj,6,theta_t(l2),1)

      l2=l2+nb+1

      do i=0,nb
         theta_t(l2)=uta(i)
         l2=l2+1
      enddo

      ! i for temperture, j for velocity, j should change to k
      ! should use utuj, in utuj, velocity uses j index
      call mxm(utuj,(nb+1)**2,alphaj,6,theta_t(l2),1)
      do j=0,nb
      do i=0,nb
         theta_t(l2)=theta_t(l2)+utua_wol(1+i+(nb+1)*j)
         l2=l2+1
      enddo
      enddo

      do i=1,nres_t
         if (nio.eq.0) write (6,*) theta_t(i),'theta_t'
      enddo

      if (nio.eq.0) write (6,*) 'exitng set_theta_uns'

      return
      end
c-----------------------------------------------------------------------
      subroutine cres_uns(aa,bb)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real aa(1:nres_u,1:nres_u)
      real bb(1:nres_t,1:nres_t)

      call set_theta_uns

      res_uu=0.
      res_tt=0.
      res=0.

      do j=1,nres_u
      do i=1,nres_u
         res_uu=res_uu+aa(i,j)*theta_u(i)*theta_u(j)
      enddo
      enddo
      if (nid.eq.0)write(6,*)'velocity dual norm',sqrt(res_uu)

      do j=1,nres_t
      do i=1,nres_t
         res_tt=res_tt+bb(i,j)*theta_t(i)*theta_t(j)
      enddo
      enddo
      if (nid.eq.0)write(6,*)'temp dual norm',sqrt(res_tt)

      res=sqrt(res_uu+res_tt)

      if (res.le.0) call exitti('negative semidefinite dual norm$',n)

      if (nid.eq.0) write (6,*) 'dual norm:',res

      return
      end
c-----------------------------------------------------------------------
      subroutine c_rieszrd_uns

      include 'SIZE'
      include 'TOTAL'
      include 'ORTHOT'
      include 'CTIMER'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /screi/ wk1(lt),wk2(lt),wk3(lt),wk4(lt),wk5(lt),wk6(lt)

      real work(lx2*ly2*lz2*lelt)
      real coef(0:nb)
      real coef2(0:(nb+1)**2-1)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside c_rieszrd_uns'

      call nekgsync

c     call set_theta_uns
      call set_betaj
      call set_alphaj

      call opzero(res_u(1,1),res_u(1,2),res_u(1,ldim))
      call rzero(res_t,n)

      ! use res_u and res_t for the velocity and temperature residual
      ! setup res_u for velocity representator
      l1=1
      call mxm(uj,nb+1,betaj,6,coef,1)
      do i=0,nb
         ifield=1
         call opcopy(wk1,wk2,wk3,ub(1,i),vb(1,i),wb(1,i))
         call opcolv(wk1,wk2,wk3,bm1)
         call opchsgn(wk1,wk2,wk3)

         call cfill(wk4,coef(i),n)
         call add2col2(res_u(1,1),wk1,wk4,n)
         call add2col2(res_u(1,2),wk2,wk4,n)
         if (ldim.eq.3) then
            call add2col2(res_u(1,ldim),wk3,wk4,n)
         endif
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_1'

      do i=0,nb
         ifield=1
         call opcopy(wk4,wk5,wk6,ub(1,i),vb(1,i),wb(1,i))
         call axhelm(wk1,wk4,ones,zeros,1,1)
         call axhelm(wk2,wk5,ones,zeros,1,2)
         if (ldim.eq.3) then
            call axhelm(wk3,wk6,ones,zeros,1,3)
         endif
         call opcmult(wk1,wk2,wk3,param(2))
         call opchsgn(wk1,wk2,wk3)

         coef(i)=ua(i)
         call cfill(wk4,coef(i),n)
         call add2col2(res_u(1,1),wk1,wk4,n)
         call add2col2(res_u(1,2),wk2,wk4,n)
         if (ldim.eq.3) then
            call add2col2(res_u(1,ldim),wk3,wk4,n)
         endif
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_2'

      call mxm(utj,(nb+1),alphaj,6,coef,1)
      do i=0,nb
         ifield=1
         call opcopy(wk1,wk2,wk3,tb(1,i),zeros,zeros)
         call opcolv(wk1,wk2,wk3,bm1)
         call opchsgn(wk1,wk2,wk3)

         coef(i)=-sin(bu_angle)*ad_ra*(coef(i)+uta_wol(i))
         call cfill(wk4,coef(i),n)
         call add2col2(res_u(1,1),wk1,wk4,n)
         call add2col2(res_u(1,2),wk2,wk4,n)
         if (ldim.eq.3) then
            call add2col2(res_u(1,ldim),wk3,wk4,n)
         endif
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_3'

      call mxm(utj,(nb+1),alphaj,6,coef,1)
      do i=0,nb
         ifield=1
         call opcopy(wk1,wk2,wk3,zeros,tb(1,i),zeros)
         call opcolv(wk1,wk2,wk3,bm1)
         call opchsgn(wk1,wk2,wk3)

         coef(i)=-cos(bu_angle)*ad_ra*(coef(i)+uta_wol(i))
         call cfill(wk4,coef(i),n)
         call add2col2(res_u(1,1),wk1,wk4,n)
         call add2col2(res_u(1,2),wk2,wk4,n)
         if (ldim.eq.3) then
            call add2col2(res_u(1,ldim),wk3,wk4,n)
         endif
         l1=l1+1
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_4'

      call mxm(u2j,(nb+1)**2,alphaj,6,coef2,1)
      do j=0,nb
         do i=0,nb
         ifield=1
            call convect_new(wk1,ub(1,i),.false.,
     $                       ub(1,j),vb(1,j),wb(1,j),.false.)
            call convect_new(wk2,vb(1,i),.false.,
     $                       ub(1,j),vb(1,j),wb(1,j),.false.)
            if (ldim.eq.3) then
               call convect_new(wk3,wb(1,i),.false.,
     $                          ub(1,j),vb(1,j),wb(1,j),.false.)
            endif
            call opchsgn(wk1,wk2,wk3)

            coef2(1+i+(nb+1)*j)=coef2(1+i+(nb+1)*j)
     $                         +u2a_wol(1+i+(nb+1)*j)
            call cfill(wk4,coef2(1+i+(nb+1)*j),n)
            call add2col2(res_u(1,1),wk1,wk4,n)
            call add2col2(res_u(1,2),wk2,wk4,n)
            if (ldim.eq.3) then
               call add2col2(res_u(1,ldim),wk3,wk4,n)
            endif
            l1=l1+1
         enddo
      enddo
      if (nid.eq.0) write(6,*)l1,'lres_u_5'

      ! setup res_t for temperature representator
      l2=1
      call mxm(utj,nb+1,betaj,6,coef,1)
      do i=0,nb
         ifield=2
         call copy(wk1,tb(1,i),n)
         call col2(wk1,bm1,n)
         call chsign(wk1,n)

         call cfill(wk2,coef(i),n)
         call add2col2(res_t,wk1,wk2,n)
         l2=l2+1
      enddo
      if (nid.eq.0) write(6,*)l2,'lres_t_1'
      do i=0,nb
         ifield=2
         call copy(wk2,tb(1,i),n)
         call axhelm(wk1,wk2,ones,zeros,1,1)
         call cmult(wk1,param(8),n)
         call chsign(wk1,n) 

         coef(i)=uta(i)
         call cfill(wk2,coef(i),n)
         call add2col2(res_t,wk1,wk2,n)
         l2=l2+1
      enddo
      if (nid.eq.0) write(6,*)l2,'lres_t_2'

      call mxm(utuj,(nb+1)**2,alphaj,6,coef2,1)
      do j=0,nb
         do i=0,nb
            ifield=2
            call convect_new(wk1,tb(1,i),.false.,
     $                       ub(1,j),vb(1,j),wb(1,j),.false.)
            call chsign(wk1,n) 

            coef2(1+i+(nb+1)*j)=coef2(1+i+(nb+1)*j)
     $                         +utua_wol(1+i+(nb+1)*j)
            call cfill(wk2,coef2(1+i+(nb+1)*j),n)
            call add2col2(res_t,wk1,wk2,n)
            l2=l2+1
         enddo
      enddo
      if (nid.eq.0) write(6,*)l2,'lres_t_3'


      ! solve two stokes problem, one for velocity and one for temp
      call rone(ones,n)
      call rzero(zeros,n)

      ifield=1
      tolhv=1e-8
      tolht(2)=1e-8
      call steady_stoke_solve(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
     $     work,res_u(1,1),res_u(1,2),res_u(1,ldim))
      if (nid.eq.0) write(6,*)'riesz_u completed'

      ifield=2
      ifld1 = ifield-1
      napproxt(1,ifld1) = laxtt
      call hsolve  ('TEMP',eh_t,res_t,ones,ones
     $             ,tmask(1,1,1,1,ifield-1)
     $             ,tmult(1,1,1,1,ifield-1)
     $             ,imesh,tolht(ifield),nmxh,1
     $             ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
      if (nid.eq.0) write(6,*)'riesz_t completed'

      ifield=1
      t1 = h10vip(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
     $         eh_u(1,1),eh_u(1,2),eh_u(1,ldim)) 
      t2 = wl2vip(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
     $           eh_u(1,1),eh_u(1,2),eh_u(1,ldim)) 
      ifield=2
      t3 = h10sip(eh_t,eh_t) 
      t4 = wl2sip(eh_t,eh_t)  

      res_uu=0.
      res_uu=(t1+t2)

      res_tt=0.
      res_tt=(t3+t4)

      res=0.
      res=sqrt(t1+t2+t3+t4)

      if (nid.eq.0)write(6,*)'velocity dual norm',sqrt(res_uu)
      if (nid.eq.0)write(6,*)'temp dual norm',sqrt(res_tt)

      if (res.le.0) call exitti('negative semidefinite dual norm$',n)

      if (nid.eq.0) write (6,*) 'dual norm:',res

      call nekgsync
      if (nio.eq.0) write (6,*) 'exiting c_rieszrd_uns'

      return
      end
