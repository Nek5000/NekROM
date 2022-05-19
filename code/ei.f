c-----------------------------------------------------------------------
      subroutine cres

      ! This subroutine has been merged to cdres except
      ! the computing eierr part

      ! TODO: Verify and remove

      include 'SIZE'
      include 'MOR'

      common /ccres/ cdiff(0:lb)

      parameter (lt=lx1*ly1*lz1*lelt)

      if (eqn.eq.'POIS') then
         call set_theta_poisson
      else if (eqn.eq.'HEAT') then
         call set_theta_heat
      else if (eqn.eq.'ADVE') then
         call set_theta_ad
      else if (eqn.eq.'VNS ') then
         call set_theta_ns
      endif

      res=0.

      do j=1,nres
      do i=1,nres
         res=res+mor_sigma(i,j)*mor_theta(i)*mor_theta(j)
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
         res=res+mor_sigma(i,j)*mor_theta(i)*mor_theta(j)
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
      subroutine set_sigma

      ! Compute the inner product between
      ! riesz respresentations and store it in
      ! sigma matrix

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      if (nio.eq.0) write (6,*) 'inside set_sigma'

      call nekgsync
      sigma_time=dnekclock()

      if (lei.eq.0) then
         ! do nothing if only one residual
         continue
      else if (lei.eq.1) then
         ! compute the size of the sigma matrix based on the eqn
         call set_nres

         ! Read or construct sigma in the following
         if (rmode.eq.'ON '.or.rmode.eq.'ONB') then
            call read_sigma
         else
            call set_riesz_reps
            call c_sigma
         endif
      else
         call exitti('Check your lei$',lei)
      endif

      call nekgsync

      if (nio.eq.0) write (6,*) 'sigma_time:',dnekclock()-sigma_time
      if (nio.eq.0) write (6,*) 'exiting set_sigma'

      return
      end
c-----------------------------------------------------------------------
      subroutine set_nres

      ! Compute the size of the sigma matrix
      ! based on the eqn

      ! TODO: Need to integrate the nres, nres_u and nres_t
      ! for different equations

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      if (eqn.eq.'POIS') then
         nres=nb+1
         nres=nb+1+nb
      else if (eqn.eq.'HEAT') then
         nres=(nb+1)*2+1
      else if (eqn.eq.'ADVE') then
         nres=(nb+1)*3
         if (ifrom(1)) nres=(nb+1)*2 + (nb+1)**2
      else if (eqn.eq.'VNS ') then
         nres=(nb+1)*2+(nb+1)**2
         if (ifbuoy) nres=nres+nb+1
      else if (eqn.eq.'NS  ') then
         nres_u=(nb+1)*3+(nb+1)**2
         if (ifrom(2)) nres_t=(nb+1)*2+(nb+1)**2
         if (ifbuoy) nres_u=nres_u+nb+1
      else if (eqn.eq.'SNS ') then
         nres_u=(nb+1)*2+(nb+1)**2
         if (ifrom(2)) nres_t=(nb+1)*1+(nb+1)**2
         if (ifbuoy) nres_u=nres_u+nb+1
      endif

      ! check whether allocation for sigma matrix is enough
      if (eqn.eq.'NS  '.OR.eqn.eq.'SNS ') then
         if (nres_u.gt.lres_u)
     $       call exitti('nres_u > lres_u$',nres_u)
         if (nres_u.le.0) call exitti('nres_u <= 0$',nres_u)
         if (ifrom(2)) then
            if (nres_t.gt.lres_t)
     $         call exitti('nres_t > lres_t$',nres_t)
            if (nres_t.le.0) call exitti('nres_t <= 0$',nres_t)
         endif
      else
         if (nres.gt.lres) call exitti('nres > lres$',nres)
         if (nres.le.0) call exitti('nres <= 0$',nres)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine read_sigma

      ! Read sigma matrix
      
      ! TODO: Need to resolve the sigma_u, sigma_t, and sigma
      ! inconsistency

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      real wk1(lt)

      n=lx1*ly1*lz1*nelv

      if (eqn.eq.'NS  '.OR.eqn.eq.'SNS ') then
         if (nio.eq.0) write (6,*) 'reading sigma_u...'
         call read_sigma_u_serial(sigma_u,nres_u,nres_u,
     $                             'ops/sigma_u ',lres_u,nres_u,
     $                             (lb+1),(nb+1),wk1,nid)
         if (nio.eq.0) write (6,*) 'reading sigma_t...'
         call read_sigma_t_serial(sigma_t,nres_t,nres_t,
     $                             'ops/sigma_t ',lres_t,nres_t,
     $                             (lb+1),(nb+1),wk1,nid)

      else
         if (nio.eq.0) write (6,*) 'reading sigma...'
         call read_serial(sigtmp,nres*nres,'ops/sigma '
     $                     ,mor_sigma,nid)
         l=1
         do j=1,nres
         do i=1,nres
            mor_sigma(i,j)=sigtmp(l,1)
            l=l+1
         enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine c_sigma

      ! Compute sigma matrix

      ! TODO: Need to resolve the xi_u, xi_t, and xi inconsistency
      ! TODO: Need to resolve the sigma_u, sigma_t, and sigma inconsistency

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      if (eqn.eq.'NS  '.OR.eqn.eq.'SNS ') then
         call csigma_u(sigma_u)
         call csigma_t(sigma_t)
         call dump_serial(sigma_u,(nres_u)**2,'ops/sigma_u ',nid)
         call dump_serial(sigma_t,(nres_t)**2,'ops/sigma_t ',nid)
      else
         if (ifield.eq.2) then
            do i=1,nres
            do j=1,nres
               mor_sigma(i,j)=glsc3(xi(1,i),xi(1,j),bm1,n)
            enddo
            enddo
         else if (ifield.eq.1) then
            if (.true.) then ! if using voritcity residual
               do i=1,nres
                  do j=1,nres
                     mor_sigma(i,j)=glsc3(xi_u(1,1,i),xi_u(1,1,j),bm1,n)
                  enddo
                  write (6,*) i,mor_sigma(1,i),'sigma'
               enddo
            else
               do i=1,nres
               do j=1,nres
                  mor_sigma(i,j)=
     $               op_glsc2_wt(xi_u(1,1,i),xi_u(1,2,i),xi_u(1,ldim,i),
     $                       xi_u(1,1,j),xi_u(1,2,j),xi_u(1,ldim,j),bm1)
               enddo
               enddo
            endif
         endif
         call dump_serial(sigma,(nres)**2,'ops/sigma ',nid)
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
         mor_theta(l)=ut(i)
         l=l+1
      enddo

      mor_theta(l)=-1.
      l=l+1

      do i=1,nb
         mor_theta(l)=-ut(i)
         mor_theta(l)=0.
         l=l+1
      enddo

      do i=1,l-1
         write (6,*) i,mor_theta(i),'theta'
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
      call mxm(utj,nb+1,betaj,8,mor_theta(l),1)

      l=l+nb+1
      do i=0,nb
         mor_theta(l)=uta(i)
         l=l+1
      enddo

      mor_theta(l)=-1.

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
      call mxm(utj,nb+1,betaj,8,mor_theta(l),1)

      l=l+nb+1

      call set_alphaj

      if (ifrom(1)) then
c        call mxm(uutj,(nb+1)**2,alphaj,8,mor_theta(l),1)
         call mxm(utuj,(nb+1)**2,alphaj,8,mor_theta(l),1)
      else
         call mxm(utj,nb+1,alphaj,8,mor_theta(l),1)
      endif

      if (ifrom(1)) then
         do j=0,nb
         do i=0,nb
            mor_theta(l)=mor_theta(l)+utua(i+(nb+1)*j)
            l=l+1
         enddo
         enddo
      else
         do i=0,nb
            mor_theta(l)=mor_theta(l)+uta(i)
            l=l+1
         enddo
      endif

      do i=0,nb
         mor_theta(l)=param(8)*uta(i)
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

      if (ifbuoy) call exitti('ifbuoy temporarily disabled in EI',1) 

      l=1
      call set_betaj
      call mxm(utj,nb+1,betaj,8,mor_theta(l),1)

      l=l+nb+1

      call set_alphaj
      call mxm(ut2j,(nb+1)**2,alphaj,8,mor_theta(l),1)
      do j=0,nb
      do i=0,nb
         mor_theta(l)=mor_theta(l)+u2a(1+i+(nb+1)*j)
         l=l+1
      enddo
      enddo

      do i=0,nb
         mor_theta(l)=param(2)*ua(i)
         l=l+1
      enddo

      if (ifbuoy) then
         do i=0,nb
            mor_theta(l)=ad_ra*uta(i)
            l=l+1
         enddo
      endif

      do i=1,nres
         if (nio.eq.0) write (6,*) mor_theta(i),'theta'
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

      eqn='POIS'

      ifflow=.false.
      ifheat=.true.

      param(174)=-1.

      call rom_init_params
      call rom_init_fields

c     call setgram
c     call setevec

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
         rhs(i)=glsc2(qq,tb(1,i,1),n)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_alphaj

      ! set the alphaj coefficients for the extrapolation terms

      ! i := map from the ei time-step to the time-stepper order

      include 'SIZE'
      include 'MOR'

      integer i(3)
      equivalence (its,i)

      ! ad_alpha(3,3)

      alphaj(1)=ad_alpha(3,i(1))
      alphaj(2)=ad_alpha(2,i(1))+ad_alpha(3,i(2))
      alphaj(3)=ad_alpha(1,i(1))+ad_alpha(2,i(2))+ad_alpha(3,i(3))
      alphaj(4)=ad_alpha(1,i(2))-ad_alpha(1,i(3))
      alphaj(5)=0.
      alphaj(6)=-ad_alpha(3,i(3))
      alphaj(7)=-ad_alpha(2,i(3))-ad_alpha(3,i(3))
      alphaj(8)=-1.

      call cmult(alphaj,1./(1.*(ad_nsteps-navg_step+1)),8)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_betaj

      ! set the betaj coefficients for the time-derivative terms

      ! i := map from the ei time-step to the time-stepper order

      include 'SIZE'
      include 'MOR'

      integer i(3)
      equivalence (its,i)

      ! ad_beta(4,3)

      betaj(1)=ad_beta(3+1,i(1))
      betaj(2)=ad_beta(2+1,i(1))+ad_beta(3+1,i(2))
      betaj(3)=ad_beta(1+1,i(1))+ad_beta(2+1,i(2))+ad_beta(3+1,i(3))

      betaj(4)=ad_beta(0+1,i(1))+ad_beta(1+1,i(2))
     $        +ad_beta(2+1,i(3))+ad_beta(3+1,i(3))

      betaj(5)=ad_beta(0+1,i(2))+ad_beta(1+1,i(3))
     $        +ad_beta(2+1,i(3))+ad_beta(3+1,i(3))

      betaj(6)=ad_beta(0+1,i(3))+ad_beta(1+1,i(3))+ad_beta(2+1,i(3))
      betaj(7)=ad_beta(0+1,i(3))+ad_beta(1+1,i(3))
      betaj(8)=ad_beta(0+1,i(3))

      call cmult(betaj,1./(ad_dt*(ad_nsteps-navg_step+1)),8)

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
         call copy(wk4,tb(1,i,1),n)
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
     $               ,tb(1,i,1),zeros,zeros)
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
     $               ,zeros,tb(1,i,1),zeros)
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
            call convect_new(riesz_rt(1,l2),tb(1,i,1),.false.,
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

      ! TODO: Remove and put ifdebug piece somewhere else

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

      ! This subroutine was made for debugging steady Boussinesq

      ! TODO: Move to somewhere else

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

      ! This subroutine was made for debugging steady Boussinesq

      ! TODO: Move to somewhere else

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

      ! This subroutine was made for debugging steady Boussinesq

      ! TODO: Move to somewhere else

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

      ! This subroutine was made for debugging steady Boussinesq

      ! TODO: Move to somewhere else

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
      subroutine steady_stokes_solve(ehu,ehv,ehw,ehp,rhs1,rhs2,rhs3)

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

      if (nio.eq.0) write (6,*) 'inside steady_stokes solver'

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

      if (nio.eq.0) write (6,*) 'exiting steady_stokes solver'

      return
      end
c-----------------------------------------------------------------------
      subroutine crd_divf(num_ts)

      ! make representation to be divergence free
      ! This subroutine was made for debugging steady Boussinesq

      ! TODO: Move to somewhere else

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
         call steady_stokes_solve(eh_u(1,1),eh_u(1,2),
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

c     Compute ROM residual for unsteady NS with Boussinesq

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
     $               ,tb(1,i,1),zeros,zeros)
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
     $               ,zeros,tb(1,i,1),zeros)
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
         call copy(riesz_rt(1,l2),tb(1,i,1),n)
         call col2(riesz_rt(1,l2),bm1,n)
         call chsign(riesz_rt(1,l2),n)
         l2=l2+1
      enddo
      if (nid.eq.0) write(6,*)l2,'lres_t_1'
      do i=0,nb
         ifield=2
         call copy(wk4,tb(1,i,1),n)
         call axhelm(riesz_rt(1,l2),wk4,ones,zeros,1,1)
         call cmult(riesz_rt(1,l2),param(8),n)
         call chsign(riesz_rt(1,l2),n) 
         l2=l2+1
      enddo
      if (nid.eq.0) write(6,*)l2,'lres_t_2'
      do j=0,nb
         do i=0,nb
         ifield=2
            call convect_new(riesz_rt(1,l2),tb(1,i,1),.false.,
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
      subroutine set_theta_uns

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside set_theta_uns'

      if (ifbuoy) call exitti('ifbuoy temporarily disabled in EI',1) 

      call rzero(theta_u,lres_u)
      call rzero(theta_t,lres_t)

c     call set_betaj
c     call set_alphaj
      call set_betaj_new
      call set_alphaj_new

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
            theta_u(l1)=-sin(bu_angle)*ad_ra*(theta_u(l1)+uta_ext(i))
            l1=l1+1
         enddo
         call mxm(utj,(nb+1),alphaj,6,theta_u(l1),1)
         do i=0,nb
            theta_u(l1)=-cos(bu_angle)*ad_ra*(theta_u(l1)+uta_ext(i))
            l1=l1+1
         enddo
      endif

      ! convection (done)
      ! Order of u2j does not matter since it is symmetric
      call mxm(u2j,(nb+1)**2,alphaj,6,theta_u(l1),1)
      do j=0,nb
      do i=0,nb
         theta_u(l1)=theta_u(l1)+u2a_ext(1+i+(nb+1)*j)
         l1=l1+1
      enddo
      enddo
      

      if (l1.ne.nres_u+1) then
         call exitti('theta_u term missing...$',l1)
      endif
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
         theta_t(l2)=theta_t(l2)+utua_ext(1+i+(nb+1)*j)
         l2=l2+1
      enddo
      enddo

      if (l2.ne.nres_t+1) then
         call exitti('theta_t term missing...$',l2)
      endif
      do i=1,nres_t
         if (nio.eq.0) write (6,*) theta_t(i),'theta_t'
      enddo

      if (nio.eq.0) write (6,*) 'exitng set_theta_uns'

      return
      end
c-----------------------------------------------------------------------
      subroutine cres_uns(aa,bb)

      ! This subroutine has been merged to cdres

      ! TODO: Verify and remove

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
      if (ifbuoy) call exitti('ifbuoy temporarily disabled in EI',1) 

      call nekgsync

c     call set_theta_uns
      call set_betaj_new
      call set_alphaj_new

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
         if (nid.eq.0) write(6,*)coef(i),'theta_u'
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
         if (nid.eq.0) write(6,*)coef(i),'theta_u'
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
         call opcopy(wk1,wk2,wk3,tb(1,i,1),zeros,zeros)
         call opcolv(wk1,wk2,wk3,bm1)
         call opchsgn(wk1,wk2,wk3)

         coef(i)=-sin(bu_angle)*ad_ra*(coef(i)+uta_ext(i))
         call cfill(wk4,coef(i),n)
         if (nid.eq.0) write(6,*)coef(i),'theta_u'
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
         call opcopy(wk1,wk2,wk3,zeros,tb(1,i,1),zeros)
         call opcolv(wk1,wk2,wk3,bm1)
         call opchsgn(wk1,wk2,wk3)

         coef(i)=-cos(bu_angle)*ad_ra*(coef(i)+uta_ext(i))
         call cfill(wk4,coef(i),n)
         if (nid.eq.0) write(6,*)coef(i),'theta_u'
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

            coef2(i+(nb+1)*j)=coef2(i+(nb+1)*j)
     $                         +u2a_ext(1+i+(nb+1)*j)
            if (nid.eq.0) write(6,*)coef2(i+(nb+1)*j),'theta_u'
            call cfill(wk4,coef2(i+(nb+1)*j),n)
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
         call copy(wk1,tb(1,i,1),n)
         call col2(wk1,bm1,n)
         call chsign(wk1,n)

         call cfill(wk2,coef(i),n)
         call add2col2(res_t,wk1,wk2,n)
         l2=l2+1
      enddo
      if (nid.eq.0) write(6,*)l2,'lres_t_1'
      do i=0,nb
         ifield=2
         call copy(wk2,tb(1,i,1),n)
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
            call convect_new(wk1,tb(1,i,1),.false.,
     $                       ub(1,j),vb(1,j),wb(1,j),.false.)
            call chsign(wk1,n) 

            coef2(i+(nb+1)*j)=coef2(i+(nb+1)*j)
     $                         +utua_ext(1+i+(nb+1)*j)
            call cfill(wk2,coef2(i+(nb+1)*j),n)
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
      call steady_stokes_solve(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
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
c-----------------------------------------------------------------------
      subroutine set_alphaj_new

      include 'SIZE'
      include 'MOR'

      ! ad_alpha(3,3)

      if (navg_step.eq.1) then
         alphaj(1)=ad_alpha(1,1)+ad_alpha(2,2)-1
         alphaj(2)=ad_alpha(1,2)-ad_alpha(1,3)
         alphaj(3)=0.
      elseif (navg_step.eq.2) then
         alphaj(1)=ad_alpha(2,2)-ad_alpha(1,3)-ad_alpha(2,3)
         alphaj(2)=ad_alpha(1,2)-ad_alpha(1,3)
         alphaj(3)=0.
      elseif (navg_step.ge.3) then
         alphaj(1)=-ad_alpha(1,3)-ad_alpha(2,3)
         alphaj(2)=-ad_alpha(1,3)
         alphaj(3)=0.
      endif
      alphaj(4)=-ad_alpha(3,3)
      alphaj(5)=-ad_alpha(2,3)-ad_alpha(3,3)
      alphaj(6)=-1.

      call cmult(alphaj,1./(1.*(ad_nsteps-navg_step+1)),6)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_betaj_new

      include 'SIZE'
      include 'MOR'

      ! ad_beta(4,3)

      if (navg_step.eq.1) then
         betaj(1)=ad_beta(1+1,1)+ad_beta(2+1,2)+ad_beta(3+1,3)
         betaj(2)=ad_beta(0+1,1)+ad_beta(1+1,2)
     $           +ad_beta(2+1,3)+ad_beta(3+1,3)
         betaj(3)=ad_beta(0+1,2)+ad_beta(1+1,3)
     $           +ad_beta(2+1,3)+ad_beta(3+1,3)
      elseif (navg_step.eq.2) then
         betaj(1)=ad_beta(2+1,2)+ad_beta(3+1,3)
         betaj(2)=ad_beta(1+1,2)+ad_beta(2+1,3)+ad_beta(3+1,3)
         betaj(3)=ad_beta(0+1,2)+ad_beta(1+1,3)
     $           +ad_beta(2+1,3)+ad_beta(3+1,3)
      elseif (navg_step.ge.3) then
         betaj(1)=ad_beta(3+1,3)
         betaj(2)=ad_beta(2+1,3)+ad_beta(3+1,3)
         betaj(3)=ad_beta(1+1,3)+ad_beta(2+1,3)+ad_beta(3+1,3)
      endif
      betaj(4)=ad_beta(0+1,3)+ad_beta(1+1,3)+ad_beta(2+1,3)
      betaj(5)=ad_beta(0+1,3)+ad_beta(1+1,3)
      betaj(6)=ad_beta(0+1,3)

      call cmult(betaj,1./(ad_dt*(ad_nsteps-navg_step+1)),6)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_rom_residual

      ! This subroutine is used when lei=0 and eqn = NS

      ! compute the rom residual
      ! make res_u and res_t in common block

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real coef(0:nb)
      real coef2(0:(nb+1)**2-1)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside set_rom_residual'

      call nekgsync
c     call set_theta_uns
      call set_betaj
      call set_alphaj

      if (ifrom(1)) call opzero(res_u(1,1),res_u(1,2),res_u(1,ldim))
      if (ifrom(2)) call rzero(res_t,n)

      ! setup velocity residual
      if (ifrom(1)) then
         if (iftran) then
            call resid_in_time('vel')
         endif

         call resid_in_diffusion('vel')

         if (ifbuoy) call resid_in_buoy

         if (ifadvc(1)) then
            call resid_in_advec('vel')
         endif
      endif

      ! setup temperature residual
      if (ifrom(2)) then
         if (iftran) call resid_in_time('tmp')

         call resid_in_diffusion('tmp')

         if (ifadvc(2)) then
            call resid_in_advec('tmp')
         endif

         if (ifsource) then
            call resid_in_source
         endif
      endif

      if (nio.eq.0) write (6,*) 'exiting set_rom_residual'
      return
      end
c-----------------------------------------------------------------------
      subroutine get_dual_norm

      ! compute dual norm of the residual

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real u_resid_h1_norm
      real u_resid_h10_norm
      real u_resid_l2_norm
      real t_resid_h1_norm
      real t_resid_h10_norm
      real t_resid_l2_norm
      real resid_h1_norm
      real resid_h10_norm
      real resid_l2_norm

      if (nio.eq.0) write (6,*) 'inside c_dual_norm'

      u_resid_h10_norm = 0
      u_resid_l2_norm = 0
      t_resid_h10_norm = 0
      t_resid_l2_norm = 0

      if (ifrom(1)) then
         ifield=1
         u_resid_h10_norm = h10vip(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
     $                             eh_u(1,1),eh_u(1,2),eh_u(1,ldim))
         u_resid_l2_norm = wl2vip(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
     $                            eh_u(1,1),eh_u(1,2),eh_u(1,ldim))
         u_resid_h1_norm = u_resid_h10_norm + u_resid_l2_norm
         if (nid.eq.0) write(6,*)'vel residual in h1 norm',
     $                            sqrt(u_resid_h1_norm)
         if (nid.eq.0) write(6,*)'vel residual in h10 norm',
     $                            sqrt(u_resid_h10_norm)
         if (nid.eq.0) write(6,*)'vel residual in l2 norm',
     $                            sqrt(u_resid_l2_norm)
      endif

      if (ifrom(2)) then
         ifield=2
         t_resid_h10_norm = h10sip(eh_t,eh_t)
         t_resid_l2_norm = wl2sip(eh_t,eh_t)

         t_resid_h1_norm = t_resid_h10_norm + t_resid_l2_norm
         if (nid.eq.0) write(6,*)'temp residual in h1 norm',
     $                            sqrt(t_resid_h1_norm)
         if (nid.eq.0) write(6,*)'temp residual in h10 norm',
     $                            sqrt(t_resid_h10_norm)
         if (nid.eq.0) write(6,*)'temp residual in l2 norm',
     $                            sqrt(t_resid_l2_norm)
      endif

      resid_h1_norm = sqrt(u_resid_h10_norm+u_resid_l2_norm+
     $                     t_resid_h10_norm+t_resid_l2_norm)
      resid_h10_norm = sqrt(u_resid_h10_norm+t_resid_h10_norm)
      resid_l2_norm = sqrt(u_resid_l2_norm+t_resid_l2_norm)

      if (resid_h1_norm.le.0) call
     $             exitti('negative semidefinite dual norm$',n)

      if (nid.eq.0) write (6,*) 'residual in h1 norm:',resid_h1_norm
      if (nid.eq.0) write (6,*) 'residual in h10 norm:',resid_h10_norm
      if (nid.eq.0) write (6,*) 'residual in l2 norm:',resid_l2_norm

      call nekgsync
      if (nio.eq.0) write (6,*) 'exiting c_dual_norm'

      return
      end
c-----------------------------------------------------------------------
      subroutine resid_in_time(msg)

      ! compute the rom residual in time derivative term

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /screi/ wk1(lt),wk2(lt),wk3(lt),wk4(lt)

      real coef(0:nb)
      character*3 msg

      if (nio.eq.0) write (6,*) 'inside resid_in_time'

      n=lx1*ly1*lz1*nelv

      if (msg.eq.'vel') then
         call mxm(uj,nb+1,betaj,8,coef,1)
         do i=0,nb
            ifield=1
            call opcopy(wk1,wk2,wk3,ub(1,i),vb(1,i),wb(1,i))
            call opcolv(wk1,wk2,wk3,bm1)
            call opchsgn(wk1,wk2,wk3)

            call cfill(wk4,coef(i),n)
            if (nid.eq.0) write(6,*)coef(i),'theta_u'
            call add2col2(res_u(1,1),wk1,wk4,n)
            call add2col2(res_u(1,2),wk2,wk4,n)
            if (ldim.eq.3) then
               call add2col2(res_u(1,ldim),wk3,wk4,n)
            endif
         enddo
      elseif (msg.eq.'tmp') then
         call mxm(utj,nb+1,betaj,8,coef,1)
         do i=0,nb
            ifield=2
            call copy(wk1,tb(1,i,1),n)
            call col2(wk1,bm1,n)
            call chsign(wk1,n)

            if (nid.eq.0) write(6,*)coef(i),'theta_t'
            call cfill(wk2,coef(i),n)
            call add2col2(res_t,wk1,wk2,n)
         enddo
      endif

      if (nio.eq.0) write (6,*) 'exitting resid_in_time'
      return
      end
c-----------------------------------------------------------------------
      subroutine resid_in_diffusion(msg)

      ! compute the rom residual in diffusion term

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /screi/ wk1(lt),wk2(lt),wk3(lt),wk4(lt),wk5(lt),wk6(lt)

      real coef(0:nb)
      character*3 msg

      if (nio.eq.0) write (6,*) 'inside resid_diffusion'

      n=lx1*ly1*lz1*nelv

      if (msg.eq.'vel') then
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
            if (nid.eq.0) write(6,*)coef(i),'theta_u'
            call add2col2(res_u(1,1),wk1,wk4,n)
            call add2col2(res_u(1,2),wk2,wk4,n)
            if (ldim.eq.3) then
               call add2col2(res_u(1,ldim),wk3,wk4,n)
            endif
         enddo
      elseif (msg.eq.'tmp') then
         do i=0,nb
            ifield=2
            call copy(wk2,tb(1,i,1),n)
            call axhelm(wk1,wk2,ones,zeros,1,1)
            call cmult(wk1,param(8),n)
            call chsign(wk1,n)

            coef(i)=uta(i)
            if (nid.eq.0) write(6,*)coef(i),'theta_t'
            call cfill(wk2,coef(i),n)
            call add2col2(res_t,wk1,wk2,n)
         enddo
         if (ifhelm) then
            if(nid.eq.0) write(6,*)'helm is on'
            do i=0,nb
               ifield=2
               call copy(wk1,tb(1,i,1),n)
               call col2(wk1,bm1,n)
               call chsign(wk1,n)

               coef(i) = uta(i) * ad_mu
               if (nid.eq.0) write(6,*)coef(i),'theta_t'
               call cfill(wk2,coef(i),n)
               call add2col2(res_t,wk1,wk2,n)
            enddo
         endif
      endif

      if (nio.eq.0) write (6,*) 'exitting resid_diffusion'
      return
      end
c-----------------------------------------------------------------------
      subroutine resid_in_advec(msg)

      ! compute the rom residual in advection term

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /screi/ wk1(lt),wk2(lt),wk3(lt),wk4(lt)

      real coef2(0:(nb+1)**2-1)
      character*3 msg

      if (nio.eq.0) write (6,*) 'inside resid_advec'

      n=lx1*ly1*lz1*nelv

      if (msg.eq.'vel') then
         call mxm(u2j,(nb+1)**2,alphaj,8,coef2,1)
         do j=0,nb
            do i=0,nb
            ifield=1
               call convect_new(wk1,ub(1,i),.false.,
     $                          ub(1,j),vb(1,j),wb(1,j),.false.)
               call convect_new(wk2,vb(1,i),.false.,
     $                          ub(1,j),vb(1,j),wb(1,j),.false.)
               if (ldim.eq.3) then
                  call convect_new(wk3,wb(1,i),.false.,
     $                             ub(1,j),vb(1,j),wb(1,j),.false.)
               endif
               call opchsgn(wk1,wk2,wk3)

               coef2(i+(nb+1)*j)=coef2(i+(nb+1)*j)
     $                            +u2a(1+i+(nb+1)*j)
               if (nid.eq.0) write(6,*)coef2(i+(nb+1)*j),'theta_u'
               call cfill(wk4,coef2(i+(nb+1)*j),n)
               call add2col2(res_u(1,1),wk1,wk4,n)
               call add2col2(res_u(1,2),wk2,wk4,n)
               if (ldim.eq.3) then
                  call add2col2(res_u(1,ldim),wk3,wk4,n)
               endif
            enddo
         enddo
      elseif (msg.eq.'tmp') then
         call mxm(utuj,(nb+1)**2,alphaj,8,coef2,1)
         do j=0,nb
            do i=0,nb
               ifield=2
               call convect_new(wk1,tb(1,i,1),.false.,
     $                          ub(1,j),vb(1,j),wb(1,j),.false.)
               call chsign(wk1,n)

               coef2(i+(nb+1)*j)=coef2(i+(nb+1)*j)
     $                            +utua(1+i+(nb+1)*j)
               call cfill(wk2,coef2(i+(nb+1)*j),n)
               if (nid.eq.0) write(6,*)coef2(i+(nb+1)*j),'theta_t'
               call add2col2(res_t,wk1,wk2,n)
            enddo
         enddo
      endif

      if (nio.eq.0) write (6,*) 'exitting resid_advec'
      return
      end
c-----------------------------------------------------------------------
      subroutine resid_in_buoy

      ! compute the rom residual in buoyancy term

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /screi/ wk1(lt),wk2(lt),wk3(lt),wk4(lt)

      real coef(0:nb)
      character*3 msg

      if (nio.eq.0) write (6,*) 'inside resid_buoy'
      if (ifbuoy) call exitti('ifbuoy temporarily disabled in EI',1) 

      n=lx1*ly1*lz1*nelv

      call mxm(utj,(nb+1),alphaj,8,coef,1)
      do i=0,nb
         ifield=1
         call opcopy(wk1,wk2,wk3,tb(1,i,1),zeros,zeros)
         call opcolv(wk1,wk2,wk3,bm1)
         call opchsgn(wk1,wk2,wk3)

         coef(i)=-sin(bu_angle)*ad_ra*(coef(i)+uta(i))
         call cfill(wk4,coef(i),n)
         if (nid.eq.0) write(6,*)coef(i),'theta_u'
         call add2col2(res_u(1,1),wk1,wk4,n)
         call add2col2(res_u(1,2),wk2,wk4,n)
         if (ldim.eq.3) then
            call add2col2(res_u(1,ldim),wk3,wk4,n)
         endif
      enddo

      call mxm(utj,(nb+1),alphaj,8,coef,1)
      do i=0,nb
         ifield=1
         call opcopy(wk1,wk2,wk3,zeros,tb(1,i,1),zeros)
         call opcolv(wk1,wk2,wk3,bm1)
         call opchsgn(wk1,wk2,wk3)

         coef(i)=-cos(bu_angle)*ad_ra*(coef(i)+uta(i))
         call cfill(wk4,coef(i),n)
         if (nid.eq.0) write(6,*)coef(i),'theta_u'
         call add2col2(res_u(1,1),wk1,wk4,n)
         call add2col2(res_u(1,2),wk2,wk4,n)
         if (ldim.eq.3) then
            call add2col2(res_u(1,ldim),wk3,wk4,n)
         endif
      enddo

      if (nio.eq.0) write (6,*) 'exitting resid_buoy'
      return
      end
c-----------------------------------------------------------------------
      subroutine resid_in_source

      ! compute the rom residual in source term

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /screi/ wk1(lt),wk2(lt),wk3(lt),wk4(lt)

      real coef(0:nb)
      real t1, t2, t3, t4
      character*3 msg

      if (nio.eq.0) write (6,*) 'inside resid_in_source'

      n=lx1*ly1*lz1*nelv

      t1 = glsc2(alphaj,rqj,6)
      t2 = t1 + rqa
      write(6,*)'t2',t2
      call add2(res_t,qq,n)

      if (ifsrct) then
         t3 = glsc2(alphaj,rqtj,6)
         t4 = t3 + rqta
         write(6,*)'t4',t4,t3,rqta
         call add2s2(res_t,qqxyz,t4,n)
      endif

      if (nio.eq.0) write (6,*) 'exitting resid_in_source'
      return
      end
c-----------------------------------------------------------------------
      subroutine cdres

      ! Computes the dual norm of the residual
      ! which is the norm of its riesz representation(s)

      include 'SIZE'
      include 'MOR'

      if (lei.eq.0) then
c        call c_rieszrd_uns
         call set_rom_residual
         call set_riesz_one
c        call get_dual_norm
      else
         if (eqn.eq.'POIS') then
            call set_theta_poisson
         else if (eqn.eq.'HEAT') then
            call set_theta_heat
         else if (eqn.eq.'ADVE') then
            call set_theta_ad
         else if (eqn.eq.'VNS ') then
            call set_theta_ns
         else if (eqn.eq.'NS  ') then
            call set_theta_uns
         else if (eqn.eq.'SNS ') then
            call exitti('no supprot for Steady NS ROM$',eqn)
         endif
      endif

      call c_riesz_norm

      return
      end
