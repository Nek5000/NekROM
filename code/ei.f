c-----------------------------------------------------------------------
      subroutine cres

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /eires/ xi(lt,lres),theta(lres),sigma(lres,lres)
      common /eiivar/ nres
      common /eivar/ res

      if (eqn.eq.'POI') then
         call set_theta_poisson
      else if (eqn.eq.'HEA') then
         call set_theta_heat
      endif

      res=0.

      do j=1,nres
      do i=1,nres
         res=res+sigma(i,j)*theta(i)*theta(j)
      enddo
      enddo

      if (res.le.0) call exitti('negative semidefinite residual$',n)

      res=sqrt(res)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_xi_poisson

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /eires/ xi(lt,lres),theta(lres),sigma(lres,lres)
      common /eiivar/ nres

      n=lx1*ly1*lz1*nelv

      if (ifield.eq.1) then
         call exitti('(set_ritz_a) ifield.eq.1 not supported...$',nb)
      else
         if (ips.eq.'L2 ') then
            do i=1,nb
               call axhelm(xi(1,i),tb(1,i),ones,zeros,1,1)
               call binv1(xi(1,i))
            enddo
            call copy(xi(1,nb+1),qq,n)
            call binv1(xi(1,nb+1))
         else
            do i=0,nb
               call copy(xi(1,i),tb(1,i),n)
            enddo
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

      common /eires/ xi(lt,lres),theta(lt,lres),sigma(lres,lres)
      common /eiivar/ nres

      n=lx1*ly1*lz1*nelv

      l=1
      if (ifield.eq.1) then
         call exitti('(set_ritz_a) ifield.eq.1 not supported...$',nb)
      else
         if (ips.eq.'L2 ') then
            do i=0,nb
               call copy(xi(l,i),tb(1,i),n)
               l=l+1
            enddo
            do i=0,nb
               call axhelm(xi(l,i),tb(1,i),ones,zeros,1,1)
               call binv1(xi(l,i))
               l=l+1
            enddo
            call copy(xi(l,nb+2),bqr,n)
            call binv1(xi(l,nb+2))
            l=l+1
         else
            do i=0,nb
               call copy(xi(1,i),tb(1,i),n)
            enddo
         endif
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
      common /eires/ xi(lt,lres),theta(lres),sigma(lres,lres)
      common /eiivar/ nres

      n=lx1*ly1*lz1*nelv

      if (eqn.eq.'POI') then
         nres=nb+1
      else if (eqn.eq.'HEA') then
         nres=(nb+1)*2+1
      endif

      if (nres.gt.lres) call exitti('nres > lres$',nres)
      if (nres.le.0) call exitti('nres <= 0$',nres)

      if (eqn.eq.'POI') then
         call set_xi_poisson
      else if (eqn.eq.'HEA') then
         call set_xi_heat
      endif

      do i=1,nres
      do j=1,nres
         sigma(i,j)=glsc3(xi(1,i),xi(1,j),bm1,n)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_theta_poisson

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /eires/ xi(lt,lres),theta(lres),sigma(lres,lres)
      common /eiivar/ nres

      n=lx1*ly1*lz1*nelv

      do i=1,nb
         theta(i)=ut(i,1)
      enddo

      theta(nb+1)=-1.

      return
      end
c-----------------------------------------------------------------------
      subroutine set_theta_heat

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /eires/ xi(lt,lres),theta(lres),sigma(lres,lres),betaj(6)
      common /eiivar/ nres

      n=lx1*ly1*lz1*nelv

      l=1
      call set_betaj(betaj)
      call mxm(uj,lub+1,betaj,6,theta(l),1)
      l=l+lub+1
      do i=0,nb
         theta(l)=uta(i)
         l=l+1
      enddo

      theta(nb+1)=-1.
      l=l+1

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_poisson

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /eires/ xi(lt,lres),theta(lres),sigma(lres,lres)
      common /eiivar/ nres

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

      nres=nb+1
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

      common /scrrhs/ tmp(0:nb)

      real rhs(nb)

      n=lx1*ly1*lz1*nelv

      do i=1,nb
         rhs(i)=glsc2(bqr,tb(1,i),n)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_betaj(betaj)

      include 'SIZE'
      include 'MOR'

      real betaj(6)

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
