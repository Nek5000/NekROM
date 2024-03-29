      subroutine IPM(rhs,vv,helm,invhelm,amax,amin,adis,bpar,bstep,
     $               tol_box,ifdiag)

      ! Interrior-point method (IPM)
      ! Transform the inequality constrained optimization problem
      ! into unconstrained optimization problem with barrier method
      ! Quasi-Newton is used to solve the unconstrained opt problem

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      common /scripm/ tlncount

      real vv(nb),rhs(nb)
      real helm(nb,nb),invhelm(nb,nb)
      real amax(nb),amin(nb),adis(nb)

      real uu(nb)
      real ngf,ysk,norm_step
      real bpar,par,tol_box
      integer bstep,tlncount,qstep

      logical ifdiag

      call copy(uu,vv,nb)

      par = bpar 

      tcopt_time=dnekclock()

      do k=1,bstep

         tquasi_time=dnekclock()
         call quasi_newton(uu,qstep,ngf,norm_step,ysk,
     $                     helm,invhelm,rhs,amax,amin,
     $                     par,tol_box,ifdiag)
         quasi_time=quasi_time+dnekclock()-tquasi_time

         if (mod(ad_step,ad_iostep).eq.0) then
            if (nio.eq.0) write (6,*) 'lnconst_ana'
            call cpod_ana(uu,par,qstep,tlncount,ngf,norm_step
     $      ,ysk)
         endif
         par = par*0.1

      enddo

      copt_time=copt_time+dnekclock()-tcopt_time

      call copy(rhs,uu,nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine quasi_newton(uu,qstep,ngf,norm_step,ysk,
     $                        helm,invhelm,rhs,amax,amin,
     $                        par,tol_box,ifdiag)

      ! BFGS Method

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      common /scripm/ tlncount

      real vv(nb),rhs(nb)
      real helm(nb,nb),invhelm(nb,nb)
      real amax(nb),amin(nb),adis(nb)

      real sk(nb,nb),yk(nb,nb)
      real qgo(nb),qngradf(nb)
      real qns(nb),qny(nb)
      real tmp(nb),uu(nb)
      real qnf,ngf,ysk
      real par,tol_box
      real norm_s,norm_step,norm_uo

      ! parameter for barrier function
      integer chekbc,qstep
      integer lncount,tlncount

      logical ifdiag

      chekbc = 0
      tlncount = 0
      ! use helm from BDF3/EXT3 as intial approximation
      call comp_qngradf(uu,rhs,helm,qngradf,amax,amin,
     $                  par,ifdiag,barr_func,nb)

      norm_uo = vlamax(uu,nb)

      do j=1,nb
         if (j.eq.1) then
            call copy(qns,qngradf,nb)
            if (ifdiag) then
               call col2(qns,invhelm,nb)
            else 
               call mxm(invhelm,nb,qngradf,nb,qns,1)
            endif
            call chsign(qns,nb)
         endif

         tlnsrch_time=dnekclock()
         call backtrackr(uu,qns,rhs,qngradf,helm,invhelm,1e-4,0.5,
     $                   amax,amin,par,chekbc,lncount,
     $                   tol_box,ifdiag)
         lnsrch_time=lnsrch_time+dnekclock()-tlnsrch_time
         tlncount = tlncount + lncount

         ! store qns
         call copy(sk(1,j),qns,nb)

         ! store old qn-gradf
         call copy(qgo,qngradf,nb) 

         ! update qn-gradf
         tcompgf_time=dnekclock()
         call comp_qngradf(uu,rhs,helm,qngradf,amax,amin,
     $                     par,ifdiag,barr_func,nb)
         compgf_time=compgf_time+dnekclock()-tcompgf_time
         call sub3(qny,qngradf,qgo,nb) 

         ysk = vlsc2(qny,qns,nb)

         norm_s = vlamax(qns,nb)
         norm_step = norm_s/norm_uo

         ! compute norm of gradf
         ngf = vlamax(qngradf,nb)

         ! store qny 
         call copy(yk(1,j),qny,nb)

         tinvhm_time=dnekclock()
         call invH_multiply(qns,invhelm,sk,yk,qngradf,j,ifdiag,nb)
         invhm_time=invhm_time+dnekclock()-tinvhm_time

         if (ngf .lt. 1e-6 .OR. ysk .lt. 1e-6 .OR. norm_step .lt.
     $      1e-6  ) then 
            exit
         endif

      enddo

      qstep=j

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_qngradf(uu,rhs,helm,s,amax,amin,
     $                        bpar,ifdiag,barr_func,nb)
      
      real helm(nb,nb)
      real uu(nb),rhs(nb),s(nb)
      real amax(nb),amin(nb) 
      real tmp1(nb),tmp2(nb),tmp3(nb),tmp4(nb)
      real denum(nb)
      real bpar,mpar,pert
      logical ifdiag
      integer nb,barr_func

      pert = 1e-2

      if (barr_func.eq.1) then ! use logarithmic as barrier function
         
         pert = 1e-2

         call sub3(tmp1,amax,uu,nb)  
         call sub3(tmp2,uu,amin,nb)  
   
         ! add perturbation 
         call cadd(tmp1,pert,nb)
         call cadd(tmp2,pert,nb)
   
         call sub3(tmp3,tmp1,tmp2,nb)
         call col3(denum,tmp1,tmp2,nb)
         call invcol2(tmp3,denum,nb)

         mpar = -1.0*bpar
         call add3s12(s,rhs,tmp3,-1.0,mpar,nb)

         if (ifdiag) then 
            call copy(tmp4,uu,nb)
            call col2(tmp4,helm,nb)
         else
            call mxm(helm,nb,uu,nb,tmp4,1)   
         endif
         call add2(s,tmp4,nb)

      else ! use inverse function as barrier function
         

         pert = -1e-2

         call sub3(tmp1,uu,amax,nb)  
         call sub3(tmp2,amin,uu,nb)  

         ! add perturbation
         call cadd(tmp1,pert,nb)
         call cadd(tmp2,pert,nb)

         call vsq(tmp1,nb)
         call vsq(tmp2,nb)

         call invcol1(tmp1,nb)
         call invcol1(tmp2,nb)
         call sub3(tmp3,tmp2,tmp1,nb)

c        call sub3(tmp1,amax,uu,nb)  
c        call sub3(tmp2,amin,uu,nb)  

c        ! add perturbation
c        call cadd(tmp1,1e-2,nb)
c        call cadd(tmp2,-1e-2,nb)

c        call sub3(tmp3,tmp1,tmp2,nb)
c        call add3(tmp4,tmp1,tmp2,nb)
c        call col2(tmp3,tmp4,nb)

c        call vsq(tmp1,nb)
c        call vsq(tmp2,nb)

c        call col3(denum,tmp1,tmp2,nb)
c        call invcol2(tmp3,denum,nb)

         mpar = -1.0*bpar
         call add3s12(s,rhs,tmp3,-1.0,mpar,nb)
   
         if (ifdiag) then 
            call copy(tmp4,uu,nb)
            call col2(tmp4,helm,nb)
         else
            call mxm(helm,nb,uu,nb,tmp4,1)   
         endif
         call add2(s,tmp4,nb)

      endif

      return 
      end
c-----------------------------------------------------------------------
      subroutine comp_qnf(uu,rhs,helm,invhelm,qnf,amax,amin,
     $                    bpar,ifdiag,barr_func,nb)
      
      real uu(nb),rhs(nb)
      real helm(nb,nb),invhelm(nb,nb)
      real amax(nb),amin(nb)
      real qnf,bpar
      real tmp1(nb),tmp2(nb),tmp3(nb)
      real tmp4(nb),tmp5(nb),tmp6(nb)
      real denum(nb)
      real bar1,bar2,pert
      real term1,term2,term3,term4
      integer nb,barr_func

      logical ifdiag

      ! evaluate quasi-newton f
      if (ifdiag) then
         ! 0.5*coef'*H*coef
         term1 = 0.5 * vlsc3(uu,helm,uu,nb)
         term2 = vlsc2(uu,rhs,nb) ! coef'*rhs
      else 
         ! 0.5*coef'*H*coef
         call mxm(helm,nb,uu,nb,tmp6,1)
         term1 = 0.5 * vlsc2(tmp6,uu,nb)
         term2 = vlsc2(uu,rhs,nb) ! coef'*rhs
      endif

      if (barr_func.eq.1) then ! use logarithmetic as barrier function

         pert=1e-2
         ! barrier term
         call sub3(tmp1,amax,uu,nb)  
         call sub3(tmp2,uu,amin,nb)  
   
         ! add perturbation
         call cadd(tmp1,pert,nb)
         call cadd(tmp2,pert,nb)
   
         do i=1,nb
            tmp3(i) = log(tmp1(i))
            tmp4(i) = log(tmp2(i))
         enddo

         bar1 = vlsum(tmp3,nb)
         bar2 = vlsum(tmp4,nb)
         term4 = bpar*(bar1+bar2)

      else 

         pert=-1e-2
         call sub3(tmp1,uu,amax,nb)  
         call sub3(tmp2,amin,uu,nb)  

         ! add perturbation
         call cadd(tmp1,pert,nb)
         call cadd(tmp2,pert,nb)

         call invers2(tmp3,tmp1,nb)
         call invers2(tmp4,tmp2,nb)

         bar1 = vlsum(tmp3,nb)
         bar2 = vlsum(tmp4,nb)
         term4 = bpar*(bar1+bar2)

      endif

      qnf = term1 - term2 - term4

      return
      end
c-----------------------------------------------------------------------
      subroutine Hessian_update(B,s,y,nb)

      real B(nb,nb)
      real rk1up(nb,nb),rk2up(nb,nb)
      real w1(nb,nb),w2(nb,nb)
      real s(nb),y(nb),w3(nb)
      real ys,sBs
      
      ! s_k * s_k^T               
      call mxm(s,nb,s,1,w1,nb)

      ! s_k * s_k^T * B_k
      call mxm(w1,nb,B,nb,w2,nb)

      ! B_k * s_k * s_k^T * B_k 
      call mxm(B,nb,w2,nb,rk1up,nb)

      ! s_k^T * B_k * s_k 
      call mxm(B,nb,s,nb,w3,1)

      sBs = vlsc2(s,w3,nb)

      call cmult(rk1up,-1.0/sBs,nb*nb)

      ! second rank-one update
      ! y_k * y_k^T               
      call mxm(y,nb,y,1,rk2up,nb)
      ! y_k^T * s_k
      ys = vlsc2(y,s,nb)

      call cmult(rk2up,1.0/ys,nb*nb)

      call add4(B(1,1),B(1,1),rk1up(1,1),rk2up(1,1),nb*nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine backtrackr(uu,s,rhs,qngradf,helm,invhelm,sigmab,facb,
     $            amax,amin,bpar,chekbc,counter,tol_box,ifdiag)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real rhs(nb),s(nb)
      real uuo(nb),uu(nb)
      real helm(nb,nb),invhelm(nb,nb)
      real Jfk(nb),Jfk1(nb),qngradf(nb)
      real amax(nb),amin(nb)
      real fk,fk1
      real Jfks,Jfks1
      real sigmab,facb,alphak
      real bpar,minalpha,tol_box

      integer chekbc,counter,nb
      logical cond1,cond2,cond3,ifdiag

c     alphak = 1.0
      counter = 0

      tcompf_time=dnekclock()
      call comp_qnf(uu,rhs,helm,invhelm,fk,amax,amin,
     $              bpar,ifdiag,barr_func,nb) ! get old f
      compf_time=compf_time+dnekclock()-tcompf_time

      call copy(Jfk,qngradf,nb)

      call findminalpha(minalpha,s,uu,amax,amin,nb)

      call copy(uuo,uu,nb)

      if ((minalpha-1.0).gt.1e-10) then
         alphak = 1.0
      else
         alphak = minalpha
      endif

      call add2s2(uu,s,alphak,nb)

      call check_box(chekbc,uu,amax,amin,tol_box,nb)

      tcompf_time=dnekclock()
      call comp_qnf(uu,rhs,helm,invhelm,fk1,amax,amin,
     $              bpar,ifdiag,barr_func,nb) ! get new f
      compf_time=compf_time+dnekclock()-tcompf_time

      tcompgf_time=dnekclock()
      call comp_qngradf(uu,rhs,helm,Jfk1,amax,amin,
     $                  bpar,ifdiag,barr_func,nb)
      compgf_time=compgf_time+dnekclock()-tcompgf_time

      Jfks  = vlsc2(Jfk,s,nb)   
      Jfks1 = vlsc2(Jfk1,s,nb)   

      cond1 = (fk1-(fk+sigmab*alphak*Jfks)).gt.1e-10
      cond2 = (Jfks1-(0.9*Jfks)).lt.1e-10

      do while (cond1.OR.(chekbc.eq.1))
         counter = counter + 1

         alphak = alphak * facb

         call add3s2(uu,uuo,s,1.0,alphak,nb)

         call check_box(chekbc,uu,amax,amin,tol_box,nb)

         tcompf_time=dnekclock()
         call comp_qnf(uu,rhs,helm,invhelm,fk1,amax,amin,
     $                 bpar,ifdiag,barr_func,nb)
         compf_time=compf_time+dnekclock()-tcompf_time

         tcompgf_time=dnekclock()
         call comp_qngradf(uu,rhs,helm,Jfk1,amax,amin,
     $                     bpar,ifdiag,barr_func,nb)
         compgf_time=compgf_time+dnekclock()-tcompgf_time

         Jfks1 = vlsc2(Jfk1,s,nb)   

         cond1 = (fk1-(fk+sigmab*alphak*Jfks)).gt.1e-10
         cond2 = (Jfks1-(0.9*Jfks)).lt.1e-10

         cond3 = (alphak-minalpha).lt.1e-10
c        if ((cond3.AND..not.cond1).OR.alphak.lt.1e-8) then
         if ((.not.cond1).OR.alphak.lt.1e-8) then
            exit
         endif
      enddo

      call cmult(s,alphak,nb)

      if (mod(ad_step,ad_iostep).eq.0) then
         if (nio.eq.0) write(6,1)
     $         ad_step,bpar,' #lnsrch:',counter,' step length:',
     $         alphak,minalpha,chekbc,cond1,cond2
      endif

    1 format (i6,1p1e16.8,a9,i3,a13,1p2e16.8,i3,' ',2L2) 

      return
      end
c-----------------------------------------------------------------------
      subroutine cpod_ana(uu,par,qstep,tlncount,
     $                    ngf,qndf,ysk)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real uu(nb)
      real par
      real ngf,qndf,ysk
      integer qstep,tlncount

      if (nio.eq.0) then
         write (6,2)'ad_step:',ad_step,scopt,par,
     $            qstep,tlncount,
     $            ngf,qndf,ysk
         if (ad_step.eq.ad_nsteps) then
            do j=1,nb
               write(6,*) j,uu(j),'final'
            enddo
         else
            do j=1,nb
               write(6,*) j,uu(j)
            enddo
         endif
      endif
    2 format (a8,i6,1x,a5,1p1e16.8,i3,i5,1p4e16.8)  
      return
      end
c-----------------------------------------------------------------------
      subroutine invH_multiply(qnsol,invh0,sk,yk,qnd,qnstep,ifdiag,nb)

      include 'SIZE'
      include 'TOTAL'

      real qnsol(nb),invh0(nb)
      real sk(nb,nb),yk(nb,nb),qnd(nb)
      integer qnstep,nb

      real qnrho(nb),qnalpha(nb),qnbeta(nb)
      real tmp(nb)
      real qnfact(nb)
      real work(nb)
      logical ifdiag

      call copy(qnsol,qnd,nb)

      ! compute right product
      do i=qnstep,1,-1
         qnrho(i) = vlsc2(yk(1,i),sk(1,i),nb)
         qnalpha(i) = vlsc2(sk(1,i),qnsol,nb)/qnrho(i)
         call add2s2(qnsol,yk(1,i),-qnalpha(i),nb)
      enddo

      ! compute center
      if (ifdiag) then
         call col2(qnsol,invh0,nb)
      else
         call copy(tmp,qnsol,nb)
         call mxm(invh0,nb,tmp,nb,qnsol,1)
      endif

      ! compute left product
      do i=1,qnstep
         qnbeta(i) = vlsc2(yk(1,i),qnsol,nb)/qnrho(i)
         qnfact(i) = (qnalpha(i)-qnbeta(i))
         call add2s2(qnsol,sk(1,i),qnfact(i),nb)
      enddo
      call chsign(qnsol,nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine findminalpha(minalpha,s,uu,amax,amin,nb)

      real minalpha
      real s(nb),uu(nb),amax(nb),amin(nb)
      real alphai(nb)
      integer nb

      do ii=1,nb
         if (s(ii).le.1e-10) then
            alphai(ii) = abs((uu(ii)-amin(ii))/s(ii))
         else
            alphai(ii) = ((amax(ii)-uu(ii))/s(ii))
         endif
      enddo

      minalpha = vlmin(alphai,nb)
      return
      end
c-----------------------------------------------------------------------
      subroutine invHessian_update(B,s,y,nb)

      real B(nb,nb)
      real s(nb),y(nb)
      real rk1up(nb,nb)
      real rk2up(nb,nb)
      real w1(nb,nb),w2(nb,nb)
      real w3(nb,nb),w4(nb,nb),w5(nb,nb)
      real ss(nb,nb),ys,ys2
      real hy(nb),yhy,nomina

      ! y_k^T * s_k
      ys = vlsc2(y,s,nb)
      ysinv = 1./ys

      ys2 = ys*ys
      ys2inv = 1./ys2

      ! y_k * s_k^T
      call mxm(y,nb,s,1,w1,nb)
      ! s_k * y_k^T
      call mxm(s,nb,y,1,w2,nb)

      call mxm(B,nb,w1,nb,w3,nb)
      call cmult(w3(1,1),ysinv,nb*nb)

      call mxm(w2,nb,B,nb,w4,nb)
      call cmult(w4(1,1),ysinv,nb*nb)

      call mxm(B,nb,y,nb,hy,1)
      yhy = vlsc2(y,hy,nb)

      nomina=ys+yhy

      ! second rank-one update
      ! s_k * s_k^T               
      call mxm(s,nb,s,1,ss,nb)

      call cmult(ss(1,1),ys2inv,nb*nb)
      call cmult(ss(1,1),nomina,nb*nb)

      call sub2(B(1,1),w3(1,1),nb*nb)
      call sub2(B(1,1),w4(1,1),nb*nb)
      call add2(B(1,1),ss(1,1),nb*nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine BFGS(rhs,vv,helm,invhelm,amax,amin,adis,
     $   bpar,bstep,ifdiag)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real helm(nb,nb),invhelm(nb,nb),B_qn(nb,nb)
      real qgo(nb),qngradf(nb)
      real qns(nb),qny(nb)
      real uu(nb),vv(nb),rhs(nb)
      real amax(nb),amin(nb),adis(nb)
      real tmp(nb,nb)
      real ngf,ysk
      real bpar,par
      real norm_s,norm_step,norm_uo

      real invB_qn(nb,nb)
      real copt_work(nb*nb)

      ! parameter for barrier function
      integer par_step,jmax,bstep,chekbc
      integer lncount
      logical ifdiag

      call copy(uu,vv,nb)

      jmax = 0

      par = bpar 
      par_step = bstep 

      ! BFGS method with barrier function starts
      tcopt_time=dnekclock()
      do k=1,par_step

         chekbc = 0

         ! use helm from BDF3/EXT3 as intial approximation
         call copy(B_qn(1,1),helm(1,1),nb*nb)
         call comp_qngradf(uu,rhs,helm,qngradf,amax,amin,par,
     $                     ifdiag,barr_func,nb)

         if (isolve.eq.5) then
            call copy(invB_qn(1,1),invhelm(1,1),nb*nb)
            call dgetri(nb,invB_qn,nb,ipiv,copt_work,nb*nb,info)
         endif
         norm_uo = vlamax(uu,nb)

         ! compute quasi-Newton step
         tquasi_time=dnekclock()
         do j=1,100

            if (isolve.eq.5) then
               ONE = 1.
               ZERO = 0.
               call dgemv('N',nb,nb,ONE,invB_qn,nb,qngradf,1,ZERO,qns,1)
               call chsign(qns,nb)
            else
               call copy(tmp(1,1),B_qn(1,1),nb*nb)
               call dgetrf(nb,nb,tmp,nb,ipiv,info)
               call copy(qns,qngradf,nb)
               call chsign(qns,nb)
               call dgetrs('N',nb,1,tmp,nb,ipiv,qns,nb,info)
            endif

            tlnsrch_time=dnekclock()
            call backtrackr(uu,qns,rhs,qngradf,helm,invhelm,1e-4,0.5,
     $                  amax,amin,par,chekbc,lncount,tol_box,ifdiag)
            lnsrch_time=lnsrch_time+dnekclock()-tlnsrch_time

            norm_s = vlamax(qns,nb)
            norm_step = norm_s/norm_uo

            ! store old qn-gradf
            call copy(qgo,qngradf,nb) 
            ! update qn-gradf
            call comp_qngradf(uu,rhs,helm,qngradf,amax,amin,par,
     $      ifdiag,barr_func,nb)
            call sub3(qny,qngradf,qgo,nb) 

            ! compute curvature condition
            ysk = vlsc2(qny,qns,nb)

            if (isolve.eq.5) then
               call invHessian_update(invB_qn,qns,qny,nb)
            else
               call Hessian_update(B_qn,qns,qny,nb)
            endif

            ! compute norm of gradf
            ngf = vlamax(qngradf,nb)

            ! reset chekbc 
            chekbc = 0
            
            jmax = max(j,jmax)

            if (ngf .lt. 1e-6 .OR. ysk .lt. 1e-6 .OR. norm_step .lt.
     $      1e-6  ) then 
               if (ysk .lt. 1e-10) then 
                  if (nio.eq.0) then
                  ! curvature condition not satisfied
                  write(6,*)'WARNING:decrease barrseq or increase barr0'
                  endif
               endif
               exit
            endif

         enddo
         quasi_time=quasi_time+dnekclock()-tquasi_time

         if (mod(ad_step,ad_iostep).eq.0) then
            if (nio.eq.0) write (6,*) 'lnconst_ana'
            call cpod_ana(uu,par,j,lncount,ngf,norm_step
     $                   ,ysk)
         endif
         par = par*0.1

      enddo
      copt_time=copt_time+dnekclock()-tcopt_time

      call copy(rhs,uu,nb)

      return
      end
