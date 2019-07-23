      subroutine BFGS(rhs,helm,invhelm,amax,amin,adis,bpar,bstep)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real B_qn(nb,nb),helm(nb,nb),invhelm(nb,nb)
      real tmp(nb,nb)
      real qgo(nb),qngradf(nb),ngf
      real qnf
      real ww(nb),pert
      real uu(nb),rhs(nb)
      real amax(nb),amin(nb),adis(nb)
      real bpar,par
      real norm_s,norm_step
      real norm_uo
      real ysk

      ! parameter for barrier function
      real bctol
      integer par_step,bflag,bstep

      ! parameter for quasi-newton
      integer uHcount,jamx,chekbc

      ! parameter for lnsrch
      integer lncount
      lncount = 0

      call copy(uu,u(1,1),nb)

      jmax = 0

      bflag = 1 
      par = bpar 
      par_step = bstep 
      bctol = 1e-8

      ! BFGS method with barrier function starts
      do k=1,par_step

         chekbc = 0
         uHcount = 0

         ! use helm from BDF3/EXT3 as intial approximation
         call copy(B_qn(1,1),helm(1,1),nb*nb)
         call comp_qnf(uu,rhs,helm,invhelm,qnf,amax,amin,par,bflag)
         call comp_qngradf(uu,rhs,helm,qngradf,amax,amin,par,bflag)

         norm_uo = glamax(uu,nb)

c        compute quasi-Newton step
         do j=1,100

            call copy(tmp(1,1),B_qn(1,1),nb*nb)
            call dgetrf(nb,nb,tmp,lub,ipiv,info)
            call copy(qns,qngradf,nb)
            call chsign(qns,nb)
            call dgetrs('N',nb,1,tmp,lub,ipiv,qns,nb,info)

            if (isolve.eq.1) then
               call backtrackr(uu,qns,rhs,helm,invhelm,1e-4,0.5,
     $                     amax,amin,bctol,bflag,par,chekbc,lncount)
            elseif (isolve.eq.2) then      
               call add2(uu,qns,nb)
               ! check the boundary 
               do ii=1,nb
                  if ((uu(ii)-amax(ii)).ge.bctol) then
                     chekbc = 1
                     uu(ii) = amax(ii) - 0.1*adis(ii)
                  elseif ((amin(ii)-uu(ii)).ge.bctol) then
                     chekbc = 1
                     uu(ii) = amin(ii) + 0.1*adis(ii)
                  endif
               enddo
            endif

            norm_s = glamax(qns,nb)
            norm_step = norm_s/norm_uo

            call copy(qgo,qngradf,nb) ! store old qn-gradf
            call comp_qngradf(uu,rhs,helm,qngradf,amax,amin,
     $                  par,bflag) ! update qn-gradf
            call sub3(qny,qngradf,qgo,nb) 

            ysk = glsc2(qny,qns,nb)

            ! update approximate Hessian by two rank-one update if chekbc = 0
            if (chekbc .ne. 1) then
               uHcount = uHcount + 1
               call Hessian_update(B_qn,qns,qny,nb)
            endif

            ! compute H^{-1} norm of gradf
            call copy(ww,qngradf,nb)
            ngf = glamax(qngradf,nb)

            call comp_qnf(uu,rhs,helm,invhelm,qnf,amax,amin,par,bflag) ! update qn-f

            ! reset chekbc 
            chekbc = 0
            
            jmax = max(j,jmax)

            if (ngf .lt. 1e-6 .OR. norm_step .lt. 1e-10) then 
               exit
            endif

         enddo

         if (mod(ad_step,ad_iostep).eq.0) then
            if (nio.eq.0) write (6,*) 'lnconst_ana'
            call cpod_ana(uu,par,j,uHcount,lncount,ngf,norm_step
     $      ,norm_s,ysk)
         endif
         par = par*0.1
      enddo
      call copy(rhs,uu,nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine BFGS_new(rhs,vv,helm,invhelm,amax,amin,adis,bpar,bstep)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real helm(nb,nb),invhelm(nb,nb)
      real tmp(nb,nb)
      real qgo(nb),qngradf(nb),ngf
      real qnf
      real ww(nb),pert
      real uu(nb),vv(nb),rhs(nb)
      real amax(nb),amin(nb),adis(nb)
      real bpar,par
      real sk(nb,50),yk(nb,50)

      ! parameter for barrier function
      integer par_step,jmax,bflag,bstep
      integer chekbc ! flag for checking boundary
      integer uHcount,lncount
      real bctol
      real norm_s,norm_step,norm_uo
      real ysk

      call copy(uu,vv,nb)

      bctol = 1e-12
      jmax = 0

      bflag = 1
      par = bpar 
      par_step = bstep 

   
      tcopt_time=dnekclock()
      ! BFGS method with barrier function starts
      do k=1,par_step

         chekbc = 0
         uHcount = 0

         ! use helm from BDF3/EXT3 as intial approximation
         call comp_qnf(uu,rhs,helm,invhelm,qnf,amax,amin,par,bflag)
         call comp_qngradf(uu,rhs,helm,qngradf,amax,amin,par,bflag)

         norm_uo = glamax(uu,nb)

         tquasi_time=dnekclock()
c        compute quasi-Newton step
         do j=1,50

            if (isolve.eq.1.OR.isolve.eq.2.OR.isolve.eq.3) then
               if (j.eq.1) then
                  call copy(qns,qngradf,nb)
                  call chsign(qns,nb)
                  call dgetrs('N',nb,1,invhelm,lub,ipiv,qns,nb,info)
               endif

               tlnsrch_time=dnekclock()
               call backtrackr(uu,qns,rhs,helm,invhelm,1e-4,0.5,
     $                     amax,amin,bctol,bflag,par,chekbc,lncount)
               lnsrch_time=lnsrch_time+dnekclock()-tlnsrch_time
               ! store qns
               call copy(sk(1,j),qns,nb)

               call copy(qgo,qngradf,nb) ! store old qn-gradf
               call comp_qngradf(uu,rhs,helm,qngradf,amax,amin,
     $                     par,bflag) ! update qn-gradf
               call sub3(qny,qngradf,qgo,nb) 

               ysk = glsc2(qny,qns,nb)

               norm_s = glamax(qns,nb)
               norm_step = norm_s/norm_uo

               ! compute H^{-1} norm of gradf
               call copy(ww,qngradf,nb)
               ngf = glamax(qngradf,nb)

               call comp_qnf(uu,rhs,helm,invhelm,qnf,amax,amin,
     $               par,bflag) ! update qn-f

               jmax = max(j,jmax)

               ! store qny 
               call copy(yk(1,j),qny,nb)
               call invH_multiply(qns,invhelm,sk,yk,qngradf,j)
            endif

            if (ngf .lt. 1e-6 .OR. ysk .lt. 1e-10 .OR. norm_step .lt.
     $      1e-10  ) then 
               exit
            endif

      ! update solution
         enddo

         quasi_time=quasi_time+dnekclock()-tquasi_time

         if (mod(ad_step,ad_iostep).eq.0) then
            if (nio.eq.0) write (6,*) 'lnconst_ana'
            call cpod_ana(uu,par,j,uHcount,lncount,ngf,norm_step
     $      ,norm_s,ysk)
         endif
         par = par*0.1

      enddo

      copt_time=copt_time+dnekclock()-tcopt_time

      call copy(rhs,uu,nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_qngradf(uu,rhs,helm,s,amax,amin,bpar,barr_func)
      
      include 'SIZE'
      include 'MOR'

      real helm(nb,nb)
      real uu(nb),rhs(nb),s(nb)
      real amax(nb),amin(nb) 
      real tmp1(nb),tmp2(nb),tmp3(nb),tmp4(nb)
      real bpar,mpar,pert
      integer barr_func

      if (barr_func .eq. 1) then ! use logarithmic as barrier function

         call sub3(tmp1,uu,amax,nb)  
         call sub3(tmp2,uu,amin,nb)  
   
         ! add perturbation 
         call cadd(tmp1,-1e-2,nb)
         call cadd(tmp2,1e-2,nb)
   
         call invcol1(tmp1,nb)
         call invcol1(tmp2,nb)
         call add3(tmp3,tmp1,tmp2,nb)

         mpar = -1.0*bpar
         call add3s12(s,rhs,tmp3,-1.0,mpar,nb)
   
         ONE = 1.
         ZERO= 0.
         call dgemv('N',nb,nb,ONE,helm,nb,uu,1,ZERO,tmp4,1)
         call add2(s,tmp4,nb)
c        call sub2(tmp4,s,nb)

      else ! use inverse function as barrier function

         call sub3(tmp1,uu,amax,nb)  
         call sub3(tmp2,amin,uu,nb)  
         call vsq(tmp1,nb)
         call vsq(tmp2,nb)
         call invcol1(tmp1,nb)
         call invcol1(tmp2,nb)
         call sub3(tmp3,tmp2,tmp1,nb)

         mpar = -1.0*bpar
         call add3s12(s,rhs,tmp3,-1.0,mpar,nb)
   
         ONE = 1.
         ZERO= 0.
         call dgemv('N',nb,nb,ONE,helm,nb,uu,1,ZERO,tmp4,1)
         call add2(s,tmp4,nb)

      endif

      return 
      end
c-----------------------------------------------------------------------
      subroutine comp_qnf(uu,rhs,helm,invhelm,qnf,amax,amin,bpar,
     $                     barr_func)
      
      include 'SIZE'
      include 'MOR'

      real tmp1(nb),tmp2(nb),tmp3(nb)
      real tmp4(nb),tmp5(nb),tmp6(nb)
      real term1,term2,term3,term4
      real bar1,bar2 ! bar1 and bar2 are the barrier function for two constrains
      real uu(nb), rhs(nb)
      real amax(nb), amin(nb)
      real helm(nb,nb), invhelm(nb,nb)
      real qnf
      real bpar
      integer barr_func 

      barr_func = 1

      ! evaluate quasi-newton f

      ONE = 1.
      ZERO= 0.

      call dgemv('N',nb,nb,ONE,helm,nb,uu,1,ZERO,tmp6,1) ! H*coef
      term1 = 0.5 * glsc2(tmp6,uu,nb) ! coef'*H*coef

      term2 = glsc2(uu,rhs,nb) ! coef'*rhs

      ! 0.5*rhs'*inv(H)*rhs
      call copy(tmp5,rhs,nb)
      call dgetrs('N',nb,1,invhelm,lub,ipiv,tmp5,nb,info)
      term3 = 0.5 * glsc2(rhs,tmp5,nb)

      if (barr_func .eq. 1) then ! use logarithmetic as barrier function
         ! barrier term
         call sub3(tmp1,amax,uu,nb)  
         call sub3(tmp2,uu,amin,nb)  
   
         ! add perturbation
         call cadd(tmp1,1e-2,nb)
         call cadd(tmp2,1e-2,nb)
   
         do i=1,nb
            tmp3(i) = log(tmp1(i))
            tmp4(i) = log(tmp2(i))
         enddo
      else 
         call sub3(tmp1,uu,amax,nb)  
         call sub3(tmp2,amin,uu,nb)  
   
         do i=1,nb
            tmp3(i) = 1./tmp1(i)
            tmp4(i) = 1./tmp2(i)
         enddo
      endif

      bar1 = vlsum(tmp3,nb)
      bar2 = vlsum(tmp4,nb)
      term4 = bpar*(bar1+bar2)

      qnf = term1 - term2 + term3 - term4

      return
      end
c-----------------------------------------------------------------------
      subroutine Hessian_update(B,s,y,nb)

      real B(nb,nb)
      real s(nb),y(nb)
      real w1(nb,nb),w2(nb,nb),w3(nb,nb)
      real w4(nb),w5(nb),w6(nb,nb),w7(nb,nb)
      real yy(nb,nb),ys,sBs
      
      ! s_k * s_k^T               
      call mxm(s,nb,s,1,w1,nb)

      ! s_k * s_k^T * B_k
      call mxm(w1,nb,B,nb,w2,nb)

      ! B_k * s_k * s_k^T * B_k 
      call mxm(B,nb,w2,nb,w3,nb)

      ! s_k^T * B_k * s_k 
      call mxm(B,nb,s,nb,w4,1)

      sBs = glsc2(s,w4,nb)

      ! second rank-one update
      ! y_k * y_k^T               
      call mxm(y,nb,y,1,yy,nb)
      ! y_k^T * s_k
      ys = glsc2(y,s,nb)

      call cmult(w3,-1.0/sBs,nb*nb)
      call cmult(yy,1.0/ys,nb*nb)

      call add4(B(1,1),B(1,1),w3(1,1),yy(1,1),nb*nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine backtrackr(uu,s,rhs,helm,invhelm,sigmab,facb,
     $            amax,amin,bctol,bflag,bpar,chekbc,counter)

      include 'SIZE'
      include 'MOR'
      include 'TOTAL'

      real rhs(nb),s(nb)
      real uuo(nb),uu(nb)
      real helm(nb,nb),invhelm(nb,nb)
      real Jfk(nb)
      real amax(nb),amin(nb)
      real fk,fk1
      real Jfks
      real sigmab,facb,alphak
      real bctol,bpar,minalpha

      integer chekbc,counter
      integer bflag
      integer countbc
      logical cond1

      alphak = 1.0
      chekbc = 0
      counter = 0
      countbc = 0

      call comp_qnf(uu,rhs,helm,invhelm,fk,amax,amin,bpar,bflag) ! get old f
      call comp_qngradf(uu,rhs,helm,Jfk,amax,amin,bpar,bflag)

      call findminalpha(minalpha,s,uu,amax,amin)

      call copy(uuo,uu,nb)
      call add2s1(uu,s,alphak,nb)

      do ii=1,nb
         if ((uu(ii)-amax(ii)).ge.bctol) then
            chekbc = 1
            countbc = countbc + 1
         elseif ((amin(ii)-uu(ii)).ge.bctol) then
            chekbc = 1
            countbc = countbc + 1
         endif
      enddo

      call comp_qnf(uu,rhs,helm,invhelm,fk1,amax,amin,bpar,bflag) ! get new f
      Jfks = vlsc2(Jfk,s,nb)   

      cond1 = fk1 .gt. (fk+sigmab*alphak*Jfks)

      do while (cond1.OR.(chekbc.eq.1))
         counter = counter + 1
         alphak = alphak * facb
         call add3s2(uu,uuo,s,1.0,alphak,nb)

         chekbc = 0
         countbc = 0
         do ii=1,nb
            if ((uu(ii)-amax(ii)).ge.bctol) then
               chekbc = 1
               countbc = countbc + 1
            elseif ((amin(ii)-uu(ii)).ge.bctol) then
               chekbc = 1
               countbc = countbc + 1
            endif
         enddo

         call comp_qnf(uu,rhs,helm,invhelm,fk1,amax,amin,bpar,bflag)

         cond1 = fk1 .gt. (fk+sigmab*alphak*Jfks)
         
c        if (alphak < minalpha .AND. (fk1.gt.(fk+sigmab*alphak*Jfks))) 
c    $   then
         if ((alphak < minalpha.AND..not.cond1).OR.alphak < 1e-8) then
            exit
         endif
      enddo

      call cmult(s,alphak,nb)

      if (mod(ad_step,ad_iostep).eq.0) then
         if (nio.eq.0) write(6,*)
     $         ad_step,'# lnsrch:',counter,'alpha',alphak,
     $         minalpha,chekbc,cond1
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cpod_ana(uu,par,qstep,uhcount,lncount,ngf,qndf,norm_s,
     $   ysk)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real uu(nb)
      real par
      real ngf,qndf,norm_s,yks
      integer qstep 
      integer uhcount,lncount

      if (nio.eq.0) then
         write (6,*)'ad_step:',ad_step,ad_iostep,par,qstep,uhcount,
     $            lncount,ngf,qndf,norm_s,ysk
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
      
      return
      end
c-----------------------------------------------------------------------
      subroutine invH_multiply(qnsol,invh0,sk,yk,qnd,qnstep)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'            

      real invh0(nb,nb)
      real sk(nb,50),yk(nb,50),qnd(nb),qnsol(nb)
      real qnrho(50),qnalpha(50),qnbeta(50)
      real qnfact(50)
      real work(nb)
      integer qnstep

      call copy(qnsol,qnd,nb)
      call chsign(qnsol,nb)

      ! compute right product
      do i=qnstep,1,-1
         qnrho(i) = glsc2(yk(1,i),sk(1,i),nb)
         qnalpha(i) = glsc2(sk(1,i),qnsol,nb)/qnrho(i)
         call add2s2(qnsol,yk(1,i),-qnalpha(i),nb)
      enddo

      ! compute center
      ONE = 1.
      ZERO= 0.
      call dgetrs('N',nb,1,invh0,lub,ipiv,qnsol,nb,info)

      ! compute left product
      do i=1,qnstep
         qnbeta(i) = glsc2(yk(1,i),qnsol,nb)/qnrho(i)
         qnfact(i) = (qnalpha(i)-qnbeta(i))
         call add2s2(qnsol,sk(1,i),qnfact(i),nb)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine findminalpha(minalpha,s,uu,amax,amin)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real minalpha
      real s(nb),uu(nb),amax(nb),amin(nb)
      real alphai(nb)

      do ii=1,nb
         if (s(ii).le.0) then
            alphai(ii) = abs((uu(ii)-amin(ii))/s(ii))
         else
            alphai(ii) = ((amax(ii)-uu(ii))/s(ii))
         endif
      enddo

      minalpha = vlmin(alphai,nb)
      return
      end
