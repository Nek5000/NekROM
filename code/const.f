      subroutine BFGS(rhs,helm,invhelm,amax,amin,adis,bpar,bstep)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real B_qn(nb,nb),helm(nb,nb),invhelm(nb,nb)
      real tmp(nb,nb)
      real qgo(nb),qngradf(nb),ngf
      real fo,qnf,qndf
      real ww(nb),pert
      real uu(nb),rhs(nb)
      real amax(nb),amin(nb), adis(nb)
      real bpar,par
      real alphak

      ! parameter for barrier function
      integer par_step,jmax,bflag,bstep
      integer chekbc ! flag for checking boundary
      integer uHcount
      real bctol

      call copy(uu,u(1,1),nb)

      bctol = 1e-8
      jmax = 0

      bflag = 1 
      par = bpar 
      par_step = bstep 
   
      ! BFGS method with barrier function starts
      do k=1,par_step

         chekbc = 0
         uHcount = 0

c        use helm from BDF3/EXT3 as intial approximation
         call copy(B_qn(1,1),helm(1,1),nb*nb)

         call comp_qnf(uu,rhs,helm,invhelm,qnf,amax,amin,par,bflag)
         call comp_qngradf(uu,rhs,helm,qngradf,amax,amin,par,bflag)

c        compute quasi-Newton step
         do j=1,100

            call copy(tmp(1,1),B_qn(1,1),nb*nb)
            call dgetrf(nb,nb,tmp,lub,ipiv,info)
            call copy(qns,qngradf,nb)
            call chsign(qns,nb)
            call dgetrs('N',nb,1,tmp,lub,ipiv,qns,nb,info)

            if (isolve.eq.1) then
               call backtrackr(uu,qns,rhs,helm,invhelm,1e-2,0.5,alphak,amax,
     $                     amin,bctol,bflag,par)
            elseif (isolve.eq.2) then      
               call add2(uu,qns,nb)
            endif

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

            call copy(qgo,qngradf,nb) ! store old qn-gradf
            call comp_qngradf(uu,rhs,helm,qngradf,amax,amin,par,bflag) ! update qn-gradf
            call sub3(qny,qngradf,qgo,nb) 

            ! update approximate Hessian by two rank-one update if chekbc = 0
            if (chekbc .ne. 1) then
               uHcount = uHcount + 1
               call Hessian_update(B_qn,qns,qny,nb)
            endif

            ! compute H^{-1} norm of gradf
            call copy(ww,qngradf,nb)
            call dgetrs('N',nb,1,invhelm,lub,ipiv,ww,nb,info)
            ngf = glsc2(ww,qngradf,nb)
            ngf = sqrt(ngf)

            fo = qnf      ! store old qn-f
            call comp_qnf(uu,rhs,helm,invhelm,qnf,amax,amin,par,bflag) ! update qn-f
            qndf = abs(qnf-fo)/abs(fo) 

            if (mod(ad_step,ad_iostep).eq.0) then
               if (nio.eq.0) write (6,*) 'lnconst_ana'
               call cpod_ana(uu,par,j,uHcount,ngf,qndf)
            endif

            ! reset chekbc 
            chekbc = 0
            
            jmax = max(j,jmax)

            if (ngf .lt. 1e-4 .OR. qndf .lt. 1e-6  ) then 
               exit
            endif

      ! update solution
         enddo
         par = par*0.1
      enddo
      call copy(rhs,uu,nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_qngradf(uu,rhs,helm,s,amax,amin,bpar,barr_func)
      
      include 'SIZE'
      include 'MOR'

      real tmp1(nb),tmp2(nb),tmp3(nb),tmp4(nb)
      real helm(nb,nb)
      real mpar
      real uu(nb), rhs(nb), s(nb)
      real amax(nb), amin(nb) 
      real bpar
      integer barr_func


!      if (nio.eq.0) write (6,*) 'inside com_qngradf'

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

!      if (nio.eq.0) write (6,*) 'exiting com_qngradf'

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

!      if (nio.eq.0) write (6,*) 'inside com_qnf'

c     evaluate quasi-newton f
      ONE = 1.
      ZERO= 0.

      call dgemv('N',nb,nb,ONE,helm,nb,uu,1,ZERO,tmp6,1) ! H*coef
      term1 = 0.5 * glsc2(tmp6,uu,nb) ! coef'*H*coef

      term2 = glsc2(uu,rhs,nb) ! coef'*rhs

      ! 0.5*rhs'*inv(H)*rhs
      call copy(tmp5,rhs,nb)
      call dgetrs('N',nb,1,invhelm,lub,ipiv,tmp5,nb,info)
c     call dgemv('N',nb,nb,ONE,invhelm,nb,rhs,1,ZERO,tmp5,1)
c     call solve(tmp5,invhelm,1,nb,nb,irv,icv)
      term3 = 0.5 * glsc2(rhs,tmp5,nb)

c     call dgemv('N',nb,nb,ONE,helm,nb,tmp5,1,ZERO,work,1)

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

!      if (nio.eq.0) write (,*) 'exitting com_qnf'

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
      subroutine invHessian_update(B,s,y,nb)

      real B(nb,nb)
      real s(nb), y(nb)
      real w1(nb,nb), w2(nb,nb)
      real w3(nb,nb), w4(nb,nb), w5(nb,nb)
      real ss(nb,nb), ys, sBs
      real sds

      ! y_k * s_k^T
      call mxm(y,nb,s,1,w1,nb)
      ! s_k * y_k^T
      call mxm(s,nb,y,1,w2,nb)
      ! y_k^T * s_k
      ys = glsc2(y,s,nb)

      call mxm(B,nb,w1,nb,w3,nb)
      call cmult(w3,-1.0/ys,nb*nb)
      do j=1,nb
      do i=1,nb
      write(6,*)i,j,w3(i,j)
      enddo
      enddo

      call mxm(w2,nb,B,nb,w4,nb)
      call cmult(w4,-1.0/ys,nb*nb)

      call mxm(w4,nb,w1,nb,w5,nb)
      call cmult(w5,1.0/(ys**2),nb*nb)

      ! second rank-one update
      ! s_k * s_k^T               
      call mxm(s,nb,s,1,ss,nb)
      ! s_k^T * s_k               
      sds = glsc2(s,s,nb)
      call cmult(ss,1.0/sds,nb*nb)

      do ii=1,nb
         call add4(B(1,ii),w3(1,ii),w4(1,ii)
     $            ,w5(1,ii),nb)
      enddo

      call add2(B(1,1),ss(1,1),nb*nb)
      call exitt0

      return
      end
c-----------------------------------------------------------------------
      subroutine backtrackr(uu,s,rhs,helm,invhelm,sigmab,facb,alphak,amax,
     $            amin,bctol,bflag,bpar)

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
      real bctol,bpar
      integer chekbc,counter
      integer bflag

      alphak = 1.0
      chekbc = 1
      counter = 0

      call comp_qnf(uu,rhs,helm,invhelm,fk,amax,amin,bpar,bflag) ! get old f
      call comp_qngradf(uu,rhs,helm,Jfk,amax,amin,bpar,bflag)

      call copy(uuo,uu,nb)
      call add2s1(uu,s,alphak,nb)

      call comp_qnf(uu,rhs,helm,invhelm,fk1,amax,amin,bpar,bflag) ! get new f
      Jfks = vlsc2(Jfk,s,nb)

c     do while ((fk1 > fk + sigmab * alphak * Jfks) .OR. (chekbc.eq.1))
      do while ((chekbc.neq.0) .and. (fk1 .gt. fk + sigmab * alphak * Jfks))
         counter = counter + 1
         alphak = alphak * facb
         call add3s2(uu,uuo,s,1.0,alphak,nb)

         chekbc = 0
         do ii=1,nb
            if ((uu(ii)-amax(ii)).ge.bctol) then
               chekbc = 1
            elseif ((amin(ii)-uu(ii)).ge.bctol) then
               chekbc = 1
            endif
         enddo

         call comp_qnf(uu,rhs,helm,invhelm,fk1,amax,amin,bpar,bflag)
         
         if (alphak < 1e-4) then
            if (mod(ad_step,ad_iostep).eq.0) then
               if (nio.eq.0) write(6,*)
     $         '# lnsrch:',counter,'alpha',alphak
            endif
            exit
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cpod_ana(uu,par,qstep,uhcount,ngf,qndf)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real uu(nb)
      real par
      real ngf, qndf
      integer qstep 
      integer uhcount

      if (nio.eq.0) then
         write (6,*)'ad_step:',ad_step,ad_iostep,par,qstep,uhcount,
     $            ngf,qndf
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
c     subroutine invH_multiply(invh0,s,y,d)

c     include 'MOR'            

c     real invh0(nb,nb)
c     real s(nb),y(nb),d(nb)

c     return
c     end
