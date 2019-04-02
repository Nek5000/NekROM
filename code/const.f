      subroutine BFGS_freeze(rhs)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real B_qn(nb,nb)
      real qgo(nb), qngradf(nb), ngf
      real fo,qnf,qndf
      real yy(nb,nb),ys,sBs
      real ww(nb), pert
      real uu(nb), rhs(nb)

      integer par_step, jmax
      integer chekbc ! flag for checking boundary

      real tmp(nb,nb),tmp1(nb,nb),tmp2(nb,nb),tmp3(nb,nb)
      real tmp4(nb),tmp5(nb),tmp6(nb,nb),tmp7(nb,nb)

c      if (nio.eq.0) write (6,*) 'inside BFGS_freeze'
      jmax = 0

      call copy(uu,u(1,1),nb)
   
      ! parameter for barrier function
      par = 1e-3 
      par_step = 4
   
      ! invhelm for computing qnf
      if (ad_step.le.3) then 
         call copy(invhelm(1,1),helm(1,1),nb*nb)
         call lu(invhelm,nb,nb,irv,icv)
      endif

      ! BFGS method with barrier function starts
      do k=1,par_step

         chekbc = 0

c        use helm from BDF3/EXT3 as intial approximation
         call copy(B_qn(1,1),helm(1,1),nb*nb)
!         call copy(B_qn(1,1),invhelm(1,1),nb*nb)

         call comp_qnf(uu,rhs,qnf)
         call comp_qngradf(uu,rhs,qngradf)

c        compute quasi-Newton step
         do j=1,500

            call copy(tmp3(1,1),B_qn(1,1),nb*nb)
            call lu(tmp3,nb,nb,irv,icv)
            call copy(qns,qngradf,nb)
            call chsign(qns,nb)
            call solve(qns,tmp3,1,nb,nb,irv,icv)

c            if (j .eq. 1) then
c               call copy(qns,qngradf,nb)
c               call solve(qns,B_qn,1,nb,nb,ir,ic)
c            else     
c               call mxm(B_qn,nb,qngradf,nb,qns,1)
c            endif

c            call chsign(qns,nb)
            call add2(uu,qns,nb)

            ! check the boundary 
            do ii=1,nb
               if ((uu(ii)-sample_max(ii)).ge.1e-8) then
                  chekbc = 1
                  uu(ii) = sample_max(ii) - 0.1*sam_dis(ii)
               elseif ((sample_min(ii)-uu(ii)).ge.1e-8) then
                  chekbc = 1
                  uu(ii) = sample_min(ii) + 0.1*sam_dis(ii)
               endif
            enddo

            call copy(qgo,qngradf,nb) ! store old qn-gradf
            call comp_qngradf(uu,rhs,qngradf)        ! update qn-gradf
            call sub3(qny,qngradf,qgo,nb) 

            ! update approximate Hessian by two rank-one update if chekbc = 0
            ! first rank-one update
            if (chekbc .NE. 1) then  

               call Hessian_update(B_qn,qns,qny,nb)
!               call invHessian_update(B_qn,qns,qny,nb)

            endif   

            call copy(ww,qngradf,nb)
            call solve(ww,invhelm,1,nb,nb,irv,icv)

            ngf = glsc2(ww,qngradf,nb)
            ngf = sqrt(ngf)

            fo = qnf      ! store old qn-f
            call comp_qnf(uu,rhs,qnf) ! update qn-f
            qndf = abs(qnf-fo)/abs(fo) 
c            write(6,*)'f and old f',j,qnf,fo,qndf,ngf

            ! reset chekbc 
            chekbc = 0
            
            jmax = max(j,jmax)

            if (ngf .lt. 1e-4 .OR. qndf .lt. 1e-6  ) then
               exit
            endif

c     update solution
         enddo
         par = par*0.1

      enddo

      if (mod(ad_step,ad_iostep).eq.0) then
         if (nio.eq.0) then
            write(6,*)'ad_step, par, jmax, ngf, qndf:',
     $               ad_step,par,jmax,ngf,qndf 
         endif
      endif

      call copy(rhs,uu,nb)

c      if (nio.eq.0) write (6,*) 'exitting BFGS_freeze'

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_qngradf(uu,rhs,s)
      
      include 'SIZE'
      include 'MOR'

      real tmp1(nb),tmp2(nb),tmp3(nb),tmp4(nb)
      real mpar
      real uu(nb), rhs(nb), s(nb)
      integer barr_func

      barr_func = 1

!      if (nio.eq.0) write (6,*) 'inside com_qngradf'

      if (barr_func .eq. 1) then ! use logarithmic as barrier function

         call sub3(tmp1,uu,sample_max,nb)  
         call sub3(tmp2,uu,sample_min,nb)  
   
         ! add perturbation 
         call cadd(tmp1,-1e-2,nb)
         call cadd(tmp2,1e-2,nb)
   
         call invcol1(tmp1,nb)
         call invcol1(tmp2,nb)
         call add3(tmp3,tmp1,tmp2,nb)

         mpar = -1.0*par
         call add3s12(s,rhs,tmp3,-1.0,mpar,nb)
   
         ONE = 1.
         ZERO= 0.
         call dgemv('N',nb,nb,ONE,helm,nb,uu,1,ZERO,tmp4,1)
         call add2(s,tmp4,nb)


      else ! use inverse function as barrier function

         call sub3(tmp1,uu,sample_max,nb)  
         call sub3(tmp2,sample_min,uu,nb)  
         call vsq(tmp1,nb)
         call vsq(tmp2,nb)
         call invcol1(tmp1,nb)
         call invcol1(tmp2,nb)
         call sub3(tmp3,tmp2,tmp1,nb)

         mpar = -1.0*par
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
      subroutine comp_qnf(uu,rhs,qnf)
      
      include 'SIZE'
      include 'MOR'

      real tmp1(nb),tmp2(nb),tmp3(nb)
      real tmp4(nb),tmp5(nb),tmp6(nb)
      real term1,term2,term3,term4
      real bar1,bar2 ! bar1 and bar2 are the barrier function for two constrains
      real uu(nb), rhs(nb)
      real qnf
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
      call solve(tmp5,invhelm,1,nb,nb,irv,icv)
      term3 = 0.5 * glsc2(rhs,tmp5,nb)

      if (barr_func .eq. 1) then ! use logarithmetic as barrier function

         ! barrier term
         call sub3(tmp1,sample_max,uu,nb)  
         call sub3(tmp2,uu,sample_min,nb)  
   
         ! add perturbation
         call cadd(tmp1,1e-2,nb)
         call cadd(tmp2,1e-2,nb)
   
         do i=1,nb
            tmp3(i) = log(tmp1(i))
            tmp4(i) = log(tmp2(i))
         enddo

      else 

         call sub3(tmp1,uu,sample_max,nb)  
         call sub3(tmp2,sample_min,uu,nb)  
   
         do i=1,nb
            tmp3(i) = 1./tmp1(i)
            tmp4(i) = 1./tmp2(i)
         enddo

      endif

      bar1 = vlsum(tmp3,nb)
      bar2 = vlsum(tmp4,nb)
      term4 = par*(bar1+bar2)

      qnf = term1 - term2 + term3 - term4

!      if (nio.eq.0) write (,*) 'exitting com_qnf'

      return
      end
c-----------------------------------------------------------------------
      subroutine Hessian_update(B,s,y,nb)

      real B(nb,nb)
      real s(nb), y(nb)
      real w1(nb,nb),w2(nb,nb),w3(nb,nb)
      real w4(nb),w5(nb),w6(nb,nb),w7(nb,nb)
      real yy(nb,nb),ys,sBs
      
c      if (nio.eq.0) write(6,*) 'inside Hessian_update'

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
!      do ii=1,nb
!         call cmult(w3(1,ii),-1.0/sBs,nb)
!         call cmult(yy(1,ii),1.0/ys,nb)
!      enddo
!
!      do ii=1,nb
!         call add4(B(1,ii),B(1,ii),w3(1,ii)
!     $            ,yy(1,ii),nb)
!      enddo

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
      
c      if (nio.eq.0) write(6,*) 'inside invHessian_update'

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
      subroutine BFGS(rhs)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real B_qn(nb,nb)
      real qgo(nb), qngradf(nb), ngf
      real fo,qnf,qndf
      real yy(nb,nb),ys,sBs
      real ww(nb), pert
      real uu(nb), rhs(nb)
      real alphak

      integer par_step, jmax
      integer chekbc ! flag for checking boundary

      real tmp(nb,nb),tmp1(nb,nb),tmp2(nb,nb),tmp3(nb,nb)
      real tmp4(nb),tmp5(nb),tmp6(nb,nb),tmp7(nb,nb)

c      if (nio.eq.0) write (6,*) 'inside BFGS'

      jmax = 0

      call copy(uu,u(1,1),nb)
   
      ! parameter for barrier function
      par = 1e-3 
      par_step = 4
   
      ! invhelm for computing qnf
      if (ad_step.le.3) then 
         call copy(invhelm(1,1),helm(1,1),nb*nb)
         call lu(invhelm,nb,nb,irv,icv)
      endif

      ! BFGS method with barrier function starts
      do k=1,par_step

         chekbc = 0

c        use helm from BDF3/EXT3 as intial approximation
         call copy(B_qn(1,1),helm(1,1),nb*nb)

         call comp_qnf(uu,rhs,qnf)
         call comp_qngradf(uu,rhs,qngradf)

c        compute quasi-Newton step
         do j=1,100

            call copy(tmp3(1,1),B_qn(1,1),nb*nb)
            call lu(tmp3,nb,nb,irv,icv)
            call copy(qns,qngradf,nb)
            call chsign(qns,nb)
            call solve(qns,tmp3,1,nb,nb,irv,icv)

c            call add2(uu,qns,nb)
            call backtrackr(uu,qns,rhs,1e-2,0.5,alphak)

            ! check the boundary 
            do ii=1,nb
               if ((uu(ii)-sample_max(ii)).ge.1e-8) then
                  chekbc = 1
                  uu(ii) = sample_max(ii) - 0.1*sam_dis(ii)
               elseif ((sample_min(ii)-uu(ii)).ge.1e-8) then
                  chekbc = 1
                  uu(ii) = sample_min(ii) + 0.1*sam_dis(ii)
               endif
            enddo

            call copy(qgo,qngradf,nb) ! store old qn-gradf
            call comp_qngradf(uu,rhs,qngradf)        ! update qn-gradf
            call sub3(qny,qngradf,qgo,nb) 

            ! update approximate Hessian by two rank-one update if chekbc = 0
            ! first rank-one update
            call Hessian_update(B_qn,qns,qny,nb)

            call copy(ww,qngradf,nb)
            call solve(ww,invhelm,1,nb,nb,irv,icv)

            ngf = glsc2(ww,qngradf,nb)
            ngf = sqrt(ngf)

            fo = qnf      ! store old qn-f
            call comp_qnf(uu,rhs,qnf) ! update qn-f
            qndf = abs(qnf-fo)/abs(fo) 
c            write(6,*)'f and old f',j,qnf,fo,qndf,ngf

            ! reset chekbc 
            chekbc = 0
            
            jmax = max(j,jmax)

            if (ngf .lt. 1e-4 .OR. qndf .lt. 1e-6  ) then 
               exit
            endif

c     update solution
         enddo
         par = par*0.1
      enddo

      if (mod(ad_step,ad_iostep).eq.0) then
         if (nio.eq.0) then
            write(6,*)'ad_step, par, jmax, ngf, qndf:',
     $               ad_step,par,jmax,ngf,qndf 
         endif
      endif

      call copy(rhs,uu,nb)

c      if (nio.eq.0) write (6,*) 'exitting BFGS'

      return
      end
c-----------------------------------------------------------------------
      subroutine backtrackr(uu,s,rhs,sigmab,facb,alphak)

      include 'SIZE'
      include 'MOR'
      include 'TOTAL'

      real rhs(nb), s(nb)
      real uuo(nb), uu(nb)
      real Jfk(nb)
      real fk, fk1
      real Jfks
      integer chekbc, counter
      real sigmab, facb, alphak

      alphak = 1.0
      chekbc = 1
      counter = 0

      call comp_qnf(uu,rhs,fk) ! get old f
      call comp_qngradf(uu,rhs,Jfk)

      call copy(uuo,uu,nb)
      call add2s1(uu,s,alphak,nb)

      call comp_qnf(uu,rhs,fk1) ! get new f
      Jfks = vlsc2(Jfk,s,nb)

      do while ((fk1 > fk + sigmab * alphak * Jfks) .OR. (chekbc.eq.1))
         counter = counter + 1
         alphak = alphak * facb
         call add3s2(uu,uuo,s,1.0,alphak,nb)

         chekbc = 0
         do ii=1,nb
            if ((uu(ii)-sample_max(ii)).ge.1e-8) then
               chekbc = 1
            elseif ((sample_min(ii)-uu(ii)).ge.1e-8) then
               chekbc = 1
            endif
         enddo

         call comp_qnf(uu,rhs,fk1)

         if (alphak < 1e-4) then
            exit
         endif
      enddo

      return
      end
