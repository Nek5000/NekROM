      subroutine BFGS_freeze

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real B_qn(nb,nb)
      real go(nb),fo,qndf
      real tmp(nb,nb),tmp1(nb,nb),tmp2(nb,nb),tmp3(nb,nb)
      real tmp4(nb),tmp5(nb),tmp6(nb,nb),tmp7(nb,nb)
      real yy(nb,nb),ys,sBs
      real ww(nb), ngf, pert
      integer par_step
      integer chekbc ! flag for checking boundary

      if (nio.eq.0) write (6,*) 'inside BFGS_freeze'
   
c     parameter for barrier function
      par = 1e-3 
      par_step = 4
   
c     invhelm for computing qnf
      if (ad_step.le.3) then 
         call copy(invhelm(1,1),helm(1,1),nb*nb)
         call lu(invhelm,nb,nb,ir,ic)
      endif

c     BFGS method with barrier function starts
      do k=1,par_step

         chekbc = 0

c        use helm from BDF3/EXT3 as intial approximation
         call copy(B_qn(1,1),helm(1,1),nb*nb)
!         call copy(B_qn(1,1),invhelm(1,1),nb*nb)

         call comp_qnf
         call comp_qngradf

c        compute quasi-Newton step
         do j=1,500

            call copy(tmp3(1,1),B_qn(1,1),nb*nb)
            call lu(tmp3,nb,nb,ir,ic)
            call copy(qns,qngradf,nb)
            call chsign(qns,nb)
            call solve(qns,tmp3,1,nb,nb,ir,ic)
c            if (j .eq. 1) then
c               call copy(qns,qngradf,nb)
c               call solve(qns,B_qn,1,nb,nb,ir,ic)
c            else     
c               call mxm(B_qn,nb,qngradf,nb,qns,1)
c            endif

c            call chsign(qns,nb)
            call add2(u(1,1),qns,nb)

            ! check the boundary 
            do ii=1,nb
               if ((u(ii,1)-sample_max(ii)).ge.1e-8) then
                  chekbc = 1
                  u(ii,1) = sample_max(ii) - 0.1*sam_dis(ii)
               elseif ((sample_min(ii)-u(ii,1)).ge.1e-8) then
                  chekbc = 1
                  u(ii,1) = sample_min(ii) + 0.1*sam_dis(ii)
               endif
            enddo

            call copy(go,qngradf,nb) ! store old qn-gradf
            call comp_qngradf        ! update qn-gradf
            call sub3(qny,qngradf,go,nb) 

            ! update approximate Hessian by two rank-one update if chekbc = 0
            ! first rank-one update
            if (chekbc .NE. 1) then  

               call Hessian_update(B_qn,qns,qny,nb)
!               call invHessian_update(B_qn,qns,qny,nb)

            endif   

            call copy(ww,qngradf,nb)
            call solve(ww,invhelm,1,nb,nb,ir,ic)

            ngf = glsc2(ww,qngradf,nb)
            ngf = sqrt(ngf)

            fo = qnf      ! store old qn-f
            call comp_qnf ! update qn-f
            qndf = abs(qnf-fo)/abs(fo) 
c            write(6,*)'f and old f',j,qnf,fo,qndf,ngf

            ! reset chekbc 
            chekbc = 0

            if (ngf .lt. 1e-4 .OR. qndf .lt. 1e-6  ) goto 900

c     update solution
         enddo
  900    write(6,*)'criterion reached, number of iteration:'
     $              ,ad_step,j,par,ngf,qndf 
         par = par*0.1
c         par = 2.**(-k)
      enddo

      if (nio.eq.0) write (6,*) 'exitting BFGS_freeze'

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_qngradf
      
      include 'SIZE'
      include 'MOR'

      real tmp1(nb),tmp2(nb),tmp3(nb),tmp4(nb)
      real mpar
      integer barr_func

      barr_func = 1

!      if (nio.eq.0) write (6,*) 'inside com_qngradf'

      if (barr_func .eq. 1) then ! use logarithmic as barrier function

         call sub3(tmp1,u(1,1),sample_max,nb)  
         call sub3(tmp2,u(1,1),sample_min,nb)  
   
         ! add perturbation 
         call cadd(tmp1,-1e-2,nb)
         call cadd(tmp2,1e-2,nb)
   
         call invcol1(tmp1,nb)
         call invcol1(tmp2,nb)
         call add3(tmp3,tmp1,tmp2,nb)

         mpar = -1.0*par
         call add3s12(qngradf,opt_rhs(1),tmp3,-1.0,mpar,nb)
   
         ONE = 1.
         ZERO= 0.
         call dgemv('N',nb,nb,ONE,helm,nb,u(1,1),1,ZERO,tmp4,1)
         call add2(qngradf,tmp4,nb)

      else ! use inverse function as barrier function

         call sub3(tmp1,u(1,1),sample_max,nb)  
         call sub3(tmp2,sample_min,u(1,1),nb)  
         call vsq(tmp1,nb)
         call vsq(tmp2,nb)
         call invcol1(tmp1,nb)
         call invcol1(tmp2,nb)
         call sub3(tmp3,tmp2,tmp1,nb)

         mpar = -1.0*par
         call add3s12(qngradf,opt_rhs(1),tmp3,-1.0,mpar,nb)
   
         ONE = 1.
         ZERO= 0.
         call dgemv('N',nb,nb,ONE,helm,nb,u(1,1),1,ZERO,tmp4,1)
         call add2(qngradf,tmp4,nb)

      endif

!      if (nio.eq.0) write (6,*) 'exiting com_qngradf'

      return 
      end
c-----------------------------------------------------------------------
      subroutine comp_qnf
      
      include 'SIZE'
      include 'MOR'

      real tmp1(nb),tmp2(nb),tmp3(nb)
      real tmp4(nb),tmp5(nb),tmp6(nb)
      real term1,term2,term3,term4
      real bar1,bar2 ! bar1 and bar2 are the barrier function for two constrains
      integer barr_func

      barr_func = 1

!      if (nio.eq.0) write (6,*) 'inside com_qnf'

c     evaluate quasi-newton f
      ONE = 1.
      ZERO= 0.

      call dgemv('N',nb,nb,ONE,helm,nb,u(1,1),1,ZERO,tmp6,1) ! H*coef
      term1 = 0.5 * glsc2(tmp6,u(1,1),nb) ! coef'*H*coef

      term2 = glsc2(u(1,1),opt_rhs(1),nb) ! coef'*rhs

      ! 0.5*rhs'*inv(H)*rhs
      call copy(tmp5,opt_rhs(1),nb)
      call solve(tmp5,invhelm,1,nb,nb,ir,ic)
      term3 = 0.5 * glsc2(opt_rhs(1),tmp5,nb)

      if (barr_func .eq. 1) then ! use logarithmetic as barrier function

         ! barrier term
         call sub3(tmp1,sample_max,u(1,1),nb)  
         call sub3(tmp2,u(1,1),sample_min,nb)  
   
         ! add perturbation
         call cadd(tmp1,1e-2,nb)
         call cadd(tmp2,1e-2,nb)
   
         do i=1,nb
            tmp3(i) = log(tmp1(i))
            tmp4(i) = log(tmp2(i))
         enddo

      else 

         call sub3(tmp1,u(1,1),sample_max,nb)  
         call sub3(tmp2,sample_min,u(1,1),nb)  
   
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
