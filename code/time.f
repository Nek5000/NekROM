c-----------------------------------------------------------------------
      subroutine rom_step_heat

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

c     Matrices and vectors for advance
      real tmp(0:nb),tmat(nb,nb+1),rhs(0:nb)
      real coef(1:nb), e0(0:nb)

      common /scrk3/ work(lt)
      common /scrk1/ t1(lt),t2(lt),t3(lt)

c     Variable for vorticity
      real vort(lt,3)

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      if (ad_step.eq.1) then
         step_time = 0.
      endif

      last_time = dnekclock()

      n  = lx1*ly1*lz1*nelt

      ad_step=max(ad_step,1)

      count = min(ad_step,3)

      call rzero(e0,nb+1)
      e0(0) = 1

c     write (6,*) 'ad_beta',count,ad_beta(1,count)

      if (ad_step.le.3) then
         call cmult2(flu,b,ad_beta(1,count)/ad_dt,nb*nb)
         call add2s2(flu,a,1/ad_re,nb*nb)
      endif

      ONE = 1.
      ZERO= 0.

      call rzero(rhs,nb+1)

      call mxm(u,nb+1,ad_beta(2,count),3,tmp,1)
c     call mxm(b0,nb+1,tmp,nb+1,rhs,1)
      call mxm(b,nb,tmp(1),nb,rhs(1),1)
      call cmult(rhs,-1.0/ad_dt,nb+1)
c     do i=0,nb
c        write (6,*) 'rhs',i,rhs(i)
c     enddo

      s=-1.0/ad_re

c     call add2s2(rhs,a0,s,nb+1) ! not working...
      do i=0,nb
         rhs(i)=rhs(i)+s*a0(i,0)
c        write (6,*) 'rhs',i,rhs(i)
      enddo

c     call copy(conv(1,3),conv(1,2),nb)
c     call copy(conv(1,2),conv(1,1),nb)

c     call evalc(conv)

c     call mxm(conv,nb,ad_alpha(1,count),3,tmp(1),1)

c     call sub2(rhs,tmp,nb+1)

      if (ad_step.le.3) call lu(flu,nb,nb,ir,ic)

      call solve(rhs(1),flu,1,nb,nb,ir,ic)
      call shiftu(rhs(1))

c     call drago
c     call comp_rms

      step_time=step_time+dnekclock()-last_time

      if (ifheat) call recon(vx,vy,vz,u)

      if (ifravg) then
         call recon(vx,vy,vz,u)
         call avg_all
      endif

      if (mod(ad_step,ad_iostep).eq.0) then

!        This output is to make sure the ceof matches with matlab code

         if (nio.eq.0) then
            write (6,*)'ad_step:',ad_step,ad_iostep,npp,nid,step_time
            if (ad_step.eq.ad_nsteps) then
               do j=1,nb
                  write(6,*) j,u(j,1),'final'
               enddo
               write (6,*) 'step_time: ',step_time
            else
               do j=1,nb
                  write(6,*) j,u(j,1)
               enddo
            endif
         endif

         if (ifdump) then
            idump=ad_step/ad_iostep
            if (.not.ifravg) call recon(vx,vy,vz,u)

            ! compute the vorticity of the ROM reconstructed field
            call opcopy(t1,t2,t3,vx,vy,vz)
            call comp_vort3(vort,work1,work2,t1,t2,t3)
            ifto = .true. ! turn on temp in fld file
            call outpost(vx,vy,vz,pr,vort,'rom')
            if (nio.eq.0) write (6,*) 'inside ifdump'
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_step

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

c     Matrices and vectors for advance
      real tmp(0:nb),tmat(nb,nb+1),rhs(0:nb)
      real coef(1:nb), e0(0:nb)

      common /scrk3/ work(lt)
      common /scrk1/ t1(lt),t2(lt),t3(lt)

c     Variable for vorticity
      real vort(lt,3)

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      if (ad_step.eq.1) then
         step_time = 0.
         call rzero(ua,nb+1)
         call rzero(u2a,(nb+1)**2)
      endif

      last_time = dnekclock()

      n  = lx1*ly1*lz1*nelt

      count = min0(ad_step,3)

      call rzero(e0,nb+1)
      e0(0) = 1

      if (ad_step.le.3) then
         call cmult2(flu,b,ad_beta(1,count)/ad_dt,nb*nb)
         call add2s2(flu,a,1/ad_re,nb*nb)
      endif

      ONE = 1.
      ZERO= 0.

      call mxm(u,nb+1,ad_beta(2,count),3,tmp,1)
c     call mxm(b0,nb+1,tmp,nb+1,rhs,1)
      call mxm(b,nb,tmp(1),nb,rhs(1),1)

      call cmult(rhs,-1.0/ad_dt,nb+1)

      s=-1.0/ad_re

c     call add2s2(rhs,a0,s,nb+1) ! not working...
      do i=0,nb
         rhs(i)=rhs(i)+s*a0(i,0)
      enddo

      call copy(conv(1,3),conv(1,2),nb)
      call copy(conv(1,2),conv(1,1),nb)

      call evalc(conv)

      call mxm(conv,nb,ad_alpha(1,count),3,tmp(1),1)

      call sub2(rhs,tmp,nb+1)

      if (ad_step.le.3) call lu(flu,nb,nb,ir,ic)

      call solve(rhs(1),flu,1,nb,nb,ir,ic)

      call shiftu(rhs(1))
      call add2(ua,u,nb+1)

      do j=0,nb
      do i=0,nb
         u2a(i,j)=u2a(i,j)+u(i,1)*u(j,1)
      enddo
      enddo

      if (ad_step.eq.3) call copy(uj,u,3*(nb+1))
      if (ad_step.eq.ad_nsteps) then
         call copy(uj(0,4),u,3*(nb+1))
         do k=1,6
         do j=0,nb
         do i=0,nb
            u2j(i,j,k)=uj(i,k)*uj(j,k)
         enddo
         enddo
         enddo
         rinsteps=1./real(ad_nsteps)
         call cmult(ua,rinsteps,nb+1)
         call cmult(u2a,rinsteps,(nb+1)**2)
      endif

      call comp_drag
c     call comp_rms ! old

      step_time=step_time+dnekclock()-last_time

      if (ifheat) call recon(vx,vy,vz,u)

      if (ifravg) then
         call recon(vx,vy,vz,u)
         call avg_all
      endif

      if (mod(ad_step,ad_iostep).eq.0) then

!        This output is to make sure the ceof matches with matlab code

         if (nio.eq.0) then
            write (6,*)'ad_step:',ad_step,ad_iostep,npp,nid,step_time
            if (ad_step.eq.ad_nsteps) then
               do j=1,nb
                  write(6,*) j,u(j,1),'final'
               enddo
               write (6,*) 'step_time: ',step_time
            else
               do j=1,nb
                  write(6,*) j,u(j,1)
               enddo
            endif
         endif

         if (.true.) then
            idump=ad_step/ad_iostep
            if (.not.ifravg) call recon(vx,vy,vz,u)

            ! compute the vorticity of the ROM reconstructed field
            call opcopy(t1,t2,t3,vx,vy,vz)
            call comp_vort3(vort,work1,work2,t1,t2,t3)
            ifto = .true. ! turn on temp in fld file
            call outpost(vx,vy,vz,pr,vort,'rom')
            if (nio.eq.0) write (6,*) 'inside ifdump'
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_BDF_coef(ad_alpha,ad_beta)

      real ad_alpha(3,3), ad_beta(4,3)

      call rzero(ad_alpha,3*3)
      call rzero(ad_beta,3*4)

      ad_beta(1,1) = 1.
      ad_beta(2,1) = -1.

      ad_beta(1,2) = 1.5
      ad_beta(2,2) = -2
      ad_beta(3,2) = 0.5

      ad_beta(1,3) = 11./6
      ad_beta(2,3) = -3
      ad_beta(3,3) = 1.5
      ad_beta(4,3) = -1./3.

      ad_alpha(1,1)=1

      ad_alpha(1,2)=2
      ad_alpha(2,2)=-1

      ad_alpha(1,3)=3
      ad_alpha(2,3)=-3
      ad_alpha(3,3)=1

      return
      end
c-----------------------------------------------------------------------
      subroutine evalc(cu)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real cu(nb)

      integer icalld
      save    icalld
      data    icalld /0/

      common /scrk4/ work(lx1*ly1*lz1*lelt)

      if (icalld.eq.0) then
         evalc_time=0.
         icalld=1
      endif

      stime=dnekclock()

      l=1

      call rzero(cu,nb)

      do l=1,ncloc
         i=icloc(1,l)
         j=icloc(2,l)
         k=icloc(3,l)
         cu(i)=cu(i)+clocal(l)*u(j,1)*u(k,1)
      enddo

      call gop(cu,work,'+  ',nb)

      evalc_time=evalc_time+dnekclock()-stime

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_const
c This subroutine is solving rom with constrains
c The subroutine is based on BFGS method with barrier function
      
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

c     Matrices and vectors for advance
      real tmp(0:nb),tmat(nb,nb+1)
      real coef(1:nb)

      common /scrk3/ work(lt)
      common /scrk1/ t1(lt),t2(lt),t3(lt)

c     Variable for vorticity
      real vort(lt,3)

c     Working arrays for LU

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      n  = lx1*ly1*lz1*nelt

      time=time+ad_dt

      count = min0(ad_step,3)

      if (ad_step.le.3) then
         call cmult2(helm,b,ad_beta(1,count)/ad_dt,nb*nb)
         call add2s2(helm,a,1/ad_re,nb*nb)
      endif

      ONE = 1.
      ZERO= 0.

      call mxm(u,nb+1,ad_beta(2,count),3,tmp,1)
      call mxm(b,nb,tmp(1),nb,opt_rhs(1),1)

      call cmult(opt_rhs,-1.0/ad_dt,nb+1)

      s=-1.0/ad_re

c     call add2s2(rhs,a0,s,nb+1) ! not working...
      do i=0,nb
         opt_rhs(i)=opt_rhs(i)+s*a0(i,0)
      enddo
c      call add2s2(rhs,a0(1,0),-1/ad_re,nb)

      call copy(conv(1,3),conv(1,2),nb)
      call copy(conv(1,2),conv(1,1),nb)

      if (np.eq.1) then
         call mxm(c,nb*(nb+1),u,nb+1,tmat,1)
         call mxm(tmat,nb,u,nb+1,conv,1)
      else
         call evalc(conv)
      endif

      call mxm(conv,nb,ad_alpha(1,count),3,tmp(1),1)

      call sub2(opt_rhs,tmp,nb+1)

      call copy(u(1,3),u(1,2),nb)
      call copy(u(1,2),u(1,1),nb)

      call opt_const

      if (mod(ad_step,ad_iostep).eq.0) then

!        This output is to make sure the ceof matches with matlab code

         if (nio.eq.0) then
            write (6,*)'ad_step:',ad_step,ad_iostep,npp,nid
            if (ad_step.eq.ad_nsteps) then
               do j=1,nb
                  write(6,*) j,u(j,1),'final'
               enddo
            else
               do j=1,nb
                  write(6,*) j,u(j,1)
               enddo
            endif
         endif

         if (ifdump) then
            call opzero(vx,vy,vz)
            do j=1,nb
               call opadds(vx,vy,vz,ub(1,j),vb(1,j),wb(1,j),coef(j),n,2)
            enddo
            call opadd2  (vx,vy,vz,ub,vb,wb)

            ! compute the vorticity of the ROM reconstructed field
            call opcopy(t1,t2,t3,vx,vy,vz)
            call comp_vort3(vort,work1,work2,t1,t2,t3)
            ifto = .true. ! turn on temp in fld file
            call copy(t,vort,n)

            call outpost (vx,vy,vz,pr,t,'rom')
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine opt_const

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

c     parameter for barrier function
c     it should start from value greater than one and decrease
      real B_qn(nb,nb)
      real go(nb),fo,qndf
      real tmp(nb,nb),tmp1(nb,nb),tmp2(nb,nb),tmp3(nb,nb)
      real tmp4(nb),tmp5(nb),tmp6(nb,nb),tmp7(nb,nb)
      real yy(nb,nb),ys,sBs
      integer par_step

      if (nio.eq.0) write (6,*) 'inside opt_const'

      par_step = 1
c      par = 1.0 
      par = 0.01 


c     invhelm for computing qnf
c     not changing during BFGS
c     only changes during the time stepper
      if (ad_step.le.3) then 
         call copy(invhelm(1,1),helm(1,1),nb*nb)
         call lu(invhelm,nb,nb,ir,ic)
      endif

c     BFGS method with barrier function starts
      do k=1,par_step

c     use helm from BDF3/EXT3 as intial approximation
         do i=1,nb
            call copy(B_qn(1,i),helm(1,i),nb)
         enddo

         call comp_qnf
         call comp_qngradf

c     compute quasi-Newton step
         do j=1,500

            call copy(tmp3(1,1),B_qn(1,1),nb*nb)
            call lu(tmp3,nb,nb,ir,ic)
            call copy(qns,qngradf,nb)
            call chsign(qns,nb)
            call solve(qns,tmp3,1,nb,nb,ir,ic)
            call add2(u(1,1),qns,nb)

            ! check whether solution exceed the boundary
            do ii=1,nb
               if ((u(ii,1)-sample_max(ii)).ge.1e-10) then
                  u(ii,1) = 0.9 * sample_max(ii)
               elseif ((sample_min(ii)-u(ii,1)).ge.1e-10) then
                  u(ii,1) = 0.9 * sample_min(ii)
               endif
            enddo

            call copy(go,qngradf,nb) ! store old qn-gradf
            call comp_qngradf        ! update qn-gradf
            call sub3(qny,qngradf,go,nb)

            ! update approximate Hessian by two rank-one update
            ! first rank-one update
            ! outer product: s_k * s_k^T               
            call mxm(qns,nb,qns,1,tmp,nb)
             
            ! s_k * s_k^T * B_k
            call mxm(tmp,nb,B_qn,nb,tmp1,nb)

            ! B_k * s_k * s_k^T * B_k 
            call mxm(B_qn,nb,tmp1,nb,tmp2,nb)

            ! s_k^T * B_k * s_k 
            call mxm(B_qn,nb,qns,nb,tmp5,1)
            sBs = glsc2(qns,tmp5,nb)

            ! second rank-one update
            ! outer product: y_k * y_k^T               
            call mxm(qny,nb,qny,1,yy,nb)

            ys = glsc2(qny,qns,nb)

            do ii=1,nb
               call cmult(tmp2(1,ii),-1.0/sBs,nb)
               call cmult(yy(1,ii),1.0/ys,nb)
            enddo

            do ii=1,nb
               call add4(B_qn(1,ii),B_qn(1,ii),tmp2(1,ii),yy(1,ii),nb)
            enddo

            fo = qnf      ! store old qn-f
            call comp_qnf ! update qn-f
            qndf = abs(qnf-fo) 
            write(6,*)'f and old f',qnf,fo,qndf

            if (qndf .lt. 1e-10) goto 900

c     update solution
         enddo
  900    write(6,*)'criterion reached, number of iteration:',j,par 
         par = par*0.1

      enddo
      if (nio.eq.0) write (6,*) 'exitting opt_const'

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_qngradf
      
      include 'SIZE'
      include 'MOR'
      real tmp1(nb),tmp2(nb),tmp3(nb),tmp4(nb)
      real mpar

!      if (nio.eq.0) write (6,*) 'inside com_qngradf'

      call sub3(tmp1,u(1,1),sample_max,nb)  
      call sub3(tmp2,u(1,1),sample_min,nb)  
      call invcol1(tmp1,nb)
      call invcol1(tmp2,nb)

      call add3(tmp3,tmp1,tmp2,nb)

      mpar = -1.0*par

      call add3s12(qngradf,opt_rhs(1),tmp3,-1.0,mpar,nb)

      ONE = 1.
      ZERO= 0.
      call dgemv( 'N',nb,nb,ONE,helm,nb,u(1,1),1,ZERO,tmp4,1)
      call add2(qngradf,tmp4,nb)


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

c     bar1 and bar2 are the barrier function for two constrains
      real bar1,bar2

!      if (nio.eq.0) write (6,*) 'inside com_qnf'

c     evaluate quasi-newton f

c     term1 represents 0.5*x'*H*x
      ONE = 1.
      ZERO= 0.
      call dgemv( 'N',nb,nb,ONE,helm,nb,u(1,1),1,ZERO,tmp6,1)
      term1 = 0.5 * glsc2(tmp6,u(1,1),nb)
      write(6,*)'term1',term1

c     term2 represents x'*g
      term2 = glsc2(u(1,1),opt_rhs(1),nb)
      write(6,*)'term2',term2

c     term3 represents 0.5*g'*inv(H)*g
      call copy(tmp5,opt_rhs(1),nb)
      call solve(tmp5,invhelm,1,nb,nb,ir,ic)
      term3 = 0.5 * glsc2(opt_rhs(1),tmp5,nb)
      write(6,*)'term3',term3

c     barrier term
      call sub3(tmp1,sample_max,u(1,1),nb)  
      call sub3(tmp2,u(1,1),sample_min,nb)  

c     currently can only come up with this way to compute log for an array
      do i=1,nb
         tmp3(i) = log(tmp1(i))
         tmp4(i) = log(tmp2(i))
      enddo

      bar1 = vlsum(tmp3,nb)
      bar2 = vlsum(tmp4,nb)
      term4 = par*(bar1+bar2)
      write(6,*)'term4',term4

      qnf = term1 - term2 + term3 - term4

!      if (nio.eq.0) write (,*) 'exitting com_qnf'

      return
      end
c-----------------------------------------------------------------------
      subroutine evalc_const(cu)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real cu(nb)

      integer icalld
      save    icalld
      data    icalld /0/

      common /scrk4/ work(lx1*ly1*lz1*lelt)

      if (icalld.eq.0) then
         evalc_time=0.
         icalld=1
      endif

      stime=dnekclock()

      l=1

      call rzero(cu,nb)

      do l=1,ncloc
         i=icloc(1,l)
         j=icloc(2,l)
         k=icloc(3,l)
         cu(i)=cu(i)+clocal(l)*u(j,1)*u(k,1)
         write (6,*) i,j,k,u(j,1),u(k,1)
      enddo

      call exitt0

      call gop(cu,work,'+  ',nb)

      evalc_time=evalc_time+dnekclock()-stime

      return
      end
c-----------------------------------------------------------------------
