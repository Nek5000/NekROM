c-----------------------------------------------------------------------
      subroutine rom_step_heat

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

c     Matrices and vectors for advance
      real tmp(0:nb),tmat(nb,nb+1),rhs(0:nb)
      real coef(1:nb), e0(0:nb)

      common /scrkrstep/ t1(lt),t2(lt),t3(lt),work(lt)

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

      s=-1.0/ad_re

c     call add2s2(rhs,a0,s,nb+1) ! not working...
      do i=0,nb
         rhs(i)=rhs(i)+s*a0(i,0)
c        write (6,*) 'rhs',i,rhs(i)
      enddo

      call copy(conv(1,3),conv(1,2),nb)
      call copy(conv(1,2),conv(1,1),nb)

      call evalc(conv)

      call mxm(conv,nb,ad_alpha(1,count),3,tmp(1),1)
c     do i=0,nb
c        write (6,*) 'evalc',tmp(i)
c     enddo

      call sub2(rhs,tmp,nb+1)

      if (ad_step.le.3) call lu(flu,nb,nb,ir,ic)

      call solve(rhs(1),flu,1,nb,nb,ir,ic)
      call shiftu(rhs(1))

      if (ifheat) call recon(vx,vy,vz,u)

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

         if (rmode.eq.'ALL'.or.rmode.eq.'ONB') then
            idump=ad_step/ad_iostep
            call recon(vx,vy,vz,u)

            ! compute the vorticity of the ROM reconstructed field
            call opcopy(t1,t2,t3,vx,vy,vz)
            call comp_vort3(vort,work1,work2,t1,t2,t3)
            ifto = .true. ! turn on temp in fld file
            call outpost(vx,vy,vz,pr,vort,'rom')
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_setup_heat

      include 'SIZE'
      include 'MOR'

      if (nio.eq.0) write (6,*) 'inside rom_setup'

      setup_start=dnekclock()

      call rom_init_params
      call rom_init_fields

      call setbases_heat

      call seta
      call setb
      call setc_const
      call setu
c     call setops

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ')
     $   call dump_all

c     call qoisetup

      setup_end=dnekclock()

      if (nio.eq.0) write (6,*) 'exiting rom_setup'
      if (nio.eq.0) write (6,*) 'setup_time:', setup_end-setup_start

      return
      end
c-----------------------------------------------------------------------
      subroutine setc_const

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real cux(lt),cuy(lt),cuz(lt)

      common /scrk1/ t1(lt),binv(lt),wk1(lt),wk2(lt),wk3(lt)
      common /scrcwk/ wk(lcloc)

      conv_time=dnekclock()

      call invers2(binv,bm1,lx1*ly1*lz1*nelv)
      call rone(binv,lx1*ly1*lz1*nelv)

      if (nio.eq.0) write (6,*) 'inside setc'

      l=0
      n=lx1*ly1*lz1*nelv

      do k=0,0
         call setcnv_c(ub(1,k),vb(1,k),wb(1,k))
         do j=0,nb
            call setcnv_u(ub(1,j),vb(1,j),wb(1,j))
            call ccu(cux,cuy,cuz)
            do i=1,nb
               l=l+1
               clocal(l) = 
     $            op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),cux,cuy,cuz,binv)
               icloc(1,l) = i
               icloc(2,l) = j
               icloc(3,l) = k
            enddo
         enddo
      enddo
      ncloc=l

      if (nio.eq.0) write (6,*) 'exiting setc'

      return
      end
c-----------------------------------------------------------------------
      subroutine setbases_heat

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      n=lx1*ly1*lz1*nelv

      do ib=0,nb
         call rzero(ub(1,ib),n)
         call rzero(vb(1,ib),n)
         call rzero(wb(1,ib),n)
      enddo

      one=1.
      pi=4.*atan(one)

      do i=1,n
         x=xm1(i,1,1,1)
         y=ym1(i,1,1,1)
         ub(i,0)=1.
         k=2
         ub(i,1)=sin(k*pi*x)*sin(k*pi*y)
         if (nb.ge.2) ub(i,2)=cos(k*pi*x)*sin(k*pi*y)
         if (nb.ge.3) ub(i,3)=sin(k*pi*x)*cos(k*pi*y)
         if (nb.ge.4) ub(i,4)=cos(k*pi*x)*cos(k*pi*y)
         if (nb.ge.4) vb(i,0)=0.5
      enddo

      return
      end
c-----------------------------------------------------------------------
