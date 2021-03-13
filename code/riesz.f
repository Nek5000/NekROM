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
               call exitti('Buoyancy in EI disabled for now...l',1)
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
      subroutine set_riesz_rep

      ! compute the riesz representation (only one)
      ! solve one stoke problem and one poisson problem

      include 'SIZE'
      include 'TOTAL'
      include 'ORTHOT'
      include 'CTIMER'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real work(lx2*ly2*lz2*lelt)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside set_riesz_rep'

      call rone(ones,n)
      call rzero(zeros,n)

      if (ifrom(1)) then
         ifield=1
         tolhv=1e-8

         call steady_stoke_solve(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
     $            work,res_u(1,1),res_u(1,2),res_u(1,ldim))
      endif

      if (ifrom(2)) then
         ifield=2
         tolht(2)=1e-8
         ifld1 = ifield-1
         napproxt(1,ifld1) = laxtt
         ! check with imesh
         imesh  = 1

         call hsolve  ('TEMP',eh_t,res_t,ones,ones
     $            ,tmask(1,1,1,1,ifield-1)
     $            ,tmult(1,1,1,1,ifield-1)
     $            ,imesh,tolht(ifield),nmxh,1
     $            ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
      endif

      call nekgsync

      if (nio.eq.0) write (6,*) 'exiting set_riesz_rep'

      return
      end
c-----------------------------------------------------------------------
