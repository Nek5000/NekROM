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
            call set_xi_a(xi(1,l),tb(1,1,1),ones,1,nb,2)
            l=l+nb
            call set_xi_b(xi(1,l),qq,1,1,2)
            l=l+1
            do i=1,nb
               call set_gradn(wk,tb(1,i,1))
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
            call set_xi_b(xi(1,l),tb,1,nb+1,2)
            l=l+nb+1
            call set_xi_a(xi(1,l),tb,ones,1,nb+1,2)
            l=l+nb+1
            call set_xi_b(xi(1,l),qq,1,1,2)
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
            call set_xi_b(xi(1,l),tb,1,nb+1,2)
            l=l+nb+1
            if (ifrom(1)) then
               if (ifaxis) then
                  call push_op(vx,vy,vz)
                  do j=0,nb
                     call opcopy(vx,vy,vz,ub(1,j),vb(1,j),wb(1,j))
                     do i=0,nb
                        call conv1d(xi(1,l),tb(1,i,1))
                        l=l+1
                     enddo
                  enddo
                  call pop_op(vx,vy,vz)
               else
                  call set_xi_c(xi(1,l),ub,vb,wb,tb,1,nb+1,nb+1,2)
                  l=l+(nb+1)**2
               endif
            else
               call set_xi_c(xi(1,l),ub,vb,wb,tb,1,1,nb+1,2)
               l=l+nb+1
            endif
            call set_xi_a(xi(1,l),tb,ones,1,nb+1,2)
            l=l+nb+1
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
                  call opcolv(wk1,wk2,wk3,tb(1,i,1))
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

c     Compute the non divergence-free Riesz representations for steady
c     Boussinesq incompressible NSE
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
      subroutine set_sns_divfrr

c     Compute the divergence-free
c     Riesz representations for steady Boussinesq
c     incompressible NSE (quadratically nonlinear elliptic problem)

      include 'SIZE'
      include 'TOTAL'
      include 'ORTHOT'
      include 'CTIMER'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      real work(lx2*ly2*lz2*lelt)

      common /screi/ wk1(lt),wk2(lt),wk3(lt),wk4(lt),wk5(lt)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside set_sns_divfrr'

      call rone(ones,n)
      call rzero(zeros,n)

      l1=1
      l2=1
      do i=0,nb
         ifield=1
         tolhv=1e-8
         tolht(2)=1e-8
         call steady_stokes_solve(xi_u(1,1,l1),xi_u(1,2,l1),
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
         call steady_stokes_solve(xi_u(1,1,l1),xi_u(1,2,l1),
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
         call steady_stokes_solve(xi_u(1,1,l1),xi_u(1,2,l1),
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
            call steady_stokes_solve(xi_u(1,1,l1),xi_u(1,2,l1),
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

      if (nio.eq.0) write (6,*) 'exiting set_sns_divfrr'

      return
      end
c-----------------------------------------------------------------------
      subroutine set_uns_divfrr

c     Compute the divergence free
c     Riesz representators for unsteady Boussinesq
c     incompressible NSE (quadratically nonlinear elliptic problem)

      include 'SIZE'
      include 'TOTAL'
      include 'ORTHOT'
      include 'CTIMER'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      real work(lx2*ly2*lz2*lelt)

      common /screi/ wk1(lt),wk2(lt),wk3(lt),wk4(lt),wk5(lt)

      n=lx1*ly1*lz1*nelv

      if (nio.eq.0) write (6,*) 'inside set_uns_divfrr'

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
         call steady_stokes_solve(xi_u(1,1,l1),xi_u(1,2,l1),
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
         call steady_stokes_solve(xi_u(1,1,l1),xi_u(1,2,l1),
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
         call steady_stokes_solve(xi_u(1,1,l1),xi_u(1,2,l1),
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
         call steady_stokes_solve(xi_u(1,1,l1),xi_u(1,2,l1),
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
            call steady_stokes_solve(xi_u(1,1,l1),xi_u(1,2,l1),
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
      if (nio.eq.0) write (6,*) 'set_uns_divfrr_time:',
     $             dnekclock()-srru_start

      return
      end
c-----------------------------------------------------------------------
      subroutine set_riesz_reps

      ! compute the Riesz representation pieces
      ! solve multiple stoke problems and / or poisson problems

      ! TODO: Need to resolve xi, xi_u, xi_t inconsistency

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      if (nio.eq.0) write (6,*) 'inside set_riesz_reps'

      if (eqn.eq.'POIS') then
         ifield=2
         call set_xi_poisson
      else if (eqn.eq.'HEAT') then
         ifield=2
         call set_xi_heat
      else if (eqn.eq.'ADVE') then
         ifield=2
         call set_xi_ad
      else if (eqn.eq.'VNS ') then
         ifield=1
         call set_xi_ns
      else if (eqn.eq.'NS  ') then
         call set_residual_unsteady
         call set_uNS_divfrr
      else if (eqn.eq.'SNS ') then
         call set_residual
         call set_sNS_divfrr
      endif

      call nekgsync

      if (nio.eq.0) write (6,*) 'exiting set_riesz_reps'

      return
      end
c-----------------------------------------------------------------------
      subroutine set_riesz_one

      ! compute the Riesz representation (only one)
      ! solve one stoke problem and / or one poisson problem

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

         call steady_stokes_solve(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
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
      subroutine set_xi_a(xi_a,sb,dfld,mdim,nxi,jfld)

      ! set Riesz representation corresponding to diffusion term

      ! xi_a := output Riesz representations
      ! sb   := input fields
      ! dfld := diffusivity fields
      ! mdim := # of components of each field in sb
      ! nxi  := # of fields in sb

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)
      real xi_a(lt,mdim,nxi),sb(lt,mdim,nxi),dfld(lt,mdim)

c     jfld=1
c     if (mdim.eq.1) jfld=2

      do ix=1,nxi
      do idim=1,mdim
         call binva(
     $      xi_a(1,idim,ix),sb(1,idim,ix),dfld(1,idim),jfld,idim)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_xi_b(xi_b,sb,mdim,nxi,jfld)

      ! set Riesz representation corresponding to mass term

      ! xi_b := output Riesz representations
      ! sb   := input fields
      ! mdim := # of components of each field in sb
      ! nxi  := # of fields in sb

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)
      real xi_b(lt,mdim,nxi),sb(lt,mdim,nxi)

c     jfld=1
c     if (mdim.eq.1) jfld=2

c     nel=nelv
c     if (jfld.eq.2) nel=nelt

      do ix=1,nxi
      do idim=1,mdim
         call binvb(xi_b(1,idim,ix),sb(1,idim,ix),jfld)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_xi_c(xi_c,cxb,cyb,czb,sb,mdim,mxi,nxi,jfld)

      ! set Riesz representation corresponding to advection term

      ! xi_c        := output Riesz representations
      ! cxb,cyb,czb := input advection field components
      ! sb          := input fields
      ! mdim        := # of components of each field in sb
      ! mxi         := # of fields in cb
      ! nxi         := # of fields in sb

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)
      real xi_c(lt,mdim,nxi,mxi),cxb(lt,mxi),cyb(lt,mxi),czb(lt,mxi),
     $   sb(lt,mdim,nxi)

      do jx=1,mxi
      do ix=1,nxi
      do idim=1,mdim
         call binvc(xi_c(1,idim,ix,jx),cxb(1,jx),cyb(1,jx),czb(1,jx),
     $     sb(1,idim,ix),.false.,.false.,jfld)
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine binva(bia,s,dfld,jfld,idir)

      ! apply local B^{-1} A to a scalar field

      ! s    := input scalar field
      ! dfld := diffusivity field
      ! jfld := field # of s
      ! idir := direction of s
      ! bia  := B^{-1} A s

      include 'SIZE'
      include 'MASS'     ! dep: bm1
      include 'MOR'      ! dep: zeros
      include 'TSTEP'    ! dep: nelfld
      include 'PARALLEL' ! dep: nelg,nelgv

      parameter (lt=lx1*ly1*lz1*lelt)
      real bia(lt,1),s(lt,1),dfld(lt,1)

      nel=nelfld(jfld)
      melg=nelg(jfld)

      imsh=2
      if (melg.eq.nelgv) imsh=1 

      mfld=1
      if (jfld.eq.1) mfld=ldim

      do i=1,mfld
         call axhelm(bia(1,i),s(1,i),dfld(1,i),zeros,imsh,idir)
      enddo

      call binv(bia,jfld)

      return
      end
c-----------------------------------------------------------------------
      subroutine binvb(bib,s,jfld)

      ! apply local B^{-1} B to a scalar field

      ! s    := input scalar field
      ! jfld := field # of s
      ! bib  := B^{-1} B s

      include 'SIZE'
c     include 'MASS'
      include 'TSTEP' ! dep: nelfld

      parameter (lt=lx1*ly1*lz1*lelt)
      real bib(lt,1),s(lt,1)

      nel=nelfld(jfld)

      mfld=1
      if (jfld.eq.1) mfld=ldim

      do i=1,mfld
         call col3(bib(1,i),s(1,i),bm1,lx1*ly1*lz1*nel)
      enddo

      call binv(bib,jfld)

      return
      end
c-----------------------------------------------------------------------
      subroutine binvc(bic,ux,uy,uz,s,ifcf,ifuf,jfld)
      
      ! apply local B^{-1} C to a scalar field

      ! ux,uy,uz := input vector field
      ! s        := input scalar field
      ! ifcf     := true when ux,uy,uz is on the fine mesh
      ! ifuf     := true when s is on the fine mesh
      ! bic      := B^{-1} C(ux,uy,uz) s
      ! jfld     := field # of s

      include 'SIZE'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      real bic(lt,1),ux(lt),uy(lt),uz(lt),s(lt,1)
      logical ifcf,ifuf

      mfld=1
      if (jfld.eq.1) mfld=ldim

      do i=1,mfld
         call convect_new(bic(1,i),s(1,i),ifuf,ux,uy,uz,ifcf)
      enddo

      call binv(bic,jfld)

      return
      end
c-----------------------------------------------------------------------
      subroutine binv(fldi,ifld)

      ! compute fldi = B^-1 * fldi

      ! fldi: input field, can be a vector if ifld==1 
      ! ifld: field number associated with fldi

      include 'SIZE'
      include 'MASS'  ! dep: binvm1, bintm1
      include 'SOLN'  ! dep: tmask
      include 'TSTEP' ! dep: ifield, nelfld

      real fldi(lx1*ly1*lz1*lelt,1)

      jfield=ifield
      ifield=ifld

      nel=nelfld(ifield)
      n=lx1*ly1*lz1*nel

      if (ifield.eq.1) then
         call opmask(fldi(1,1),fldi(1,2),fldi(1,3))
         call opdssum(fldi(1,1),fldi(1,2),fldi(1,3))
      else
         call col2(fldi,tmask(1,1,1,1,ifield-1),n)
         call dssum(fldi,lx1,ly1,lz1)
      endif

      if (nel.eq.nelv) then
         call col2(fldi(1,1),binvm1,n)
         if (ifield.eq.1) then
            call col2(fldi(1,2),binvm1,n)
            if (ldim.eq.3) call col2(fldi(1,3),binvm1,n)
         endif
      else if (nel.eq.nelt) then
         call col2(fldi(1,1),bintm1,n)
      endif

      ifield=jfield

      return
      end
c-----------------------------------------------------------------------
      subroutine c_riesz_norm

      ! compute the norm of the riesz
      ! TODO: use ips to determine the norm

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real u_resid_h1_norm,u_resid_h10_norm,u_resid_l2_norm
      real t_resid_h1_norm,t_resid_h10_norm,t_resid_l2_norm
      real resid_h1_norm,resid_h10_norm,resid_l2_norm

      if (nio.eq.0) write (6,*) 'inside c_riesz_norm'

      if (lei.eq.0) then

         u_resid_h10_norm = 0
         u_resid_l2_norm = 0
         t_resid_h10_norm = 0
         t_resid_l2_norm = 0

         if (ifrom(1)) then
            ifield=1
            u_resid_h10_norm = h10vip(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
     $                                eh_u(1,1),eh_u(1,2),eh_u(1,ldim))
            u_resid_l2_norm = wl2vip(eh_u(1,1),eh_u(1,2),eh_u(1,ldim),
     $                               eh_u(1,1),eh_u(1,2),eh_u(1,ldim))
            u_resid_h1_norm = u_resid_h10_norm + u_resid_l2_norm
            if (nid.eq.0) write(6,*)'vel residual in h1 norm',
     $                               sqrt(u_resid_h1_norm)
            if (nid.eq.0) write(6,*)'vel residual in h10 norm',
     $                               sqrt(u_resid_h10_norm)
            if (nid.eq.0) write(6,*)'vel residual in l2 norm',
     $                               sqrt(u_resid_l2_norm)
         endif

         if (ifrom(2)) then
            ifield=2
            t_resid_h10_norm = h10sip(eh_t,eh_t)
            t_resid_l2_norm = wl2sip(eh_t,eh_t)

            t_resid_h1_norm = t_resid_h10_norm + t_resid_l2_norm
            if (nid.eq.0) write(6,*)'temp residual in h1 norm',
     $                               sqrt(t_resid_h1_norm)
            if (nid.eq.0) write(6,*)'temp residual in h10 norm',
     $                               sqrt(t_resid_h10_norm)
            if (nid.eq.0) write(6,*)'temp residual in l2 norm',
     $                               sqrt(t_resid_l2_norm)
         endif

         resid_h1_norm = sqrt(u_resid_h10_norm+u_resid_l2_norm+
     $                        t_resid_h10_norm+t_resid_l2_norm)
         resid_h10_norm = sqrt(u_resid_h10_norm+t_resid_h10_norm)
         resid_l2_norm = sqrt(u_resid_l2_norm+t_resid_l2_norm)

         if (resid_h1_norm.le.0) call
     $                exitti('negative semidefinite dual norm$',n)

         if (nid.eq.0) write (6,*) 'residual in h1 norm:'
     $                             ,resid_h1_norm
         if (nid.eq.0) write (6,*) 'residual in h10 norm:'
     $                             ,resid_h10_norm
         if (nid.eq.0) write (6,*) 'residual in l2 norm:'
     $                             ,resid_l2_norm
      else
         if (eqn.eq.'NS  ') then
            res_uu=0.
            res_tt=0.
            res=0.

            do j=1,nres_u
            do i=1,nres_u
               res_uu=res_uu
     $               +sigma_u(i+(nres_u)*(j-1))*theta_u(i)*theta_u(j)
            enddo
            enddo
            if (nid.eq.0)write(6,*)'velocity dual norm',sqrt(res_uu)
            if (ifrom(2)) then
               do j=1,nres_t
               do i=1,nres_t
                  res_tt=res_tt
     $                  +sigma_t(i+(nres_t)*(j-1))*theta_t(i)*theta_t(j)
               enddo
               enddo
               if (nid.eq.0)write(6,*)'temp dual norm',sqrt(res_tt)
            endif
            res=sqrt(res_uu+res_tt)
         else
            ! This block of code is to support when eqn .ne. NS
            res=0.

            do j=1,nres
            do i=1,nres
               res=res+mor_sigma(i,j)*mor_theta(i)*mor_theta(j)
            enddo
            enddo

            res=sqrt(res)
         endif
         if (res.le.0)
     $      call exitti('negative semidefinite dual norm$',nb)

         if (nid.eq.0) write (6,*) 'dual norm:',res
      endif

      call nekgsync
      if (nio.eq.0) write (6,*) 'exiting c_riesz_norm'
      return
      end
c-----------------------------------------------------------------------
      subroutine ainva(aia,s,dfld,jfld,idir)

      ! apply local A^{-1} A to a scalar field

      ! s    := input scalar field
      ! dfld := diffusivity field
      ! jfld := field # of s
      ! idir := direction of s
      ! aia  := A^{-1} A s

      include 'SIZE'
      include 'MASS'     ! dep: bm1
      include 'MOR'      ! dep: zeros
      include 'TSTEP'    ! dep: nelfld
      include 'PARALLEL' ! dep: nelg,nelgv

      parameter (lt=lx1*ly1*lz1*lelt)
      real aia(lt,1),s(lt,1),dfld(lt,1)

      nel=nelfld(jfld)
      melg=nelg(jfld)

      imsh=2
      if (melg.eq.nelgv) imsh=1 

      mfld=1
      if (jfld.eq.1) mfld=ldim

      do i=1,mfld
         call axhelm(aia(1,i),s(1,i),dfld(1,i),zeros,imsh,idir)
      enddo

      call ainv(aia,jfld)

      return
      end
c-----------------------------------------------------------------------
      subroutine ainvb(aib,s,jfld)

      ! apply local A^{-1} B to a scalar field

      ! s    := input scalar field
      ! jfld := field # of s
      ! aib  := A^{-1} B s

      include 'SIZE'
c     include 'MASS'
      include 'TSTEP' ! dep: nelfld

      parameter (lt=lx1*ly1*lz1*lelt)
      real aib(lt,1),s(lt,1)

      nel=nelfld(jfld)

      mfld=1
      if (jfld.eq.1) mfld=ldim

      do i=1,mfld
         call col3(aib(1,i),s(1,i),bm1,lx1*ly1*lz1*nel)
      enddo

      call ainv(aib,jfld)

      return
      end
c-----------------------------------------------------------------------
      subroutine ainvc(aic,ux,uy,uz,s,ifcf,ifuf,jfld)
      
      ! apply local A^{-1} C to a scalar field

      ! ux,uy,uz := input vector field
      ! s        := input scalar field
      ! ifcf     := true when ux,uy,uz is on the fine mesh
      ! ifuf     := true when s is on the fine mesh
      ! aic      := A^{-1} C(ux,uy,uz) s
      ! jfld     := field # of s

      include 'SIZE'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      real aic(lt,1),ux(lt),uy(lt),uz(lt),s(lt,1)
      logical ifcf,ifuf

      mfld=1
      if (jfld.eq.1) mfld=ldim

      do i=1,mfld
         call convect_new(aic(1,i),s(1,i),ifuf,ux,uy,uz,ifcf)
      enddo

      call ainv(aic,jfld)

      return
      end
c-----------------------------------------------------------------------
      subroutine ainv(fldi,ifld)

      ! compute fldi = A^-1 * fldi

      ! fldi: input field, can be a vector if ifld==1 
      ! ifld: field number associated with fldi

      include 'SIZE'
      include 'MASS'   ! dep: binvm1, bintm1
      include 'SOLN'   ! dep: v1mask, v2mask, v3mask, tmask
      include 'TSTEP'  ! dep: ifield, nelfld
      include 'INPUT'  ! dep: ifaxis, ifaziv
      include 'ORTHOT' ! dep: napproxt, name4t
      include 'VPROJ'  ! dep: vproj

      real fldi(lx1*ly1*lz1*lelt,1)
      real ifld1
      real h1(1),h2(1)
      real o1(lx1*ly1*lz1*lelt)
      real rie_tol, rie_it

      jfield=ifield
      ifield=ifld

      ifld1 = ifield-1

      nel=nelfld(ifield)
      n=lx1*ly1*lz1*nel

      if1=ifield-1
      write(name4t,1) if1-1
    1 format('PS',i2)
      if(ifield.eq.2) write(name4t,'(A4)') 'TEMP'

      isd = 1
      if (ifaxis.and.ifaziv.and.ifield.eq.2) isd = 2

      call rone (h1,n)
      call rzero(h2,n)

      ! set the tolerance = 1e-10 for all fields
      rie_tol = 1e-10
      ! set the max iter = 500
      rie_it = 500

      ! initial guess = 0
      if (ifield.eq.1) then
         imesh = 1
         call rzero(o1,n)
         call hsolve ('VELX',o1,fldi(1,1),h1,h2,v1mask,vmult
     $                      ,imesh,rie_tol,rie_it,1
     $                      ,vproj(1,1),0,binvm1)
         call copy(fldi(1,1),o1,n)

         call rzero(o1,n)
         call hsolve ('VELY',o1,fldi(1,2),h1,h2,v2mask,vmult
     $                      ,imesh,rie_tol,rie_it,2
     $                      ,vproj(1,2),0,binvm1)
         call copy(fldi(1,2),o1,n)

         if (if3d) then 
            call rzero(o1,n)
            call hsolve ('VELZ',o1,fldi(1,3),h1,h2,v3mask,vmult
     $                         ,imesh,rie_tol,rie_it,3
     $                         ,vproj(1,3),0,binvm1)
            call copy(fldi(1,3),o1,n)
         endif
      else
         call rzero(o1,n)
         if(nel.eq.nelv) then
            imesh = 1
            call hsolve  (name4t,o1,fldi,h1,h2
     $                    ,tmask(1,1,1,1,ifield-1)
     $                    ,tmult(1,1,1,1,ifield-1)
     $                    ,imesh,rie_tol,rie_it,1
     $                    ,approxt(1,0,ifld1),0,binvm1)
         else
            imesh = 2
            call hsolve  (name4t,o1,fldi,h1,h2
     $                    ,tmask(1,1,1,1,ifield-1)
     $                    ,tmult(1,1,1,1,ifield-1)
     $                    ,imesh,rie_tol,rie_it,1
     $                    ,approxt(1,0,ifld1),0,bintm1)
         endif
         call copy(fldi,o1,n)
      endif

      ifield=jfield

      return
      end
c-----------------------------------------------------------------------
