c-----------------------------------------------------------------------
      subroutine setcnv_c(cx,cy,cz)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ltd=lxd*lyd*lzd*lelt)

      common /convect/ c1v(ltd),c2v(ltd),c3v(ltd),
     $                 u1v(ltd),u2v(ltd),u3v(ltd)

      real cx(lt),cy(lt),cz(lt)

      call set_convect_new(c1v,c2v,c3v,cx,cy,cz)

      return
      end
c-----------------------------------------------------------------------
      subroutine setcnv_u(ux,uy,uz)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ltd=lxd*lyd*lzd*lelt)

      common /convect/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      real ux(lt), uy(lt), uz(lt)

      call intp_rstd_all(u1v,ux,nelv)
      call intp_rstd_all(u2v,uy,nelv)
      if (ldim.eq.3) call intp_rstd_all(u3v,uz,nelv)

      return
      end
c-----------------------------------------------------------------------
      subroutine cct(ct,i,j) ! compute C(u) * t set by setcnv

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'MOR'

      parameter(lt=lx1*ly1*lz1*lelt)

      real ct(lt)

      if (ifaxis) then
         n=lx1*ly1*lz1*nelv
         call push_op(vx,vy,vz)
         call opcopy(vx,vy,vz,ub(1,i),vb(1,i),wb(1,i))
         call conv1d(ct,tb(1,j))
         call col2(ct,bm1,n)
         call pop_op(vx,vy,vz)
      else
         call convect_new(ct,u1v(1,j),.true.,
     $                    c1v(1,i),c2v(1,i),c3v(1,i),.true.)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ccu(cu1,cu2,cu3,i,j) ! compute C(u) * u set by setcnv

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MOR'
      include 'MASS'

      parameter(lt=lx1*ly1*lz1*lelt)

      common /scrwk/ wk1(lt),wk2(lt),wk3(lt),wk4(lt),wk5(lt),wk6(lt)

      real cu1(lt),cu2(lt),cu3(lt)

      if (ifaxis) then
         n=lx1*ly1*lz1*nelv
         call push_op(vx,vy,vz)
         call opcopy(vx,vy,vz,ub(1,i),vb(1,i),wb(1,i))
         call conv1d(cu1,ub(1,j))
         call conv1d(cu2,vb(1,j))
         call col2(cu1,bm1,n)
         call col2(cu2,bm1,n)
         call pop_op(vx,vy,vz)
      else
         call convect_new(cu1,u1v(1,j),.true.,
     $                    c1v(1,i),c2v(1,i),c3v(1,i),.true.)
         call convect_new(cu2,u2v(1,j),.true.,
     $                    c1v(1,i),c2v(1,i),c3v(1,i),.true.)
         if (ldim.eq.3) call convect_new(cu3,u3v(1,j),.true.,
     $                                c1v(1,i),c2v(1,i),c3v(1,i),.true.)
      endif

      if (ifcdrag.and.ifield.eq.1.and..not.ifread) then
         call opcopy(wk4,wk5,wk6,cu1,cu2,cu3)
         call opbinv1(wk1,wk2,wk3,wk4,wk5,wk6,1.)
         call cint(fd2(1,j,i),wk1,wk2,wk3)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ccu_new(cu1,cu2,cu3,i,j)

      include 'SIZE'
      include 'MOR'

      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(ltd=lxd*lyd*lzd*lelt)

      common /scrccunew/ t1(lt),t2(lt),t3(lt)
      common /scrwk/ wk1(lt),wk2(lt),wk3(lt),wk4(lt),wk5(lt),wk6(lt)

      real cu1(lt),cu2(lt),cu3(lt)

      call convect_new(cu1,u1v(1,j),.true.,
     $                 c1v(1,i),c2v(1,i),c3v(1,i),.true.)
      call convect_new(cu2,u2v(1,j),.true.,
     $                 c1v(1,i),c2v(1,i),c3v(1,i),.true.)
      if (ldim.eq.3) call convect_new(cu3,u3v(1,j),.true.,
     $                                c1v(1,i),c2v(1,i),c3v(1,i),.true.)

      if (ifcdrag.and.ifield.eq.1) then
         call opcopy(wk4,wk5,wk6,cu1,cu2,cu3)
         call opbinv1(wk1,wk2,wk3,wk4,wk5,wk6,1.)
         call cint(fd2(1,i,j),wk1,wk2,wk3)
      endif

      if (i.ne.j) then
         call convect_new(t1,u1v(1,i),.true.,
     $                    c1v(1,j),c2v(1,j),c3v(1,j),.true.)
         call convect_new(t2,u2v(1,i),.true.,
     $                    c1v(1,j),c2v(1,j),c3v(1,j),.true.)
         if (ldim.eq.3) call convect_new(t3,u3v(1,i),.true.,
     $                                c1v(1,j),c2v(1,j),c3v(1,j),.true.)
         call opadd2(cu1,cu2,cu3,t1,t2,t3)
         if (ifcdrag.and.ifield.eq.1) then
            call opcopy(wk4,wk5,wk6,t1,t2,t3)
            call opbinv1(wk1,wk2,wk3,wk4,wk5,wk6,1.)
            call cint(fd2(1,j,j),wk1,wk2,wk3)
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine intp_rstd_all(uf,u,nel)

      include 'SIZE'
      include 'INPUT'

      parameter (lxyz1=lx1*ly1*lz1)
      parameter (lxyzd=lxd*lyd*lzd)

      real uf(lxyzd,lelt), u(lxyz1,lelt)

      do i=1,nel
         call intp_rstd(uf(1,i),u(1,i),lx1,lxd,if3d,0) ! 0 --> forward
      enddo

      return
      end
c-----------------------------------------------------------------------
