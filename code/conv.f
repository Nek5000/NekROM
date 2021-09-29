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

      common /convect/ c1v(ltd),c2v(ltd),c3v(ltd),
     $                 u1v(ltd),u2v(ltd),u3v(ltd)

      real ux(lt),uy(lt),uz(lt)

      call intp_rstd_all(u1v,ux,nelv)
      call intp_rstd_all(u2v,uy,nelv)
      if (ldim.eq.3) call intp_rstd_all(u3v,uz,nelv)

      return
      end
c-----------------------------------------------------------------------
      subroutine setcnv_u1(u)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ltd=lxd*lyd*lzd*lelt)

      common /convect/ c1v(ltd),c2v(ltd),c3v(ltd),
     $                 u1v(ltd),u2v(ltd),u3v(ltd)

      real u(lt)

      call intp_rstd_all(u1v,u,nelv)

      return
      end
c-----------------------------------------------------------------------
      subroutine cc(ct,mdim) ! compute C(u) * t set by setcnv

      ! compute c * u set by setcnv_c and set_cnv_u

      ! ct   := output convected field
      ! mdim := dimension of ct and u field

      include 'SIZE' ! dep: lt,ltd

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ltd=lxd*lyd*lzd*lelt)

      common /convect/ c1v(ltd),c2v(ltd),c3v(ltd),uv(ltd,3)

      real ct(lt,mdim)

      do idim=1,mdim
         call convect_new(
     $      ct(1,idim),uv(1,idim),.true.,c1v,c2v,c3v,.true.)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine convect_axis(ct,mdim,ux,uy,uz,s)

      ! convect vector/scalar field in axisymmetric formulation

      ! ct         := advection field
      ! mdim       := dimension of input field s and output field ct
      ! <ux,uy,uz> := advecting velocity field
      ! s          := advected field

      include 'SIZE' ! dep: lt
      include 'SOLN' ! dep: vx,vy,vz
      include 'MASS' ! dep: bm1
      include 'MOR'  ! dep: ub,vb,wb,tb

      parameter (lt=lx1*ly1*lz1*lelt)

      real ct(lt,mdim),ux(lt),uy(lt),uz(lt),s(lt,mdim)

      n=lx1*ly1*lz1*nelv

      call push_op(vx,vy,vz)

      call opcopy(vx,vy,vz,ux,uy,uz)

      do idim=1,mdim
         call conv1d(ct(1,idim),s(1,idim))
         call col2(ct(1,idim),bm1,n)
      enddo

      call pop_op(vx,vy,vz)

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
      subroutine cpart(ic1,ic2,jc1,jc2,kc1,kc2,nloc,nb,np,ip)

      call cpart_helper(npi,npj,npk,np)

      l=1

      do i=1,npi
         if (i.eq.1) then
            ic=0
         else
            ic=ic+ni
         endif
         ni=(nb+1)/npi+max(min(i-npi+(nb+1)-((nb+1)/npi)*npi,1),0)
         do j=1,npj
            if (j.eq.1) then
               jc=0
            else
               jc=jc+nj
            endif
            nj=(nb+1)/npj+max(min(j-npj+(nb+1)-((nb+1)/npj)*npj,1),0)
            do k=1,npk
               if (k.eq.1) then
                  kc=1
               else
                  kc=kc+nk
               endif
               nk=nb/npk+max(min(k-npk+nb-(nb/npk)*npk,1),0)
               if (ip.eq.l) then
                  ic1=ic
                  ic2=ic+ni-1
                  jc1=jc
                  jc2=jc+nj-1
                  kc1=kc
                  kc2=kc+nk-1
                  nloc=(ic2-ic1+1)*(jc2-jc1+1)*(kc2-kc1+1)
               endif
               l=l+1
            enddo
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cpart_helper(ni,nj,nk,np)

      ntmp=np

      nk=ceiling((1.*ntmp)**(1./3.))

      do while ((ntmp/nk)*nk.ne.ntmp)
         nk=nk+1
      enddo

      ntmp=ntmp/nk

      nj=ceiling((1.*ntmp)**(1./2.))

      do while ((ntmp/nj)*nj.ne.ntmp)
         nj=nj+1
      enddo

      ntmp=ntmp/nj

      ni=ceiling((1.*ntmp)**(1./1.))

      if (ni*nj*nk.ne.np) call exitti('ni*nj*nk != np$',ni*nj*nk)

      return
      end
c-----------------------------------------------------------------------
      subroutine setc_local(cl,cel,ic1,ic2,jc1,jc2,kc1,kc2,ic,jc,kc)

      real cl(ic1:ic2,jc1:jc2,kc1:kc2)

      if (ic.ge.ic1.and.ic.le.ic2 .and.
     $    jc.ge.jc1.and.jc.le.jc2 .and.
     $    kc.ge.kc1.and.kc.le.kc2) cl(ic,jc,kc)=cel

      return
      end
c-----------------------------------------------------------------------
