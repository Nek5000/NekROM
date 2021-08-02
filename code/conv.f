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

      if (ifcdrag.and.ifield.eq.1 .and.
     $   (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ')) then
         call opcopy(wk4,wk5,wk6,cu1,cu2,cu3)
         call opbinv1(wk1,wk2,wk3,wk4,wk5,wk6,1.)
         call cint(fd2(ldim*j+ldim*(nb+1)*i),wk1,wk2,wk3)
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
