c-----------------------------------------------------------------------
      subroutine cint(fd,px,py,pz)

      ! fd: integral
      ! px,py,pz: gradient field

      ! computes the contour integral on object 1 given the gradient field

      include 'SIZE'
      include 'SOLN'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /pdrag/ wsck(lt,2,lk),abveck(2,ldim,lk),wk1(lt)

      real px(lt),py(lt),pz(lt)
      real fd(ldim)

c     call setup_pdrag(px,py,pz)
      call cint_helper(px,py,pz)

      n=lx1*ly1*lz1*nelt

      if (ldim.eq.3) 
     $   call exitti('ldim.eq.3 not supported in cpdrag$',n)

      fd(1)=0.
      fd(2)=0.

      one=1.
      twopi=8.*atan(one)

      do k=1,lk
         ak=glsc2(wsck(1,2,k),g,n) ! cosine coefficients
         bk=glsc2(wsck(1,1,k),g,n) ! sine coefficients

         fd(1)=fd(1)-(ak*abveck(1,1,k)+bk*abveck(2,1,k))/real(2*k)
         fd(2)=fd(2)+(ak*abveck(1,2,k)+bk*abveck(2,2,k))/real(2*k)
c        write (6,*) 'ak pdrag',k,ak,bk,fd(1),fd(2)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cint_helper(px,py,pz)

      include 'SIZE'
      include 'TOTAL'
      include 'LMOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /pdrag/ wsck(lt,2,lk),abveck(2,ldim,lk)

      integer icalld
      save    icalld
      data    icalld /0/

      real px(lt),py(lt),pz(lt)
      logical ifinit

      n=lx1*ly1*lz1*nelt
      twopi=2.*pi

      ifinit=.false.

      if (icalld.eq.0) then
         ifinit=.true.
         call rzero(abveck,2*ldim*lk)
         icalld=1
      endif

      do iobj=1,nobj
         memtot=nmember(iobj)
         do mem=1,memtot
            ieg=object(iobj,mem,1)
            ifc=object(iobj,mem,2)
            if (gllnid(ieg).eq.nid) then
               if (ifinit) call cint_helper_init(ifc,gllel(ieg))
               call cint_cdgds(px,py,pz,ifc,gllel(ieg))
            endif
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cint_helper_init(if,ie)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      common /pdrag/ wsck(lx1,ly1,lz1,lelt,2,lk),abveck(2,ldim,lk)

      call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,if)

      twopi=2.*pi

      if (ldim.eq.3) call exitti('ldim.eq.3 not supported (pdrag)$',lt)

      ia=1
      do iz=kz1,kz2
      do iy=ky1,ky2
      do ix=kx1,kx2
         stwopi=sobj(ix,iy,iz,ie)*twopi
         aa=area(ia,1,if,ie)
         ux=unx(ia,1,if,ie)
         uy=uny(ia,1,if,ie)
         do k=1,lk
            wsk=sin(stwopi*k)*aa
            wck=cos(stwopi*k)*aa

            abveck(1,1,k)=abveck(1,1,k)+ux*wsk/real(k)
            abveck(1,2,k)=abveck(1,2,k)+uy*wsk

            abveck(2,1,k)=abveck(2,1,k)+ux*wck
            abveck(2,2,k)=abveck(2,2,k)+uy*wck

            wsck(ix,iy,iz,ie,1,k)=2.*wsk/real(k*pi)
            wsck(ix,iy,iz,ie,2,k)=2.*wck/real(k*pi)
         enddo
         ia=ia+1
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cint_cdgds(px,py,pz,if,ie)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real px(lx1,ly1,lz1,lelt),
     $     py(lx1,ly1,lz1,lelt),
     $     pz(lx1,ly1,lz1,lelt)

      call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,if)

      twopi=2.*pi

      if (ldim.eq.3) call exitti('ldim.eq.3 not supported (pdrag)$',lt)

      ia=1
      do iz=kz1,kz2
      do iy=ky1,ky2
      do ix=kx1,kx2
         tx=t1x(ia,1,if,ie)
         ty=t1y(ia,1,if,ie)
         tt=sqrt(tx*tx+ty*ty)
         g(ix,iy,iz,ie)=-tx*px(ix,iy,iz,ie)/tt
     $                 -ty*py(ix,iy,iz,ie)/tt
         ia=ia+1
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cvdrag_setup

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk1/ ux(lt),uy(lt),uz(lt)

      if (ifcdrag) then
         call opcopy(ux,uy,uz,vx,vy,vz)

         do i=0,nb
            nio = -1
            call opcopy(vx,vy,vz,ub(1,i),vb(1,i),wb(1,i))
            call torque_calc(1.,x0,.true.,.false.)
            rdgx(i)=dragvx(1)
            rdgy(i)=dragvy(1)
            rdgz(i)=dragvz(1)
            nio = nid
            if (nio.eq.0) then
               write (6,*) i,rdgx(i),rdgy(i),rdgz(i),'dragi'
            endif
         enddo

         call opcopy(vx,vy,vz,ux,uy,uz)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cdrag

      include 'SIZE'
      include 'TSTEP'
      include 'MOR'

      real tmp(0:nb)

      if (ifcdrag) then
         if (nio.eq.0) then
            vdx=vlsc2(rdgx,u,nb+1)
            vdy=vlsc2(rdgy,u,nb+1)
            pdx1=0.
            pdx2=0.
            pdx3=0.
            pdy1=0.
            pdy2=0.
            pdy3=0.
            call mxm(u,nb+1,ad_beta(2,count),3,tmp,1)
            do j=0,nb
               pdx1=pdx1+fd1(1,j)*tmp(j)
               pdy1=pdy1+fd1(2,j)*tmp(j)
               do i=0,nb
                  pdx2=pdx2+fd2(1,i,j)*u(j,1)*u(i,1)
                  pdy2=pdy2+fd2(2,i,j)*u(j,1)*u(i,1)
               enddo
               pdx3=pdx3+fd3(1,j)*u(j,1)
               pdy3=pdy3+fd3(2,j)*u(j,1)
            enddo
            dx=pdx1+pdx2+pdx3+vdx
            dy=pdy1+pdy2+pdy3+vdy
c           if (nio.eq.0) write (6,1) time,vdx,pdx1,pdx2,pdx3,dx,'dragx'
c           if (nio.eq.0) write (6,1) time,vdy,pdy1,pdy2,pdy3,dy,'dragy'
            if (nio.eq.0) write (6,2) time,vdx,pdx3,dx,'dragx'
            if (nio.eq.0) write (6,2) time,vdy,pdy3,dy,'dragy'
            if (ldim.eq.3) then
               dz=vlsc2(rdgz,u,nb+1)
               write (6,*) ad_step*dt,vdz,'dragz'
            endif
         endif
      endif

    1 format (1p6e16.8,2x,a)
    2 format (1p4e16.8,2x,a)

      return
      end
c-----------------------------------------------------------------------
