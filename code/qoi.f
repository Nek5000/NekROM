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
         call push_op(vx,vy,vz)

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

         call pop_op(vx,vy,vz)
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
      subroutine cnuss_setup

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /cnus1/ tbn(0:nb,0:nb),tbd(0:nb),tsa(0:nb)
      common /cnus2/ qwall(0:nb)
      common /scrk0/ tx(lt),ty(lt),tz(lt)

      if (inus.eq.1) then
         do j=0,nb
            do i=0,nb
               call ctbulk_num(tbn(i,j),ub(1,i),vb(1,i),wb(1,i),tb(1,j))
            enddo
            call ctbulk_den(tbd(j),ub(1,j),vb(1,j),wb(1,j))
            call ctsurf(tsa(j),tb(1,j))
         enddo

         do i=0,nb
            if (nio.eq.0) write (6,*) i,tsa(i),'tsa'
         enddo

         do j=0,nb
         do i=0,nb
            if (nio.eq.0) write (6,*) i,j,tbn(i,j),'tbn'
         enddo
         enddo

         do i=0,nb
            if (nio.eq.0) write (6,*) i,tbd(i),'tbd'
         enddo
      else if (inus.eq.2) then
         do i=0,nb
            call gradm1(tx,ty,tz,tb(1,i))

            eps=1.e-6
            ta=0.
            a=0.
            do ie=1,nelt
            do ifc=1,ldim*2
               call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,ifc)
               if (cbc(ifc,ie,2).eq.'t  ') then
                  x1=xm1(kx1,ky1,kz1,ie)
                  x2=xm1(kx2,ky2,kz2,ie)
                  z1=zm1(kx1,ky1,kz1,ie)
                  z2=zm1(kx2,ky2,kz2,ie)
                  xa=.5*(x1+x2)
                  za=.5*(z1+z2)
                  if (ifaxis.and.xa.gt.-eps) then
                     ta=ta+facint_v(tx,area,ifc,ie)
                     a=a+facint_v(ones,area,ifc,ie)
                  endif
                  if (ldim.eq.3.and.za.lt.eps) then
                     ta=ta-facint_v(tz,area,ifc,ie)
                     a=a+facint_v(ones,area,ifc,ie)
                  endif
               endif
            enddo
            enddo

            ta=glsum(ta,1)
            a=glsum(a,1)
            qwall(i)=ta/a
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ctsurf(tsurf,tt)

      include 'SIZE'
      include 'TOTAL'

      call savg(tsurf,a_surf,tt,1,'W  ')  ! tbar on wall

      return
      end
c-----------------------------------------------------------------------
      subroutine ctbulk_den(tbulk_den,uu,vv,ww)

      include 'SIZE'
      include 'TOTAL'

      common /nusvars/ diam

      parameter (lt=lx1*ly1*lz1*lelt)
      real uu(lt),vv(lt),ww(lt),tt(lt)

      integer icalld,idir
      save    icalld,idir
      data    icalld /0/
      data    idir /0/

      n=lx1*ly1*lz1*nelv

      if (icalld.eq.0) then
         ! determine stream-wise direction
         qu=glsc2(vx,bm1,n)
         qv=glsc2(vy,bm1,n)
         qw=glsc2(vz,bm1,n)
         idir=3
         if (abs(qu).gt.max(abs(qv),abs(qw))) then
            idir=1
         else if (abs(qv).gt.max(abs(qu),abs(qw))) then
            idir=2
         endif
         if (ldim.eq.2) then
            if (idir.eq.1) diam=glmax(ym1,n)-glmin(ym1,n)
            if (idir.eq.2) diam=glmax(xm1,n)-glmin(xm1,n)
         else
            if (idir.eq.1) then
               diam=glmax(ym1,n)+glmax(zm1,n)-glmin(ym1,n)-glmin(zm1,n)
            else if (idir.eq.2) then
               diam=glmax(xm1,n)+glmax(zm1,n)-glmin(xm1,n)-glmin(zm1,n)
            else
               diam=glmax(xm1,n)+glmax(ym1,n)-glmin(xm1,n)-glmin(ym1,n)
            endif
         endif
         if (nio.eq.0) write (6,*) qu,qv,qw,idir,diam,'bulk_vars'
         icalld=1
      endif

      if (idir.eq.1) then
         tbulk_den=glsc2(uu,bm1,n)
      else if (idir.eq.2) then
         tbulk_den=glsc2(vv,bm1,n)
      else
         tbulk_den=glsc2(ww,bm1,n)
      endif

    1 format (1p4e13.4,i4,'bulk_vars')

      return
      end
c-----------------------------------------------------------------------
      subroutine ctbulk_num(tbulk_num,uu,vv,ww,tt)

      include 'SIZE'
      include 'TOTAL'

      common /nusvars/ diam

      parameter (lt=lx1*ly1*lz1*lelt)
      real uu(lt),vv(lt),ww(lt),tt(lt)

      integer icalld,idir
      save    icalld,idir
      data    icalld /0/
      data    idir /0/

      n=lx1*ly1*lz1*nelv

      if (icalld.eq.0) then
         ! determine stream-wise direction
         qu=glsc2(vx,bm1,n)
         qv=glsc2(vy,bm1,n)
         qw=glsc2(vz,bm1,n)
         idir=3
         if (abs(qu).gt.max(abs(qv),abs(qw))) then
            idir=1
         else if (abs(qv).gt.max(abs(qu),abs(qw))) then
            idir=2
         endif
         if (ldim.eq.2) then
            if (idir.eq.1) diam=glmax(ym1,n)-glmin(ym1,n)
            if (idir.eq.2) diam=glmax(xm1,n)-glmin(xm1,n)
         else
            if (idir.eq.1) then
               diam=glmax(ym1,n)+glmax(zm1,n)-glmin(ym1,n)-glmin(zm1,n)
            else if (idir.eq.2) then
               diam=glmax(xm1,n)+glmax(zm1,n)-glmin(xm1,n)-glmin(zm1,n)
            else
               diam=glmax(xm1,n)+glmax(ym1,n)-glmin(xm1,n)-glmin(ym1,n)
            endif
         endif
         if (nio.eq.0) write (6,*) qu,qv,qw,idir,diam,'bulk_vars'
         icalld=1
      endif

      if (idir.eq.1) then
         tbulk_num=glsc3(uu,tt,bm1,n)
      else if (idir.eq.2) then
         tbulk_num=glsc3(vv,tt,bm1,n)
      else
         tbulk_num=glsc3(ww,tt,bm1,n)
      endif

    1 format (1p4e13.4,i4,'bulk_vars')

      return
      end
c-----------------------------------------------------------------------
      subroutine ctbulk(tbulk,uu,vv,ww,tt)

      include 'SIZE'
      include 'TOTAL'

      common /nusvars/ diam

      parameter (lt=lx1*ly1*lz1*lelt)
      real uu(lt),vv(lt),ww(lt),tt(lt)

      integer icalld,idir
      save    icalld,idir
      data    icalld /0/
      data    idir /0/

      n=lx1*ly1*lz1*nelv

      if (icalld.eq.0) then
         ! determine stream-wise direction
         qu=glsc2(vx,bm1,n)
         qv=glsc2(vy,bm1,n)
         qw=glsc2(vz,bm1,n)
         idir=3
         if (abs(qu).gt.max(abs(qv),abs(qw))) then
            idir=1
         else if (abs(qv).gt.max(abs(qu),abs(qw))) then
            idir=2
         endif
         if (ldim.eq.2) then
            if (idir.eq.1) diam=glmax(ym1,n)-glmin(ym1,n)
            if (idir.eq.2) diam=glmax(xm1,n)-glmin(xm1,n)
         else
            if (idir.eq.1) then
               diam=glmax(ym1,n)+glmax(zm1,n)-glmin(ym1,n)-glmin(zm1,n)
            else if (idir.eq.2) then
               diam=glmax(xm1,n)+glmax(zm1,n)-glmin(xm1,n)-glmin(zm1,n)
            else
               diam=glmax(xm1,n)+glmax(ym1,n)-glmin(xm1,n)-glmin(ym1,n)
            endif
         endif
         if (nio.eq.0) write (6,*) qu,qv,qw,idir,diam,'bulk_vars'
         icalld=1
      endif

      if (idir.eq.1) then
         tbulk=glsc3(uu,tt,bm1,n)/glsc2(uu,bm1,n)
      else if (idir.eq.2) then
         tbulk=glsc3(vv,tt,bm1,n)/glsc2(vv,bm1,n)
      else
         tbulk=glsc3(ww,tt,bm1,n)/glsc2(ww,bm1,n)
      endif

    1 format (1p4e13.4,i4,'bulk_vars')

      return
      end
c-----------------------------------------------------------------------
      subroutine cnuss

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      common /cnus1/ tbn(0:nb,0:nb),tbd(0:nb),tsa(0:nb)
      common /cnus2/ qwall(0:nb)
      common /nusvars/ diam

      parameter (lt=lx1*ly1*lz1*lelt)

      nv=nx1*ny1*nz1*nelv

      if (inus.eq.1) then
         qsurf=1.

         rhocp=param(7)
         cond=param(8)

         tbulk_num=0.
         tbulk_den=0.
         twall=0.

         do j=0,nb
            do i=0,nb
               tbulk_num=tbulk_num+tbn(i,j)*u(i,1)*ut(j,1)
            enddo
            tbulk_den=tbulk_den+tbd(j)*u(j,1)
            twall=twall+tsa(j)*ut(j,1)
         enddo

         tbulk=tbulk_num/tbulk_den
         h=(twall-tbulk)
         if (h.gt.0) h=qsurf/h
         rnus=diam*h/cond

         if (nio.eq.0) write (6,1) istep,time,twall,tbulk,rnus
      else if (inus.eq.2) then
         rnus=vlsc2(qwall,ut,nb+1)
         if (nio.eq.0) write (6,*) istep,time,rnus,'nus'
      endif

    1 format (i10,1p1e16.8,1p2e14.6,1p1e16.8,' fluxes')

      return
      end
c-----------------------------------------------------------------------
      subroutine savg(s_bar,a_surf,s,ifld,cb3)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real s(lx1,ly1,lz1,lelt)
      character*3 cb3

      integer e,f

      s_bar  = 0
      s_surf = 0
      a_surf = 0
      n_surf = 0

      do e=1,nelt
      do f=1,2*ndim
         if (cbc(f,e,ifld).eq.cb3) then
            call facint3(s1,a1,s,area,f,e) ! integrate a() on face f
            s_surf = s_surf + s1
            a_surf = a_surf + a1
            n_surf = n_surf + 1
         endif
      enddo
      enddo

      s_surf = glsum(s_surf,1)   ! sum across all processors
      a_surf = glsum(a_surf,1)
      n_surf = iglsum(n_surf,1)

      if (a_surf.gt.0) s_bar  = s_surf/a_surf

      return
      end
c-----------------------------------------------------------------------
      subroutine facint3(s1,a1,s,area,f,e) ! integrate a() on face f

      include 'SIZE'
      include 'TOPOL'

      real s(lx1,ly1,lz1,lelv),area(lx1,lz1,6,lelv)

      integer e,f

      call dsset(nx1,ny1,nz1) !     Set up counters
      iface  = eface1(f)
      js1    = skpdat(1,iface)
      jf1    = skpdat(2,iface)
      jskip1 = skpdat(3,iface)
      js2    = skpdat(4,iface)
      jf2    = skpdat(5,iface)
      jskip2 = skpdat(6,iface)

      i =0
      s1=0
      a1=0
      do 100 j2=js2,jf2,jskip2
      do 100 j1=js1,jf1,jskip1
         i  = i+1
         s1 = s1 + area(i,1,f,e)*s(j1,j2,1,e)
         a1 = a1 + area(i,1,f,e)
  100 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine neumann(a1,a2)

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      real a1(lx1,ly1,lz1,lelt),a2(lx1,ly1,lz1,lelt)

      do ie=1,nel
      do if=1,2*ldim
         if (cbc(if,ie,2).eq.'f  ') then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,if)
            l=1
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
               a2(ix,iy,iz,ie)=a2(ix,iy,iz,ie)*area(l,1,if,ie)
               l=l+1
            enddo
            enddo
            enddo
         endif
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
