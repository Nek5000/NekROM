#define INTP_NMAX 4000
#define LPART INTP_NMAX /* max number of particles per MPI rank */
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
c        write (6,*) 'ak pdrag',k,ak,bk,fd(1),fd(2),lk
      enddo
      call gop(fd(1),wk1,'+  ',1)
      call gop(fd(2),wk1,'+  ',1)

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

      common /scrk1/ ux(lt),uy(lt),uz(lt),w(lub+1)

      if (ifcdrag) then
         if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
            call read_serial(rdgx,nb+1,'qoi/rdgx ',w,nid)
            call read_serial(rdgy,nb+1,'qoi/rdgy ',w,nid)
            if (ldim.eq.3) call read_serial(rdgz,nb+1,'qoi/rdgz ',w,nid)
         else
            call push_op(vx,vy,vz)

            do i=0,nb
               call opcopy(vx,vy,vz,ub(1,i),vb(1,i),wb(1,i))
               nio = -1
               call torque_calc(1.,x0,.true.,.false.)
               rdgx(i)=dragvx(1)*ad_re
               rdgy(i)=dragvy(1)*ad_re
               rdgz(i)=dragvz(1)*ad_re
               nio = nid
               if (nio.eq.0) then
                  write (6,*) i,rdgx(i),rdgy(i),rdgz(i),ad_re,'dragi'
               endif
            enddo

            call pop_op(vx,vy,vz)
         endif
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
            vdx=vlsc2(rdgx,u,nb+1)/ad_re
            vdy=vlsc2(rdgy,u,nb+1)/ad_re
            pdx1=0.
            pdx2=0.
            pdx3=0.
            pdy1=0.
            pdy2=0.
            pdy3=0.
            call mxm(u,nb+1,ad_beta(2,count),3,tmp,1)
            do j=0,nb
               pdx1=pdx1+fd1(1+ldim*j)*tmp(j)
               pdy1=pdy1+fd1(2+ldim*j)*tmp(j)
               do i=0,nb
                  pdx2=pdx2+fd2(1+ldim*i+ldim*(nb+1)*j)*u(j)*u(i)
                  pdy2=pdy2+fd2(2+ldim*i+ldim*(nb+1)*j)*u(j)*u(i)
               enddo
               pdx3=pdx3+fd3(1+ldim*j)*u(j)
               pdy3=pdy3+fd3(2+ldim*j)*u(j)
            enddo
            pdx3=pdx3/ad_re
            pdy3=pdy3/ad_re
            dx=pdx1+pdx2+pdx3+vdx
            dy=pdy1+pdy2+pdy3+vdy
            if (nio.eq.0) write (6,2) time,vdx,pdx3,dx,'dragx'
            if (nio.eq.0) write (6,2) time,vdy,pdy3,dy,'dragy'
            if (ldim.eq.3) then
               dz=vlsc2(rdgz,u,nb+1)/ad_re
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

      common /cnus1/ tbn(0:lb,0:lb),tbd(0:lb),tsa(0:lb)
      common /cnus2/ qwall(0:lb)
      common /scrk0/ tx(lx1,ly1,lz1,lelt),
     $               ty(lx1,ly1,lz1,lelt),
     $               tz(lx1,ly1,lz1,lelt)

      if (inus.eq.1) then
         if (rmode.eq.'ON '.or.rmode.eq.'ONB') then
            call read_mat_serial(tbn,nb+1,nb+1,
     $         'qoi/tbn ',mb+1,mb+1,rtmp1,nid)
            call read_serial(tbd,nb+1,'qoi/tbd ',rtmp1,nid)
            call read_serial(tsa,nb+1,'qoi/tsa ',rtmp1,nid)
         else
            do j=0,nb
               do i=0,nb
                  call ctbulk_num(tbn(i+j*(nb+1),0),
     $               ub(1,i),vb(1,i),wb(1,i),tb(1,j,1))
               enddo
               call ctbulk_den(tbd(j),ub(1,j),vb(1,j),wb(1,j))
               call ctsurf(tsa(j),tb(1,j,1))
            enddo

            call dump_serial(tsa,nb+1,'qoi/tsa ',nid)

            do i=0,nb
               if (nio.eq.0) write (6,*) i,tsa(i),'tsa'
            enddo

            call dump_serial(tbn,(nb+1)**2,'qoi/tbn ',nid)

            do j=0,nb
            do i=0,nb
               if (nio.eq.0) write (6,*) i,j,tbn(i+j*(nb+1),0),'tbn'
            enddo
            enddo

            call dump_serial(tbd,nb+1,'qoi/tbd ',nid)

            do i=0,nb
               if (nio.eq.0) write (6,*) i,tbd(i),'tbd'
            enddo
         endif
      else if (inus.eq.2) then
         do i=0,nb
            call gradm1(tx,ty,tz,tb(1,i,1))

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
                  if (.not.ifaxis.and.xa.gt.(1.-eps)) then
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
      else if (inus.eq.3) then
         tbn(0,0)=2.
         do j=0,nb
            call ctsurf3(tsa(j),tb(1,j,1))
         enddo

         do i=0,nb
            if (nio.eq.0) write (6,*) i,tsa(i),'tsa'
         enddo
      else if (inus.eq.4) then
         call rone(ones,lx1*ly1*lz1*nelt)
         if (rmode.ne.'ON ') then
         do i=0,nb
            call gradm1(tx,ty,tz,tb(1,i,1))

            eps=1.e-3
            ta=0.
            a=0.
            do ie=1,nelt
            do ifc=1,ldim*2
               call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,ifc)
               if (cbc(ifc,ie,2).eq.'t  ') then
                  x1=xm1(kx1,ky1,kz1,ie)
                  x2=xm1(kx2,ky2,kz2,ie)
                  y1=ym1(kx1,ky1,kz1,ie)
                  y2=ym1(kx2,ky2,kz2,ie)
                  z1=zm1(kx1,ky1,kz1,ie)
                  z2=zm1(kx2,ky2,kz2,ie)
                  xa=.5*(x1+x2)
                  ya=.5*(y1+y2)
                  za=.5*(z1+z2)
                  if (ya.gt.(1.-eps)) then
                     ta=ta+facint_v(ty,area,ifc,ie)
                     a=a+facint_v(ones,area,ifc,ie)
                  endif
               endif
            enddo
            enddo

            ta=glsum(ta,1)
            a=glsum(a,1)
            qwall(i)=ta/a
         enddo
         call dump_serial(qwall,nb+1,'qoi/qwall ',nid)
         else
         call read_serial(qwall,nb+1,'qoi/qwall ',rtmp1,nid)
         endif
      else if (inus.eq.5) then
         call rone(ones,lx1*ly1*lz1*nelt)
         if (rmode.ne.'ON ') then
         do i=0,nb
            call gradm1(tx,ty,tz,tb(1,i,1))

            eps=1.e-3
            ta=0.
            a=0.
            do ie=1,nelt
            do ifc=1,ldim*2
               call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,ifc)
               if (cbc(ifc,ie,2).eq.'t  ') then
                  x1=xm1(kx1,ky1,kz1,ie)
                  x2=xm1(kx2,ky2,kz2,ie)
                  y1=ym1(kx1,ky1,kz1,ie)
                  y2=ym1(kx2,ky2,kz2,ie)
                  z1=zm1(kx1,ky1,kz1,ie)
                  z2=zm1(kx2,ky2,kz2,ie)
                  xa=.5*(x1+x2)
                  ya=.5*(y1+y2)
                  za=.5*(z1+z2)
                  if (xa.gt.(0.5-eps)) then
                     ta=ta+facint_v(tx,area,ifc,ie)
                     a=a+facint_v(ones,area,ifc,ie)
                  endif
               endif
            enddo
            enddo

            ta=glsum(ta,1)
            a=glsum(a,1)
            qwall(i)=ta/2
            if(nio.eq.0) write(6,*)qwall(i),'surface'
         enddo
         call dump_serial(qwall,nb+1,'qoi/qwall ',nid)
         else
         call read_serial(qwall,nb+1,'qoi/qwall ',rtmp1,nid)
         endif
      else if (inus.eq.6) then
         iobj=1

         do i=0,nb
            call gradm1(tx,ty,tz,tb(1,i,1))
            a=0.
            s=0.

            do imem=1,nmember(iobj)
               ieg=object(iobj,imem,1)
               ifc=object(iobj,imem,2)
               ie=gllel(ieg)
               call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,ifc)
               l=1

               do iz=kz1,kz2
               do iy=ky1,ky2
               do ix=kx1,kx2
                  a=a+area(l,1,ifc,ie)
                  s=s+area(l,1,ifc,ie)*(unx(l,1,ifc,ie)*tx(ix,iy,iz,ie)+
     $                                  uny(l,1,ifc,ie)*ty(ix,iy,iz,ie)+
     $                                  unz(l,1,ifc,ie)*tz(ix,iy,iz,ie))
                  l=l+1
               enddo
               enddo
               enddo
            enddo

            s=glsum(s,1)
            a=glsum(a,1)
            qwall(i)=s/a

            call dump_serial(qwall,nb+1,'qoi/qwall ',nid)
         enddo
      endif

c     if (iftflux) then
c        do j=0,nbavg
c        do i=0,nbavg
c           if (idirf.eq.1) then
c              call calc_tbulk(tbulkn(i,j),tbulkd(i,j),tb(1,i),ub(1,j))
c           else if (idirf.eq.2) then
c              call calc_tbulk(tbulkn(i,j),tbulkd(i,j),tb(1,i),vb(1,j))
c           else if (idirf.eq.3) then
c              call calc_tbulk(tbulkn(i,j),tbulkd(i,j),tb(1,i),wb(1,j))
c           endif
c        enddo
c        enddo

c        do i=0,nbavg
c           call calc_tsurf(tsurf(i),tb(1,i))
c           call calc_tmean(ttmean(i),tb(1,i))
c           call calc_sterm(st0(i),tb(1,i))
c           write (6,*) i,st0(i),'st0'
c        enddo
c     endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ctsurf(tsurf,tt)

      include 'SIZE'
      include 'TOTAL'

      common /nusselt/ fpmask(lx1,ly1,lz1,lelt),ffmask(2*ldim,lelt)

      call savg(tsurf,a_surf,tt,2,'f  ')  ! tbar on wall

      return
      end
c-----------------------------------------------------------------------
      subroutine ctsurf3(tsurf,s)

      include 'SIZE'
      include 'TOTAL'

      real s(lx1,ly1,lz1,lelt)

      tsurf=0.
      a=0.
      eps=1.e-3

      do ie=1,nelt
      do ifc=1,2*ldim
         call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,ifc)
         l=0
         ya=.5*(ym1(kx1,ky1,kz1,ie)+ym1(kx2,ky2,kz2,ie))
         if (ya.gt.(1.-eps).or.ya.lt.(-1.+eps)) then
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
               l=l+1
               x=xm1(ix,iy,iz,ie)
               if (x.gt.(-1.5-eps).and.x.lt.(2.5+eps)) then
                  a=a+area(l,1,ifc,ie)
                  tsurf=tsurf+area(l,1,ifc,ie)*s(ix,iy,iz,ie)
               endif
            enddo
            enddo
            enddo
         endif
      enddo
      enddo

      tsurf=glsum(tsurf,1)/glsum(a,1)

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

      common /cnus1/ tbn(0:lb,0:lb),tbd(0:lb),tsa(0:lb)
      common /cnus2/ qwall(0:lb)
      common /nusvars/ diam

      parameter (lt=lx1*ly1*lz1*lelt)

      nv=nx1*ny1*nz1*nelv

      if (inus.eq.1) then
         qsurf=1.
         diam=1.0

         rhocp=param(7)
         cond=param(8)

         tbulk_num=0.
         tbulk_den=0.
         twall=0.

         do j=0,nb
            do i=0,nb
               tbulk_num=tbulk_num+tbn(i,j)*u(i)*ut(j)
            enddo
            tbulk_den=tbulk_den+tbd(j)*u(j)
            twall=twall+tsa(j)*ut(j)
         enddo

         tbulk=tbulk_num/tbulk_den
         h=(twall-tbulk)
         if (h.gt.0) h=qsurf/h
         rnus=diam*h/cond

         if (nio.eq.0) write (6,1) ad_step,time,twall,tbulk,rnus
      else if (inus.eq.2) then
         rnus=vlsc2(qwall,ut,nb+1)
         if (nio.eq.0) write (6,*) ad_step,time,rnus,'nus'
      else if (inus.eq.3) then
         tbulk=tbn(0,0)
         twall=0.
         do j=0,nb
            twall=twall+tsa(j)*ut(j)
         enddo
         h=(twall-tbulk)
         diam=2.
         qsurf=1.
         if (h.gt.0) h=qsurf/h
         rnus=ad_pe*diam*h
         if (nio.eq.0) write (6,1) ad_step,time,twall,tbulk,rnus
      else if (inus.eq.4) then
         rnus=vlsc2(qwall,ut,nb+1)
         if (nio.eq.0) write (6,*) ad_step,time,rnus,'nus'
      else if (inus.eq.5) then
         rnus=vlsc2(qwall,ut,nb+1)
         if (nio.eq.0) write (6,*) ad_step,time,rnus,'nus'
      else if (inus.eq.6) then
         rnus=vlsc2(qwall,ut,nb+1)
         if (nio.eq.0) write (6,*) ad_step,time,rnus,'nus'
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
      subroutine cubar_setup

      include 'SIZE'
      include 'MASS'
      include 'MOR'

      common /morubar/ uu(0:lb),vv(0:lb),ww(0:lb)

      n=lx1*ly1*lz1*nelv

      do i=0,nb
         uu(i)=glsc2(ub(1,i),bm1,n)
         vv(i)=glsc2(vb(1,i),bm1,n)
         if (ldim.eq.3) ww(i)=glsc2(wb(1,i),bm1,n)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cubar

      include 'SIZE'
      include 'MASS'
      include 'TSTEP'
      include 'MOR'

      common /morubar/ uu(0:lb),vv(0:lb),ww(0:lb)

      n=lx1*ly1*lz1*nelv

      vol=glsum(bm1,n)
      ubar=vlsc2(uu,u,nb+1)/vol
      vbar=vlsc2(vv,v,nb+1)/vol
      if (ldim.eq.3) wbar=vlsc2(ww,w,nb+1)/vol

      if (nio.eq.0) write (6,*) ad_step,time,ubar,'ubar'

      return
      end
c-----------------------------------------------------------------------
      subroutine find_nelx(nelx)
      include 'SIZE'
      include 'TOTAL'
      integer e,f

      n=lx1*ly1*lz1*nelt

      ymin = glmin(ym1,n)
      ymax = glmax(ym1,n)

      dymin = ymax-ymin
      do e=1,nelt
         dy    = ym1(1,ny1,1,e)-ym1(1,1,1,e)
         dymin = min(dymin,dy)
      enddo
      dymin = glmin(dymin,1) ! Min across all processors

      nelx=0
      ytop = ymin+dymin
      do e=1,nelt
         ymid = ym1(1,2,1,e)
         if (ymid.lt.ytop) nelx = nelx+1
      enddo
      nelx = iglsum(nelx,1) ! Sum across all processors

      return
      end
c-----------------------------------------------------------------------
      subroutine find_nely(nely)
      include 'SIZE'
      include 'TOTAL'
      integer e,f

      n=lx1*ly1*lz1*nelt

      xmin = glmin(xm1,n)
      xmax = glmax(xm1,n)

      dxmin = xmax-xmin
      do e=1,nelt
         dx    = xm1(nx1,1,1,e)-xm1(1,1,1,e)
         dxmin = min(dxmin,dx)
      enddo
      dxmin = glmin(dxmin,1) ! Min across all processors

      nely=0
      xtop = xmin+dxmin
      do e=1,nelt
         xmid = xm1(2,1,1,e)
         if (xmid.lt.xtop) nely = nely+1
      enddo
      nely = iglsum(nely,1) ! Sum across all processors

      return
      end
c-----------------------------------------------------------------------
      subroutine gfldi(fi,f,x,y,z,n,mdim)

      ! generic field interpolation

      ! fi      := interpolated field
      ! f       := original field
      ! <x,y,z> := interpolation points
      ! n       := number of interpolation points
      ! mdim    := dimension of f

      include 'SIZE'
      include 'TOTAL'

      real fi(n,mdim),f(lx1,ly1,lz1,lelt,mdim),x(n),y(n),z(n)

      integer ih
      save    ih
      data    ih /0/

      common /rwk_intp/ rwk(INTP_NMAX,ldim+1)
      common /iwk_intp/ iwk(INTP_NMAX,3)

      if (n.gt.INTP_NMAX) call exitti ('n > INTP_NMAX in gfldi!$',n)

      if (ih.eq.0) call interp_setup(ih,0.0,0,nelt)

      call interp_nfld(fi,f,mdim,x,y,z,n,iwk,rwk,INTP_NMAX,.true.,ih)

      return
      end
c-----------------------------------------------------------------------
      subroutine interp_s(si,s,xyz,n)

      ! scalar field interpolation

      ! si   := interpolated field
      ! s    := original scalar field
      ! xyz  := interpolation points
      ! n    := number of interpolation points

      include 'SIZE'

      common /rwk_intp_pts/ x(INTP_NMAX),y(INTP_NMAX),z(INTP_NMAX)
      real si(n),s(lx1,ly1,lz1,lelt),xyz(ldim,n)

      call splitvec(x,y,z,xyz,ldim,n)
      call gfldi(si,s,x,y,z,n,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine interp_v(vi,v,xyz,n)

      ! vector field interpolation

      ! vi   := interpolated field
      ! v    := original vector field
      ! xyz  := interpolation points
      ! n    := number of interpolation points

      include 'SIZE'

      common /rwk_intp_pts/ x(INTP_NMAX),y(INTP_NMAX),z(INTP_NMAX)
      real vi(n,ldim),v(lx1,ly1,lz1,lelt,ldim),xyz(ldim,n)

      call splitvec(x,y,z,xyz,ldim,n)
      call gfldi(vi,v,x,y,z,n,ldim)

      return
      end
c-----------------------------------------------------------------------
      subroutine average_in_y ! Y averages yield functions of x
      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter(lt=lx1*ly1*lz1*lelt)
      common /myavg/ ubar(lt),vbar(lt),wbar(lt),tbar(lt)
      common /myavi/ igs_x,igs_y,igs_z,igs_xy

      parameter (lprof=201)
      common /mypts/ xyz(ldim,lprof),scalar(lprof)

      integer icalld,nelx,nely,nprofile
      save    icalld,nelx,nely,nprofile
      data    icalld /0/

      character*12 filename
      integer iunit

      n  = lx1*ly1*lz1*nelt
      nv = lx1*ly1*lz1*nelv

      if (icalld.eq.0) then
         icalld = 1

         call find_nely(nely)
         nelx = nelgt/nely
         call gtpp_gs_setup(igs_y,nelx,nely,1,2) ! Contract in (x & y)

         call domain_size(xmn,xmx,ymn,ymx,zmn,zmx)
         nprofile = lprof
         dx = (xmx-xmn)/(nprofile-1)
         ymid = (ymn+ymx)/2
         do i=1,nprofile
            xyz(1,i)=xmn + dx*(i-1) ! Pertub off of centerline (better for interpolation)
            xyz(2,i)=ymid + .0001
         enddo
         if (nid.gt.0) nprofile=0  !!! ONLY Node 0 initiates the interpolation!

      endif

        do i=0,nb
c       call planar_avg(tbar,tb(1,i) ,igs_y) ! Contract in x+y: work=f(z)
        call interp_t(scalar,xyz,nprofile,tbar)
        
        write(filename,'("q",I0,".dat")')i
        open(unit=iunit+100,file=filename,status='unknown')
        if (nid.eq.0) write(iunit+100,*) ' '
        if (nid.eq.0) write(iunit+100,20)(time,xyz(1,k),
     $  scalar(k),k=1,nprofile)
        close(unit=iunit+100)
        enddo
   20   format(1p3e15.6)  ! time, z , tbar


      return
      end

c-----------------------------------------------------------------------
      subroutine average_in_xy ! X-Y averages yield functions of Z
      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter(lt=lx1*ly1*lz1*lelt)
      common /myavg/ ubar(lt),vbar(lt),wbar(lt),tbar(lt)
      common /myavi/ igs_x,igs_y,igs_z,igs_xy

      parameter (lprof=201)
      common /mypts/ xyz(ldim,lprof),scalar(lprof)

      integer icalld,nelxy,nprofile
      save    icalld,nelxy,nprofile
      data    icalld /0/

      character*12 filename
      integer iunit

      n  = lx1*ly1*lz1*nelt
      nv = lx1*ly1*lz1*nelv

      if (icalld.eq.0) then
         icalld = 1

         call find_nelx(nelx)
         nely = nelgt/nelx
         call gtpp_gs_setup(igs_x,nelx,nely,1,1) ! Contract in (x & y)

         call domain_size(xmn,xmx,ymn,ymx,zmn,zmx)
         nprofile = lprof
         dy = (ymx-ymn)/(nprofile-1)
         xmid = (xmn+xmx)/2
         do i=1,nprofile
            xyz(1,i)=xmid + .0001  ! Pertub off of centerline (better for interpolation)
            xyz(2,i)=ymn + dy*(i-1)
         enddo
         if (nid.gt.0) nprofile=0  !!! ONLY Node 0 initiates the interpolation!
      endif

        do i=0,nb
c       call planar_avg(tbar,tb(1,i) ,igs_x) ! Contract in x+y: work=f(z)
        call interp_t(scalar,xyz,nprofile,tbar)
        
        write(filename,'("q",I0,".dat")')i
        open(unit=iunit+100,file=filename,status='unknown')
        if (nid.eq.0) write(iunit+100,*) ' '
        if (nid.eq.0) write(iunit+100,20)(time,xyz(2,k),
     $  scalar(k),k=1,nprofile)
        close(unit=iunit+100)
        enddo
   20   format(1p3e15.6)  ! time, z , tbar

      return
      end
c-----------------------------------------------------------------------
      subroutine calc_tbulk(tbulkn,tbulkd,tt,vv)

      include 'SIZE'
      include 'TOTAL'

      real tt(1)
      real vv(1)

      nv=lx1*ly1*lz1*nelv

      tbulkn=glsc3(tt,vv,bm1,nv)
      tbulkd=glsc2(vv,bm1,nv)

      return
      end
c-----------------------------------------------------------------------
      subroutine calc_tsurf(tsurf,tt)

      include 'SIZE'
      include 'TOTAL'

      real tt(lx1,ly1,lz1,lelt)

      nt=lx1*ly1*lz1*nelt

      nbavg=lb

      s=0.
      a=0.

      do ie=1,nelt
      do ifc=1,2*ldim
         call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,ifc)
         if (cbc(ifc,ie,2).eq.'f  ') then
            l=0
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
               l=l+1
               a=a+area(l,1,ifc,ie)
               s=s+area(l,1,ifc,ie)*tt(ix,iy,iz,ie)
            enddo
            enddo
            enddo
         endif
      enddo
      enddo

      s=glsum(s,1)
      a=glsum(a,1)

      tsurf=s/a

      return
      end
c-----------------------------------------------------------------------
      subroutine calc_sterm(s,tf)

      include 'SIZE'
      include 'TOTAL'

      real tf(lx1,ly1,lz1,lelt)

      a=0.
      s=0.

      do ie=1,nelt
      do ifc=1,2*ldim
         if (cbc(ifc,ie,2).eq.'f  ') then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,ifc)
            l=0
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
               l=l+1
               a=a+area(l,1,ifc,ie)
               s=s+area(l,1,ifc,ie)*tf(ix,iy,iz,ie)
            enddo
            enddo
            enddo
         endif
      enddo
      enddo

      s=glsum(s,1)
      a=glsum(a,1)

      return
      end
c-----------------------------------------------------------------------
