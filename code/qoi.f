c-----------------------------------------------------------------------
      subroutine comp_pdrag(fd,px,py,pz)

      ! fd: drag vector
      ! px,py,pz: pressure gradient

      ! computes the pressure drag on object 1 given the pressure gradient

      include 'SIZE'
      include 'SOLN'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /pdrag/ wsck(lt,2,lk),abveck(2,ldim,lk),wk1(lt)

      real px(lt),py(lt),pz(lt)
      real fd(ldim)

      call setup_pdrag(px,py,pz)

      do k=1,lk
         write (6,*) abveck(1,1,k),abveck(1,2,k),'avec'
      enddo

      do k=1,lk
         write (6,*) abveck(2,1,k),abveck(2,2,k),'bvec'
      enddo

      n=lx1*ly1*lz1*nelt

      if (ldim.eq.3) 
     $   call exitti('ldim.eq.3 not supported in comp_pdrag$',n)

      fd(1)=0.
      fd(2)=0.

      one=1.
      twopi=8.*atan(one)

      do k=1,lk
         ak=glsc2(wsck(1,2,k),g,n) ! cosine coefficients
         bk=glsc2(wsck(1,1,k),g,n) ! sine coefficients

         fd(1)=fd(1)+(ak*abveck(1,1,k)+bk*abveck(2,1,k))/(twopi*k)
         fd(2)=fd(2)+(ak*abveck(1,2,k)+bk*abveck(2,2,k))/(twopi*k)
         write (6,*) 'ak pdrag',k,ak,bk,fd(1),fd(2)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setup_pdrag(px,py,pz)

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
               if (ifinit) call setup_pdrag_init(ifc,gllel(ieg))
               call setup_pdrag_helper(px,py,pz,ifc,gllel(ieg))
            endif
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setup_pdrag_init(if,ie)

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
      subroutine setup_pdrag_helper(px,py,pz,if,ie)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real px(lx1,ly1,lz1,lelt),
     $     py(lx1,ly1,lz1,lelt),
     $     pz(lx1,ly1,lz1,lelt)

      call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,if)
c     write (6,*) 'ks',kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,if

c     call outpost(px,py,pz,pr,t,'ppp')

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
