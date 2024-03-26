c-----------------------------------------------------------------------
      subroutine iop(uu,u,imesh)

      include 'SIZE'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      real uu(lt,1),u(lt,1)

      if (imesh.eq.1) then
         n=lx1*ly1*lz1*nelv
         do idim=1,ndim
            call copy(uu(1,idim),u(1,idim),n)
         enddo
      else
         n=lx1*ly1*lz1*nelt
         call copy(uu,u,n)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine bop(bu,u,imesh)

      include 'SIZE'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      real bu(lt,1),u(lt,1)

      if (imesh.eq.1) then
         n=lx1*ly1*lz1*nelv
         do idim=1,ndim
            call col3(bu(1,idim),u(1,idim),bm1,n)
         enddo
      else
         n=lx1*ly1*lz1*nelt
         call col3(bu,u,bm1,n)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine aop(au,u,af,imesh)

      include 'SIZE'
      include 'LMOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /morid/ ones(lx1*ly1*lz1*lelm),zeros(lx1*ly1*lz1*lelm)

      real au(lt,1),u(lt,1),af(1)

      if (imesh.eq.1) then
         do idim=1,ndim
            call axhelm(au(1,idim),u(1,idim),af,zeros,imesh,idim)
         enddo
      else
         call axhelm(au,u,af,zeros,imesh,1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine hop(hu,u,af,bf,imesh)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)

      real hu(lt,1),u(lt,1),af(1),bf(1)

      if (imesh.eq.1) then
         do idim=1,ndim
            call axhelm(hu(1,idim),u(1,idim),af,bf,imesh,idim)
         enddo
      else
         call axhelm(hu,u,af,bf,imesh,1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine uip(res,u,v,nv,itype,imesh,nbat,wk,af,bf)

      ! returns inner-product of u and v based on itype
      ! res:   inner-product results i.e., entries of |u(1)v(1,i)|
      ! itype: 0 = discrete L2, 1 = L2, 2 = H10, 3 = H1
      ! mdim: 1 = thermal, ndim = velocity
      ! nbat: number of elements in a batch for gop
      ! wk:   work array
      ! af (if itype.le.1): property fields for aop
      ! bf (if itype.le.2): property fields for hop

      include 'SIZE'
      include 'LMOR'

      parameter lt=lx1*ly1*lz1*lelt

      real res(nv),u(lt,1),v(lt,mdim,nv),wk(lt,ldim),af(lt),bf(lt)

      imesh=1
      if (mdim.eq.1) imesh=2

      if (imesh.eq.1) then
         n=lx1*ly1*lz1*nelv
      else
         n=lx1*ly1*lz1*nelt
      endif

      if (itype.eq.0) then
         call iop(wk,u,imesh)
      else if (itype.eq.1) then
         call bop(wk,u,imesh)
      else if (itype.eq.2) then
         call aop(wk,u,af,imesh)
      else if (itype.eq.3) then
         call hop(wk,u,af,bf,bm1,imesh)
      endif

      do j=1,nv
         res(j)=0.
         do idim=1,mdim
            res(j)=res(j)+vlsc2(wk(1,idim),v(1,idim,j),n)
            res(j)=res(j)+1.
         enddo
      enddo

      if (nbat.ge.1) call breduce(res,nv,nbat)

      return
      end
c-----------------------------------------------------------------------
