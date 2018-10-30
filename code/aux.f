c-----------------------------------------------------------------------
      subroutine factor3(mq,mp,mr,m)

      integer dmin,d

      n=m
      l=nint(real(n)**(1/3))

      dmin=n
      imin=-1

      do i=1,n
          d=abs(n-i**3)
          if (d.lt.dmin.and.mod(n,i).eq.0) then
              dmin=d
              imin=i
          endif
      enddo

      mp=imin
      n=n/mp

      dmin=n
      imin=-1

      do i=1,n
          d=abs(n-i*i)
          if (d.lt.dmin.and.mod(n,i).eq.0) then
              dmin=d
              imin=i
          endif
      enddo

      mq=imin
      mr=n/mq

      if (nio.eq.0) write (6,*) 'mp,mq,mr,mp*mq*mr',mp,mq,mr,mp*mq*mr

      return
      end
c-----------------------------------------------------------------------
      subroutine setpart(mps,mp,n)

      integer mps(mp)

      do i=0,mp-1
         mps(i+1)=n/mp+max(mod(n,mp)-i,0)/max(mod(n,mp)-i,1)
         mps(i+1)=mps(i+1)+mps(max(i,1))*max(i,0)/max(i,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setpart3(mps,mqs,mrs,mp,mq,mr,nb)

      integer mps(mp),mqs(mq),mrs(mr)

      call setpart(mps,mp,nb)
      call setpart(mqs,mq,nb+1)
      call setpart(mrs,mr,nb+1)

      return
      end
c-----------------------------------------------------------------------
      function i2p(i,mps)

      integer mps(1)

      i2p=0

      do while (mps(i2p+1).ne.0)
         i2p=i2p+1
         if (i.le.mps(i2p)) return
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ijk2pqr(ip,iq,ir,i,j,k,mps,mqs,mrs)

      real mps(1),mqs(1),mrs(1)

      ip=i2p(i,mps)
      iq=i2p(j+1,mqs)
      ir=i2p(k+1,mrs)

      return
      end
c-----------------------------------------------------------------------
      function ijk2pid(i,j,k,mps,mqs,mrs,mp,mq,mr)

      real mps(1),mqs(1),mrs(1)

      call ijk2pqr(ip,iq,ir,i,j,k,mps,mqs,mrs)
      ijk2pid=(ip-1)+mp*(iq-1)+mp*mq*(ir-1)

      return
      end
c-----------------------------------------------------------------------
      subroutine setrange(mps,mqs,mrs,mp,mq,mr)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      integer mps(1),mqs(1),mrs(1)

      ip=mod(nid,mp)
      iq=mod(nid/mp,mq)
      irr=   (nid/mp)/mq

      i0=mps(max(ip,1))*max(ip,0)/max(ip,1)+1
      i1=mps(ip+1)

      j0=mqs(max(iq,1))*max(iq,0)/max(iq,1)
      j1=mqs(iq+1)-1

      k0=mrs(max(irr,1))*max(irr,0)/max(irr,1)
      k1=mrs(irr+1)-1

      return
      end
c-----------------------------------------------------------------------
      subroutine ijk2l(l,i,j,k)
      
      include 'SIZE'
      include 'MOR'

      il=i-i0
      jl=j-j0
      kl=k-k0

      l=il+jl*(i1-i0+1)+kl*(i1-i0+1)*(j1-j0+1)+1

      return
      end
c-----------------------------------------------------------------------
      subroutine opadd3 (a1,a2,a3,b1,b2,b3,c1,c2,c3)

      include 'SIZE'

      real a1(1),a2(1),a3(1),b1(1),b2(1),b3(1)
      real c1(1),c2(1),c3(1)

      ntot1=lx1*ly1*lz1*nelv
      call add2(a1,b1,c1,ntot1)
      call add2(a2,b2,c2,ntot1)
      if (ldim.eq.3) call add2(a3,b3,c3,ntot1)

      return
      end
c-----------------------------------------------------------------------
      subroutine reconstruct(ux,uy,uz)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ux(lt),uy(lt),uz(lt)

      n=lx1*ly1*lz1*nelv

      call opzero(ux,uy,uz)

      do i=0,nb
         call opadds(ux,uy,uz,ub(1,i),vb(1,i),wb(1,i),u(i,1),n,2)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine fom_analysis

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrns/ t1(lt),t2(lt),t3(lt)
      common /ctracker/ cmax(0:nb), cmin(0:nb)

      if (istep.eq.0) then
         call rom_init

         call gengram
         call genevec
         call genbases
         call genops

         do i=0,nb
            cmax(i) = -1e10
            cmin(i) =  1e10
         enddo

         call rom_setup
      endif

      u(0,1) = 1.

      call opsub3(t1,t2,t3,vx,vy,vz,ub(1,0),vb(1,0),wb(1,0))

      if (ifl2) then
         call wl2proj(u(1,1),t1,t2,t3)
      else
         call h10proj(u(1,1),t1,t2,t3)
      endif

      do i=0,nb
         if (u(i,1).lt.cmin(i)) cmin(i)=u(i,1)
         if (u(i,1).gt.cmax(i)) cmax(i)=u(i,1)
      enddo

      if (mod(istep,100).eq.0) then
         n=istep/100

         if (nio.eq.0) then
            do i=0,nb
               write (6,*) n,time,u(i,1),'cmax'
               write (6,*) n,time,u(i,1),'coef'
               write (6,*) n,time,u(i,1),'cmin'
            enddo
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
