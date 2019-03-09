c-----------------------------------------------------------------------
      function csig(g1,g2,g3,h1,h2)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrsig/ r1(lt),r2(lt),r3(lt)
      real g1(1),g2(1),g3(1),h1(1),h2(1)

      call ophinv(r1,r2,r3,g1,g2,g3,h1,h2,tolh,nmxhi)

      csig=h10prod(r1,r2,r3,r1,r2,r3,h1,h2)

      return
      end
c-----------------------------------------------------------------------
      subroutine sets() ! set sigmas

      include 'SIZE'
      include 'MOR'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrsets/ g1(lt),g2(lt),g3(lt),h1(lt),h2(lt)

      n=lx1*ly1*lz1*nelt

      call rone(h1,n)
      call rzero(h2,n)

      do j=0,nb
         call col3(g1,ub(1,j),bm1,n)
         call col3(g2,vb(1,j),bm1,n)
         if (ldim.eq.3) call col3(g3,wb(1,j),bm1,n)
         sigb(j)=csig(g1,g2,g3,h1,h2)

         call axhelm(g1,ub(1,j),h1,h2,1,1)
         call axhelm(g2,vb(1,j),h1,h2,1,1)
         if (ldim.eq.3) call axhelm(g3,wb(1,j),h1,h2,1,1)
         siga(j)=csig(g1,g2,g3,h1,h2)

         call setcnv_c(ub(1,j),vb(1,j),wb(1,j))
         do i=0,nb
            call setcnv_u(ub(1,i),vb(1,i),wb(1,i))
            call ccu(g1,g2,g3)
            sigc(i,j)=csig(g1,g2,g3,h1,h2)
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine sett() ! set thetas

      include 'SIZE'
      include 'MOR'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrsets/ g1(lt),g2(lt),g3(lt),h1(lt),h2(lt)

      esum=0.

      n=lx1*ly1*lz1*nelt

      call rone(h1,n)
      call rzero(h2,n)

      rinstep=1./real(ad_nsteps)

      call rzero(thb,nb+1)
      call add2s2(thb,uj(0,1),-5./6,nb+1)
      call add2s2(thb,uj(0,2), 1./6,nb+1)
      call add2s2(thb,uj(0,3),-1./3,nb+1)
      call add2s2(thb,uj(0,4), 1./3,nb+1)
      call add2s2(thb,uj(0,5),-7./6,nb+1)
      call add2s2(thb,uj(0,6),11./6,nb+1)
      s=ad_dt*rinstep
      call cmult(thb,s,nb+1)

      s=1./ad_re
      call cmult2(tha,ua,s,nb+1)

      call copy(thc,u2a,(nb+1)**2)

      call add2s2(thc,u2j(0,0,2),-1.*rinstep,(nb+1)**2)
      call add2s2(thc,u2j(0,0,5),-1.*rinstep,(nb+1)**2)
      call add2s2(thc,u2j(0,0,6), 2.*rinstep,(nb+1)**2)

      return
      end
c-----------------------------------------------------------------------
      function eest() ! compute error estimate

      include 'SIZE'
      include 'MOR'

      eest=0.

      do j=0,nb
         eest=eest+sigb(j)*thb(j)**2
         eest=eest+siga(j)*tha(j)**2
         do i=0,nb
            eest=eest+sigc(i,j)*thc(i,j)**2
         enddo
      enddo

      eest=sqrt(eest)

      return
      end
c-----------------------------------------------------------------------
