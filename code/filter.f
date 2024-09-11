c-----------------------------------------------------------------------
      subroutine pod_proj(uu,r1,nn,msg)

      ! Apply smoothing function to vector uu

      ! uu:= output/input, vector to be filtered
      ! r1:= input, number of modes to be filtered
      ! nn:= input, length of vector uu
      ! msg:= input, used as an indicator for different filter
      ! functions. Currently, we support msg = step, linear, parabo, and
      ! cubic.

      real uu(nn)
      real a1,a2,a3,a4
      integer r1,nn,ncut
      character*6  msg

      if (msg.eq.'step  ') then
         call rzero(uu(nn-r1+1),r1)
      elseif (msg.eq.'linear') then
         ncut = nn-r1
         do i=ncut+1,nn
            uu(i) = (-uu(ncut)/r1)*(i-nn)
         enddo
      elseif (msg.eq.'parabo') then
         ncut = nn-r1
         do i=ncut+1,nn
            uu(i) = (-uu(ncut)/r1**2)*(i-ncut)**2 + uu(ncut)
         enddo
      elseif (msg.eq.'cubic ') then
         ncut = nn-r1
         a1 = 2
         a2 = -3*(ncut+nn)
         a3 = 6*ncut*nn
         a4 = nn**3-3*nn**2*ncut
         do i=ncut+1,nn
            uu(i) = (a1*i**3+a2*i**2+a3*i+a4)*uu(ncut)/(ncut+nn)**3
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine pod_df(uu)

      ! Apply differential filter (DF) to vector uu

      ! uu := vector that will be filtered

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      integer icalld
      save    icalld
      data    icalld /0/

      real tmp(nb)
      real uu(nb)

      if (icalld.eq.0) then
         call setdf(dfops,au,bu,rdft*rdft)
         call dgetrf(nb,nb,dfops,nb,ipiv,info)
         icalld=1
      endif

      if (rdft.le.1e-12) then
      else
         ! setup right-handed side
         call mxm(bu,nb,uu,nb,tmp,1)
         ! solve for uu
         call dgetrs('N',nb,1,dfops,nb,ipiv,tmp,nb,info)
         call copy(uu,tmp,nb)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setdf(flu,a,b,ad_diff)

      ! Set up differential filter (DF) operator B + ad_diff * A

      ! flu := differential filter operator
      ! a := ROM stiffness matrix
      ! b := ROM mass matrix
      ! ad_diff := radius of the DF

      include 'SIZE'
      include 'MOR'

      real flu(nb*nb),a(nb*nb),b(nb*nb)
      real ad_diff
      
      call cmult2(flu,a,ad_diff,nb*nb)
      call add2(flu,b,nb*nb)
         
      return
      end
c-----------------------------------------------------------------------
      subroutine set_les_imp(fles1,fles2)

      ! set implicit les matrices

      ! fles1 := filtering matrix for velocity
      ! fles2 := filtering matrix for temperature

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /mytmp/ dfld(lt),t1(lt)

      real fles1(0:nb,nb),fles2(0:nb,nb)

      nv=lx1*ly1*lz1*nelv
      nt=lx1*ly1*lz1*nelv

      call push_sol(vx,vy,vz,pr,t)

      do i=1,nb
         if (ifrom(1)) then
            call opcopy(vx,vy,vz,ub(1,i),vb(1,i),wb(1,i))
            call opadd2(vx,vy,vz,ub,vb,wb)
         endif

         if (ifrom(2)) then
            call copy(t,tb(1,i,1),nt)
            call add2(t,tb,nt)
         endif

         call q_filter(1.)
         if (ifrom(1)) call pv2b(fles1(0,i),vx,vy,vz,ub,vb,wb)
         if (ifrom(2)) call ps2b(fles2(0,i),t,tb)
      enddo

      do i=1,nb
         if (ifrom(1)) fles1(i,i)=fles1(i,i)-1.
         if (ifrom(2)) fles2(i,i)=fles2(i,i)-1.
      enddo

      call pop_sol(vx,vy,vz,pr,t)

      return
      end
c-----------------------------------------------------------------------
      subroutine apply_les_imp(uu,tt,sig,fles1,fles2,tmp)

      ! apply implicit les filter to coefficients

      ! uu    := velocity coefficients
      ! tt    := temperature coefficients
      ! sig   := filter strength
      ! fles1 := filter matrix for velocity
      ! fles2 := filter matrix for temperature
      ! tmp   := work array

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /mytmp/ dfld(lt),t1(lt)

      real uu(0:nb),tt(0:nb),fles1(0:nb,nb),fles2(0:nb,nb),tmp(0:nb)

      do i=1,nb
         if (ifrom(1)) then
            call mxm(fles1,nb+1,uu(1),nb,tmp,1)
            tmp(0)=0.
            call add2s1(tmp,uu,-sig,nb+1)
         endif
         if (ifrom(2)) then
            call mxm(fles2,nb+1,tt(1),nb,tmp,1)
            tmp(0)=0.
            call add2s1(tmp,tt,-sig,nb+1)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine evalnut(nut,u1,u2,u3)

      ! evaluate eddy viscosity

      ! nut        := eddy viscosity field
      ! <u1,u2,u3> := velocity field

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      real nut(lx1,ly1,lz1,lelt),u1(lt),u2(lt),u3(lt)

      common /scrns/ du(lx1,ly1,lz1,lelt),
     $   t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)
      common /scruz/ t7(lt),t8(lt),t9(lt)
      common /cfldx/ dri(lx1),dsi(ly1),dti(lz1)
c
      integer icalld
      save    icalld
      data    icalld /0/

      if (icalld.eq.0) then
         icalld=1
         call getdr(dri,zgm1(1,1),lx1)
         call getdr(dsi,zgm1(1,2),ly1)
         if (if3d) call getdr(dti,zgm1(1,3),lz1)
      endif

      n=lx1*ly1*lz1*nelv

      call gradm1(t1,t2,t3,u1)
      call gradm1(t4,t5,t6,u2)
      if (ldim.eq.3) call gradm1(t7,t8,t9,u3)

      do i=1,n
         if (ldim.eq.2) then
            du(i,1,1,1)=sqrt(t1(i)*t1(i)+t2(i)*t2(i)+t3(i)*t3(i)
     $                +t4(i)*t4(i)+t5(i)*t5(i)+t6(i)*t6(i))
         else
            du(i,1,1,1)=sqrt(t1(i)*t1(i)+t2(i)*t2(i)+t3(i)*t3(i)
     $                +t4(i)*t4(i)+t5(i)*t5(i)+t6(i)*t6(i)
     $                +t7(i)*t7(i)+t8(i)*t8(i)+t9(i)*t9(i))
         endif
      enddo

      sq2=sqrt(2.)
      cs=.2

      do ie=1,nelt
         l=0
         do iz=1,lz1
         do iy=1,ly1
         do ix=1,lx1
            l=l+1
            hk=min(min(dri(ix),dsi(iy)),dti(iz))*jacmi(l,ie)
            nut(ix,iy,iz,ie)=sq2*du(ix,iy,iz,ie)*(cs*hk)**2
         enddo
         enddo
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
