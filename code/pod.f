c-----------------------------------------------------------------------
      subroutine setbases

      ! set ub,vb,wb,tb in /morbasis/ with pod basis

      include 'SIZE'
      include 'MOR'

      if (nio.eq.0) write (6,*) 'inside setbases'

      call nekgsync
      bas_time=dnekclock()

      if (rmode.ne.'ON ') ifrecon=.true.

      if (rmode.eq.'ONB'.or.rmode.eq.'CP ') then
         call loadbases
      else if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ') then
         if (ifrom(1)) then
            call pod(uvwb(1,1,1),evec(1,1,1),eval(1,1),ug(1,1,1),
     $         us0,ldim,ips,nb,ns,ifpb)
            do ib=1,nb
               call opcopy(ub(1,ib),vb(1,ib),wb(1,ib),
     $            uvwb(1,1,ib),uvwb(1,2,ib),uvwb(1,ldim,ib))
            enddo
            if (.not.ifcomb.and.ifpb) call vnorm(ub,vb,wb)
         else
            call opcopy(ub,vb,wb,uic,vic,wic)
         endif
         if (ifrom(2)) then
            call pod(tb(1,1),evec(1,1,2),eval(1,2),ug(1,1,2),
     $         ts0,1,ips,nb,ns,ifpb)
            if (.not.ifcomb.and.ifpb) call snorm(tb)
         endif

         if (ifcomb.and.ifpb) call cnorm(ub,vb,wb,tb)
      endif

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ') then
         call dump_gram
         call dump_bas
      endif

      call nekgsync
      if (nio.eq.0) write (6,*) 'bas_time:',dnekclock()-bas_time
      if (nio.eq.0) write (6,*) 'exiting setbases'

      return
      end
c-----------------------------------------------------------------------
      subroutine ps2k(ck,ux,uub)

      ! set snapshot coefficients for a given scalar basis

      ! ck  := coefficients
      ! ux  := snapshots
      ! uub := basis functions

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ck(0:nb,ls),ux(lt,ls),uub(lt,0:nb)

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ') then
         do i=1,ns
            if (nio.eq.0) write (6,*) 'ps2k: ',i,'/',ns
            nio=-1
            call ps2b(ck(0,i),ux(1,i),uub)
            nio=nid
         enddo
      else
         ! implement read here
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine pv2k(ck,usnap,uub,vvb,wwb)

      ! set snapshot coefficients for a given vector basis

      ! ck           := coefficients
      ! usnap        := snapshots
      ! uub,vvb,wwb, := x,y,z components of basis functions

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrgg/ uu(lt),vv(lt),ww(lt)

      real ck(0:nb,ls),usnap(lt,ldim,ls),
     $     uub(lt,0:nb),vvb(lt,0:nb),wwb(lt,0:nb)

      n=lx1*ly1*lz1*nelt

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ') then
c        if (ips.eq.'H10') then
c           do j=1,ns
c              call axhelm(uu,usnap(1,1,j),ones,zeros,1,1)
c              call axhelm(vv,usnap(1,2,j),ones,zeros,1,2)
c              if (ldim.eq.3)
c    $            call axhelm(ww,usnap(1,ldim,j),ones,zeros,1,3)
c              ck(0,j)=1.
c              do i=1,nb
c                 ck(i,j)=glsc2(uu,uub(1,i),n)+glsc2(vv,vvb(1,i),n)
c                 if (ldim.eq.3) ck(i,j)=ck(i,j)+glsc2(ww,wwb(1,i),n)
c              enddo
c           enddo
c        else
         do i=1,ns
            if (nio.eq.0) write (6,*) 'pv2k: ',i,'/',ns
            nio=-1
            call pv2b(ck(0,i),usnap(1,1,i),usnap(1,2,i),usnap(1,ldim,i),
     $           uub,vvb,wwb)
            nio=nid
         enddo
c        endif
      else
         ! implement read here
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ps2b(coef,tt,sb)

      ! get coordinates of a scalar field for a given basis

      ! ck  := coordinates of <ux> in <uub>
      ! ux  := FOM scalar field
      ! uub := basis functions

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real coef(0:nb),tt(lt),sb(lt,0:nb)

      if (nio.eq.0) write (6,*) 'inside ps2b'

      n=lx1*ly1*lz1*nelt

      coef(0) = 1.
      if (nio.eq.0) write (6,1) coef(0),coef(0),1.

      do i=1,nb
         ww=sip(sb(1,i),sb(1,i))
         vv=sip(sb(1,i),tt)
         coef(i) = vv/ww
         if (nio.eq.0) write (6,1) coef(i),vv,ww,ips
      enddo

      if (nio.eq.0) write (6,*) 'exiting ps2b'

    1 format(' coef',1p3e16.8,1x,a3)

      return
      end
c-----------------------------------------------------------------------
      subroutine pv2b(coef,ux,uy,uz,uub,vvb,wwb)

      ! get coordinates of a vector field for a given basis

      ! coef        := coordinates of <ux,uy,uz> in <uub,vvb,wwb>
      ! ux,uy,uz    := x,y,z components of FOM field
      ! uub,vvb,wwb := x,y,z components of basis functions

      include 'SIZE'
      include 'MOR'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ux(lt),uy(lt),uz(lt),uub(lt,0:nb),vvb(lt,0:nb),wwb(lt,0:nb)

      real coef(0:nb)

      if (nio.eq.0) write (6,*) 'inside pv2b'

      n=lx1*ly1*lz1*nelt

      coef(0) = 1.
      if (nio.eq.0) write (6,1) 0,coef(0),coef(0),1.

      jfield=ifield
      ifield=1

      do i=1,nb
         ww=vip(uub(1,i),vvb(1,i),wwb(1,i),uub(1,i),vvb(1,i),wwb(1,i))
         vv=vip(uub(1,i),vvb(1,i),wwb(1,i),ux,uy,uz)
         coef(i) = vv/ww
         if (nio.eq.0) write (6,1) i,coef(i),vv,ww,ips
      enddo

      ifield=jfield

      if (nio.eq.0) write (6,*) 'exiting pv2b'

    1 format('coef',i8,1p3e16.8,1x,a3)

      return
      end
c-----------------------------------------------------------------------
      subroutine pc2b(cfu,cft,ux,uy,uz,tt,uub,vvb,wwb,ttb)

      ! get coordinates of a combined field for a given basis

      ! cfu & cft       := coords of <ux,uy,uz,tt> in <uub,vvb,wwb,ttb>
      ! ux,uy,uz,tt     := vx,vy,vz,t components of combined FOM field
      ! uub,vvb,wwb,ttb := vx,vy,vz,t components of basis functions

      include 'SIZE'
      include 'MOR'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ux(lt),uy(lt),uz(lt),tt(lt)
      real uub(lt,0:nb),vvb(lt,0:nb),wwb(lt,0:nb),ttb(lt,0:nb)
      real cfu(0:nb),cft(0:nb)

      if (nio.eq.0) write (6,*) 'inside pc2b'

      n=lx1*ly1*lz1*nelt

      cfu(0) = 1.
      cft(0) = 1.
      if (nio.eq.0) write (6,1) 0,cfu(0),cfu(0),1.

      jfield=ifield
      ifield=1

      do i=1,nb
         ww=cip(uub(1,i),vvb(1,i),wwb(1,i),ttb(1,i),
     $          uub(1,i),vvb(1,i),wwb(1,i),ttb(1,i))
         vv=cip(uub(1,i),vvb(1,i),wwb(1,i),ttb(1,i),ux,uy,uz,tt)
         cfu(i) = vv/ww
         cft(i) = cfu(i)
         if (nio.eq.0) write (6,1) i,cfu(i),vv,ww,ips
      enddo

      ifield=jfield

      if (nio.eq.0) write (6,*) 'exiting pc2b'

    1 format('coef',i8,1p3e16.8,1x,a3)

      return
      end
c-----------------------------------------------------------------------
      function sip(t1,t2)

      ! return inner-product of scalar fields

      ! t1,t2 := scalar fields
      ! ips   := inner product type
      !          (L2 = L_2, H10 = H^1_0, HLM = Helmholtz)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt)

      if (ips.eq.'L2 ') then
         sip=wl2sip(t1,t2)
      else if (ips.eq.'H10') then
         sip=h10sip(t1,t2)
      else if (ips.eq.'HLM') then
         sip=hlmsip(t1,t2)
      else
         if (nid.eq.0) write (6,*) 'ips: ',ips
         call exitti('did not provide supported inner product space$',1)
      endif

      return
      end
c-----------------------------------------------------------------------
      function vip(t1,t2,t3,t4,t5,t6)

      ! return inner-product of vector fields

      ! t1,t2,t3 := x,y,z components of field 1
      ! t4,t5,t6 := x,y,z components of field 2
      ! ips   := inner product type
      !          (L2 = L_2, H10 = H^1_0, HLM = Helmholtz)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)

      if (ips.eq.'L2 ') then
         vip=wl2vip(t1,t2,t3,t4,t5,t6)
      else if (ips.eq.'H10') then
         vip=h10vip(t1,t2,t3,t4,t5,t6)
      else if (ips.eq.'HLM') then
         vip=hlmvip(t1,t2,t3,t4,t5,t6)
      else
         if (nid.eq.0) write (6,*) 'ips: ',ips
         call exitti('did not provide supported inner product space$',1)
      endif

      return
      end
c-----------------------------------------------------------------------
      function cip(t1,t2,t3,t4,t5,t6,t7,t8)

      ! return inner-product of vector fields

      ! t1,t2,t3,t4 := vx,vy,vz,t components of field 1
      ! t5,t6,t7,t8 := vx,vy,vz,t components of field 2

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt),t7(lt),t8(lt)

      cip=podrat*vip(t1,t2,t3,t5,t6,t7)+(1.-podrat)*sip(t4,t8)

      return
      end
c-----------------------------------------------------------------------
      function h10sip_vd(t1,t2,vd)

      ! return inner-product of scalar fields using the H^1_0
      ! inner-product with an arbitrary diffusivity field

      ! t1,t2 := scalar fields
      ! vd := variable diffusivity

      include 'SIZE'
      include 'SOLN'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),vd(lt)

      common /scrip/ t3(lt)

      call axhelm(t3,t1,vd,zeros,1,1)
      h10sip_vd=glsc2(t3,t2,lx1*ly1*lz1*nelt)

      return
      end
c-----------------------------------------------------------------------
      function h10sip(t1,t2)

      ! return inner-product of scalar fields using the H^1_0
      ! inner-product

      ! t1,t2 := scalar fields

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt)

      common /scrip/ t3(lt)

      call axhelm(t3,t1,ones,zeros,1,1)
      h10sip=glsc2(t3,t2,lx1*ly1*lz1*nelt)

      return
      end
c-----------------------------------------------------------------------
      function h10vip(t1,t2,t3,t4,t5,t6)

      ! return inner-product of vector fields using the H^1_0
      ! inner-product

      ! t1,t2,t3; t4,t5,t6 := x,y,z components of vector fields

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)

      common /scrip/ t7(lt),t8(lt),t9(lt)

      n=lx1*ly1*lz1*nelt

      call axhelm(t7,t1,ones,zeros,1,1)
      h10vip=glsc2(t7,t4,n)

      call axhelm(t8,t2,ones,zeros,1,2)
      h10vip=h10vip+glsc2(t8,t5,n)

      if (ldim.eq.3) then
         call axhelm(t9,t3,ones,zeros,1,3)
         h10vip=h10vip+glsc2(t9,t6,n)
      endif

      return
      end
c-----------------------------------------------------------------------
      function h10vip_vd(t1,t2,t3,t4,t5,t6,vd)

      ! return inner-product of vector fields using the H^1_0
      ! inner-product with arbitrary diffusivity fields

      ! t1,t2,t3; t4,t5,t6 := x,y,z components of vector fields
      ! vd := variable diffusivity in multiple dimensions

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt),vd(lt,ldim)

      common /scrip/ t7(lt),t8(lt),t9(lt)

      n=lx1*ly1*lz1*nelt

      call axhelm(t7,t1,vd(1,1),zeros,1,1)
      h10vip_vd=glsc2(t7,t4,n)

      call axhelm(t8,t2,vd(1,2),zeros,1,2)
      h10vip_vd=h10vip_vd+glsc2(t8,t5,n)

      if (ldim.eq.3) then
         call axhelm(t9,t3,vd(1,3),zeros,1,3)
         h10vip_vd=h10vip_vd+glsc2(t9,t6,n)
      endif

      return
      end
c-----------------------------------------------------------------------
      function hlmsip(t1,t2)

      ! return inner-product of scalar fields using the Helmholtz
      ! inner-product

      ! t1,t2 := scalar fields

      include 'SIZE'
      include 'MASS'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt)

      hlmsip=h10sip(t1,t2)/ad_re+wl2sip(t1,t2)*ad_beta(1,3)/ad_dt

      return
      end
c-----------------------------------------------------------------------
      function hlmvip(t1,t2,t3,t4,t5,t6)

      ! return inner-product of vector fields using the Helmholtz
      ! inner-product

      ! t1,t2,t3; t4,t5,t6 := x,y,z components of vector fields

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)

      hlmvip=h10vip(t1,t2,t3,t4,t5,t6)/ad_re
     $      +wl2vip(t1,t2,t3,t4,t5,t6)*ad_beta(1,3)/ad_dt

      return
      end
c-----------------------------------------------------------------------
      function wl2sip(t1,t2)

      ! return inner-product of scalar fields using the L^2
      ! inner-product

      ! t1,t2 := scalar fields

      include 'SIZE'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt)

      n=lx1*ly1*lz1*nelt

      wl2sip = glsc3(t1,t2,bm1,n)

      return
      end
c-----------------------------------------------------------------------
      function wl2vip(t1,t2,t3,t4,t5,t6)

      ! return inner-product of vector fields using the L^2
      ! inner-product

      ! t1,t2,t3; t4,t5,t6 := x,y,z components of vector fields

      include 'SIZE'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)

      wl2vip=op_glsc2_wt(t1,t2,t3,t4,t5,t6,bm1)

      return
      end
c-----------------------------------------------------------------------
      subroutine setgram

      ! set the Gramian according to ifpod flags

      include 'SIZE'
      include 'MOR'
      include 'SOLN'

      call nekgsync
      sg_start=dnekclock()

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ') then
         jfield=ifield
         ifield=1
         if (ifpod(1)) call gengram(ug(1,1,1),us0,ns,ldim,ips)
         ifield=2
         if (ifpod(2)) call gengram(ug(1,1,2),ts0,ns,1,ips)
         ifield=jfield
      endif

      if (ifcomb) then
         podrati=1.-podrat
         call add3s2(ug(1,1,1),ug(1,1,1),ug(1,1,2),podrat,podrati,ns*ns)
         call copy(ug(1,1,2),ug(1,1,1),ns*ns)
      endif

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ')
     $   call dump_gram

      call nekgsync
      if (nio.eq.0) write (6,*) 'gram_time:',dnekclock()-sg_start

      return
      end
c-----------------------------------------------------------------------
      subroutine setevec

      ! set the eigenvectors based on the generated Gramians

      include 'SIZE'
      include 'MOR'

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ') then
         if (nb > 1) then 
c           do i=0,ldimt1
c              if (ifpod(i)) then 
c                 call genevec(evec(1,1,i),eval(1,i),ug(1,1,i),i)
c                 call cnmax(eval(1,i),fname,i)
c              endif
c           enddo
            if (ifpod(1)) then 
               call genevec(evec(1,1,1),eval(1,1),ug(1,1,1),ns,nb,1)
               call cnmax(eval(1,1),'ops/enlu ',1)
               call cenpm(eval(1,1),'ops/enru ',1)
            endif
            if (ifpod(2)) then 
               call genevec(evec(1,1,2),eval(1,2),ug(1,1,2),ns,nb,2)
               call cnmax(eval(1,2),'ops/enlt ',2)
            endif
         else
            if(nid.eq.0) write(6,*)'nb = 1, no POD'
            do i=0,ldimt1
               if (ifpod(i)) then 
                  call rzero(evec(1,1,i),ls*lb)
                  evec(1,1,i) = 1
               endif
            enddo
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine hlmgg(gram,s,ms,mdim)

      ! set the Gramian based on the Helmholtz inner-product

      ! gram := Gramian
      ! s    := snapshots
      ! ms   := number of snapshots
      ! mdim := vector dimension

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrgg/ uu(lt),vv(lt),ww(lt)

      real gram(ms,ms)
      real s(lt,mdim,ms)

      if (nio.eq.0) write (6,*) 'inside gengram HLM'

      n=lx1*ly1*lz1*nelt

      s2=ad_beta(1,3)/ad_dt

      if (ifield.eq.1) then
         s1=1./ad_re
      else
         s1=1./ad_pe
      endif

      do j=1,ms
         call axhelm(uu,s(1,1,j),ones,zeros,1,1)
         if (mdim.ge.2) call axhelm(vv,s(1,2,j),ones,zeros,1,2)
         if (mdim.eq.3) call axhelm(ww,s(1,3,j),ones,zeros,1,3)
         do i=j,ms ! Form the Gramian, U=U_K^T A U_K using H^1_0 Norm
            gram(i,j)=s1*glsc2(uu,s(1,1,i),n)
     $               +s2*glsc3(s(1,1,i),s(1,1,j),bm1,n)
            if (mdim.ge.2) then
               gram(i,j)=gram(i,j)+s1*glsc2(vv,s(1,2,i),n)
     $                            +s2*glsc3(s(1,2,i),s(1,2,j),bm1,n)
            endif
            if (mdim.eq.3) then
               gram(i,j)=gram(i,j)+s1*glsc2(ww,s(1,3,i),n)
     $                            +s2*glsc3(s(1,3,i),s(1,3,j),bm1,n)
            endif
            if (i.ne.j) gram(j,i)=gram(i,j)
         enddo
         if (nio.eq.0) write(6,1) j,gram(1,j),'HLM'
      enddo

      if (nio.eq.0) write (6,*) 'exiting gengram'
    1 format (' gram',i5,1p1e16.6,2x,a3)

      return
      end
c-----------------------------------------------------------------------
      subroutine h10gg(gram,s,ms,mdim)

      ! set the Gramian based on the H^1_0 inner-product

      ! gram := Gramian
      ! s    := snapshots
      ! ms   := number of snapshots
      ! mdim := vector dimension

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrgg/ uu(lt),vv(lt),ww(lt)

      real gram(ms,ms)
      real s(lt,mdim,ms)

      if (nio.eq.0) write (6,*) 'inside gengram H10'

      n=lx1*ly1*lz1*nelt

      do j=1,ms
         call axhelm(uu,s(1,1,j),ones,zeros,1,1)
         if (mdim.ge.2) call axhelm(vv,s(1,2,j),ones,zeros,1,2)
         if (mdim.eq.3) call axhelm(ww,s(1,3,j),ones,zeros,1,3)
         do i=j,ms ! Form the Gramian, U=U_K^T A U_K using H^1_0 Norm
            gram(i,j)=glsc2(uu,s(1,1,i),n)
            if (mdim.ge.2) then
               gram(i,j)=gram(i,j)+glsc2(vv,s(1,2,i),n)
            endif
            if (mdim.eq.3) then
               gram(i,j)=gram(i,j)+glsc2(ww,s(1,3,i),n)
            endif
            if (i.ne.j) gram(j,i)=gram(i,j)
         enddo
         if (nio.eq.0) write(6,1) j,gram(1,j),'H10'
      enddo

      if (nio.eq.0) write (6,*) 'exiting gengram H10'
    1 format (' gram',i5,1p1e16.6,2x,a3)

      return
      end
c-----------------------------------------------------------------------
      subroutine wl2gg(gram,s,ms,mdim)

      ! set the Gramian based on the L^2 inner-product
      ! (duplicate, will be deprecated)

      ! gram := Gramian
      ! s    := snapshots
      ! ms   := number of snapshots
      ! mdim := vector dimension

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrgg/ uu(lt),vv(lt),ww(lt)

      real gram(ms,ms)
      real s(lt,mdim,ms)

      if (nio.eq.0) write (6,*) 'inside gengram'

      n=lx1*ly1*lz1*nelt

      do j=1,ms
         do i=j,ms ! Form the Gramian, U=U_K^T A U_K using H^1_0 Norm
            if (mdim.eq.1) then
               gram(i,j)=sip(s(1,1,i),s(1,1,j))
            else
               gram(i,j)=vip(s(1,1,i),s(1,2,i),s(1,ldim,i),
     $                       s(1,1,j),s(1,2,j),s(1,ldim,j))
            endif
            if (i.ne.j) gram(j,i)=gram(i,j)
         enddo
         if (nio.eq.0) write(6,1) j,gram(1,j),'L2 '
      enddo

      if (nio.eq.0) write (6,*) 'exiting gengram'
    1 format (' gram',i5,1p1e16.6,2x,a3)

      return
      end
c-----------------------------------------------------------------------
      subroutine gengraml2(gram,s,ms,mdim)

      ! set the Gramian based on the L^2 inner-product

      ! gram := Gramian
      ! s    := snapshots
      ! ms   := number of snapshots
      ! mdim := vector dimension

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrgram/ uw(lt),vw(lt),ww(lt)
      real gram(ms,ms)
      real s(lt,mdim,ms)

      if (nio.eq.0) write (6,*) 'inside gengraml2'

      n=lx1*ly1*lz1*nelv

      do j=1,ms ! Form the Gramian, U=U_K^T A U_K using L2 Norm
      do i=j,ms
         gram(i,j)=glsc3(s(1,1,i),s(1,1,j),bm1,n)
         if (mdim.ge.2)
     $      gram(i,j)=gram(i,j)+glsc3(s(1,2,i),s(1,2,j),bm1,n)
         if (mdim.ge.3)
     $      gram(i,j)=gram(i,j)+glsc3(s(1,3,i),s(1,3,j),bm1,n)
         if (i.ne.j) gram(j,i)=gram(i,j)
      enddo
         if (nio.eq.0) write (6,1) j,gram(1,j)
      enddo

      if (nio.eq.0) write (6,*) 'exiting gengraml2'

    1 format (' gram',i5,' ',1p1e16.6)

      return
      end
c-----------------------------------------------------------------------
      subroutine gengram(gram,s,ms,mdim,cips)

      ! set the Gramian based on the inner-product set by ips

      ! gram := Gramian
      ! s    := snapshots
      ! ms   := number of snapshots
      ! mdim := vector dimension
      ! cips := inner-product space specifier

      real gram(1),s(1)
      character*3 cips

      if (cips.eq.'L2 ') then
         call gengraml2(gram,s,ms,mdim)
      else if (cips.eq.'H10') then
         call h10gg(gram,s,ms,mdim)
      else if (cips.eq.'HLM') then
         call hlmgg(gram,s,ms,mdim)
      else
         if (nid.eq.0) write (6,*) 'unsupported ips in gengram'
         call exitti('failed in gengram, exiting...$',1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine genevec(vec,val,gg,ms,mb,ifld)

      ! set the eigenvectors based on the given Gramian

      ! vec  := eigenvectors
      ! val  := eigenvalues
      ! gg   := Gramian
      ! ms   := number of snapshots
      ! mb   := number of basis
      ! ifld := field number

      include 'SIZE'
      include 'TSTEP'
      include 'LMOR'
      
      common /scrgvec/ gc(ls*ls),wk(ls*ls*2)

      real gg(ms,ms),vec(ms,mb),val(ms)

      if (nio.eq.0) write (6,*) 'inside genevec'

      call nekgsync
      eig_time=dnekclock()

      call copy(gc,gg,ms*ms)

      call regularev(gc,val,ms,wk)
      call copy(wk,gc,ms*ms)

      call nekgsync
      eval_time=dnekclock()

      do l = 1,mb
         call copy(vec(1,l),wk(1+(ms-l)*ms),ms) ! reverse order of wk
      enddo

      do i=1,ns
         if (nio.eq.0) write (6,'(i5,1p1e16.6,3x,a,i1)')
     $      i,val(ns-i+1),'eval',ifld
      enddo

      call nekgsync
      evec_time=dnekclock()

      if (nio.eq.0) write (6,*) 'eval_time:',eval_time-eig_time
      if (nio.eq.0) write (6,*) 'evec_time:',evec_time-eval_time

      if (nio.eq.0) write (6,*) 'exiting genevec'

      return
      end
c-----------------------------------------------------------------------
      subroutine vnorm(uub,vvb,wwb)

      ! normalizes vector field

      ! uub,vvb,wwb := x,y,z components of vector field

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real uub(lt,0:nb),vvb(lt,0:nb),wwb(lt,0:nb)

      jfield=ifield
      ifield=1
      nio=-1
      do i=1,nb
         p=vip(uub(1,i),vvb(1,i),wwb(1,i),uub(1,i),vvb(1,i),wwb(1,i))
         s=1./sqrt(p)
         call opcmult(uub(1,i),vvb(1,i),wwb(1,i),s)
      enddo
      nio=nid
      ifield=jfield

      return
      end
c-----------------------------------------------------------------------
      subroutine cnorm(uub,vvb,wwb,ttb)

      ! normalizes combined field

      ! uub,vvb,wwb,ttb := vx,vy,vz,t components of vector field

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real uub(lt,0:nb),vvb(lt,0:nb),wwb(lt,0:nb),ttb(lt,0:nb)

      jfield=ifield
      ifield=1
      nio=-1
      do i=1,nb
         p=cip(uub(1,i),vvb(1,i),wwb(1,i),ttb(1,i),
     $         uub(1,i),vvb(1,i),wwb(1,i),ttb(1,i))
         s=1./sqrt(p)
         call opcmult(uub(1,i),vvb(1,i),wwb(1,i),s)
         call cmult(ttb(1,i),s,lx1*ly1*lz1*nelt)
      enddo
      nio=nid
      ifield=jfield

      return
      end
c-----------------------------------------------------------------------
      subroutine snorm(ssb)

      ! normalizes scalar field

      ! ssb := scalar field

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ssb(lt,0:nb)

      nio=-1
      do i=1,nb
         p=sip(ssb(1,i),ssb(1,i))
         s=1./sqrt(p)
         call cmult(ssb(1,i),s,lx1*ly1*lz1*nelt)
      enddo
      nio=nid

      return
      end
c-----------------------------------------------------------------------
      subroutine h10pv2b(coef,ux,uy,uz,uub,vvb,wwb)

      ! get coordinates of a vector field for a H^1_0 orthogonal basis

      ! coef        := coordinates of <ux,uy,uz> in <uub,vvb,wwb>
      ! ux,uy,uz    := x,y,z components of FOM field
      ! uub,vvb,wwb := x,y,z components of basis functions

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ux(lt),uy(lt),uz(lt),uub(lt,0:nb),vvb(lt,0:nb),wwb(lt,0:nb)

      common /scrp/ t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)

      real coef(0:nb)

      if (nio.eq.0) write (6,*) 'inside h10vpv2b'

      n=lx1*ly1*lz1*nelt

      call opcopy(t1,t2,t3,ux,uy,uz)

      coef(0) = 1.
      if (nio.eq.0) write (6,1) coef(0),coef(0),1.

      do i=1,nb
         call axhelm(t4,uub(1,i),ones,zeros,1,1)
         call axhelm(t5,vvb(1,i),ones,zeros,1,2)

         ww = glsc2(t4,uub(1,i),n)+glsc2(t5,vvb(1,i),n)
         vv = glsc2(t4,t1,n)+glsc2(t5,t2,n)

         if (ldim.eq.3) then
            call axhelm(t6,wwb(1,i),ones,zeros,1,3)
            ww = ww + glsc2(t6,wwb(1,i),n)
            vv = vv + glsc2(t6,t3,n)
         endif

         coef(i) = vv/ww
         if (nio.eq.0) write (6,1) coef(i),vv,ww
      enddo

      if (nio.eq.0) write (6,*) 'exiting h10pv2b'

    1 format(' h10coef',1p3e16.8)

      return
      end
c-----------------------------------------------------------------------
      subroutine hlmpv2b(coef,ux,uy,uz,uub,vvb,wwb)

      ! get coordinates of a vector field for a Helmholtz orthogonal basis

      ! coef        := coordinates of <ux,uy,uz> in <uub,vvb,wwb>
      ! ux,uy,uz    := x,y,z components of FOM field
      ! uub,vvb,wwb := x,y,z components of basis functions

      include 'SIZE'
      include 'MASS'
      include 'INPUT'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ux(lt),uy(lt),uz(lt),uub(lt,0:nb),vvb(lt,0:nb),wwb(lt,0:nb)

      common /scrp/ t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)

      real coef(0:nb)

      if (nio.eq.0) write (6,*) 'inside hlmvpv2b'

      n=lx1*ly1*lz1*nelt

      call opcopy(t1,t2,t3,ux,uy,uz)

      coef(0) = 1.
      if (nio.eq.0) write (6,1) coef(0),coef(0),1.

      s1=1./ad_re
      s2=ad_beta(1,3)/ad_dt

      do i=1,nb
         call axhelm(t4,uub(1,i),ones,zeros,1,1)
         call axhelm(t5,vvb(1,i),ones,zeros,1,2)

         ww=s1*(glsc2(t4,uub(1,i),n)+glsc2(t5,vvb(1,i),n))
         vv=s1*(glsc2(t4,t1,n)+glsc2(t5,t2,n))

         vv=vv+s2*(glsc3(t1,uub(1,i),bm1,n)+glsc3(t2,vvb(1,i),bm1,n))
         ww=ww+s2*(glsc3(uub(1,i),uub(1,i),bm1,n)
     $            +glsc3(vvb(1,i),vvb(1,i),bm1,n))

         if (ldim.eq.3) then
            call axhelm(t6,wwb(1,i),ones,zeros,1,3)
            ww=ww+s1*glsc2(t6,wwb(1,i),n)
            vv=vv+s1*glsc2(t6,t3,n)
            vv=vv+s2*glsc3(t3,wwb(1,i),bm1,n)
            ww=ww+s2*glsc3(wwb(1,i),wwb(1,i),bm1,n)
         endif

         coef(i) = vv/ww
         if (nio.eq.0) write (6,1) coef(i),vv,ww
      enddo

      if (nio.eq.0) write (6,*) 'exiting hlmpv2b'

    1 format(' hlmcoef',1p3e16.8)

      return
      end
c-----------------------------------------------------------------------
      subroutine regularev(a,lam,n,wk)
 
c     Solve the eigenvalue problem  A x = lam x
c
c     A -- symmetric matrix
c
c     "SIZE" is included here only to deduce WDSIZE, the working
c     precision, in bytes, so as to know whether dsygv or ssygv
c     should be called.

      include 'SIZE'
      include 'PARALLEL'

      real a(n,n),lam(n),wk(n,n)
      real aa(100)

      call copy(aa,a,100)

      call dsyev('V','U',n,a,n,lam,wk,2*n*n,info)

      if (info.ne.0) then
         if (nid.eq.0) then
            call outmat2(aa ,n,n,n,'aa  ')
            call outmat2(a  ,n,n,n,'Aeig')
            call outmat2(lam,1,n,n,'Deig')
         endif

         ninf = n-info
         write(6,*) 'Error in regularev, info=',info,n,ninf
         call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cnmax(val,fname,ifld)

      ! compute maximum number of POD modes to be used

      ! val   := eigenvalues based on the generated Gramians
      ! fname := filename that will be created
      ! ifld  := field number

      include 'SIZE'
      include 'TSTEP'
      include 'MOR'

      integer icalld
      save    icalld
      data    icalld /0/
      
      real val(ls)
      real total,ena(ls),enl(ls)
      character*128 fname

      if (nio.eq.0) write (6,*) 'inside cnmax'

      call nekgsync

      if (icalld.eq.0) then
         icalld=1
      endif

c     total=vlsum(val,ls)
      total=0
      do i=1,ns
         total=total+val(ns-i+1) 
         if (nio.eq.0) write(6,*)i,total,val(ns-i+1),'total'
      enddo
      ena(1) = val(ns)
      enl(1) = sqrt(total-ena(1))/sqrt(total)
      do i=2,ns
        ena(i)=ena(i-1)+val(ns-i+1)
        enl(i) = sqrt(total-ena(i))/sqrt(total)
        if (enl(i).le.1e-4) then 
        if (icalld.eq.1) then
           if (nio.eq.0) write(6,*)i,enl(i),'cnmax Nmax for field',ifld
           icalld=2
        endif
        endif
      enddo
      
      call dump_serial(enl,ns,fname,nid)

      call nekgsync

      icalld=0

      if (nio.eq.0) write (6,*) 'exiting cnmax'

      return
      end
c-----------------------------------------------------------------------
      subroutine cenpm(val,fname,ifld)

      ! compute the percentage that each POD mode represents in the
      ! total averaged energy. The energy is defined based on the
      ! inner-product set by ips

      ! val   := eigenvalues based on the generated Gramians
      ! fname := filename that will be created
      ! ifld  := field number

      include 'SIZE'
      include 'TSTEP'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      integer icalld
      save    icalld
      data    icalld /0/
      
      real val(ls)
      real total,enr(ls)
      character*128 fname

      if (nio.eq.0) write (6,*) 'inside cenpm'

      call nekgsync

      if (icalld.eq.0) then
         icalld=1
      endif

      total=0
      do i=1,ns
         total=total+val(ns-i+1) 
         if (nio.eq.0) write(6,*)i,total,val(ns-i+1),'total'
      enddo
      do i=1,ns
        enr(i)=val(ns-i+1)/total
c       if (nio.eq.0) write(6,*)i,enr(i),'Nmax for field',ifld
        if (enr(i).le.1e-3) then 
        if (icalld.eq.1) then
           if (nio.eq.0) write(6,*)i,enr(i),'cenpm Nmax for field',ifld
           icalld=2
        endif
        endif
      enddo
      
      call dump_serial(enr,ns,fname,nid)

      call nekgsync

      icalld=0

      if (nio.eq.0) write (6,*) 'exiting cenpm'

      return
      end
c-----------------------------------------------------------------------
      subroutine setbases_p(nocp)

      ! Compute the reduced-bases for p_greedy
      ! First ncop basis are from pbas.list
      ! ncop+1 to nb constructed in this subroutine

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real u0(lt,3)
      integer nocp

      if (nio.eq.0) write (6,*) 'inside setbases_p'

      call nekgsync
      bas_time=dnekclock()

      n=lx1*ly1*lz1*nelt

      if (rmode.ne.'ON ') ifrecon=.true.

      if (rmode.eq.'ONB'.or.rmode.eq.'CP ') then
         call loadbases
      else if (rmode.eq.'ALL'.or.rmode.eq.'OFF') then
         do i=nocp+1,nb
            call rzero(ub(1,i),n)
            call rzero(vb(1,i),n)
            if (ldim.eq.3) call rzero(wb(1,i),n)
         enddo

            ! ub, vb, wb, are the modes
            if (ifrom(1)) then
              do j=1,ns
              do i=1,nb-nocp
                 call opadds(ub(1,i+nocp),vb(1,i+nocp),wb(1,i+nocp),
     $             us0(1,1,j),us0(1,2,j),us0(1,ldim,j),evec(j,i,1),n,2)
              enddo
              enddo

              call vnorm(ub,vb,wb)
            else
               call opcopy(ub,vb,wb,uic,vic,wic)
            endif

            if (ifrom(2)) then
               do i=1,nb-nocp
                  call rzero(tb(1,i+nocp),n)
                  do j=1,ns
                     call add2s2(tb(1,i+nocp),ts0(1,j),evec(j,i,2),n)
                  enddo
               enddo
               call snorm(tb)
            endif
         endif

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF') call dump_bas

      call nekgsync
      if (nio.eq.0) write (6,*) 'bas_time:',dnekclock()-bas_time

      if (nio.eq.0) write (6,*) 'exiting setbases_p'

      return
      end
c-----------------------------------------------------------------------
      subroutine pod(basis,evec,eval,gram,snaps,mdim,cips,nb,ns,ifpod)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)

      real basis(lt,mdim,nb),snaps(lt,mdim,ns),
     $   evec(ns,ns),eval(ns),gram(ns,ns)

      character*3 cips
      logical ifpod

      n=lx1*ly1*lz1*nelt

      call gengram(gram,snaps,ns,mdim,cips)

      if (ifpod) then
         call genevec(evec,eval,gram,ns,nb,mdim)
         do i=1,mdim
            call dgemm('N','N',n,nb,ns,1.,
     $         snaps(1,i,1),lt*mdim,evec,ns,0.,basis(1,i,1),lt*mdim)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
