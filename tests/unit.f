c-----------------------------------------------------------------------
      ! Unit Tests
c-----------------------------------------------------------------------
      subroutine gramian_l2_unit
      call gramian_unit(.true.)
      return
      end
c-----------------------------------------------------------------------
      subroutine gramian_h10_unit
      call gramian_unit(.false.)
      return
      end
c-----------------------------------------------------------------------
      subroutine gramian_unit(iflag)

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      common /scrtest/ wk(ls*ls)

      parameter (lt=lx1*ly1*lz1*lelt)

      logical iflag
      real vv(ls,ls)

      param(171) = 1.
      if (iflag) param(171) = 0.
      param(172) = 1.
      param(173) = 0.

      call rom_setup
      call gengram(ug,us0,ns,ldim,ips)

      call read_serial(vv,ls*ls,'tops/gu ',wk,nid)

      iexit=0

      s1=0.
      s2=0.
      s3=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do j=1,ls
      do i=1,ls
         s1=s1+(ug(i,j)-ug(j,i))**2
         s2=s2+(ug(i,j)-vv(i,j))**2
         s3=s3+vv(i,j)**2
         if (nio.eq.0) write (6,*) 'gram',i,j,ug(i,j),vv(i,j)
      enddo
      enddo

      esym=sqrt(s1/s3)
      edif=sqrt(s2/s3)

      if (nio.eq.0) write (6,*) 'esym',esym,s1,s3
      if (nio.eq.0) write (6,*) 'edif',edif,s2,s3

      if (ips.eq.'L2 '.and.esym.gt.1e-16) iexit=iexit+1
      if (ips.eq.'H10'.and.esym.gt.1e-14) iexit=iexit+1
      if (edif.gt.5.e-15) iexit=iexit+2

      call exitm(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine ic_l2_unit
      call ic_unit(.true.)
      return
      end
c-----------------------------------------------------------------------
      subroutine ic_h10_unit
      call ic_unit(.false.)
      return
      end
c-----------------------------------------------------------------------
      subroutine ic_unit(iflag)

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MOR'

      common /scrtest/ wk(lb+1),u0(0:lb)

      logical iflag

      param(171) = 1.
      if (iflag) param(171) = 0.
      param(172) = 1.
      param(173) = 0.

      call rom_setup
      call read_serial(u0,nb+1,'tops/u ',wk,nid)

      s1=0.
      s2=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do i=0,nb
         write (6,*) i,u(i),u0(i),'ic'
      enddo

      do i=1,nb
         s1=s1+(u0(i)-u(i))**2
         s2=s2+u0(i)**2
      enddo

      edif=sqrt(s1/s2)
      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      iexit=1
      if (edif.lt.8.e-13) iexit=0

      call exitm(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine a0_l2_unit
      call a0_unit(.true.)
      return
      end
c-----------------------------------------------------------------------
      subroutine a0_h10_unit
      call a0_unit(.false.)
      return
      end
c-----------------------------------------------------------------------
      subroutine a0_unit(iflag)

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      logical iflag

      param(171) = 1.
      if (iflag) param(171) = 0.
      param(172) = 1.
      param(173) = 0.

      call rom_setup
      call a0_unit_helper(au0)

      return
      end
c-----------------------------------------------------------------------
      subroutine a0_unit_helper(a0)

      include 'SIZE'
      include 'MOR'

      common /scrtest/ wk(lb+1,lb+1),aa(0:lb,0:lb)

      real a0(0:nb,0:nb)
      call read_serial(aa,(nb+1)**2,'tops/au ',wk,nid)

      iexit=0

      s1=0.
      s2=0.
      s3=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do j=0,nb
      do i=0,nb
         if (nio.eq.0) write (6,*) i,j,a0(i,j),aa(i,j),'a'
         s1=s1+(aa(i,j)-a0(i,j))**2
         s2=s2+(a0(i,j)-a0(j,i))**2
         s3=s3+aa(i,j)**2
      enddo
      enddo

      edif=sqrt(s1/s3)
      esym=sqrt(s2/s3)

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s3
      if (nio.eq.0) write (6,*) 'esym',esym,s2,s3

      if (edif.gt.5.e-11) iexit=iexit+1
      if (esym.gt.1.e-15) iexit=iexit+2

      s1=0.
      s2=0.

      do j=1,nb
      do i=1,nb
         if (i.ne.j) s1=s1+a0(i,j)**2
         s2=s2+a0(i,j)**2
      enddo
      enddo

      edia=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'edia',edia,s1,s2
      if (ips.eq.'H10'.and.edia.gt.1.e-11) iexit=iexit+4

      s1=0.

      do i=1,nb
         s1=s1+(a0(i,i)-1.)**2
      enddo

      euni=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'euni',euni,s1,s2
      if (ips.eq.'H10'.and.euni.gt.1.e-13) iexit=iexit+8

      call exitm(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine b0_l2_unit
      call b0_unit(.true.)
      return
      end
c-----------------------------------------------------------------------
      subroutine b0_h10_unit
      call b0_unit(.false.)
      return
      end
c-----------------------------------------------------------------------
      subroutine b0_unit(iflag)

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      logical iflag

      param(171) = 1.
      if (iflag) param(171) = 0.
      param(172) = 1.
      param(173) = 0.

      call rom_setup
      call b0_unit_helper(bu0)

      return
      end
c-----------------------------------------------------------------------
      subroutine b0_unit_helper(b0)

      include 'SIZE'
      include 'SOLN'
      include 'MOR'

      common /scrtest/ wk(lb+1,lb+1),bb(0:lb,0:lb)

      real b0(0:nb,0:nb)

      call read_serial(bb,(nb+1)**2,'tops/bu ',wk,nid)

      iexit=0

      s1=0.
      s2=0.
      s3=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do j=0,nb
      do i=0,nb
         if (nio.eq.0) write (6,*) i,j,b0(i,j),bb(i,j),'b'
         s1=s1+(bb(i,j)-b0(i,j))**2
         s2=s2+(b0(i,j)-b0(j,i))**2
         s3=s3+bb(i,j)**2
      enddo
      enddo

      edif=sqrt(s1/s3)
      esym=sqrt(s2/s3)

      if (edif.gt.4.e-13) iexit=iexit+1
      if (esym.gt.1.e-16) iexit=iexit+2

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s3
      if (nio.eq.0) write (6,*) 'esym',esym,s2,s3

      s1=0.
      s2=0.

      do j=1,nb
      do i=1,nb
         if (i.ne.j) s1=s1+b0(i,j)**2
         s2=s2+bb(i,j)**2
      enddo
      enddo

      edia=sqrt(s1/s2)

      if (ips.eq.'L2 '.and.edia.gt.1.e-9) iexit=iexit+4
      if (nio.eq.0) write (6,*) 'edia',edia,s1,s2

      s1=0.

      do i=1,nb
         s1=s1+(b0(i,i)-1.)**2
      enddo

      euni=sqrt(s1/s2)

      if (ips.eq.'L2 '.and.euni.gt.1.e-14) iexit=iexit+8
      if (nio.eq.0) write (6,*) 'euni',euni,s1,s2

      call exitm(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine c0_l2_unit
      call c0_unit(.true.)
      return
      end
c-----------------------------------------------------------------------
      subroutine c0_h10_unit
      call c0_unit(.false.)
      return
      end
c-----------------------------------------------------------------------
      subroutine c0_unit(iflag)

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      logical iflag

      param(171) = 1.
      if (iflag) param(171) = 0.
      param(172) = 1.
      param(173) = 0.

      call rom_setup
      call c0_unit_helper(cul)

      return
      end
c-----------------------------------------------------------------------
      subroutine c0_unit_helper(c,iflag)

      include 'SIZE'
      include 'SOLN'
      include 'MOR'

      common /scrtest/ wk(lb+1,lb+1,lb+1),cc(lb,0:lb,0:lb)

      real c(nb,0:nb,0:nb)

      call read_serial(cc,nb*(nb+1)**2,'tops/cu ',wk,nid)

      iexit=0

      s1=0.
      s2=0.
      s3=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do k=0,nb
      do j=0,nb
      do i=1,nb
         s1=s1+(cc(i,j,k)-c(i,j,k))**2
         s3=s3+cc(i,j,k)**2
         write (6,*) i,j,k,c(i,j,k),cc(i,j,k),'c'
      enddo
      enddo
      enddo

      do k=1,nb
      do j=1,nb
      do i=1,nb
         s2=s2+(c(i,j,k)+c(j,i,k))**2
      enddo
      enddo
      enddo

      edif=sqrt(s1/s3)
      if (edif.gt.2.e-10) iexit=iexit+1
      if (nio.eq.0) write (6,*) 'edif',edif,s1,s3

      eskew=sqrt(s2/s3)
c     if (eskew.gt.1.e-16) iexit=iexit+2
      if (nio.eq.0) write (6,*) 'eskew',eskew,s2,s3

      call exitm(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine cp_l2_unit
      call cp_unit(.true.)
      return
      end
c-----------------------------------------------------------------------
      subroutine cp_h10_unit
      call cp_unit(.false.)
      return
      end
c-----------------------------------------------------------------------
      subroutine cp_unit(iflag)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrtens_norm/ norm_c,norm_c0

      logical iflag
      real wk(lb*ltr),wk3(lt)
      real cp_a(lb*ltr),cp_b(lb*ltr),cp_c(lb*ltr),cp_w(ltr)
      real cl0(lcglo)
      real norm_c,norm_c0

      param(171) = 1.
      if (iflag) param(171) = 0.
      param(172) = 1.
      param(173) = 0.
      nb = 20

      call rom_setup

      ntr = 80
      ifcore=.true.
      if (nio.eq.0) write (6,*) 'reading initial A2...'
      call read_mat_serial(cub,nb,ntr,'./tcp/cub0 ',nb,ntr,wk3,nid)

      if (nio.eq.0) write (6,*) 'reading initial A3...'
      call read_mat_serial(cuc,nb,ntr,'./tcp/cuc0 ',nb,ntr,wk3,nid)

      local_size = (ic2-ic1+1)*(jc2-jc1+1)*(kc2-kc1+1)
      norm_c = vlsc2(cul,cul,local_size)
      if (ifcore) then
         call set_c_slice(cuj0,cu0k,cul,ic1,ic2,jc1,jc2,kc1,kc2,nb)
         call set_c_core(cl0,cul,ic1,ic2,jc1,jc2,kc1,kc2,nb)
         local_size = (ic2-ic1+1)*(jc2-jc1)*(kc2-kc1)
         norm_c0 = vlsc2(cl0,cl0,local_size)
      endif
      call als_core(cua,cub,cuc,cp_uw,cl0,
     $              ic1,ic2,jc1+1,jc2,kc1+1,kc2,nb,ntr)

      call read_serial(cp_a,nb*ntr,'./tcp/cua ',wk,nid)
      call read_serial(cp_b,nb*ntr,'./tcp/cub ',wk,nid)
      call read_serial(cp_c,nb*ntr,'./tcp/cuc ',wk,nid)
      call read_serial(cp_w,ntr   ,'./tcp/cuw ',wk,nid)

      iexit=0

      s1=0.
      s2=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do j=1,ntr
      do i=1,nb
         s1=s1+(cua(i+(j-1)*nb)-cp_a(i+(j-1)*nb))**2
         s2=s2+cp_a(i+(j-1)*nb)**2
         if (nio.eq.0) write (6,*) 'cp a factor',
     $                 i,j,cua(i+(j-1)*nb),cp_a(i+(j-1)*nb)
      enddo
      enddo

      edif=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      if (edif.gt.9.e-11) iexit=iexit+1

      s1=0.
      s2=0.

      do j=1,ntr
      do i=1,nb
         s1=s1+(cub(i+(j-1)*nb)-cp_b(i+(j-1)*nb))**2
         s2=s2+cp_b(i+(j-1)*nb)**2
         if (nio.eq.0) write (6,*) 'cp b factor',
     $                 i,j,cub(i+(j-1)*nb),cp_b(i+(j-1)*nb)
      enddo
      enddo

      edif=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      if (edif.gt.9.e-11) iexit=iexit+2

      s1=0.
      s2=0.

      do j=1,ntr
      do i=1,nb
         s1=s1+(cuc(i+(j-1)*nb)-cp_c(i+(j-1)*nb))**2
         s2=s2+cp_c(i+(j-1)*nb)**2
         if (nio.eq.0) write (6,*) 'cp c factor',
     $                 i,j,cuc(i+(j-1)*nb),cp_c(i+(j-1)*nb)
      enddo
      enddo

      edif=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      if (edif.gt.9.e-11) iexit=iexit+4

      s1=0.
      s2=0.

      do j=1,ntr
         s1=s1+(cp_uw(j)-cp_w(j))**2
         s2=s2+cp_w(j)**2
         if (nio.eq.0) write (6,*) 'cp w factor',
     $                 j,cp_uw(j),cp_w(j)
      enddo

      edif=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      if (edif.gt.9.e-11) iexit=iexit+8

      call exitm(iexit)
      
      return
      end
c-----------------------------------------------------------------------
      subroutine evalc_l2_unit
      call evalc_unit(.true.)
      return
      end
c-----------------------------------------------------------------------
      subroutine evalc_h10_unit
      call evalc_unit(.false.)
      return
      end
c-----------------------------------------------------------------------
      subroutine evalc_unit(iflag)

      include 'SIZE'
      include 'SOLN'
      include 'MOR'

      parameter (llb=100)
      common /scrtest_evalc/ cu(llb),wk(llb),
     $   c_ref(llb*(llb+1)**2),u_ref(llb+1),cu_ref(llb),
     $   tmp(llb*(llb+1)),c(llb*(llb+1)**2)

      call srand(123)

      do i=1,llb*(llb+1)**2
         c_ref(i)=rand()
      enddo

      do i=1,llb+1
         u_ref(i)=rand()
      enddo

      dl2max=0.
      cl2max=0.
      el2max=0.

      iexit=0.

      do mb=1,llb
         nb=mb
         do mp=1,32
            call rzero(cu,mb)
            do ip=1,mp
               call cpart(kc1,kc2,jc1,jc2,ic1,ic2,ncloc,mb,mp,ip)

               if (mp.eq.1) then
                  call rzero(cu_ref,mb)
                  call evalc(cu_ref,tmp,c_ref,u_ref,u_ref)
               endif

               do k=0,nb
               do j=0,mb
               do i=1,mb
                  cel=c_ref(i+j*mb+k*mb*(mb+1))
                  call setc_local(c,cel,ic1,ic2,jc1,jc2,kc1,kc2,i,j,k)
               enddo
               enddo
               enddo

               call rzero(wk,mb)
               call evalc(wk,tmp,c,u_ref,u_ref)

               call add2(cu,wk,mb)
            enddo
            call sub2(cu,cu_ref,mb)
            dl2=sqrt(vlsc2(cu,cu,mb))
            cl2=sqrt(vlsc2(cu_ref,cu_ref,mb))
            el2=dl2/cl2
            cl2max=max(cl2,cl2max)
            dl2max=max(dl2,dl2max)
            el2max=max(el2,el2max)
            write (6,*) mb,mp,dl2,cl2,el2,'error'
         enddo
      enddo

      write (6,*) dl2max,cl2max,el2max,'errormax'

      if (el2max.gt.3.5e-15) iexit=iexit+1

      call exitm(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine exitm(iexit)

      open (unit=10,file='ecode')
      write (10,*) iexit
      close (unit=10)

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine cp_skew_l2_unit
      call cp_skew_unit(.true.)
      return
      end
c-----------------------------------------------------------------------
      subroutine cp_skew_h10_unit
      call cp_skew_unit(.false.)
      return
      end
c-----------------------------------------------------------------------
      subroutine cp_skew_unit(iflag)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrtens_norm/ norm_c,norm_c0

      logical iflag
      real wk(lb*ltr),wk3(lt)
      real cp_a(lb*ltr),cp_b(lb*ltr),cp_c(lb*ltr),cp_w(ltr)
      real cl0(lcglo)
      real norm_c,norm_c0

      param(171) = 1.
      if (iflag) param(171) = 0.
      param(172) = 1.
      param(173) = 0.
      nb = 20

      call rom_setup

      ntr = 80
      ifcore=.true.
      ifquad=.true.
      if (nio.eq.0) write (6,*) 'reading initial A1...'
      call read_mat_serial(cua,nb,ntr,'./tcp/cua_r80skew ',
     $                     nb,ntr,wk3,nid)

      if (nio.eq.0) write (6,*) 'reading initial A2...'
      call read_mat_serial(cub,nb,ntr,'./tcp/cub_r80skew ',
     $                     nb,ntr,wk3,nid)

      if (ifquad) then
         call force_skew(cul,ic1,ic2,jc1,jc2,kc1,kc2)
      endif
      local_size = (ic2-ic1+1)*(jc2-jc1+1)*(kc2-kc1+1)
      norm_c = vlsc2(cul,cul,local_size)
      if (ifcore) then
         call set_c_slice(cuj0,cu0k,cul,ic1,ic2,jc1,jc2,kc1,kc2,nb)
         call set_c_core(cl0,cul,ic1,ic2,jc1,jc2,kc1,kc2,nb)
         local_size = (ic2-ic1+1)*(jc2-jc1)*(kc2-kc1)
         norm_c0 = vlsc2(cl0,cl0,local_size)
      endif
      call als_quad(cua,cub,cuc,cp_uw,cl0,
     $              ic1,ic2,jc1+1,jc2,kc1+1,kc2,nb,ntr)

      call read_serial(cp_a,nb*ntr,'./tcp/cua_skew ',wk,nid)
      call read_serial(cp_b,nb*ntr,'./tcp/cub_skew ',wk,nid)
      call read_serial(cp_c,nb*ntr,'./tcp/cuc_skew ',wk,nid)
      call read_serial(cp_w,ntr   ,'./tcp/cuw_skew ',wk,nid)

      iexit=0

      s1=0.
      s2=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do j=1,ntr
      do i=1,nb
         s1=s1+(cua(i+(j-1)*nb)-cp_a(i+(j-1)*nb))**2
         s2=s2+cp_a(i+(j-1)*nb)**2
         if (nio.eq.0) write (6,*) 'cp a factor',
     $                 i,j,cua(i+(j-1)*nb),cp_a(i+(j-1)*nb)
      enddo
      enddo

      edif=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      if (edif.gt.2.e-11) iexit=iexit+1

      s1=0.
      s2=0.

      do j=1,ntr
      do i=1,nb
         s1=s1+(cub(i+(j-1)*nb)-cp_b(i+(j-1)*nb))**2
         s2=s2+cp_b(i+(j-1)*nb)**2
         if (nio.eq.0) write (6,*) 'cp b factor',
     $                 i,j,cub(i+(j-1)*nb),cp_b(i+(j-1)*nb)
      enddo
      enddo

      edif=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      if (edif.gt.2.e-11) iexit=iexit+2

      s1=0.
      s2=0.

      do j=1,ntr
      do i=1,nb
         s1=s1+(cuc(i+(j-1)*nb)-cp_c(i+(j-1)*nb))**2
         s2=s2+cp_c(i+(j-1)*nb)**2
         if (nio.eq.0) write (6,*) 'cp c factor',
     $                 i,j,cuc(i+(j-1)*nb),cp_c(i+(j-1)*nb)
      enddo
      enddo

      edif=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      if (edif.gt.2.e-11) iexit=iexit+4

      s1=0.
      s2=0.

      do j=1,ntr
         s1=s1+(cp_uw(j)-cp_w(j))**2
         s2=s2+cp_w(j)**2
         if (nio.eq.0) write (6,*) 'cp w factor',
     $                 j,cp_uw(j),cp_w(j)
      enddo

      edif=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      if (edif.gt.3.e-11) iexit=iexit+8

      call exitm(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine rbf_wt_unit

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      common /scrtest/ wk((lb+1)*ls)

      real vv(ls*lb)

      ifedvs = .true.

      anch(1) = 100000
      anch(2) = 750000
      anch(3) = 400000
      rbf_sigma(1) = 650000
      rbf_sigma(2) = 650000
      rbf_sigma(3) = 650000

      call read_serial(edk,ls*(lb+1),'tops/edk ',wk,nid)
      call read_serial(vv,ls*lb,'tops/rbfwt ',wk,nid)
      call c_rbfwt(rbfwt,rbf_sigma,edk,anch,ls,lb)

      iexit=0

      s1=0.
      s2=0.

      if (nio.eq.0) write (6,*) 'live | data'

      do j=1,ls
      do i=1,lb
         s1=s1+(rbfwt(i+(j-1)*lb)-vv(i+(j-1)*lb))**2
         s2=s2+(vv(i+(j-1)*lb))**2
         if (nio.eq.0) write(6,*)'rbfwt',i+(j-1)*lb,
     $                 rbfwt(i+(j-1)*lb),vv(i+(j-1)*lb)
      enddo
      enddo

      edif=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      if (edif.gt.3.e-11) iexit=iexit+1

      call exitm(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine rf_unit

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /romdbas/ tmp(lt,ldimt)

      character*128 fname1

      fname1='file.list '

      nt=lx1*ly1*lz1*nelt

      nb = lb
      ns = ls
      nskip = 0

      do i=0,ldimt1
         ifreads(i)=.true.
      enddo

      call read_fields(us0,prs,ts0,ns,ls,nskip,ifreads,
     $     timek,fname1,.true.)

      iftmp=ifxyo
      ifxyo=.true.
      do i=1,ns
         time=i
         itmp=i
         do j=1,ldimt
            call copy(tmp(1,j),ts0(1,i,j),nt)
         enddo
         param(66)=-1
         call outpost2(us0(1,1,i),us0(1,2,i),us0(1,ldim,i),prs(1,i),
     $        tmp,ldimt,'us0')
      enddo
      ifxyo=.false.

      return
      end
