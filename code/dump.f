c-----------------------------------------------------------------------
      subroutine dump_global(a,n,fname,wk1,wk2,nid)

      ! dump a distributed real array to a file

      ! a       := real data array
      ! n       := number of local entries to dump
      ! fname   := file name
      ! wk1,wk2 := work arrays
      ! nid     := processor id

      real a(n),wk1(1),wk2(1)

      character*128 fname
      character*128 fntrunc

      if (nid.eq.0) then
         call blank(fntrunc,128)
         len=ltruncr(fname,128)
         call chcopy(fntrunc,fname,len)
      endif

      call dump_global_helper(a,n,fntrunc,wk1,wk2,nid)

      return
      end
c-----------------------------------------------------------------------
      subroutine idump_serial(a,n,fname,nid)

      ! dump an integer array to a file

      ! a       := integer data array
      ! n       := number of entries to dump
      ! fname   := file name
      ! nid     := processor id

      integer a(n)

      character*128 fname
      character*128 fntrunc

      if (nid.eq.0) then
         call blank(fntrunc,128)

         len=ltruncr(fname,128)
         call chcopy(fntrunc,fname,len)

         call idump_serial_helper(a,n,fntrunc)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_serial(a,n,fname,nid)

      ! dump a real array to a file

      ! a       := integer data array
      ! n       := number of entries to dump
      ! fname   := file name
      ! nid     := processor id

      real a(n)

      character*128 fname
      character*128 fntrunc

      if (nid.eq.0) then
         call blank(fntrunc,128)

         len=ltruncr(fname,128)
         call chcopy(fntrunc,fname,len)

         call dump_serial_helper(a,n,fntrunc)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_global_helper(a,n,fname,wk1,wk2,nid)

      ! helper routine to dump_global

      ! a       := integer data array
      ! n       := number of local entries to dump
      ! fname   := file name
      ! wk1,wk2 := work arrays
      ! nid     := processor id

      real a(n),wk1(1),wk2(1)
      integer iwk(1)

      character*128 fname

      if (nid.eq.0) open (unit=12,file=fname)

      iwk(1)=n
      nmax=iglmax(iwk,1)

      iwk(1)=nid
      ipmax=iglmax(iwk,0)

      do ip=0,ipmax
         if (nid.eq.ip) then
            call copy(wk1,a,nmax)
            iwk(1)=n
         else
            call rzero(wk1,nmax)
            iwk(1)=0
         endif

         iwk(1)=iglmax(iwk,1)

         call gop(wk1,wk2,'+  ',nmax)

         if (nid.eq.0) then
            do i=1,iwk(1)
               write (12,*) wk1(i)
            enddo
         endif
      enddo

      if (nid.eq.0) close (unit=12)

      return
      end
c-----------------------------------------------------------------------
      subroutine idump_serial_helper(a,n,fname)

      ! helper routine for idump_serial

      ! a       := integer data array
      ! n       := number of entries to dump
      ! fname   := file name

      integer a(n)

      character*128 fname

      open (unit=12,file=fname)

      do i=1,n
         write (12,*) a(i)
      enddo

      close (unit=12)

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_serial_helper(a,n,fname)

      ! helper routine for dump_serial

      ! a       := real data array
      ! n       := number of entries to dump
      ! fname   := file name

      real a(n)

      character*128 fname

      open (unit=12,file=fname)

      do i=1,n
         write (12,1) a(i)
      enddo

      close (unit=12)
    1 format(1pe24.16)

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_all

      ! dump 'all' operators and data

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /dumpglobal/ wk1(lcloc),wk2(lcloc)

      logical iftmp1,iftmp2,iftmp3

      call nekgsync
      dump_time=dnekclock()

      if (ifpod(1)) then
         call dump_serial(au0,(nb+1)**2,'ops/au ',nid)
         call dump_serial(bu0,(nb+1)**2,'ops/bu ',nid)
         call dump_serial(u,(nb+1)*3,'ops/u ',nid)
         call dump_serial(uk,ns*(nb+1),'ops/uk ',nid)
         call dump_serial(umin,nb,'ops/umin ',nid)
         call dump_serial(umax,nb,'ops/umax ',nid)
         call dump_serial(timek,ns,'ops/timek ',nid)
         call dump_global(cul,ncloc,'ops/cu ',wk1,wk2,nid)

         if (ifcdrag) then
            call dump_serial(rdgx,nb+1,'qoi/rdgx ',nid)
            call dump_serial(rdgy,nb+1,'qoi/rdgy ',nid)
            if (ldim.eq.3) call dump_serial(rdgz,nb+1,'qoi/rdgz ',nid)

            call dump_serial(fd1,ldim*(nb+1),'qoi/fd1 ',nid)
            call dump_serial(fd2,ldim*(nb+1)**2,'qoi/fd2 ',nid)
            call dump_serial(fd3,ldim*(nb+1),'qoi/fd3 ',nid)
         endif
      endif

      if (ifpod(2)) then
         call dump_serial(at0,(nb+1)**2,'ops/at ',nid)
         call dump_serial(bt0,(nb+1)**2,'ops/bt ',nid)
         call dump_serial(ut,(nb+1)*3,'ops/t ',nid)
         call dump_serial(tk,ns*(nb+1),'ops/tk ',nid)
         call dump_serial(tmin,nb,'ops/tmin ',nid)
         call dump_serial(tmax,nb,'ops/tmax ',nid)
         if (.not.ifpod(1))
     $      call dump_serial(timek,ns,'ops/timek ',nid)
         call dump_global(ctl,ncloc,'ops/ct ',wk1,wk2,nid)
      endif

      if (ifforce)  call dump_serial(rf,nb,'ops/rf ',nid)
      if (ifsource) call dump_serial(rq,nb,'ops/rq ',nid)

      if (ifei) then
         l=1
         do j=1,nres
         do i=1,nres
            sigtmp(l,1)=mor_sigma(i,j)
            l=l+1
         enddo
         enddo
         call dump_serial(sigtmp,nres*nres,'ops/sigma ',nid)
      endif

      ttmp=time
      itmp=istep

      iftmp1=ifxyo
      iftmp2=ifpo
      iftmp3=ifto

      call nekgsync
      dbas_time=dnekclock()

      ifto=ifrom(2)

      do i=0,nb
         time=i
         itmp=i
         ifxyo=(i.eq.0)
         call outpost(ub(1,i),vb(1,i),wb(1,i),pb(1,i),tb(1,i,1),'bas')
      enddo

      istep=itmp
      time=ttmp

      ifxyo=iftmp1
      ifpo=iftmp2
      ifto=iftmp3

      call nekgsync
      done_time=dnekclock()
      if (nio.eq.0) write (6,*) 'dbas_time:',done_time-dbas_time
      if (nio.eq.0) write (6,*) 'dump_time:',done_time-dump_time

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_ops

      ! dump core operators (c out disabled)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /dumpglobal/ wk1(lcloc),wk2(lcloc)

      logical iftmp1,iftmp2,iftmp3

      call nekgsync
      dops_time=dnekclock()

      open (unit=10,file='ops/ips')
      if (nio.eq.0) write (10,*) ips
      close (unit=10)

      if (ifrom(1)) then
         call dump_serial(au0,(nb+1)**2,'ops/au ',nid)
         if (rmode.eq.'AEQ')
     $      call dump_serial(aue,nb*(nb+1)**2,'ops/aue ',nid)
         call dump_serial(bu0,(nb+1)**2,'ops/bu ',nid)
c        call dump_global(cul,ncloc,'ops/cu ',wk1,wk2,nid)
      endif

      if (ifrom(2)) then
         call dump_serial(at0,(nb+1)**2,'ops/at ',nid)
         if (rmode.eq.'AEQ')
     $      call dump_serial(ate,nb*(nb+1)**2,'ops/ate ',nid)
         call dump_serial(bt0,(nb+1)**2,'ops/bt ',nid)
         call dump_serial(st0,nb+1,'ops/st ',nid)
c        call dump_global(ctl,ncloc,'ops/ct ',wk1,wk2,nid)
      endif

      call nekgsync
      if (nio.eq.0) write (6,*) 'dops_time:',dnekclock()-dops_time

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_bas

      ! dump pod basis

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /dumpglobal/ wk1(lcloc),wk2(lcloc)
      common /romdbas/ tmp(lt,ldimt)

      logical iftmp1,iftmp2,iftmp3

      call nekgsync
      dbas_time=dnekclock()

      n=lx1*ly1*lz1*lelt
      ttmp=time
      itmp=istep

      iftmp1=ifxyo
      iftmp2=ifpo
      iftmp3=ifto

      ifpo=.false.
      ifto=ifrom(2)

      call nekgsync
      dbas_time=dnekclock()

      do i=0,nb
         time=i
         itmp=i
         ifxyo=(i.eq.0)
         do j=1,npscal+1
            call copy(tmp(1,j),tb(1,i,j),n)
         enddo
         call outpost2(ub(1,i),vb(1,i),wb(1,i),pb(1,i),tmp,ldimt,'bas')
      enddo

      istep=itmp
      time=ttmp

      ifxyo=iftmp1
      ifpo=iftmp2
      ifto=iftmp3

      rtmp1(1,1)=nb*1.
      call dump_serial(rtmp1(1,1),1,'ops/nb ',nid)

      rtmp1(1,1)=nbo*1.
      call dump_serial(rtmp1(1,1),1,'ops/nbo ',nid)

      call nekgsync
      if (nio.eq.0) write (6,*) 'dbas_time:',dnekclock()-dbas_time

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_cbas
      ! Calculates and dumps the POD basis of the convection snapshots
      !
      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      logical iftmp,iftmp2
      integer istep0
      parameter (ltd=lxd*lyd*lzd)
      common /romfine/ xfine(ltd,lelt),yfine(ltd,lelt),zfine(ltd,lelt),
     $                 ufine1(ltd,lelt),ufine2(ltd,lelt),
     $                 ufine3(ltd,lelt)

      integer pind(1)
      integer pmat(1,1) 

      ! Compute convection field for each snapshot and store in snapt
      call evalcflds(snapt,us0,us0,ldim,ns,.false.)

      iftmp=ifxyo
      iftmp2=ifpo

      ifpo=.false.
      ttmp=time
      istep0=istep

      ! Dump the convection snapshots if enabled
      if (ifdumpnls .and. .not.ifdumpfine) then
         do i=1,ns
            time=i
            istep=i
            ifxyo=(i.eq.1)
            call outpost(snapt(1,1,i),snapt(1,2,i),snapt(1,ldim,i),
     $                 pr,t,'csn')
         enddo
      endif

      call pod(uvwbnl,eval2,ug,snapt,ldim,ips,nbnl,ns,ifpb,
     $         'ops/guc  ',nbat)

      ! B-normalize the POD basis
      ! Could modify vnorm_ to handle this, but do it manually for now 
      do i=1,nbnl
        p=vip(uvwbnl(1,1,i),uvwbnl(1,2,i),uvwbnl(1,ldim,i),
     $        uvwbnl(1,1,i),uvwbnl(1,2,i),uvwbnl(1,ldim,i))
        s=1./sqrt(p)
        call opcmult(uvwbnl(1,1,i),uvwbnl(1,2,i),uvwbnl(1,ldim,i),s)
      enddo

      if (ifdeim) call dump_deim_inds

      if (ifdumpfine) then
         call intp_rstd_all(xfine,xm1,nelv)
         call intp_rstd_all(yfine,ym1,nelv)
         if (if3d) then
            call intp_rstd_all(zfine,zm1,nelv)
         else
            call rzero(zfine(1,1),ltd*nelv)
         endif

         if (ifdumpnls) then
            do i=1,ns
               time=i
               istep=i
               ifxyo=(i.eq.1)
               call convect_new(snapt(1,1,i),us0(1,1,i),.false.,
     $                          us0(1,1,i),us0(1,2,i),us0(1,ldim,i),
     $                          .false.)
               call copy(ufine1(1,1),ufine(1,1),ltd*nelv)

               call convect_new(snapt(1,2,i),us0(1,2,i),.false.,
     $                          us0(1,1,i),us0(1,2,i),us0(1,ldim,i),
     $                          .false.)
               call copy(ufine2(1,1),ufine(1,1),ltd*nelv)

               if (if3d) then
                  call convect_new(snapt(1,ldim,i),us0(1,ldim,i),
     $                             .false.,us0(1,1,i),us0(1,2,i),
     $                             us0(1,ldim,i),.false.)
                  call copy(ufine3(1,1),ufine(1,1),ltd*nelv)
               else
                  call rzero(ufine3(1,1),ltd*nelv)
               endif

               call dump_fine_vecfile('csn',xfine,yfine,zfine,
     $            ufine1,ufine2,ufine3,ifxyo)
            enddo
         endif
      endif

      ! Dump the convection basis
      if (ifdumpfine) then
         do i=1,nbnl
            time=i
            istep=i
            ifxyo=(i.eq.1)
            call intp_rstd_all(ufine1,uvwbnl(1,1,i),nelv)
            call intp_rstd_all(ufine2,uvwbnl(1,2,i),nelv)
            if (if3d) then
               call intp_rstd_all(ufine3,uvwbnl(1,ldim,i),nelv)
            else
               call rzero(ufine3(1,1),ltd*nelv)
            endif
            call dump_fine_vecfile('cba',xfine,yfine,zfine,
     $         ufine1,ufine2,ufine3,ifxyo)
         enddo
      else
         do i=1,nbnl
            time=i
            istep=i
            ! The temperature field isn't correct, but doesn't matter right now
            ! since temperature and pressure are not supported.
            ifxyo=(i.eq.1)
            call outpost2(uvwbnl(1,1,i),uvwbnl(1,2,i),uvwbnl(1,ldim,i),
     $                   pb(1,0),tb(1,0,1),ldimt,'cba')
         enddo
      endif

      istep=istep0
      time=ttmp

      ifxyo=iftmp
      ifpo=iftmp2

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_deim_inds
      ! Select a greedy DEIM point set from the nonlinear POD basis and
      ! emit the operator bundle consumed by the embedded DEIM runtime.
      ! The selection itself stays deterministic: each rank contributes a
      ! local residual maximizer in rank order, and the global maximizer is
      ! retained for the next greedy step.

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      integer lt
      parameter (lt=lx1*ly1*lz1*lelt)

      integer i,j,k,ip,comp,inode,info,nsel,isnap
      integer nnode_local,nstack_local,local_row,local_best_row
      integer local_best_comp
      integer global_row,offset,best_ip,best_comp,best_inode
      integer iuofs,ivofs,iwofs,iuxofs,iuyofs,iuzofs
      integer cand_info(3),work_info(3),ipivl(lbnl)
      integer nrowpack
      real    amat(lbnl,lbnl),rhs(lbnl),coeff(lbnl)
      real    sel_rows(lbnl,lbnl)
      real    loc_rowvals(lbnl),sel_rowvals(lbnl)
      real    cand_pack(lbnl+1),work_pack(lbnl+1)
      real    row_pack(6*(lub+1)),row_work(6*(lub+1))
      real    loc_best_val,glob_best_val,rowval
      real    snapcoef(lbnl),scale,rtmp(1)
      real    gx(lt),gy(lt),gz(lt)
      real    uadv(lt,ldim,1),tadv(lt,ldim,1)
      real    cf1(lt,ldim,1),cf2(lt,ldim,1)
      integer itmp(1)

      if (.not.ifdeim) return
      if (nbnl.le.0) return

      nsel = min(nbnl,lbnl)
      if (nsel.gt.ndeim_max) call exitti('ndeim_max too small$',nsel)
      if (if3d) then
         nrowpack = 6*(nb+1)
      else
         nrowpack = 4*(nb+1)
      endif
      if (nrowpack.gt.6*(lub+1)) call exitti('row pack too small$',nb)

      nnode_local = lx1*ly1*lz1*nelv
      nstack_local = ldim*nnode_local

      iuofs = 1
      ivofs = 1 + (nb+1)
      if (if3d) then
         iwofs = 1 + 2*(nb+1)
         iuxofs = 1 + 3*(nb+1)
         iuyofs = 1 + 4*(nb+1)
         iuzofs = 1 + 5*(nb+1)
      else
         iuxofs = 1 + 2*(nb+1)
         iuyofs = 1 + 3*(nb+1)
      endif

      call izero(deim_inds,ndeim_max)
      call izero(deim_inds_os,ndeim_max)
      call izero(deim_eval_inds,ndeim_max)
      call rzero(deim_eval_weights,ndeim_max)
      call rzero(deim_nl_bas_p_eval,ndeim_max*lbnl)
      call rzero(deim_u_p,ndeim_max*(lub+1))
      call rzero(deim_v_p,ndeim_max*(lub+1))
      call rzero(deim_w_p,ndeim_max*(lub+1))
      call rzero(deim_ux_p,ndeim_max*(lub+1))
      call rzero(deim_uy_p,ndeim_max*(lub+1))
      call rzero(deim_uz_p,ndeim_max*(lub+1))
      call rzero(deim_proj_mat,lub*lbnl)
      call rzero(deim_zmc,lub*(lub+1))
      call rzero(deim_Ainv,lbnl*lbnl)
      call rzero(deim_interp_mat,lbnl*ndeim_max)
      call rzero(sel_rows,lbnl*lbnl)
      call rzero(coeff,lbnl)
      call rzero(rhs,lbnl)
      call rzero(amat,lbnl*lbnl)
      call rzero(loc_rowvals,lbnl)
      call rzero(sel_rowvals,lbnl)

      do k=1,nsel
         best_ip = -1
         best_comp = 0
         best_inode = 0

         if (k.gt.1) then
            call rzero(amat,lbnl*lbnl)
            call rzero(rhs,lbnl)
            do i=1,k-1
               rhs(i) = sel_rows(k,i)
               do j=1,k-1
                  amat(i,j) = sel_rows(j,i)
               enddo
            enddo

            call izero(ipivl,lbnl)
            call dgetrf(k-1,k-1,amat,lbnl,ipivl,info)
            if (info.ne.0) call exitti(
     $           'DEIM selector factorization$',info)
            call dgetrs('N',k-1,1,amat,lbnl,ipivl,rhs,lbnl,info)
            if (info.ne.0) call exitti(
     $           'DEIM selector solve$',info)
            call rzero(coeff,lbnl)
            do i=1,k-1
               coeff(i) = rhs(i)
            enddo
         endif

         loc_best_val = -1.0e30
         local_best_row = 0
         local_best_comp = 0
         call rzero(loc_rowvals,lbnl)
         local_row = 0

         do comp=1,ldim
            do inode=1,nnode_local
               local_row = local_row + 1
               rowval = uvwbnl(inode,comp,k)
               do j=1,k-1
                  rowval = rowval - uvwbnl(inode,comp,j)*coeff(j)
               enddo
               rowval = abs(rowval)

               if (rowval.gt.loc_best_val) then
                  loc_best_val = rowval
                  local_best_row = local_row
                  local_best_comp = comp
                  do i=1,nsel
                     loc_rowvals(i) = uvwbnl(inode,comp,i)
                  enddo
               endif
            enddo
         enddo

         cand_pack(1) = loc_best_val
         do i=1,nsel
            cand_pack(i+1) = loc_rowvals(i)
         enddo
         cand_info(1) = local_best_row
         cand_info(2) = nstack_local
         cand_info(3) = local_best_comp

         glob_best_val = -1.0e30
         global_row = 0
         offset = 0

         do ip=0,np-1
            if (nid.ne.ip) then
               call rzero(cand_pack,nsel+1)
               call izero(cand_info,3)
            endif

            call gop(cand_pack,work_pack,'+  ',nsel+1)
            call igop(cand_info,work_info,'+  ',3)

            if (cand_pack(1).gt.glob_best_val) then
               glob_best_val = cand_pack(1)
               global_row = offset + cand_info(1)
               best_ip = ip
               best_comp = cand_info(3)
               best_inode = cand_info(1) - (cand_info(3)-1)*nnode_local
               do i=1,nsel
                  sel_rows(i,k) = cand_pack(i+1)
               enddo
            endif

            offset = offset + cand_info(2)
         enddo

         deim_inds(k) = global_row
         deim_eval_inds(k) = global_row
         deim_eval_weights(k) = 1.0
         do i=1,nsel
            deim_nl_bas_p_eval(k,i) = sel_rows(i,k)
         enddo

         call rzero(row_pack,6*(lub+1))
         if (nid.eq.best_ip) then
            do j=0,nb
               row_pack(iuofs+j) = ub(best_inode,j)
               row_pack(ivofs+j) = vb(best_inode,j)
               if (if3d) then
                  row_pack(iwofs+j) = wb(best_inode,j)
               endif
            enddo

            if (best_comp.eq.1) then
               do j=0,nb
                  call gradm1(gx,gy,gz,ub(1,j))
                  row_pack(iuxofs+j) = gx(best_inode)
                  row_pack(iuyofs+j) = gy(best_inode)
                  if (if3d) then
                     row_pack(iuzofs+j) = gz(best_inode)
                  endif
               enddo
            else if (best_comp.eq.2) then
               do j=0,nb
                  call gradm1(gx,gy,gz,vb(1,j))
                  row_pack(iuxofs+j) = gx(best_inode)
                  row_pack(iuyofs+j) = gy(best_inode)
                  if (if3d) then
                     row_pack(iuzofs+j) = gz(best_inode)
                  endif
               enddo
            else
               do j=0,nb
                  call gradm1(gx,gy,gz,wb(1,j))
                  row_pack(iuxofs+j) = gx(best_inode)
                  row_pack(iuyofs+j) = gy(best_inode)
                  if (if3d) then
                     row_pack(iuzofs+j) = gz(best_inode)
                  endif
               enddo
            endif
         endif

         call gop(row_pack,row_work,'+  ',nrowpack)

         do j=0,nb
            deim_u_p(k,j) = row_pack(iuofs+j)
            deim_v_p(k,j) = row_pack(ivofs+j)
            if (if3d) then
               deim_w_p(k,j) = row_pack(iwofs+j)
            endif
            deim_ux_p(k,j) = row_pack(iuxofs+j)
            deim_uy_p(k,j) = row_pack(iuyofs+j)
            if (if3d) then
               deim_uz_p(k,j) = row_pack(iuzofs+j)
            endif
         enddo

         if (nid.eq.0) write (6,*) 'DEIM selector step',k,'row',global_row
      enddo

      ndeim_pts = nsel
      ndeim_pts_eval = nsel
      ndeim_pts_os = 0

      call rzero(amat,lbnl*lbnl)
      do i=1,nsel
         do j=1,nsel
            do k=1,nsel
               amat(i,j) = amat(i,j) + sel_rows(i,k)*sel_rows(j,k)
            enddo
         enddo
      enddo

      call rzero(deim_Ainv,lbnl*lbnl)
      do i=1,nsel
         deim_Ainv(i,i) = 1.0
      enddo

      call izero(ipivl,lbnl)
      call dgetrf(nsel,nsel,amat,lbnl,ipivl,info)
      if (info.ne.0) call exitti('DEIM inverse factorization$',info)
      call dgetrs('N',nsel,nsel,amat,lbnl,ipivl,deim_Ainv,lbnl,info)
      if (info.ne.0) call exitti('DEIM inverse solve$',info)

      do i=1,nsel
         do j=1,nsel
            deim_interp_mat(i,j) = 0.0
            do k=1,nsel
               deim_interp_mat(i,j) =
     $            deim_interp_mat(i,j) + deim_Ainv(i,k)*sel_rows(k,j)
            enddo
         enddo
      enddo

      call rzero(deim_proj_mat,lub*lbnl)
      do i=1,nb
         do j=1,nsel
            deim_proj_mat(i,j) = op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),
     $         uvwbnl(1,1,j),uvwbnl(1,2,j),uvwbnl(1,ldim,j),bm1)
         enddo
      enddo

      call rzero(cf1,lt*ldim)
      call rzero(cf2,lt*ldim)
      do j=0,nb
         call opcopy(uadv(1,1,1),uadv(1,2,1),uadv(1,ldim,1),
     $      ub(1,j),vb(1,j),wb(1,j))
         call opcopy(tadv(1,1,1),tadv(1,2,1),tadv(1,ldim,1),
     $      ub(1,0),vb(1,0),wb(1,0))
         call evalcflds(cf1,uadv,tadv,ldim,1,.false.)

         call opcopy(uadv(1,1,1),uadv(1,2,1),uadv(1,ldim,1),
     $      ub(1,0),vb(1,0),wb(1,0))
         call opcopy(tadv(1,1,1),tadv(1,2,1),tadv(1,ldim,1),
     $      ub(1,j),vb(1,j),wb(1,j))
         call evalcflds(cf2,uadv,tadv,ldim,1,.false.)

         do i=1,ldim
            call add2(cf1(1,i,1),cf2(1,i,1),nnode_local)
         enddo

         do i=1,nb
            deim_zmc(i,j) = op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),
     $         cf1(1,1,1),cf1(1,2,1),cf1(1,ldim,1),bm1)
         enddo
      enddo

      do i=1,nb
         deim_zmc(i,0) = 0.5*deim_zmc(i,0)
      enddo

      call rzero(deim_eval_weights,ndeim_max)
      do i=1,nsel
         deim_eval_weights(i) = 1.0
      enddo

      itmp(1) = ndeim_pts
      call idump_serial(itmp,1,'ops/deim_npts ',nid)
      itmp(1) = ndeim_pts_os
      call idump_serial(itmp,1,'ops/deim_npts_os ',nid)
      itmp(1) = ndeim_pts_eval
      call idump_serial(itmp,1,'ops/deim_npts_eval ',nid)
      call idump_serial(deim_inds,nsel,'ops/deim_inds ',nid)
      call idump_serial(deim_eval_inds,nsel,'ops/deim_eval_inds ',nid)
      call dump_serial(deim_eval_weights,nsel,'ops/deim_eval_weights ',
     $   nid)
      call dump_serial(deim_u_p,nsel*(nb+1),'ops/deim_u_p ',nid)
      call dump_serial(deim_v_p,nsel*(nb+1),'ops/deim_v_p ',nid)
      if (if3d) then
         call dump_serial(deim_w_p,nsel*(nb+1),'ops/deim_w_p ',nid)
      endif
      call dump_serial(deim_ux_p,nsel*(nb+1),'ops/deim_ux_p ',nid)
      call dump_serial(deim_uy_p,nsel*(nb+1),'ops/deim_uy_p ',nid)
      if (if3d) then
         call dump_serial(deim_uz_p,nsel*(nb+1),'ops/deim_uz_p ',nid)
      endif
      call dump_serial(deim_nl_bas_p_eval,nsel*nsel,
     $   'ops/deim_nl_bas_p_eval ',nid)
      call dump_serial(deim_proj_mat,nb*nsel,'ops/deim_proj_mat ',nid)
      call dump_serial(deim_zmc,nb*(nb+1),'ops/deim_zmc ',nid)
      call dump_serial(deim_Ainv,nsel*nsel,'ops/deim_Ainv ',nid)
      call dump_serial(deim_interp_mat,nsel*nsel,
     $   'ops/deim_interp_mat ',nid)

      if (deimmode.eq.'MCLSDEIM') then
         call rzero(deim_mu,lbnl)
         call rzero(deim_tau,lbnl*lbnl)
         call rzero(deim_A_tau_inv,lbnl*lbnl)
         call rzero(snapcoef,lbnl)

         do isnap=1,ns
            do i=1,nsel
               snapcoef(i)=0.
               do comp=1,ldim
                  snapcoef(i)=snapcoef(i)+glsc2(
     $               uvwbnl(1,comp,i),snapt(1,comp,isnap),nnode_local)
               enddo
               deim_mu(i)=deim_mu(i)+snapcoef(i)
            enddo
         enddo

         if (ns.gt.0) then
            scale=1./real(ns)
            call cmult(deim_mu,scale,nsel)
         endif

         call rzero(amat,lbnl*lbnl)
         do isnap=1,ns
            do i=1,nsel
               snapcoef(i)=0.
               do comp=1,ldim
                  snapcoef(i)=snapcoef(i)+glsc2(
     $               uvwbnl(1,comp,i),snapt(1,comp,isnap),nnode_local)
               enddo
               snapcoef(i)=snapcoef(i)-deim_mu(i)
            enddo

            do i=1,nsel
            do j=1,nsel
               amat(i,j)=amat(i,j)+snapcoef(i)*snapcoef(j)
            enddo
            enddo
         enddo

         if (ns.gt.1) then
            scale=1./real(ns-1)
            call cmult(amat,scale,nsel*nsel)
         else
            call rzero(amat,nsel*nsel)
         endif

         do i=1,nsel
            amat(i,i)=amat(i,i)+1.e-15
         enddo

         call izero(ipivl,lbnl)
         call dgetrf(nsel,nsel,amat,lbnl,ipivl,info)
         if (info.ne.0) call exitti(
     $        'MCLSDEIM tau factorization$',info)
         call rzero(deim_tau,lbnl*lbnl)
         do i=1,nsel
            deim_tau(i,i)=1.
         enddo
         call dgetrs('N',nsel,nsel,amat,lbnl,ipivl,deim_tau,lbnl,info)
         if (info.ne.0) call exitti('MCLSDEIM tau solve$',info)

         call rzero(amat,lbnl*lbnl)
         do i=1,nsel
         do j=1,nsel
            do k=1,nsel
               amat(i,j)=amat(i,j)+deim_nl_bas_p_eval(k,i)
     $            *deim_nl_bas_p_eval(k,j)
            enddo
            amat(i,j)=amat(i,j)+deim_alpha*deim_tau(i,j)
         enddo
         enddo

         call izero(ipivl,lbnl)
         call dgetrf(nsel,nsel,amat,lbnl,ipivl,info)
         if (info.ne.0) call exitti(
     $        'MCLSDEIM A_tau factorization$',info)
         call rzero(deim_A_tau_inv,lbnl*lbnl)
         do i=1,nsel
            deim_A_tau_inv(i,i)=1.
         enddo
         call dgetrs('N',nsel,nsel,amat,lbnl,ipivl,deim_A_tau_inv,
     $      lbnl,info)
         if (info.ne.0) call exitti('MCLSDEIM A_tau solve$',info)

         rtmp(1)=deim_alpha
         call dump_serial(rtmp,1,'ops/deim_alpha ',nid)
         call dump_serial(deim_mu,nsel,'ops/deim_mu ',nid)
         call dump_serial(deim_tau,nsel*nsel,'ops/deim_tau ',nid)
         call dump_serial(deim_A_tau_inv,nsel*nsel,
     $      'ops/deim_A_tau_inv ',nid)
      endif

      if (nid.eq.0) write (6,*) 'DEIM selector complete, points:',nsel

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_misc
      ! Dump miscellaneous items
      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /dumpglobal/ wk1(lcloc),wk2(lcloc)

      call nekgsync
      dmisc_time=dnekclock()

      if (ifrom(1)) then
         call dump_serial(u,(nb+1)*3,'ops/u ',nid)
         call dump_serial(uk,ns*(nb+1),'ops/uk ',nid)
         call dump_serial(umin,nb,'ops/umin ',nid)
         call dump_serial(umax,nb,'ops/umax ',nid)
         call dump_serial(timek,ns,'ops/timek ',nid)

         if (ifcdrag) then
            call dump_serial(rdgx,nb+1,'qoi/rdgx ',nid)
            call dump_serial(rdgy,nb+1,'qoi/rdgy ',nid)
            if (ldim.eq.3) call dump_serial(rdgz,nb+1,'qoi/rdgz ',nid)

            call dump_serial(fd1,ldim*(nb+1),'qoi/fd1 ',nid)
            call dump_serial(fd2,ldim*(nb+1)**2,'qoi/fd2 ',nid)
            call dump_serial(fd3,ldim*(nb+1),'qoi/fd3 ',nid)
         endif
      endif

      if (ifrom(2)) then
         call dump_serial(ut,(nb+1)*3,'ops/t ',nid)
         call dump_serial(tk,ns*(nb+1),'ops/tk ',nid)
         call dump_serial(tmin,nb,'ops/tmin ',nid)
         call dump_serial(tmax,nb,'ops/tmax ',nid)
         if (.not.ifpod(1))
     $      call dump_serial(timek,ns,'ops/timek ',nid)
      endif

      if (ifforce)  call dump_serial(rf,nb,'ops/rf ',nid)
      if (ifsource) call dump_serial(rq,nb,'ops/rq ',nid)
      if (ifbuoy)   call dump_serial(but0,(nb+1)**2,'ops/but ',nid)
      if (ifedvs)   call dump_serial(edk,ns*(nb+1),'ops/edk ',nid)

      if (ifedvs) then
         call dump_serial(rbfwt,ns*nb,'ops/rbfwt ',nid)
      endif

      if (ifei) then
         l=1
         do j=1,nres
         do i=1,nres
            sigtmp(l,1)=mor_sigma(i,j)
            l=l+1
         enddo
         enddo
         call dump_serial(sigtmp,nres*nres,'ops/sigma ',nid)
      endif

      call nekgsync
      if (nio.eq.0) write (6,*) 'dmisc_time:',dnekclock()-dmisc_time

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_snaps

      ! dump velocity and temperature snapshots

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      logical iftmp,iftmp2

      iftmp=ifxyo
      iftmp2=ifpo

      ifxyo=.true.
      ifpo=.false.

      do i=1,ns
         call outpost(us0(1,1,i),us0(1,2,i),us0(1,ldim,i),
     $      pr,ts0(1,i,1),'sna')
         ifxyo=.false.
      enddo

      ifxyo=iftmp
      ifpo=iftmp2

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_fine_vecfile(prefix,x,y,z,u,v,w,ifcoord)

      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'

      character*3 prefix
      logical ifcoord

      real x(1),y(1),z(1),u(1),v(1),w(1)

      integer ierr,nout,nxo0,nyo0,nzo0

      call blank(rdcode1,10)
      if (ifcoord) then
         rdcode1(1)='X'
         rdcode1(2)='U'
      else
         rdcode1(1)='U'
      endif

      nxo0=nxo
      nyo0=nyo
      nzo0=nzo

      nxo=lxd
      nyo=lyd
      nzo=lzd

      ierr=0
      if (nid.eq.pid0) call mfo_open_files(prefix,ierr)
      call err_chk(ierr,'Error opening file in dump_fine_vecfile. $')

      call mfo_write_hdr(rdcode1)

      nout=nelt
      if (ifcoord) then
         call mfo_outv_fine(x,y,z,nout,lxd,lyd,lzd)
      endif
      call mfo_outv_fine(u,v,w,nout,lxd,lyd,lzd)

      if (nid.eq.pid0) then
      if (ifmpiio) then
            call byte_close_mpi(ifh_mbyte,ierr)
         else
            call byte_close(ierr)
         endif
      endif
      call err_chk(ierr,'Error closing file in dump_fine_vecfile. $')

      nxo=nxo0
      nyo=nyo0
      nzo=nzo0

      return
      end
c-----------------------------------------------------------------------
      subroutine mfo_outv_fine(u,v,w,nel,mx,my,mz)   ! output a vector field

      include 'SIZE'
      include 'INPUT'
      include 'RESTART'

      real u(mx*my*mz,1),v(mx*my*mz,1),w(mx*my*mz,1)

      common /SCRNF/ u4(2+lxd*lxd*lzd*6*lelt)
      real*4         u4
      real*8         u8(1+lxd*lxd*lzd*3*lelt)
      equivalence    (u4,u8)

      integer e
      integer cnt

      call nekgsync() ! clear outstanding message queues.
      if(mx.gt.lxd .or. my.gt.lyd .or. mz.gt.lzd) then
        if(nid.eq.0) write(6,*) 'ABORT: fine output buffer too small'
        call exitt
      endif

      nxyz  = mx*my*mz
      lrecv = 8 + 8*(lelt*nxyz*ldim)   ! recv buffer size (u4)
      lsend = 8 + wdsizo*(nel*nxyz*ldim)
      idum  = 1
      ierr  = 0

      if (nid.eq.pid0) then
         cnt = 0
         j = 0
         if (wdsizo.eq.4) then             ! 32-bit output
             do iel = 1,nel
               if(out_mask(iel).ne.0) then
                 call copyx4   (u4(j+1),u(1,iel),nxyz)
                 j = j + nxyz
                 call copyx4   (u4(j+1),v(1,iel),nxyz)
                 j = j + nxyz
                 if(if3d) then
                   call copyx4 (u4(j+1),w(1,iel),nxyz)
                   j = j + nxyz
                 endif
                 cnt = cnt + 1
               endif
             enddo
         else
             do iel = 1,nel
               if(out_mask(iel).ne.0) then
                 call copy     (u8(j+1),u(1,iel),nxyz)
                 j = j + nxyz
                 call copy     (u8(j+1),v(1,iel),nxyz)
                 j = j + nxyz
                 if(if3d) then
                   call copy   (u8(j+1),w(1,iel),nxyz)
                   j = j + nxyz
                 endif
                 cnt = cnt + 1
               endif
             enddo
         endif
         nout = wdsizo/4 * ldim*cnt * nxyz
         if(ierr.eq.0) then
           if(ifmpiio) then
             call byte_write_mpi(u4,nout,-1,ifh_mbyte,ierr)
           else
             call byte_write(u4,nout,ierr)          ! u4 :=: u8
           endif
         endif

         ! write out the data of my childs
         do k=pid0+1,pid1
            mtype = k
            call csend(mtype,idum,4,k,0)           ! handshake
            call crecv(mtype,u4,lrecv)
            nout = wdsizo/4 * ldim*nxyz * u8(1)

            if (wdsizo.eq.4.and.ierr.eq.0) then
               if(ifmpiio) then
                 call byte_write_mpi(u4(3),nout,-1,ifh_mbyte,ierr)
               else
                 call byte_write(u4(3),nout,ierr)
               endif
            elseif(ierr.eq.0) then
               if(ifmpiio) then
                 call byte_write_mpi(u8(2),nout,-1,ifh_mbyte,ierr)
               else
                 call byte_write(u8(2),nout,ierr)
               endif
            endif
         enddo
      else
         cnt = 0
         if (wdsizo.eq.4) then             ! 32-bit output
             j = 2
             do iel = 1,nel
               if(out_mask(iel).ne.0) then
                 call copyx4   (u4(j+1),u(1,iel),nxyz)
                 j = j + nxyz
                 call copyx4   (u4(j+1),v(1,iel),nxyz)
                 j = j + nxyz
                 if(if3d) then
                   call copyx4 (u4(j+1),w(1,iel),nxyz)
                   j = j + nxyz
                 endif
                 cnt = cnt + 1
               endif
             enddo
         else
             j = 1
             do iel = 1,nel
               if(out_mask(iel).ne.0) then
                 call copy     (u8(j+1),u(1,iel),nxyz)
                 j = j + nxyz
                 call copy     (u8(j+1),v(1,iel),nxyz)
                 j = j + nxyz
                 if(if3d) then
                   call copy   (u8(j+1),w(1,iel),nxyz)
                   j = j + nxyz
                 endif
                 cnt = cnt + 1
               endif
             enddo
         endif
         u8(1) = cnt

         mtype = nid
         call crecv(mtype,idum,4)            ! hand-shake
         call csend(mtype,u4,lsend,pid0,0)     ! u4 :=: u8
      endif

      call err_chk(ierr,'Error writing data in mfo_outv_fine. $')
      return
      end
c-----------------------------------------------------------------------
