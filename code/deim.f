c-----------------------------------------------------------------------
      subroutine dump_deim_inds_impl
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
      integer cand_info(3),work_info(3),ipivl(lbnl_eff)
      integer irl(lbnl_eff),icl(lbnl_eff)
      integer nrowpack
      real    amat(lbnl_eff,lbnl_eff),rhs(lbnl_eff),coeff(lbnl_eff)
      real    sel_rows(lbnl_eff,lbnl_eff)
      real    loc_rowvals(lbnl_eff),sel_rowvals(lbnl_eff)
      real    cand_pack(lbnl_eff+1),work_pack(lbnl_eff+1)
      real    row_pack(6*(lub+1)),row_work(6*(lub+1))
      real    loc_best_val,glob_best_val,rowval
      real    snapcoef(lbnl_eff),scale
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
      call rzero(deim_nl_bas_p_eval,ndeim_max*lbnl_eff)
      call rzero(deim_u_p,ndeim_max*(lub+1))
      call rzero(deim_v_p,ndeim_max*(lub+1))
      call rzero(deim_w_p,ndeim_max*(lub+1))
      call rzero(deim_ux_p,ndeim_max*(lub+1))
      call rzero(deim_uy_p,ndeim_max*(lub+1))
      call rzero(deim_uz_p,ndeim_max*(lub+1))
      call rzero(deim_proj_mat,lub*lbnl_eff)
      call rzero(deim_zmc,lub*(lub+1))
      call rzero(deim_Ainv,lbnl_eff*lbnl_eff)
      call rzero(deim_interp_mat,lbnl_eff*ndeim_max)
      call rzero(sel_rows,lbnl_eff*lbnl_eff)
      call rzero(coeff,lbnl_eff)
      call rzero(rhs,lbnl_eff)
      call rzero(amat,lbnl_eff*lbnl_eff)
      call rzero(loc_rowvals,lbnl_eff)
      call rzero(sel_rowvals,lbnl_eff)

      do k=1,nsel
         best_ip = -1
         best_comp = 0
         best_inode = 0

         if (k.gt.1) then
            call rzero(amat,lbnl_eff*lbnl_eff)
            call rzero(rhs,lbnl_eff)
            do i=1,k-1
               rhs(i) = sel_rows(k,i)
               do j=1,k-1
                  amat(i,j) = sel_rows(j,i)
               enddo
            enddo

            call izero(ipivl,lbnl_eff)
            call dgetrf(k-1,k-1,amat,lbnl_eff,ipivl,info)
            if (info.ne.0) call exitti(
     $           'DEIM selector factorization$',info)
            call dgetrs('N',k-1,1,amat,lbnl_eff,ipivl,rhs,lbnl_eff,info)
            if (info.ne.0) call exitti(
     $           'DEIM selector solve$',info)
            call rzero(coeff,lbnl_eff)
            do i=1,k-1
               coeff(i) = rhs(i)
            enddo
         endif

         loc_best_val = -1.0e30
         local_best_row = 0
         local_best_comp = 0
         call rzero(loc_rowvals,lbnl_eff)
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

      call rzero(amat,lbnl_eff*lbnl_eff)
      do i=1,nsel
         do j=1,nsel
            do k=1,nsel
               amat(i,j) = amat(i,j) + sel_rows(i,k)*sel_rows(j,k)
            enddo
         enddo
      enddo

      call rzero(deim_Ainv,lbnl_eff*lbnl_eff)
      do i=1,nsel
         deim_Ainv(i,i) = 1.0
      enddo

      call izero(ipivl,lbnl_eff)
      call dgetrf(nsel,nsel,amat,lbnl_eff,ipivl,info)
      if (info.ne.0) call exitti('DEIM inverse factorization$',info)
      call dgetrs('N',nsel,nsel,amat,lbnl_eff,ipivl,deim_Ainv,lbnl_eff,
     $   info)
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

      call rzero(deim_proj_mat,lub*lbnl_eff)
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
      call dump_mat_serial(deim_u_p,ndeim_max,lub+1,
     $   'ops/deim_u_p ',nsel,nb+1,nid)
      call dump_mat_serial(deim_v_p,ndeim_max,lub+1,
     $   'ops/deim_v_p ',nsel,nb+1,nid)
      if (if3d) then
         call dump_mat_serial(deim_w_p,ndeim_max,lub+1,
     $      'ops/deim_w_p ',nsel,nb+1,nid)
      endif
      call dump_mat_serial(deim_ux_p,ndeim_max,lub+1,
     $   'ops/deim_ux_p ',nsel,nb+1,nid)
      call dump_mat_serial(deim_uy_p,ndeim_max,lub+1,
     $   'ops/deim_uy_p ',nsel,nb+1,nid)
      if (if3d) then
         call dump_mat_serial(deim_uz_p,ndeim_max,lub+1,
     $      'ops/deim_uz_p ',nsel,nb+1,nid)
      endif
      call dump_mat_serial(deim_nl_bas_p_eval,ndeim_max,lbnl_eff,
     $   'ops/deim_nl_bas_p_eval ',nsel,nsel,nid)
      call dump_mat_serial(deim_proj_mat,lub,lbnl_eff,
     $   'ops/deim_proj_mat ',nb,nsel,nid)
      call dump_mat_serial(deim_zmc,lub,lub+1,'ops/deim_zmc ',
     $   nb,nb+1,nid)
      call dump_mat_serial(deim_Ainv,lbnl_eff,lbnl_eff,
     $   'ops/deim_Ainv ',nsel,nsel,nid)
      call dump_mat_serial(deim_interp_mat,lbnl_eff,ndeim_max,
     $   'ops/deim_interp_mat ',nsel,nsel,nid)

      if (deimmode.eq.'MCLSDEIM') then
         call rzero(deim_mu,lbnl_eff)
         call rzero(deim_tau,lbnl_eff*lbnl_eff)
         call rzero(deim_A_tau_inv,lbnl_eff*lbnl_eff)
         call rzero(snapcoef,lbnl_eff)

         do isnap=1,ns
            do i=1,nsel
               snapcoef(i)=0.
               snapcoef(i)=op_glsc2_wt(uvwbnl(1,1,i),uvwbnl(1,2,i),
     $            uvwbnl(1,ldim,i),snapt(1,1,isnap),snapt(1,2,isnap),
     $            snapt(1,ldim,isnap),bm1)
               deim_mu(i)=deim_mu(i)+snapcoef(i)
            enddo
         enddo

         if (ns.gt.0) then
            scale=1./real(ns)
            call cmult(deim_mu,scale,nsel)
         endif

         call rzero(amat,lbnl_eff*lbnl_eff)
         do isnap=1,ns
            do i=1,nsel
               snapcoef(i)=0.
               snapcoef(i)=op_glsc2_wt(uvwbnl(1,1,i),uvwbnl(1,2,i),
     $            uvwbnl(1,ldim,i),snapt(1,1,isnap),snapt(1,2,isnap),
     $            snapt(1,ldim,isnap),bm1)
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
            call cmult(amat,scale,lbnl_eff*lbnl_eff)
         else
            call rzero(amat,lbnl_eff*lbnl_eff)
         endif

         do i=1,nsel
            amat(i,i)=amat(i,i)+1.e-15
         enddo

         call izero(irl,lbnl_eff)
         call izero(icl,lbnl_eff)
         call lu(amat,nsel,lbnl_eff,irl,icl)
         call rzero(deim_tau,lbnl_eff*lbnl_eff)
         do i=1,nsel
            deim_tau(i,i)=1.
         enddo
         call solve(deim_tau,amat,nsel,nsel,lbnl_eff,irl,icl)

         call rzero(amat,lbnl_eff*lbnl_eff)
         do i=1,nsel
         do j=1,nsel
            do k=1,nsel
               amat(i,j)=amat(i,j)+deim_nl_bas_p_eval(k,i)
     $            *deim_nl_bas_p_eval(k,j)
            enddo
            amat(i,j)=amat(i,j)+deim_alpha*deim_tau(i,j)
         enddo
         enddo

         call izero(irl,lbnl_eff)
         call izero(icl,lbnl_eff)
         call lu(amat,nsel,lbnl_eff,irl,icl)
         call rzero(deim_A_tau_inv,lbnl_eff*lbnl_eff)
         do i=1,nsel
            deim_A_tau_inv(i,i)=1.
         enddo
         call solve(deim_A_tau_inv,amat,nsel,nsel,lbnl_eff,irl,icl)

         call dump_serial(deim_mu,nsel,'ops/deim_mu ',nid)
         call dump_mat_serial(deim_tau,lbnl_eff,lbnl_eff,
     $      'ops/deim_tau ',nsel,nsel,nid)
         call dump_mat_serial(deim_A_tau_inv,lbnl_eff,lbnl_eff,
     $      'ops/deim_A_tau_inv ',nsel,nsel,nid)
         call dump_serial(deim_alpha,1,'ops/deim_alpha ',nid)
      endif

      if (nid.eq.0) write (6,*) 'DEIM selector complete, points:',nsel

      return
      end
c-----------------------------------------------------------------------
      subroutine setdeim_impl

      ! load DEIM-family artifacts

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      logical have_alpha,have_mu,have_tau,have_ainv,
     $        have_w_p,have_uz_p

      integer ndeimwrk
      parameter (ndeimwrk=ndeim_max*(lub+1+lbnl_eff))

      real rwk(ndeimwrk)
      integer iwk(ndeimwrk)
      integer itmp(1)

      if (nio.eq.0) write (6,*) 'inside setdeim'

      if (.not.ifdeim) return
      if (nbnl.le.0) call exitti('nbnl <= 0$',nbnl)

      call iread_serial(itmp,1,'ops/deim_npts ',iwk,nid)
      ndeim_pts=itmp(1)
      call iread_serial(itmp,1,'ops/deim_npts_os ',iwk,nid)
      ndeim_pts_os=itmp(1)
      call iread_serial(itmp,1,'ops/deim_npts_eval ',iwk,nid)
      ndeim_pts_eval=itmp(1)

      if (ndeim_pts.le.0) call exitti('ndeim_pts <= 0$',ndeim_pts)
      if (ndeim_pts.gt.ndeim_max) then
         call exitti('ndeim_pts > ldeim$',ndeim_pts)
      endif
      if (ndeim_pts_os.gt.ndeim_max) then
         call exitti('ndeim_pts_os > ldeim$',ndeim_pts_os)
      endif
      if (ndeim_pts_eval.le.0) call exitti(
     $   'ndeim_pts_eval <= 0$',ndeim_pts_eval)
      if (ndeim_pts_eval.gt.ndeim_max) then
         call exitti('ndeim_pts_eval > ldeim$',ndeim_pts_eval)
      endif

      call iread_serial(deim_inds,ndeim_pts,'ops/deim_inds ',iwk,nid)
      if (ndeim_pts_os.gt.0) then
         call iread_serial(deim_inds_os,ndeim_pts_os,
     $      'ops/deim_inds_os ',iwk,nid)
      endif
      call iread_serial(deim_eval_inds,ndeim_pts_eval,
     $   'ops/deim_eval_inds ',iwk,nid)

      call read_serial(deim_eval_weights,ndeim_pts_eval,
     $   'ops/deim_eval_weights ',rwk,nid)
      call read_mat_serial(deim_u_p,ndeim_max,lub+1,
     $   'ops/deim_u_p ',ndeim_pts_eval,nb+1,rwk,nid)
      call read_mat_serial(deim_v_p,ndeim_max,lub+1,
     $   'ops/deim_v_p ',ndeim_pts_eval,nb+1,rwk,nid)
      if (if3d) then
         inquire (file='ops/deim_w_p',exist=have_w_p)
         inquire (file='ops/deim_uz_p',exist=have_uz_p)
         if ((.not.have_w_p).or.(.not.have_uz_p)) then
            call exitti('missing 3D DEIM artifacts$',1)
         endif
         call read_mat_serial(deim_w_p,ndeim_max,lub+1,
     $      'ops/deim_w_p ',ndeim_pts_eval,nb+1,rwk,nid)
      else
         call rzero(deim_w_p,ndeim_max*(lub+1))
      endif
      call read_mat_serial(deim_ux_p,ndeim_max,lub+1,
     $   'ops/deim_ux_p ',ndeim_pts_eval,nb+1,rwk,nid)
      call read_mat_serial(deim_uy_p,ndeim_max,lub+1,
     $   'ops/deim_uy_p ',ndeim_pts_eval,nb+1,rwk,nid)
      if (if3d) then
         call read_mat_serial(deim_uz_p,ndeim_max,lub+1,
     $      'ops/deim_uz_p ',ndeim_pts_eval,nb+1,rwk,nid)
      else
         call rzero(deim_uz_p,ndeim_max*(lub+1))
      endif
      call read_mat_serial(deim_nl_bas_p_eval,ndeim_max,lbnl_eff,
     $   'ops/deim_nl_bas_p_eval ',ndeim_pts_eval,nbnl,rwk,nid)
      call read_mat_serial(deim_proj_mat,lub,lbnl_eff,
     $   'ops/deim_proj_mat ',nb,nbnl,rwk,nid)
      call read_mat_serial(deim_zmc,lub,lub+1,'ops/deim_zmc ',
     $   nb,nb+1,rwk,nid)
      call read_mat_serial(deim_Ainv,lbnl_eff,lbnl_eff,
     $   'ops/deim_Ainv ',nbnl,nbnl,rwk,nid)
      call read_mat_serial(deim_interp_mat,lbnl_eff,ndeim_max,
     $   'ops/deim_interp_mat ',nbnl,ndeim_pts_eval,rwk,nid)

      inquire (file='ops/deim_mu',exist=have_mu)
      if (have_mu) then
         call read_serial(deim_mu,nbnl,'ops/deim_mu ',rwk,nid)
      else
         call rzero(deim_mu,nbnl)
      endif

      inquire (file='ops/deim_tau',exist=have_tau)
      if (have_tau) then
         call read_mat_serial(deim_tau,lbnl_eff,lbnl_eff,
     $      'ops/deim_tau ',nbnl,nbnl,rwk,nid)
      else
         call rzero(deim_tau,nbnl*nbnl)
      endif

      inquire (file='ops/deim_A_tau_inv',exist=have_ainv)
      if (have_ainv) then
         call read_mat_serial(deim_A_tau_inv,lbnl_eff,lbnl_eff,
     $      'ops/deim_A_tau_inv ',
     $      nbnl,nbnl,rwk,nid)
      else
         call rzero(deim_A_tau_inv,nbnl*nbnl)
      endif

      inquire (file='ops/deim_alpha',exist=have_alpha)
      if (have_alpha) then
         call read_serial(deim_alpha,1,'ops/deim_alpha ',rwk,nid)
      endif

      if (deimmode.eq.'MCLSDEIM') then
         if ((.not.have_mu).or.(.not.have_tau).or.(.not.have_ainv)) then
            call exitti('missing MCLSDEIM artifacts$',1)
         endif
      endif

      call nekgsync
      if (nio.eq.0) write (6,*) 'deim setup complete'

      return
      end
c-----------------------------------------------------------------------
      subroutine evalc_deim_impl(cu,uu)

      ! Compute the DEIM-family convection term for the velocity equation.

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real cu(nb)
      real uu(0:nb)

      real fraw(ndeim_max),c_hat(lbnl_eff),bvec(lbnl_eff),tmpb(lbnl_eff)
      real up,vp,wp,ux,uy,uz,lambda,s1,s2

      integer i,j

      call rzero(cu,nb)
      call rzero(fraw,ndeim_max)
      call rzero(c_hat,lbnl_eff)
      call rzero(bvec,lbnl_eff)
      call rzero(tmpb,lbnl_eff)

      if (.not.ifdeim) return

      if (ndeim_pts_eval.le.0.or.nbnl.le.0) then
         call exitti('invalid DEIM setup$',ndeim_pts_eval)
      endif

      do i=1,ndeim_pts_eval
         up=0.
         vp=0.
         wp=0.
         ux=0.
         uy=0.
         uz=0.
         do j=1,nb
            up=up+deim_u_p(i,j)*uu(j)
            vp=vp+deim_v_p(i,j)*uu(j)
            if (if3d) then
               wp=wp+deim_w_p(i,j)*uu(j)
            endif
            ux=ux+deim_ux_p(i,j)*uu(j)
            uy=uy+deim_uy_p(i,j)*uu(j)
            if (if3d) then
               uz=uz+deim_uz_p(i,j)*uu(j)
            endif
         enddo
         fraw(i)=(up*ux+vp*uy+wp*uz)*deim_eval_weights(i)
      enddo

      if (deimmode.eq.'MCLSDEIM') then
         do i=1,nbnl
            c_hat(i)=0.
            do j=1,ndeim_pts_eval
               c_hat(i)=c_hat(i)+deim_nl_bas_p_eval(j,i)*fraw(j)
            enddo
            do j=1,nbnl
               c_hat(i)=c_hat(i)+deim_alpha*deim_tau(i,j)*deim_mu(j)
            enddo
         enddo

         do i=1,nbnl
            tmpb(i)=0.
            do j=1,nbnl
               tmpb(i)=tmpb(i)+deim_A_tau_inv(i,j)*c_hat(j)
            enddo
         enddo

         do i=1,nbnl
            c_hat(i)=tmpb(i)
         enddo

         do i=1,nbnl
            bvec(i)=0.
            do j=1,nb
               bvec(i)=bvec(i)+deim_proj_mat(j,i)*uu(j)
            enddo
         enddo

         do i=1,nbnl
            tmpb(i)=0.
            do j=1,nbnl
               tmpb(i)=tmpb(i)+deim_A_tau_inv(i,j)*bvec(j)
            enddo
         enddo

         s1=0.
         s2=0.
         do i=1,nbnl
            s1=s1+bvec(i)*c_hat(i)
            s2=s2+bvec(i)*tmpb(i)
         enddo

         if (abs(s2).gt.0.) then
            lambda=s1/s2
            do i=1,nbnl
               c_hat(i)=c_hat(i)-lambda*tmpb(i)
            enddo
         endif
      else
         do i=1,nbnl
            do j=1,ndeim_pts_eval
               c_hat(i)=c_hat(i)+deim_interp_mat(i,j)*fraw(j)
            enddo
         enddo

         if (deimmode.eq.'CLSDEIM') then
            do i=1,nbnl
               do j=1,nb
                  bvec(i)=bvec(i)+deim_proj_mat(j,i)*uu(j)
               enddo
            enddo

            do i=1,nbnl
               do j=1,nbnl
                  tmpb(i)=tmpb(i)+deim_Ainv(i,j)*bvec(j)
               enddo
            enddo

            s1=0.
            s2=0.
            do i=1,nbnl
               s1=s1+bvec(i)*c_hat(i)
               s2=s2+bvec(i)*tmpb(i)
            enddo

            if (abs(s2).gt.0.) then
               lambda=s1/s2
               do i=1,nbnl
                  c_hat(i)=c_hat(i)-lambda*tmpb(i)
               enddo
            endif
         endif
      endif

      do i=1,nb
         do j=1,nbnl
            cu(i)=cu(i)+deim_proj_mat(i,j)*c_hat(j)
         enddo
         do j=0,nb
            cu(i)=cu(i)+deim_zmc(i,j)*uu(j)
         enddo
      enddo

      return
      end
