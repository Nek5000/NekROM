c-----------------------------------------------------------------------
      subroutine read_serial(a,n,fname,wk,nid)

      ! read in array

      ! a     := read target array
      ! n     := number of items
      ! fname := file name
      ! wk    := work array
      ! nid   := id of core

      character*128 fname
      character*128 fntrunc

      real a(n),wk(n)

      if (nid.le.0) then
         call blank(fntrunc,128)
         len=ltruncr(fname,128)
         call chcopy(fntrunc,fname,len)
         call read_serial_helper(a,n,fntrunc)
      else
         call rzero(a,n)
      endif

      if (nid.ge.0) call gop(a,wk,'+  ',n)

      return
      end
c-----------------------------------------------------------------------
      subroutine read_serial_helper(a,n,fname)

      ! core reading routine

      ! a     := read target array
      ! n     := number of items
      ! fname := file name

      real a(n)
      character*128 fname

      open (unit=12,file=fname)
      read (12,*) (a(i),i=1,n)
      close (unit=12)

      return
      end
c-----------------------------------------------------------------------
      subroutine read_mat_serial(a,n1,n2,fname,m1,m2,wk,nid)

      ! read in matrix

      ! a     := read target matrix
      ! n1    := number of rows of a
      ! n2    := number of columns of a
      ! fname := file name
      ! m1    := number of rows of fname
      ! m2    := number of columns of fname
      ! wk    := work array
      ! nid   := id of core

      character*128 fname
      character*128 fntrunc

      real a(n1,n2),wk(n1*n2)

      if (nid.eq.0) then
         call blank(fntrunc,128)
         len=ltruncr(fname,128)
         call chcopy(fntrunc,fname,len)
         call read_mat_serial_helper(a,n1,n2,fntrunc,m1,m2)
      else
         call rzero(a,n1*n2)
      endif

      call gop(a,wk,'+  ',n1*n2)

      return
      end
c-----------------------------------------------------------------------
      subroutine read_mat_serial_helper(a,n1,n2,fname,m1,m2)

      ! matrix core reader

      ! a     := read target matrix
      ! n1    := number of rows of a
      ! n2    := number of columns of a
      ! fname := file name
      ! m1    := number of rows of fname
      ! m2    := number of columns of fname

      real a(n1,n2)
      character*128 fname

      open (unit=12,file=fname)

      do j=1,m2
      do i=1,m1
         read (12,*) b
         if (j.le.n2.and.i.le.n1) a(i,j)=b
      enddo
      enddo

      close (unit=12)

      return
      end
c-----------------------------------------------------------------------
      subroutine loadbases

      ! set POD bases ub,vb,wb,pb,tb according to parameter flags

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      character*128 fname
      character*1   fn1(128)
      character*5   fnum

      logical ifexist

      equivalence (fname,fn1)

      if (nio.eq.0) write (6,*) 'inside loadbases'

      n=lx1*ly1*lz1*nelt
      n2=lx2*ly2*lz2*nelt

      call push_sol(vx,vy,vz,pr,t)

      inquire (file='bas.list',exist=ifexist)

      if (ifexist) then
         nn=nb+1
         do i=0,ldimt1
            ifreads(i)=ifrom(i)
         enddo

         call read_fields(
     $      us0,prs,ts0,nn,ls,0,ifreads,tk,'bas.list ',.false.)

         do i=0,nb
            if (ifrom(0)) call copy(pb(1,i),prs(1,i+1),n2)
            if (ifrom(1)) call opcopy(ub(1,i),vb(1,i),wb(1,i),
     $                        us0(1,1,i+1),us0(1,2,i+1),us0(1,ldim,i+1))
            if (ifrom(2)) call copy(tb(1,i,1),ts0(1,i+1,1),n)
         enddo
         if (nn.lt.nb) call exitti(
     $   'number of files in bas.list fewer than nb$',nb-nn)
      else
         do i=0,nb 
            len=ltrunc(session,132)
            call chcopy(fn1,'bas',3)
            call chcopy(fn1(4),session,len)
            call chcopy(fn1(4+len),'0.f',3)
            write (fnum,'(i5.5)') i+1
            call chcopy(fn1(7+len),fnum,5)

            inquire (file=fname,exist=ifexist)
            if (.not.ifexist)
     $        call exitti('missing basis file, exiting...$',i+1)

            call restart_filen(fname,11+len)
            if (ifrom(0)) call copy(pb(1,i),pr,n2)
            if (ifrom(1)) call opcopy(ub(1,i),vb(1,i),wb(1,i),vx,vy,vz)
            if (ifrom(2)) call copy(tb(1,i,1),t,n)
         enddo
      endif

    1 continue

      call pop_sol(vx,vy,vz,pr,t)

      if (nio.eq.0) write (6,*) 'exiting loadbases'

      return
      end
c-----------------------------------------------------------------------
      subroutine read_fields(usave,psave,tsave,ns,ls,nskp,ifread,
     $                       tk,fn,ifa)

      ! Reads and stores field files in given arrays

      ! usave,psave,tsave := velocity, pressure, temperature storage
      ! ns                := number of fields to save
      ! ls                := number of fields in (u,p,t)save
      ! nskp              := skipping interval
      ! ifread            := flags for reading each field
      ! tk                := time at each snapshot
      ! fn                := file name of the list file
      ! ifa               := to average or not to average

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (lt2=lx2*ly2*lz2*lelt)

      real usave(lt,ldim,1),psave(lt2,1),tsave(lt,ls,ldimt)
      real tk(1)

      character*128 fn
      character*128 fnlint

      logical ifa,ifread(0:ldimt1)

      common /scrk2/ t4(lt),t5(lt),t6(lt)

      call nekgsync
      rf_time=dnekclock()

      ierr = 0
      call lints(fnlint,fn,128)
      if (nid.eq.0) open(77,file=fnlint,status='old',err=199)
      ierr = iglmax(ierr,1)
      if (ierr.gt.0) goto 199

      n = lx1*ly1*lz1*nelt
      n2 = lx2*ly2*lz2*nelv

      call push_sol(vx,vy,vz,pr,t)

      if (ifa) then
         call zero_sol(uavg,vavg,wavg,pavg,tavg)
         call zero_sol(urms,vrms,wrms,prms,trms)
         call opzero(vwms,wums,uvms)
      endif

      call opcopy(t4,t5,t6,xm1,ym1,zm1)

      if (nio.eq.0) write (6,*) 'ns:',ns

      icount = 0
      do ip=1,ns
         do i=1,nskp+1
            call blank(initc,127)
            initc(1) = 'done '
            if (nid.eq.0) read(77,127,end=998) initc(1)
  998       call bcast(initc,127)
  127       format(a127)
         enddo
         if (nio.eq.0) write (6,*) ip,' '

         if (indx1(initc,'done ',5).eq.0) then ! We're not done
            icount = icount+1
            nfiles = 1
            ttmp=time
            call restart(nfiles)  ! Note -- time is reset.
            tk(icount)=time
            time=ttmp

            if (ifa) then
               call add_sol(uavg,vavg,wavg,pavg,tavg,vx,vy,vz,pr,t)
               call add2col2(urms,vx,vx,n)
               call add2col2(vrms,vy,vy,n)
               if (ldim.eq.3) call add2col2(wrms,vz,vz,n)
               call add2col2(trms,t,t,n)
               call add2col2(vwms,vx,t,n)
               call add2col2(wums,vy,t,n)
               if (ldim.eq.3) call add2col2(uvms,vz,t,n)
            endif

            if (ifread(0)) call copy(psave(1,ip),pr,n2)
            if (ifread(1))
     $         call opcopy(usave(1,1,ip),usave(1,2,ip),usave(1,ldim,ip),
     $                    vx,vy,vz)

            do j=1,ldimt
               idx=j+1
               if (ifread(idx)) call copy(tsave(1,ip,j),t(1,1,1,1,j),n)
            enddo
         else
            goto 999
         endif
      enddo

  999 continue  ! clean up averages
      if (nid.eq.0) close(77)

      ns = icount ! Actual number of files read

      call pop_sol(vx,vy,vz,pr,t)

      if (ifa) then
         s=1./real(ns)

         call scale_sol(uavg,vavg,wavg,pavg,tavg,s)
         call scale_sol(urms,vrms,wrms,prms,trms,s)
         call opcmult(vwms,wums,uvms,s,n)
      endif

      call opcopy(xm1,ym1,zm1,t4,t5,t6)

      call nekgsync
      if (nio.eq.0) write (6,*) 'rf_time:',dnekclock()-rf_time

      return

  199 continue ! exception handle for file not found

      call nekgsync
      if (nio.eq.0) write (6,*) 'rf_time:',dnekclock()-rf_time

      ierr = 1
      if (nid.eq.0) ierr = iglmax(ierr,1)
      write (6,*) fnlint
      call exitti('read_fields did not find list file.$',ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine restart_filen(fname127,nch)

      ! restart by specifying the number of file name characters

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      character*127 fname127

      parameter (lt=lx1*ly1*lz1*lelt)

      character*1 s1(127)
      equivalence (s1,initc) ! equivalence to initial condition

      logical iffexist

      call blank(initc,127)
      call chcopy(initc,fname127,nch)

      iblank = indx1(initc,' ',1)-1
      if (nio.eq.0) write(6,1) ipass,(s1(k),k=1,iblank)
    1 format(i8,'Reading: ',127a1)

      inquire (file=initc(1),exist=iffexist)

      if (iffexist) then
         nfiles = 1
         ttmp=time
         call restart(nfiles)
         time=ttmp
      else
         if (nio.eq.0) write (6,*) initc(1),'did not exist...'
         call exitt0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine restart_file(fname127)

      ! restart w/o specifying the number of file name characters

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      character*127 fname127

      character*1 s1(127)
      equivalence (s1,initc) ! equivalence to initial condition

      logical iffexist

      call blank(initc,127)
      nch=ltruncr(fname127,127)
      if (nio.eq.0) write (6,*) 'nch=',nch

      call chcopy(initc,fname127,nch)

      iblank = indx1(initc,' ',1)-1
      if (nio.eq.0) write(6,1) ipass,(s1(k),k=1,iblank)
    1 format(i8,'Reading: ',127a1)

      inquire (file=initc(1),exist=iffexist)

      if (iffexist) then
         nfiles = 1
         ttmp=time
         call restart(nfiles)
         time=ttmp
      else
         if (nio.eq.0) write (6,*) fname127,'did not exist...'
         call exitt0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine restart_time_file(fname127)

      ! restart and change 'time' to the value in the file

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      character*127 fname127

      character*1 s1(127)
      equivalence (s1,initc) ! equivalence to initial condition

      logical iffexist

      call blank(initc,127)
      nch=ltruncr(fname127,127)
      write (6,*) 'nch=',nch

      call chcopy(initc,fname127,nch)

      iblank = indx1(initc,' ',1)-1
      if (nio.eq.0) write(6,1) ipass,(s1(k),k=1,iblank)
    1 format(i8,'Reading: ',127a1)

      inquire (file=initc(1),exist=iffexist)

      if (iffexist) then
         nfiles = 1
         call restart(nfiles)
      else
         if (nio.eq.0) write (6,*) fname127,'did not exist...'
         call exitt0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine read_sigma_u_serial(a,n1,n2,fname,m1,m2,m3,m4,wk,nid)

      ! TODO: add description

      character*128 fname
      character*128 fntrunc

      real a(n1,n2),wk(n1*n2)

      if (nid.eq.0) then
         call blank(fntrunc,128)
         len=ltruncr(fname,128)
         call chcopy(fntrunc,fname,len)
         call read_sigma_u_serial_helper(a,n1,n2,fntrunc,m1,m2,m3,m4)
      else
         call rzero(a,n1*n2)
      endif

      call gop(a,wk,'+  ',n1*n2)

      return
      end
c-----------------------------------------------------------------------
      subroutine read_sigma_u_serial_helper(a,n1,n2,fname,m1,m2,m3,m4)

      ! TODO: add description

      real a(n1,n2)
      character*128 fname

      open (unit=12,file=fname)

      l1=1
      l2=1
      do j=1,m3
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      write(6,*)'l1,l2',l1,l2
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      write(6,*)'l1,l2',l1,l2
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      write(6,*)'l1,l2',l1,l2
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      write(6,*)'l1,l2',l1,l2
      do k=1,m3
         do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
         enddo
      enddo
      write(6,*)'l1,l2',l1,l2
      if (j.le.m4) l2=l2+1
      l1=1
      enddo

      do j=1,m3
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      do k=1,m3
         do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
         enddo
      enddo
      if (j.le.m4) l2=l2+1
      l1=1
      enddo

      do j=1,m3
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      do k=1,m3
         do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
         enddo
      enddo
      if (j.le.m4) l2=l2+1
      l1=1
      enddo

      do j=1,m3
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      do k=1,m3
         do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
         enddo
      enddo
      if (j.le.m4) l2=l2+1
      l1=1
      enddo

      do j=1,m3
         do s=1,m3
            do i=1,m3
               read (12,*) b
               if (s.le.m4.and.i.le.m4) then
                  a(l1,l2)=b
                  l1=l1+1
               endif
            enddo
            do i=1,m3
               read (12,*) b
               if (s.le.m4.and.i.le.m4) then
                  a(l1,l2)=b
                  l1=l1+1
               endif
            enddo
            do i=1,m3
               read (12,*) b
               if (s.le.m4.and.i.le.m4) then
                  a(l1,l2)=b
                  l1=l1+1
               endif
            enddo
            do i=1,m3
               read (12,*) b
               if (s.le.m4.and.i.le.m4) then
                  a(l1,l2)=b
                  l1=l1+1
               endif
            enddo
            do k=1,m3
               do i=1,m3
               read (12,*) b
               if (s.le.m4.and.i.le.m4) then
                  a(l1,l2)=b
                  l1=l1+1
               endif
               enddo
            enddo
            if (s.le.m4) l2=l2+1
            l1=1
         enddo
      enddo

      close (unit=12)

      return
      end
c-----------------------------------------------------------------------
      subroutine read_sigma_t_serial(a,n1,n2,fname,m1,m2,m3,m4,wk,nid)

      ! TODO: add description

      character*128 fname
      character*128 fntrunc

      real a(n1,n2),wk(n1*n2)

      if (nid.eq.0) then
         call blank(fntrunc,128)
         len=ltruncr(fname,128)
         call chcopy(fntrunc,fname,len)
         call read_sigma_t_serial_helper(a,n1,n2,fntrunc,m1,m2,m3,m4)
      else
         call rzero(a,n1*n2)
      endif

      call gop(a,wk,'+  ',n1*n2)

      return
      end
c-----------------------------------------------------------------------
      subroutine read_sigma_t_serial_helper(a,n1,n2,fname,m1,m2,m3,m4)

      ! TODO: add description

      real a(n1,n2)
      character*128 fname

      open (unit=12,file=fname)

      l1=1
      l2=1
      do j=1,m3
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      write(6,*)'l1,l2',l1,l2
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      write(6,*)'l1,l2',l1,l2
      do k=1,m3
         do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
         enddo
      enddo
      write(6,*)'l1,l2',l1,l2
      if (j.le.m4) l2=l2+1
      l1=1
      enddo

      do j=1,m3
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
      enddo
      do k=1,m3
         do i=1,m3
         read (12,*) b
         if (j.le.m4.and.i.le.m4) then
            a(l1,l2)=b
            l1=l1+1
         endif
         enddo
      enddo
      if (j.le.m4) l2=l2+1
      l1=1
      enddo

      do j=1,m3
         do s=1,m3
            do i=1,m3
               read (12,*) b
               if (s.le.m4.and.i.le.m4) then
                  a(l1,l2)=b
                  l1=l1+1
               endif
            enddo
            do i=1,m3
               read (12,*) b
               if (s.le.m4.and.i.le.m4) then
                  a(l1,l2)=b
                  l1=l1+1
               endif
            enddo
            do k=1,m3
               do i=1,m3
               read (12,*) b
               if (s.le.m4.and.i.le.m4) then
                  a(l1,l2)=b
                  l1=l1+1
               endif
               enddo
            enddo
            if (s.le.m4) l2=l2+1
            l1=1
         enddo
      enddo

      close (unit=12)

      return
      end
c-----------------------------------------------------------------------
      subroutine loadpbases(nsave)

      ! This subroutine is for p-greedy. It reads
      ! in files in pbas.list and return nsave for how
      ! many files it has read in.

      ! This subroutine should be called before calling
      ! subroutine projtoprerb

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      character*128 fname
      character*1   fn1(128)
      character*5   fnum

      logical ifexist

      equivalence (fname,fn1)

      if (nio.eq.0) write (6,*) 'inside loadpbases'

      n=lx1*ly1*lz1*nelt
      n2=lx2*ly2*lz2*nelt

      call push_sol(vx,vy,vz,pr,t)

      inquire (file='pbas.list',exist=ifexist)

      if (ifexist) then
         nn=nb
         nsu=1
         nsp=1
         nst=1
         if (ifrom(0)) nsp=nn
         if (ifrom(1)) nsu=nn
         if (ifrom(2)) nst=nn

         call get_p_rb(nsave,nsu,nsp,nst,timek,'pbas.list ')
         if (nn.lt.nb) call exitti(
     $   'number of files in pbas.list fewer than nb$',nb-nn)
      else
         do i=0,nb
            len=ltrunc(session,132)
            call chcopy(fn1,'pbas',4)
            call chcopy(fn1(5),session,len)
            call chcopy(fn1(5+len),'0.f',3)
            write (fnum,'(i5.5)') i+1
            call chcopy(fn1(8+len),fnum,5)

            inquire (file=fname,exist=ifexist)
            if (.not.ifexist)
     $        call exitti('missing pbas file, exiting...$',i+1)

            call restart_filen(fname,11+len)
            if (ifrom(0)) call copy(pb(1,i),pr,n2)
            if (ifrom(1)) call opcopy(ub(1,i),vb(1,i),wb(1,i),vx,vy,vz)
            if (ifrom(2)) call copy(tb(1,i,1),t,n)
         enddo
      endif

    1 continue

      call pop_sol(vx,vy,vz,pr,t)

      if (nio.eq.0) write (6,*) 'exiting loadpbases'

      return
      end
c-----------------------------------------------------------------------
      subroutine get_p_rb(nocp,nsu,nsp,nst,ttk,fn)

c     This routine reads files specificed in fname

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'
      include 'AVG'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (lt2=lx2*ly2*lz2*lelt)

      real ttk(nsu)
      integer nocp

      character*128 fn
      character*128 fnlint

      common /scrk2/ t4(lt),t5(lt),t6(lt)

      call nekgsync
      gsf_time=dnekclock()

      ierr = 0
      call lints(fnlint,fn,128)
      if (nid.eq.0) open(77,file=fnlint,status='old',err=199)
      ierr = iglmax(ierr,1)
      if (ierr.gt.0) goto 199

      n = lx1*ly1*lz1*nelt
      n2= lx2*ly2*lz2*nelt

      call push_sol(vx,vy,vz,pr,t)
      call opcopy(t4,t5,t6,xm1,ym1,zm1)

      nsave=max(max(nsu,nsp),nst)

      if (nio.eq.0) write (6,*) 'nsu,nsp,nst:',nsu,nsp,nst

      icount = 0
      do ipass=1,nsave
         call blank(initc,127)
         initc(1) = 'done '
         if (nid.eq.0) read(77,127,end=998) initc(1)
  998    call bcast(initc,127)
  127    format(a127)
         if (nio.eq.0) write (6,*) ipass,' '

         if (indx1(initc,'done ',5).eq.0) then ! We're not done
            icount = icount+1
            nfiles = 1
            ttmp=time
            call restart(nfiles)  ! Note -- time is reset.
            ttk(icount)=time
            time=ttmp

            ip=ipass
            if (icount.le.nsu)
     $         call opcopy(ub(1,ip),vb(1,ip),wb(1,ip),
     $                    vx,vy,vz)

            if (icount.le.nsp) call copy(pb(1,ip),pr,n2)
            if (icount.le.nst) call copy(tb(1,ip,1),t,n)
         else
            goto 999
         endif
      enddo

  999 continue  ! clean up averages
      if (nid.eq.0) close(77)

      nsave = icount ! Actual number of files read
      nocp = nsave

      call pop_sol(vx,vy,vz,pr,t)

      call nekgsync
      if (nio.eq.0) write (6,*) 'gsf_time:',dnekclock()-gsf_time

      return

  199 continue ! exception handle for file not found

      call nekgsync
      if (nio.eq.0) write (6,*) 'gsf_time:',dnekclock()-gsf_time

      ierr = 1
      if (nid.eq.0) ierr = iglmax(ierr,1)
      write (6,*) fnlint
      call exitti('get_p_rb did not find list file.$',ierr)

      return
      end
c-----------------------------------------------------------------------
