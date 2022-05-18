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
         do j=1,ldimt
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
      subroutine dump_misc

      ! dump miscellaneous items

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
