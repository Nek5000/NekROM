c-----------------------------------------------------------------------
      subroutine read_serial(a,n,fname,wk,nid)

      character*128 fname
      character*128 fntrunc

      real a(n)

      if (nid.eq.0) then
         call blank(fntrunc,128)
         len=ltruncr(fname,128)
         call chcopy(fntrunc,fname,len)
         call read_serial_helper(a,n,fntrunc)
      else
         call rzero(a,n)
      endif

      call gop(a,wk,'+  ',n)

      return
      end
c-----------------------------------------------------------------------
      subroutine read_serial_helper(a,n,fname)

      real a(n)
      character*128 fname

      open (unit=12,file=fname)
      read (12,*) (a(i),i=1,n)
      close (unit=12)

      return
      end
c-----------------------------------------------------------------------
      subroutine loadbases

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk1/ tmp1(lt),tmp2(lt),tmp3(lt),
     $               tmp4(lt),tmp5(lt),tmp6(lt,3)

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
         call get_saved_fields(us,ps,ts,nn,'bas.list ')
         do i=0,nb
            call copy_sol(ub(1,i),vb(1,i),wb(1,i),pb(1,i),tb(1,i,1),
     $       us(1,1,i+1),us(1,2,i+1),us(1,ldim,i+1),ps(1,i+1),ts(1,i,1))
         enddo
         if (nn.lt.nb) call exitti(
     $   'number of files in bas.list fewer than nb',nb-nn)
      else
         do i=0,nb 
            len=ltrunc(session,132)
            call chcopy(fn1,'bas',3)
            call chcopy(fn1(4),session,len)
            call chcopy(fn1(4+len),'0.f',3)
            write (fnum,'(i5.5)') i+1
            call chcopy(fn1(7+len),fnum,5)

            call restart_filen(fname,11+len)
            call opcopy(ub(1,i),vb(1,i),wb(1,i),vx,vy,vz)
         enddo
      endif

      call pop_sol(vx,vy,vz,pr,t)

      if (nio.eq.0) write (6,*) 'exiting loadbases'

      return
      end
c-----------------------------------------------------------------------
      subroutine get_saved_fields(usave,psave,tsave,nsave,fname)

c     This routine reads files specificed in fname

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (lt2=lx2*ly2*lz2*lelt)

      real usave(lt,ldim,nsave),psave(lt2,nsave),tsave(lt,ldimt,nsave)
      character*128 fname
      character*128 fnlint

      common /scrk5/ uu(lt),vv(lt),ww(lt),t1(lt),t2(lt,ldimt),t3(lt)

      ierr = 0
      call lints(fnlint,fname,128)
      if (nid.eq.0) open(77,file=fnlint,status='old',err=199)
      ierr = iglmax(ierr,1)
      if (ierr.gt.0) goto 199

      n = lx1*ly1*lz1*nelt
      n2= lx2*ly2*lz2*nelt

      call push_sol(vx,vy,vz,pr,t)
      call zero_sol(uavg,vavg,wavg,pavg,tavg)

      icount = 0
      do ipass=1,nsave
         call blank(initc,127)
         initc(1) = 'done '
         if (nid.eq.0) read(77,127,end=998) initc(1)
  998    call bcast(initc,127)
  127    format(a127)
         if (nio.eq.0) write (6,*) ipass,' '

         if (indx1(initc,'done ',5).eq.0) then ! We're not done
            nfiles = 1
            ttmp=time
            call restart(nfiles)  ! Note -- time is reset.
            time=ttmp

            ip=ipass
            call add_sol(uavg,vavg,wavg,pavg,tavg,vx,vy,vz,pr,t)
            call copy_sol(usave(1,1,ip),usave(1,2,ip),usave(1,ldim,ip),
     $                    psave(1,ip),tsave(1,1,ip),vx,vy,vz,pr,t)
            icount = icount+1
         else
            goto 999
         endif
      enddo

  999 continue  ! clean up averages
      if (nid.eq.0) close(77)

      nsave = icount ! Actual number of files read

      call pop_sol(vx,vy,vz,pr,t)

      s=1./real(nsave)
      call scale_sol(uavg,vavg,wavg,pavg,tavg,s)

      return

  199 continue ! exception handle for file not found
      ierr = 1
      if (nid.eq.0) ierr = iglmax(ierr,1)
      write (6,*) fnlint
      call exitti('get_saved_fields did not find list file.$',ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine restart_filen(fname127,nch)

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
