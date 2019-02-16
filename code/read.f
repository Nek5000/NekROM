c-----------------------------------------------------------------------
      subroutine read_serial(a,n,fname,wk,nid)

      character*128 fname
      character*128 fntrunc

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
      subroutine loadbases(ub,vb,wb,nb)

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk1/ tmp1(lt),tmp2(lt),tmp3(lt),tmp4(lt),tmp5(lt)
      real ub(lt,0:nb), vb(lt,0:nb), wb(lt,0:nb)

      character*128 fname
      character*1   fn1(128)
      character*5   fnum

      equivalence (fname,fn1)

      if (nio.eq.0) write (6,*) 'inside loadbases'

      n=lx1*ly1*lz1*nelt
      n2=lx2*ly2*lz2*nelt

      call opcopy(tmp1,tmp2,tmp3,vx,vy,vz)
      call copy(tmp4,pr,n2)
      call copy(tmp5,t,n)

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

      call opcopy(vx,vy,vz,tmp1,tmp2,tmp3)
      call copy(pr,tmp4,n2)
      call copy(t,tmp5,n)

      if (nio.eq.0) write (6,*) 'exiting loadbases'

      return
      end
c-----------------------------------------------------------------------
      subroutine get_saved_fields(usave,vsave,wsave,nsave,u0,ifvort)

c     This routine reads files specificed in file.list

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      parameter (lt=lx1*ly1*lz1*lelt)
      real usave(lt,nsave),vsave(lt,nsave),wsave(lt,nsave)
      real u0(lt,3) ! Initial condtion
      logical ifvort

      common /scrk5/ uu(lt),vv(lt),ww(lt),t1(lt),t2(lt),t3(lt)

      ierr = 0
      if (nid.eq.0) open(77,file='file.list',status='old',err=199)
      ierr = iglmax(ierr,1)
      if (ierr.gt.0) goto 199
      n = lx1*ly1*lz1*nelt
      n2= lx2*ly2*lz2*nelt

      call opcopy(uu,vv,ww,vx,vy,vz)

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
            call restart(nfiles)  ! Note -- time is reset.

!           Usave = U_snapshot - U_0:

            if (ifvort) then
               call comp_vort3(t1,t2,t3,vx,vy,vz)
               call sub3(usave(1,ipass),t1,u0,n)
            else
               call opsub3(usave(1,ipass),vsave(1,ipass),wsave(1,ipass),
     $                     vx,vy,vz,u0(1,1),u0(1,2),u0(1,3))
            endif
            icount = icount+1
         else
            goto 999
         endif

      enddo

      call opcopy(vx,vy,vz,uu,vv,ww)

  999 continue  ! clean up averages
      if (nid.eq.0) close(77)

      nsave = icount ! Actual number of files read

      return

  199 continue ! exception handle for file not found
      ierr = 1
      if (nid.eq.0) ierr = iglmax(ierr,1)
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
         call restart(nfiles)
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
         call restart(nfiles)
      else
         if (nio.eq.0) write (6,*) fname127,'did not exist...'
         call exitt0
      endif

      return
      end
c-----------------------------------------------------------------------
