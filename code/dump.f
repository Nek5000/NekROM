c-----------------------------------------------------------------------
      subroutine dump_global(a,n,fname,wk1,wk2,nid)

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
      subroutine dump_serial(a,n,fname,nid)

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

      real a(n),wk1(1),wk2(1)
      integer iwk(1)

      character*128 fname

      if (nid.eq.0) open (unit=12,file=fname)

      iwk(1)=n
      nmax=iglmax(iwk,1)

      iwk(1)=nid
      ipmax=iglmax(iwk,1)

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
      subroutine dump_serial_helper(a,n,fname)

      real a(n)

      character*128 fname

      open (unit=12,file=fname)

      do i=1,n
         write (12,*) a(i)
      enddo

      close (unit=12)

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_all

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      common /dumpglobal/ wk1(lcloc),wk2(lcloc)

      call dump_serial(ug(1,1,1),ls*ls,'ops/gv ',nid)
      call dump_serial(av0,(nb+1)**2,'ops/av ',nid)
      call dump_serial(bv0,(nb+1)**2,'ops/bv ',nid)
      call dump_serial(u,(nb+1)*3,'ops/v ',nid)
      call dump_global(cvl,ncloc,'ops/cv ',wk1,wk2,nid)


      ttmp=time
      itmp=istep
      do i=0,nb
         time=i
         itmp=i
         call outpost(ub(1,i),vb(1,i),wb(1,i),pr,t,'bas')
      enddo
      istep=itmp
      time=ttmp

      return
      end
c-----------------------------------------------------------------------
