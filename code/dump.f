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

      parameter (lt=lx1*ly1*lz1*lelt)

      common /dumpglobal/ wk1(lcloc),wk2(lcloc)

      if (ifpod(1)) then
         call dump_serial(ug(1,1,1),ls*ls,'ops/gu ',nid)
         call dump_serial(au0,(nb+1)**2,'ops/au ',nid)
         call dump_serial(bu0,(nb+1)**2,'ops/bu ',nid)
         call dump_serial(u,(nb+1)*3,'ops/u ',nid)
         call dump_serial(uk,ns*(nb+1),'ops/uk ',nid)
         call dump_serial(umin,nb,'ops/umin ',nid)
         call dump_serial(umax,nb,'ops/umax ',nid)
         call dump_serial(timek,ns,'ops/timek ',nid)
         call dump_global(cul,ncloc,'ops/cu ',wk1,wk2,nid)
      endif

      if (ifpod(2)) then
         call dump_serial(ug(1,1,2),ls*ls,'ops/gt ',nid)
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

      if (ifforce) call dump_serial(rf,nb,'ops/rf ',nid)
      if (ifsource) call dump_serial(rq,nb,'ops/rq ',nid)

      if (ifbuoy) then
         call dump_serial(but0,(nb+1)**2,'ops/but ',nid)
      endif

      if (ifei) then
         l=1
         do j=1,nres
         do i=1,nres
            sigtmp(l,1)=sigma(i,j)
            l=l+1
         enddo
         enddo
         call dump_serial(sigtmp,nres*nres,'ops/sigma ',nid)
      endif

      ttmp=time
      itmp=istep
      do i=0,nb
         time=i
         itmp=i
         call outpost(ub(1,i),vb(1,i),wb(1,i),pb(1,i),tb(1,i),'bas')
      enddo
      istep=itmp
      time=ttmp

      return
      end
c-----------------------------------------------------------------------
