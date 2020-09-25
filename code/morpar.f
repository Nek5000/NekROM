c-----------------------------------------------------------------------
      subroutine readat_mor
C
C     Read in run parameters from .mor file
C
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'

      if (nid.eq.0) call mor_read(ierr)
      call bcast(ierr,isize)
      if (ierr .ne. 0) call exitt
      call bcastparam

      return
      end
c-----------------------------------------------------------------------
      subroutine mor_read(ierr)
c
c     parse .mor file and set run parameters
c
      include 'SIZE'
      include 'INPUT'
      include 'ADJOINT'
      include 'RESTART'
      include 'PARALLEL'
      include 'CTIMER'
      include 'TSTEP'
      include 'MOR'

      call set_morfle

      character*132 c_out,txt,txt2
      character*3 chartmp

      call finiparser_load(morfle,ierr)
      if (ierr.ne.0) return

      call mor_verify(ierr)
      if (ierr.ne.0) return

      ! general

      call finiparser_getstring(c_out,'general:mode',ifnd)
      if (ifnd.eq.1) then
         call capit(c_out,132)
         if (index(c_out,'ALL').eq.1) then
            rmode='ALL'
         else if (index(c_out,'ONB').eq.1) then
            rmode='ONB'
         else if (index(c_out,'ON').eq.1) then
            rmode='ON'
         else if (index(c_out,'OFF').eq.1) then
            rmode='OFF'
         else if (index(c_out,'CP').eq.1) then
            rmode='CP '
            max_tr=ltr
         else
            write (6,*) 'invalid option for general:mode ',c_out
            ierr=ierr+1
         endif
      endif

      if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
         open (unit=10,file='ops/ips')
         read (10,*) chartmp
         close (unit=10)
         if (chartmp.ne.ips) then
            write (6,*) 'online ips does not match offline ips',nb
            ierr=ierr+1
         endif
      endif

      ifrecon=(rmode.ne.'ON '.and.rmode.ne.'CP ')

      call finiparser_getstring(c_out,'general:field',ifnd)
      if (ifnd.eq.1) then
         call capit(c_out,132)
         if (index(c_out,'v').eq.1) then
            ifpod(1)=.true.
         else if (index(c_out,'t').eq.1) then
            rmode='ONB'
         else if (index(c_out,'vt').eq.1) then
            ifpod(1)=.true.
            ifpod(2)=ifheat
         else
            write (6,*) 'invalid option for general:field ',c_out
            ierr=ierr+1
         endif
      endif

      ifrom(1)=ifpod(1)
      ifpod(1)=ifpod(1).or.ifrom(2)

      call finiparser_getdbl(d_out,'general:nb',ifnd)
      if (ifnd.eq.1) nb=min(nint(d_out),lb)
      if (nb.eq.0) nb=lb

      if (ifavg0.and.(nb.eq.ls)) then
         write (6,*) 'nb == ls results in linear dependent bases',nb
         ierr=ierr+1
      endif

      if (nb.gt.ls) then
         write (6,*) 'nb > ls is undefined configuration',nb
         ierr=ierr+1
      endif

      call finiparser_getbool(i_out,'general:ei',ifnd)
      if (ifnd.eq.1) ifei=i_out.eq.1

      call finiparser_getdbl(d_out,'general:avginit',ifnd)
      if (ifnd.eq.1) navg_step=nint(max(1.,d_out))

      call finiparser_getdbl(d_out,'general:rktol',ifnd)
      if (ifnd.eq.1) rktol=d_out
      if (rktol.lt.0.) rktol=10.**rktol

      call finiparser_getbool(i_out,'general:tneubc',ifnd)
      if (ifnd.eq.1) iftneu=i_out.eq.1

      call finiparser_getdbl(d_out,'general:nplay',ifnd)
      if (ifnd.eq.1) nplay=nint(d_out)

      if (nplay.ne.0) then
         nplay=max(nplay,0)
         ifplay=.true.
      else
         ifplay=.false.
      endif

      ! POD

      call finiparser_getstring(c_out,'pod:type',ifnd)
      if (ifnd.eq.1) then
         call capit(c_out,132)
         if (index(c_out,'L2').eq.1) then
            ips='L2 '
         else if (index(c_out,'H10').eq.1) then
            ips='H10'
         else if (index(c_out,'HLM').eq.1) then
            ips='HLM'
         else
            write (6,*) 'invalid option for pod:type ',c_out
            ierr=ierr+1
         endif
      endif

      call finiparser_getstring(c_out,'pod:mode0',ifnd)
      if (ifnd.eq.1) then
         call capit(c_out,132)
         if (index(c_out,'AVG').eq.1) then
            ifavg0=.true.
         else if (index(c_out,'STATE').eq.1) then
            ifavg0=.false.
         else
            write (6,*) 'invalid option for pod:mode0 ,c_out
            ierr=ierr+1
         endif
      endif

      ! QOI

      call finiparser_getdbl(d_out,'qoi:freq',ifnd)
      if (ifnd.eq.1) ad_qstep=nint(d_out)+ad_iostep*max(1-nint(d_out),0)

      call finiparser_getbool(i_out,'qoi:drag',ifnd)
      if (ifnd.eq.1) ifcdrag=i_out.eq.1

      call finiparser_getbool(i_out,'qoi:tke',ifnd)
      if (ifnd.eq.1) ifctke=i_out.eq.1

      call finiparser_getdbl(d_out,'qoi:nus',ifnd)
      if (ifnd.eq.1) inus=min(max(nint(d_out),0),5)

      ! COPT

      call finiparser_getstring(c_out,'copt:mode',ifnd)
      if (ifnd.eq.1) then
         call capit(c_out,132)
         if (index(c_out,'ON').eq.1) then
            isolve=1
         else if (index(c_out,'OFF').eq.1) then
            isolve=0
         else if (index(c_out,'SELECT').eq.1) then
            isolve=2
         else
            write (6,*) 'invalid option for copt:mode ',c_out
            ierr=ierr+1
         endif
      endif

      call finiparser_getstring(c_out,'copt:field',ifnd)
      if (ifnd.eq.1) then
         call capit(c_out,132)
         if (index(c_out,'v').eq.1) then
            icopt=2
         else if (index(c_out,'t').eq.1) then
            icopt=1
         else if (index(c_out,'vt').eq.1) then
            icopt=0
         else
            write (6,*) 'invalid option for copt:field ',c_out
            ierr=ierr+1
         endif
      endif

      call finiparser_getstring(c_out,'copt:barrier',ifnd)
      if (ifnd.eq.1) then
         call capit(c_out,132)
         if (index(c_out,'LOG').eq.1) then
            barr_func=1.
         else if (index(c_out,'INV').eq.1) then
            barr_func=0.
         else
            write (6,*) 'invalid option for copt:barrier ',c_out
            ierr=ierr+1
         endif
      endif

      call finiparser_getdbl(d_out,'copt:boxtol',ifnd)
      if (ifnd.eq.1) box_tol=d_out

      call finiparser_getdbl(d_out,'copt:vpar0',ifnd)
      if (ifnd.eq.1) ubarr0=d_out

      call finiparser_getdbl(d_out,'copt:vnloop',ifnd)
      if (ifnd.eq.1) ubarrseq=nint(d_out)

      call finiparser_getdbl(d_out,'copt:tpar0',ifnd)
      if (ifnd.eq.1) tbarr0=d_out

      call finiparser_getdbl(d_out,'copt:tnloop',ifnd)
      if (ifnd.eq.1) tbarrseq=nint(d_out)

      ! Fast

      call finiparser_getbool(i_out,'fast:ceval',ifnd)
      if (ifnd.eq.1) iffastc=i_out.eq.1

      call finiparser_getbool(i_out,'fast:heval',ifnd)
      if (ifnd.eq.1) iffasth=i_out.eq.1.and.ips.eq.'HLM'

      call finiparser_getbool(i_out,'forcing:body',ifnd)
      if (ifnd.eq.1) ifforce=i_out.eq.1

      call finiparser_getbool(i_out,'forcing:source',ifnd)
      if (ifnd.eq.1) ifsource=i_out.eq.1

      call finiparser_getbool(i_out,'forcing:buoyancy',ifnd)
      if (ifnd.eq.1) ifbuoy=i_out.eq.1

      ! Filter

      call finiparser_getstring(c_out,'filter:location',ifnd)
      if (ifnd.eq.1) then
         call capit(c_out,132)
         if (index(c_out,'NONE').eq.1) then
            rfilter='STD'
         else if (index(c_out,'CONV').eq.1) then
            rfilter='LER'
         else if (index(c_out,'POST').eq.1) then
            rfilter='EF '
         else
            write (6,*) 'invalid option for filter:location ',c_out
            ierr=ierr+1
         endif
      endif

      call finiparser_getstring(c_out,'filter:type',ifnd)
      if (ifnd.eq.1) then
         call capit(c_out,132)
         if (index(c_out,'TFUNC').eq.1) then
            call finiparser_getdbl(d_out,'filter:modes',ifnd)
            if (ifnd.eq.1) then
               rbf=d_out
               if (d_out.ge.0) then
                  rbf=min(nint(d_out),nb-1)
               else
                  rbf=nint(nb*min(1.,-d_out))
               endif
            else
               write (6,*) 'transfer filter needs a filter:modes value'
               ierr=ierr+1
            endif
         else if (index(c_out,'DIFF').eq.1) then
            call finiparser_getdbl(d_out,'filter:radius',ifnd)
            if (ifnd.eq.1) then
               rdft=d_out
            else
               write (6,*) 'diff. filter needs a filter:radius value'
               ierr=ierr+1
            endif
         else
            write (6,*) 'invalid option for filter:type ',c_out
            ierr=ierr+1
         endif
      endif

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF') then
         rtmp1(1,1)=nb*1.
         call dump_serial(rtmp1(1,1),1,'ops/nb ',nid)
      else
         call read_serial(rtmp1(1,1),1,'ops/nb ',b,nid)
         mb=rtmp1(1,1)
         if (mb.lt.nb) then
            write (6,*) 'mb less than nb... ',mb
         ierr=ierr+1
      endif

100   if (ierr.eq.0) call finiparser_dump()

      return
      end
c-----------------------------------------------------------------------
      subroutine bcastrpar
C
C     Broadcast mor parameters to all processors
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      INCLUDE 'RESTART'
      INCLUDE 'PARALLEL'
      INCLUDE 'CTIMER'
      INCLUDE 'ADJOINT'
      INCLUDE 'CVODE'

c     call bcast(loglevel, isize)
c     call bcast(optlevel, isize)

c     call bcast(param , 200*wdsize)
c     call bcast(uparam, 20*wdsize)

c     call bcast(filterType, wdsize)

c     call bcast(atol ,  (ldimt1+1)*wdsize)
c     call bcast(restol, (ldimt1+1)*wdsize)

c     call bcast(ifchar , lsize)
c     call bcast(iftran  , lsize)
c     call bcast(ifflow  , lsize)
c     call bcast(ifheat  , lsize)     
c     call bcast(iflomach, lsize)
c     call bcast(ifstrs  , lsize)
c     call bcast(ifmvbd  , lsize)
c     call bcast(ifusermv, lsize)
c     call bcast(ifdp0dt,  lsize)
c     call bcast(ifaxis  , lsize)
c     call bcast(ifcyclic, lsize)
c     call bcast(ifmhd   , lsize)
c     call bcast(ifuservp, lsize)
c     call bcast(ifpert, lsize)
c     call bcast(ifbase, lsize)
c     call bcast(ifmoab, lsize)
c     call bcast(ifaziv, lsize)
c     call bcast(ifadj , lsize)

c     call bcast(ifadvc ,  ldimt1*lsize)
c     call bcast(ifdiff ,  ldimt1*lsize)
c     call bcast(ifdeal ,  ldimt1*lsize)
c     call bcast(iffilter, ldimt1*lsize)

c     call bcast(idpss    ,  ldimt*isize)
c     call bcast(iftmsh   , (ldimt1+1)*lsize)
c     call bcast(ifprojfld, (ldimt1+1)*lsize)

c     call bcast(cpfld, 3*ldimt1*wdsize)

c     call bcast(ifxyo , lsize)
c     call bcast(ifvo  , lsize)
c     call bcast(ifpo  , lsize)
c     call bcast(ifto  , lsize)
c     call bcast(ifpsco, ldimt1*lsize)

c     call bcast(initc, 15*132*csize) 

c     call bcast(timeioe,sizeof(timeioe))

      return
      END
c-----------------------------------------------------------------------
      subroutine mor_verify(ierr)

      INCLUDE 'MORDICT'
      
      character*132  key
      character*1024 val

      character*132 txt
      character*1   tx1(132)
      equivalence   (tx1,txt)

      ierr = 0

      call finiparser_getDictEntries(n)
      do i = 1,n
         call finiparser_getPair(key,val,i,ifnd)
         call capit(key,132)

         is = index(key,'_') ! ignore user keys
         if (is.eq.1) goto 10

         do j = 1,mordict_nkeys ! do we find the key in the mar-dictionary
            if (index(mordictkey(j),key).eq.1) goto 10
         enddo
         write (6,*) 'ERROR: .mor file contains unknown key ', key
         ierr = ierr + 1
   10 enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_morfle

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'MOR'

      character*132 name
      character*1 sess1(132),path1(132),nam1(132)
      equivalence (session,sess1)
      equivalence (path,path1)
      equivalence (name,nam1)
      character*1 mor(4)
      character*4 mor4
      equivalence (mor,mor4)
      character*78  string

      character*132 morfle data  /'.mor'/

      len = ltrunc(path,132)
      if (indx1(path1(len),'/',1).lt.1) then
         call chcopy(path1(len+1),'/',1)
      endif

      call blank(morfle,132)
      call blank(name,132)

      ls=ltrunc(session,132)
      lpp=ltrunc(path,132)
      lsp=ls+lpp

      call chcopy(nam1(    1),path1,lpp)
      call chcopy(nam1(lpp+1),sess1,ls )
      l1 = lpp+ls+1
      ln = lpp+ls+4

      call chcopy(nam1  (l1),mor , 4)
      call chcopy(morfle    ,nam1,ln)

      return
      end
c-----------------------------------------------------------------------
