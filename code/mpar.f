c-----------------------------------------------------------------------
      subroutine mor_set_params_par

      ! Read in run parameters from .mor file

      include 'SIZE'
      include 'PARALLEL'

      if (nid.eq.0) call mpar_read(ierr)
      call bcast(ierr,isize)
      if (ierr.ne.0) call exitt

      call bcastmpar

      return
      end
c-----------------------------------------------------------------------
      subroutine mpar_read(ierr)

      ! parse .mor file and set run parameters

      include 'SIZE'
      include 'INPUT'
      include 'ADJOINT'
      include 'RESTART'
      include 'PARALLEL'
      include 'CTIMER'
      include 'TSTEP'
      include 'MOR'

      character*132 c_out,txt,txt2,morfle
      character*3 chartmp

      call set_morfle(morfle)
      write (6,*) 'morfle: ',morfle

      call finiparser_load(morfle,ierr)
      if (ierr.ne.0) return

      call mpar_verify(ierr)
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
         else if (index(c_out,'AEQ').eq.1) then
            rmode='AEQ'
            ifpb=.false.
         else
            write (6,*) 'invalid option for general:mode ',c_out
            ierr=ierr+1
         endif
      endif

      ifrecon=rmode.ne.'ON '.and.rmode.ne.'CP '

      call finiparser_getstring(c_out,'general:field',ifnd)
      if (ifnd.eq.1) then
         call capit(c_out,132)

         if (index(c_out,'VT').eq.1) then
            ifpod(1)=.true.
            ifpod(2)=ifheat
         elseif (index(c_out,'V').eq.1) then
            ifpod(1)=.true.
         else if (index(c_out,'T').eq.1) then
            ifpod(2)=.true.
         else
            write (6,*) 'invalid option for general:field ',c_out
            ierr=ierr+1
         endif
      endif

      ifrom(1)=ifpod(1)
      ifpod(1)=ifpod(1).or.ifrom(2)
      ifrom(2)=ifpod(2)

      call finiparser_getdbl(d_out,'general:nb',ifnd)
      if (ifnd.eq.1) nb=min(nint(d_out),lb)
      if (nb.eq.0) nb=lb
      nbo=nb

      if (ifavg0.and.(nb.eq.ls)) then
         write (6,*) 'nb == ls results in linear dependent bases',nb
         ierr=ierr+1
      endif

      if (nb.gt.ls) then
         write (6,*) 'nb > ls is undefined configuration',nb
         ierr=ierr+1
      endif

      call finiparser_getdbl(d_out,'general:ns',ifnd)
      if (ifnd.eq.1) ns=min(nint(d_out),ls)
      if (ns.eq.0) ns=ls

      call finiparser_getdbl(d_out,'general:nskip',ifnd)
      nskip=0
      if (ifnd.eq.1) nskip=nint(d_out)

      call finiparser_getdbl(d_out,'general:avginit',ifnd)
      if (ifnd.eq.1) navg_step=nint(max(1.,d_out))

      do i=1,3
         its(i)=min(navg_step+i-1,3)
      enddo

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

      call finiparser_getbool(i_out,'general:decoupled',ifnd)
      if (ifnd.eq.1) ifdecpl=i_out.eq.1

      call finiparser_getbool(i_out,'general:cflow',ifnd)
      if (ifnd.eq.1) ifcflow=i_out.eq.1

      call finiparser_getbool(i_out,'general:eddy_vis',ifnd)
      if (ifnd.eq.1) then
         ifedvs=i_out.eq.1
         if (ldimt.lt.4) then
            write(6,*) 'recompile with ldimt=4!'
            ierr=ierr+1
         else
            ifrom(5)=ifedvs
            ifpod(5)=ifrom(5)
            ifpod(3)=ifrom(5)
            ifpod(4)=ifrom(5)
            ifrom(3)=ifrom(5)
            ifrom(4)=ifrom(5)
         endif
      endif

      ibuoy=0

      call finiparser_getdbl(d_out,'buoyancy:magnitude',ifnd)
      if (ifnd.eq.1) then
         gmag=d_out
         ibuoy=ibuoy+1
      endif

      call finiparser_getdbl(d_out,'buoyancy:angle',ifnd)
      if (ifnd.eq.1) then
         gangle=d_out
         ibuoy=ibuoy+2
         one=1.
         pi=4.*atan(one)
         gx=gmag*cos(gangle*pi/180)
         gy=gmag*sin(gangle*pi/180)
      endif
      
      call finiparser_getdbl(d_out,'buoyancy:gx',ifnd)
      if (ifnd.eq.1) then
         gx=d_out
         ibuoy=ibuoy+4
      endif

      call finiparser_getdbl(d_out,'buoyancy:gy',ifnd)
      if (ifnd.eq.1) then
         gy=d_out
         ibuoy=ibuoy+8
      endif

      call finiparser_getdbl(d_out,'buoyancy:gz',ifnd)
      if (ifnd.eq.1) then
         gz=d_out
         ibuoy=ibuoy+16
      endif

      if (ldim.eq.2) then
         if (ibuoy.ne.12.and.ibuoy.ne.3.and.ibuoy.ne.0) then
            write (6,*) 'invalid options for buoyancy(2D)'
            ierr=ierr+1
         endif
      endif

      if (ldim.eq.3) then
         if (ibuoy.ne.28.and.ibuoy.ne.0) then
            write (6,*) 'invalid options for buoyancy(3D)'
            ierr=ierr+1
         endif
      endif

      if (ibuoy.ne.0) ifbuoy=.true.

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
         else if (index(c_out,'OFF').eq.1) then
            ips='L2 '
            ifpb=.false.
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
            write (6,*) 'invalid option for pod:mode0 ',c_out
            ierr=ierr+1
         endif
      endif

      call finiparser_getbool(i_out,'pod:combined',ifnd)
      if (ifnd.eq.1) ifcomb=i_out.eq.1

      call finiparser_getdbl(d_out,'pod:ratio',ifnd)
      if (ifnd.eq.1) podrat=d_out

      call finiparser_getdbl(d_out,'pod:augment',ifnd)
      if (ifnd.eq.1) iaug=nint(d_out)

      ! QOI

      call finiparser_getdbl(d_out,'qoi:freq',ifnd)
      if (ifnd.eq.1) ad_qstep=nint(d_out)+ad_iostep*max(1-nint(d_out),0)

      call finiparser_getbool(i_out,'qoi:drag',ifnd)
      if (ifnd.eq.1) ifcdrag=i_out.eq.1

      call finiparser_getbool(i_out,'qoi:tke',ifnd)
      if (ifnd.eq.1) ifctke=i_out.eq.1

      call finiparser_getdbl(d_out,'qoi:nu',ifnd)
      if (ifnd.eq.1) inus=min(max(nint(d_out),0),6)

      call finiparser_getdbl(d_out,'qoi:nintp',ifnd)
      if (ifnd.eq.1) nintp=nint(d_out)

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
         if (index(c_out,'VT').eq.1) then
            icopt=0
         else if (index(c_out,'T').eq.1) then
            icopt=1
         else if (index(c_out,'V').eq.1) then
            icopt=2
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
            cfloc='NONE'
         else if (index(c_out,'CONV').eq.1) then
            cfloc='CONV'
         else if (index(c_out,'POST').eq.1) then
            cfloc='POST'
         else
            write (6,*) 'invalid option for filter:location ',c_out
            ierr=ierr+1
         endif
      endif

      call finiparser_getstring(c_out,'filter:type',ifnd)
      if (ifnd.eq.1) then
         call capit(c_out,132)
         if (index(c_out,'TFUNC').eq.1) then
            cftype='TFUN'
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
            cftype='DIFF'
            call finiparser_getdbl(d_out,'filter:radius',ifnd)
            if (ifnd.eq.1) then
               rdft=d_out
            else
               write (6,*) 'diff. filter needs a filter:radius value'
               ierr=ierr+1
            endif
         else if (index(c_out,'POLY').eq.1) then
            cftype='POLY'
            call finiparser_getdbl(d_out,'filter:radius',ifnd)
            if (ifnd.eq.1) then
               rdft=d_out
            else
               write (6,*) 'poly. filter needs a filter:radius value'
               ierr=ierr+1
            endif
         else if (index(c_out,'SMAG').eq.1) then
            cftype='SMAG'
            write (6,*) 'smag. filter type not supported ',c_out
            ierr=ierr+1
         else
            cftype='INVA'
            write (6,*) 'invalid option for filter:type ',c_out
            ierr=ierr+1
         endif
      endif

      ! EI

      call finiparser_getbool(i_out,'ei:mode',ifnd)
      if (ifnd.eq.1) ifei=i_out.eq.1

      call finiparser_getstring(c_out,'ei:eqn',ifnd)
      if (ifnd.eq.1) then
         call capit(c_out,132)
         if (index(c_out,'NS').eq.1) then
            eqn='NS  '
         else if (index(c_out,'SNS').eq.1) then
            eqn='SNS '
         else if (index(c_out,'VNS').eq.1) then
            eqn='VNS '
         else if (index(c_out,'VSNS').eq.1) then
            eqn='VSNS'
         else if (index(c_out,'POIS').eq.1) then
            eqn='POIS'
         else if (index(c_out,'HEAT').eq.1) then
            eqn='HEAT'
         else if (index(c_out,'ADVE').eq.1) then
            eqn='ADVE'
         else
            write (6,*) 'invalid option for ei:equation',c_out
            ierr=ierr+1
         endif
      endif

      ! TENSOR DECOMP
      call finiparser_getstring(c_out,'tendec:mode',ifnd)
      if (ifnd.eq.1) then
         call capit(c_out,132)
         if (index(c_out,'CP').eq.1) then
            ifcp=.true.
            if (ltr.lt.1) then
               write (6,*) 'increase CP rank ltr',ltr
               ierr=ierr+1
            endif
         else
            write (6,*) 'invalid option for tendec:mode ',c_out
            ierr=ierr+1
         endif
      endif

      call finiparser_getdbl(d_out,'tendec:rank',ifnd)
      if (ifnd.eq.1) ntr=min(nint(d_out),ltr)
      if (ntr.eq.0) ntr=ltr

      if (ntr.gt.ltr) then
         write (6,*) 'ntr > ltr is undefined configuration',ntr
         ierr=ierr+1
      endif

      call finiparser_getbool(i_out,'tendec:core',ifnd)
      if (ifnd.eq.1) ifcore=i_out.eq.1

      call finiparser_getbool(i_out,'tendec:skew',ifnd)
      if (ifnd.eq.1) ifquad=i_out.eq.1

      if (rmode.eq.'ON'.or.rmode.eq.'ONB') then
         call read_serial(rtmp1(1,1),1,'ops/nb ',rtmp2,-1)
         mb=nint(rtmp1(1,1))
         if (mb.lt.nb) then
            write (6,*) 'mb less than nb... ',mb
            ierr=ierr+1
         endif

         call read_serial(rtmp1(1,1),1,'ops/ns ',rtmp2,-1)
         nns=nint(rtmp1(1,1))
         if (nns.lt.ns) then
            write (6,*) 'ms less than ns... ',nns
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

      if (ierr.eq.0) call finiparser_dump()

      return
      end
c-----------------------------------------------------------------------
      subroutine bcastmpar

      ! Broadcast mor parameters to all processors

      include 'SIZE'
      include 'PARALLEL'
      include 'MOR'

      ! characters

      call bcast(rmode,csize*3)
      call bcast(eqn,csize*4)
      call bcast(ips,csize*3)
      call bcast(cfloc,csize*4)
      call bcast(cftype,csize*4)

      ! integers

      call bcast(ntr,isize)
      call bcast(nb,isize)
      call bcast(nbo,isize)
      call bcast(mb,isize)
      call bcast(ns,isize)
      call bcast(nns,isize)
      call bcast(nskip,isize)
      call bcast(navg_step,isize)
      call bcast(nplay,isize)
      call bcast(ad_qstep,isize)
      call bcast(inus,isize)
      call bcast(isolve,isize)
      call bcast(icopt,isize)
      call bcast(barr_func,isize)
      call bcast(ubarrseq,isize)
      call bcast(tbarrseq,isize)
      call bcast(nintp,isize)
      call bcast(iaug,isize)

      ! reals

      call bcast(rktol,wdsize)
      call bcast(box_tol,wdsize)
      call bcast(ubarr0,wdsize)
      call bcast(tbarr0,wdsize)
      call bcast(rbf,wdsize)
      call bcast(rdft,wdsize)
      call bcast(gx,wdsize)
      call bcast(gy,wdsize)
      call bcast(gz,wdsize)
      call bcast(podrat,wdsize)

      ! logicals

      call bcast(ifrecon,lsize)

      do i=0,ldimt1
         call bcast(ifpod(i),lsize)
         call bcast(ifrom(i),lsize)
      enddo

      call bcast(ifcomb,lsize)

      call bcast(ifei,lsize)
      call bcast(iftneu,lsize)
      call bcast(ifplay,lsize)
      call bcast(ifavg0,lsize)

      call bcast(ifcdrag,lsize)
      call bcast(ifctke,lsize)

      call bcast(iffastc,lsize)
      call bcast(iffasth,lsize)

      call bcast(ifforce,lsize)
      call bcast(ifsource,lsize)
      call bcast(ifbuoy,lsize)

      call bcast(ifpb,lsize)
      call bcast(ifdecpl,lsize)

      call bcast(ifcflow,lsize)

      call bcast(ifcp,lsize)
      call bcast(ifcore,lsize)
      call bcast(ifquad,lsize)
      call bcast(ifedvs,lsize)

      return
      END
c-----------------------------------------------------------------------
      subroutine mpar_verify(ierr)

      ! verifies the keys of the .mor file

      include 'MORDICT'

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
      subroutine set_morfle(morfle)

      ! set morfle to $session.mor

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      character*132 name
      character*1 sess1(132),path1(132),nam1(132)
      equivalence (session,sess1)
      equivalence (path,path1)
      equivalence (name,nam1)
      character*1 mor(4)
      character*4 mor4
      equivalence (mor,mor4)
      character*78  string

      character*132 morfle
      data mor4  /'.mor'/

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
