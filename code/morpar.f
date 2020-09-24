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
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'ADJOINT'
      INCLUDE 'RESTART'
      INCLUDE 'PARALLEL'
      INCLUDE 'CTIMER'
      INCLUDE 'TSTEP'

      character*132 c_out,txt,txt2

      call finiparser_load(morfle,ierr)
      if (ierr .ne. 0) return

      call mor_verify(ierr)
      if (ierr .ne. 0) return

      ! set parameters
      call finiparser_getDbl(d_out,'general:loglevel',ifnd)
      if (ifnd .eq. 1) loglevel=d_out

      call finiparser_getDbl(d_out,'general:optlevel',ifnd)
      if (ifnd .eq. 1) optlevel=d_out

      call finiparser_getString(c_out,'general:stopAt',ifnd)
      if (ifnd .eq. 1) then
         call capit(c_out,132)
         if (index(c_out,'ENDTIME') .eq. 1) then
            call finiparser_getDbl(d_out,'general:endTime',ifnd)
            if (ifnd .eq. 1) then
               param(10) = d_out
            else
               write (6,*) 'general:endTime'
               write (6,*) 'is required for general:stopAt = endTime!'
               goto 999
            endif
         else if (index(c_out,'NUMSTEPS') .eq. 1) then
            call finiparser_getDbl(d_out,'general:numSteps',ifnd)
            if (ifnd .eq. 1) then
               param(11) = d_out 
            else
               write (6,*) 'general:numSteps'
               write (6,*) 'is required for general:stopAt = numSteps!'
               goto 999
            endif
         else
            write (6,*) 'value: ',trim(c_out)
            write (6,*) 'is invalid for general:stopAt!'
            goto 999
         endif
      else
         call finiparser_getDbl(d_out,'general:numSteps',ifnd)
         if (ifnd .eq. 1) then 
            param(11) = d_out 
         else
            write (6,*) 'general:numSteps not found!'
            goto 999
         endif
      endif

      call finiparser_getDbl(d_out,'general:dt',ifnd)
      if (ifnd .eq. 1) then
         param(12) = d_out
      endif

      param(12) = -1*abs(param(12))
      call finiparser_getBool(i_out,'general:variableDt',ifnd)
      if (ifnd .eq. 1) then
         if (i_out .eq. 1) then
            param(12) = abs(param(12)) 
            call finiparser_getDbl(d_out,'general:targetCFL',ifnd)
            if (ifnd .eq. 1) then
               param(26) = d_out
            else
               write (6,*) 'general:targetCFL'
               write (6,*) 'is required for general:variableDt!'
               goto 999
            endif
         endif
      endif

      call finiparser_getDbl(d_out,'general:writeInterval',ifnd)
      if (ifnd .eq. 1) param(15) = d_out
      call finiparser_getString(c_out,'general:writeControl',ifnd)
      if (ifnd .eq. 1) then
         call capit(c_out,132)
         if (index(c_out,'RUNTIME') .eq. 1) then
            param(14) = d_out
            param(15) = 0
         else if (index(c_out,'TIMESTEP') .eq. 1) then
            param(14) = 0
            param(15) = d_out
         else if (index(c_out,'ELAPSEDTIME') .eq. 1) then
            param(14) = d_out
            param(15) = 0
            timeioe = 1
         else
            write (6,*) 'value: ',trim(c_out)
            write (6,*) 'is invalid for general:writeControl!'
            goto 999
         endif
      endif

      call finiparser_getDbl(d_out,'pressure:residualTol',ifnd)
      if (ifnd .eq. 1) param(21) = d_out 
      call finiparser_getDbl(d_out,'velocity:residualTol',ifnd)
      if (ifnd .eq. 1) then
        restol(1) = d_out 
        param(22) = d_out
      endif

      call finiparser_find(i_out,'temperature',ifnd)
      if (ifnd .eq. 1) then
        ifheat = .true.
        ifto   = .true.
        idpss(1) = 0 ! Helmholtz is default
      endif

      j = 0
      do i = 1,99
         write (txt,"('scalar',i2.2)") i
         call finiparser_find(i_out,txt,ifnd)
         if (ifnd .eq. 1) j = j + 1
      enddo
      if (j.gt.ldimt-1) then
         write (6,*) 'found more scalars than specified in SIZE!' 
         goto 999
      endif

      j = 0
      do i = 1,ldimt-1
         write (txt,"('scalar',i2.2)") i
         call finiparser_find(i_out,txt,ifnd)
         if (ifnd .eq. 1) then
            j = j + 1 
            ifpsco(i) = .true.
            idpss(i+1) = 0 ! Helmholtz is default
         endif
      enddo
      param(23) = j ! number of scalars 
      n = param(23)
 
      is = 2
      if (ifheat) is = 1 

      do i = is,n+1

      if (i.eq.1) then
        txt = 'temperature'
      else
        write (txt,"('scalar',i2.2)") i-1
      endif

      call finiparser_getString(c_out,trim(txt)//':solver',ifnd)
      call capit(c_out,132)
      if (ifnd .eq. 1) then 
        if (index(c_out,'CVODE') .eq. 1) then
          idpss(i) = 1
          call finiparser_getDbl(d_out,trim(txt)//':absoluteTol',ifnd)
          if (ifnd .eq. 1) then
             atol(i+1) = d_out 
          else
             write (6,*) trim(txt) // ':absoluteTol' 
             write (6,*) 'is required for ',trim(txt)//':solver = CVODE'
             goto 999
          endif
        else if (index(c_out,'HELM') .eq. 1) then
          continue
        else if (index(c_out,'NONE') .eq. 1) then
          idpss(i) = -1
        else
           write (6,*) 'value: ',trim(c_out)
           write (6,*) 'is invalid for ',trim(txt)//':solver!' 
           goto 999
        endif
      else
        call finiparser_getDbl(d_out,trim(txt)//':residualTol',ifnd)
        if (ifnd .eq. 1) then
           restol(i+1) = d_out 
        endif
      endif
      call finiparser_getBool(i_out,trim(txt)//':residualProj',ifnd)
      if (ifnd .eq. 1) then
         ifprojfld(i+1) = .false.
         if (i_out .eq. 1) ifprojfld(i+1) = .true.
      endif
 
      enddo

      call finiparser_getDbl(d_out,'cvode:absoluteTol',ifnd)
      if (ifnd .eq. 1) param(162) = d_out 
      call finiparser_getDbl(d_out,'cvode:relativeTol',ifnd)
      if (ifnd .eq. 1) param(163) = d_out 
      call finiparser_getDbl(d_out,'cvode:dtmax',ifnd)
      if (ifnd .eq. 1) param(164) = d_out 
      call finiparser_getDbl(d_out,'cvode:DQJincrementFactor',ifnd)
      if (ifnd .eq. 1) param(165) = d_out 
      call finiparser_getDbl(d_out,'cvode:ratioLNLtol',ifnd)
      if (ifnd .eq. 1) param(166) = d_out 

      call finiparser_getString(c_out,'cvode:preconditioner',ifnd)
      call capit(c_out,132)
      if (ifnd .eq. 1) then 
        if (index(c_out,'USER') .eq. 1) then
           param(167) = 1
        else if (index(c_out,'NONE') .eq. 1) then
           param(167) = 0
        else
           write (6,*) 'value: ',trim(c_out)
           write (6,*) 'is invalid for cvode:preconditioner!' 
           goto 999
        endif
      endif

      j = 0
      do i = 1,ldimt
         if (idpss(i).ge.0) j = j + 1
      enddo
      if (j .ge. 1) then ! we have to solve for temp and/or ps
         ifheat = .true.  
      else
         ifheat = .false.
      endif

      call finiparser_getDbl(d_out,'magnetic:viscosity',ifnd)
      if (ifnd .eq. 1) param(29) = d_out 
      if (param(29).lt.0.0) param(29) = -1.0/param(29)

      call finiparser_getDbl(d_out,'mesh:numberOfBCFields',ifnd)
      if (ifnd .eq. 1) param(32) = int(d_out)

      call finiparser_getDbl(d_out,'mesh:firstBCFieldIndex',ifnd)
      if (ifnd .eq. 1) param(33) = int(d_out)

      call finiparser_getString(c_out,'pressure:preconditioner',ifnd)
      if (ifnd .eq. 1) then 
         call capit(c_out,132)
         if (index(c_out,'SEMG_XXT') .eq. 1) then
            param(40) = 0
        else if (index(c_out,'SEMG_AMG_HYPRE') .eq. 1) then
            param(40) = 2
         else if (index(c_out,'SEMG_AMG') .eq. 1) then
            param(40) = 1
         else if (index(c_out,'FEM_AMG_HYPRE') .eq. 1) then
            param(40) = 3
         else
           write (6,*) 'value: ',trim(c_out)
           write (6,*) 'is invalid for pressure:preconditioner!'
           goto 999
         endif
      endif 

      call finiparser_getBool(i_out,'general:writeDoublePrecision',ifnd)
      if (ifnd .eq. 1 .and. i_out .eq. 1) param(63) = 1 

      call finiparser_getDbl(d_out,'general:writeNFiles',ifnd)
      if (ifnd .eq. 1) param(65) = int(d_out) 

      call finiparser_getBool(i_out,'velocity:residualProj',ifnd)
      if (ifnd .eq. 1) then
        ifprojfld(1) = .false.
        if (i_out .eq. 1) ifprojfld(1) = .true. 
      endif

      call finiparser_getBool(i_out,'pressure:residualProj',ifnd)
      if (ifnd .eq. 1) then
        param(95) = 0
        if (i_out .eq. 1) param(95) = 5 
      endif

      call finiparser_getBool(i_out,'general:dealiasing',ifnd)
      if (ifnd .eq. 1 .and. i_out .eq. 0) param(99) = -1 

c     filtering parameters
      call finiparser_getString(c_out,'general:filtering',ifnd)
      if (ifnd .eq. 1) then
c        stabilization type: none, explicit or hpfrt    
         call capit(c_out,132)
         if (index(c_out,'NONE') .eq. 1) then
            filterType = 0
            goto 101
         else if (index(c_out,'EXPLICIT') .eq. 1) then
            filterType = 1
            call ltrue(iffilter,size(iffilter))
         else if (index(c_out,'HPFRT') .eq. 1) then
            filterType = 2
            call ltrue(iffilter,size(iffilter))
         else
           write (6,*) 'value: ',c_out
           write (6,*) 'is invalid for general:filtering!'
           goto 999
         endif
         call finiparser_getDbl(d_out,'general:filterWeight',ifnd)
         if (ifnd .eq. 1) then
            param(103) = d_out 
         else
            write (6,*) 'general:filterWeight'
            write (6,*) 'is required for general:filtering!'
            goto 999
         endif
         call finiparser_getDbl(d_out,'general:filterCutoffRatio',ifnd)
         if (ifnd .eq. 1) then
            dtmp = anint(lx1*(1.0 - d_out)) 
            param(101) = max(dtmp-1,0.0)
            if (abs(1.0 - d_out).lt.0.01) filterType = 0
         else
            write (6,*) 'general:filterCutoffRatio'
            write (6,*) 'is required for general:filtering!'
            goto 999
         endif
 101     continue 
      endif

      call finiparser_getString(c_out,'cvode:mode',ifnd)
      call capit(c_out,132)
      if (index(c_out,'NORMAL') .eq. 1) param(160) = 1
      if (index(c_out,'NORMAL_TSTOP' ) .eq. 1) param(160) = 3
 
      do i = 1,20
         call blank(txt,132)
         write (txt,"('general:userParam',i2.2)") i
         call finiparser_getDbl(d_out,txt,ifnd)
         if (ifnd .eq. 1) uparam(i) = d_out
      enddo

c set logical flags
      call finiparser_getString(c_out,'general:timeStepper',ifnd)
      if (ifnd .eq. 1) then
        call capit(c_out,132)

        if (index(c_out,'BDF1') .eq. 1) then
           param(27) = 1 
        else if (index(c_out,'BDF2') .eq. 1) then
           param(27) = 2 
        else if (index(c_out,'BDF3') .eq. 1) then
           param(27) = 3 
        else
           write (6,*) 'value: ',trim(c_out)
           write (6,*) 'is invalid for general:timeStepper!'
           goto 999
        endif
      endif

      call finiparser_getString(c_out,'general:extrapolation',ifnd)
      if (ifnd .eq. 1) then
        call capit(c_out,132)
        if (index(c_out,'OIFS') .eq. 1) then
           ifchar = .true.

           call finiparser_getDbl(d_out,'general:targetCFL',ifnd)
           if (ifnd .eq. 1) then
              param(26) = d_out
           else
              write (6,*) 'general:targetCFL'
              write (6,*) 'is required for general:extrapolation!'
              goto 999
           endif
        else if (index(c_out,'STANDARD') .eq. 1) then
           continue
        else
           write (6,*) 'value: ',trim(c_out)
           write (6,*) 'is invalid for general:extrapolation!'
           goto 999
        endif
      endif

      call finiparser_find(i_out,'velocity',ifnd)
      if (ifnd .eq. 1) then
        ifflow = .true.
        ifvo   = .true.
        ifpo   = .true.
      endif

      call finiparser_getString(c_out,'mesh:motion',ifnd)
      if (ifnd .eq. 1) then
       call capit(c_out,132)
       if (index(c_out,'ELASTICITY') .eq. 1) then
          ifmvbd = .true.
          call finiparser_getDbl(d_out,'mesh:viscosity',ifnd)
          if (ifnd .eq. 1) param(47) = d_out
          call finiparser_getBool(i_out,'mesh:residualProj',ifnd)
          if (ifnd .eq. 1) then
             ifprojfld(0) = .false.
             if (i_out .eq. 1) ifprojfld(0) = .true. 
          endif
        else if (index(c_out,'USER') .eq. 1) then
          ifmvbd = .true.
          ifusermv = .true.
        else if (index(c_out,'NONE') .eq. 1) then
          continue
        else
          write (6,*) 'value: ',trim(c_out)
          write (6,*) 'is invalid for mesh:motion!'
          goto 999
        endif
      endif
      call finiparser_getDbl(d_out,'mesh:residualTol',ifnd)
      if (ifnd .eq. 1) restol(0) = d_out 

      call finiparser_getBool(i_out,'problemType:axiSymmetry',ifnd)
      if (ifnd .eq. 1) then
        ifaxis = .false.
        if (i_out .eq. 1) ifaxis = .true.
      endif

      call finiparser_getBool(i_out,'problemType:swirl',ifnd)
      if (ifnd .eq. 1) then
        ifaziv = .false.
        if (i_out .eq. 1) ifaziv = .true.
      endif

      call finiparser_getBool(i_out,'problemType:cyclicBoundaries',ifnd)
      if (ifnd.eq.1) then
        ifcyclic = .false.
        if (i_out .eq. 1) ifcyclic = .true.
      endif

      call finiparser_getString(c_out,'problemType:equation',ifnd)
      call capit(c_out,132)
      if (index(c_out,'STEADYSTOKES').eq.1) then
         iftran = .false.
      else if (index(c_out,'INCOMPNS').eq.1) then
         continue
      else if (index(c_out,'LOWMACHNS').eq.1) then
         iflomach = .true.
      else if (index(c_out,'INCOMPLINNS').eq.1 .or.
     $         index(c_out,'INCOMPLINADJNS').eq.1) then
         ifpert = .true.
         if (index(c_out,'INCOMPLINADJNS').eq.1) ifadj  = .true.
         call finiparser_getDbl
     $        (d_out,'problemType:numberOfPerturbations',ifnd)
         if (ifnd .eq. 1) then
            param(31) = int(d_out) 
         else
            param(31) = 1 
         endif
         call finiparser_getBool
     $        (i_out,'problemType:solveBaseFlow',ifnd)
         if (ifnd .eq. 1) then
            ifbase = .false.
            if (i_out .eq. 1) ifbase = .true.
         else
            write (6,*) 'problemType:solveBaseFlow'
            write (6,*) 'is required for ', trim(c_out) 
            goto 999
         endif
      else if (index(c_out,'COMPNS') .eq. 1) then
#ifdef CMTNEK
         continue
#else
         write (6,*) 'value: ',trim(c_out)
         write (6,*) 'not supported for problemType:equation!'
         write (6,*) 'Recompile with CMTNEK ...'
         goto 999
#endif
      else if (index(c_out,'INCOMPMHD') .eq. 1) then
         write (6,*) 'value: ',trim(c_out)
         write (6,*) 'not yet supported for problemType:equation!'
         goto 999
      endif

      call finiparser_getBool(i_out,
     &                        'problemType:stressFormulation',ifnd)
      if (ifnd .eq. 1) then
        ifstrs = .false.
        if (i_out .eq. 1) ifstrs = .true.
      endif

      call finiparser_getBool(i_out,
     &                        'problemType:variableProperties',ifnd)
      if (ifnd .eq. 1) then
        ifuservp = .false.
        if (i_out .eq. 1) ifuservp = .true.
      endif

      call finiparser_getBool(i_out,'problemType:dp0dt',ifnd)
      if (ifnd .eq. 1) then
        if (i_out .eq. 1) ifdp0dt = .true.
      endif

      call finiparser_getBool(i_out,'cvode:stiff',ifnd)
      if (ifnd .eq. 1) then
        param(161) = 1 ! AM
        if (i_out .eq. 1) param(161) = 2 !BDF
      endif

c set advection
      call finiparser_getBool(i_out,'velocity:advection',ifnd)
      if (ifnd .eq. 1) then
        ifadvc(1) = .false.
        if (i_out .eq. 1) ifadvc(1) = .true.
      endif
 
      call finiparser_getBool(i_out,'temperature:advection',ifnd)
      if (ifnd .eq. 1) then
        ifadvc(2) = .false.
        if (i_out .eq. 1) ifadvc(2) = .true.
      endif

      do i = 1,ldimt-1
         write (txt,"('scalar',i2.2,a)") i,':advection'
         call finiparser_getBool(i_out,txt,ifnd)
         if (ifnd .eq. 1) then
           ifadvc(i+2) = .false.
           if (i_out .eq. 1) ifadvc(i+2) = .true.
         endif
      enddo

c set mesh-field mapping
      call finiparser_getBool(i_out,'temperature:conjugateHeatTransfer',
     &                        ifnd)
      if (ifnd .eq. 1) then
        if (i_out .eq. 1) then
          iftmsh(0) = .true.
          iftmsh(2) = .true.
        endif
      endif

      do i = 1,ldimt-1
         write (txt,"('scalar',i2.2,a)") i,':conjugateHeatTransfer'
         call finiparser_getBool(i_out,txt,ifnd)
         if (ifnd .eq. 1) then
           iftmsh(i+2) = .false.
           if (i_out .eq. 1) iftmsh(i+2) = .true.
         endif
      enddo

c set output flags
      call finiparser_getBool(i_out,'mesh:writeToFieldFile',ifnd) 
      if (ifnd .eq. 1) then
        ifxyo = .false.
        if (i_out .eq. 1) ifxyo = .true.
      endif

      call finiparser_getBool(i_out,'velocity:writeToFieldFile',ifnd)
      if (ifnd .eq. 1) then
        ifvo = .false.
        if (i_out .eq. 1) ifvo = .true.
      endif

      call finiparser_getBool(i_out,'pressure:writeToFieldFile',ifnd)
      if (ifnd .eq. 1) then
        ifpo = .false.
        if (i_out .eq. 1) ifpo = .true.
      endif

      call finiparser_getBool(i_out,'temperature:writeToFieldFile',ifnd)
      if (ifnd .eq. 1) then
        ifto = .false.
        if (i_out .eq. 1) ifto = .true.
      endif

      do i = 1,ldimt-1
         write (txt,"('scalar',i2.2,a)") i,':writeToFieldFile'
         call finiparser_getBool(i_out,txt,ifnd)
         if (ifnd .eq. 1) then
           ifpsco(i) = .false.
           if (i_out .eq. 1) ifpsco(i) = .true.
         endif
      enddo

c set properties
      call finiparser_getDbl(d_out,'velocity:viscosity',ifnd)
      if (ifnd .eq. 1) cpfld(1,1) = d_out 
      if (cpfld(1,1) .lt.0.0) cpfld(1,1)  = -1.0/cpfld(1,1)
      call finiparser_getDbl(d_out,'velocity:density',ifnd)
      if (ifnd .eq. 1) cpfld(1,2) = d_out 

      call finiparser_getDbl(d_out,'temperature:conductivity',ifnd)
      if (ifnd .eq. 1) cpfld(2,1) = d_out 
      if (cpfld(2,1) .lt.0.0) cpfld(2,1)  = -1.0/cpfld(2,1)
      call finiparser_getDbl(d_out,'temperature:rhoCp',ifnd)
      if (ifnd .eq. 1) cpfld(2,2) = d_out 

      do i = 1,ldimt-1
         write (txt,"('scalar',i2.2,a)") i,':diffusivity'
         call finiparser_getDbl(d_out,txt,ifnd)
         if (ifnd .eq. 1) cpfld(2+i,1) = d_out 
         if (cpfld(2+i,1) .lt.0.0) cpfld(2+i,1)  = -1.0/cpfld(2+i,1)
         write (txt,"('scalar',i2.2,a)") i,':density'
         call finiparser_getDbl(d_out,txt,ifnd)
         if (ifnd .eq. 1) cpfld(2+i,2) = d_out 
      enddo

c set restart options
      call finiparser_findTokens('general:startFrom', ',' , ifnd)
      do i = 1,min(ifnd,15)
         call finiparser_getToken(initc(i),i)
         if (index(initc(i),'0') .eq. 1) call blank(initc(i),132)
      enddo


100   if (ierr.eq.0) call finiparser_dump()
      return

c error handling
 999  continue
      ierr = 1
      goto 100

      end
c-----------------------------------------------------------------------
      subroutine bcastrpar
C
C     Broadcast run parameters to all processors
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

         do j = 1,MORDICT_NKEYS ! do we find the key in the par-dictionary 
            if (index(pardictkey(j),key).eq.1) goto 10      

            is = index(key,'SCALAR')
            if (is .eq. 1) then
              call chcopy(txt,key,132)
              call chcopy(tx1(is+6),'%%',2) 
              if (index(pardictkey(j),txt).eq.1) goto 10
            endif
         enddo
         write (6,*) 'ERROR: Par file contains unknown key ', key
         ierr = ierr + 1
   10 enddo

      return
      end
