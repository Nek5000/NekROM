c-----------------------------------------------------------------------
      subroutine get_saved_fields(usave,vsave,wsave,nsave,u0)

c     This routine reads files specificed in file.list

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      parameter (lt=lx1*ly1*lz1*lelt)
      real usave(lt,nsave),vsave(lt,nsave),wsave(lt,nsave)
      real u0(lt,3) ! Initial condtion


      ierr = 0
      if (nid.eq.0) open(77,file='file.list',status='old',err=199)
      ierr = iglmax(ierr,1)
      if (ierr.gt.0) goto 199
      n = lx1*ly1*lz1*nelt
      n2= lx2*ly2*lz2*nelt

      icount = 0
      do ipass=1,nsave

         call blank(initc,127)
         initc(1) = 'done '
         if (nid.eq.0) read(77,127,end=998) initc(1)
  998    call bcast(initc,127)
  127    format(a127)

         if (indx1(initc,'done ',5).eq.0) then ! We're not done
            nfiles = 1
            call restart(nfiles)  ! Note -- time is reset.

!           Usave = U_snapshot - U_stokes:

            call opsub3 (usave(1,ipass),vsave(1,ipass),wsave(1,ipass)
     $                  ,vx,vy,vz,u0(1,1),u0(1,2),u0(1,3))

            icount = icount+1
         else
            goto 999
         endif

      enddo

  999 continue  ! clean up averages
      if (nid.eq.0) close(77)

      nsave = icount ! Actual number of files read

      return

  199 continue ! exception handle for file not found
      ierr = 1
      if (nid.eq.0) ierr = iglmax(ierr,1)
      call exitti('Auto averager did not find list file.$',ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine genmodes

      include 'SIZE'
      include 'TOTAL'
      include 'POD'

      parameter (lt=lx1*ly1*lz1*lelt)

      real usave(lt,ms),vsave(lt,ms),wsave(lt,ms)
      real uu(ms,ms),Identity(ms,ms),eig(ms),eigv(ms,ms),w(ms,ms)
      real uw(lt),vw(lt),ww(lt),h1(lt),h2(lt)
      real u0(lt,3)
      real u0r(ms)

      real evec(ms,nb)

      character(len=10) fname

      if (nio.eq.0) write (6,*) 'inside genmodes'

      n  = lx1*ly1*lz1*nelt
      ns = ms

      call rzero(vz,n)
      call rzero(wb,n)

      call opcopy(u0(1,1),u0(1,2),u0(1,3),ub(1,0),vb(1,0),wb(1,0))

      call get_saved_fields(usave,vsave,wsave,ns,u0)

      call rone (h1,n)
      call rzero(h2,n)
      call rzero(Identity,ms*ms)

      do j=1,ns                    ! Form the Gramian, U=U_K^T A U_K
         call axhelm(uw,usave(1,j),h1,h2,1,1)
         call axhelm(vw,vsave(1,j),h1,h2,1,1)
         if (ldim.eq.3) call axhelm(ww,wsave(1,j),h1,h2,1,1)
         do i=1,ns
            uu(i,j) = glsc2(usave(1,i),uw,n)+glsc2(vsave(1,i),vw,n)
            if (ldim.eq.3) uu(i,j) = uu(i,j)+glsc2(wsave(1,i),ww,n)
         enddo
         identity(j,j) = 1
         if (nio.eq.0) write(6,*) j,uu(1,j),' uu'
      enddo

      call generalev(uu,Identity,eig,ms,w)
      call copy(eigv,uu,ms*ms)
      eig = eig(ms:1:-1)

      nvecs = nb
      if (nio.eq.0) write(6,*)'number of mode:',nb

      do l = 1,nvecs
         call copy(evec(1,l),eigv(1,ms-l+1),ms)
      enddo

      ONE = 1.
      ZERO= 0.

      ! ub, vb, wb, are the modes
      call dgemm( 'N','N',n,nb,ms,ONE,usave,lt,evec,ms,ZERO,ub(1,1),lt)
      call dgemm( 'N','N',n,nb,ms,ONE,vsave,lt,evec,ms,ZERO,vb(1,1),lt)
      if (ldim.eq.3)
     $call dgemm( 'N','N',n,nb,ms,ONE,wsave,lt,evec,ms,ZERO,wb(1,1),lt)

      do i=0,nb
         call outpost(ub(1,i),vb(1,i),wb(1,i),pr,t,'bas')
      enddo

      if (nio.eq.0) write (6,*) 'exiting genmodes'

      return
      end
c-----------------------------------------------------------------------
