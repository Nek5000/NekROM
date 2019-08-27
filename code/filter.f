c-----------------------------------------------------------------------
      subroutine pod_proj(uu,r1)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real uu(nb)
      integer r1

c     if (nid.eq.0) write(6,*)'r1:',r1,'nb-r1',nb-r1

c     write(6,*)'quick check before'
c     do ii=1,nb
c     write(6,*)ii,uu(ii),u(ii,1)
c     enddo

      call cfill(uu(r1+1),0.0,nb-r1)

c     write(6,*)'quick check after'
c     do ii=1,nb
c     write(6,*)ii,uu(ii),u(ii,1)
c     enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine pod_df(uu)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      integer icalld
      save    icalld
      data    icalld /0/

      real dfops(nb,nb)
      real tmp(nb)

      if (icalld.eq.0) then
         call setdf(dfops,au,bu,rdft*rdft)
         call dgetrf(nb,nb,dfops,lub,ipiv,info)
      endif

      call mxm(bu,nb,uu,nb,tmp,1)
      call dgetrs('N',nb,1,dfops,lub,ipiv,tmp,nb,info)
      call copy(uu,tmp,nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine setdf(flu,a,b,ad_diff)

      include 'SIZE'
      include 'MOR'


      real flu(nb,nb),a(nb,nb),b(nb,nb)
      real ad_diff

      call cmult2(flu,a,ad_diff,nb*nb)
      call add2(flu,b,nb*nb)
         
      return
      end
c-----------------------------------------------------------------------
