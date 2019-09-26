c-----------------------------------------------------------------------
      subroutine pod_proj(uu,r1)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real uu(nb)
      integer r1

      call cfill(uu((nb-r1+1):nb),0.0,r1)

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

      real tmp(nb)

      if (icalld.eq.0) then
         call setdf(dfops,au,bu,rdft*rdft)
         call dgetrf(nb,nb,dfops,nb,ipiv,info)
         icalld=1
      endif

      if (rdft.le.1e-12) then
      else
         call mxm(bu,nb,uu,nb,tmp,1)
         call dgetrs('N',nb,1,dfops,nb,ipiv,tmp,nb,info)
         call copy(uu,tmp,nb)
      endif

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
