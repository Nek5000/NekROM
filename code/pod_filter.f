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
      subroutine pod_df(uu,delta)

      return
      end
c-----------------------------------------------------------------------
