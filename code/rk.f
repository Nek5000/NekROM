c-----------------------------------------------------------------------
      subroutine rk_step(un,ue,u,t,dt,k,wk,n)

      common /rk_btab/ a(6,6),b(6),c(6),d(6),bd(6)
      common /rk_ivar/ ns
      common /rk_rvar/ tol

      real k(nb,6),wk(n)
      real u(n),un(n),ue(n)

      call evf(t,u,k(1,1))
      call rk_helper(u,t,dt,k,wk,n)
      call mxm(k,n,b,ns,un,1)
      call add2s1(un,u,dt,n)

      call mxm(k,n,bd,ns,ue,1)
      call cmult(ue,dt,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine rk_helper(u,t,dt,k,wk,n)

      common /rk_btab/ a(6,6),b(6),c(6),d(6),bd(6)
      common /rk_ivar/ ns

      real u(n),k(n,6),wk(n)

      do i=2,ns
         call mxm(k,n,a(1,i),i,wk,1)
         call add2s1(wk,u,dt,n)
         call evf(t+c(i)*dt,wk,k(1,i),1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rk1_setup

      common /rk_btab/ a(6,6),b(6),c(6),d(6),bd(6)
      common /rk_ivar/ ns
      common /rk_rvar/ tol

      c(1)=0.0
      b(1)=1.0
      ns=1
      tol=0.

      return
      end
c-----------------------------------------------------------------------
      subroutine rkmp_setup

      common /rk_btab/ a(6,6),b(6),c(6),d(6),bd(6)
      common /rk_ivar/ ns
      common /rk_rvar/ tol

      a(1,2)= 1.0 / 2.0

      c(1)= 0.0
      c(2)= 1.0 / 2.0

      b(1)= 0.0
      b(2)= 1.0

      ns=2
      tol=0.

      return
      end
c-----------------------------------------------------------------------
      subroutine rk4_setup

      common /rk_btab/ a(6,6),b(6),c(6),d(6),bd(6)
      common /rk_ivar/ ns
      common /rk_rvar/ tol

      a(1,2)= 1.0 / 2.0

      a(1,3)= 0.0
      a(2,3)= 1.0 / 2.0

      a(1,4)= 0.0
      a(2,4)= 0.0
      a(3,4)= 1.0

      c(1)= 0.0
      c(2)= 1.0 / 2.0
      c(3)= 1.0 / 2.0
      c(4)= 1.0

      b(1)= 1.0/6.0
      b(2)= 1.0/3.0
      b(3)= 1.0/3.0
      b(4)= 1.0/6.0

      ns=4
      tol=0.

      return
      end
c-----------------------------------------------------------------------
      subroutine rk38_setup

      common /rk_btab/ a(6,6),b(6),c(6),d(6),bd(6)
      common /rk_ivar/ ns
      common /rk_rvar/ tol

      a(1,2)= 1.0 / 3.0

      a(1,3)= -1.0 / 3.0
      a(2,3)=  1.0

      a(1,4)=  1.0
      a(2,4)= -1.0
      a(3,4)=  1.0

      c(1)= 0.0
      c(2)= 1.0 / 3.0
      c(3)= 2.0 / 3.0
      c(4)= 1.0

      b(1)= 1.0/8.0
      b(2)= 3.0/8.0
      b(3)= 3.0/8.0
      b(4)= 1.0/8.0

      ns=4
      tol=0.

      return
      end
c-----------------------------------------------------------------------
      subroutine rkck_setup

      common /rk_btab/ a(6,6),b(6),c(6),d(6),bd(6)
      common /rk_ivar/ ns
      common /rk_rvar/ tol

      c(1)= 0.0
      c(2)= 1.0 /  5.0
      c(3)= 3.0 / 10.0
      c(4)= 3.0 /  5.0
      c(5)= 1.0
      c(6)= 7.0 /  8.0

      a(1,2)= 1.0 / 5.0

      a(1,3)= 3.0 / 40.0
      a(2,3)= 9.0 / 40.0

      a(1,4)=  3.0 / 10.0
      a(2,4)= -9.0 / 10.0
      a(3,4)=  6.0 /  5.0

      a(1,5)= -11.0 / 54.0
      a(2,5)=   5.0 / 2.0
      a(3,5)= -70.0 / 27.0
      a(4,5)=  35.0 / 27.0

      a(1,6)=  1631.0 /  55296.0
      a(2,6)=   175.0 /    512.0
      a(3,6)=   575.0 /  13824.0
      a(4,6)= 44275.0 / 110592.0
      a(5,6)=   253.0 /   4096.0

      b(1)=  37.0 /  378.0
      b(2)=   0.0
      b(3)= 250.0 /  621.0
      b(4)= 125.0 /  594.0
      b(5)=   0.0
      b(6)= 512.0 / 1771.0

      d(1)=  2825.0  / 27648.0
      d(2)=     0.0
      d(3)= 18575.0  / 48384.0
      d(4)= 13525.0  / 55296.0
      d(5)=   277.0  / 14336.0
      d(6)=     1.0  /     4.0

      call sub3(bd,b,d,6)

      ns=6
      tol=0.

      return
      end
c-----------------------------------------------------------------------
