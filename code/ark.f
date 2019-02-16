c-----------------------------------------------------------------------
      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs) 
      
      INTEGER n,NMAX
      REAL h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50) Set to the maximum number of functions. USESderivs

c     Given values for n variables y and their derivatives dydx known at x, use the fifth-order Cash-Karp Runge-Kutta method to advance the solution over an interval h and return the incremented variables as yout. Also return an estimate of the local truncation er- ror in yout using the embedded fourth-order method. The user supplies the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
      INTEGER i
      REAL ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),
      ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51, B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3, DC4,DC5,DC6
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40., B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5, B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512., B63=575./13824.,B64=44275./110592.,B65=253./4096., C1=37./378.,C3=250./621.,C4=125./594.,C6=512./1771., DC1=C1-2825./27648.,DC3=C3-18575./48384., DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25)
      
c     The routine rkqs calls the routine rkck to take a Cash-Karp Runge-Kutta step:
      do i=1,n !First step.
         ytemp(i)=y(i)+B21*h*dydx(i)
      enddo

      call derivs(x+A2*h,ytemp,ak2) !Second step.

      do i=1,n
      ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
      enddo

      call derivs(x+A3*h,ytemp,ak3) !Third step.
      
      do i=1,n
         ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
      enddo

c     Sample page from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
c     Copyright (C) 1986-1992 by Cambridge University Press. Programs Copyright (C) 1986-1992 by Numerical Recipes Software.
c     Permission is granted for internet users to make one paper copy for their own personal use. Further reproduction, or any copying of machine- readable files (including this one) to any server computer, is strictly prohibited. To order Numerical Recipes books, diskettes, or CDROMs visit website http://www.nr.com or call 1-800-872-7423 (North America only), or send email to trade@cup.cam.ac.uk (outside North America).

      call derivs(x+A4*h,ytemp,ak4) !Fourth step.
      do i=1,n
         ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+ B54*ak4(i))
      enddo

      call derivs(x+A5*h,ytemp,ak5) !Fifth step.

      do i=1,n
         ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)
     $           +B64*ak4(i)+B65*ak5(i))
      enddo

      call derivs(x+A6*h,ytemp,ak6) !Sixth step.

      do i=1,n !Accumulate increments with proper weights.
         yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+ C6*ak6(i))
      enddo

      do i=1,n !Estimate error as difference between fourth and fifth order methods.
         yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i) +DC6*ak6(i))
      enddo

      return
      END
c-----------------------------------------------------------------------
      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)

      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.e-30)
c     Runge-Kutta driver with adaptive stepsize control. Integrate the starting values ystart(1:nvar) from x1 to x2 with accuracy eps, storing intermediate results in the common block /path/.
c     h1 should be set as a guessed first stepsize, hmin as the minimum allowed stepsize (can
c     be zero). On output nok and nbad are the number of good and bad (but retried and
c     Sample page from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
c     Copyright (C) 1986-1992 by Cambridge University Press. Programs Copyright (C) 1986-1992 by Numerical Recipes Software.
c     Permission is granted for internet users to make one paper copy for their own personal use. Further reproduction, or any copying of machine- readable files (including this one) to any server computer, is strictly prohibited. To order Numerical Recipes books, diskettes, or CDROMs visit website http://www.nr.com or call 1-800-872-7423 (North America only), or send email to trade@cup.cam.ac.uk (outside North America).
c      *
c     16.2 Adaptive Stepsize Control for Runge-Kutta 715
c     fixed) steps taken, and ystart is replaced by values at the end of the integration interval. derivs is the user-supplied subroutine for calculating the right-hand side derivative, while rkqs is the name of the stepper routine to be used. /path/ contains its own information about how often an intermediate value is to be stored.
      INTEGER i,kmax,kount,nstp
      REAL dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
      yp(NMAX,KMAXX),yscal(NMAX) COMMON /path/ kmax,kount,dxsav,xp,yp

c     User storage for intermediate results. Preset dxsav and kmax. x=x1

      h=sign(h1,x2-x1) nok=0
      nbad=0
      kount=0

      do i=1,nvar
         y(i)=ystart(i)
      enddo

      if (kmax.gt.0) xsav=x-2.*dxsav

c     Assures storage of first step. Take at most MAXSTP steps.
      do 16 nstp=1,MAXSTP
         call derivs(x,y,dydx)
         do 12 i=1,nvar
c     Scaling used to monitor accuracy. This general-purpose choice can be modified if need be.
            yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
         enddo 12
         if (kmax.gt.0) then
            if (abs(x-xsav).gt.abs(dxsav)) then
               if (kount.lt.kmax-1) then kount=kount+1
                  xp(kount)=x
                  do i=1,nvar
                     yp(i,kount)=y(i)
                  enddo
                  xsav=x
               endif
            endif
         endif
         if ((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
         call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
         if (hdid.eq.h) then
            nok=nok+1 else
            nbad=nbad+1
         endif

         if ((x-x2)*(x2-x1).ge.0.) then
            do i=1,nvar
               ystart(i)=y(i)
            enddo
            if (kmax.ne.0) then
               kount=kount+1
               xp(kount)=x
               do i=1,nvar
                  yp(i,kount)=y(i)
               enddo
            endif
            return
         endif
         if (abs(hnext).lt.hmin)
     $       pause ’stepsize smaller than minimum in odeint’
         h=hnext
      enddo
      pause ’too many steps in odeint’
      return
      END
c-----------------------------------------------------------------------
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      
      INTEGER n,NMAX
      REAL eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50) !Maximum number of equations.
C USESderivs,rkck
c     Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy and adjust stepsize. Input are the dependent variable vector y(1:n) and its derivative dydx(1:n) at the starting value of the independent variable x. Also input are the stepsize to be attempted htry, the required accuracy eps, and the vector yscal(1:n) against which the error is scaled. On output, y and x are replaced by their new values, hdid is the stepsize that was actually accomplished, and hnext is the estimated next stepsize. derivs is the user-supplied subroutine that computes the right-hand side derivatives.
      INTEGER i

      REAL errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,
     $     PSHRNK,ERRCON

      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)

c     The value ERRCON equals (5/SAFETY)**(1/PGROW), see use below.
      h=htry !Set stepsize to the initial trial value.

    1 call rkck(y,dydx,n,x,h,ytemp,yerr,derivs) !Take a step.

      Permission is granted for internet users to make one paper copy for their own personal use. Further reproduction, or any copying of machine- readable files (including this one) to any server computer, is strictly prohibited. To order Numerical Recipes books, diskettes, or CDROMs visit website http://www.nr.com or call 1-800-872-7423 (North America only), or send email to trade@cup.cam.ac.uk (outside North America).

      errmax=0. !Evaluate accuracy.
      do i=1,n
         errmax=max(errmax,abs(yerr(i)/yscal(i)))
      enddo
      errmax=errmax/eps !Scale relative to required tolerance. 
      if (errmax.gt.1.) then !Truncation error too large, reduce stepsize.
         htemp=SAFETY*h*(errmax**PSHRNK)
         h=sign(max(abs(htemp),0.1*abs(h)),h) !No more than a factor of 10. 
         xnew=x+h
         if (xnew.eq.x) pause ’stepsize underflow in rkqs’
         goto 1
      else
         if (errmax.gt.errcon) then
            hnext=safety*h*(errmax**pgrow)
         else ! No more than a factor of 5 increase.
            hnext=5.*h
         endif
         hdid=h
         x=x+h
         do i=1,n
            y(i)=ytemp(i)
         enddo
         return
      endif
      END
c-----------------------------------------------------------------------
