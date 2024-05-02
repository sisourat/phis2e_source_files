
        SUBROUTINE CCHG(A,B,Z,CHG)
!
!       ===================================================
!       Purpose: Compute confluent hypergeometric function
!                M(a,b,z) with real parameters a, b and a
!                complex argument z
!       Input :  a --- Parameter
!                b --- Parameter
!                z --- Complex argument
!       Output:  CHG --- M(a,b,z)
!       Routine called: GAMMA for computing gamma function
!       ===================================================
!
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX *16 (C,Z)
        CHW=0.0D0
        PI=3.141592653589793D0
        CI=(0.0D0,1.0D0)
        A0=A
        A1=A
        Z0=Z
        IF (B.EQ.0.0.OR.B.EQ.-INT(ABS(B))) THEN
           CHG=(1.0D+300,0.0D0)
        ELSE IF (A.EQ.0.0D0.OR.CDABS(Z).LT.1.0D-15) THEN
           CHG=(1.0D0,0.0D0)
        ELSE IF (A.EQ.-1.0D0) THEN
           CHG=1.0D0-Z/B
        ELSE IF (A.EQ.B) THEN
           CHG=CDEXP(Z)
        ELSE IF (A-B.EQ.1.0D0) THEN
           CHG=(1.0D0+Z/B)*CDEXP(Z)
        ELSE IF (A.EQ.1.0D0.AND.B.EQ.2.0D0) THEN
           CHG=(CDEXP(Z)-1.0D0)/Z
        ELSE IF (A.EQ.INT(A).AND.A.LT.0.0D0) THEN
           M=INT(-A)
           CR=(1.0D0,0.0D0)
           CHG=(1.0D0,0.0D0)
           DO 10 K=1,M
              CR=CR*(A+K-1.0D0)/K/(B+K-1.0D0)*Z
10            CHG=CHG+CR
        ELSE
           X0=REAL(Z)
           IF (X0.LT.0.0D0) THEN
              A=B-A
              A0=A
              Z=-Z
           ENDIF
           IF (A.LT.2.0D0) NL=0
           IF (A.GE.2.0D0) THEN
              NL=1
              LA=INT(A)
              A=A-LA-1.0D0
           ENDIF
           DO 30 N=0,NL
              IF (A0.GE.2.0D0) A=A+1.0D0
              IF (CDABS(Z).LT.20.0D0+ABS(B).OR.A.LT.0.0D0) THEN
                 CHG=(1.0D0,0.0D0)
                 CRG=(1.0D0,0.0D0)
                 DO 15 J=1,500
                    CRG=CRG*(A+J-1.0D0)/(J*(B+J-1.0D0))*Z
                    CHG=CHG+CRG
                    IF (CDABS((CHG-CHW)/CHG).LT.1.D-15) GO TO 25
                    CHW=CHG
15               CONTINUE
              ELSE
                 CALL GAMMA(A,G1)
                 CALL GAMMA(B,G2)
                 BA=B-A
                 CALL GAMMA(BA,G3)
                 CS1=(1.0D0,0.0D0)
                 CS2=(1.0D0,0.0D0)
                 CR1=(1.0D0,0.0D0)
                 CR2=(1.0D0,0.0D0)
                 DO 20 I=1,8
                    CR1=-CR1*(A+I-1.0D0)*(A-B+I)/(Z*I)
                    CR2=CR2*(B-A+I-1.0D0)*(I-A)/(Z*I)
                    CS1=CS1+CR1
20                  CS2=CS2+CR2
                 X=REAL(Z)
                 Y=DIMAG(Z)
                 IF (X.EQ.0.0.AND.Y.GE.0.0) THEN
                    PHI=0.5D0*PI
                 ELSE IF (X.EQ.0.0.AND.Y.LE.0.0) THEN
                    PHI=-0.5D0*PI
                 ELSE
                    PHI=DATAN(Y/X)
                 ENDIF
                 IF (PHI.GT.-0.5*PI.AND.PHI.LT.1.5*PI) NS=1
                 IF (PHI.GT.-1.5*PI.AND.PHI.LE.-0.5*PI) NS=-1
                 CFAC=CDEXP(NS*CI*PI*A)
                 IF (Y.EQ.0.0D0) CFAC=DCOS(PI*A)
                 CHG1=G2/G3*Z**(-A)*CFAC*CS1
                 CHG2=G2/G1*CDEXP(Z)*Z**(A-B)*CS2
                 CHG=CHG1+CHG2
              ENDIF
25            IF (N.EQ.0) CY0=CHG
              IF (N.EQ.1) CY1=CHG
30         CONTINUE
           IF (A0.GE.2.0D0) THEN
              DO 35 I=1,LA-1
                 CHG=((2.0D0*A-B+Z)*CY1+(B-A)*CY0)/A
                 CY0=CY1
                 CY1=CHG
35               A=A+1.0D0
           ENDIF
           IF (X0.LT.0.0D0) CHG=CHG*CDEXP(-Z)
        ENDIF
        A=A1
        Z=Z0
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)
!
!       ==================================================
!       Purpose: Compute gamma function ‚(x)
!       Input :  x  --- Argument of ‚(x)
!                       ( x is not equal to 0,-1,-2,˙˙˙)
!       Output:  GA --- ‚(x)
!       ==================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,  &
               -0.6558780715202538D0, -0.420026350340952D-1,  &
               0.1665386113822915D0,-.421977345555443D-1,  &
               -.96219715278770D-2, .72189432466630D-2,  &
               -.11651675918591D-2, -.2152416741149D-3,  &
               .1280502823882D-3, -.201348547807D-4,  &
               -.12504934821D-5, .11330272320D-5,  &
               -.2056338417D-6, .61160950D-8,  &
               .50020075D-8, -.11812746D-8,  &
               .1043427D-9, .77823D-11,  &
               -.36968D-11, .51D-12,  &
               -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      double complex function i1cerf(n,z)
!
! compute the n'th derivative of Erf(z)
!
      implicit none
!      double complex  i1cerf
      double complex z,cv,cold,z2,ci
      double precision pi,b,eps,epsr
      integer n,nn,i
      parameter ( eps = 1d-12  )
      parameter ( ci = dcmplx(0d0,1d0) )
      parameter ( pi =  3.1415926535897932d0  )
      if (n.eq.0) then
! Erf
       z2=z*z
       if (abs(z).lt.5) then
        b=0.5d0
        cold=0
        epsr=eps
        i1cerf=1
        cv=i1cerf
        i=0
1        i=i+1
         cv=cv*z2/(b+i)
         i1cerf=i1cerf+cv
        if (abs(real(cv))+abs(imag(cv)).gt.epsr*  &
          (abs(real(i1cerf))+abs(imag(i1cerf)))) goto 1
        if (cold.eq.0) then
         cold=i1cerf
         i1cerf=0
         epsr=1d-3
         goto 1
        endif
        i1cerf=(cold+i1cerf)*2*z*exp(-z2)/sqrt(pi)
       else
        z2=1/(2*z*z)
        cv=-z2
        cold=cv
        i=1
2        i=i+1
         cv=-cv*(2*i-1)*z2
         cold=cold+cv
         if (2*i-1.gt.abs(1/z2)) goto 3
         if (i.eq.300) then
          write(*,*) 'not converged'
          goto 3
         endif
        if (abs(real(cv))+abs(imag(cv)).gt.eps* &
          (abs(real(cold))+abs(imag(cold)))) goto 2
3       i1cerf=1-(1+cold)*exp(-z*z)/(z*sqrt(pi))
       endif
      else
! d Erf /dz
       nn=n-1
       i1cerf=1
       if (nn.gt.0) then
        cold=i1cerf
        i1cerf=2*z
        if (nn.gt.1) then
         do i=2,nn
          cv=i1cerf
          i1cerf=2*z*i1cerf-2*(i-1)*cold
          cold=cv
         enddo
        endif
       endif
       i1cerf=i1cerf*(-1)**nn*2*exp(-z*z)/sqrt(pi)
      endif
      return
      end function i1cerf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
