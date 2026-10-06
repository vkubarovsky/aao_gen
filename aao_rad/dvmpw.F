c
c*******************Valery Model****************************************
c
      subroutine DVMPW(COST,W,Q2, Phi_g,E,heli,MESONMASS,
     *                   sigma0,S_T,S_L,S_TT,S_LT,S_LTP)
C
c  dsigma/dQ2 dphi dCosTheta* dW for ep-->ep pi0
C
C  exc pi0 x-section
c
c input:
c COST = Cos(Theta*), Theta* is the angle (pi0-gamma* or  p-p')in the CM system
c W is W
c xb,Q2 x and Q^2 (GeV^2)
c Phi_g angle in the photon frame (radians)
c E energy of the electron in GeV
c heli electron helicity -1.0 or +1.0
c
c del2=t (negative GeV^2)           NEGATIVE !!!!
C
C MESONMASS is the mass of the pi0 or eta.
c

      IMPLICIT NONE
      REAL    COST,W,del2,xb,Q2, Phi_g,E,mesonmass
      INTEGER heli
      REAL    Mp, mele, pi
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)

      REAL     sigma0, S_T, S_L, S_LT, S_TT,S_LTP
      REAL     SIGMA_T, SIGMA_L, SIGMA_LT, SIGMA_TT,SIGMA_LTP
      REAL     EPSILON, FLUXW, SIGMA_TOT
      EXTERNAL EPSILON,SIGMA_T,SIGMA_L, SIGMA_LT, SIGMA_TT
      EXTERNAL SIGMA_LTP
      REAL     JACG,JACR
      LOGICAL  CHECK_KINE
      REAL S_TOT
      DOUBLE PRECISION eps
C      INCLUDE 'pi0eta.par'
c
      CALL XSINIT(MESONMASS)
c     
      sigma0 = 0.0
      S_T = 0.0
      S_L = 0.0
      S_LT = 0.0
      S_TT = 0.0
      S_LTP = 0.0
      IF(.NOT.CHECK_KINE(COST,W,Q2,E,del2,xb,JACG,JACR)) RETURN 
      EPS=EPSILON(XB,q2,e)

      

      S_T  =  0.001/(2*pi)*SIGMA_T     (COST,W,Q2,E)
      S_L  =  0.001/(2*pi)*SIGMA_L     (COST,W,Q2,E)
      S_LT =  0.001/(2*pi)*2*SIGMA_LT  (COST,W,Q2,E)
      S_TT =  0.001/(2*pi)*SIGMA_TT    (COST,W,Q2,E)
      S_LTP = 0.001/(2*pi)*2*SIGMA_LTP (COST,W,Q2,E)

c
c     AO convention.  S_LT and S_LTP above already carry the factor 2 that
c     the two OFF-DIAGONAL structure functions have in the AO/MAID basis
c     (cf. xsection.F: vlt = sqrt(eps*(1+eps)/2), "divided by 2 ... to
c     match AO formalism").  The matching kinematic coefficient is therefore
c     sqrt(eps*(1+-eps)/2) and NOT the canonical sqrt(2*eps*(1+-eps)):
c     using the canonical one here applies a canonical prefactor to an
c     AO-basis quantity and doubles the cos(phi) interference term.
c     Measured on the RG-A aao_norad production (2026-09-22): cos(phi)
c     harmonic 0.222 vs the canonical 0.097 at Q2 = 1.0-1.6 GeV2, while
c     cos(2phi) agrees (-0.472 vs -0.467); sigma_TT is diagonal and is
c     not doubled.  EXCLURAD never hit this because it ignores sigma0 and
c     builds its own Born from the structure functions with the AO
c     coefficient (exclurad.F:407).  The radiative branch of aao_rad is
c     likewise already correct (aao_rad.F:977).
      sigma0 = S_T   + EPS*S_L +         
     *   EPS                  *S_TT  * COS(2*PHI_G)   + 
     * SQRT(EPS*(1.+EPS)/2.)*S_LT  * COS(PHI_G)     +
     * HELI*SQRT(EPS*(1.-EPS)/2.)*S_LTP * SIN(PHI_G)     
       IF(sigma0 .LT. 0.0) sigma0=0.0     
c
c       S_TOT =FLUXW(xb,Q2,E)/(2.*PI)*(
c     *                       S_T   + EPS*S_L          + 
c     * EPS                  *S_TT  * COS(2*PHI_G)       + 
c     * SQRT(2.*EPS*(1.+EPS))*S_LT  * COS(PHI_G)     +
c     * HELI*SQRT(2.*EPS*(1.-EPS))*S_LTP * SIN(PHI_G)     
c     * )
      
c       IF(DVMPW.LT.0.0) DVMPW=0.
c       PRINT *,'DVMPW: S_tot ',S_TOT,COST,W,Q2,-del2,xb,phi_g,mesonmass,
c     *  eps,FLUXW(xb,Q2,E)

c       print *,'S_T ',  S_T,S_L,S_T   + EPS*S_L 
c       print *,'S_TT  ',  S_TT, EPS     *S_TT  * COS(2*PHI_G) 
c       print *, 'S_LT ',S_LT,SQRT(2.*EPS*(1.+EPS))*S_LT  * COS(PHI_G) 

      RETURN
      END

C====================================================================================c

      REAL FUNCTION SIGMA_T(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,del2,x,Q2,E
      REAL          JACG,JACR
      LOGICAL       CHECK_KINE
      EXTERNAL      XSIGMA_T
      REAL          XSIGMA_T
      
      IF(CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) THEN
        SIGMA_T=JACR*XSIGMA_T(del2,x,Q2,E)
      ELSE
        SIGMA_T=0.0
      ENDIF
      RETURN
      END
C=============================================================      
      REAL FUNCTION SIGMA_L(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,del2,x,Q2,E
      REAL          JACG,JACR
      LOGICAL       CHECK_KINE
      EXTERNAL      XSIGMA_L
      REAL          XSIGMA_L
      
      IF(CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) THEN
        SIGMA_L=JACR*XSIGMA_L(del2,x,Q2,E)
      ELSE
        SIGMA_L=0.0
      ENDIF
      RETURN
      END
      
C====================================================================================c

      REAL FUNCTION SIGMA_U(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,Q2,E
      REAL          del2,x,JACG,JACR
      EXTERNAL      SIGMA_T,SIGMA_L,EPSILON
      REAL          SIGMA_T,SIGMA_L,EPSILON
      LOGICAL       CHECK_KINE

      SIGMA_U=0.
      IF(.NOT.CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) RETURN
      SIGMA_U=          SIGMA_T(COST,W,Q2,E)+
     * EPSILON(X,q2,e)*SIGMA_L(COST,W,Q2,E)

      RETURN
      END
C=============================================================      
      REAL FUNCTION SIGMA_TT(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,del2,x,Q2,E
      REAL          JACG,JACR
      LOGICAL       CHECK_KINE
      EXTERNAL      XSIGMA_TT
      REAL          XSIGMA_TT
      
      IF(CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) THEN
        SIGMA_TT=JACR*XSIGMA_TT(del2,x,Q2,E)
      ELSE
        SIGMA_TT=0.0
      ENDIF
      RETURN
      END
      
C=============================================================      
      REAL FUNCTION SIGMA_LT(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,del2,x,Q2,E
      REAL          JACG,JACR
      LOGICAL       CHECK_KINE
      EXTERNAL      XSIGMA_LT
      REAL          XSIGMA_LT
      
      IF(CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) THEN
        SIGMA_LT=JACR*XSIGMA_LT(del2,x,Q2,E)
      ELSE
        SIGMA_LT=0.0
      ENDIF
      RETURN
      END
      
C====================================================================================c

       REAL FUNCTION SIGMA_LTP(COST,W,Q2,E)
       IMPLICIT NONE
       REAL                    COST,W,Q2,E

       SIGMA_LTP=0.

       RETURN
       END

C====================================================================================c

      LOGICAL FUNCTION CHECK_KINE(COST,W,Q2,E,del2,xb,JACG,JACR)
      IMPLICIT NONE
      REAL COST,W,del2,xb,Q2,E,JACG,JACR
      real Mp, mele, pi,mpi0
      real fluxw, rxs
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)
            

      real nu,W2,qmod,E1cm,P1cm,E2cm,P2cm,del2max,del2min
      real  xmin1,xmax1
      double precision y, e1, epsilon

      REAL AM(2)
      DATA AM/0.134976,0.547300/
      INTEGER K
      COMMON/PIETA/K

      MPI0=AM(K)

      CHECK_KINE=.FALSE.
      W2=W*W
      XB=Q2/(W2+Q2-Mp*Mp)
      nu  = (W2+Q2-Mp*Mp)/(2*Mp)
      qmod = sqrt(nu**2 + Q2)
      IF(W         .LE. Mp+Mpi0)              RETURN
      IF(W2        .GE. -Q2+2*Mp*E+Mp*Mp)     RETURN
      IF(Q2/(4*E*(E-NU)) .GE. 1.0)            RETURN
      IF(XB        .LE. Q2/(2*Mp*E))          RETURN
      IF(ABS(COST) .GE. 1.0)                  RETURN
      xmin1 = Q2/(2.0*Mp*E)
      xmax1 = 1.0

      E1cm = Mp*(Mp + nu)/W
      E2cm = (W2 + Mp**2-Mpi0**2)/(2D0*W)
      IF(E1cm.LE.Mp .OR. E2cm.LE. Mp) RETURN
      P1cm = Mp*qmod/W
      P2cm = SQRT(E2CM**2 - Mp**2)
      del2max = 2.0*(Mp**2 - E1cm*E2cm - P1cm*P2cm)
      del2min = 2.0*(Mp**2 - E1cm*E2cm + P1cm*P2cm)
      del2=del2min-2.*P1cm*P2cm*(1-COST)
      IF( xb.le.xmin1 .or. xb.gt.xmax1 )         return   !    x  out of range
      IF( del2.ge.del2min .or. del2.le.del2max ) return   ! delta out of range

      y=Q2/(2*Mp*xb*E)
      e1=(y*xb*Mp)**2/Q2
      EPSILON=(1.0-y-e1)/(1-y+y**2/2+e1)

      IF(EPSILON.LE.0. .OR.EPSILON .GE. 1.0)         RETURN

C      jacobian
C      JAC=4*P1cm*P2CM*W*XB/(W2+Q2-Mp*Mp)
      JACG=2*W*XB/(W2+Q2-Mp*Mp)
      JACR=2*P1cm*P2CM

      CHECK_KINE=.TRUE.
c      PRINT *,'Check_kine: ', W,xB,nu,qmod,E1cm,E2cm,del2,y,e1
      RETURN
      END


c*****************The end of Valery Model******************************
c
