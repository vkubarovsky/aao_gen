c-----------------------------------------------------------------------
c     check_lt_convention -- are the two Born treatments consistent?
c
c     dvmpw returns BOTH a cross section (sigma0) and the structure
c     functions S_T, S_L, S_TT, S_LT, S_LTP.  Two places in aao then use
c     them, and they must agree:
c
c       sigma0   as dvmpw builds it internally   (used by aao_norad to
c                sample events, and as the Born of aao_rad)
c       born_AO  as aao_rad.F:977 rebuilds it from the exported
c                structure functions, in the AO convention
c
c     S_LT and S_LTP carry the factor 2 of the AO/MAID basis, so the
c     matching coefficient is sqrt(eps*(1+-eps)/2) -- see xsection.F:22,
c     "vlt and vltp are divided by 2 ... to match AO formalism".  Using
c     the canonical sqrt(2*eps*(1+-eps)) on an AO-basis quantity doubles
c     the cos(phi) interference; that was the state of the code before
c     c83ce8e, and it is what this program detects.
c
c     Build (needs nothing but the two model files):
c       gfortran -fno-automatic -ffixed-line-length-none \
c         tools/check_lt_convention.F aao_rad/dvmpw.F aao_rad/dvmpx.F \
c         -o tools/check_lt_convention
c
c     Run (defaults are the point of the October 2026 aao/exclurad
c     comparison; override as  check_lt_convention E0 Q2 xB |t|):
c       tools/check_lt_convention
c       tools/check_lt_convention 5.75 1.125 0.1372 0.12
c
c     Expected output: ratio = 1.000000 at every phi.  Before the fix it
c     ran 0.801 / 1.000 / 1.987 at phi = 0 / 90 / 180 -- exactly 1 at 90
c     and 270, where cos(phi) = 0 kills the LT term.
c-----------------------------------------------------------------------
      PROGRAM CHECKLT
      IMPLICIT NONE
      REAL COST,W,Q2,PHI,E,MESONMASS,TABS,XB
      REAL SIGMA0,S_T,S_L,S_TT,S_LT,S_LTP
      REAL AMP2,AMU2,EPICM,PPICM,QVCM,NUCM,W2,EPS,BORNAO,PI
      REAL RATIO,WORST,EPSILON
      EXTERNAL EPSILON
      INTEGER HELI,I,IARGC
      CHARACTER*32 ARG
      PARAMETER (PI=3.1415926536)

      E    = 10.6
      Q2   = 1.5
      XB   = 0.30
      TABS = 0.30
      IF(IARGC().GE.1) THEN
        CALL GETARG(1,ARG)
        READ(ARG,*) E
      ENDIF
      IF(IARGC().GE.2) THEN
        CALL GETARG(2,ARG)
        READ(ARG,*) Q2
      ENDIF
      IF(IARGC().GE.3) THEN
        CALL GETARG(3,ARG)
        READ(ARG,*) XB
      ENDIF
      IF(IARGC().GE.4) THEN
        CALL GETARG(4,ARG)
        READ(ARG,*) TABS
      ENDIF

      MESONMASS = 0.134976
      HELI = 0
      AMP2 = 0.93827**2
      AMU2 = MESONMASS**2
      W2   = AMP2 + Q2*(1.0/XB - 1.0)
      W    = SQRT(W2)
      EPICM = 0.5*(W2 + AMU2 - AMP2)/W
      PPICM = SQRT(EPICM**2 - AMU2)
      QVCM  = SQRT(((W2 + Q2 + AMP2)/2.0/W)**2 - AMP2)
      NUCM  = (W2 - AMP2 - Q2)/(2*W)
      COST  = (Q2 - AMU2 + 2*NUCM*EPICM - TABS)/(2*QVCM*PPICM)
      IF(ABS(COST).GT.1.0) THEN
        WRITE(*,'(A,F8.4)') 'unphysical point: cos(theta*) = ',COST
        STOP 1
      ENDIF
      EPS = EPSILON(XB,Q2,E)

      WRITE(*,'(A,F6.2,A,F6.3,A,F6.3,A,F6.3)')
     *  'E0 = ',E,'   Q2 = ',Q2,'   xB = ',XB,'   |t| = ',TABS
      WRITE(*,'(A,F8.4,A,F8.4,A,F8.4)')
     *  'W  = ',W,'   cos(theta*) = ',COST,'   eps = ',EPS
      WRITE(*,'(A)') ' '
      WRITE(*,'(A)')
     *  '  phi [deg]     sigma0_dvmpw         born_AO      ratio'

      WORST = 0.0
      DO I = 0,24
        PHI = I*15.0
        CALL DVMPW(COST,W,Q2,PHI*PI/180.,E,HELI,MESONMASS,
     *             SIGMA0,S_T,S_L,S_TT,S_LT,S_LTP)
        BORNAO = S_T + EPS*S_L
     *         + EPS*S_TT*COS(2*PHI*PI/180.)
     *         + SQRT(EPS*(1.+EPS)/2.)*S_LT*COS(PHI*PI/180.)
        RATIO = BORNAO/SIGMA0
        WORST = MAX(WORST,ABS(RATIO-1.0))
        WRITE(*,'(F10.1,2E18.7,F11.6)') PHI,SIGMA0,BORNAO,RATIO
      ENDDO

      WRITE(*,'(A)') ' '
      WRITE(*,'(A,E10.3)') 'worst |ratio - 1| = ',WORST
      IF(WORST.GT.1.0E-4) THEN
        WRITE(*,'(A)') 'FAIL: the two Born treatments disagree.'
        STOP 1
      ENDIF
      WRITE(*,'(A)') 'PASS: the two Born treatments agree.'
      END
