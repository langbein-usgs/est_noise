      SUBROUTINE SRAND(ISEED)
c from https://www.cisl.ucar.edu/zine/96/spring/articles/3.random-6.html
C
C  This subroutine sets the integer seed to be used with the
C  companion RAND function to the value of ISEED.  A flag is 
C  set to indicate that the sequence of pseudo-random numbers 
C  for the specified seed should start from the beginning.
C
      COMMON /SEED/JSEED,IFRST
C
      JSEED = ISEED
      IFRST = 0
C
      RETURN
      END
      REAL FUNCTION RAND()
C
C  This function returns a pseudo-random number for each invocation.
C  It is a FORTRAN 77 adaptation of the "Integer Version 2" minimal 
C  standard number generator whose Pascal code appears in the article:
C
C     Park, Steven K. and Miller, Keith W., "Random Number Generators: 
C     Good Ones are Hard to Find", Communications of the ACM, 
C     October, 1988.
C
      PARAMETER (MPLIER=16807,MODLUS=2147483647,MOBYMP=127773,
     +           MOMDMP=2836)
C
      COMMON  /SEED/JSEED,IFRST
      INTEGER HVLUE, LVLUE, TESTV, NEXTN
      SAVE    NEXTN
C
      IF (IFRST .EQ. 0) THEN
        NEXTN = JSEED
        IFRST = 1
      ENDIF
C
      HVLUE = NEXTN / MOBYMP
      LVLUE = MOD(NEXTN, MOBYMP)
      TESTV = MPLIER*LVLUE - MOMDMP*HVLUE
      IF (TESTV .GT. 0) THEN
        NEXTN = TESTV
      ELSE
        NEXTN = TESTV + MODLUS
      ENDIF
      RAND = REAL(NEXTN)/REAL(MODLUS)
C
      RETURN
      END
      BLOCKDATA RANDBD
      COMMON /SEED/JSEED,IFRST
C
      DATA JSEED,IFRST/123456789,0/
C
      END
      FUNCTION ZBQLNOR(MU,SIGMA)
c From http://www.ucl.ac.uk/~ucakarc/work/software/randgen.f
*
*       Returns a random number Normally distributed with mean
*       MU and standard deviation |SIGMA|, using the Box-Muller
*       algorithm
*
      DOUBLE PRECISION THETA,R,ZBQLNOR,ZBQLU01,PI,MU,SIGMA
      DOUBLE PRECISION SPARE
      INTEGER STATUS
      SAVE STATUS,SPARE,PI
      DATA STATUS /-1/
      IF (STATUS.EQ.-1) PI = 4.0D0*DATAN(1.0D0)

      IF (STATUS.LE.0) THEN
c       THETA = 2.0D0*PI*ZBQLU01(0.0D0)
       THETA = 2.0D0*PI*dble(rand())
c       R = DSQRT( -2.0D0*DLOG(ZBQLU01(0.0D0)) )
       R=DSQRT (-2.0d+0*dlog(dble(rand()) ) )
       ZBQLNOR = (R*DCOS(THETA))
       SPARE = (R*DSIN(THETA))
       STATUS = 1
      ELSE
       ZBQLNOR = SPARE
       STATUS = 0
      ENDIF
      
      ZBQLNOR = MU + (SIGMA*ZBQLNOR)

      END
