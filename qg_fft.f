C******************************************************************************|
C qg_fft.f, the FFT package for diablo.                               VERSION 0.3e
C
C This file isolates all calls to the FFTW package (available at: www.fftw.org)
C These wrapper routines were written by T. Bewley (spring 2001).
c
c I have slightly modified this file so that at the time being it only performs
c 2D ffts.  Hopefully I will not do too much damage! AFT (December '06)
c
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C The arrangement of the significant real numbers in the arrays (denoted by +)
C in physical space, in Fourier space, and in Fourier space after packing are
C shown below for the 2D (X-Z) plane.  The third direction (Y) is handled in
C an identical matter as the Z direction shown here.
C
C      oooooooooooooooooo         oooooooooooooooooo         oooooooooooooooooo
C      oooooooooooooooooo         oooooooooooooooooo         oooooooooooooooooo
C NZ-1 ++++++++++++++++oo     -1  ++++++++++++oooooo         oooooooooooooooooo
C      ++++++++++++++++oo     -2  ++++++++++++oooooo         oooooooooooooooooo
C      ++++++++++++++++oo     -3  ++++++++++++oooooo         oooooooooooooooooo
C      ++++++++++++++++oo         ++++++++++++oooooo         oooooooooooooooooo
C      ++++++++++++++++oo    -NKZ ++++++++++++oooooo         oooooooooooooooooo
C      ++++++++++++++++oo         oooooooooooooooooo     -1  ++++++++++++oooooo
C      ++++++++++++++++oo         oooooooooooooooooo     -2  ++++++++++++oooooo
C      ++++++++++++++++oo         oooooooooooooooooo     -3  ++++++++++++oooooo
C      ++++++++++++++++oo         oooooooooooooooooo         ++++++++++++oooooo
C      ++++++++++++++++oo         oooooooooooooooooo    -NKZ ++++++++++++oooooo
C      ++++++++++++++++oo     NKZ ++++++++++++oooooo     NKZ ++++++++++++oooooo
C      ++++++++++++++++oo         ++++++++++++oooooo         ++++++++++++oooooo
C   3  ++++++++++++++++oo      3  ++++++++++++oooooo      3  ++++++++++++oooooo
C   2  ++++++++++++++++oo      2  ++++++++++++oooooo      2  ++++++++++++oooooo
C   1  ++++++++++++++++oo      1  ++++++++++++oooooo      1  ++++++++++++oooooo
C   0  ++++++++++++++++oo      0  +o++++++++++oooooo      0  +o++++++++++oooooo
C      ^^^^           ^           ^ ^ ^     ^                ^ ^ ^     ^
C      0123           NX-1        0 1 2     NKX              0 1 2     NKX
C
C       PHYSICAL SPACE              FOURIER SPACE         FOURIER SPACE (PACKED)
C
C After the Real->Fourier transform, the significant coefficients are put next
C to each other in the array, so a loop such as
C
C        DO K=0,TNKZ           [where TNKZ = 2*NKZ = 2*(NZ/3) ]
C          DO I=0,NKX          [where  NKX = NX/3             ]
C            CP(I,K,J)= ...
C          END DO
C        END DO
C
C includes all the Fourier coefficients of interest.  The subsequent loops in
C Fourier space just work on these coefficients in the matrix.
C  
C Before a Fourier->Real transform, the significant coefficients are unpacked
C and the higher wavenumbers are SET TO ZERO before the inverse transform.
C This has the effect of doing the required dealiasing.
C
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_FFT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'qgbt_diablo.h'
      INTEGER I,J,K

      INTEGER         FFTW_FORWARD,      FFTW_BACKWARD,
     *                FFTW_ESTIMATE,     FFTW_MEASURE,
     *                FFTW_OUT_OF_PLACE, FFTW_IN_PLACE,
     *                FFTW_USE_WISDOM,   FFTW_THREADSAFE
      PARAMETER(      FFTW_FORWARD=-1,      FFTW_BACKWARD=1,
     *                FFTW_ESTIMATE=0,      FFTW_MEASURE=1,
     *                FFTW_OUT_OF_PLACE=0,  FFTW_IN_PLACE=8,
     *                FFTW_USE_WISDOM=16,   FFTW_THREADSAFE=128 )

      WRITE(6,*) 'Initializing FFTW package.'
      PI = 4. * ATAN(1.0)
      CI = CMPLX(0.0,1.0)
      EPS= 0.000000001

      IF (NUM_PER_DIR .GT. 0) THEN
        CALL RFFTWND_F77_CREATE_PLAN(FFTW_X_TO_F_PLAN, 1, NX,
     *        FFTW_FORWARD,  FFTW_MEASURE + FFTW_IN_PLACE )
        CALL RFFTWND_F77_CREATE_PLAN(FFTW_X_TO_P_PLAN, 1, NX,
     *        FFTW_BACKWARD, FFTW_MEASURE + FFTW_IN_PLACE )
        NKX=NX/2
        RNX=1.0*NX
        DO I=0,NKX
          KX(I)=I*(2.*PI)/LX
          KX2(I)=KX(I)*KX(I)
          CIKX(I)=CI*KX(I)
        END DO
      END IF

      IF (NUM_PER_DIR .GT. 1) THEN
        CALL FFTWND_F77_CREATE_PLAN(FFTW_Y_TO_F_PLAN, 1, NY,
     *       FFTW_FORWARD,  FFTW_MEASURE + FFTW_IN_PLACE )
        CALL FFTWND_F77_CREATE_PLAN(FFTW_Y_TO_P_PLAN, 1, NY,
     *       FFTW_BACKWARD, FFTW_MEASURE + FFTW_IN_PLACE )
        NKY=NY/2
        TNKY=2*NKY
        RNY=1.0*NY
        DO J=0,NKY
          KY(J)=J*(2.*PI)/LY
        END DO
        DO J=1,NKY
          KY(TNKY+1-J)=-J*(2.*PI)/LY
        END DO
        DO J=0,TNKY
          KY2(J)=KY(J)*KY(J)
          CIKY(J)=CI*KY(J)
        END DO
      END IF

      IF (NUM_PER_DIR .GT. 2) THEN
        CALL FFTWND_F77_CREATE_PLAN(FFTW_Y_TO_F_PLAN, 1, NY,
     *       FFTW_FORWARD,  FFTW_MEASURE + FFTW_IN_PLACE )
        CALL FFTWND_F77_CREATE_PLAN(FFTW_Y_TO_P_PLAN, 1, NY,
     *       FFTW_BACKWARD, FFTW_MEASURE + FFTW_IN_PLACE )
        NKY=NY/2
        TNKY=2*NKY
        RNY=1.0*NY
        DO J=0,NKY
          KY(J)=J*(2.*PI)/LY
        END DO
        DO J=1,NKY
          KY(TNKY+1-J)=-J*(2.*PI)/LY
        END DO
        DO J=0,TNKY
          KY2(J)=KY(J)*KY(J)
          CIKY(J)=CI*KY(K)
        END DO
      END IF

      DO J=0,TNKY
         DO I=0,NKX
            KSQD(I,J)=KX2(I)+KY2(J)
         END DO
      END DO
	write(6,*) ksqd(0,1) , ksqd(1,0)
      KSQD(0,0)=0.0
! 	stop

      DO I=0,NKX
        DO J=0,NY
          CYX_PLANE(J,I)=CMPLX(0.0,0.0)
        END DO
      END DO
c      DO K=0,TNKZ
c        DO J=0,NY
c          CYZ_PLANE(J,K)=CMPLX(0.0,0.0)
c        END DO
c      END DO

      WRITE(6,*) 'FFTW package initialized.'

      RETURN
      END

C******************************************************************************|
C-------------> The transform routines for the duct flow follow. <-------------|
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_X_TO_FOURIER(U,CU,JMIN,JMAX,idx)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'qgbt_diablo.h'
      REAL*8     U (0:NX+1,0:NY+1,0:2)
      COMPLEX*16 CU(0:NX/2,0:NY+1,0:2)
      INTEGER JMIN, JMAX, I, J, idx

C Looping over the planes of interest, simply perform a real -> complex
C transform in place in the big storage array, scaling appropriately.

c      DO J=JMIN,JMAX
       CALL RFFTWND_F77_REAL_TO_COMPLEX(FFTW_X_TO_F_PLAN,(JMAX-JMIN+1),
     *    U(0,JMIN,idx), 1, NX+2, CU(0,JMIN,idx), 1, NX/2+1)
        DO J=JMIN,JMAX
          DO I=0,NKX
            CU(I,J,idx)=CU(I,J,idx)/RNX
          END DO
        END DO
c      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_X_TO_PHYSICAL(CPP,P,JMIN,JMAX,idx)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) planes to physical space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'qgbt_diablo.h'
      REAL*8     P (0:NX+1,0:NY+1,0:2)
      COMPLEX*16 CPP(0:NX/2,0:NY+1,0:2)
      INTEGER JMIN, JMAX, I, J, idx

C Looping over the planes of interest, simply set the higher wavenumbers to
C zero and then perform a complex -> real transform in place in the big
C storage array.

c      DO K=KMIN,KMAX
        DO J=JMIN,JMAX
          DO I=NKX+1,NX/2
            CPP(I,J,idx)=0.
          END DO
        END DO
       CALL RFFTWND_F77_COMPLEX_TO_REAL(FFTW_X_TO_P_PLAN,(JMAX-JMIN+1),
     *    CPP(0,JMIN,idx), 1, NX/2+1, P(0,JMIN,idx), 1, NX+2)
c      END DO

      RETURN
      END

C******************************************************************************|
C-----------> The transform routines for the channel flow follow. <------------|
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XY_TO_FOURIER(U,CU,idx)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 2 directions)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'qgbt_diablo.h'
      REAL*8     U (0:NX+1,0:NY+1,0:2)
      COMPLEX*16 CU(0:NX/2,0:NY+1,0:2)
      INTEGER I, J, idx

C Transform in the X direction using FFT_X_TO_FOURIER, put the data into the 
C CZX_PLANE temporary storage variable, perform a complex -> complex transform 
C in the z direction, then put the data back into the big storage array, packing 
C the data towards K=0 and scaling appropriately.

c      DO J=JMIN,JMAX
        CALL FFT_X_TO_FOURIER(U,CU,0,NYM,idx)
        DO I=0,NKX
          DO J=0,NYM
            CYX_PLANE(J,I)=CU(I,J,idx)
          END DO
        END DO        
        CALL FFTWND_F77(FFTW_Y_TO_F_PLAN, NKX+1,
     *    CYX_PLANE(0,0), 1, NY+1, CYX_PLANE(0,0), 1, NY+1)
        DO I=0,NKX
          DO J=0,NKY
            CU(I,J,idx)=CYX_PLANE(J,I)/RNY
          END DO
          DO J=1,NKY
            CU(I,NKY+J,idx)=CYX_PLANE(NYM-NKY+J,I)/RNY
          END DO
        END DO
c      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XY_TO_PHYSICAL(CU,U,idx)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 2 directions) planes.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'qgbt_diablo.h'
      REAL*8     U (0:NX+1,0:NY+1,0:2)
      COMPLEX*16 CU(0:NX/2,0:NY+1,0:2)
      INTEGER I, J, idx

C Unpack data into the CYX_PLANE temporary storage variable (setting higher 
C wavenumbers to zero), perform a complex -> complex transform in the z 
C direction, then put the data back into the big storage array and transform in 
C the X direction using FFT_X_TO_PHYSICAL.

c      DO J=JMIN,JMAX
        DO I=0,NKX
          DO J=0,NKY
            CYX_PLANE(J,I)=CU(I,J,idx)
          END DO
          DO J=NKY+1,NYM-NKY
            CYX_PLANE(J,I)=CMPLX(0.0,0.0)
          END DO
          DO J=1,NKY
            CYX_PLANE(NYM-NKY+J,I)=CU(I,NKY+J,idx)
          END DO
        END DO
        CALL FFTWND_F77(FFTW_Y_TO_P_PLAN, NKX+1,
     *    CYX_PLANE(0,0), 1, NY+1, CYX_PLANE(0,0), 1, NY+1)
        DO I=0,NKX
          DO J=0,NYM
            CU(I,J,idx)=CYX_PLANE(J,I)
          END DO
        END DO        
        CALL FFT_X_TO_PHYSICAL(CU,U,0,NYM,idx)
c      END DO

      RETURN
      END



C******************************************************************************|
C---------> A test subroutine for debugging the FFT wrappers follow. <---------|
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE TEST_FFT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'qgbt_diablo.h'
      INTEGER I, J

      WRITE(6,*) 
      WRITE(6,*) 'Now testing 1D FFT routine.  Fourier-space data:'
      CUG(3,0,0)=CMPLX(1.0,0.0)
      DO I=0,NKX
        WRITE(6,'(2F8.4)') CUG(I,0,0)
      END DO
      CALL FFT_X_TO_PHYSICAL(CUG,UG,0,0,0)
      WRITE(6,*) 'After Fourier->Physical transform:'
      DO I=0,NXM
        WRITE(6,'(F8.4)') UG(I,0,0)
      END DO
      CALL FFT_X_TO_FOURIER(UG,CUG,0,0,0)
      WRITE(6,*) 'After Physical->Fourier transform:'
      DO I=0,NKX
        WRITE(6,'(2F8.4)') CUG(I,0,0)
      END DO

      WRITE(6,*)
      WRITE(6,*) 'NOW TESTING 2D FFT ROUTINE.  Fourier-space data:'
      CVG(0,0,0)=CMPLX(1.0,0.0)
      DO J=0,TNKY
        WRITE(6,'(10F8.4)') (CVG(I,J,0),I=0,NKX)
      END DO
      CALL FFT_XY_TO_PHYSICAL(CVG,VG,0)
      WRITE(6,*) 'After Fourier->Physical transform:'
      DO J=0,NYM
        WRITE(6,'(10F8.4)') (VG(I,J,0),I=0,NXM)
      END DO
      CALL FFT_XY_TO_FOURIER(VG,CVG,0)
      WRITE(6,*) 'After Physical->Fourier transform:'
      DO J=0,TNKY
        WRITE(6,'(10F8.4)') (CVG(I,J,0),I=0,NKX)
      END DO

      WRITE(6,*)
      STOP 'FFT Test Complete.'

      RETURN
      END






