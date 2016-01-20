#include "definesld.com"
MODULE INI_3D_MODULE

USE KINDS
USE PARMSLD
USE CONSTLD
USE WORKCOM
USE CONST3D
USE UTILS
USE BOUND
USE UPDATE_TRACER_MODULE
USE Z_COORD
USE DOMAIN_DECOMPOSITION
#if defined (DIFFUSION)
USE TURB_3D_MODULE
#if defined (LSM)
USE LAND_MODULE
#endif
#endif
USE PETSCSYS

IMPLICIT NONE
PRIVATE

! PUBLIC MEMBER FUNCTIONS
PUBLIC ::     &
   INI_3D,     &
   PTSPEC  
   
CONTAINS

!=======================================================================
   SUBROUTINE INI_3D (ITT)
!=======================================================================

!     SET UP I.C. 

      INTEGER (KIND=INT_KIND), INTENT(IN) :: ITT

      INTEGER (KIND=INT_KIND) :: &
         RDHOUR  
      INTEGER (KIND=INT_KIND) :: &
         I, J, K,  & ! DO LOOP INDICES FOR ZONAL, MERIDIONAL AND VERTICAL DIMENSIONS
         N,        & ! TIME LEVEL INDEX
         I00, J00, & !
         IQ,JQ 
      REAL (KIND=DBL_KIND) ::    &
         SCALEQ=0._DBL_KIND, & 
         RKAPPA,             &
         RKAPPAI,            &
         GAMFAC, A0, Z00, XR, YR, ZR, DIST,TMP
      REAL (KIND=DBL_KIND), PARAMETER :: TH00 = 300._DBL_KIND      
      REAL (KIND=DBL_KIND), PARAMETER :: RGAS = 287.04     
      REAL (KIND=DBL_KIND), PARAMETER :: RAD1 = -2.    
      REAL (KIND=DBL_KIND), PARAMETER :: SECDAY = 86400.    
      REAL (KIND=DBL_KIND) :: PZERO = 1000._DBL_KIND      
      REAL (KIND=DBL_KIND) :: P_INT(NK2),PI_INT(NK2) 
      REAL (KIND=INT_KIND) :: TEMP(NK3)

      INTEGER NU1, NT1, K1, K2, IT
      LOGICAL LF
      DATA NU1/33/,NT1/33/
      REAL(KIND=DBL_KIND) :: RBCP, T8, P8, QV8, QC8, QI8, DUMMY, GDZBCP, CPBR, &
                             ES1, ALPHAW, PBARMB
      REAL (KIND=DBL_KIND) :: T(NK3),H(NK3),TV(NK3),PILOG(NK3),ALPHA(NK3)  &
                             ,QCBAR(NK3),QVSBAR(NK3)
      REAL (KIND=DBL_KIND) :: ZT1(33),PI1(33),TV1(33),PILOG1(33)
      REAL (KIND=DBL_KIND) :: P1(33),T1(33),QV1(33),QC1(33)
      REAL (KIND=DBL_KIND) :: U1(33),V1(33)
      REAL (KIND=DBL_KIND) :: Q1OBS(20),Q2OBS(20),PQ1(20)

      REAL (KIND=DBL_KIND) :: SST2
      REAL (KIND=REAL_KIND) :: SST3(MI_GLOB)
      REAL (KIND=REAL_KIND):: TGOUT(MI_GLOB,MJ_GLOB)
      CHARACTER*100 :: PTH
      INTEGER :: ii,num
      

      DATA P1 /1006.7, 993.8, 977.9, 959.1, 937.7, 913.7, 887.4,   &
                858.9, 828.5, 796.4, 762.8, 728.0, 692.2, 655.7,   &
                618.6, 581.3, 543.9, 506.8, 470.2, 434.1, 398.9,   &
                364.6, 331.4, 299.4, 268.7, 239.5, 211.8, 185.8,   &
                161.6, 139.4, 119.2, 101.4,  85.8 /
      DATA T1/300.41,299.95,299.27,298.40,297.60,296.81,295.85,294.69 &
             ,293.38,291.91,290.29,288.72,286.88,284.90,282.62,280.14 &
             ,277.51,274.55,271.49,268.14,264.55,260.73,256.46,251.93 &
             ,246.90,241.30,234.76,227.68,219.85,211.77,204.24,196.16 &
             ,196.03/
      DATA QV1/22.4,22.0,21.6,20.8,20.1,19.3,18.6,18.0,17.1,16.2,15.3 &
              ,14.2,13.3,12.1,11.0,9.8,8.6,7.5,6.5,5.4,4.5,3.6,2.8,2.0 &
              ,1.4,0.8,0.4,0.2,0.1,0.0,0.0,0.0,0.0/
      DATA QC1/33*0./
     

      CASENAME= 'GATE_PHASE_III'

      CZ2 = ( DZ - DZ1 ) / ( DZ * ( DOMAIN - DZ ) )
      CZ1 = 1. - CZ2 * DOMAIN

      CALL COORDS_2D ( CZ1, CZ2, DZ, ZB )

!*****************************************************
!     1. OBTAIN 1D THERMODYNAMICAL PROFILE
!*****************************************************

      SST = 300
      SST2 = 303      !HIGHER SST

      DO K = 1, NU1
      U1(K) = 0
      V1(K) = 0
      ENDDO

      DO 100 K = 1, NT1
      T1(K) = T1(K)
      QV1(K)= QV1(K)*1E-3
  100 CONTINUE

      RBCP = RGAS / CP
      DO 110 K = 1, NT1
      PI1(K) = ( P1(K) / PZERO ) ** RBCP
  110 CONTINUE

!     ADJUST T1, QV1, AND QC1. UNITS : P1, MB.

      DO 120 K = 1, NT1
      T8  = T1(K)
      P8  = P1(K)
      QV8 = QV1(K)
      QC8 = 0.D0
      QI8 = 0.D0
      CALL SATURATION ( 1, 1, K, T8, P8, QV8, QC8, QI8, DUMMY )
      QV1(K) = QV8
  120 CONTINUE

      DO 130 K = 1, NT1
      PILOG1(K) = LOG ( PI1(K) )
      TV1(K) = T1(K) * ( 1. + 0.608 * QV1(K) - QC1(K) )
  130 CONTINUE

        ZT1(1) = ZZ(1)
          T(1) = T1(1)
      QVBAR(1) = QV1(1)
      DO 140 K = 2, NT1
      ZT1(K) = ZT1(K-1)-CP/(2.*GRAV)         &
              *(PILOG1(K)-PILOG1(K-1))*(TV1(K)+TV1(K-1))
  140 CONTINUE

!     CALCULATE CONSTANT :  CALCULATE PIBAR AND PBAR USING THE HYDROSTATIC EQ.

      GDZBCP = 2. * GRAV * DZ / CP
      CPBR = CP / RGAS
      RBCP = 1. / CPBR

      PSFC  = P1(1)
      PSFC  = PSFC * 100.
      PZERO = PZERO * 100.

      PISFC = ( PSFC / PZERO ) ** RBCP

      PBAR(1)  = PSFC
      PIBAR(1) = PISFC
      PILOG(1) = LOG ( PIBAR(1) )

! =====================     INTERPLATE TO MODEL LEVEL    ===================

      DO 150 K = 2, NK3
      K1 = INDEXR ( ZT(K), NT1, ZT1, LF )
      K2 = K1 + 1
      T(K) = FINTRP ( 1, ZT(K), ZT1(K1), T1(K1), ZT1(K2), T1(K2) )
      QVBAR(K) = FINTRP ( 1, ZT(K), ZT1(K1), QV1(K1), ZT1(K2), QV1(K2) )
  150 CONTINUE

      DO 162 K = 1, NK3
      TV(K) = T(K) * ( 1. + 0.608 * QVBAR(K) )
  162 CONTINUE

      PILOG(2)=PILOG(1)-GDZBCP/(TV(1)+TV(2))     &
              *(ZT(2)-ZZ(1))/DZ
      DO 164 K = 3, NK3
      PILOG(K)=PILOG(K-1)-GDZBCP/(TV(K-1)+TV(K))/FNZ(K-1)
  164 CONTINUE

      DO 166 K = 2, NK3
      PIBAR(K)= EXP ( PILOG(K) )
      PBAR(K) = PZERO * PIBAR(K) ** CPBR
  166 CONTINUE

      DO 168 K = 1, NK3
      T8 = T(K)
      P8 = PBAR(K) * 0.01
      QV8 = QVBAR(K)
      QC8 = 0.D0
      QI8 = 0.D0
      CALL SATURATION ( 1, 1, K, T8, P8, QV8, QC8, QI8, DUMMY )
      QVBAR(K) = QV8
      ES1 = ES ( T(K) )
      QVSBAR(K) = ES1 * 0.622 / ( PBAR(K) - ES1 )
      THBAR(K) = T(K) / PIBAR(K)
      H(K) = QVBAR(K) / QVSBAR(K)
  168 CONTINUE

      IF (MY_TASK .EQ. 0 ) PRINT*,(PIBAR(K),PBAR(K),K,K=1,NK3)
!     PROFILES OF RHO, RHOZ, AND GAMMA

      GAMFAC = HLF / CP

      DO 170 K = 2, NK3
      GAMMA(K) = GAMFAC / PIBAR(K)
      ALPHA(K) = RGAS * TV(K) / PBAR(K)
  170 CONTINUE

      DO 180 K = 2, NK2
      ALPHAW = ( ALPHA(K) + ALPHA(K+1) ) / 2.
      RHOZ(K) = 1. / ALPHAW
  180 CONTINUE

!     DENSITY
      RHOZ(1) = PSFC / ( RGAS * TV(1) )
      DO 190 K = 2, NK3
      RHO(K) = 1. / ALPHA(K)
  190 CONTINUE
!     DENSITY FOR K=1/2
      RHO(1) = RHOZ(1)


      DO 101 K=1,NK3
      ZW(K)=ZZ(K)  ! THIS LINE ADDED WITH REMOVAL OF EQUIVALENCE IN CONSTLD
      ZU(K)=ZT(K)  ! THIS LINE ADDED WITH REMOVAL OF EQUIVALENCE IN CONSTLD
      FNU(K)=FNT(K)  ! THIS LINE ADDED WITH REMOVAL OF EQUIVALENCE IN CONSTLD
      FNW(K)=FNZ(K)  ! THIS LINE ADDED WITH REMOVAL OF EQUIVALENCE IN CONSTLD
      RHOT(K)=RHO(K)  ! THIS LINE ADDED WITH REMOVAL OF EQUIVALENCE IN CONSTLD
      RHOU(K)=RHO(K)  ! THIS LINE ADDED WITH REMOVAL OF EQUIVALENCE IN CONSTLD
  101 CONTINUE  
      DO 102 K=1,NK2
      RHOW(K)=RHOZ(K)  ! THIS LINE ADDED WITH REMOVAL OF EQUIVALENCE IN CONSTLD
  102 CONTINUE

      DO 200 K = 2, NK3
      K1 = INDEXR(PBAR(K)*0.01,20,P1,LF)
      K2 = K1 + 1
      Q1LS (K) =  &
        FINTRP(1,PBAR(K)*0.01,P1(K1),Q1OBS(K1),P1(K2),Q1OBS(K2))
      Q2LS (K) =  &
        FINTRP(1,PBAR(K)*0.01,P1(K1),Q2OBS(K1),P1(K2),Q2OBS(K2))
  200 CONTINUE

!     CHANGE UNIT FROM K/DAY TO K/S (Q1LS) AND G/G/S (Q2LS) AND MULTIPLY
!     THEM BY SCALEQ TO GET THE ACTUALLY IMPOSED Q1 & Q2.

      IF ( Q1Q2 ) SCALEQ = SCALE

      DO 210 K = 1, NK3
      Q1LS (K) = -1.* SCALEQ * Q1LS (K) / PIBAR(K) / SECDAY
      Q2LS (K) = SCALEQ * Q2LS (K) / (GAMFAC * SECDAY)
  210 CONTINUE

#if !defined (RADCODE)
      DO 212 K = 1, NK3
      Q1LS (K) = Q1LS (K) + RAD1 / PIBAR(K) / SECDAY
  212 CONTINUE
#endif

      DO 220 K = 2, NK2
      PBARMB = PBAR(K) * 0.01
      K1 = INDEXR(PBARMB,NU1,P1,LF)
      K2 = K1 + 1
      UG(K) = FINTRP(1,PBARMB,P1(K1),U1(K1),P1(K2),U1(K2))
      VG(K) = FINTRP(1,PBARMB,P1(K1),V1(K1),P1(K2),V1(K2))
  220 CONTINUE 


      PTH='/data5/jsrdoz/11.MODEL/54.VVM/RUNDATA/sst_5400.bin'
      OPEN(UNIT=99,FILE=pth,FORM='unformatted', &
      STATUS='OLD',access='direct',recl=1)
      num=1
      do k=1,MI_GLOB
      READ(99,rec=num) SST3(k)
      num=num+1
      enddo
      CLOSE(99)


!     OCEAN POINTS

      PTH='/data5/jsrdoz/11.MODEL/54.VVM/RUNDATA/sst_5400.bin'
      OPEN(UNIT=99,FILE=PTH,FORM='unformatted', &
      STATUS='OLD',access='direct',recl=1)
      num=1
      do k=1,MI_GLOB
      READ(99,rec=num) SST3(k)
      num=num+1
      enddo
      CLOSE(99)

      DO 300 J = MJM,MJP
      DO 300 I = MIM,MIP
      GWET(I,J)  = -1.
      IQ = NI_SBDM*MI1+I
      JQ = NJ_SBDM*MJ1+J
      TG(I,J) = SST3(IQ) + 273.15

      ZROUGH(I,J)= ZRSEA
  300 CONTINUE

      LOCEAN = .TRUE.
      PSFX = .FALSE.

      WRITE(95,*) '************* INITIAL PROFILE ****************'
      WRITE(95,401)
      WRITE(95,402) (ZT(K),PBAR(K),T(K),THBAR(K),QVBAR(K),UG(K),  &
                     VG(K),K=2,NK2)
      WRITE(95,*) 
      WRITE(95,*) 
      WRITE(95,*) '************* LS FORCING **********************'
      WRITE(95,403)
      WRITE(95,404) (ZT(K),Q1LS(K),Q2LS(K),K=2,NK2)


  401 FORMAT (7X,'Z(M)',2X,'PBAR(MB)',4X,'T(K)',4X,'THBAR(K)',2X,  &
       'QVBAR(G/KG)',2X,'U(M/S)'//)

  402 FORMAT ( 0PF10.0,-2PF10.1,0PF10.2,0PF10.2,3PF10.1,0PF10.1)

  403 FORMAT (7X,'Z(M)',3X,'Q1LS(K/DAY)',1X,'Q2LS(K/DAY)'//)
  404 FORMAT (0PF10.0,0PE12.3,0PE12.3)

      IF(MY_TASK == 0) THEN
      OPEN(UNIT=99,FILE='BAR.DAT',FORM='UNFORMATTED',STATUS='UNKNOWN', &
      ACCESS='DIRECT',RECL=NK3)
      TEMP=THBAR
      WRITE(99,REC=1) TEMP
      TEMP=QVBAR
      WRITE(99,REC=2) TEMP
      TEMP=PBAR
      WRITE(99,REC=3) TEMP
      TEMP=PIBAR
      WRITE(99,REC=4) TEMP
      TEMP=RHO
      WRITE(99,REC=5) TEMP
      CLOSE(99)
      ENDIF


!*********************************************************
!     2. SET UP THE THERMODYNAMIC VARIABLES FOR 3D-FRAME
!*********************************************************
      DO 500 K = 1, NK3
      DO 500 J = MJM,MJP
      DO 500 I = MIM,MIP
      TH3D(I,J,K) = THBAR(K)
      QV3D(I,J,K) = QVBAR(K)
      QC3D(I,J,K) = 0.
      QI3D(I,J,K) = 0.
      QR3D(I,J,K) = 0.
      QS3D(I,J,K) = 0.
      QG3D(I,J,K) = 0.
  500 CONTINUE


      DO 515 K = 1, NK2
      DO 515 J = 1, MJ1
      DO 515 I = 1, MI1
      DO 515 N = 1, 2
      FTH3D(I,J,K,N)=0.
      FQV3D(I,J,K,N)=0.
      FQC3D(I,J,K,N)=0.
      FQI3D(I,J,K,N)=0.
      FQR3D(I,J,K,N)=0.
      FQS3D(I,J,K,N)=0.
      FQG3D(I,J,K,N)=0.
      FZX(I,J,K,N)=0.
      FZY(I,J,K,N)=0.
  515 CONTINUE

      DO 517 N = 1, 2
      DO 517 J = 1, MJ1
      DO 517 I = 1, MI1
      FZTOP(I,J,N)=0.
  517 CONTINUE

      DO 518 N = 1, 2
      FUTMN(N)=0.
      FVTMN(N)=0.
  518 CONTINUE

!*******************************************
!     3. SET UP CIRCULATION
!*******************************************

      DO 600 K = 1, NK2
      DO 600 J = MJM,MJP
      DO 600 I = MIM,MIP
      Z3DX(I,J,K)=0.
      Z3DY(I,J,K)=0.
      W3D(I,J,K) =0.
  600 CONTINUE

      DO 610 K = 1, NK3
      DO 610 J = MJM,MJP
      DO 610 I = MIM,MIP
      Z3DZ(I,J,K)=0.
      U3DX(I,J,K)=0
      U3DY(I,J,K)=0
  610 CONTINUE

      DO 660 K = 2, NK1
      DO 660 J = 1, MJ1
      DO 660 I = 1, MI1
      Z3DY(I,J,K)=(W3D(I+1,J,K)-W3D(I,J,K))/DX    &
                 -(U3DX(I,J,K+1)-U3DX(I,J,K))*FNZ(K)/DZ
      Z3DX(I,J,K)=(W3D(I,J+1,K)-W3D(I,J,K))/DYNEW    &
                 -(U3DY(I,J,K+1)-U3DY(I,J,K))*FNZ(K)/DZ
  660 CONTINUE
      CALL BOUND_3D_XY


#if defined (AREAMN)
!-------------------
!     INSTANT NUDGING
!     OF THE AREA MEAN OF HORIZONTAL WIND TO THE INITIAL VALUE.

      CALL XYAVG2(Z3DX,NHALO,MI1,MJ1,NK2,2,NK1,Z3DX0)
      CALL XYAVG2(Z3DY,NHALO,MI1,MJ1,NK2,2,NK1,Z3DY0)
      CALL XYAVG1(Z3DZ(:,:,NK2), NHALO, MI1, MJ1, Z3DZ0)
!-------------------
#endif

      DO 700 K = 1, NK3
      DO 700 J = 1, MJ1
      DO 700 I = 1, MI1
      DEFXY(I,J,K)=(U3DY(I+1,J,K)-U3DY(I,J,K))/DX   &
                  +(U3DX(I,J+1,K)-U3DX(I,J,K))/DYNEW
  700 CONTINUE
      DO 710 K = 1, NK2
      DO 710 J = 1, MJ1 
      DO 710 I = 1, MI1 
      DEFXZ(I,J,K)=(W3D(I+1,J,K)-W3D(I,J,K))/DX   &
                  +(U3DX(I,J,K+1)-U3DX(I,J,K))*FNZ(K)/DZ
      DEFYZ(I,J,K)=(W3D(I,J+1,K)-W3D(I,J,K))/DYNEW   &
                  +(U3DY(I,J,K+1)-U3DY(I,J,K))*FNZ(K)/DZ
  710 CONTINUE

#if defined (PERIODIC)
      CALL BOUND_ARB (NK3,DEFXY)
      CALL BOUND_ARB (NK2,DEFXZ)
      CALL BOUND_ARB (NK2,DEFYZ)
#endif

      CALL XYAVG1(U3DX(:,:,NK2), NHALO, MI1, MJ1, UTMN)
      CALL XYAVG1(U3DY(:,:,NK2), NHALO, MI1, MJ1, VTMN)
      UTMN0 = UTMN
      VTMN0 = VTMN
      
      IF(ITT == 0) CALL UPDATE_TRACER_INIT
      
   END SUBROUTINE INI_3D


      SUBROUTINE PTSPEC
      REAL(KIND=DBL_KIND) :: &
      HEIGHT,L,H,HT,TMP
      INTEGER(KIND=INT_KIND) :: &
      I,J,K,NN,MIDI,MIDJ
      REAL(KIND=INT_KIND) :: &
      TWTOPO(MI_GLOB,MJ_GLOB) 
      CHARACTER(25) FILENAME,FILENAME2

      ITYPEU = 1
      ITYPEV = 1
      ITYPEW = 1

#if defined (TOPOGRAPHY) 

!      DATA L/75000./,H/4000./ !WIDTH AND HEIGHT OF BELL SHAPED MOUNTAIN. 

      WRITE(FILENAME,'(A12,I3.3,A1,I3.3,A4)') &
      'RUNDATA/TOPO',NI_SBDM+1,'_',NJ_SBDM+1,'.DAT'

      PRINT*,FILENAME

      OPEN(99,FILE=FILENAME,FORM='UNFORMATTED',STATUS='OLD')
      READ(99) ((HX(I,J),I=1,MI1),J=1,MJ1)
      CLOSE(99)

      TMP = MAXVAL(HX)
      CALL MPI_ALLREDUCE(TMP,H,1,MPI_DBL_KIND, MPI_MAX, MPI_COMM_WORLD, MPIERR)
      MAXTOPO=INT(H) + 1

      CALL BOUND_ARB(1,HX)

      DO I=1,MI1
      DO J=1,MJ1
      IF(HX(I+1,J)-HX(I,J) .GT. 0) HXU(I,J)=HX(I+1,J)
      IF(HX(I,J+1)-HX(I,J) .GT. 0) HXV(I,J)=HX(I,J+1)
      ENDDO
      ENDDO

      DO J=1,MJ1
      DO I=1,MI1
      NN=INT(HX(I,J))
      IF(NN .NE. 0) THEN
      DO K=2,NN
      ITYPEU(I,J,K)=0
      ITYPEV(I,J,K)=0
      ITYPEW(I,J,K)=0
      ENDDO
      ENDIF
      ENDDO
      ENDDO

      CALL BOUND_INT_ARB(NK2,ITYPEW)

      DO K=2,NK2
      DO J=MJM,MJP
      DO I=MIM,MIP
      NN=ITYPEW(I,J,K)
      IF(NN .EQ. 0) THEN
      ITYPEU(I-1,J,K)=0
      ITYPEV(I,J-1,K)=0
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      CALL BOUND_INT_ARB(NK2,ITYPEU)
      CALL BOUND_INT_ARB(NK2,ITYPEV)

      DO J=MJM,MJP
      DO I=MIM,MIP
      IF(HX(I,J) .EQ. 0.) HX(I,J)=1.
      ENDDO
      ENDDO


#else
      HX=1
      MAXTOPO=2
#endif

#if defined (DIFFUSION)
      CALL TURB_3D_INI
#if defined (LSM)
      CALL LAND_INI
#endif
#endif

      ENDSUBROUTINE PTSPEC


!^^^^^ENDIF FOR #IF !DEFINED (RESTART)


END MODULE INI_3D_MODULE
