      module integrate
      contains
C Authors: Harold Youngren (guppy@maine.rr.com), Mark Drela (drela@mit.edu)
C Program Management: Brad Tousley, Paul Eremenko (eremenko@alum.mit.edu)
C
C=========================================================================

      SUBROUTINE LAMP( X1, R1, X2, R2,  XF, RF, 
     &                 UG1, VG1, UG2, VG2, 
     &                 US1, VS1, US2, VS2 )
C-----------------------------------------------------------------------
C
C     Computes the velocities at XF,RF induced by a vortex + source 
C     sheet "lampshade"  extending from X1,R1 to X2,R2.  
C
C      Vortex sheet density = gamma  (positive counterclockwise)
C      Source sheet density = sigma
C
C     Both densities are assumed to be linear in meridional 
C     arc length over the panel.
C
C
C  Input:
C  ------
C     X1,R1   x,r of one   lampshade edge
C     X2,R2   x,r of other lampshade edge
C     XF,RF   x,r of the field point
C
C  Output: 
C  -------
C     UG1,VG1   x,r velocities at XF,RF for unit gamma at X1,R1
C     UG2,VG2   x,r velocities at XF,RF for unit gamma at X2,R2
C     US1,VS1   x,r velocities at XF,RF for unit sigma at X1,R1
C     US2,VS2   x,r velocities at XF,RF for unit sigma at X2,R2
C
C     Total x,r velocities for given endpoint sheet densities are:
C
C       U = UG1*gamma1 + UG2*gamma2 + US1*sigma1 + US2*sigma2
C       V = VG1*gamma1 + VG2*gamma2 + VS1*sigma1 + VS2*sigma2
C
C-----------------------------------------------------------------------
C
C     Integrates point vortex/source ring velocities over 
C     the lampshade using a Romberg sequence applied to the 
C     simple midpoint rule:
C
C        |       *       |  ->  I1  -   I21  -   I321  -   I4321
C                                   /        /         /     .
C        |   *       *   |  ->  I2  -   I32  -   I432        .
C                                   /        /     .       8th order
C        | *   *   *   * |  ->  I3  -   I43        .       
C                                   /    .       6th order
C        |* * * * * * * *|  ->  I4       .        
C                .               .      4th order
C                .               .      
C                               2nd order
C               etc             
C
C  The first column I1,I2... are the 2nd-order integral approximations
C  computed with the stock midpoint rule.  The subsequent columns
C  are the Richardson extrapolations to higher-order accuracy.
C
C
C  Algorithm for four Romberg stages:
C
C    I21 = (4*I2 - I1) / 3           | extrapolation of 2nd-order results
C    I32 = (4*I3 - I2) / 3           |
C    I43 = (4*I4 - I3) / 3           |
C 
C    I321 = (16*I32 - I21) / 15      | extrapolation of 4th-order results
C    I432 = (16*I43 - I32) / 15      |
C
C    I4321 = (64*I432 - I321) / 63   | extrapolation of 6th-order results
C
C
C  Quantities stored in the NROMX arrays after each IROM stage:
C
C    IROM  =   1      2       3       4     
C            ----   ----    -----   ------
C    U(1)  =  I1     I21     I321    I4321  
C    U(2)  =         I2      I32     I432   
C    U(3)  =                 I3      I43    
C    U(4)  =                         I4      
C
C-----------------------------------------------------------------------
C
C---- NROMX = max number of Romberg stages
C-      This limits the integration resolution and cost.
C-      The influence of this panel on a field point 
C-      which is closer than  O(panel_length)/2**NROMX 
C-      will not be accurately represented.
      PARAMETER (NROMX=10)
C
      DIMENSION UG1I(NROMX), VG1I(NROMX),
     &          UG2I(NROMX), VG2I(NROMX),
     &          US1I(NROMX), VS1I(NROMX),
     &          US2I(NROMX), VS2I(NROMX)
C
      PARAMETER(PI=3.14159265358979)
C
C
C---- Romberg convergence tolerance
      DATA ROMTOL / 1.0E-6 /
ccc      DATA ROMTOL / 1.0E-12 /
C
C---- reference length for convergence tolerance 
      REFL = 0.5*(R1 + R2)
C
C---- evaluate integrals over  0..t..1  on increasingly fine grids
      DO 100 IROM=1, NROMX
       NT = 2**IROM / 2
C
       UG1I(IROM) = 0.
       VG1I(IROM) = 0.
       UG2I(IROM) = 0.
       VG2I(IROM) = 0.
       US1I(IROM) = 0.
       VS1I(IROM) = 0.
       US2I(IROM) = 0.
       VS2I(IROM) = 0.
C
C------ visit the midpoints of each of the NT intervals
       DO 10 IT=1, NT
         T = (FLOAT(IT) - 0.5) / FLOAT(NT)
         TB = 1.0 - T
C
         DT = 1.0/FLOAT(NT)
C
         XT = X1*TB + X2*T
         RT = R1*TB + R2*T
C
C-------- get induced velocities for vortex,source ring at XT,RT
         CALL RING( XT, RT, XF, RF, UGT, VGT, UST, VST )
C
C-------- accumulate the separate unit-gamma, unit-sigma integrals
         UG1I(IROM) = UG1I(IROM) + DT*UGT*TB
         VG1I(IROM) = VG1I(IROM) + DT*VGT*TB
         UG2I(IROM) = UG2I(IROM) + DT*UGT*T
         VG2I(IROM) = VG2I(IROM) + DT*VGT*T
         US1I(IROM) = US1I(IROM) + DT*UST*TB
         VS1I(IROM) = VS1I(IROM) + DT*VST*TB
         US2I(IROM) = US2I(IROM) + DT*UST*T
         VS2I(IROM) = VS2I(IROM) + DT*VST*T
 10     CONTINUE
C
C------ Romberg sequence using all previous grid results
       DO KROM = IROM, 2, -1
C-------- weight needed to cancel lowest-order error terms in KROM level
         W = 2.0 ** (2*(IROM-KROM+1))
C
C-------- put Richardson extrapolation for KROM level into KROM-1 level
         UG1I(KROM-1) = (W*UG1I(KROM) - UG1I(KROM-1)) / (W-1.0)
         VG1I(KROM-1) = (W*VG1I(KROM) - VG1I(KROM-1)) / (W-1.0)
         UG2I(KROM-1) = (W*UG2I(KROM) - UG2I(KROM-1)) / (W-1.0)
         VG2I(KROM-1) = (W*VG2I(KROM) - VG2I(KROM-1)) / (W-1.0)
         US1I(KROM-1) = (W*US1I(KROM) - US1I(KROM-1)) / (W-1.0)
         VS1I(KROM-1) = (W*VS1I(KROM) - VS1I(KROM-1)) / (W-1.0)
         US2I(KROM-1) = (W*US2I(KROM) - US2I(KROM-1)) / (W-1.0)
         VS2I(KROM-1) = (W*VS2I(KROM) - VS2I(KROM-1)) / (W-1.0)
       ENDDO
C
       IF(IROM.GT.1) THEN
C------- compare the best-current and best-previous integrals
        ERRUG1 = UG1I(1) - UG1I(2)
        ERRVG1 = VG1I(1) - VG1I(2)
        ERRUG2 = UG2I(1) - UG2I(2)
        ERRVG2 = VG2I(1) - VG2I(2)
        ERRUS1 = US1I(1) - US1I(2)
        ERRVS1 = VS1I(1) - VS1I(2)
        ERRUS2 = US2I(1) - US2I(2)
        ERRVS2 = VS2I(1) - VS2I(2)
C
         ERR = MAX(
     &     ABS(ERRUG1),
     &     ABS(ERRVG1),
     &     ABS(ERRUG2),
     &     ABS(ERRVG2),
     &     ABS(ERRUS1),
     &     ABS(ERRVS1),
     &     ABS(ERRUS2),
     &     ABS(ERRVS2) )
        
c         write(13,1200) 
c     &     ABS(ERRUG1),
c     &     ABS(ERRVG1),
c     &     ABS(ERRUG2),
c     &     ABS(ERRVG2),
c     &     ABS(ERRUS1),
c     &     ABS(ERRVS1),
c     &     ABS(ERRUS2),
c     &     ABS(ERRVS2)
c 1200    format(8e16.9)

        IF(ERR*REFL .LT. ROMTOL) GO TO 101
       ENDIF
 100  CONTINUE
      WRITE(*,*) 'LAMP: Romberg convergence failed.  Error =', ERR
C
 101  CONTINUE

ccc      write(*,*) IROM, ERR
C
C---- return best final results
      DELSQ = (X1-X2)**2 + (R1-R2)**2
      DELS = SQRT(DELSQ)
C
      UG1 = UG1I(1)*DELS
      VG1 = VG1I(1)*DELS
      UG2 = UG2I(1)*DELS
      VG2 = VG2I(1)*DELS
      US1 = US1I(1)*DELS
      VS1 = VS1I(1)*DELS
      US2 = US2I(1)*DELS
      VS2 = VS2I(1)*DELS
C
      RETURN
      END ! LAMP



      SUBROUTINE LAMPC( X1, R1, X2, R2,
     &                  UG1, VG1, UG2, VG2, 
     &                  US1, VS1, US2, VS2)
C     &                  UGs1,VGs1,UGs2,VGs2,
C     &                  UGa1,VGa1,UGa2,VGa2)
C---------------------------------------------------------
C     Same as LAMP, but the field point is assumed 
C     to be at the lampshade-panel midpoint.
C
C     The 1/r and log(r) singularities in the integrands 
C     are removed from the numerical integration,
C     and are computed analytically.
C
C     The induced velocites returned by this routine 
C     are the average of the two values on each side 
C     of the sheet.  The sheet jumps are not included.
C---------------------------------------------------------
C
C---- max number of Romberg integration stages
      PARAMETER (NROMX=9)
C
      DIMENSION UG1I(NROMX), VG1I(NROMX),
     &          UG2I(NROMX), VG2I(NROMX),
     &          US1I(NROMX), VS1I(NROMX),
     &          US2I(NROMX), VS2I(NROMX)
C
      PARAMETER(PI=3.14159265358979)
C
C
C---- Romberg convergence tolerance (actual error may be much less than this)
      !DATA ROMTOL / 1.0E-6 /
      DATA ROMTOL / 1.0E-12 /
C
C---- lampshade meridional length**2
      DELSQ = (X1-X2)**2 + (R1-R2)**2
C
C---- reference length for convergence tolerance 
      REFL = 0.5*(R1 + R2)
C
C
C---- field point is assumed to be at midpoint
      XF = 0.5*(X1 + X2)
      RF = 0.5*(R1 + R2)
C
C---- evaluate integrals on increasingly fine grids 
C-     (start with two intervals to avoid landing right on the midpoint)
      DO 100 IROM=1, NROMX
       NT = 2**IROM
C
       UG1I(IROM) = 0.
       VG1I(IROM) = 0.
       UG2I(IROM) = 0.
       VG2I(IROM) = 0.
       US1I(IROM) = 0.
       VS1I(IROM) = 0.
       US2I(IROM) = 0.
       VS2I(IROM) = 0.
C
C------ visit the midpoints of each of the NT intervals
       DO 10 IT=1, NT
         T = (FLOAT(IT) - 0.5) / FLOAT(NT)
         PRINT*, "T: ", T
         TB = 1.0 - T
C
         DT = 1.0/FLOAT(NT)
C
         XT = X1*TB + X2*T
         RT = R1*TB + R2*T
C
         CALL RING( XT, RT, XF, RF, UGT, VGT, UST, VST )
C
C-------- singular parts of velocities in the limit  XT,RT -> XF,RF
         DSQ = (XT-XF)**2 + (RT-RF)**2
         UGA =  (RT-RF)/(4.0*PI*DSQ) ! vortex axial singular part from E term
     &        - 0.5*LOG(DSQ/(64.0*RF**2)) / (8.0*PI*RF) ! vortex axial singular part from K term
         VGA = -(XT-XF)/(4.0*PI*DSQ) ! vortex radial singular part from E term
         USA = -(XT-XF)/(4.0*PI*DSQ) ! source axial singular part from E term
         VSA = -(RT-RF)/(4.0*PI*DSQ) ! source radial singular part from E term
     &        - 0.5*LOG(DSQ/      RF**2 ) / (8.0*PI*RF) ! source radial singular part from K term (where did 64 go?)
C
C-------- accumulate integrals, with singular parts (at t=0.5) removed
         UG1I(IROM) = UG1I(IROM) + DT*(UGT*TB - UGA)
         VG1I(IROM) = VG1I(IROM) + DT*(VGT*TB - VGA)
         UG2I(IROM) = UG2I(IROM) + DT*(UGT*T  - UGA)
         VG2I(IROM) = VG2I(IROM) + DT*(VGT*T  - VGA)
         US1I(IROM) = US1I(IROM) + DT*(UST*TB - USA)
         VS1I(IROM) = VS1I(IROM) + DT*(VST*TB - VSA)
         US2I(IROM) = US2I(IROM) + DT*(UST*T  - USA)
         VS2I(IROM) = VS2I(IROM) + DT*(VST*T  - VSA)
 10     CONTINUE
C
C------ Romberg sequence using all previous grid results
       DO KROM = IROM, 2, -1
C-------- weight needed to cancel lowest-order error terms in KROM level
         W = 2.0**(2*(IROM-KROM+1))
         WM1 = W - 1.0
C
C-------- put Richardson extrapolation for KROM level into KROM-1 level
         UG1I(KROM-1) = (W*UG1I(KROM) - UG1I(KROM-1)) / WM1
         VG1I(KROM-1) = (W*VG1I(KROM) - VG1I(KROM-1)) / WM1
         UG2I(KROM-1) = (W*UG2I(KROM) - UG2I(KROM-1)) / WM1
         VG2I(KROM-1) = (W*VG2I(KROM) - VG2I(KROM-1)) / WM1
         US1I(KROM-1) = (W*US1I(KROM) - US1I(KROM-1)) / WM1
         VS1I(KROM-1) = (W*VS1I(KROM) - VS1I(KROM-1)) / WM1
         US2I(KROM-1) = (W*US2I(KROM) - US2I(KROM-1)) / WM1
         VS2I(KROM-1) = (W*VS2I(KROM) - VS2I(KROM-1)) / WM1
       ENDDO
C
       IF(IROM.GT.1) THEN
C------- compare the best-current and best-previous integrals
        ERRUG1 = UG1I(1) - UG1I(2)
        ERRVG1 = VG1I(1) - VG1I(2)
        ERRUG2 = UG2I(1) - UG2I(2)
        ERRVG2 = VG2I(1) - VG2I(2)
        ERRUS1 = US1I(1) - US1I(2)
        ERRVS1 = VS1I(1) - VS1I(2)
        ERRUS2 = US2I(1) - US2I(2)
        ERRVS2 = VS2I(1) - VS2I(2)
C
        ERR = MAX(
     &     ABS(ERRUG1),
     &     ABS(ERRVG1),
     &     ABS(ERRUG2),
     &     ABS(ERRVG2),
     &     ABS(ERRUS1),
     &     ABS(ERRVS1),
     &     ABS(ERRUS2),
     &     ABS(ERRVS2) )
C 
        IF(ERR*REFL .LT. ROMTOL) GO TO 101
       ENDIF
 100  CONTINUE
      WRITE(*,*) 'LAMPC: Romberg convergence failed.  Error =', ERR
C
 101  CONTINUE
C
      DELSQ = (X1-X2)**2 + (R1-R2)**2
      DELS = SQRT(DELSQ)
C
C---- analytically-integrated singular parts which were removed
      UGAI = (1.0 + LOG(16.0*RF/DELS)) / (4.0*PI*RF)
      VGAI = 0.
      USAI = 0.
      VSAI = (1.0 + LOG( 2.0*RF/DELS)) / (4.0*PI*RF)
      ! open(unit=77,file="singular_components.jl")
        write(77,*) '# Analytic Components'
        write(77,*) 'VGa = [',UGAI,';',VGAI,']'
        write(77,*) 'VSa = [',USAI,';',VSAI,']'
        write(77,*) 'VGmsi = [',UG1I(1),';',VG1I(1),']'
        write(77,*) 'VGmsip1 = [',UG2I(1),';',VG2I(1),']'
        write(77,*) 'VSmsi = [',US1I(1),';',VS1I(1),']'
        write(77,*) 'VSmsip1 = [',US2I(1),';',VS2I(1),']'     

C
C---- return final results, with removed parts added back on
      UG1 = (UG1I(1) + UGAI*0.5)*DELS
      VG1 = (VG1I(1) + VGAI*0.5)*DELS
      UG2 = (UG2I(1) + UGAI*0.5)*DELS
      VG2 = (VG2I(1) + VGAI*0.5)*DELS
      US1 = (US1I(1) + USAI*0.5)*DELS
      VS1 = (VS1I(1) + VSAI*0.5)*DELS
      US2 = (US2I(1) + USAI*0.5)*DELS
      VS2 = (VS2I(1) + VSAI*0.5)*DELS
C
      close(77)
      RETURN
      END ! LAMPC

      SUBROUTINE RING( XV, RV, XF, RF, UX, UR, SX, SR )
C-----------------------------------------------------------------------
C     Computes the velocities induced by a vortex/source ring 
C     located at XV with radius RV with unit circulation 
C     and unit source/length density.
C
C     Adapted from routines provided by J. Kerwin.
C-----------------------------------------------------------------------
C  Input:
C     XV,RV  x,r of the ring
C     XF,RF  x,r of the field point
C
C  Output:
C     UX,UR  velocity at XF,RF for unit circulation (+ counterclockwise)
C     SX,SR  velocity at XF,RF for unit source/perimeter
C-----------------------------------------------------------------------
      PARAMETER(PI=3.14159265358979)
C
      IF(RV.LE.0.0) THEN
C------ zero-radius ring
          UX = 0.
          UR = 0.
          SX = 0.
          SR = 0.
          RETURN
      ENDIF
C
C---- this fails if R=1 and X=0  (on the ring itself)
      R = RF/RV
      X = (XV-XF)/RV
C
      IF(R.EQ.1.0 .AND. X.EQ.0.0) THEN
      UX = 0.
      UR = 0.
      SX = 0.
      SR = 0.
      RETURN
      ENDIF
C 
      IF (RF .EQ. 0.0) THEN
C----- Control Point on the axis
      UX = 1.0 / SQRT(1.0+X**2)**3 / (2.0*RV)
      UR = 0.0
      SX =  -X / SQRT(1.0+X**2)**3 / (2.0*RV)
      SR = 0.0
C
      ELSE
C----- Control Point not on X-axis
      XRP = X**2 + (1.0+R)**2
      XRM = X**2 + (1.0-R)**2
C
      SRP = SQRT(XRP)
C
      AK = XRM/XRP
      CALL ELLEK(AK,ELE,ELK)
      PRINT*, 'E: ', ELE
      PRINT*, 'K: ', ELK
C
      F = 2.0/XRM
C
      UX = ( 1.0/ SRP   )*(ELK - ELE*(1.0 + F*(R-1.0)))/(2.0*PI*RV)
      UR = (   X/(SRP*R))*(ELK - ELE*(1.0 + F* R     ))/(2.0*PI*RV) ! where is the negative out front?
C
      SX =  (  X/ SRP   )*(    - ELE*       F         )/(2.0*PI*RV) ! I don't have a negative for this one.
      SR =  (1.0/(SRP*R))*(ELK - ELE*(1.0 + F*(R-R*R)))/(2.0*PI*RV)
C
ccC----- streamfunction due to vortex
cc       PSI = ((1.0 - 2.0*XRP*R)*ELK - ELE)*RV / (2.0*PI*SQRT(XRP))
      ENDIF
C
      RETURN
      END


      SUBROUTINE ELLEK(AK,ELE,ELK)
C-----------------------------------------------------------------------
C              ELLIPTIC FUNCTIONS ROUTINE
C
C     Adapted from routines provided by J. Kerwin.
C-----------------------------------------------------------------------
C   Input
C     AK     elliptic-integral argument
C
C   Output
C     ELE    complete elliptic integral of the second kind
C     ELK    complete elliptic integral of the first  kind
C_______________________________________________________________________
C
      ALK = -LOG(AK)
C
      ELE = 1.00000000000
     &    +(0.44325141463
     &    +(0.06260601220
     &    +(0.04757383546
     &    + 0.01736506451*AK)*AK)*AK)*AK
     &  +( (0.24998368310
     &    +(0.09200180037
     &    +(0.04069697526
     &    + 0.00526449639*AK)*AK)*AK)*AK )*ALK
C          
      ELK = 1.38629436112
     &    +(0.09666344259
     &    +(0.03590092383
     &    +(0.03742563713
     &    + 0.01451196212*AK)*AK)*AK)*AK
     &  +(  0.50000000000
     &    +(0.12498593597
     &    +(0.06880248576
     &    +(0.03328355346
     &    + 0.00441787012*AK)*AK)*AK)*AK )*ALK
C
      RETURN
      END   

      end module integrate