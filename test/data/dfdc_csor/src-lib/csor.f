      module csor
      contains
C=========================================================================

      SUBROUTINE CIRCRLX(BGX, BGAM, BGMAG, DBGOLD, NRC)

          DIMENSION BGX(NRC), BGAM(NRC), DBGOLD(NRC)
C
C---- Check for rational relaxation factors based on BGAM changes
      RLXB = 0.5
      IRMAX = 0
      DBGMAX = 0.0
      DO IR = 1, NRC
       DBG = (BGX(IR)-BGAM(IR))
       IF(ABS(DBG).GT.ABS(DBGMAX)) THEN
         DBGMAX = DBG
         IRMAX = IR
       ENDIF
       IF(BGMAG.NE.0.0) THEN
         FDBG = DBG/BGMAG
         IF(FDBG*RLXB.LT.-0.2) RLXB = -0.2/FDBG
         IF(FDBG*RLXB.GT. 0.4) RLXB =  0.4/FDBG
       ENDIF
      END DO
C ----- TODO: save RLXB value
      write(99,*) 'RLXB = ', RLXB
C
C---- Update blade circulation using CSOR
      write(99,*) 'RLXBG = ['
      DO IR = 1, NRC
        DBG = (BGX(IR)-BGAM(IR))
        RLXBG = 0.5*RLXB
        IF(DBG*DBGOLD(IR).LT.0.0) RLXBG = 0.6*RLXB
        write(99,*) RLXBG
        DBGOLD(IR) = DBG*RLXBG
        BGAM(IR) = BGAM(IR) + RLXBG*DBG
      END DO
      write(99,*) ']'
C---------TODO: save RLXBG as well a BGAM
      write(99,*) 'BGAMupdate = [', BGAM, ']'
      RETURN
      END

      SUBROUTINE GAMWRLX(GAMTH, GTH, DGOLD, NPTOT)

          DIMENSION GAMTH(NPTOT), GTH(NPTOT), DGOLD(NPTOT)

C---- Generate GTH estimate for updated circulations
         RLX = 0.5
         IPMAX = 0
         DGTHMAX = 0.0
         write(99,*) 'RLXG = ['
         DO IP = 1, NPTOT
           DG = GAMTH(IP) - GTH(IP)
           IF(ABS(DG).GT.ABS(DGTHMAX)) THEN
             DGTHMAX = DG
             IPMAX = IP
           ENDIF
           RLXG = RLX
           IF(DG*DGOLD(IP).LT.0.0) RLXG = 0.6*RLX
           IF(DG*DGOLD(IP).GT.0.0) RLXG = 1.2*RLX
C--------TODO: save RLXG value
           write(99,*) RLXG
           DGOLD(IP) = DG*RLXG
C---- Update GTH wake gamma using CSOR
           GTH(IP) = GTH(IP) + RLXG*DG
         END DO
         write(99,*) ']'
C ------ TODO: save DGTHMAX value
         write(99,*) 'DGTHMAX = ', DGTHMAX
C--------TODO: save GTH values
         write(99,*) 'GTHupdate = [', GTH(:), ']'
C
      RETURN
      END

      end module csor
