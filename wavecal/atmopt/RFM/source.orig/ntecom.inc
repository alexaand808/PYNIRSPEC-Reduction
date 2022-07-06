C NTECOM.INC -------------------------------------------------------------------
      COMMON / NTECOM / NNTE, MNTE, IDGNTE, IDXNTE, ISONTE, IGQNTE, 
     &                  QPRNTE,
     &                  ENGNTE, TEMNTE, QFNNTE
C
C VERSION
C     09-JAN-08  AD  Add IDXNTE
C     09-AUG-00  AD  Comment change only
C     04-APR-97  AD  Add QPRNTE(MAXNTE)
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Common variables associated with RFM non-LTE calculations
C     Loaded by ATMNTE, NTETEM, Read by NTECLC.
C
      INTEGER NNTE           ! No. of nonLTE vib.temp profiles
      INTEGER MNTE           ! No. of nonLTE VT profiles loaded via ATM files
      INTEGER IDGNTE(MAXNTE) ! Nonlte gas ids (NTEGAS)
      INTEGER IDXNTE(MAXNTE) ! Encoded IDG+ISO+IGQ as single integer
      INTEGER ISONTE(MAXNTE) ! Nonlte isotope ids (NTEISO)
      INTEGER IGQNTE(MAXNTE) ! Global Quantum index of NTE levels (NTENTE)
      LOGICAL QPRNTE(MAXNTE) ! T=use QFNNTE profile, F=use LTE calc. of Qtot.
      REAL    ENGNTE(MAXNTE) ! Nonlte NTE level energy [cm-1]
      REAL    TEMNTE(MAXATM,MAXNTE) ! Nonlte vibrational-kineter temps (NTETEM)
      REAL    QFNNTE(MAXATM,MAXNTE) ! Nonlte vibrational partition functions
C
C ------------------------------------------------------------------------------
