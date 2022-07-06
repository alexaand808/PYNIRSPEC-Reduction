C RECLEN.INC ------------------------------------------------------------------
C
C VERSION
C     10-APR-02  AD  Documentation change only.
C     18-JAN-00  AD  Documentation change only.
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Machine-dependent record length for accessing HITRAN binary files.
C     Use  88 for IBC RISC, Sun, Fujitsu f95 linux; 
C          22 for for Compaq Alpha, DEC VAX VMS (longwords).
C
       INTEGER       RECLEN            ! RECL for open
         PARAMETER ( RECLEN = 88 )     ! 22 = Alpha, VMS
C
C -----------------------------------------------------------------------------