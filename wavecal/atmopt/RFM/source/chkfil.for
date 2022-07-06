      LOGICAL FUNCTION CHKFIL ( LUN, FILNAM )
C
C VERSION
C     23-APR-12  AD  Original.
C
C DESCRIPTION
C     Check if argument represents an existing filename
C     General purpose RFM module
C
      IMPLICIT NONE 
C
C ARGUMENTS
      INTEGER       LUN    !  I  Spare LUN to be used for opening/closing
      CHARACTER*(*) FILNAM !  I  File name to be tested
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      OPEN ( UNIT=LUN, FILE=FILNAM, STATUS='OLD', ERR=100 )
      CLOSE ( LUN ) 
      CHKFIL = .TRUE.
      RETURN
C
 100  CHKFIL = .FALSE.
C
      END 
