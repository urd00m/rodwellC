      program main 
      IMPLICIT DOUBLE PRECISION (A-H,K-M,O-Z)
      character PRNTR*20
      integer i
      PRNTR = "fortrangoto.txt"
      OPEN(unit=10,FILE=PRNTR,STATUS="NEW")

      DO I=1,10
         WRITE(10, 3030) I
3030     format(1x, I1)
        GOTO 2
1       WRITE(10, 3000)
3000    format(1x, "Back into")
      enddo 

2     WRITE(10, 3030) I
      GOTO 1




      END