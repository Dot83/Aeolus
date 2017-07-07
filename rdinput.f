c subroutine for reading in observations of BDs

      SUBROUTINE rdinput(fmeas,xin,sigmain,tin)

      IMPLICIT NONE 
      INCLUDE 'const_incl'

      INTEGER i
      REAL(8), DIMENSION(obsmax) :: xin,sigmain,tin
      CHARACTER fmeas*80,header*40

      xin=0.
      sigmain=0.
      tin=0.

      OPEN(unit=9,file=trim(fmeas),err=100)
      READ(9,*) header
      READ(9,*) 

      DO i=1,obsmax
       READ(9,*,end=200) xin(i),sigmain(i),tin(i)
      ENDDO




 200  continue

      CLOSE(9)

  !    print*,'readin:',xin

 !      tin=tin*2!-1.

      RETURN

 100  print*,'error reading in observations'
      STOP
 300  print*,'cannot open file?'
      STOP


      END
