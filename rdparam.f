c subroutine to read in the first guess parameters for the MCMC 

       SUBROUTINE rdparam(fmeas,Nmax,inc,per,ld,nspots,lon,lat,
     .                   ss,fc,outpt)!Irotsym,Crotsym)


       IMPLICIT NONE
       INCLUDE 'const_incl'

       INTEGER nspots(9),i,Nmax,j

       DOUBLE PRECISION inc,per,Irotsym,Crotsym,ld

       DOUBLE PRECISION lat(1:N*9),lon(1:N*9),ss(1:N*9),
     .               fc(1:N*9)
       CHARACTER header*46, fmeas*80,outpt*80

!      OPEN(unit=9,file='mcmcin.dat',err=100)


       lat=0.
       lon=0.
       ss=0.
       fc=0.

       OPEN(unit=9,file=
     .  '/extra/tkaralidi/MCinputs/STORMS_INS/2m1821_1.dat',
     .      err=100) 
       READ(9,*) header
       READ(9,*) 


       READ(9,*) fmeas              !file with observations to compare with

       READ(9,*) outpt        !fileaname where output will be written
   
       READ(9,*) Nmax               !max chain length

       READ(9,*) inc                !inclination
       READ(9,*) per                !rotational period
       READ(9,*) ld                !rotational period

       DO i=1,9
       READ(9,*) nspots(i)             !number of spots
       ENDDO

       DO j=1,9

       DO i=1,nspots(j)
       READ(9,*) lon((j-1)*N+i)             !longitudes of spots (degrees)
       ENDDO 

       DO i=1,nspots(j)
       READ(9,*) lat((j-1)*N+i)             !latitudes of spot (degrees)
       ENDDO 

       DO i=1,nspots(j)
       READ(9,*) ss((j-1)*N+i)             !angular size of spot (degrees)
       ENDDO 

       DO i=1,nspots(j)
       READ(9,*) fc((j-1)*N+i)             !flux ratio Fspot/Fbackground
       ENDDO 

       ENDDO

       CLOSE(9)

 



       RETURN

100   write(*,*) 'ERROR READING input file'
       STOP

       END
