c subroutine for writting out MCMC fit models

      SUBROUTINE wrout(nmax,outpt,in,test,ld,per,fcz,inc,nspots,sp_lat,
     .           sp_lon,sp_ss,sp_fc,xguess,chiout,update_test)

      IMPLICIT NONE
      INCLUDE 'const_incl'

      INTEGER test,in,nspots(9),j,nmax,cntrl,update_test
      REAL(8), DIMENSION(3) :: fcz
      INTEGER, PARAMETER :: iun=4,jun=5
      DOUBLE PRECISION chiout,ld
      REAL(8), DIMENSION(N*9) :: sp_lon,sp_lat,sp_ss,sp_fc
      REAL(8), DIMENSION(ndata) :: xguess
      DOUBLE PRECISION per,inc

      CHARACTER outpt*80 

c open standard output file       
!      IF (in.eq.1) THEN 

!        OPEN(unit=iun,file=trim(outpt)//'.dat',err=999)
!        WRITE(iun,400)
!        OPEN(unit=jun,file=trim(outpt)//'_f.dat',err=998)
!        WRITE(jun,500)
!      ENDIF

!      IF (nmax*0.1.gt.10000) THEN
!          cntrl=int(nmax*0.1)
!      ELSE
!          cntrl=10000
!      ENDIF

!      IF (in.ge.cntrl) THEN

!       print*,'!',trim(outpt),'!'


!        stop

       OPEN(iun,file=trim(outpt)//'.dat',status='old',
     .      position='append',err=997)


!       OPEN(jun,file=trim(outpt)//'_f.dat',status='old',
!     .      position='append',err=996)

       DO j=1,9

         WRITE(iun,404) in,update_test,test,ld,inc,per,fcz(1), 
     .     fcz(2),fcz(3),nspots(j),
     .     sp_lon((j-1)*N+1),sp_lon((j-1)*N+2),sp_lon(j-1)*N+(3),
     .     sp_lon((j-1)*N+4),sp_lon((j-1)*N+5),sp_lon((j-1)*N+6),
     .     sp_lat((j-1)*N+1),sp_lat((j-1)*N+2),sp_lat(j-1)*N+(3),
     .     sp_lat((j-1)*N+4),sp_lat((j-1)*N+5),sp_lat((j-1)*N+6),
     .     sp_ss((j-1)*N+1),sp_ss((j-1)*N+2),sp_ss((j-1)*N+3),
     .     sp_ss((j-1)*N+4),sp_ss((j-1)*N+5),sp_ss((j-1)*N+6),
     .     sp_fc((j-1)*N+1),sp_fc((j-1)*N+2),sp_fc((j-1)*N+3),
     .     sp_fc((j-1)*N+4),sp_fc((j-1)*N+5),sp_fc((j-1)*N+6),
     .            chiout

       ENDDO

!      ENDIF




!       IF (in.gt.5000) THEN
!       WRITE(jun,520) in,chiout
!       DO j=1,30
!        WRITE(jun,510) xguess((j-1)*10+1:j*10)
!       ENDDO
!       ENDIF


!       WRITE(iun,*) '    '

       CLOSE(iun)

!       CLOSE (jun)

      RETURN

 400  FORMAT('# N  Y/N  L.D. Inc Per Spots Lon  Lat SS FC ')
! 401  FORMAT(I7.7,2x,I1.1,2x,F8.3,2x,F8.3,2x,I3.3,2x,F8.3,2x,F8.3,2x,
!     .       F8.2,2x,F8.2,2x,E16.8)
! 402  FORMAT(I7.7,2x,I1.1,2x,F8.3,2x,F8.3,2x,I3.3,2x,F8.3,2x,F8.3,2x,
!     .       F8.3,2x,F8.3,2x,F8.2,2x,F8.2,2x,F8.2,2x,F8.2,2x,E16.8)
! 403  FORMAT(I7.7,2x,I1.1,2x,F8.3,2x,F8.3,2x,I3.3,2x,F8.3,2x,F8.3,2x,
!     .      F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.2,2x,F8.2,2x,F8.2,2x,
!     .      F8.2,2x,F8.2,2x,F8.2,2x,E16.8)
 404  FORMAT(I9.9,2x,I3.3,2x,I1.1,2x,F8.3,2x,f8.3,2x,F8.3,2x,F8.3,2x,
     .  F8.3,2x,F8.3,2x,I3.3,2x,
     .  F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,
     .  F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,
     .  F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,
     .  F8.3,2x,F8.3,2xE16.8)
! 405  FORMAT(I7.7,2x,I1.1,2x,F8.3,2x,F8.3,2x,I3.3,2x,F8.3,2x,F8.3,2x,
!     .       F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,
!     .       F8.3,2x,F8.2,2x,F8.2,2x,F8.2,2x,F8.2,2x,F8.2,2x,F8.2,2x,
!     .       F8.2,2x,F8.2,2x,F8.2,2x,F8.2,2x,E16.8)


! 500  FORMAT('Output fluxes')
! 520  FORMAT (I7.7,2x,E16.8)
! 510  FORMAT(F10.8,2x,F10.8,2x,F10.8,2x,F10.8,2x,F10.8,2x,F10.8,2x,
!     .       F10.8,2x,F10.8,2x,F10.8,2x,F10.8)

 999  WRITE(*,*) 'ERROR OPENING OUTPUT FILE! wrout.f',in 
      STOP
 997  WRITE(*,*) 'ERROR RE-OPENING FILE! wrout.f'
      STOP

 998  WRITE(*,*) 'ERROR OPENING SECONDARY FLUX FILE wrout.f'
 996  WRITE(*,*) 'ERROR RE-OPENING FLUXES output; wrout.f '
      STOP


      END
