
       PROGRAM MCMC
c-------------------------------------------------------------------------
c MCMC for mapping of the observed Brown Dwarf signal. 
c-------------------------------------------------------------------------

       IMPLICIT NONE
       INCLUDE 'const_incl'

       INTEGER in,Nsp,Nmax,test,loc1
       INTEGER update_test,ict
       DOUBLE PRECISION a1,a2,a3,a4,a5,b1,b2,b3,b4,b5

       INTEGER nspots(9),cntrl,nspots_n(9) !fail_c,pass_c,
       DOUBLE PRECISION u !Io,Irotsystem,Isp,Csp,Isp2,Csp2,

       REAL(8), DIMENSION(obsmax) :: xin,sigmain,tin

!       DOUBLE PRECISION lat(1:N),lon(1:N),ss(1:N),fc(1:N),
       REAL(8), DIMENSION (N*9) :: lat,lon,ss,fc
       REAL(8), DIMENSION(ndata) :: xguess,xn,xo,
     .                  tguess,tn

       REAL(8), DIMENSION(3) :: fcz,fcz_n
       DOUBLE PRECISION inc,per,per_n,inc_n,chiout,
     .       ld,ld_n                    !,Irotsym,Crotsym

       CHARACTER fmeas*80,outpt*80 

       REAL(8), DIMENSION(N*9) :: lat_n,ss_n,fc_n,lon_n
       INTEGER, PARAMETER :: cun=8

 
c------------------------------------------------------------------------
c-------------------------------------------------------------------------

!initialize arrays!
       inc=0.
       per=0.
       lon=0.
       lat=0.
       ss=0.
       fc=0.
       xguess=0.
       tguess=0.
       xn=0.
       tn=0.
       xo=0.
       xin=0.
       tin=0.
       sigmain=0.
       lat_n=0.
       lon_n=0.
       ss_n=0.
       fc_n=0.


       a1=0.
       a2=0.
       a3=0.
       a4=0.
       a5=0.
       b1=0.
       b2=0.
       b3=0. 
       b4=0.
       b5=0.


       in=1 
!       in=4582209

       call rdparam(fmeas,Nmax,inc,per,ld,nspots,lon,lat,ss,fc,
     .                    outpt)!Irotsym,Crotsym)


!        fmeas='2m1324_v54.dat'

!       print*,'mcmc:',lat(1),fc(1)

       call rdinput(fmeas,xin,sigmain,tin)


!       sigmain=sigmain*2.!/3.
!       sigmain=0.004  !4 !0.01
!       inc=89.-inc

       tin=tin-tin(1)
       ict=count(tin.gt.0)


       call BDsignal(cntrl,inc,tin,ld,per,fcz,nspots,lon,lat,ss,fc,  !Irotsym,Crotsym, 
     .                     xguess,tguess)


 500  continue     
      


        print*,'in:',in

        per_n=per
        inc_n=inc
        fcz_n=fcz
        nspots_n=nspots
        lat_n=lat
        lon_n=lon
        ss_n=ss
        fc_n=fc
        ld_n=ld

       call gmh_on(xguess,per_n,fcz_n,inc_n,ld_n,nspots_n,tin,lat_n,     !Irotsym,Crotsym,
     .                lon_n,ss_n,fc_n,xn,tn,update_test)


!       print*,xguess
!       stop

       call randnum(u)

       call control(xin,tin,sigmain,xguess,tguess,xn,tn,u,
     .              test,chiout)


!---------------------------------------------------
      IF (in.lt.int(nmax*0.1)) THEN 
        IF (test.eq.1) THEN
         lon=lon_n
         lat=lat_n
         ld=ld_n
         ss=ss_n
         fc=fc_n
         nspots=nspots_n
         per=per_n
         inc=inc_n
         xguess=xn
        ENDIF
      ENDIF



!      print*,nspots,lon(1),lat(1)

!      STOP

!      IF (in.ge.100) THEN 
      IF (in.ge.int(nmax*0.1)) THEN 



!      IF (mod(in,2000).eq.0) THEN

!       OPEN(unit=cun,file='conv_'//trim(outpt)//'.dat',status='old',
!     .      position='append',err=994) 

!        WRITE(cun,*) a1/(a1+b1),a2/(a2+b2),a3/(a3+b3),a4/(a4+b4),
!!     .         a5/(a5+b5),
!     .         (a1+a2+a3+a4)/(a1+a2+a3+a4+b1+b2+b3+b4)
!       CLOSE(cun)
     
!----------------
!      IF (in.eq.int(nmax/2.)) THEN 

!       IF (a1/(a1+b1).gt.0.6) THEN 
!          print*,'control convergence 1'
!        !  STOP
!       ENDIF
!       IF (a2/(a2+b2).gt.0.6) THEN 
!          print*,'control convergence 2'
!         ! STOP
!       ENDIF

!       IF (a3/(a3+b3).gt.0.6) THEN 
!          print*,'control convergence 3'
!          !STOP
!       ENDIF
!       IF (a4/(a4+b4).gt.0.6) THEN 
!          print*,'control convergence 4'
!         ! STOP
!       ENDIF
!       IF (a5/(a5+b5).gt.0.6) THEN 
!          print*,'control convergence 5'
!         ! STOP
!       ENDIF
!      ENDIF 
!___________________
!      ENDIF



!       print*,test,inc,nspots,lon(1),lat(1)
      


       IF (test.eq.0) THEN 

       call wrout(Nmax,outpt,in,test,ld,per,inc,nspots,lat,lon,ss,
     .            fc,xguess,chiout,update_test) 

       ENDIF

      IF (test.eq.1) THEN

        call wrout(Nmax,outpt,in,test,ld_n,per_n,inc_n,nspots_n,lat_n,
     .            lon_n,ss_n,fc_n,xn,chiout,update_test)

         lon=lon_n
         lat=lat_n
         ss=ss_n
         fc=fc_n
         nspots=nspots_n
         per=per_n
         inc=inc_n
         xguess=xn

       ENDIF


      ENDIF

!---------------------------------------------------


       IF (in.le.Nmax) THEN 
         in=in+1
         GO TO 500

       ELSE 

         print*,'end of MCMC; goodbye!' !,fail_c,pass_c

       ENDIF 
      

       RETURN




 505   FORMAT('Convergence 1 2 3 4 5 total')
 994   WRITE(*,*) 'Error opening conv. contr. file @ mcmc.f'


       END









