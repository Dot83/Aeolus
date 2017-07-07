!
       SUBROUTINE gmh_on(xguess,per,incin,ldin,inspots,tin,latin,nrot,
     .                lonin,ssin,fcin,fmod,tmod,update_test)


       IMPLICIT NONE

       INCLUDE 'const_incl'

       INTEGER i,k,guess_spot,nspots,update_test, nrot,
     .       cnt,cntrl,itst, rot_test, 
     .       itstlim,upd, PCA!,rossby

       DOUBLE PRECISION beta,qb,qjump,r,r2,lim1,lim2,lim3,lim4
       REAL(8), DIMENSION(obsmax) ::tin

       DOUBLE PRECISION pi
       PARAMETER (pi=3.14159265359)

       DOUBLE PRECISION per,inc,incin,ld,ldin!,maxfc,minfc
       REAL(8), DIMENSION(N*9) :: latin,ssin,fcin,lonin,lttm,
     .       lntm,sstm,fctm 
       REAL(8), DIMENSION(N) ::  lat,lon,ss,fc,rb
       REAL(8), DIMENSION(N) :: lonintm,ssintm,fcintm,latintm
     
       INTEGER inspots(9),nsptm(9)

       REAL(8), DIMENSION(ndata) :: fmod,tmod,ch2,xguess

       DOUBLE PRECISION beta_inc,beta_nsp,beta_ll,beta_s,
     .         beta_f,beta_i,beta_ld
       DOUBLE PRECISION guess_in,guess_lat!,

       DOUBLE PRECISION new_rule(4)

       
       
! We use a gaussian candidate transition function q(x|xguess)~exp(-(x-xguess)^2/2β)/(sqrt(2π)β)
! where the step β can vary per parameter we change
! the optimization of the steps (βi) for every parameter i, has to be done 
! manually after running a couple of chains and controling the acceptance rates

 77    continue


      upd=0
      itstlim=150

 38    upd=upd+1

      IF (upd.gt.250) THEN
       print*,'stuck'
       print*,update_test,nspots,lonin,latin
       STOP
      ENDIF

      lat=0.
      lon=0.
      ss=0.
      fc=0.
      itst=0
      inc=0.
      ld=0.
      latintm=0
      lonintm=0.
      ssintm=0.
      fcintm=0.


      beta_ld=1.8
      beta_i=1.48
      beta_nsp=2.6 
      beta_ll=1.24
      beta_s=3.4
      beta_f=1.32

      beta_ld=1.8! 1.8 !apo 1.4  H: 1.2  J:1.4
      beta_i=1.6!1.4!1.6! 1.4
      beta_nsp=2.6 
      beta_ll=1.24
      beta_s=3.4
      beta_f=1.32


      lim1=0.01
      lim2=0.01
      lim3=0.01
      lim4=0.01

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c If PCA results show one spot then PCA=1 
       PCA=1
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c Define the maximum contrast ratio you allow the code to fit (*100)
!       maxfc=240.           
!       minfc=101.
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



!c the parameters should be updated at random in sets we split: inclination+# of spots
!c latitude+longitude & contrast+size 
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!!!!!!TEMP FIX FOR lat=[-90 +90] :

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       call ftst(r)
       rot_test=int(r*9)+1


       call ftst(r) 
       update_test=int(r*(3*nspots+3))+1

!
!       IF (update_test.eq.1) print*,'oops'
!       IF (update_test.eq.1) GOTO 2

 !     stop

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------


      IF (update_test.eq.1) THEN

 11     inc=incin

       call ftst(r) 
 !!      inc=int(r*22.)*2.
 !!      inc=89.-inc

 !!!      inc=int(r*89)*2.-89.
 !      inc=int(r*60)*2.-60.

        inc=int(r*40)*2.
 !!!      inc=int(r*50)*2.-42.

!       print*,'in:',r,inc
!        stop

!       IF (inc.gt.90) inc=90.
!       IF (inc.gt.90.) print*,'error inc:gmhon.f'
!       IF (inc.gt.90) STOP


!add an extra control to make sure the new random point is not *exactly* the same as the initial one
     
       IF (inc.eq.incin) GO TO 11
 
       call BDsignal(cntrl,inc,tin,ldin,per,nspots,lonin,
     .         latin,ssin,fc,Fmod,tmod)

       IF (cntrl.eq.1) GO TO 11

!c calculate the transition probability to the new state:

       beta=beta_i

       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))
   
       IF (qjump.lt.lim4) THEN     !0.1

        itst=itst+1

        IF (itst.eq.itstlim) GO TO 38

          GO TO 11

        ELSE

          CONTINUE

        ENDIF 

        incin=inc
!        print*,incin
!        stop

      ENDIF





      IF (update_test.eq.2) THEN

 10      ld=ldin

       call ftst(r) 
       ld=(r*0.79)+0.01


!add an extra control to make sure the new random point is not *exactly* the same as the initial one
     
       IF (ld.eq.ldin) GO TO 10
 
       call BDsignal(cntrl,incin,tin,ld,per,nspots,lonin,
     .         latin,ssin,fc,Fmod,tmod)

       IF (cntrl.eq.1) GO TO 10

!c calculate the transition probability to the new state:

       beta=beta_ld

       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))
   
       IF (qjump.lt.lim4) THEN     !0.1

        itst=itst+1

        IF (itst.eq.itstlim) GO TO 38

          GO TO 10

        ELSE

          CONTINUE

        ENDIF 

        ldin=ld
!        print*,incin
!        stop

      ENDIF



!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------

      IF (update_test.eq.3) THEN

       nspots=inspots(rot_test)

        rb=0.

 22   call ftst(r)           !call random number generator    
      guess_spot=int(r*N)+1                          !set the # spot guess for the model


      latintm=latin((rot_test-1)*N+1:N*rot_test)
      lonintm=lonin((rot_test-1)*N+1:N*rot_test)
      ssintm=ssin((rot_test-1)*N+1:N*rot_test)
      fcintm=fcin((rot_test-1)*N+1:N*rot_test)
      
       lttm=latin
       lntm=lonin
       sstm=ssin
       fctm=fcin
       nsptm=inspots

!add an extra control to make sure the new random point is not *exactly* the same as the initial one
      IF (guess_spot.eq.nspots) GO TO 22

!c If number of spots changes we need to change the available lat,lon,ss,fc :

       IF (guess_spot.lt.nspots) THEN
     
        lat=latintm
        lon=lonintm
        ss=ssintm
        fc=fcintm

        lat(guess_spot+1:nspots)=0.
        lon(guess_spot+1:nspots)=0.
        ss(guess_spot+1:nspots)=0.
        fc(guess_spot+1:nspots)=0.
      
       ENDIF

       IF (guess_spot.gt.nspots) THEN
     
        lat(1:nspots)=latintm((rot_test-1)*N+1:(rot_test-1)*N+nspots)
        lon(1:nspots)=lonintm((rot_test-1)*N+1:(rot_test-1)*N+nspots)
        ss(1:nspots)=ssintm((rot_test-1)*N+1:(rot_test-1)*N+nspots)
        fc(1:nspots)=fcintm((rot_test-1)*N+1:(rot_test-1)*N+nspots)

       DO k=nspots+1,guess_spot
        call ftstb(rb,guess_spot-nspots)                           !call random number generator    
  
!        guess_lat=(int(rb(k-nspots)*45.)+1)*2. !(int(rb(k-nspots)*10.)+1)*9. !45.)+1)*2.    
        guess_lat=-60.+(int(rb(k-nspots)*74.)+1)*2.          
        lat(k)=guess_lat 

!        lat(1)=latin(1)
!        lat(2)=latin(2)
!        lat(3)=latin(3)
!        lat(4)=latin(4)


        call ftstb(rb,guess_spot-nspots)                       !call random number generator    
        lon(k)=int(rb(k)*180.)*2. ! int(rb(k)*360.)            !set the longitude guess for the model (between 0 and 360)

        call ftstb(rb,guess_spot-nspots)                       !call random number generator    
        ss(k)=int(rb(k-nspots)*40)+1.                          !use rn to pick a random size for the spot betwen 1 and 40 deg 
        IF (ss(k).lt.1.)   ss(k)=1.
        IF (ss(k).gt.41.)  ss(k)=41.


          call ftstb(rb,guess_spot-nspots)                                     !call random number generator    
          


!!          fc(k)=(int(rb(k-nspots)*maxfc)+1)/100.+30/100.               !use rn to pick a random contrast from 0.01 to maxfc/100.
!!          fc(k)=(int(rb(k-nspots)*(maxfc+1-minfc))+minfc)/100.

! this is what it was before: 
!          fc(k)=(int(rb(k-nspots)*(200+1-60))+60)/100.
           fc(k)=(int(rb(k-nspots)*(200+1-60))+60)/100.

          IF (fc(k).eq.1.) fc(k)=1.01


       ENDDO

       ENDIF

        nsptm(rot_test)=guess_spot
        lttm((rot_test-1)*N+1:N*rot_test)=lat
        lntm((rot_test-1)*N+1:N*rot_test)=lon
        sstm((rot_test-1)*N+1:N*rot_test)=ss
        fctm((rot_test-1)*N+1:N*rot_test)=fc
      

!now that we fixed the number of spots and their properties we check the signal 

       call BDsignal(cntrl,incin,tin,ldin,per,nsptm,lntm,lttm,
     .               sstm,fctm,Fmod,tmod)

       IF (cntrl.eq.1) GO TO 22

!c calculate the transition probability to the new state:

       beta=beta_nsp


       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-(Fmod-xguess)**2./
     .              (2*beta**2.),mask=fmod.ne.0.))
    

      IF (qjump.lt.lim1) THEN

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 22

      ELSE

       CONTINUE

      ENDIF 
 

      inspots=nsptm
      latin=lttm
      lonin=lntm
      ssin=sstm
      fcin=fctm

      ENDIF

!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------

      IF (update_test.eq.4) THEN

      rb=0.
             

 331  lonintm(1)=lonin((rot_test-1)*N+1)
      latintm(1)=latin((rot_test-1)*N+1)

      lttm=latin
      lntm=lonin

       call ftst(r)               
       lat(1)=-60.+(int(r*74.)+1)*2 !(int(r*45.)+1)*2.

       call ftst(r)  
       lon(1)=int(r*180.)*2. 



!        print*,'in:',latin(1),fcin(1)


!        print*,lat(1),maxfc,minfc,fcin(1)

!         stop

!add an extra control to make sure the new random point is not *exactly* the same as the initial one

       IF ((lat(1).eq.latintm(1)).AND.(lon(1).eq.lonintm(1))) GO TO 331


       lttm((rot_test-1)*N+1)=lat(1)
       lntm((rot_test-1)*N+1)=lon(1)


       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lntm,lttm,
     .               ssin,fcin,Fmod,tmod)

       IF (cntrl.eq.1) GO TO 331

!c calculate the transition probability to the new state:

       beta=beta_ll


       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))
       
!c call a random number-> if qjump>ran accept if not retry
       call ftst(r2)
       DO i=1,ndata
        call r8_normal_ab( xguess(i) ,beta, int(r2*100), ch2(i))

       ENDDO


      IF (qjump.lt.lim2) THEN  !0.01

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 331

      ELSE

       CONTINUE

      ENDIF 

       latin=lttm
       lonin=lntm
       

      ENDIF 

!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------
      IF (update_test.eq.5) THEN
      
      ssintm(1)=ssin((rot_test-1)*N+1)
      itst=0
      sstm=ssin

 441  call ftst(r)  
      ss(1)=int(r*60.)+1  
      IF (ss(1).lt.1.)   ss(1)=1.
      IF (ss(1).gt.61.)  ss(1)=61.



!add an extra control to make sure the new random point is not *exactly* the same as the initial one
 
       IF (ss(1).eq.ssintm(1)) GO TO 441

       sstm((rot_test-1)*N+1)=ss(1)

       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lonin,
     .         latin,sstm,fcin,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 441

!c calculate the transition probability to the new state:

       beta=beta_s 

       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))

c call a random number-> if qjump>ran accept if not retry
       call ftst(r2)
       DO i=1,ndata
        call r8_normal_ab(xguess(i) ,beta, int(r2*100), ch2(i))
       ENDDO


      IF (qjump.lt.lim3) THEN   !0.11

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 441

      ELSE

       CONTINUE

      ENDIF 

      ssin=sstm


      ENDIF


!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------
      IF (update_test.eq.6) THEN
 
 551     fcintm(1)=fcin((rot_test-1)*N+1)
         fctm=fcin
!          lat(1)=-28.
     
       call ftst(r) 
!       fc(1)=(int(r*maxfc)+1)/100.+0.3



!       fc(1)=(int(r*(maxfc+1-minfc))+minfc)/100.
       fc(1)=(int(r*(240+1-60))+60)/100.
       IF (fc(1).eq.1.) fc(1)=1.01
!         print*,lat(1),fc(1)

!         stop


      IF (ANY(fc.gt.240.)) print*,'error fc(1):gmhon.f',fc
      IF (ANY(fc.gt.240.)) STOP


!add an extra control to make sure the new random point is not *exactly* the same as the initial one

       IF (fc(1).eq.fcintm(1)) GO TO 551
       fctm((rot_test-1)*N+1)=fc(1)

       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lonin,
     .         latin,ssin,fctm,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 551

!c calculate the transition probability to the new state:

       beta=beta_f

       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))
   

      IF (qjump.lt.lim4) THEN     !0.1

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 551

      ELSE

       CONTINUE

      ENDIF 

      fcin=fctm

      ENDIF


!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------

      IF (update_test.eq.7) THEN

      rb=0.
       
 332  lonintm(2)=lonin((rot_test-1)*N+2)
      latintm(2)=latin((rot_test-1)*N+2)

      lttm=latin
      lntm=lonin

      call ftst(r)               
      lat(2)=-60.+(int(r*74.)+1)*2.   
 !     lat(2)=(int(r*45.)+1)*2.!-90.+(int(r*90.)+1)*2. !(int(r*45.)+1)*2.


!      lat(2)=latin(2)

      call ftst(r)  
      lon(2)=int(r*180.)*2. 


!        lon(2)=lonin(2)




!add an extra control to make sure the new random point is not *exactly* the same as the initial one
!
       IF ((lat(2).eq.latintm(2)).AND.(lon(2).eq.lonintm(2))) GO TO 332

       lttm((rot_test-1)*N+2)=lat(2)
       lntm((rot_test-1)*N+2)=lon(2)


       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lntm,lttm,
     .               ssin,fcin,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 332

!c calculate the transition probability to the new state:

       beta=beta_ll


       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))

!       qjump=qjump*(1./sqrt(2.*pi*(10.**2.)))
!       qjump=qjump*exp(-(lat(2)**2.)/(2.*10.**2.))
       
!c call a random number-> if qjump>ran accept if not retry
       call ftst(r2)
       DO i=1,ndata
        call r8_normal_ab( xguess(i) ,beta, int(r2*100), ch2(i))
       ENDDO


      IF (qjump.lt.lim2) THEN  !0.01

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 332

      ELSE

       CONTINUE

      ENDIF 

       latin=lttm
       lonin=lntm

      ENDIF 

!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------
      IF (update_test.eq.8) THEN
      
      ssintm(2)=ssin((rot_test-1)*N+2)
      itst=0
      sstm=ssin

 442  call ftst(r)  
      ss(2)=int(r*60.)+1  
      IF (ss(2).lt.1.)   ss(2)=1.
      IF (ss(2).gt.61.)  ss(2)=61.

!add an extra control to make sure the new random point is not *exactly* the same as the initial one

       IF (ss(2).eq.ssintm(2)) GO TO 442

       sstm((rot_test-1)*N+2)=ss(2)

       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lonin,
     .         latin,sstm,fcin,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 442

!c calculate the transition probability to the new state:

       beta=beta_s

       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))

c call a random number-> if qjump>ran accept if not retry
       call ftst(r2)
       DO i=1,ndata
        call r8_normal_ab(xguess(i) ,beta, int(r2*100), ch2(i))
       ENDDO


      IF (qjump.lt.lim3) THEN   !0.11

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 442

      ELSE

       CONTINUE

      ENDIF 

      ssin=sstm

      ENDIF


!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------
      IF (update_test.eq.9) THEN
 
 552     fcintm(2)=fcin((rot_test-1)*N+2)
         fctm=fcin

       call ftst(r) 
!       fc(2)=(int(r*maxfc)+1)/100.
!       fc(2)=(int(r*(maxfc+1-minfc))+minfc)/100.
       fc(2)=(int(r*(240+1-60))+60)/100.
       IF (fc(2).eq.1.) fc(2)=1.01


      IF (ANY(fc.gt.240.)) print*,'error fc(2):gmhon.f',fc
      IF (ANY(fc.gt.240.)) STOP


!add an extra control to make sure the new random point is not *exactly* the same as the initial one

       IF (fc(2).eq.fcintm(2)) GO TO 552

       fctm((rot_test-1)*N+2)=fc(2)

       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lonin,
     .         latin,ssin,fctm,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 552

!c calculate the transition probability to the new state:

       beta=beta_f

       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))
   

      IF (qjump.lt.lim4) THEN     !0.1

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 552

      ELSE

       CONTINUE

      ENDIF 

      fcin=fctm

      ENDIF

!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------


      IF (update_test.eq.10) THEN

      rb=0.
       
 333  lonintm(3)=lonin((rot_test-1)*N+3)
      latintm(3)=latin((rot_test-1)*N+3)

      lttm=latin
      lntm=lonin

      call ftst(r)               
      lat(3)=-60.+(int(r*74.)+1)*2.   !-90.+(int(r*90.)+1)*2. !
!      lat(3)=(int(r*45.)+1)*2. !(int(r*10.)+1)*9. !45.)+1)*2. 

!       lat(3)=latin(3)

      call ftst(r)  
      lon(3)=int(r*180.)*2. !int(r*360.)


!add an extra control to make sure the new random point is not *exactly* the same as the initial one

       IF ((lat(3).eq.latintm(3)).AND.(lon(3).eq.lonintm(3))) GO TO 333


       lttm((rot_test-1)*N+3)=lat(3)
       lntm((rot_test-1)*N+3)=lon(3)


       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lntm,lttm,
     .               ssin,fcin,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 333

!c calculate the transition probability to the new state:

       beta=beta_ll


       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))
       
!c call a random number-> if qjump>ran accept if not retry
       call ftst(r2)
       DO i=1,ndata
        call r8_normal_ab( xguess(i) ,beta, int(r2*100), ch2(i))

       ENDDO


      IF (qjump.lt.lim2) THEN  !0.01

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 333

      ELSE

       CONTINUE

      ENDIF 

       latin=lttm
       lonin=lntm

      ENDIF 

!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------
      IF (update_test.eq.11) THEN
      
      ssintm(3)=ssin((rot_test-1)*N+3)
      itst=0
      sstm=ssin


 443  call ftst(r)  
      ss(3)=int(r*60.)+1  
      IF (ss(3).lt.1.)   ss(3)=1.
      IF (ss(3).gt.61.)  ss(3)=61.


!add an extra control to make sure the new random point is not *exactly* the same as the initial one

       IF (ss(3).eq.ssintm(3)) GO TO 443
 
       sstm((rot_test-1)*N+3)=ss(3)

       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lonin,
     .         latin,sstm,fcin,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 443

!c calculate the transition probability to the new state:

       beta=beta_s

       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))

c call a random number-> if qjump>ran accept if not retry
       call ftst(r2)
       DO i=1,ndata
        call r8_normal_ab(xguess(i) ,beta, int(r2*100), ch2(i))
       ENDDO


      IF (qjump.lt.lim3) THEN   !0.11

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 443

      ELSE

       CONTINUE

      ENDIF 

      ssin=sstm

      ENDIF


!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------
      IF (update_test.eq.12) THEN
 
 553     fcintm(3)=fcin((rot_test-1)*N+3)
         fctm=fcin


     
       call ftst(r) 
!       fc(3)=(int(r*maxfc)+1)/100.
  !     fc(3)=(int(r*(maxfc+1-minfc))+minfc)/100.
       fc(3)=(int(r*(200+1-60))+60)/100. 
       IF (fc(3).eq.1.) fc(3)=1.01

      IF (ANY(fc.gt.240.)) print*,'error fc(3):gmhon.f',fc
      IF (ANY(fc.gt.240.)) STOP


!add an extra control to make sure the new random point is not *exactly* the same as the initial one
 
       IF (fc(3).eq.fcintm(3)) GO TO 553

       fctm((rot_test-1)*N+3)=fc(3)

       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lonin,
     .         latin,ssin,fctm,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 553

!c calculate the transition probability to the new state:

       beta=beta_f

       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))
   

      IF (qjump.lt.lim4) THEN     !0.1

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 553

      ELSE

       CONTINUE

      ENDIF 

      fcin=fctm

      ENDIF

!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------

      IF (update_test.eq.13) THEN

      rb=0.
       
 334  lonintm(4)=lonin((rot_test-1)*N+4)
      latintm(4)=latin((rot_test-1)*N+4)

      lttm=latin
      lntm=lonin

      call ftst(r)               
      lat(4)=-60.+(int(r*74.)+1)*2.   !-90.+(int(r*90.)+1)*2. 
!(int(r*45.)+1)*2.       
!      lat(4)=(int(r*45.)+1)*2.


!       lat(4)=latin(4)

      call ftst(r)  
      lon(4)=int(r*180.)*2. 



!add an extra control to make sure the new random point is not *exactly* the same as the initial one
       IF ((lat(4).eq.latintm(4)).AND.(lon(4).eq.lonintm(4))) GO TO 334

       lttm((rot_test-1)*N+4)=lat(4)
       lntm((rot_test-1)*N+4)=lon(4)


       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lntm,lttm,
     .               ssin,fcin,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 334

!c calculate the transition probability to the new state:

       beta=beta_ll


       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))
       
!c call a random number-> if qjump>ran accept if not retry
       call ftst(r2)
       DO i=1,ndata
        call r8_normal_ab( xguess(i) ,beta, int(r2*100), ch2(i))

       ENDDO


      IF (qjump.lt.lim2) THEN  !0.01

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 334

      ELSE

       CONTINUE

      ENDIF 

       latin=lttm
       lonin=lntm

      ENDIF 

!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------
      IF (update_test.eq.14) THEN
      
      ssintm(4)=ssin((rot_test-1)*N+4)
      itst=0
      sstm=ssin


 444  call ftst(r)  
      ss(4)=int(r*60.)+1  
      IF (ss(4).lt.1.)   ss(4)=1.
      IF (ss(4).gt.61.)  ss(4)=61.

!add an extra control to make sure the new random point is not *exactly* the same as the initial one

       IF (ss(4).eq.ssintm(4)) GO TO 444
 
       sstm((rot_test-1)*N+4)=ss(4)

       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lonin,
     .         latin,sstm,fcin,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 444

!c calculate the transition probability to the new state:

       beta=beta_s

       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))

!c call a random number-> if qjump>ran accept if not retry
       call ftst(r2)
       DO i=1,ndata
        call r8_normal_ab(xguess(i) ,beta, int(r2*100), ch2(i))
       ENDDO


      IF (qjump.lt.lim3) THEN   !0.11

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 444

      ELSE

       CONTINUE

      ENDIF 

       ssin=sstm

      ENDIF


!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------
      IF (update_test.eq.15) THEN
 
 554    fcintm(4)=fcin((rot_test-1)*N+4)
         fctm=fcin

    
       call ftst(r) 
!       fc(4)=(int(r*maxfc)+1)/100.
       fc(4)=(int(r*(240+1-60))+60)/100.
       IF (fc(4).eq.1.) fc(4)=1.01

      IF (ANY(fc.gt.240.)) print*,'error fc(4):gmhon.f',fc
      IF (ANY(fc.gt.240.)) STOP


!add an extra control to make sure the new random point is not *exactly* the same as the initial one

       IF (fc(4).eq.fcintm(4)) GO TO 554

       fctm((rot_test-1)*N+4)=fc(4)
       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lonin,
     .         latin,ssin,fctm,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 554

!c calculate the transition probability to the new state:

       beta=beta_f

       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))
   

      IF (qjump.lt.lim4) THEN     !0.1

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 554

      ELSE

       CONTINUE

      ENDIF 

      fcin=fctm

      ENDIF

!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------


      IF (update_test.eq.16) THEN

      rb=0.
       
 335  lonintm(5)=lonin((rot_test-1)*N+5)
      latintm(5)=latin((rot_test-1)*N+5)

      lttm=latin
      lntm=lonin

      call ftst(r)               
      lat(5)=-60.+(int(r*74.)+1)*2.   !-90.+(int(r*90.)+1)*2. 

!       lat(5)=(int(r*45.)+1)*2. 

      call ftst(r)  
      lon(5)=int(r*180.)*2. 



!add an extra control to make sure the new random point is not *exactly* the same as the initial one
       IF ((lat(5).eq.latintm(5)).AND.(lon(5).eq.lonintm(5))) GO TO 335

       lttm((rot_test-1)*N+5)=lat(5)
       lntm((rot_test-1)*N+5)=lon(5)

       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lntm,lttm,
     .               ssin,fcin,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 335

!c calculate the transition probability to the new state:

       beta=beta_ll


       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))
       
!c call a random number-> if qjump>ran accept if not retry
       call ftst(r2)
       DO i=1,ndata
        call r8_normal_ab( xguess(i) ,beta, int(r2*100), ch2(i))

       ENDDO


      IF (qjump.lt.lim2) THEN  !0.01

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 335

      ELSE

       CONTINUE

      ENDIF 

       latin=lttm
       lonin=lntm

      ENDIF 

!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------
      IF (update_test.eq.17) THEN
      
      ssintm(5)=ssin((rot_test-1)*N+5)
      itst=0
      sstm=ssin


 445  call ftst(r)  
      ss(5)=int(r*60.)+1  
      IF (ss(5).lt.1.)   ss(5)=1.
      IF (ss(5).gt.61.)  ss(5)=61.

!add an extra control to make sure the new random point is not *exactly* the same as the initial one

       IF (ss(5).eq.ssintm(5)) GO TO 445
 
       sstm((rot_test-1)*N+5)=ss(5)

       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lonin,
     .         latin,sstm,fcin,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 445

!c calculate the transition probability to the new state:

       beta=beta_s

       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))

c call a random number-> if qjump>ran accept if not retry
       call ftst(r2)
       DO i=1,ndata
        call r8_normal_ab(xguess(i) ,beta, int(r2*100), ch2(i))
       ENDDO


      IF (qjump.lt.lim3) THEN   !0.11

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 445

      ELSE

       CONTINUE

      ENDIF 

       ssin=sstm

      ENDIF


!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------
      IF (update_test.eq.18) THEN
 
 555     fcintm(5)=fcin((rot_test-1)*N+5)
         fctm=fcin

    
       call ftst(r) 
!       fc(5)=(int(r*maxfc)+1)/100.
       fc(5)=(int(r*(240+1-60))+60)/100.
       IF (fc(5).eq.1.) fc(5)=1.01

      IF (ANY(fc.gt.240.)) print*,'error fc(5):gmhon.f',fc
      IF (ANY(fc.gt.240.)) STOP


!add an extra control to make sure the new random point is not *exactly* the same as the initial one

       IF (fc(5).eq.fcintm(5)) GO TO 555

       fctm((rot_test-1)*N+5)=fc(5)

       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lonin,
     .         latin,ssin,fctm,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 555

!c calculate the transition probability to the new state:

       beta=beta_f

       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))
   

      IF (qjump.lt.lim4) THEN     !0.1

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 555

      ELSE

       CONTINUE

      ENDIF 

      fcin=fctm

      ENDIF

!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------


      IF (update_test.eq.19) THEN

      rb=0.
       
 336  lonintm(6)=lonin((rot_test-1)*N+6)
      latintm(6)=latin((rot_test-1)*N+6)

      lttm=latin
      lntm=lonin

      call ftst(r)               
      lat(6)=-60.+(int(r*74.)+1)*2.   !-90.+(int(r*90.)+1)*2. 
!(int(r*45.)+1)*2. 
!       lat(6)=(int(r*45.)+1)*2.

      call ftst(r)  
      lon(6)=int(r*180.)*2. 



!add an extra control to make sure the new random point is not *exactly* the same as the initial one
       IF ((lat(6).eq.latintm(6)).AND.(lon(6).eq.lonintm(6))) GO TO 336

       lttm((rot_test-1)*N+6)=lat(6)
       lntm((rot_test-1)*N+6)=lon(6)


       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lntm,lttm,
     .               ssin,fcin,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 336

!c calculate the transition probability to the new state:

       beta=beta_ll


       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))
       
!c call a random number-> if qjump>ran accept if not retry
       call ftst(r2)
       DO i=1,ndata
        call r8_normal_ab( xguess(i) ,beta, int(r2*100), ch2(i))

       ENDDO


      IF (qjump.lt.lim2) THEN  !0.01

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 336

      ELSE

       CONTINUE

      ENDIF 

       latin=lttm
       lonin=lntm

      ENDIF 

!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------
      IF (update_test.eq.20) THEN
      
      ssintm(6)=ssin((rot_test-1)*N+6)
      itst=0
      sstm=ssin


 446  call ftst(r)  
      ss(6)=int(r*60.)+1  
      IF (ss(6).lt.1.)   ss(6)=1.
      IF (ss(6).gt.61.)  ss(6)=61.

!add an extra control to make sure the new random point is not *exactly* the same as the initial one

       IF (ss(6).eq.ssintm(6)) GO TO 446
 
       sstm((rot_test-1)*N+6)=ss(6)

       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lonin,
     .         latin,sstm,fcin,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 446

!c calculate the transition probability to the new state:

       beta=beta_s

       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))

c call a random number-> if qjump>ran accept if not retry
       call ftst(r2)
       DO i=1,ndata
        call r8_normal_ab(xguess(i) ,beta, int(r2*100), ch2(i))
       ENDDO


      IF (qjump.lt.lim3) THEN   !0.11

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 446

      ELSE

       CONTINUE

      ENDIF 

      ssin=sstm

      ENDIF


!c----------------------------------------------------------------------------
!c----------------------------------------------------------------------------
      IF (update_test.eq.21) THEN
 
 556    fcintm(6)=fcin((rot_test-1)*N+6)
         fctm=fcin


       call ftst(r) 
!       fc(6)=(int(r*maxfc)+1)/100.
       fc(6)=(int(r*(240+1-60))+60)/100.
       IF (fc(6).eq.1.) fc(6)=1.01

      IF (ANY(fc.gt.240.)) print*,'error fc(6):gmhon.f',fc
      IF (ANY(fc.gt.240.)) STOP


!add an extra control to make sure the new random point is not *exactly* the same as the initial one

       IF (fc(6).eq.fcintm(6)) GO TO 556

       fctm((rot_test-1)*N+6)=fc(6)

       call BDsignal(cntrl,incin,tin,ldin,per,nspots,lonin,
     .         latin,ssin,fctm,Fmod,tmod)


       IF (cntrl.eq.1) GO TO 556

!c calculate the transition probability to the new state:

       beta=beta_f

       qb=(1./sqrt(2.*pi*(beta**2.)))

       qjump=qb*exp(SUM(-((Fmod-xguess)**2./(2*beta**2.))))
   

      IF (qjump.lt.lim4) THEN     !0.1

       itst=itst+1

       IF (itst.eq.itstlim) GO TO 38

       GO TO 556

      ELSE

       CONTINUE

      ENDIF 

      fcin=fctm

      ENDIF



!c------------------------------------------------------
!c----------------------------------------------------------------------------


c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      IF (PCA.eq.1) THEN 

      IF (update_test.eq.3) fcin((rot_test-1)*N+2:(rot_test-1)*N+
     .                  nspots)=fcin((rot_test-1)*N+1) !1.3
      IF (update_test.eq.6) fcin((rot_test-1)*N+2:(rot_test-1)*N+
     .                  nspots)=fcin((rot_test-1)*N+1)
      IF (update_test.eq.9) fcin((rot_test-1)*N+1)=
     .                  fcin((rot_test-1)*N+2)
      IF (update_test.eq.9) fcin((rot_test-1)*N+3:(rot_test-1)*N+
     .                  nspots)=fcin((rot_test-1)*N+2)
      IF (update_test.eq.12) fcin((rot_test-1)*N+1:(rot_test-1)*N+2)=
     .                  fcin((rot_test-1)*N+3)
      IF (update_test.eq.12) fcin((rot_test-1)*N+4:(rot_test-1)*N+
     .                  nspots)=fcin((rot_test-1)*N+3)
      IF (update_test.eq.15) fcin((rot_test-1)*N+1:(rot_test-1)*N+3)=
     .                  fcin((rot_test-1)*N+4)
      IF (update_test.eq.15) fcin((rot_test-1)*N+5:(rot_test-1)*N+
     .                  nspots)=fcin((rot_test-1)*N+4)
      IF (update_test.eq.18) fcin((rot_test-1)*N+1:(rot_test-1)*N+4)=
     .                  fcin((rot_test-1)*N+5)
      IF (update_test.eq.18) fcin((rot_test-1)*N+6)=
     .                  fcin((rot_test-1)*N+5)
      IF (update_test.eq.21) fcin((rot_test-1)*N+1:(rot_test-1)*N+5)=
     .                  fcin((rot_test-1)*N+6)
 
      ENDIF

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!       incin=64.


 !     print*,'#sp:',cnt,nspots!,update_test
!          stop

      RETURN



! 999   write(*,*) 'error reading grid file!'
! 10    format(f10.4,2x,f12.8)


       END SUBROUTINE gmh_on

