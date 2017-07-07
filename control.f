c Subroutine to control whether the step is approved or not 

      SUBROUTINE control(xin,tin,sigmain,xguess,tguess,xn,tn,u,!xo,
     .                   test,chiout) !(xguess,tguess,xin,tin,sigmain,xn,tn,u,xo,test)

      implicit none
      INCLUDE 'const_incl'

      INTEGER test,i,nel,p,i1,i2,lim
      DOUBLE PRECISION u,cnt,chi2x,chi2xn,chiout
       
      REAL(8), DIMENSION(obsmax) :: xin,sigmain,tin
      REAL(8), DIMENSION(ndata) :: xguess,xn,xo,tguess,tn
      REAL(8), DIMENSION(ndata) :: tup,fup

      logical, dimension(obsmax) :: elmask
      logical, dimension(ndata) :: elmaskb

       chi2x=0.
       chi2xn=0.
       fup=0.
       tup=0.
       
c to compare observations with the model we need to match their phases 
c by bracketing the normalized period  



c find the number of elements in the periods array

!      elmask=tin.ne.0. !accept the elements where we have an observation (period measurement)

!      nel = COUNT(elmask) !cound how many elements that is
      
!      nel=nel+1  !to correct for the P=0 element that was not taken into account due to the mask


!      elmaskb=tguess.ne.0    !find how many elements the model has
!      lim=COUNT(elmaskb)

!bracket the observations at the corresponding model phases

!       DO i=1,lim+1
    
!       call  brack(tguess(i),tin,lim,ndata,i1,i2)

!       IF ((i2-i1).gt.1) THEN
!         p=int(i1+((i2-i1)/2.))
!       ELSE
!         p=i1
!       ENDIF

!         tup(i)=tin(p)         !fix the order of the observations to match those of the model
!         fup(i)=xin(p)


!      ENDDO


       DO i=1,N
        chi2x=chi2x+(xguess(i)-xin(i))**2./sigmain(i)**2.
       ENDDO

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ repeat for the other model
!      tup=0.
!      fup=0.
!      elmaskb=tn.ne.0
!      lim=COUNT(elmaskb)

!      DO i=1,lim+1!nel
       
!       call  brack(tn(i),tin,lim,ndata,i1,i2)!(tn(i),tguess,lim,ndata,i1,i2)
       
!       IF ((i2-i1).gt.1) THEN
!         p=int(i1+((i2-i1)/2.))
!       ELSE
!         p=i1
!       ENDIF

!         tup(i)=tin(p)         
!         fup(i)=xin(p)

!      ENDDO



       DO i=1,N
        chi2xn=chi2xn+(xn(i)-xin(i))**2./sigmain(i)**2.
       ENDDO

!compare models and check which one you accept

       cnt=MIN(exp((chi2x-chi2xn)/2.),1.)

!       IF (u.le.cnt) xo=xn
!       IF (u.gt.cnt) xo=xguess

       IF (u.le.cnt) test=1
       IF (u.gt.cnt) test=0
       
       IF (u.le.cnt) chiout=chi2xn
       IF (u.gt.cnt) chiout=chi2x
       


       RETURN

       END



!       DO i=1,N

!        chi2x=chi2x+(xguess(i)-xin(i))**2./sigmain(i)**2.
!        chi2xn=chi2xn+(xn(i)-xin(i))**2./sigmain(i)**2.

!       ENDDO
