c Subroutine to create the target's light curve


       subroutine BDsignal(cntrl,inc,ld, per,fcz,nspots,ilon,ilat,iss,ifc,     
     .                     Fmodout,tmodout)


       IMPLICIT NONE

       INCLUDE 'const_incl' 

       INTEGER nspots,cntrl,imp
       DOUBLE PRECISION inc,per,ld,el_par,sz,tsm

!we initially define the amount of pixels the planet shall be split in:
       INTEGER, PARAMETER :: nlon=361
       INTEGER, PARAMETER :: nlat=181, data=1000

       
       INTEGER :: i, m,k, alpha,j,chk,corr_fc,inc_shft,nspotsm(10)
!     note that exceptions are raised when a spot is at 0deg/360deg longitude and to be plotted
!     correctly the spot occupies 2 slots (for the 2 edges of the map; 0/360)
       REAL(8), DIMENSION(N*20) :: lat,ss,fc,lon,lat_inv
       REAL(8), DIMENSION(N*9) ::ilat,ilon,iss,ifc
       REAL(8), PARAMETER :: degs = 0.0174533D0, pi=3.14159D0
       REAL(8), DIMENSION(3) :: fcz,phs
       
       REAL(8), DIMENSION(5000) :: Fmodout,tmodout
       REAL(8), DIMENSION(nlat,nlon) :: mp,mp0
       REAL(8), DIMENSION (nlat,nlat) :: mp_int, theta,mp_tmp1,mp_tmp2
      
       REAL(8), DIMENSION(nlat,nlon*10) :: mp_tot,sin1,sin2,sin3
       REAL(8), DIMENSION(nlat) :: lats
       REAL(8), DIMENSION(nlon) ::lons
       REAL(8), DIMENSION (701) :: lc
       REAL(8), DIMENSION(nlat,nlon) :: lats2,lons2
       REAL(8), DIMENSION(nlon*10) :: sin1t1,sin1t2,sin1t3
 
       
       el_par=1.
       phs(1)=0.*degs
       phs(2)=0.*degs
       phs(3)=0.*degs
       ilat(1)=lat(1)
       ilat(2)=lat(13)
       ilat(3)=lat(25)
       ilat(4)=lat(37)
       ilat(5)=lat(49)
       ilat(6)=lat(61)
       ilat(7)=lat(73)
       ilat(8)=lat(91)
       ilat(9)=lat(97)
       ilon(1)=lon(1)
       ilon(2)=lon(13)
       ilon(3)=lon(25)
       ilon(4)=lon(37)
       ilon(5)=lon(49)
       ilon(6)=lon(61)
       ilon(7)=lon(73)
       ilon(8)=lon(91)
       ilon(9)=lon(97)
       iss(1)=ss(1)
       iss(2)=ss(13)
       iss(3)=ss(25)
       iss(4)=ss(37)
       iss(5)=ss(49)
       iss(6)=ss(61)
       iss(7)=ss(73)
       iss(8)=ss(91)
       iss(9)=ss(97)
       ifc(1)=fc(1)
       ifc(2)=fc(13)
       ifc(3)=fc(25)
       ifc(4)=fc(37)
       ifc(5)=fc(49)
       ifc(6)=fc(61)
       ifc(7)=fc(73)
       ifc(8)=fc(91)
       ifc(9)=fc(97)







       DO i=1,nspots
          IF ((lon(i)-ss(i).lt.0).or.(lon(i)+ss(i).gt.360)) THEN
             lon(i+6)=360.
             lat(i+6)=lat(i)
          IF (lon(i)-ss(i).lt.0) THEN
             ss(i+6)=lon(i)-ss(i)
          ENDIF
          IF (lon(i)+ss(i).gt.360) THEN
             ss(i+6)=lon(i)-360.
          ENDIF
             fc(i+6)=fc(i)
          ENDIF   
       ENDDO   
       
       
!set all necessary parameters to their zero-state values
      lc=0.
      lons(1)=0.
      lats(1)=-90.
      mp=3.54410D-8
      corr_fc=0.
      mp_tot=0.
      mp_tmp1=0.
      mp_tmp2=0.
      fmodout=0.
      tmodout=0.
      
!read-in the file with projection angles theta (interpolated version!)
       OPEN(unit=2,file='000.thetaB',err=100) 
       DO imp=1,181,1 
          READ(2,*,err=100) theta(imp,1:181) 
       END DO
       CLOSE(2)

       theta=transpose(theta)
       
!make the longitude-latitude arrays
      DO i=2,nlon
         lons(i)=lons(i-1)+1
      ENDDO
      DO i=2,nlat
         lats(i)=lats(i-1)+1
      ENDDO
      DO i=1,nlon
        DO j=1,nlat
           lons2(j,i)=lons(i)
           lats2(j,i)=lats(j)
        ENDDO
      ENDDO   



      DO i=1,nlon*10 
       sin1t1(i)=sin(2*pi*i/((12.57/per)*360.)+phs(1))*fcz(1)
       sin1t2(i)=sin(2*pi*i/((14.22/per)*360.)+phs(2))*fcz(2)
       sin1t3(i)=sin(2*pi*i/((6.2/per)*360.)+phs(3))*fcz(3)
      ENDDO
      
       sin1t1=sin1t1+abs(minval(sin1t1))
       sin1t2=sin1t2+abs(minval(sin1t2))
       sin1t3=sin1t3+abs(minval(sin1t3))

       tsm=sum(mp)/nlon/nlat


       DO i=1,nlon*10
             sin1(110:130,i)=sin1t1(i)*0.1+tsm
             sin2(50:70,i)=sin1t2(i)*0.1+tsm
             sin3(80:100,i)=sin1t3(i)*0.1+tsm
       ENDDO

       
!     make map

       mp0=mp

 !##########################################################################      
      DO m=1,10

         mp=mp0
         
      DO k=1+12*(m-1),nspotsm(m)+12*(m-1)
         
         WHERE ((((lons2-lon(k))/ss(k))**2.+ ((lats2-lat(k))/
     .    (el_par*ss(k)))**2./(cos(lats2*degs))).lt.1.)

         mp=fc(k)*mp

        END WHERE
         
      ENDDO   

      
      
      mp=cshift(mp,180,2) 

!     shift map for the inclination of the target

      DO i=1,nlat
         IF (lats(i).eq.inc) THEN
            chk=i
            exit
         ENDIF
      ENDDO
       
!      calculate how many pixels the shift is 
      inc_shft=-91+chk


!     calculate the latitude of spots in the 'dark side'
      DO i=1,nspots
       lat_inv(i)=lat(i)+inc
      ENDDO
       
      DO i=1,nspots
         IF (lat_inv(i).ge.90.-inc) THEN
            corr_fc=1
            exit
         ENDIF
      ENDDO


 
        mp_tot(1:181,(m-1)*360+1:m*360)=mp(1:181,1:361)
      ENDDO


!##########################################################################
      
!   add belts        
      mp_tot(110:130,:)=sin1(110:130,:)
      mp_tot(50:70,:)=sin2(50:70,:)
      mp_tot(80:100,:)=sin1(80:100,:)

     
!     start the disk-integration process
      DO alpha=225,3900
         mp_int(1:181,1:181)=mp_tot(1:181,alpha*2:alpha*2+180)
 
         IF (corr_fc.eq.1) THEN
            
            mp_int=cshift(mp_int,inc_shft,1)
           
            mp_tmp1(1:inc_shft+1,1:181)=
     .             mp_tot(181-inc_shft:181,alpha*2-181+1:alpha*2+1)

           
          DO i=1,inc_shft
             mp_tmp2(i,1:181)=mp_tmp1(inc_shft+1-i,1:181)
          ENDDO

          DO i=1,181
             mp_tmp1(1:inc_shft,i)=mp_tmp2(1:inc_shft,182-i)
          ENDDO


          mp_int(181-inc_shft:181,1:181)=mp_tmp1(1:inc_shft,1:181)
          
         ENDIF

         
        mp_int=mp_int*(1-ld*(1-cos(theta*degs)))
        
        Fmodout(alpha-224)=sum(mp_int*cos(theta*degs))
        tmodout(alpha-224)=((alpha-225.)/180.)*per
        
        
      ENDDO


         fmodout(1:1500)=fmodout(1:1500)/(sum(fmodout(1:1500)/1500.))

         DO i=1,1500
               print*, tmodout(i), fmodout(i)
         ENDDO

         
       RETURN
      
*-------------------------------------------------------------------------------
*     Error messages:
*-------------------------------------------------------------------------------

 404   FORMAT(F23.17,2x,F23.17,2x,F23.17,2x,F23.17,2x,F23.17,2x,
     .      F23.17,2x,F23.17,2x,F23.17,2x,F23.17,2x,F23.17)
 405   FORMAT(4x)      
 100   WRITE(*,*) 'ERROR READING input file: thetas'
       STOP
       
 997   WRITE(*,*) 'ERROR OPENING FILE wrout.f'
      STOP
*-------------------------------------------------------------------------------


!      RETURN
      END
