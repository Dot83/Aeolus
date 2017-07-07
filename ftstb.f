c random number generator

      SUBROUTINE ftstb(k,n)
      REAL :: r(5,5) 
      INTEGER n
      DOUBLE PRECISION k(1:n),m(1:15)

      k=0.
      m=0.

      CALL init_random_seed()         ! see example of RANDOM_SEED
      CALL RANDOM_NUMBER(r)

      m=(/r(1,1),r(1,2),r(1,3),r(1,4),r(1,5),r(2,1),r(2,2),r(2,3),
     .   r(2,4),r(2,5),r(3,1),r(3,2),r(3,3),r(3,4),r(3,5) /)
      
      k(1:n)=m(1:n)


      RETURN

      end subroutine
