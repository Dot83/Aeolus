c random number generator

      SUBROUTINE ftst(k)
      REAL :: r(5,5) 
      DOUBLE PRECISION k

      CALL init_random_seed()         ! see example of RANDOM_SEED
      CALL RANDOM_NUMBER(r)

      k=r(1,1)

      RETURN

      end subroutine
