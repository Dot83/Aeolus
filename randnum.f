c Create the random number in the 0,1 space to control the validity of the step 

       SUBROUTINE randnum(k) 

       implicit none

       DOUBLE PRECISION r(5,5)
       DOUBLE PRECISION k

       call init_random_seed() 
       call random_number(r)

       k=r(1,1)

       RETURN

       END
