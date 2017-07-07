      SUBROUTINE random_flat(n,r)

      INTEGER i, n
      DOUBLE PRECISION seeda, r(n)

      call random_seed()
      DO i=1,n
        call random_number(seeda)
        r(i)=seeda
      END DO


      RETURN
      END 
