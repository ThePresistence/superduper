      integer function idamaxp(n,dx) 
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision dmax
      DOUBLE PRECISION, DIMENSION (:), POINTER :: dx
      integer i,ix,n
c
      idamaxp = 1
c
      dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamaxp = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
