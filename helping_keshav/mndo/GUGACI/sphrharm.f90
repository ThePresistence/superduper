Module sphrharm

Public

Double Precision, Parameter :: Fac1 = 0.86602540378443864676D0

Type CartComp
  Integer          :: k, j
  Double Precision :: c
End Type CartComp

Type RealSphrHarm
  Integer                      :: n
  Type(CartComp), Dimension(3) :: cart
End Type RealSphrHarm

Type(RealSphrHarm), Dimension(0:2, 5) :: Ylm


Contains


Subroutine InitYlm

  Ylm(0,1)%n       = 1                       ! s
  Ylm(0,1)%cart(1) = CartComp(0, 0,  1.0D0)

  Ylm(1,1)%n       = 1                       ! px
  Ylm(1,1)%cart(1) = CartComp(0, 0,  1.0D0)

  Ylm(1,2)%n       = 1                       ! py
  Ylm(1,2)%cart(1) = CartComp(0, 1,  1.0D0)

  Ylm(1,3)%n       = 1                       ! pz
  Ylm(1,3)%cart(1) = CartComp(1, 0,  1.0D0)

  Ylm(2,1)%n       = 2                       ! d(+2) = dx2-y2
  Ylm(2,1)%cart(1) = CartComp(0, 0,  Fac1) 
  Ylm(2,1)%cart(2) = CartComp(0, 2, -Fac1) 

  Ylm(2,2)%n       = 1                       ! d(+1) = dxz
  Ylm(2,2)%cart(1) = CartComp(1, 0,  1.0D0)

  Ylm(2,3)%n       = 3                       ! d(0)  = dz2
  Ylm(2,3)%cart(1) = CartComp(0, 0, -0.5D0)
  Ylm(2,3)%cart(2) = CartComp(0, 2, -0.5D0)
  Ylm(2,3)%cart(3) = CartComp(2, 0,  1.0D0)

  Ylm(2,4)%n       = 1                       ! d(-1) = dyz
  Ylm(2,4)%cart(1) = CartComp(1, 1,  1.0D0)

  Ylm(2,5)%n       = 1                       ! d(-2) = dxy
  Ylm(2,5)%cart(1) = CartComp(0, 1,  1.0D0)

End Subroutine InitYlm


End Module sphrharm
