Subroutine gpudavliu
  Use gugaglobal, only: Stdout
  Implicit None
  Write(Stdout,'(1X,"CUDA implementation of Davidson diagonalizer not available.")')
  Stop 'GPUDAVLIU'
End Subroutine gpudavliu
