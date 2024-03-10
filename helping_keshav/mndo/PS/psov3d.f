C     ******************************************************************
C
C     NMR 3-center integrals over Slater AOs: Low-level routines.
C
C     ******************************************************************
C
C     Definition of the expansion coefficients for Sharma's expansion 
C     at the displaced center (Phys. Rev. A, 13, 517 (1976)).
C
      BLOCK DATA PS3SHS
C
C   Define values of the expansion coefficients for PS3SHR.
C
C   Coded by: Serge Pachkovsky (sort of)
C
C   NOTE(WT): Coefficients are now computed as needed.
C             Preferred by PS (10/3/98).
C
C   Assisted by: Mathematica
C
C   Defined common blocks:
C
C      PS3SHD - Values of the expansion coefficients.
C         COEFF - packed 5D array of precomputed expansion
C                 coefficients. First two indices (Expansion
C                 term order LS and orbital quantum number L)
C                 are packed using the IBASE array. Remaining
C                 three (NU, S and magnetic quantum number M)
C                 are represented by rectangular array. Valid
C                 range of the indices is:
C                     NU: 0 to L+LS-S (although space is reserved
C                                      for L+LS+1 elements)
C                     S:  0 to L+LS
C                     M:  0 to L
C                     L:  0 to MAXL
C                     LS: 0 to MAXLS
C         IBASE - index into COEFF for given LS,L pair
C
C   Local storage:
C
C   Module logic:
C
C   Possible optimizations:
C
C   Bugs:
C
C      Approximately half of the space occupied by COEFF is wasted.
C      However, it makes life so much simpler...
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=3)
      PARAMETER (MAXLS=40)
C     PARAMETER (LCOEFF=274700)
C
      COMMON 
     ./PS3SHE/ IBASE(0:MAXLS,0:MAXL)
C    ./PS3SHD/ COEFF(LCOEFF)
C     SAVE /PS3SHD/
      SAVE /PS3SHE/
C    Coefficients which have space reserved for 'em
      DATA IBASE(0,0) /-1/
      DATA IBASE(1,0) /-2/
      DATA IBASE(2,0) /-6/
      DATA IBASE(3,0) /-15/
      DATA IBASE(4,0) /-31/
      DATA IBASE(5,0) /-56/
      DATA IBASE(6,0) /-92/
      DATA IBASE(7,0) /-141/
      DATA IBASE(8,0) /-205/
      DATA IBASE(9,0) /-286/
      DATA IBASE(10,0) /-386/
      DATA IBASE(11,0) /-507/
      DATA IBASE(12,0) /-651/
      DATA IBASE(13,0) /-820/
      DATA IBASE(14,0) /-1016/
      DATA IBASE(15,0) /-1241/
      DATA IBASE(16,0) /-1497/
      DATA IBASE(17,0) /-1786/
      DATA IBASE(18,0) /-2110/
      DATA IBASE(19,0) /-2471/
      DATA IBASE(20,0) /-2871/
      DATA IBASE(21,0) /-3312/
      DATA IBASE(22,0) /-3796/
      DATA IBASE(23,0) /-4325/
      DATA IBASE(24,0) /-4901/
      DATA IBASE(25,0) /-5526/
      DATA IBASE(26,0) /-6202/
      DATA IBASE(27,0) /-6931/
      DATA IBASE(28,0) /-7715/
      DATA IBASE(29,0) /-8556/
      DATA IBASE(30,0) /-9456/
      DATA IBASE(31,0) /-10417/
      DATA IBASE(32,0) /-11441/
      DATA IBASE(33,0) /-12530/
      DATA IBASE(34,0) /-13686/
      DATA IBASE(35,0) /-14911/
      DATA IBASE(36,0) /-16207/
      DATA IBASE(37,0) /-17576/
      DATA IBASE(38,0) /-19020/
      DATA IBASE(39,0) /-20541/
      DATA IBASE(40,0) /-22141/
      DATA IBASE(0,1) /-23822/
      DATA IBASE(1,1) /-23830/
      DATA IBASE(2,1) /-23848/
      DATA IBASE(3,1) /-23880/
      DATA IBASE(4,1) /-23930/
      DATA IBASE(5,1) /-24002/
      DATA IBASE(6,1) /-24100/
      DATA IBASE(7,1) /-24228/
      DATA IBASE(8,1) /-24390/
      DATA IBASE(9,1) /-24590/
      DATA IBASE(10,1) /-24832/
      DATA IBASE(11,1) /-25120/
      DATA IBASE(12,1) /-25458/
      DATA IBASE(13,1) /-25850/
      DATA IBASE(14,1) /-26300/
      DATA IBASE(15,1) /-26812/
      DATA IBASE(16,1) /-27390/
      DATA IBASE(17,1) /-28038/
      DATA IBASE(18,1) /-28760/
      DATA IBASE(19,1) /-29560/
      DATA IBASE(20,1) /-30442/
      DATA IBASE(21,1) /-31410/
      DATA IBASE(22,1) /-32468/
      DATA IBASE(23,1) /-33620/
      DATA IBASE(24,1) /-34870/
      DATA IBASE(25,1) /-36222/
      DATA IBASE(26,1) /-37680/
      DATA IBASE(27,1) /-39248/
      DATA IBASE(28,1) /-40930/
      DATA IBASE(29,1) /-42730/
      DATA IBASE(30,1) /-44652/
      DATA IBASE(31,1) /-46700/
      DATA IBASE(32,1) /-48878/
      DATA IBASE(33,1) /-51190/
      DATA IBASE(34,1) /-53640/
      DATA IBASE(35,1) /-56232/
      DATA IBASE(36,1) /-58970/
      DATA IBASE(37,1) /-61858/
      DATA IBASE(38,1) /-64900/
      DATA IBASE(39,1) /-68100/
      DATA IBASE(40,1) /-71462/
      DATA IBASE(0,2) /-74990/
      DATA IBASE(1,2) /-75017/
      DATA IBASE(2,2) /-75065/
      DATA IBASE(3,2) /-75140/
      DATA IBASE(4,2) /-75248/
      DATA IBASE(5,2) /-75395/
      DATA IBASE(6,2) /-75587/
      DATA IBASE(7,2) /-75830/
      DATA IBASE(8,2) /-76130/
      DATA IBASE(9,2) /-76493/
      DATA IBASE(10,2) /-76925/
      DATA IBASE(11,2) /-77432/
      DATA IBASE(12,2) /-78020/
      DATA IBASE(13,2) /-78695/
      DATA IBASE(14,2) /-79463/
      DATA IBASE(15,2) /-80330/
      DATA IBASE(16,2) /-81302/
      DATA IBASE(17,2) /-82385/
      DATA IBASE(18,2) /-83585/
      DATA IBASE(19,2) /-84908/
      DATA IBASE(20,2) /-86360/
      DATA IBASE(21,2) /-87947/
      DATA IBASE(22,2) /-89675/
      DATA IBASE(23,2) /-91550/
      DATA IBASE(24,2) /-93578/
      DATA IBASE(25,2) /-95765/
      DATA IBASE(26,2) /-98117/
      DATA IBASE(27,2) /-100640/
      DATA IBASE(28,2) /-103340/
      DATA IBASE(29,2) /-106223/
      DATA IBASE(30,2) /-109295/
      DATA IBASE(31,2) /-112562/
      DATA IBASE(32,2) /-116030/
      DATA IBASE(33,2) /-119705/
      DATA IBASE(34,2) /-123593/
      DATA IBASE(35,2) /-127700/
      DATA IBASE(36,2) /-132032/
      DATA IBASE(37,2) /-136595/
      DATA IBASE(38,2) /-141395/
      DATA IBASE(39,2) /-146438/
      DATA IBASE(40,2) /-151730/
      DATA IBASE(0,3) /-157277/
      DATA IBASE(1,3) /-157341/
      DATA IBASE(2,3) /-157441/
      DATA IBASE(3,3) /-157585/
      DATA IBASE(4,3) /-157781/
      DATA IBASE(5,3) /-158037/
      DATA IBASE(6,3) /-158361/
      DATA IBASE(7,3) /-158761/
      DATA IBASE(8,3) /-159245/
      DATA IBASE(9,3) /-159821/
      DATA IBASE(10,3) /-160497/
      DATA IBASE(11,3) /-161281/
      DATA IBASE(12,3) /-162181/
      DATA IBASE(13,3) /-163205/
      DATA IBASE(14,3) /-164361/
      DATA IBASE(15,3) /-165657/
      DATA IBASE(16,3) /-167101/
      DATA IBASE(17,3) /-168701/
      DATA IBASE(18,3) /-170465/
      DATA IBASE(19,3) /-172401/
      DATA IBASE(20,3) /-174517/
      DATA IBASE(21,3) /-176821/
      DATA IBASE(22,3) /-179321/
      DATA IBASE(23,3) /-182025/
      DATA IBASE(24,3) /-184941/
      DATA IBASE(25,3) /-188077/
      DATA IBASE(26,3) /-191441/
      DATA IBASE(27,3) /-195041/
      DATA IBASE(28,3) /-198885/
      DATA IBASE(29,3) /-202981/
      DATA IBASE(30,3) /-207337/
      DATA IBASE(31,3) /-211961/
      DATA IBASE(32,3) /-216861/
      DATA IBASE(33,3) /-222045/
      DATA IBASE(34,3) /-227521/
      DATA IBASE(35,3) /-233297/
      DATA IBASE(36,3) /-239381/
      DATA IBASE(37,3) /-245781/
      DATA IBASE(38,3) /-252505/
      DATA IBASE(39,3) /-259561/
      DATA IBASE(40,3) /-266957/
C    Precomputed coefficients (do we really need 'em?)
      END
