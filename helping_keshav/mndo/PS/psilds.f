C     ******************************************************************
      SUBROUTINE PSILDS(R,ARO,BRO,AD,IMOD,W,W1,W2)
C
C  This subroutine is not intended to be human-readable.
C  It was generated by a Mathematica program (makeints.m).
C  General version of 14 May 1997.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARO(9), BRO(9), AD(5)
      DIMENSION W(461),W1(461),W2(461)
      PARAMETER (CN1=-0.5D0)
      PARAMETER (CN2=0.5D0)
      PARAMETER (CN3=2.D0)
      PARAMETER (CN4=0.25D0)
      PARAMETER (CN5=1.4142135623730950488016887242096980785697D0)
      PARAMETER (CN6=-1.4142135623730950488016887242096980785697D0)
      PARAMETER (CN7=1.1547005383792515290182975610039149112952D0)
      PARAMETER (CN8=-1.3333333333333333333333333333333333333333D0)
      PARAMETER (CN9=0.6666666666666666666666666666666666666667D0)
      PARAMETER (CN10=1.3333333333333333333333333333333333333333D0)
      PARAMETER (CN11=3.D0)
      PARAMETER (CN12=-1.5D0)
      PARAMETER (CN13=1.5D0)
      PARAMETER (CN14=0.75D0)
      PARAMETER (CN15=-0.25D0)
      T1 = R**2
      T42 = SQRT(T1 + (ARO(1) + BRO(1))**2)
      W(1) = 1/T42
      T44 = SQRT(T1 + (ARO(2) + BRO(1))**2)
      W(2) = 1/T44
      T3 = (ARO(3) + BRO(1))**2
      T2 = -R
      T46 = T2 - AD(1)
      T4 = T46**2
      T47 = SQRT(T3 + T4)
      T48 = T2 + AD(1)
      T5 = T48**2
      T49 = SQRT(T3 + T5)
      W(3) = CN1/T47 + CN2/T49
      T52 = SQRT(T1 + (ARO(5) + BRO(1))**2)
      T12 = 1/T52
      T7 = (ARO(4) + BRO(1))**2
      T54 = SQRT(T1 + T7)
      T13 = CN1/T54
      T8 = CN3*AD(2)**2
      T53 = SQRT(T1 + T7 + T8)
      T10 = T12 + T13 + CN2/T53
      W(4) = T10
      W(5) = T10
      T61 = T2 + CN5*AD(2)
      T15 = T61**2
      T62 = SQRT(T15 + T7)
      T63 = T2 + CN6*AD(2)
      T17 = T63**2
      T64 = SQRT(T17 + T7)
      W(6) = T12 + T13 + CN4/T62 + CN4/T64
      T18 = (ARO(8) + BRO(1))**2
      T70 = T2 - AD(4)
      T19 = T70**2
      T71 = SQRT(T18 + T19)
      T72 = T2 + AD(4)
      T20 = T72**2
      T73 = SQRT(T18 + T20)
      T22 = AD(4)**2
      T74 = SQRT(T1 + T18 + T22)
      W(7) = CN7*(CN4/T71 + CN4/T73 + CN1/T74)
      T23 = (ARO(7) + BRO(1))**2
      T78 = T2 - AD(3)
      T24 = T78**2
      T79 = SQRT(T23 + T24)
      T80 = T2 + AD(3)
      T25 = T80**2
      T81 = SQRT(T23 + T25)
      T27 = CN1/T79 + CN2/T81
      W(8) = T27
      W(9) = CN7*T27
      W(10) = T27
      T29 = (ARO(9) + BRO(1))**2
      T86 = T2 - AD(5)
      T30 = T86**2
      T87 = SQRT(T29 + T30)
      T88 = T2 + AD(5)
      T31 = T88**2
      T89 = SQRT(T29 + T31)
      T33 = AD(5)**2
      T90 = SQRT(T1 + T29 + T33)
      T34 = CN4/T87 + CN4/T89 + CN1/T90
      T91 = SQRT(T1 + (ARO(6) + BRO(1))**2)
      T35 = 1/T91
      T40 = CN8*T34 + T35
      W(11) = T40
      T38 = CN9*T34 + T35
      W(12) = T38
      W(13) = CN10*T34 + T35
      W(14) = T38
      W(15) = T40
      T43 = SQRT(T1 + (ARO(1) + BRO(2))**2)
      W(16) = 1/T43
      T45 = SQRT(T1 + (ARO(2) + BRO(2))**2)
      W(17) = 1/T45
      T6 = (ARO(3) + BRO(2))**2
      T50 = SQRT(T4 + T6)
      T51 = SQRT(T5 + T6)
      W(18) = CN1/T50 + CN2/T51
      T55 = SQRT(T1 + (ARO(5) + BRO(2))**2)
      T14 = 1/T55
      T9 = (ARO(4) + BRO(2))**2
      T57 = SQRT(T1 + T9)
      T16 = CN1/T57
      T56 = SQRT(T1 + T8 + T9)
      T11 = T14 + T16 + CN2/T56
      W(19) = T11
      W(20) = T11
      T67 = SQRT(T15 + T9)
      T68 = SQRT(T17 + T9)
      W(21) = T14 + T16 + CN4/T67 + CN4/T68
      T21 = (ARO(8) + BRO(2))**2
      T75 = SQRT(T19 + T21)
      T76 = SQRT(T20 + T21)
      T77 = SQRT(T1 + T21 + T22)
      W(22) = CN7*(CN4/T75 + CN4/T76 + CN1/T77)
      T26 = (ARO(7) + BRO(2))**2
      T82 = SQRT(T24 + T26)
      T83 = SQRT(T25 + T26)
      T28 = CN1/T82 + CN2/T83
      W(23) = T28
      W(24) = CN7*T28
      W(25) = T28
      T32 = (ARO(9) + BRO(2))**2
      T92 = SQRT(T30 + T32)
      T93 = SQRT(T31 + T32)
      T94 = SQRT(T1 + T32 + T33)
      T36 = CN4/T92 + CN4/T93 + CN1/T94
      T95 = SQRT(T1 + (ARO(6) + BRO(2))**2)
      T37 = 1/T95
      T41 = CN8*T36 + T37
      W(26) = T41
      T39 = CN9*T36 + T37
      W(27) = T39
      W(28) = CN10*T36 + T37
      W(29) = T39
      W(30) = T41
      IF( IMOD.LE.0 ) RETURN
C
C   Expressions for first derivatives
C
      T104 = T42**2
      T105 = T104*T42
      T106 = 1/T105
      W1(1) = -(R*T106)
      T110 = T44**2
      T111 = T110*T44
      T112 = 1/T111
      W1(2) = -(R*T112)
      T116 = T47**2
      T117 = T116*T47
      T118 = 1/T117
      T119 = T49**2
      T120 = T119*T49
      T121 = 1/T120
      W1(3) = CN1*T118*T46 + CN2*T121*T48
      T130 = T53**2
      T131 = T130*T53
      T136 = 1/T131
      T133 = T52**2
      T134 = T133*T52
      T135 = 1/T134
      T60 = -(R*T135)
      T128 = T54**2
      T129 = T128*T54
      T132 = 1/T129
      T65 = CN2*R*T132
      T58 = CN1*R*T136 + T60 + T65
      W1(4) = T58
      W1(5) = T58
      T154 = T62**2
      T155 = T154*T62
      T156 = 1/T155
      T150 = T64**2
      T151 = T150*T64
      T152 = 1/T151
      W1(6) = T60 + CN4*T156*T61 + CN4*T152*T63 + T65
      T168 = T74**2
      T169 = T168*T74
      T174 = 1/T169
      T166 = T71**2
      T167 = T166*T71
      T170 = 1/T167
      T171 = T73**2
      T172 = T171*T73
      T173 = 1/T172
      W1(7) = CN7*(CN2*R*T174 + CN4*T170*T70 + CN4*T173*T72)
      T184 = T79**2
      T185 = T184*T79
      T186 = 1/T185
      T187 = T81**2
      T188 = T187*T81
      T189 = 1/T188
      T84 = CN1*T186*T78 + CN2*T189*T80
      W1(8) = T84
      W1(9) = CN7*T84
      W1(10) = T84
      T202 = T90**2
      T203 = T202*T90
      T208 = 1/T203
      T200 = T87**2
      T201 = T200*T87
      T204 = 1/T201
      T205 = T89**2
      T206 = T205*T89
      T207 = 1/T206
      T96 = CN2*R*T208 + CN4*T204*T86 + CN4*T207*T88
      T198 = T91**2
      T199 = T198*T91
      T209 = 1/T199
      T97 = -(R*T209)
      T102 = CN8*T96 + T97
      W1(11) = T102
      T100 = CN9*T96 + T97
      W1(12) = T100
      W1(13) = CN10*T96 + T97
      W1(14) = T100
      W1(15) = T102
      T107 = T43**2
      T108 = T107*T43
      T109 = 1/T108
      W1(16) = -(R*T109)
      T113 = T45**2
      T114 = T113*T45
      T115 = 1/T114
      W1(17) = -(R*T115)
      T122 = T50**2
      T123 = T122*T50
      T124 = 1/T123
      T125 = T51**2
      T126 = T125*T51
      T127 = 1/T126
      W1(18) = CN1*T124*T46 + CN2*T127*T48
      T139 = T56**2
      T140 = T139*T56
      T145 = 1/T140
      T142 = T55**2
      T143 = T142*T55
      T144 = 1/T143
      T66 = -(R*T144)
      T137 = T57**2
      T138 = T137*T57
      T141 = 1/T138
      T69 = CN2*R*T141
      T59 = CN1*R*T145 + T66 + T69
      W1(19) = T59
      W1(20) = T59
      T163 = T67**2
      T164 = T163*T67
      T165 = 1/T164
      T159 = T68**2
      T160 = T159*T68
      T161 = 1/T160
      W1(21) = CN4*T165*T61 + CN4*T161*T63 + T66 + T69
      T177 = T77**2
      T178 = T177*T77
      T183 = 1/T178
      T175 = T75**2
      T176 = T175*T75
      T179 = 1/T176
      T180 = T76**2
      T181 = T180*T76
      T182 = 1/T181
      W1(22) = CN7*(CN2*R*T183 + CN4*T179*T70 + CN4*T182*T72)
      T190 = T82**2
      T191 = T190*T82
      T192 = 1/T191
      T193 = T83**2
      T194 = T193*T83
      T195 = 1/T194
      T85 = CN1*T192*T78 + CN2*T195*T80
      W1(23) = T85
      W1(24) = CN7*T85
      W1(25) = T85
      T214 = T94**2
      T215 = T214*T94
      T220 = 1/T215
      T212 = T92**2
      T213 = T212*T92
      T216 = 1/T213
      T217 = T93**2
      T218 = T217*T93
      T219 = 1/T218
      T98 = CN2*R*T220 + CN4*T216*T86 + CN4*T219*T88
      T210 = T95**2
      T211 = T210*T95
      T221 = 1/T211
      T99 = -(R*T221)
      T103 = CN8*T98 + T99
      W1(26) = T103
      T101 = CN9*T98 + T99
      W1(27) = T101
      W1(28) = CN10*T98 + T99
      W1(29) = T101
      W1(30) = T103
      IF( IMOD.LE.1 ) RETURN
C
C   Expressions for second derivatives
C
      W2(1) = CN11*T1/(T104*T105) - T106
      W2(2) = CN11*T1/(T110*T111) - T112
      W2(3) = CN2*T118 + CN1*T121 + CN12*T4/(T116*T117) + 
     1 CN13*T5/(T119*T120)
      T148 = CN12*T1/(T128*T129)
      T149 = CN2*T132
      T153 = CN11*T1/(T133*T134) - T135
      T146 = CN13*T1/(T130*T131) + CN1*T136 + T148 + T149 + T153
      W2(4) = T146
      W2(5) = T146
      W2(6) = T148 + T149 + CN15*T152 + T153 + CN14*T15/(T154*T155) + 
     1 CN15*T156 + CN14*T17/(T150*T151)
      W2(7) = CN7*(CN12*T1/(T168*T169) + CN15*T170 + CN15*T173 + 
     1 CN2*T174 + CN14*T19/(T166*T167) + CN14*T20/(T171*T172))
      T196 = CN2*T186 + CN1*T189 + CN12*T24/(T184*T185) + 
     1 CN13*T25/(T187*T188)
      W2(8) = T196
      W2(9) = CN7*T196
      W2(10) = T196
      T222 = CN11*T1/(T198*T199)
      T223 = CN12*T1/(T202*T203) + CN15*T204 + CN15*T207 + CN2*T208 + 
     1 CN14*T30/(T200*T201) + CN14*T31/(T205*T206)
      T224 = -T209
      T230 = T222 + CN8*T223 + T224
      W2(11) = T230
      T228 = T222 + CN9*T223 + T224
      W2(12) = T228
      W2(13) = T222 + CN10*T223 + T224
      W2(14) = T228
      W2(15) = T230
      W2(16) = CN11*T1/(T107*T108) - T109
      W2(17) = CN11*T1/(T113*T114) - T115
      W2(18) = CN2*T124 + CN1*T127 + CN12*T4/(T122*T123) + 
     1 CN13*T5/(T125*T126)
      T157 = CN12*T1/(T137*T138)
      T158 = CN2*T141
      T162 = CN11*T1/(T142*T143) - T144
      T147 = CN13*T1/(T139*T140) + CN1*T145 + T157 + T158 + T162
      W2(19) = T147
      W2(20) = T147
      W2(21) = T157 + T158 + CN15*T161 + T162 + CN14*T15/(T163*T164) + 
     1 CN15*T165 + CN14*T17/(T159*T160)
      W2(22) = CN7*(CN12*T1/(T177*T178) + CN15*T179 + CN15*T182 + 
     1 CN2*T183 + CN14*T19/(T175*T176) + CN14*T20/(T180*T181))
      T197 = CN2*T192 + CN1*T195 + CN12*T24/(T190*T191) + 
     1 CN13*T25/(T193*T194)
      W2(23) = T197
      W2(24) = CN7*T197
      W2(25) = T197
      T225 = CN11*T1/(T210*T211)
      T226 = CN12*T1/(T214*T215) + CN15*T216 + CN15*T219 + CN2*T220 + 
     1 CN14*T30/(T212*T213) + CN14*T31/(T217*T218)
      T227 = -T221
      T231 = T225 + CN8*T226 + T227
      W2(26) = T231
      T229 = T225 + CN9*T226 + T227
      W2(27) = T229
      W2(28) = T225 + CN10*T226 + T227
      W2(29) = T229
      W2(30) = T231
      RETURN
      END