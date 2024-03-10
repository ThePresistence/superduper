C     ******************************************************************
      SUBROUTINE PSILPX(R,ARO,AD,IMOD,W,W1,W2)
C
C  This subroutine is not intended to be human-readable.
C  It was generated by a Mathematica program (makeints.m).
C  General version of 14 May 1997.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARO(9), AD(5)
      DIMENSION W(461),W1(461),W2(461)
      PARAMETER (CN1=4.D0)
      PARAMETER (CN2=0.5D0)
      PARAMETER (CN3=-0.5D0)
      PARAMETER (CN4=2.D0)
      PARAMETER (CN5=0.25D0)
      PARAMETER (CN6=1.4142135623730950488016887242096980785697D0)
      PARAMETER (CN7=-1.4142135623730950488016887242096980785697D0)
      PARAMETER (CN8=-0.25D0)
      PARAMETER (CN9=-2.D0)
      PARAMETER (CN10=-0.125D0)
      PARAMETER (CN11=0.125D0)
      PARAMETER (CN12=0.375D0)
      PARAMETER (CN13=8.D0)
      PARAMETER (CN14=-2.828427124746190097603377448419396157139D0)
      PARAMETER (CN15=2.828427124746190097603377448419396157139D0)
      PARAMETER (CN16=0.0625D0)
      PARAMETER (CN17=-0.375D0)
      PARAMETER (CN18=-0.0625D0)
      PARAMETER (CN19=-4.D0)
      PARAMETER (CN20=3.D0)
      PARAMETER (CN21=1.5D0)
      PARAMETER (CN22=-1.5D0)
      PARAMETER (CN23=0.75D0)
      PARAMETER (CN24=-0.1875D0)
      PARAMETER (CN25=0.1875D0)
      PARAMETER (CN26=-0.75D0)
      PARAMETER (CN27=1.125D0)
      PARAMETER (CN28=-3.D0)
      T1 = R**2
      T225 = ARO(1)**2
      T97 = SQRT(T1 + CN1*T225)
      W(1) = 1/T97
      T98 = SQRT(T1 + (ARO(1) + ARO(2))**2)
      T8 = 1/T98
      W(2) = T8
      T3 = (ARO(1) + ARO(3))**2
      T2 = -R
      T44 = -AD(1)
      T101 = T2 + T44
      T9 = T101**2
      T102 = SQRT(T3 + T9)
      T33 = 1/T102
      T99 = T2 + AD(1)
      T10 = T99**2
      T100 = SQRT(T10 + T3)
      T34 = 1/T100
      W(3) = CN3*T33 + CN2*T34
      T21 = AD(2)**2
      T12 = CN4*T21
      T4 = (ARO(1) + ARO(4))**2
      T104 = SQRT(T1 + T12 + T4)
      T103 = SQRT(T1 + (ARO(1) + ARO(5))**2)
      T6 = 1/T103
      T105 = SQRT(T1 + T4)
      T7 = CN3/T105
      T5 = CN2/T104 + T6 + T7
      W(4) = T5
      W(5) = T5
      T47 = CN6*AD(2)
      T108 = T2 + T47
      T16 = T108**2
      T109 = SQRT(T16 + T4)
      T45 = CN7*AD(2)
      T110 = T2 + T45
      T18 = T110**2
      T111 = SQRT(T18 + T4)
      T83 = CN5/T109 + CN5/T111 + T6 + T7
      W(6) = T83
      W(16) = T8
      T258 = ARO(2)**2
      T114 = SQRT(T1 + CN1*T258)
      W(17) = 1/T114
      T11 = (ARO(2) + ARO(3))**2
      T116 = SQRT(T11 + T9)
      T35 = 1/T116
      T115 = SQRT(T10 + T11)
      T36 = 1/T115
      W(18) = CN3*T35 + CN2*T36
      T13 = (ARO(2) + ARO(4))**2
      T118 = SQRT(T1 + T12 + T13)
      T117 = SQRT(T1 + (ARO(2) + ARO(5))**2)
      T15 = 1/T117
      T119 = SQRT(T1 + T13)
      T17 = CN3/T119
      T14 = CN2/T118 + T15 + T17
      W(19) = T14
      W(20) = T14
      T122 = SQRT(T13 + T16)
      T123 = SQRT(T13 + T18)
      T84 = CN5/T122 + CN5/T123 + T15 + T17
      W(21) = T84
      T22 = ARO(3)**2
      T19 = CN1*T22
      T20 = AD(1)**2
      T126 = SQRT(T1 + T19 + CN1*T20)
      T125 = SQRT(T1 + T19)
      T37 = CN2/T125
      T31 = CN3/T126 + T37
      W(31) = T31
      T23 = T1 + T20
      T28 = AD(1)*AD(2)
      T24 = CN7*T28
      T54 = ARO(4)**2
      T26 = T21 + T22 + T54 + CN4*ARO(3)*ARO(4)
      T25 = R*AD(2)
      T27 = CN7*T25
      T129 = SQRT(T23 + T24 + T26 + T27)
      T72 = 1/T129
      T29 = CN6*T25
      T127 = SQRT(T23 + T24 + T26 + T29)
      T73 = 1/T127
      T30 = CN6*T28
      T131 = SQRT(T23 + T26 + T27 + T30)
      T74 = 1/T131
      T133 = SQRT(T23 + T26 + T29 + T30)
      T75 = 1/T133
      T81 = CN5*T72 + CN8*T73 + CN8*T74 + CN5*T75
      W(32) = T81
      W(41) = T31
      W(42) = T81
      W(51) = CN2*T33 + CN3*T34
      W(52) = CN2*T35 + CN3*T36
      T141 = T2 + CN9*AD(1)
      T312 = T141**2
      T142 = SQRT(T19 + T312)
      T143 = T2 + CN4*AD(1)
      T317 = T143**2
      T144 = SQRT(T19 + T317)
      W(53) = CN8/T142 + CN8/T144 + T37
      T38 = (ARO(3) + ARO(4))**2
      T40 = T12 + T38
      T146 = SQRT(T40 + T9)
      T50 = 1/T146
      T148 = SQRT(T10 + T40)
      T53 = 1/T148
      T147 = SQRT(T38 + T9)
      T48 = 1/T147
      T85 = CN8*T48
      T145 = SQRT(T10 + T38)
      T49 = 1/T145
      T86 = CN5*T49
      T39 = (ARO(3) + ARO(5))**2
      T149 = SQRT(T39 + T9)
      T51 = 1/T149
      T150 = SQRT(T10 + T39)
      T52 = 1/T150
      T89 = CN2*T51 + CN3*T52
      T68 = CN5*T50 + CN8*T53 + T85 + T86 + T89
      W(54) = T68
      W(55) = T68
      T157 = T2 + T44 + T45
      T344 = T157**2
      T158 = SQRT(T344 + T38)
      T87 = 1/T158
      T159 = T2 + T45 + AD(1)
      T348 = T159**2
      T160 = SQRT(T348 + T38)
      T88 = 1/T160
      T152 = T2 + T44 + T47
      T355 = T152**2
      T153 = SQRT(T355 + T38)
      T90 = 1/T153
      T154 = T2 + T47 + AD(1)
      T359 = T154**2
      T155 = SQRT(T359 + T38)
      T91 = 1/T155
      W(56) = T85 + T86 + CN11*T87 + CN10*T88 + T89 + CN11*T90 + 
     1 CN10*T91
      W(66) = T5
      W(67) = T14
      T42 = CN5*T48
      T43 = CN8*T49
      T46 = CN3*T51 + CN2*T52
      T41 = T42 + T43 + T46 + CN8*T50 + CN5*T53
      W(68) = T41
      T55 = CN1*T54
      T174 = SQRT(T1 + T12 + T55)
      T65 = 1/T174
      T58 = CN3*T65
      T169 = SQRT(T1 + CN1*ARO(5)**2)
      T60 = 1/T169
      T56 = (ARO(4) + ARO(5))**2
      T171 = SQRT(T1 + T56)
      T61 = -(1/T171)
      T170 = SQRT(T1 + T12 + T56)
      T62 = 1/T170
      T59 = T60 + T61 + T62
      T173 = SQRT(T1 + CN13*T21 + T55)
      T67 = 1/T173
      T172 = SQRT(T1 + T55)
      T57 = 1/T172
      T94 = CN12*T57
      T70 = T58 + T59 + CN11*T67 + T94
      W(69) = T70
      T64 = CN5*T57
      T63 = CN1*T21 + T55
      T177 = SQRT(T1 + T63)
      T66 = 1/T177
      T69 = T58 + T59 + T64 + CN5*T66
      W(70) = T69
      T76 = CN8*T65
      T185 = SQRT(T1 + CN14*T25 + T63)
      T77 = CN11/T185
      T186 = SQRT(T1 + CN15*T25 + T63)
      T78 = CN11/T186
      T182 = SQRT(T16 + T55)
      T95 = 1/T182
      T79 = CN10*T95
      T187 = SQRT(T18 + T55)
      T96 = 1/T187
      T80 = CN10*T96
      T180 = SQRT(T16 + T56)
      T92 = 1/T180
      T181 = SQRT(T18 + T56)
      T93 = 1/T181
      T71 = T60 + T61 + CN2*T62 + T64 + T76 + T77 + T78 + T79 + T80 + 
     1 CN5*T92 + CN5*T93
      W(71) = T71
      W(83) = CN16*T57 + CN10*T66 + CN16*T67
      W(87) = T5
      W(88) = T14
      W(89) = T41
      W(90) = T69
      W(91) = T70
      W(92) = T71
      T32 = CN8*T72 + CN5*T73 + CN5*T74 + CN8*T75
      W(104) = T32
      T82 = T64 + T76 + T77 + T78 + T79 + T80
      W(105) = T82
      W(114) = T32
      W(115) = T82
      W(124) = T83
      W(125) = T84
      W(126) = T42 + T43 + T46 + CN10*T87 + CN11*T88 + CN10*T90 + 
     1 CN11*T91
      W(127) = T71
      W(128) = T71
      T218 = T2 + CN15*AD(2)
      T458 = T218**2
      T219 = SQRT(T458 + T55)
      T223 = T2 + CN14*AD(2)
      T464 = T223**2
      T224 = SQRT(T464 + T55)
      W(129) = CN16/T219 + CN16/T224 + T60 + T61 + CN2*T92 + CN2*T93 + 
     1 T94 + CN8*T95 + CN8*T96
      IF( IMOD.LE.0 ) RETURN
C
C   Expressions for first derivatives
C
      T226 = T97**2
      T227 = T226*T97
      W1(1) = -(R/T227)
      T228 = T98**2
      T229 = T228*T98
      T230 = 1/T229
      T113 = -(R*T230)
      W1(2) = T113
      T231 = T102**2
      T232 = T102*T231
      T233 = 1/T232
      T136 = T101*T233
      T234 = T100**2
      T235 = T100*T234
      T236 = 1/T235
      T137 = T236*T99
      W1(3) = CN3*T136 + CN2*T137
      T242 = T103**2
      T243 = T103*T242
      T244 = 1/T243
      T107 = -(R*T244)
      T237 = T105**2
      T238 = T105*T237
      T241 = 1/T238
      T112 = CN2*R*T241
      T239 = T104**2
      T240 = T104*T239
      T245 = 1/T240
      T106 = T107 + T112 + CN3*R*T245
      W1(4) = T106
      W1(5) = T106
      T249 = T111**2
      T250 = T111*T249
      T251 = 1/T250
      T253 = T109**2
      T254 = T109*T253
      T255 = 1/T254
      T207 = T107 + T112 + CN5*T110*T251 + CN5*T108*T255
      W1(6) = T207
      W1(16) = T113
      T259 = T114**2
      T260 = T114*T259
      W1(17) = -(R/T260)
      T261 = T116**2
      T262 = T116*T261
      T263 = 1/T262
      T138 = T101*T263
      T264 = T115**2
      T265 = T115*T264
      T266 = 1/T265
      T139 = T266*T99
      W1(18) = CN3*T138 + CN2*T139
      T272 = T117**2
      T273 = T117*T272
      T274 = 1/T273
      T121 = -(R*T274)
      T267 = T119**2
      T268 = T119*T267
      T271 = 1/T268
      T124 = CN2*R*T271
      T269 = T118**2
      T270 = T118*T269
      T275 = 1/T270
      T120 = T121 + T124 + CN3*R*T275
      W1(19) = T120
      W1(20) = T120
      T279 = T123**2
      T280 = T123*T279
      T281 = 1/T280
      T283 = T122**2
      T284 = T122*T283
      T285 = 1/T284
      T208 = T121 + T124 + CN5*T110*T281 + CN5*T108*T285
      W1(21) = T208
      T286 = T125**2
      T287 = T125*T286
      T288 = 1/T287
      T140 = CN3*R*T288
      T289 = T126**2
      T290 = T126*T289
      T291 = 1/T290
      T134 = T140 + CN2*R*T291
      W1(31) = T134
      T128 = CN4*R
      T132 = T128 + T47
      T294 = T127**2
      T295 = T127*T294
      T296 = 1/T295
      T197 = T132*T296
      T302 = T133**2
      T303 = T133*T302
      T304 = 1/T303
      T198 = T132*T304
      T130 = T128 + T45
      T298 = T131**2
      T299 = T131*T298
      T300 = 1/T299
      T199 = T130*T300
      T292 = T129**2
      T293 = T129*T292
      T196 = 1/T293
      T427 = CN9*R + T47
      T205 = CN11*T197 + CN10*T198 + CN11*T199 + CN11*T196*T427
      W1(32) = T205
      W1(41) = T134
      W1(42) = T205
      W1(51) = CN2*T136 + CN3*T137
      W1(52) = CN2*T138 + CN3*T139
      T313 = T142**2
      T314 = T142*T313
      T316 = 1/T314
      T318 = T144**2
      T319 = T144*T318
      T320 = 1/T319
      W1(53) = T140 + CN8*T141*T316 + CN8*T143*T320
      T327 = T146**2
      T328 = T146*T327
      T329 = 1/T328
      T165 = T101*T329
      T330 = T148**2
      T331 = T148*T330
      T334 = 1/T331
      T166 = T334*T99
      T336 = T150**2
      T337 = T150*T336
      T338 = 1/T337
      T167 = T338*T99
      T332 = T149**2
      T333 = T149*T332
      T335 = 1/T333
      T168 = T101*T335
      T211 = CN3*T167 + CN2*T168
      T321 = T147**2
      T322 = T147*T321
      T325 = 1/T322
      T163 = T101*T325
      T214 = CN8*T163
      T323 = T145**2
      T324 = T145*T323
      T326 = 1/T324
      T164 = T326*T99
      T215 = CN5*T164
      T192 = CN5*T165 + CN8*T166 + T211 + T214 + T215
      W1(54) = T192
      W1(55) = T192
      T360 = T155**2
      T361 = T155*T360
      T362 = 1/T361
      T209 = T154*T362
      T356 = T153**2
      T357 = T153*T356
      T358 = 1/T357
      T210 = T152*T358
      T349 = T160**2
      T350 = T160*T349
      T351 = 1/T350
      T212 = T159*T351
      T345 = T158**2
      T346 = T158*T345
      T347 = 1/T346
      T213 = T157*T347
      W1(56) = CN10*T209 + CN11*T210 + T211 + CN10*T212 + CN11*T213 + 
     1 T214 + T215
      W1(66) = T106
      W1(67) = T120
      T156 = CN2*T167 + CN3*T168
      T161 = CN8*T164
      T162 = CN5*T163
      T151 = T156 + T161 + T162 + CN8*T165 + CN5*T166
      W1(68) = T151
      T384 = T170**2
      T385 = T170*T384
      T386 = 1/T385
      T179 = R*T386
      T378 = T169**2
      T379 = T169*T378
      T382 = 1/T379
      T183 = -(R*T382)
      T380 = T171**2
      T381 = T171*T380
      T383 = 1/T381
      T184 = R*T383
      T175 = -T179 + T183 + T184
      T371 = T174**2
      T372 = T174*T371
      T374 = 1/T372
      T189 = R*T374
      T178 = CN2*T189
      T375 = T173**2
      T376 = T173*T375
      T377 = 1/T376
      T190 = R*T377
      T369 = T172**2
      T370 = T172*T369
      T373 = 1/T370
      T176 = R*T373
      T222 = CN17*T176
      T194 = T175 + T178 + CN10*T190 + T222
      W1(69) = T194
      T188 = CN8*T176
      T390 = T177**2
      T391 = T177*T390
      T392 = 1/T391
      T191 = R*T392
      T193 = T175 + T178 + T188 + CN8*T191
      W1(70) = T193
      T395 = R + T45
      T396 = T185**2
      T397 = T185*T396
      T405 = 1/T397
      T200 = CN10*T395*T405
      T201 = CN5*T189
      T406 = R + T47
      T407 = T186**2
      T408 = T186*T407
      T409 = 1/T408
      T202 = CN10*T406*T409
      T402 = T182**2
      T403 = T182*T402
      T420 = 1/T403
      T220 = T108*T420
      T203 = CN10*T220
      T400 = T187**2
      T401 = T187*T400
      T404 = 1/T401
      T221 = T110*T404
      T204 = CN10*T221
      T413 = T180**2
      T414 = T180*T413
      T415 = 1/T414
      T216 = T108*T415
      T411 = T181**2
      T412 = T181*T411
      T416 = 1/T412
      T217 = T110*T416
      T195 = CN3*T179 + T183 + T184 + T188 + T200 + T201 + T202 + 
     1 T203 + T204 + CN5*T216 + CN5*T217
      W1(71) = T195
      W1(83) = CN18*T176 + CN18*T190 + CN11*T191
      W1(87) = T106
      W1(88) = T120
      W1(89) = T151
      W1(90) = T193
      W1(91) = T194
      W1(92) = T195
      T135 = CN11*T130*T196 + CN10*T197 + CN11*T198 + CN10*T199
      W1(104) = T135
      T206 = T188 + T200 + T201 + T202 + T203 + T204
      W1(105) = T206
      W1(114) = T135
      W1(115) = T206
      W1(124) = T207
      W1(125) = T208
      W1(126) = T156 + T161 + T162 + CN11*T209 + CN10*T210 + 
     1 CN11*T212 + CN10*T213
      W1(127) = T195
      W1(128) = T195
      T465 = T224**2
      T466 = T224*T465
      T467 = 1/T466
      T459 = T219**2
      T460 = T219*T459
      T471 = 1/T460
      W1(129) = T183 + T184 + CN2*T216 + CN2*T217 + CN8*T220 + 
     1 CN8*T221 + T222 + CN16*T223*T467 + CN16*T218*T471
      IF( IMOD.LE.1 ) RETURN
C
C   Expressions for second derivatives
C
      T257 = CN4*T1
      W2(1) = (CN19*T225 + T257)/(T226*T227)
      T256 = CN20*T1/(T228*T229) - T230
      W2(2) = T256
      T307 = T9/(T231*T232)
      T308 = T10/(T234*T235)
      W2(3) = CN2*T233 + CN3*T236 + CN22*T307 + CN21*T308
      T247 = CN22*T1/(T237*T238)
      T248 = CN2*T241
      T252 = CN20*T1/(T242*T243) - T244
      T246 = CN21*T1/(T239*T240) + CN3*T245 + T247 + T248 + T252
      W2(4) = T246
      W2(5) = T246
      T444 = T247 + T248 + CN23*T18/(T249*T250) + CN8*T251 + T252 + 
     1 CN23*T16/(T253*T254) + CN8*T255
      W2(6) = T444
      W2(16) = T256
      W2(17) = (T257 + CN19*T258)/(T259*T260)
      T309 = T9/(T261*T262)
      T310 = T10/(T264*T265)
      W2(18) = CN2*T263 + CN3*T266 + CN22*T309 + CN21*T310
      T277 = CN22*T1/(T267*T268)
      T278 = CN2*T271
      T282 = CN20*T1/(T272*T273) - T274
      T276 = CN21*T1/(T269*T270) + CN3*T275 + T277 + T278 + T282
      W2(19) = T276
      W2(20) = T276
      T445 = T277 + T278 + CN23*T18/(T279*T280) + CN8*T281 + T282 + 
     1 CN23*T16/(T283*T284) + CN8*T285
      W2(21) = T445
      T311 = CN21*T1/(T286*T287)
      T315 = CN3*T288
      T305 = CN22*T1/(T289*T290) + CN2*T291 + T311 + T315
      W2(31) = T305
      T428 = 1/(T292*T293)
      T301 = T132**2
      T429 = T301/(T294*T295)
      T297 = T130**2
      T430 = T297/(T298*T299)
      T431 = T301/(T302*T303)
      T442 = CN8*T196 + CN5*T296 + CN5*T300 + CN8*T304 + 
     1 CN25*T427**2*T428 + CN24*T429 + CN24*T430 + CN25*T431
      W2(32) = T442
      W2(41) = T305
      W2(42) = T442
      W2(51) = CN3*T233 + CN2*T236 + CN21*T307 + CN22*T308
      W2(52) = CN3*T263 + CN2*T266 + CN21*T309 + CN22*T310
      W2(53) = T311 + CN26*T312/(T313*T314) + T315 + CN5*T316 + 
     1 CN26*T317/(T318*T319) + CN5*T320
      T365 = T9/(T327*T328)
      T366 = T10/(T330*T331)
      T363 = T9/(T321*T322)
      T446 = CN26*T363
      T447 = CN5*T325
      T364 = T10/(T323*T324)
      T448 = CN23*T364
      T449 = CN8*T326
      T367 = T9/(T332*T333)
      T452 = CN21*T367
      T453 = CN3*T335
      T368 = T10/(T336*T337)
      T454 = CN2*T338 + CN22*T368
      T423 = CN8*T329 + CN5*T334 + CN23*T365 + CN26*T366 + T446 + 
     1 T447 + T448 + T449 + T452 + T453 + T454
      W2(54) = T423
      W2(55) = T423
      T450 = T344/(T345*T346)
      T451 = T348/(T349*T350)
      T455 = T355/(T356*T357)
      T456 = T359/(T360*T361)
      W2(56) = CN10*T347 + CN11*T351 + CN10*T358 + CN11*T362 + T446 + 
     1 T447 + T448 + T449 + CN12*T450 + CN17*T451 + T452 + T453 + 
     2 T454 + CN12*T455 + CN17*T456
      W2(66) = T246
      W2(67) = T276
      T340 = CN23*T363
      T341 = CN8*T325
      T342 = CN26*T364
      T343 = CN5*T326
      T352 = CN22*T367
      T353 = CN2*T335
      T354 = CN3*T338 + CN21*T368
      T339 = CN5*T329 + CN8*T334 + T340 + T341 + T342 + T343 + T352 + 
     1 T353 + T354 + CN26*T365 + CN23*T366
      W2(68) = T339
      T394 = T1/(T371*T372)
      T388 = CN22*T394
      T389 = CN2*T374
      T410 = T1/(T384*T385)
      T419 = CN20*T1/(T378*T379)
      T417 = -T382
      T418 = CN28*T1/(T380*T381)
      T470 = T417 + T418
      T393 = T383 - T386 + CN20*T410 + T419 + T470
      T422 = T1/(T375*T376)
      T387 = T1/(T369*T370)
      T462 = CN27*T387
      T463 = CN17*T373
      T425 = CN10*T377 + T388 + T389 + T393 + CN12*T422 + T462 + T463
      W2(69) = T425
      T398 = CN23*T387
      T399 = CN8*T373
      T421 = T1/(T390*T391)
      T424 = T388 + T389 + CN8*T392 + T393 + T398 + T399 + CN23*T421
      W2(70) = T424
      T432 = CN12*T406**2/(T407*T408)
      T461 = T18/(T400*T401)
      T433 = CN17*T461
      T434 = CN10*T409
      T457 = T16/(T402*T403)
      T435 = CN17*T457
      T436 = CN11*T420
      T437 = CN11*T404
      T438 = CN26*T394
      T439 = CN5*T374
      T440 = CN12*T395**2/(T396*T397)
      T441 = CN10*T405
      T468 = T16/(T413*T414)
      T469 = T18/(T411*T412)
      T426 = T383 + CN3*T386 + T398 + T399 + CN21*T410 + CN8*T415 + 
     1 CN8*T416 + T417 + T418 + T419 + T432 + T433 + T434 + T435 + 
     2 T436 + T437 + T438 + T439 + T440 + T441 + CN23*T468 + CN23*T469
      W2(71) = T426
      W2(83) = CN18*T373 + CN18*T377 + CN25*T387 + CN11*T392 + 
     1 CN17*T421 + CN25*T422
      W2(87) = T246
      W2(88) = T276
      W2(89) = T339
      W2(90) = T424
      W2(91) = T425
      W2(92) = T426
      T306 = CN5*T196 + CN8*T296 + CN8*T300 + CN5*T304 + 
     1 CN24*T297*T428 + CN25*T429 + CN25*T430 + CN24*T431
      W2(104) = T306
      T443 = T398 + T399 + T432 + T433 + T434 + T435 + T436 + T437 + 
     1 T438 + T439 + T440 + T441
      W2(105) = T443
      W2(114) = T306
      W2(115) = T443
      W2(124) = T444
      W2(125) = T445
      W2(126) = T340 + T341 + T342 + T343 + CN11*T347 + CN10*T351 + 
     1 T352 + T353 + T354 + CN11*T358 + CN10*T362 + CN17*T450 + 
     2 CN12*T451 + CN17*T455 + CN12*T456
      W2(127) = T426
      W2(128) = T426
      W2(129) = T383 + CN5*T404 + CN3*T415 + CN3*T416 + T419 + 
     1 CN5*T420 + CN26*T457 + CN25*T458/(T459*T460) + CN26*T461 + 
     2 T462 + T463 + CN25*T464/(T465*T466) + CN18*T467 + CN21*T468 + 
     3 CN21*T469 + T470 + CN18*T471
      RETURN
      END