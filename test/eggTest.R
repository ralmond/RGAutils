data(eggs)
eggsPlus <- cbind(as.matrix(eggs),1)

eggsQplus <- t(eggsPlus)%*%eggsPlus

demp.eggQplus <- matrix(c( 6.9964,  5.6861, 18.3292,  9.1608,
                           5.6861,  4.6246, 14.9038,  7.4489,
                          18.3292, 14.9038, 48.0368, 24.0080,
                           9.1608,  7.4489, 24.0080, 12.0000),4,4)
stopifnot(max(abs(eggsQplus-demp.eggQplus))<.0001)

matSweep(eggsQplus,1)

demp.swp1Qp <- matrix(c(-0.1429, 0.8127,   2.6198,   1.3094,
                         0.8127, 0.003329, 0.007217, 0.003708,
                         2.6198, 0.007217, 0.01756,  0.008374,
                         1.3094, 0.003708, 0.008374, 0.005174), 4,4)
stopifnot(max(abs(matSweep(eggsQplus,1)-demp.swp1Qp))<.0001)

## Original as in Dempster (1968, p152).
demp.swp12Qpo <- matrix(c(-198.5754,  244.1576, 0.8577,    0.4040,
                           244.1576, -300.4193, 2.1682,    1.1139,
                             0.8577,    2.1682, 0.001910,  0.0003347,
                             0.4040,    1.1139, 0.0003347, 0.001044),4,4)

stopifnot(max(abs(matSweep(eggsQplus,c(1,2))-demp.swp12Qpo)) < .02)


## Original as in Dempster (1968, p152).
demp.swp123Qpo <- matrix(c(-583.3215,  -729.5001,  448.5989, 0.2538,
                           -729.5001, -2761.3479, 1134.0810, 0.7342,
                            448.5989,  1134.0810, -523.0488, 0.1751,
                              0.2538,     0.7342,    0.1751, 0.0009849),4,4)
stopifnot(max(abs(demp.swp123Qpo-matSweep(eggsQplus,1:3)))<1)

stopifnot(max(matSweep(eggsQplus,1:4)+solve(eggsQplus))<.00001)

eggsQswp1234 <- matSweep(eggsQplus,1:3)
revSweep(eggsQswp1234,3)
matSweep(eggsQplus,c(1,2,4))

revSweep(eggsQswp1234,4:1)
eggsQplus

