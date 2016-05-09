### Tests for Sweep operator & friends.

## These tests seem to work with the Goodnight Version and not the
## Dempster version.
### Test taken from http://www-rci.rutgers.edu/~dhjones/APPLIED_LINEAR_STATISTICAL_MODELS%28PHD%29/LECTURES/LECTURE06/3-The%20sweep%20operator.pdf
##
A11 <- matrix(c(2,3, 4,1),2,2)
A12 <- matrix(c(5,6),2,1)
A22 <- 70

A0 <- rbind(cbind(A11,A12),c(t(A12),70))

## Expected results
A1 <- matrix(c(.5,-1.5,-2.5, 2,-5,-4, 2.5,-1.5,57.5),3,3)
A2 <- matrix(c(-.1,.3,-1.3, .4,-.2,-.8, 1.9,.3,58.7),3,3)

A0.1 <- matSweep(A0,1)
stopifnot(all(abs(A1-A0.1)<.Machine$double.eps^.5))

A0.12 <- matSweep(A0,1:2)
stopifnot(all(abs(A2-A0.12)<.Machine$double.eps^.5))


##### Example taken from Dempster 1969, p  64.

DQ <- matrix(c(19.1434,  9.0356,  9.7634, 3.2394,
                9.0356, 11.8658,  4.6232, 2.4746,
                9.7634,  4.6232, 12.2978, 3.8794,
                3.2394,  2.4746,  3.8794, 2.4604), 4,4)

DQ1 <- matrix(c(-0.0522373, 0.471995, 0.510014, 0.169218,
                 0.471995,  7.60104,  0.01492,  0.94562,
                 0.510014,  0.01492,  7.31833,  2.22726,
                 0.169218,  0.94562,  2.22726,  1.91224), 4,4)

DQ2 <- matrix(c(-0.0815463,  0.0620961, 0.509088,  0.110499,
                 0.0620961, -0.131561,  0.0019629, 0.1244067,
                 0.509088,   0.0019629, 7.318301,  2.225405,
                 0.110499,   0.1244067, 2.225405,  1.79460), 4,4)
## Art has 2.22574 for DQ2[3,4]
DQ2a <- DQ2
DQ2a[3,4] <-2.22574
DQ2a[4,3] <-2.22574

DQ3 <- matrix(c(-0.1169601,  0.0619595,  0.0695633, -0.044331,
                 0.0619595, -0.131562,   0.0002682,  0.123810,
                 0.0695633,  0.0002682, -0.136643,   0.304132,
                -0.044331,   0.123810,   0.304132,   1.117877), 4,4)
## Dempster has 1.11768 in DQ3[4,4]
DQ3a <- DQ3
DQ3a[4,4] <- 1.11768


DQ4 <- matrix(c(-0.118718,  0.066871,   0.081626, -0.039663,
                 0.066871, -0.145277,  -0.033422,  0.110774,
                 0.081626, -0.033422,  -0.219400,  0.272110,
                -0.039663,  0.110774,   0.272110, -0.894552), 4,4)
### Dempster has -0.894710 in DQ4[4,4]
DQ4a <- DQ4
DQ4a[4,4] <- -0.894710


stopifnot(all(abs(matSweep(DQ,1)-DQ1) < 1e-4))
stopifnot(all(abs(matSweep(DQ,1:2)-DQ2) < 1e-4))
stopifnot(all(abs(matSweep(DQ,1:3)-DQ3) < 1e-4))
stopifnot(all(abs(matSweep(DQ,1:4)-DQ4) < 1e-4))

stopifnot(all(abs(revSweep(DQ4,4)-DQ3) < 1e-4))
stopifnot(all(abs(revSweep(DQ3,3)-DQ2) < 1e-4))

stopifnot(all(abs(revSweep(DQ1,1)-DQ) < 1e-4))
stopifnot(all(abs(revSweep(DQ2,1:2)-DQ) < 1e-4))
stopifnot(all(abs(revSweep(DQ3,1:3)-DQ) < 1e-4))




DQT <- diag(c(19.1434, 7.60104, 7.318301, 1.117877))
DQB <- matrix(c(1,0.471995,0.510014,0.169218,
                0,1,0.0019629,0.1244067,
                0,0,1,0.304132,
                0,0,0,1),4,4)
DQA <- matrix(c(1,-0.471995,-0.509088,0.044331,
                0,1,-0.0019629,-0.123810,
                0,0,1,-0.304132,
                0,0,0,1),4,4)

oi <- orthoInvert(DQ)

stopifnot(all(abs(oi$T-DQT) < 1e-4))
stopifnot(all(abs(oi$B-DQB) < 1e-4))
stopifnot(all(abs(oi$A-DQA) < 1e-4))
stopifnot(all(abs(oi$Qinv+DQ4) < 1e-4))

