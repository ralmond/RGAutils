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
## sweep 1:4 should be the same as - the inverse.
stopifnot(all(abs(matSweep(DQ,1:4)+solve(DQ)) <1e-4))

stopifnot(all(abs(revSweep(DQ4,4)-DQ3) < 1e-4))
stopifnot(all(abs(revSweep(DQ3,3)-DQ2) < 5e-4))

stopifnot(all(abs(revSweep(DQ1,1)-DQ) < 1e-4))
stopifnot(all(abs(revSweep(DQ2,1:2)-DQ) < 1e-4))
stopifnot(all(abs(revSweep(DQ3,1:3)-DQ) < 5e-4))

stopifnot(all(abs(revSweep(matSweep(DQ,1:4),1:4)-DQ) < 1e-5))

stopifnot(max(abs(matSweep(DQ,c(1,2,4))-revSweep(DQ4,3))) < .001)




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


### Assimilation test from Dempster (1969, p. 69)

## DQ <- matrix(c(19.1434,  9.0356,  9.7634, 3.2394,
##                 9.0356, 11.8658,  4.6232, 2.4746,
##                 9.7634,  4.6232, 12.2978, 3.8794,
##                 3.2394,  2.4746,  3.8794, 2.4604), 4,4)

Q1s <- -1/DQ[1,1,drop=FALSE]
stopifnot(all(abs(Q1s-c(-0.052237)) < .0001))
Qa2.1 <- matASM(Q1s,DQ[1,2,drop=FALSE],DQ[2,2,drop=FALSE],1)

DQa2.1 <- matrix(c(-0.052237, 0.471993,
                    0.471993, 7.601060),2,2)
stopifnot(all(abs(DQa2.1 - Qa2.1) <.0001))

Q2s <- matSweep(DQa2.1,2)
DQ2s <- matrix(c(-0.081546,  0.062096,
                  0.062096, -0.131561),2,2)
stopifnot(all(abs(DQ2s - Q2s) <.0001))

Qa3.12 <- matASM(Q2s,DQ[1:2,3,drop=FALSE],DQ[3,3,drop=FALSE],1:2)
DQa3.12 <- matrix(c(-0.081546,  0.062096, 0.509084,
                     0.062096, -0.131561, 0.001965,
                     0.509084,  0.001965, 7.318324),3,3)
stopifnot(all(abs(DQa3.12 - Qa3.12) <.0001))

Q3s <- matSweep(DQa3.12,3)
DQ3s <- matrix(c(-0.116959,  0.061959,  0.069563,
                  0.061959, -0.131561,  0.000269,
                  0.069563,  0.000260, -0.136643),3,3)
stopifnot(all(abs(DQ3s - Q3s) <.0001))

Qa4.123 <- matASM(Q3s,DQ[1:3,4,drop=FALSE],DQ[4,4,drop=FALSE],1:3)
DQa4.123 <- matrix(c(-0.116959,  0.061959,  0.069563, -0.044310,
                      0.061959, -0.131562,  0.000269,  0.123809,
                      0.069563,  0.000269, -0.136643,  0.304085,
                     -0.044310,  0.123809,  0.304085,  1.117893),4,4)
stopifnot(all(abs(DQa4.123 - Qa4.123) <.0001))

Q4s <- matSweep(DQa4.123,4)
DQ4s <- matrix(c(-0.118715,  0.066866,  0.081616, -0.039637,
                  0.066866, -0.145274, -0.033409,  0.110752,
                  0.081616, -0.033409, -0.219359,  0.272016,
                 -0.039637,  0.110752,  0.272016, -0.894540),4,4)

stopifnot(all(abs(DQ4s - Q4s) <.00001))

## Trying to do tests with various sizes of Q11 and Q22.

Q2s1 <- matSweep(DQ[1:2,1:2],1)
stopifnot(all(abs(Q2s1-Qa2.1) < .000001))

Q3s1 <- matSweep(DQ[1:3,1:3],1)
Q3a2.1 <- matASM(Q2s1,DQ[1:2,3,drop=FALSE],DQ[3,3,drop=FALSE],1)
stopifnot(all(abs(Q3s1-Q3a2.1)<.000001))

Q4s1 <- matSweep(DQ,1)
Q34a2.1 <- matASM(Q2s1,DQ[1:2,3:4,drop=FALSE],DQ[3:4,3:4,drop=FALSE],1)
stopifnot(all(abs(Q4s1-Q34a2.1)<.000001))

Q4s12 <- matSweep(DQ,1:2)
Q3s12 <- matSweep(DQ[1:3,1:3],1:2)
Q4a3.12 <- matASM(Q3s12,DQ[1:3,4,drop=FALSE],DQ[4,4,drop=FALSE],1:2)
stopifnot(all(abs(Q4s12-Q4a3.12)<.000001))

Q4s13 <- matSweep(DQ,c(1,3))
Q3s13 <- matSweep(DQ[1:3,1:3],c(1,3))
Q4a3.13 <- matASM(Q3s13,DQ[1:3,4,drop=FALSE],DQ[4,4,drop=FALSE],c(1,3))
stopifnot(all(abs(Q4s13-Q4a3.13)<.000001))


##############################
## Test of Multistandardize, follows Dempster (1969, p. 71-3)

DQmst1 <- matMST(DQ,do=1)
Dmst1Q <- matrix(c(1.0000,  2.0651,  2.2314, 0.7404,
                   2.0651, 11.8658,  4.6232, 2.4746,
                   2.2314,  4.6232, 12.2978, 3.8794,
                   0.7404,  2.4746,  3.8794, 2.4604),4,4)
Dmst1C <- diag(c(0.22855,1,1,1))
stopifnot(all(abs(DQmst1$Q-Dmst1Q)<.0002),
          all(abs(DQmst1$C-Dmst1C)<.0002))


DQmst2 <- matMST(DQmst1$Q,DQmst1$C,do=2,done=DQmst1$done)
Dmst2Q <- matrix(c(1.00000, 0.00000,  2.23143, 0.74036,
                   0.00000, 1.00000,  0.00548, 0.34300,
                   2.23143, 0.00548, 12.29780, 3.87940,
                   0.74036, 0.34300,  3.87940, 2.46040),4,4)
Dmst2C <- Dmst1C
Dmst2C[2,1:2] <- c(-0.17119, 0.36270)
stopifnot(all(abs(DQmst2$Q-Dmst2Q)<.0002),
          all(abs(DQmst2$C-Dmst2C)<.0002))

DQmst3 <- matMST(DQmst2$Q,DQmst2$C,do=3,done=DQmst2$done)
Dmst3Q <- matrix(c(1.00000, 0.00000, 0.00000, 0.74036,
                   0.00000, 1.00000, 0.00000, 0.34300,
                   0.00000, 0.00000, 1.00000, 0.82262,
                   0.74036, 0.34300, 0.82262, 2.46040),4,4)
Dmst3C <- Dmst2C
Dmst3C[3,1:3] <- c(-0.18817, -0.00074, 0.36964)
stopifnot(all(abs(DQmst3$Q-Dmst3Q)<.0002),
          all(abs(DQmst3$C-Dmst3C)<.0002))


DQmst4 <- matMST(DQmst3$Q,DQmst3$C,do=4,done=DQmst3$done)
Dmst4C <- Dmst3C
Dmst4C[4,1:4] <- c(0.04190, -0.11709, -0.28759, 0.94581)
stopifnot(all(abs(DQmst4$Q-diag(4))<.0002),
          all(abs(DQmst4$C-Dmst4C)<.0002))

stopifnot(all(abs(DQmst4$C-t(solve(chol(DQ)))) < .00001))

## Checking order for partial assimilations.
DQmst13 <- matMST(DQ,do=c(1,3))
stopifnot(all(abs(DQmst13$Q[c(1,3),c(1,3)]-diag(2))<.0002),
          all(abs(DQmst13$C[-c(1,3),-c(1,3)]-diag(2))<.0002))

## Checking assimilation order
DQmst4r <- matMST(DQ,do=4:1)

## Note:  This is order depedent, as the order determines the new basis.
## Also, C will only be lower triangular, if the elimination is done
## in order.  In particular, if it is done in reverse order, C will be
## upper triangular.

DQmst4rt <- t(solve(chol(DQ[4:1,4:1]))[4:1,4:1])
stopifnot(all(abs(DQmst4r$C-DQmst4rt) <.00001))


