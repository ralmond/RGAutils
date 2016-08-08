### A set of data about the volume of Eggs from a Kitchen experiment
### done by Art Dempster (Dempster, 1968, p 151).

eggs <- data.frame(
    logLength = c(0.7659, 0.7353, 0.7416, 0.7600, 0.7861, 0.7539,
                  0.7747, 0.7718, 0.7889, 0.7659, 0.7689, 0.7478),
    logWidth  = c(0.6360, 0.6198, 0.6280, 0.6280, 0.6239, 0.6156,
                  0.6156, 0.6239, 0.6114, 0.6072, 0.6156, 0.6239),
    logVol  = c(2.031, 1.982, 1.995, 2.019, 2.031, 1.956,
                2.007, 1.995, 1.995, 1.995, 1.995, 2.007))
eggs <- as.matrix(eggs)

## Logs are all log base 10, and logVol is actually log (6/pi) V
