devtools::load_all()
library(testthat)
library(microbenchmark)
library(lineprof)

## time points are unique
set.seed(1)
n <- 2e6
bw <- c(0.1, 0.05)
xin <- round(matrix(runif(2 * n), n, 2), 12)
yin <- rnorm(n)
win <- rep(1, n)
xout1 <- seq(0, 1, length.out=60)
xout2 <- seq(0, 1, length.out=30)
# tmp <- Lwls2D2(bw, 'epan', xin, yin, win, xout1, xout2, method='plain')
# tmp1 <- Lwls2D2(bw, 'epan', xin, yin, win, xout1, xout2, method='SearchTree')
# expect_equal(tmp, tmp1)
tmp2 <- Lwls2D2(bw, 'epan', xin, yin, win, xout1, xout2, method='sort1')
tmp4 <- Lwls2D2(bw, 'epan', xin, yin, win, xout1, xout2, method='sort2')
expect_equal(tmp2, tmp4)
# tmp3 <- Lwls2D2(bw, 'epan', xin, yin, win, xout1, xout2, method='tree')
# summary(as.numeric(abs(tmp - tmp3)))

# microbenchmark(tmp <-Lwls2D2(c(bw, bw), 'epan', xin, yin, win, xout1, xout2, method='plain'), times=10L)
m2 <- microbenchmark(tmp2 <- Lwls2D2(c(bw, bw), 'epan', xin, yin, win, xout1, xout2, method='sort1'), times=5L)
m4 <- microbenchmark(tmp4 <- Lwls2D2(c(bw, bw), 'epan', xin, yin, win, xout1, xout2, method='sort2'), times=5L)
microbenchmark(tmp3 <- Lwls2D2(c(bw, bw), 'epan', xin, yin, win, xout1, xout2, method='tree'), times=1L)
expect_equal(tmp2, tmp3)

# microbenchmark(tmp1 <- Lwls2D2(c(bw, bw), 'epan', xin, yin, win, xout1, xout2), times=1L)
microbenchmark(tmp2 <- Lwls2D2(c(bw, bw), 'epan', xin, yin, win, xout1, xout2, method='sort1'), times=10L)
expect_equal(tmp, tmp2)

# expect_equal(tmp, tmp1)
# p1 <- lineprof(Lwls2D2(c(bw, bw), 'epan', xin, yin, win, xout1, xout2))
# shine(p1)

## time points are on a finite grid
set.seed(1)
n <- 1e4
bw <- 0.1
xin <- round(matrix(runif(2 * n), n, 2), 3)
yin <- rnorm(n)
win <- rep(1, n)
xout1 <- xout2 <- seq(0, 1, length.out=5)
system.time(matList <- matXinYinWin(xin, yin, win))

tmp <- invisible(with(matList, RmatLwls2d(c(bw, bw), 'epan', yMat, wMat, x1Grid, x2Grid, xout1, xout2, FALSE, FALSE)))
tmp1 <- RmullwlskCC(c(bw, bw), 'epan', t(xin), yin, win, xout1, xout2, FALSE)
test_that('Mat 2D smoother is the same as the long 2D smoother', {
  expect_equal(tmp, t(tmp1))
})

# microbenchmark(with(matList, RmatLwls2d(c(bw, bw), 'epan', yMat, wMat, x1Grid, x2Grid, xout1, xout2, FALSE, FALSE)), times=100L)
# microbenchmark(RmullwlskCC(c(bw, bw), 'epan', t(xin), yin, win, xout1, xout2, FALSE), times=100L)


## Functional setting
set.seed(1)
bw <- 0.1
pts <- seq(0, 1, length.out=100L)
mu <- rep(0, length(pts))
a <- wiener(1000, pts, sparsify=5:10)
xout1 <- xout2 <- seq(0, 1, length.out=20)
rcov <- BinRawCov(GetRawCov(a$yList, a$tList, pts, mu, 'Sparse', 0))
matList <- matXinYinWin(rcov$tPairs, rcov$meanVals, rcov$count)
tmp <- with(matList, RmatLwls2d(c(bw, bw), 'epan', yMat, wMat, x1Grid, x2Grid, xout1, xout2, FALSE, FALSE))
tmp1 <- with(rcov, RmullwlskCC(c(bw, bw), 'epan', t(tPairs), meanVals, count, xout1, xout2, FALSE))
expect_equal(tmp, t(tmp1))

# # speed test
microbenchmark(tmp <- with(matList, RmatLwls2d(c(bw, bw), 'epan', yMat, wMat, x1Grid, x2Grid, xout1, xout2, FALSE, FALSE)), times=10L)
microbenchmark(tmp1 <- with(rcov, RmullwlskCC(c(bw, bw), 'epan', t(tPairs), meanVals, count, xout1, xout2, FALSE)), times=10L)
