# test.R  (패키지 없이 base R만 사용)

set.seed(123)
n <- 1e6
x <- rnorm(n)

mean_x <- mean(x)
sd_x   <- sd(x)

cat("==== R Test Job ====\n")
cat("N           :", n, "\n")
cat("Mean        :", mean_x, "\n")
cat("SD          :", sd_x, "\n")
cat("Timestamp   :", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

# 결과를 파일로 저장 (Actions에서 artifact로 다운로드 가능)
out <- sprintf("N=%d\nMean=%.6f\nSD=%.6f\nTime=%s\n",
               n, mean_x, sd_x, format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
writeLines(out, "result.txt")

