##########Part 1##########
#(a)
library(MASS) 
library(survival)

set.seed(1)


n <- 30000  
rho <- 0.5  
gamma <- 1.5  
beta1 <- log(2) 
beta2 <- log(2) 
beta3 <- -1 

# 1. Z1, Z2, W1, W2 생성
mu <- c(0, 0)
sigma <- matrix(c(1, 0.75, 0.75, 1), ncol = 2)
Z <- mvrnorm(n, mu = mu, Sigma = sigma)
Z1 <- Z[, 1]
Z2 <- Z[, 2]

# W1은 N(0,1)에서 생성, W2는 Bernoulli(0.5)에서 생성
W1 <- rnorm(n, mean = 0, sd = 1)
W2 <- rbinom(n, 1, 0.5)

# 2. u 생성 (uniform(0, 1)에서)
u <- runif(n)

# 3. Weibull 분포의 생존함수 S(t)를 통해 t 생성
# u = 1 - S(t) = 1 - exp(-exp(β1*Z1 + β2*W1 + β3*W2) * rho * t^gamma)
# 이를 t에 대해 풀면 아래와 같이 전개 가능
linear_predictor <- beta1 * Z1 + beta2 * W1 + beta3 * W2
T <- (-log(1 - u) / (rho * exp(linear_predictor)))^(1 / gamma)

status <- rep(1, n)

data <- data.frame(T = T, status = status, Z1 = Z1, W1 = W1, W2 = W2)

# 4. Cox 비례위험 모형 적합
cox_model <- coxph(Surv(T, status) ~ Z1 + W1 + W2, data = data)
summary(cox_model)

#(b)
set.seed(1)

n <- 30000  
rho <- 0.5  
gamma <- 1.5  
beta1 <- log(2) 
beta2 <- log(2) 
beta3 <- -1 

nsim <- 100

beta1_estimates <- numeric(nsim)
beta2_estimates <- numeric(nsim)
beta3_estimates <- numeric(nsim)
se_beta1_estimates <- numeric(nsim)
se_beta2_estimates <- numeric(nsim)
se_beta3_estimates <- numeric(nsim)

# 100번 반복하여 각 추정값 저장
for (i in 1:nsim) {
  
  mu <- c(0, 0)
  sigma <- matrix(c(1, 0.75, 0.75, 1), ncol = 2)
  Z <- mvrnorm(n, mu = mu, Sigma = sigma)
  Z1 <- Z[, 1]
  Z2 <- Z[, 2]
  W1 <- rnorm(n, mean = 0, sd = 1)
  W2 <- rbinom(n, 1, 0.5)
  
  u <- runif(n)
  
  linear_predictor <- beta1 * Z1 + beta2 * W1 + beta3 * W2
  T <- (-log(1 - u) / (rho * exp(linear_predictor)))^(1 / gamma)
  status <- rep(1, n)
  data <- data.frame(T = T, status = status, Z1 = Z1, W1 = W1, W2 = W2)
  cox_model <- coxph(Surv(T, status) ~ Z1 + W1 + W2, data = data)
  
  beta1_estimates[i] <- cox_model$coef[1]
  beta2_estimates[i] <- cox_model$coef[2]
  beta3_estimates[i] <- cox_model$coef[3]
  se_beta1_estimates[i] <- summary(cox_model)$coefficients[1, 3]
  se_beta2_estimates[i] <- summary(cox_model)$coefficients[2, 3]
  se_beta3_estimates[i] <- summary(cox_model)$coefficients[3, 3]
}

# 추정된 β 값들의 표본 평균 계산
mean_beta1 <- mean(beta1_estimates)
mean_beta2 <- mean(beta2_estimates)
mean_beta3 <- mean(beta3_estimates)

# 추정된 표준 오차들의 표본 평균 계산
mean_se_beta1 <- mean(se_beta1_estimates)
mean_se_beta2 <- mean(se_beta2_estimates)
mean_se_beta3 <- mean(se_beta3_estimates)

# 추정된 β 값들의 표본 표준 편차 계산
sd_beta1 <- sd(beta1_estimates)
sd_beta2 <- sd(beta2_estimates)
sd_beta3 <- sd(beta3_estimates)

# 결과 출력
cat("Sample means of β estimates:\n")
cat("Mean β1:", mean_beta1, "\n")
cat("Mean β2:", mean_beta2, "\n")
cat("Mean β3:", mean_beta3, "\n\n")

cat("Sample means of standard error estimates:\n")
cat("Mean SE of β1:", mean_se_beta1, "\n")
cat("Mean SE of β2:", mean_se_beta2, "\n")
cat("Mean SE of β3:", mean_se_beta3, "\n\n")

cat("Sample standard deviations of β estimates:\n")
cat("SD of β1:", sd_beta1, "\n")
cat("SD of β2:", sd_beta2, "\n")
cat("SD of β3:", sd_beta3, "\n\n")

# 비교: 표준 오차의 표본 평균과 β 추정치의 표본 표준 편차
cat("Comparison of SE means and SD of estimates:\n")
cat("β1: Mean SE =", mean_se_beta1, ", SD =", sd_beta1, "\n")
cat("β2: Mean SE =", mean_se_beta2, ", SD =", sd_beta2, "\n")
cat("β3: Mean SE =", mean_se_beta3, ", SD =", sd_beta3, "\n")

#(c)
set.seed(1)

n_large <- 1e6  # 매우 큰 샘플 크기 (100만)
rho <- 0.5  
gamma <- 1.5 
beta1 <- log(2) 
beta2 <- log(2) 
beta3 <- -1  


mu <- c(0, 0)
sigma <- matrix(c(1, 0.75, 0.75, 1), ncol = 2)
Z <- mvrnorm(n_large, mu = mu, Sigma = sigma)
Z1 <- Z[, 1]
Z2 <- Z[, 2]
W1 <- rnorm(n_large, mean = 0, sd = 1)
W2 <- rbinom(n_large, 1, 0.5)

u <- runif(n_large)

linear_predictor <- beta1 * Z1 + beta2 * W1 + beta3 * W2
T <- (-log(1 - u) / (rho * exp(linear_predictor)))^(1 / gamma)

# 주어진 censoring rate에 맞는 지수 분포 매개변수를 단계적으로 찾는 함수
find_lambda_stepwise <- function(censoring_rate, T_values, lambda_values) {
  best_lambda <- lambda_values[1]
  best_diff <- Inf
  
  for (lambda in lambda_values) {
    C <- rexp(length(T_values), rate = lambda)  # 지수 분포에서 C 생성
    censoring_actual <- mean(T_values > C) 
    
    diff <- abs(censoring_actual - censoring_rate)
    
    # 가장 작은 차이를 가진 lambda 값을 기록
    if (diff < best_diff) {
      best_lambda <- lambda
      best_diff <- diff
    }
  }
  
  return(best_lambda)  # 최적의 lambda 값 반환
}


lambda_values <- seq(0.00001, 20, length.out = 1000)  

censoring_rates <- c(0.10, 0.30, 0.90, 0.95, 0.99)

lambdas <- sapply(censoring_rates, 
                  function(cr) find_lambda_stepwise(cr, T, lambda_values))

cat("Censoring rates and corresponding lambda values:\n")
for (i in 1:length(censoring_rates)) {
  cat("Censoring rate:", censoring_rates[i], " -> Lambda:", lambdas[i], "\n")
}

#(d)
set.seed(1)

n <- 30000  
rho <- 0.5 
gamma <- 1.5  
beta1 <- log(2) 
beta2 <- log(2)  
beta3 <- -1  

# lambda 값 (c에서 계산한 결과 사용)
lambdas <- c(0.04005002, 0.1801901, 2.922931, 5.205213, 17.0971)

nsim <- 100

# 결과 저장용 리스트 생성
results <- list()

# censoring rate별로 반복
for (lambda in lambdas) {
  
  # 추정된 계수와 표준 오차를 저장할 공간 생성
  beta1_estimates <- numeric(nsim)
  beta2_estimates <- numeric(nsim)
  beta3_estimates <- numeric(nsim)
  se_beta1_estimates <- numeric(nsim)
  se_beta2_estimates <- numeric(nsim)
  se_beta3_estimates <- numeric(nsim)
  
  for (i in 1:nsim) {
    mu <- c(0, 0)
    sigma <- matrix(c(1, 0.75, 0.75, 1), ncol = 2)
    Z <- mvrnorm(n, mu = mu, Sigma = sigma)
    Z1 <- Z[, 1]
    Z2 <- Z[, 2]
    W1 <- rnorm(n, mean = 0, sd = 1)
    W2 <- rbinom(n, 1, 0.5)
    
    u <- runif(n)
    
    linear_predictor <- beta1 * Z1 + beta2 * W1 + beta3 * W2
    T <- (-log(1 - u) / (rho * exp(linear_predictor)))^(1 / gamma)
    
    C <- rexp(n, rate = lambda)
    
    # X = min(T, C)와 Δ = I(T < C) 계산
    X <- pmin(T, C)
    Delta <- as.numeric(T < C)
    
    data <- data.frame(X = X, Delta = Delta, Z1 = Z1, W1 = W1, W2 = W2)
    
    cox_model <- coxph(Surv(X, Delta) ~ Z1 + W1 + W2, data = data)
    
    # 추정된 계수와 표준 오차 저장
    beta1_estimates[i] <- cox_model$coef[1]
    beta2_estimates[i] <- cox_model$coef[2]
    beta3_estimates[i] <- cox_model$coef[3]
    se_beta1_estimates[i] <- summary(cox_model)$coefficients[1, 3]
    se_beta2_estimates[i] <- summary(cox_model)$coefficients[2, 3]
    se_beta3_estimates[i] <- summary(cox_model)$coefficients[3, 3]
  }
  
  # censoring rate별 결과 저장
  results[[paste0("Lambda_", lambda)]] <- list(
    beta1 = beta1_estimates,
    beta2 = beta2_estimates,
    beta3 = beta3_estimates,
    se_beta1 = se_beta1_estimates,
    se_beta2 = se_beta2_estimates,
    se_beta3 = se_beta3_estimates
  )
}

for (name in names(results)) {
  cat("Results for", name, ":\n")
  
  beta1_mean <- mean(results[[name]]$beta1)
  beta2_mean <- mean(results[[name]]$beta2)
  beta3_mean <- mean(results[[name]]$beta3)
  se_beta1_mean <- mean(results[[name]]$se_beta1)
  se_beta2_mean <- mean(results[[name]]$se_beta2)
  se_beta3_mean <- mean(results[[name]]$se_beta3)
  
  cat("Mean beta1:", beta1_mean, "\n")
  cat("Mean beta2:", beta2_mean, "\n")
  cat("Mean beta3:", beta3_mean, "\n")
  cat("Mean SE beta1:", se_beta1_mean, "\n")
  cat("Mean SE beta2:", se_beta2_mean, "\n")
  cat("Mean SE beta3:", se_beta3_mean, "\n\n")
}

##########Part 2##########
#(a)
#i.
set.seed(1)

generate_data <- function(n, beta1, beta2, beta3, rho, gamma, lambda) {
  library(MASS)
  Sigma <- matrix(c(1, 0.75, 0.75, 1), 2, 2)
  Z <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  Z1 <- Z[, 1]
  Z2 <- Z[, 2]
  
  W1 <- rnorm(n, 0, 1)     
  W2 <- rbinom(n, 1, 0.5)    
  
  U <- runif(n)
  linear_pred <- beta1 * Z1 + beta2 * W1 + beta3 * W2
  T <- (-log(1 - U) / (lambda * exp(linear_pred)))^(1 / gamma)
  
  C <- rexp(n, rate = lambda)  
  
  X <- pmin(T, C)   #observed time     
  Delta <- as.numeric(T <= C) #status
  
  data <- data.frame(X = X, Delta = Delta, Z1 = Z1, Z2 = Z2, W1 = W1, W2 = W2)
  return(data)
}

case_cohort_sample <- function(data, subcohort_size) {
  # 무작위 하위 코호트 추출
  subcohort <- data[sample(nrow(data), subcohort_size, replace = FALSE), ]
  
  # 하위 코호트에 포함되지 않은 사건 데이터 추가
  remaining_cases <- data[data$Delta == 1 & !(rownames(data) %in% rownames(subcohort)), ]
  
  # 최종 case-cohort
  case_cohort <- rbind(subcohort, remaining_cases)
  return(case_cohort)
}

n <- 30000  
rho <- 0.5  
gamma <- 1.5  
beta1 <- log(2) 
beta2 <- log(2) 
beta3 <- -1 
lambda <- 17.0971 # 99% censoring rate

data <- generate_data(n = 30000, beta1, beta2, beta3, rho, gamma, lambda)

subcohort_size <- 100
sample_data <- case_cohort_sample(data, subcohort_size)

#ii.
num_failures <- sum(sample_data$Delta == 1)
total_sample_size <- nrow(sample_data)
cat("Number of failures:", num_failures, "\n")
cat("Case-cohort sample size:", total_sample_size, "\n")

#iii.
library(survey)

nc <- sum(data$Delta == 0)
nc_sub <- sum(sample_data$Delta == 0)

weights <- ifelse(sample_data$Delta == 1, 1, nc / nc_sub)
design <- svydesign(ids = ~1, data = sample_data, weights = ~weights)
cox_model <- svycoxph(Surv(X, Delta) ~ Z1 + W1 + W2, design)
summary(cox_model)

#iv.
nsim <- 500
results_case_cohort <- list()

for (i in 1:nsim) {
  data <- generate_data(n = 30000, beta1, beta2, beta3, rho, gamma, lambda)
  sample_data <- case_cohort_sample(data, subcohort_size)
  
  weights <- ifelse(sample_data$Delta == 1, 1, nc / nc_sub)
  design <- svydesign(ids = ~1, data = sample_data, weights = ~weights)
  cox_model <- svycoxph(Surv(X, Delta) ~ Z1 + W1 + W2, design)
  
  results_case_cohort[[i]] <- c(cox_model$coef[1], cox_model$coef[2], cox_model$coef[3])
}

#v.
nsim <- 500
subcohort_size <- 100
results_case_cohort <- list()

for (i in 1:nsim) {
  print(i)
  data <- generate_data(n = 30000, beta1, beta2, beta3, rho, gamma, lambda)
  sample_data <- case_cohort_sample(data, subcohort_size)
  
  weights <- ifelse(sample_data$Delta == 1, 1, nc / nc_sub)
  design <- svydesign(ids = ~1, data = sample_data, weights = ~weights)
  cox_model <- svycoxph(Surv(X, Delta) ~ Z1 + W1 + W2, design)
  
  num_failures <- sum(sample_data$Delta == 1)
  total_sample_size <- nrow(sample_data)
  
  results_case_cohort[[i]] <- list(coef = c(cox_model$coef[1], cox_model$coef[2], cox_model$coef[3]),
                                   num_failures = num_failures, total_sample_size = total_sample_size)
}

failure_counts <- sapply(results_case_cohort, function(x) x$num_failures)
sample_sizes <- sapply(results_case_cohort, function(x) x$total_sample_size)

cat("Average number of failures:", mean(failure_counts), "\n")
cat("Average case-cohort sample size:", mean(sample_sizes), "\n")

#vi.
beta1_estimates <- sapply(results_case_cohort, function(x) x$coef[1])
beta2_estimates <- sapply(results_case_cohort, function(x) x$coef[2])
beta3_estimates <- sapply(results_case_cohort, function(x) x$coef[3])


mean_beta1 <- mean(beta1_estimates)
mean_beta2 <- mean(beta2_estimates)
mean_beta3 <- mean(beta3_estimates)

sd_beta1 <- sd(beta1_estimates)
sd_beta2 <- sd(beta2_estimates)
sd_beta3 <- sd(beta3_estimates)

cat("Mean beta1 estimate:", mean_beta1, "\n")
cat("Mean beta2 estimate:", mean_beta2, "\n")
cat("Mean beta3 estimate:", mean_beta3, "\n")
cat("Standard deviation of beta1 estimates:", sd_beta1, "\n")
cat("Standard deviation of beta2 estimates:", sd_beta2, "\n")
cat("Standard deviation of beta3 estimates:", sd_beta3, "\n")

#vii.
subcohort_size <- 300
results_case_cohort_large <- list()

for (i in 1:nsim) {
  data <- generate_data(n = 30000, beta1, beta2, beta3, rho, gamma, lambda)
  sample_data <- case_cohort_sample(data, subcohort_size)
  
  weights <- ifelse(sample_data$Delta == 1, 1, nc / nc_sub)
  design <- svydesign(ids = ~1, data = sample_data, weights = ~weights)
  cox_model <- svycoxph(Surv(X, Delta) ~ Z1 + W1 + W2, design)
  
  results_case_cohort_large[[i]] <- c(cox_model$coef[1], 
                                      cox_model$coef[2], cox_model$coef[3])
}

beta1_estimates_large <- sapply(results_case_cohort_large, function(x) x[1])
beta2_estimates_large <- sapply(results_case_cohort_large, function(x) x[2])
beta3_estimates_large <- sapply(results_case_cohort_large, function(x) x[3])

mean_beta1_large <- mean(beta1_estimates_large)
mean_beta2_large <- mean(beta2_estimates_large)
mean_beta3_large <- mean(beta3_estimates_large)

sd_beta1_large <- sd(beta1_estimates_large)
sd_beta2_large <- sd(beta2_estimates_large)
sd_beta3_large <- sd(beta3_estimates_large)

cat("Mean beta1 estimate with subcohort size 300:", mean_beta1_large, "\n")
cat("Mean beta2 estimate with subcohort size 300:", mean_beta2_large, "\n")
cat("Mean beta3 estimate with subcohort size 300:", mean_beta3_large, "\n")
cat("Standard deviation of beta1 estimates with subcohort size 300:", sd_beta1_large, "\n")
cat("Standard deviation of beta2 estimates with subcohort size 300:", sd_beta2_large, "\n")
cat("Standard deviation of beta3 estimates with subcohort size 300:", sd_beta3_large, "\n")

#(b)
#i.
set.seed(1)

nested_case_control_sample <- function(data, controls_per_case = 1) {
  cases <- data[data$Delta == 1, ]  # 사건 데이터 선택
  controls <- data.frame()          # 대조군 데이터를 저장할 곳
  
  for (i in 1:nrow(cases)) {
    # 비사건 데이터 선택
    control_candidates <- data[data$Delta == 0 & data$X > cases$X[i], ]
    
    # 대조군을 무작위로 샘플링
    if (nrow(control_candidates) > 0) {
      sampled_controls <- control_candidates[sample(nrow(control_candidates), 
                                                    min(controls_per_case, nrow(control_candidates)), replace = FALSE), ]
      controls <- rbind(controls, sampled_controls)
    }
  }
  
  nested_sample <- rbind(cases, controls)
  return(nested_sample)
}

data <- generate_data(n = 30000, beta1, beta2, beta3, rho, gamma, lambda)
sample_data <- nested_case_control_sample(data, controls_per_case = 1)

weights <- ifelse(sample_data$Delta == 1, 1, 1)
design <- svydesign(ids = ~1, data = sample_data, weights = ~weights)
cox_model <- svycoxph(Surv(X, Delta) ~ Z1 + W1 + W2, design)

num_failures <- sum(sample_data$Delta == 1)
total_sample_size <- nrow(sample_data)
cat("Number of failures:", num_failures, "\n")
cat("Case-control sample size:", total_sample_size, "\n")
cat("Estimated beta1:", cox_model$coef[1], "\n")
cat("Estimated beta2:", cox_model$coef[2], "\n")
cat("Estimated beta3:", cox_model$coef[3], "\n")


###500번 반복
set.seed(1)

nsim <- 500
results_nested_case_control <- list()

for (i in 1:nsim) {
  data <- generate_data(n = 30000, beta1, beta2, beta3, rho, gamma, lambda)
  sample_data <- nested_case_control_sample(data, controls_per_case = 1)
  
  weights <- ifelse(sample_data$Delta == 1, 1, 1)
  design <- svydesign(ids = ~1, data = sample_data, weights = ~weights)
  cox_model <- svycoxph(Surv(X, Delta) ~ Z1 + W1 + W2, design)
  
  num_failures <- sum(sample_data$Delta == 1)
  total_sample_size <- nrow(sample_data)
  
  results_nested_case_control[[i]] <- list(coef = c(cox_model$coef[1], cox_model$coef[2], cox_model$coef[3]),
                                           num_failures = num_failures, total_sample_size = total_sample_size)
}

failure_counts <- sapply(results_nested_case_control, function(x) x$num_failures)
sample_sizes <- sapply(results_nested_case_control, function(x) x$total_sample_size)

cat("Average number of failures (nested case-control):", mean(failure_counts), "\n")
cat("Average case-control sample size (nested case-control):", mean(sample_sizes), "\n")

beta1_estimates_nested <- sapply(results_nested_case_control, function(x) x$coef[1])
beta2_estimates_nested <- sapply(results_nested_case_control, function(x) x$coef[2])
beta3_estimates_nested <- sapply(results_nested_case_control, function(x) x$coef[3])

cat("Mean beta1 estimate (nested case-control):", mean(beta1_estimates_nested), "\n")
cat("Mean beta2 estimate (nested case-control):", mean(beta2_estimates_nested), "\n")
cat("Mean beta3 estimate (nested case-control):", mean(beta3_estimates_nested), "\n")

cat("Standard deviation of beta1 estimates (nested case-control):", sd(beta1_estimates_nested), "\n")
cat("Standard deviation of beta2 estimates (nested case-control):", sd(beta2_estimates_nested), "\n")
cat("Standard deviation of beta3 estimates (nested case-control):", sd(beta3_estimates_nested), "\n")

#ii.
set.seed(1)

data <- generate_data(n = 30000, beta1, beta2, beta3, rho, gamma, lambda)
sample_data <- nested_case_control_sample(data, controls_per_case = 5)

weights <- ifelse(sample_data$Delta == 1, 1, 1)
design <- svydesign(ids = ~1, data = sample_data, weights = ~weights)
cox_model <- svycoxph(Surv(X, Delta) ~ Z1 + W1 + W2, design)

num_failures <- sum(sample_data$Delta == 1)
total_sample_size <- nrow(sample_data)
cat("Number of failures:", num_failures, "\n")
cat("Case-control sample size:", total_sample_size, "\n")
cat("Estimated beta1:", cox_model$coef[1], "\n")
cat("Estimated beta2:", cox_model$coef[2], "\n")
cat("Estimated beta3:", cox_model$coef[3], "\n")

###500번 반복
set.seed(1)
nsim <- 500
results_nested_case_control_large <- list()

for (i in 1:nsim) {
  data <- generate_data(n = 30000, beta1, beta2, beta3, rho, gamma, lambda)
  sample_data <- nested_case_control_sample(data, controls_per_case = 5)
  
  weights <- ifelse(sample_data$Delta == 1, 1, 1)
  design <- svydesign(ids = ~1, data = sample_data, weights = ~weights)
  cox_model <- svycoxph(Surv(X, Delta) ~ Z1 + W1 + W2, design)
  
  num_failures <- sum(sample_data$Delta == 1)
  total_sample_size <- nrow(sample_data)
  
  results_nested_case_control_large[[i]] <- list(coef = c(cox_model$coef[1], cox_model$coef[2], cox_model$coef[3]),
                                                 num_failures = num_failures, total_sample_size = total_sample_size)
}

failure_counts_large <- sapply(results_nested_case_control_large, function(x) x$num_failures)
sample_sizes_large <- sapply(results_nested_case_control_large, function(x) x$total_sample_size)

cat("Average number of failures (nested case-control, 5 controls):", mean(failure_counts_large), "\n")
cat("Average case-control sample size (nested case-control, 5 controls):", mean(sample_sizes_large), "\n")

beta1_estimates_nested_large <- sapply(results_nested_case_control_large, function(x) x$coef[1])
beta2_estimates_nested_large <- sapply(results_nested_case_control_large, function(x) x$coef[2])
beta3_estimates_nested_large <- sapply(results_nested_case_control_large, function(x) x$coef[3])

cat("Mean beta1 estimate (nested case-control, 5 controls):", mean(beta1_estimates_nested_large), "\n")
cat("Mean beta2 estimate (nested case-control, 5 controls):", mean(beta2_estimates_nested_large), "\n")
cat("Mean beta3 estimate (nested case-control, 5 controls):", mean(beta3_estimates_nested_large), "\n")

cat("Standard deviation of beta1 estimates (nested case-control, 5 controls):", sd(beta1_estimates_nested_large), "\n")
cat("Standard deviation of beta2 estimates (nested case-control, 5 controls):", sd(beta2_estimates_nested_large), "\n")
cat("Standard deviation of beta3 estimates (nested case-control, 5 controls):", sd(beta3_estimates_nested_large), "\n")


##########Part 3##########
#(c)
library(survey)

set.seed(1)
# 전체 데이터는 Z1 생성 후 제거
generate_data <- function(n, beta1, beta2, beta3, rho, gamma, lambda) {
  library(MASS)
  Sigma <- matrix(c(1, 0.75, 0.75, 1), 2, 2)
  Z <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  Z1 <- Z[, 1]
  Z2 <- Z[, 2]
  
  W1 <- rnorm(n, 0, 1)     
  W2 <- rbinom(n, 1, 0.5)    
  
  U <- runif(n)
  linear_pred <- beta1 * Z1 + beta2 * W1 + beta3 * W2
  T <- (-log(1 - U) / (lambda * exp(linear_pred)))^(1 / gamma)
  
  C <- rexp(n, rate = lambda)  
  
  X <- pmin(T, C)   # Observed time     
  Delta <- as.numeric(T <= C)  # Status
  
  # Z1을 제거한 전체 데이터 반환
  data <- data.frame(X = X, Delta = Delta, Z2 = Z2, W1 = W1, W2 = W2)
  return(list(data = data, Z1 = Z1))  # Z1은 별도로 반환
}

# Case-cohort 샘플에만 Z1 추가
case_cohort_sample <- function(data, Z1, subcohort_size) {
  subcohort_indices <- sample(nrow(data), subcohort_size, replace = FALSE)
  subcohort <- data[subcohort_indices, ]
  subcohort$Z1 <- Z1[subcohort_indices]  # Z1 추가
  
  remaining_cases_indices <- which(data$Delta == 1 & !(1:nrow(data) %in% subcohort_indices))
  remaining_cases <- data[remaining_cases_indices, ]
  remaining_cases$Z1 <- Z1[remaining_cases_indices]  # Z1 추가
  
  # Case-cohort 데이터 생성
  case_cohort <- rbind(subcohort, remaining_cases)
  return(case_cohort)
}

# 파라미터 설정
n <- 30000  
rho <- 0.5  
gamma <- 1.5  
beta1 <- log(2) 
beta2 <- log(2) 
beta3 <- -1 
lambda <- 17.0971  # 99% censoring rate

# 데이터 생성 및 Z1 제거
generated <- generate_data(n = n, beta1, beta2, beta3, rho, gamma, lambda)
data <- generated$data  # Z1이 없는 전체 데이터
Z1 <- generated$Z1      # Z1은 따로 보관

# Case-cohort 샘플 생성
subcohort_size <- 300
sample_data <- case_cohort_sample(data, Z1, subcohort_size)

# Case-cohort 샘플에서 비사건 수 계산
nc <- sum(data$Delta == 0)
nc_sub <- sum(sample_data$Delta == 0)

# 시뮬레이션 수행
nsim <- 3000
results_case_cohort <- list()

for (i in 1:nsim) {
  print(i)
  
  # 데이터 생성
  generated <- generate_data(n = n, beta1, beta2, beta3, rho, gamma, lambda)
  data <- generated$data
  Z1 <- generated$Z1
  
  # Case-cohort 샘플 생성
  sample_data <- case_cohort_sample(data, Z1, subcohort_size)
  
  # 가중치 계산
  weights <- ifelse(sample_data$Delta == 1, 1, nc / nc_sub)
  
  # 설계 객체 및 Cox 모델 적합
  design <- svydesign(ids = ~1, data = sample_data, weights = ~weights)
  cox_model <- svycoxph(Surv(X, Delta) ~ Z1 + W1 + W2, design)
  
  # 결과 저장
  num_failures <- sum(sample_data$Delta == 1)
  total_sample_size <- nrow(sample_data)
  
  results_case_cohort[[i]] <- list(coef = c(cox_model$coef[1], cox_model$coef[2], cox_model$coef[3]),
                                   num_failures = num_failures, total_sample_size = total_sample_size)
}

# 결과 요약
beta1_estimates <- sapply(results_case_cohort, function(x) x$coef[1])
beta2_estimates <- sapply(results_case_cohort, function(x) x$coef[2])
beta3_estimates <- sapply(results_case_cohort, function(x) x$coef[3])

mean_beta1 <- mean(beta1_estimates)
mean_beta2 <- mean(beta2_estimates)
mean_beta3 <- mean(beta3_estimates)

sd_beta1 <- sd(beta1_estimates)
sd_beta2 <- sd(beta2_estimates)
sd_beta3 <- sd(beta3_estimates)

cat("Mean beta1 estimate:", mean_beta1, "\n")
cat("Mean beta2 estimate:", mean_beta2, "\n")
cat("Mean beta3 estimate:", mean_beta3, "\n")
cat("Standard deviation of beta1 estimates:", sd_beta1, "\n")
cat("Standard deviation of beta2 estimates:", sd_beta2, "\n")
cat("Standard deviation of beta3 estimates:", sd_beta3, "\n")



#Mean beta1 estimate: 0.6989477 
#Mean beta2 estimate: 0.7008559 
#Mean beta3 estimate: -1.006281 
#Standard deviation of beta1 estimates: 0.04714039 
#Standard deviation of beta2 estimates: 0.04710016 
#Standard deviation of beta3 estimates: 0.09731998 


############################
#dfbeta
library(survey)
library(survival)
library(MASS)

# 데이터 생성 함수
generate_data <- function(n, beta1, beta2, beta3, rho, gamma, lambda) {
  Sigma <- matrix(c(1, 0.75, 0.75, 1), 2, 2)
  Z <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  Z1 <- Z[, 1]
  Z2 <- Z[, 2]
  
  W1 <- rnorm(n, 0, 1)     
  W2 <- rbinom(n, 1, 0.5)    
  
  U <- runif(n)
  linear_pred <- beta1 * Z1 + beta2 * W1 + beta3 * W2
  T <- (-log(1 - U) / (lambda * exp(linear_pred)))^(1 / gamma)
  
  C <- rexp(n, rate = lambda)  
  
  X <- pmin(T, C)   # Observed time     
  Delta <- as.numeric(T <= C)  # Status
  
  # Z1을 제거한 전체 데이터 반환
  data <- data.frame(X = X, Delta = Delta, Z2 = Z2, W1 = W1, W2 = W2)
  return(list(data = data, Z1 = Z1))  # Z1은 별도로 반환
}

# Case-cohort 샘플링 함수
case_cohort_sample <- function(data, Z1, subcohort_size) {
  subcohort_indices <- sample(nrow(data), subcohort_size, replace = FALSE)
  subcohort <- data[subcohort_indices, ]
  subcohort$Z1 <- Z1[subcohort_indices]  # Subcohort에 Z1 추가
  
  remaining_cases_indices <- which(data$Delta == 1 & !(1:nrow(data) %in% subcohort_indices))
  remaining_cases <- data[remaining_cases_indices, ]
  remaining_cases$Z1 <- Z1[remaining_cases_indices]  # 사건 데이터에 Z1 추가
  
  # 최종 Case-cohort 데이터 생성
  case_cohort <- rbind(subcohort, remaining_cases)
  return(case_cohort)
}

# 반복 설정
nsim <- 3000
n <- 300
results <- list()

# 시뮬레이션
for (i in 1:nsim) {
  print(paste("Iteration:", i))
  
  # 1. 데이터 생성
  data_gen <- generate_data(n = 30000, beta1 = log(2), beta2 = log(2), beta3 = -1, rho = 0.75, gamma = 1.5, lambda = 17.0971)
  Z1 <- data_gen$Z1  # 실제 Z1 값
  data <- data_gen$data  # 전체 데이터
  
  # 2. 하위 코호트 샘플링
  sampled_data <- case_cohort_sample(data, Z1, subcohort_size = n)
  
  # 3. 설계 객체 생성 및 Z1 예측
  sampled_data$weight <- 1
  twophase_design <- svydesign(ids = ~0, weights = ~weight, data = sampled_data)
  est_svyglm <- svyglm(Z1 ~ Z2 + W1 + W2, design = twophase_design)
  data$Z1_est <- predict(est_svyglm, newdata = data, type = "response")  # Z1 예측값
  
  # 4. Z1_final 생성 (실제 값은 그대로, 없는 경우 예측값 사용)
  data$Z1_final <- ifelse(is.na(Z1), data$Z1_est, Z1)  # Z1 존재 여부에 따라 값 선택
  
  # 5. Cox 모델 적합 및 dfbeta 계산 (Z1에만 적용)
  cox_model <- coxph(Surv(X, Delta) ~ Z1_final + W1 + W2, data = data)
  dfbeta_values <- residuals(cox_model, type = "dfbeta")
  data$dfbeta_Z1 <- dfbeta_values[, 1] + 1  # Z1에 대해서만 dfbeta 사용, 안정성을 위해 1 추가
  
  # 6. 보정된 가중치 계산 (Z1만 사용, 총합을 모집단 크기로 조정)
  pop_total <- sum(data$dfbeta_Z1, na.rm = TRUE)
  data$calibrated_weights <- (data$dfbeta_Z1 / pop_total) * 30000  # 총합을 30000으로 맞춤
  data$calibrated_weights[data$calibrated_weights <= 0] <- 1e-6  # 음수나 0 방지
  
  # 7. 보정된 Cox 모델 적합 (Z1_final만 가중치 반영)
  cox_model_calibrated <- coxph(Surv(X, Delta) ~ Z1_final + W1 + W2, data = data, weights = calibrated_weights)
  
  # 8. 결과 저장
  results[[i]] <- cox_model_calibrated$coef
}

# 결과 요약
beta1_estimates <- sapply(results, function(x) x["Z1_final"])  # Beta1 (Z1)
beta2_estimates <- sapply(results, function(x) x["W1"])        # Beta2 (W1)
beta3_estimates <- sapply(results, function(x) x["W2"])        # Beta3 (W2)

# 각 계수의 평균 계산
mean_beta1 <- mean(beta1_estimates)
mean_beta2 <- mean(beta2_estimates)
mean_beta3 <- mean(beta3_estimates)

# 각 계수의 표준 편차 계산
sd_beta1 <- sd(beta1_estimates)
sd_beta2 <- sd(beta2_estimates)
sd_beta3 <- sd(beta3_estimates)

# 결과 출력
cat("Mean beta1 estimate:", mean_beta1, "\n")
cat("Mean beta2 estimate:", mean_beta2, "\n")
cat("Mean beta3 estimate:", mean_beta3, "\n")
cat("Standard deviation of beta1 estimates:", sd_beta1, "\n")
cat("Standard deviation of beta2 estimates:", sd_beta2, "\n")
cat("Standard deviation of beta3 estimates:", sd_beta3, "\n")


#Mean beta1 estimate: 0.6935672 
#Mean beta2 estimate: 0.6936149 
#Mean beta3 estimate: -1.000317 
#Standard deviation of beta1 estimates: 0.01439908 
#Standard deviation of beta2 estimates: 0.01438206 
#Standard deviation of beta3 estimates: 0.02915825 

#(d)
# 필요한 패키지 로드
library(survey)
library(survival)
library(MASS)
library(tidyverse)

# 데이터 불러오기
data <- read.csv("C:/Data/mimic3_final.csv")
data<-data%>%filter(Glascow.coma.scale.eye.opening<=3)%>%filter(Glascow.coma.scale.verbal.response<=4)
#Glascow.coma.scale.eye.opening은 1 to 4, Glascow.coma.scale.verbal.response은 1 to 5라고 자료에 기재되어있음)
full_data<-data

# Glucose 데이터를 조작하여 subcohort에만 존재하게 변경
case_cohort_sample <- function(data, subcohort_size) {
  subcohort_indices <- sample(nrow(data), subcohort_size, replace = FALSE)
  subcohort <- data[subcohort_indices, ]
  
  remaining_cases_indices <- which(data$delta == 1 & !(1:nrow(data) %in% subcohort_indices))
  remaining_cases <- data[remaining_cases_indices, ]
  
  case_cohort_data <- rbind(subcohort, remaining_cases)
  
  # Subcohort 외에는 Glucose 값을 NA로 처리
  data$Glucose <- ifelse(1:nrow(data) %in% subcohort_indices, data$Glucose, NA)
  
  return(list(case_cohort_data = case_cohort_data, full_data = data))
}

# Glucose 값을 모든 데이터에서 사용할 수 있다고 가정 (Golden Standard), 보조변수인 Mean pressure 제외
cox_model <- coxph(Surv(futime, delta) ~ Glucose + Heart.Rate + Height +
                     Oxygen.saturation + Respiratory.rate + Temperature + Weight +
                     Glascow.coma.scale.eye.opening + Glascow.coma.scale.verbal.response, 
                   data = data)

summary(cox_model) #golden standard
golden_coef <- summary(cox_model)$coef

# Case-cohort 방식 분석
case_cohort_analysis <- function(case_cohort_data, full_data) {
  # Glucose 예측 모델 생성
  case_cohort_data$weight <- ifelse(case_cohort_data$delta == 1, 1, nrow(full_data) / nrow(case_cohort_data))
  design <- svydesign(ids = ~1, data = case_cohort_data, weights = ~weight)
  svy_model <- svyglm(Glucose ~ Mean.blood.pressure + Heart.Rate + Height + Oxygen.saturation +
                        Respiratory.rate + Temperature + Weight + Glascow.coma.scale.eye.opening +
                        Glascow.coma.scale.verbal.response, 
                      design = design)
  
  # Full data에서 Glucose 예측값 생성
  full_data$Glucose_pred <- predict(svy_model, newdata = full_data, type = "response")
  full_data$Glucose_final <- ifelse(is.na(full_data$Glucose), full_data$Glucose_pred, full_data$Glucose)
  
  # Cox 모델 적합
  cox_model <- coxph(Surv(futime, delta) ~ Glucose_final + Heart.Rate + Height + Mean.blood.pressure +
                       Oxygen.saturation + Respiratory.rate + Temperature + Weight +
                       Glascow.coma.scale.eye.opening + Glascow.coma.scale.verbal.response, 
                     data = full_data)
  return(cox_model$coef)
}


# Case-cohort 방식 분석
case_cohort_analysis <- function(case_cohort_data, full_data) {
  # 전체 데이터에서 검열된 케이스 수와 하위 코호트 내 검열된 케이스 수 계산
  nc <- sum(full_data$delta == 0)       # 전체 검열된 데이터 수
  nc_sub <- sum(case_cohort_data$delta == 0)  # 하위 코호트 내 검열된 데이터 수
  
  # 가중치 설정
  case_cohort_data$weight <- ifelse(
    case_cohort_data$delta == 1,  # 사건 케이스
    1,  # 사건 데이터의 가중치는 1
    ifelse(nc_sub > 0, nc / nc_sub, 1e-6)  # 검열 데이터의 가중치 계산
  )
  
  design <- svydesign(ids = ~1, data = case_cohort_data, weights = ~weight)
  
  # Cox 모델 적합
  cox_model <- svycoxph(Surv(futime, delta) ~ Glucose + Heart.Rate + Height +
                          Oxygen.saturation + Respiratory.rate + Temperature + Weight +
                          Glascow.coma.scale.eye.opening + Glascow.coma.scale.verbal.response, 
                        design)
  return(cox_model$coef)
}



# dfbeta 방식 분석
dfbeta_analysis <- function(case_cohort_data, full_data) {
  # 가중치 초기화
  case_cohort_data$weight <- 1
  
  twophase_design <- svydesign(ids = ~0, weights = ~weight, data = case_cohort_data)
  est_svyglm <- svyglm(Glucose ~ Mean.blood.pressure + Heart.Rate + Height +
                         Oxygen.saturation + Respiratory.rate + Temperature + Weight +
                         Glascow.coma.scale.eye.opening + Glascow.coma.scale.verbal.response, 
                       design = twophase_design)
  
  # Glucose 예측 및 보정
  full_data$Glucose_pred <- predict(est_svyglm, newdata = full_data, type = "response")
  full_data$Glucose_final <- ifelse(is.na(full_data$Glucose), full_data$Glucose_pred, full_data$Glucose)
  
  # Cox 모델 적합 후 dfbeta 계산
  cox_model <- coxph(Surv(futime, delta) ~ Glucose_final + Heart.Rate + Height +
                       Oxygen.saturation + Respiratory.rate + Temperature + Weight +
                       Glascow.coma.scale.eye.opening + Glascow.coma.scale.verbal.response, 
                     data = full_data)
  
  dfbeta_values <- residuals(cox_model, type = "dfbeta")
  full_data$dfbeta_Glucose <- dfbeta_values[, 1] + 1  # 안정성을 위해 1 추가
  
  # 보정된 가중치 생성
  total <- sum(full_data$dfbeta_Glucose, na.rm = TRUE)
  full_data$calibrated_weights <- (full_data$dfbeta_Glucose / total) * nrow(full_data)
  full_data$calibrated_weights[full_data$calibrated_weights <= 0] <- 1e-6  # 음수 또는 0 방지
  
  # 보정된 가중치를 사용한 Cox 모델 적합
  cox_model_calibrated <- coxph(Surv(futime, delta) ~ Glucose_final + Heart.Rate + Height +
                                  Oxygen.saturation + Respiratory.rate + Temperature + Weight +
                                  Glascow.coma.scale.eye.opening + Glascow.coma.scale.verbal.response, 
                                data = full_data, weights = full_data$calibrated_weights)
  
  return(cox_model_calibrated$coef)
}


set.seed(1)
subcohort_size <- 1000 

# Subcohort 데이터 생성
subcohort_result <- case_cohort_sample(data, subcohort_size)
case_cohort_data <- subcohort_result$case_cohort_data
full_data <- subcohort_result$full_data

# Golden Standard
golden_coef

# Case-cohort 방식
case_cohort_coef <- case_cohort_analysis(case_cohort_data, full_data)

# dfbeta 방식
dfbeta_coef <- dfbeta_analysis(case_cohort_data, full_data)

# 결과 비교
cat("Golden Standard Coefficients:\n")
print(golden_coef)

cat("\nCase-cohort Coefficients:\n")
print(case_cohort_coef)

cat("\nDfbeta Coefficients:\n")
print(dfbeta_coef)


#####################################################################

# 시뮬레이션 설정
nsim <- 3000
subcohort_size <- 500 #round(nrow(data) * 0.0333)  # 3.33% 비율 유지
case_cohort_results <- list()
dfbeta_results <- list()


set.seed(1)
for (i in 1:nsim) {
  print(i)
  subcohort_result <- case_cohort_sample(data, subcohort_size)
  case_cohort_data <- subcohort_result$case_cohort_data
  full_data <- subcohort_result$full_data
  
  case_cohort_results[[i]] <- case_cohort_analysis(case_cohort_data, full_data)
  
  dfbeta_results[[i]] <- dfbeta_analysis(case_cohort_data, full_data)
}

case_cohort_matrix <- do.call(rbind, case_cohort_results)
dfbeta_matrix <- do.call(rbind, dfbeta_results)

case_cohort_means <- colMeans(case_cohort_matrix)
dfbeta_means <- colMeans(dfbeta_matrix)

case_cohort_sds <- apply(case_cohort_matrix, 2, sd)
dfbeta_sds <- apply(dfbeta_matrix, 2, sd)

# 결과 출력
cat("Case-cohort 방식 추정치의 평균:\n")
print(case_cohort_means)

cat("\nCase-cohort 방식 추정치의 표준편차:\n")
print(case_cohort_sds)

cat("\nDfbeta 방식 추정치의 평균:\n")
print(dfbeta_means)

cat("\nDfbeta 방식 추정치의 표준편차:\n")
print(dfbeta_sds)


