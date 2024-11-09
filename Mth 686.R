
## MTH686 Project
# Step 1: Load the Data
data <- read.table("data-22.txt", header = FALSE)
colnames(data) <- c("t", "y")

# Question 1

# Model 1
model1 <- try(nls(y ~ alpha0 + alpha1 * exp(beta1 * t) + alpha2 * exp(beta2 * t),
                  data = data, 
                  start = list(alpha0 = 8, alpha1 = 1, beta1 = -0.1, alpha2 = 1, beta2 = -0.05),
                  control = nls.control(maxiter = 100, tol = 1e-5)), silent = TRUE)

# Model 2
model2 <- try(nls(y ~ (alpha0 + alpha1 * t) / (beta0 + beta1 * t),
                  data = data, 
                  start = list(alpha0 = 8, alpha1 = 0.1, beta0 = 1, beta1 = 0.1),
                  control = nls.control(maxiter = 100, tol = 1e-5)), silent = TRUE)

# Model 3
model3 <- try(nls(y ~ beta0 + beta1 * t + beta2 * t^2 + beta3 * t^3 + beta4 * t^4,
                  data = data, 
                  start = list(beta0 = 8, beta1 = 0.5, beta2 = 0.1, beta3 = 0.05, beta4 = 0.01),
                  control = nls.control(maxiter = 100, tol = 1e-5)), silent = TRUE)


# Check which models fit successfully
if(class(model1) == "try-error") {
  cat("Model 1 failed to converge.\n")
} else {
  print(summary(model1))
}

if(class(model2) == "try-error") {
  cat("Model 2 failed to converge.\n")
} else {
  print(summary(model2))
}

if(class(model3) == "try-error") {
  cat("Model 3 failed to converge.\n")
} else {
  print(summary(model3))
}

# Compare RSS if models converge
if(class(model1) != "try-error") {
  rss1 <- sum(residuals(model1)^2)
  cat("RSS for Model 1:", rss1, "\n")
}

if(class(model2) != "try-error") {
  rss2 <- sum(residuals(model2)^2)
  cat("RSS for Model 2:", rss2, "\n")
}

if(class(model3) != "try-error") {
  rss3 <- sum(residuals(model3)^2)
  cat("RSS for Model 3:", rss3, "\n")
}


# Question 3


# Calculate Residual Sum of Squares (RSS) for Model 3
rss3 <- sum(residuals(model3)^2)

# Calculate AIC for Model 3
aic3 <- AIC(model3)

# Calculate BIC for Model 3
n <- nrow(data)  
k <- length(coef(model3))  
bic3 <- k * log(n) + n * log(rss3 / n)

# Display the results
cat("Model 3 (Polynomial Model):\n")
cat("RSS:", rss3, "\n")
cat("AIC:", aic3, "\n")
cat("BIC:", bic3, "\n")



# Question 4

# Calculate RSS from Model 3
rss3 <- sum(residuals(model3)^2)

n <- nrow(data)
k <- length(coef(model3))

# Estimate of sigma^2
sigma_squared <- rss3 / (n - k)

# Display the result
cat("Estimated sigma^2:", sigma_squared, "\n")

#Question 5

rss3 <- sum(residuals(model3)^2)
n <- nrow(data)
k <- length(coef(model3))
sigma_squared <- rss3 / (n - k)

# Calculate the Fisher Information Matrix
library(numDeriv)
jacobian_matrix <- jacobian(function(params) {
  with(as.list(params), {
    beta0 + beta1 * data$t + beta2 * data$t^2 + beta3 * data$t^3 + beta4 * data$t^4
  })
}, coef(model3))

fisher_information <- t(jacobian_matrix) %*% jacobian_matrix / sigma_squared

# Calculate the variance-covariance matrix
cov_matrix <- sigma_squared * solve(fisher_information)

# Calculate 95% Confidence Intervals
params <- coef(model3)
se_params <- sqrt(diag(cov_matrix))
z_value <- qnorm(0.975)  # for 95% confidence

conf_intervals <- data.frame(
  Estimate = params,
  Lower = params - z_value * se_params,
  Upper = params + z_value * se_params
)

print(conf_intervals)


# Question 6


# Residuals from Model 3
residuals_model3 <- residuals(model3)

# Plot 1: Residuals vs. Fitted Values
fitted_values <- fitted(model3)
plot(fitted_values, residuals_model3, main = "Residuals vs Fitted Values",
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# Plot 2: Normal Q-Q Plot
qqnorm(residuals_model3, main = "Normal Q-Q Plot of Residuals")
qqline(residuals_model3, col = "red")

# Plot 3: Residuals vs. Time (t)
plot(data$t, residuals_model3, main = "Residuals vs Time",
     xlab = "Time (t)", ylab = "Residuals")
abline(h = 0, col = "red")




# Question 7


# Residuals from Model 3
residuals_model3 <- residuals(model3)

# Shapiro-Wilk Normality Test
shapiro_test <- shapiro.test(residuals_model3)
cat("Shapiro-Wilk Test p-value:", shapiro_test$p.value, "\n")

# Q-Q Plot for visual inspection
qqnorm(residuals_model3, main = "Normal Q-Q Plot of Residuals")
qqline(residuals_model3, col = "red")



# Question 8


# Observed data points
plot(data$t, data$y, main = "Observed Data Points and Fitted Curve",
     xlab = "Time (t)", ylab = "Observed and Fitted Values",
     pch = 16, col = "blue", cex = 0.6)

# Fitted curve from Model 3
lines(data$t, fitted(model3), col = "red", lwd = 2)

legend("topright", legend = c("Observed Data", "Fitted Curve"), 
       col = c("blue", "red"), pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2))

