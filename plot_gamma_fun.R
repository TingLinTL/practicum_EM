
# Parameters
shape_val <- 0.5 # very skewed
scale_val <- 0.5   # small scale

# Grid for x values
x_vals <- seq(0, 10, length.out = 500)   # we only need small range (most mass is near 0)

# Compute PDF
y_vals <- dgamma(x_vals, shape = shape_val, scale = scale_val)

# Plot
plot(x_vals, y_vals, type = "l", lwd = 2, col = "blue",
     main = "Gamma(0.5, 0.5) PDF",
     xlab = "u", ylab = "Density")

# Add mean line
gamma_mean <- shape_val * scale_val  # mean = 0.25
abline(v = gamma_mean, col = "red", lty = 2)
legend("topright",
       legend = c("Gamma(0.5,0.5)", "Mean (0.25)"),
       col = c("blue", "red"),
       lty = c(1,2), lwd = 2)
