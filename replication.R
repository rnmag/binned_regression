library(margins)
library(ggplot2)

# set.seed(seed = NULL)
set.seed(123)
n <- 200

# Variables are defined according to pages 8 and 9
X <- rnorm(n, 3, 1) # X for the linear simulations
X_un <- runif(n, -3, 3) # X for the non-linear simulation

epsilon <- rnorm(n, 0, 4) # Error term
x_terciles <- as.factor(findInterval(X, quantile(X, c(1/3, 2/3))) + 1) # Terciles for X
levels(x_terciles) <- c("X: Low", "X: Med", "X: High")

x_un_terciles <- as.factor(findInterval(X_un, quantile(X_un, c(1/3, 2/3))) + 1) # Terciles for X_un
levels(x_un_terciles) <- c("X_un: Low", "X_un: Med", "X_un: High")

D1 <- rbinom(n, 1, .5) # Binary treatment
D2 <- rnorm(n, 3, 1) # Continuous treatment

Y1 <- 5 - 4*X - 9*D1 + 3*D1*X + epsilon 
Y2 <- 5 - 4*X - 9*D2 + 3*D2*X + epsilon 
Y3 <- 2.5 - X_un^2 - 5*D1 + 2*D1*X_un^2 + epsilon # Non-linear with binary treatment

# Creating database
db <- data.frame(Y1, Y2, Y3, D1, D2, X, X_un, epsilon, x_terciles, x_un_terciles)
View(db)

# Figure 1a - linear, binary treatment
ggplot(data = db) +
  aes(x = X, y = Y1) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, color = "blue", size = 1.2) +
  geom_smooth(se = FALSE, color = "red", size = 1.2) + 
  facet_wrap(~D1)

# Figure 1b - non-linear, binary treatment
ggplot(data = db) +
  aes(x = X_un, y = Y3) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, color = "blue", size = 1.2) +
  geom_smooth(se = FALSE, color = "red", size = 1.2) + 
  facet_wrap(~D1)

# Figure 1c - continuous treatment
ggplot(data = db) +
  aes(x = D2, y = Y2) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, color = "blue", size = 1.2) +
  geom_smooth(se = FALSE, color = "red", size = 1.2) + 
  facet_wrap(~x_terciles)
  
# Binning estimator
# The following lines attempt to translate the wording on section 4.1 directly
# into code. The relevant excerpts are quoted for comparison.

# "First, we discretize the moderator variable X into three bins (respectively
# corresponding to the three terciles) as before and create a dummy variable for
# each bin."
dummies <- model.matrix(~ factor(x_terciles) - 1, db)
db <- cbind(db, dummies)
colnames(db)[11:13] <- c("x_low", "x_med", "x_high")

dummies_un <- model.matrix(~ factor(x_un_terciles) - 1, db)
db <- cbind(db, dummies_un)
colnames(db)[14:16] <- c("x_un_low", "x_un_med", "x_un_high")

# "Second, we pick an evaluation point within each bin, x1, x2, and x3, where we
# want to estimate the conditional marginal effect of D on Y . Typically, we choose 
# x1, x2, and x3 to be the median of X in each bin, but researchers are free to
# choose other numbers within the bins (for example, the means)."
x1 <- median(db$X[db$x_low == 1])
x2 <- median(db$X[db$x_med == 1])
x3 <- median(db$X[db$x_high == 1])

x_un1 <- median(db$X_un[db$x_un_low == 1])
x_un2 <- median(db$X_un[db$x_un_med == 1])
x_un3 <- median(db$X_un[db$x_un_high == 1])

# "Third, we estimate a model that includes interactions between the bin dummies
# G and the treatment indicator D, the bin dummies and the moderator X minus the
# evaluation points we pick (x1, x2, and x3), as well as the triple interactions."

# Discounting the median
db$centered_X <- ifelse(db$x_low == 1, X - x1,
                          ifelse(db$x_med == 1, X - x2, X - x3))

db$centered_X_un <- ifelse(db$x_un_low == 1, X_un - x_un1,
                            ifelse(db$x_un_med == 1, X_un - x_un2, X_un - x_un3))

# Binned model - linear
linear_me <- by(db, db$x_terciles, function(db) lm(Y1 ~ D1 * centered_X, db))

m_effects <- cbind(me = sapply(linear_me, function(y) coef(y)["D1"]),
                   var = sapply(linear_me, function(y) diag(vcov(y))["D1"]))

# Binned model - non-linear
non_linear_me <- by(db, db$x_un_terciles, function(db) lm(Y3 ~ D1 * I(centered_X_un^2), db))

nl_m_effects <- cbind(me = sapply(non_linear_me, function(y) coef(y)["D1"]),
                      var = sapply(non_linear_me, function(y) diag(vcov(y))["D1"]))

# Figure 2a
two_a <- cplot(lm(Y1 ~ D1 * X, db), x = "X", dx = "D1", what = "effect")
arrows(x0 = c(x1, x2, x3),  
       y0 = c(m_effects[1, 1] - m_effects[1, 2], m_effects[2, 1] - m_effects[2, 2], m_effects[3, 1] - m_effects[3, 2]),
       y1 = c(m_effects[1, 1] + m_effects[1, 2], m_effects[2, 1] + m_effects[2, 2], m_effects[3, 1] + m_effects[3, 2]),
       length = 0.05, angle = 90, code = 3, lwd = 2)
points(c(x1, x2, x3),  c(m_effects[1, 1], m_effects[2, 1], m_effects[3, 1]),
       pch = 19)
text(c(x1, x2, x3), 7, c("L", "M", "H"))


# Figure 2b
two_b <- cplot(lm(Y3 ~ D1 * I(X_un^2), db), x = "X_un", dx = "D1", what = "effect")
arrows(x0 = c(x_un1, x_un2, x_un3),  
       y0 = c(nl_m_effects[1, 1] - nl_m_effects[1, 2], nl_m_effects[2, 1] - nl_m_effects[2, 2], nl_m_effects[3, 1] - nl_m_effects[3, 2]),
       y1 = c(nl_m_effects[1, 1] + nl_m_effects[1, 2], nl_m_effects[2, 1] + nl_m_effects[2, 2], nl_m_effects[3, 1] + nl_m_effects[3, 2]),
       length = 0.05, angle = 90, code = 3, lwd = 2)
points(c(x_un1, x_un2, x_un3),  c(nl_m_effects[1, 1], nl_m_effects[2, 1], nl_m_effects[3, 1]),
       pch = 19)
abline(lm(yvals ~ xvals, data = two_b), col = "red", lwd = 2)
text(c(x_un1, x_un2, x_un3), 13, c("L", "M", "H"))
