# {Preparation}
library('PrevMap');
library('splines');

# Load data
Ghana.bndrs <- read.csv('data/Ghana_bndrs.csv');
names(Ghana.bndrs);
##001

maln <- read.csv('data/malnutrition.csv');
names(maln);
summary(maln);
##002

# Create coordinates data of sampled locations
ID.coords <- create.ID.coords(maln, ~utm_x + utm_y);
HAZ.avg <- tapply(maln$HAZ, ID.coords, mean);
coords <- unique(maln[, c('utm_x', 'utm_y')]);

# {Plots of data}
par(mfrow = c(2, 2), mar = c(4, 4, 2, 2));

# Plot locations of sampled clusters
plot(
  Ghana.bndrs,
  type = 'l',
  asp = 1,
  xlab = '',
  ylab = '',
  main = c('(a)')
); # ← Plot the border of area
points(coords, pch = 20, cex = 0.5);
# ↑ Plot the location of sampled locations

# Scatter plot of height-for-age Z-scores against age
plot(
  HAZ ~ age,
  data = maln,
  main = '(b)',
  xlab = 'Age',
  pch = 20,
  cex = 0.5
);

# Box-plots of height-for-age Z-scores to each level of education
plot(
  HAZ ~ factor(edu),
  data = maln,
  main = c('(c)'),
  xlab='Maternal education'
);

# Box-plots of height-for-age Z-scores to each score of household wealth
plot(
  HAZ ~ factor(wealth),
  data = maln,
  main = '(d)',
  xlab='Wealth index'
);
##003

# {Broken Stick}

par(mfrow = c(1,1));

max.vec <- function(val, x) {
  return(sapply(x, function(i) {
    return(max(0, i - val));
  }));
}; # ← a function that returns the vector whose element less than `val` updated to 0
lm.fit <- lm(
  HAZ ~
    age
    + I(max.vec(1, age))
    + I(max.vec(2, age))
    + edu
    + wealth,
  data = maln
); # ← create a linear model to fit HAZ to age in different range, education, and wealth
summary(lm.fit);
##004

# Plot broken stick
beta.hat <- coef(lm.fit);
##005
age.set <- seq(0, 5, length = 1000);
age.vars <- cbind(
  age.set,
  max.vec(1, age.set),
  max.vec(2, age.set)
); # ← Match age.set with its value from max.vec function
broken.sticks <- as.numeric(age.vars %*% beta.hat[2:4]);
std.errors <- sqrt(
  diag(age.vars %*% vcov(lm.fit)[2:4,2:4] %*% t(age.vars))
);
ci.95 <- cbind(
  broken.sticks - qnorm(0.975) * std.errors,
  broken.sticks + qnorm(0.975) * std.errors
);

matplot(
  age.set,
  cbind(ci.95, broken.sticks),
  type = 'l',
  xlab = 'Age (years)',
  ylab = '',
  lty = c('dashed', 'dashed', 'solid'),
  col = 1
);
##006
