# {Preparation}

# Include library
library('PrevMap');

# Read river blindness data
rb <- read.csv('data/LiberiaRemoData.csv');
summary(rb);
#001

# Read boundaries data
Liberia.bndrs <- read.csv('data/Liberia_bndrs.csv');

# Read predicion grid
names(Liberia.grid.pred);

# Calculate logistic linear model
glm.fit <- glm(
  cbind(npos, ntest - npos) ~ I(utm_x / 1000) + I(utm_y / 1000),
  data = rb,
  family = binomial
);

summary(glm.fit);
#002

# Estiamtes of the regression coefficients
beta.hat <- coef(glm.fit);
#003

# Matrix of the explanatory variables at prediction locations
D.pred <- as.matrix(cbind(1, Liberia.grid.pred / 1000));

# Linear predictor at the prediction locations
eta.hat <- D.pred %*% beta.hat;

# Covariance matrix of the regression coefficients
beta.covar <- vcov(glm.fit);
#004

# Standard errors of the linear predictor
se.eta.hat <- sqrt(diag(D.pred %*% beta.covar %*% t(D.pred)));

# Exceedance probabilities of 20% threshold
exceed.20 <- 1 - pnorm(
  -log(4),
  mean = eta.hat,
  sd = se.eta.hat
);

# Plot of the exceedance probabilities
plot(
  rasterFromXYZ(cbind(Liberia.grid.pred / 1000, exceed.20))
);
lines(Liberia.bndrs / 1000, type = 'l');

# Diagnostic Plots
check.spat <- spat.corr.diagnostic(
  npos ~ I(utm_x / 1000) + I(utm_y / 1000),
  units.m = ~ ntest,
  coords = ~ I(utm_x / 1000) + I(utm_y / 1000),
  data = rb,
  likelihood = 'Binomial',
  uvec = seq(20, 300, length = 15),
  n.sim = 1000
);
#005
