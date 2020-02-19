# {Preparing}

# Install packages
install.packages('PrevMap');
install.packages('geoR');

# Using packages
library('PrevMap');
library('geoR');
library('splancs');

# Load data
data('loaloa');

# {Explore data}

# Data columns and summary
names(loaloa);
summary(loaloa);
##001

# {Finding `K` value}

loaloa$logit <- log(
  (loaloa$NO_INF + 0.5) / (loaloa$NO_EXAM - loaloa$NO_INF + 0.5)
);
# Create logit -(1)

plot(loaloa$NO_INF, loaloa$logit);
# Try to plot
##002

attach(loaloa);
profile.kappa <- shape.matern(
  formula = logit ~ 1,
  coords = ~ LONGITUDE + LATITUDE,
  data = loaloa,
  set.kappa = seq(0.2, 1.5, length = 15),
  start.par = c(0.2, 0.05),
  coverage = 0.95
);
# help(shape.matern) = plot profile likelihood of Matern covariance function used in linear Gaussian model
##003

c(profile.kappa$lower, profile.kappa$upper);
# Confidence Interval
profile.kappa$kappa.hat;
# Value of `K` (near 1/2)
##004

# ----

# {Least-Squares Estimation of the Empirical Variogramme}
coords <- as.matrix(loaloa[, c('LONGITUDE', 'LATITUDE')]);
head(coords)
##005

vari <- variog(
  coords = coords,
  data = loaloa$logit,
  uvec = c(0, 0.1, 0.15, 0.2, 0.4, 0.8, 1.4, 1.8, 2, 2.5, 3)
);
vari.fit <- variofit(
  vari,
  ini.cov.pars = c(2, 0.2),
  cov.model = 'matern',
  fix.nugget = FALSE,
  nugget = 0,
  fix.kappa = TRUE,
  kappa = 0.5
);

par(mfrow = c(1, 2));
plot(
  coords,
  pch = 20,
  asp = 1,
  cex = 0.5,
  main = '(a)'
);
plot(vari, main = '(b)');
lines(vari.fit);
##006

vari.fit
##007

# ---

# {Linear Model with Maximum Likelihood Method}

fit.MLE <- linear.model.MLE(
  formula = logit ~ 1,
    # ↑ fitting an intercept
  coords = ~ LONGITUDE + LATITUDE,
  data = loaloa,
  start.cov.pars = c(0.2, 0.15),
    # ↑ initial value of Φ and v²=τ²/σ²
  kappa = 0.5
);
# ※Finding maximum likelihood?
summary(fit.MLE, log.cov.pars = FALSE);
##008

cp1 <- control.profile(
  rel.nugget = exp(seq(-5, 0, length = 20))
);
cp2 <- control.profile(
  rel.nugget = exp(seq(-4, 0, length = 20)),
  phi = exp(seq(-4, 4, length = 20))
);
# Set parameters profiles

lp1 <- loglik.linear.model(
  fit.MLE,
  cp1,
  plot.profile = FALSE
);
lp2 <- loglik.linear.model(
  fit.MLE,
  cp2,
  plot.profile = FALSE
);

par(mfrow = c(1, 2));
plot(
  lp1,
  type = 'l',
  log.scale = TRUE,
  xlab = expression(log(nu ^ 2)),
  ylab = 'log-likelihood',
  main = expression('Profile likelihood for ' ~ nu ^ 2)
);
plot(
  lp2,
  log.scale = TRUE,
  xlab = expression(log(phi)),
  ylab = expression(log(nu ^ 2)),
  main = expression('Profile likelihood for ' ~ nu ^ 2 ~ " and " ~ phi)
);
##009

ci0.95 <-loglik.ci(
  lp1,
  coverage = 0.95,
  plot.spline.profile = FALSE
);
##010

# ----

# {Binomial Logistic Model}

# {{Likelihood-Based}}
fit.glm <- glm(
  cbind(NO_INF, NO_EXAM - NO_INF) ~ 1,
  data = loaloa,
  family = binomial
);
par0 <- c(
  coef(fit.glm),
  vari.fit$cov.pars,
  vari.fit$nugget
);
c.mcmc <- control.mcmc.MCML(
  n.sim = 10000,
    # ↑ Iterations
  burnin = 2000,
    # ↑↓ Retaining every 8th sample after a burn-in of 2,000 values
  thin = 8,
  h = 1.65 / (nrow(loaloa) ^ (1 / 6))
    # ↑ Proposal density
);
# Set control paramters
fit.MCML1 <- binomial.logistic.MCML(
  formula = NO_INF ~ 1,
  units.m = ~ NO_EXAM,
  par0 = par0,
  coords = ~ LONGITUDE + LATITUDE,
  data = loaloa,
  control.mcmc = c.mcmc,
  kappa = 0.5,
  start.cov.pars = c(par0[3], par0[4] / par0[2])
);
##011

fit.MCML1$log.lik;
##012

# {{Repeat MCML with another values of β₀ and θ₀}}
par0 <- coef(fit.MCML1);
start <- c(par0[3], par0[4] / par0[2]);
fit.MCML2 <- binomial.logistic.MCML(
  formula = NO_INF ~ 1,
  units.m = ~ NO_EXAM,
  par0 = par0,
  coords = ~ LONGITUDE + LATITUDE,
  data = loaloa,
  control.mcmc = c.mcmc,
  kappa = 0.5,
  start.cov.pars = c(par0[3], par0[4] / par0[2])
);
fit.MCML2$log.lik;
##013

c.mcmc <- control.mcmc.MCML(
  n.sim = 65000,
  burnin = 5000,
  thin = 6,
  h = 1.65 / (nrow(loaloa) ^ (1 / 6))
);
par0 <- coef(fit.MCML2);
fit.MCML3 <- binomial.logistic.MCML(
  formula = NO_INF ~ 1,
  units.m = ~ NO_EXAM,
  par0 = par0,
  coords = ~ LONGITUDE + LATITUDE,
  data = loaloa,
  control.mcmc = c.mcmc,
  kappa = 0.5,
  start.cov.pars = c(par0[3], par0[4] / par0[2])
);
summary(fit.MCML3);
##014

# ---

## {{Summarise Predictive Distribution}}

poly <- coords[chull(coords),];
grid.pred <- gridpts(poly, xs = 0.1, ys = 0.1);
pred.MCML <- spatial.pred.binomial.MCML(
  fit.MCML3,
  grid.pred,
  control.mcmc = c.mcmc,
  type = 'marginal',
    # ↑ or 'joint'
  scale.prediction = 'prevalence',
    # ↑ or 'logit' or 'odds'
  standard.errors = TRUE,
  thresholds = 0.2,
  scale.thresholds = 'prevalence'
);

par(mfrow = c(1, 3));

plot(
  pred.MCML,
  type = 'prevalence',
  summary = 'predictions',
  zlim = c(0, 0.45),
  main = 'Prevalence - predictions\n(classical analysis)'
);
contour(
  pred.MCML,
  type = 'prevalence',
  summary = 'predictions',
  levels = c(0.05, 0.1, 0.2, 0.3),
  add = TRUE
);

plot(
  pred.MCML,
  type = 'prevalence',
  summary = 'predictions',
  zlim = c(0, 0.3),
  main = 'Prevalence - predictions\n(classical analysis)'
);
contour(
  pred.MCML,
  type = 'prevalence',
  summary = 'standard.errors',
  levels = c(0.05, 0.1, 0.15, 0.2),
  add = TRUE
);

plot(
  pred.MCML,
  summary = 'exceedance.prob',
  zlim = c(0, 1),
  main = 'Prevalence - exceedance probabilities\n(classical analysis)'
);
contour(
  pred.MCML,
  type = 'exceedance.prob',
  levels = c(0.1, 0.4, 0.5, 0.7),
  add = TRUE
);
##015

# {{Diagnostic Plots}}
par(mfrow = c(3, 3));
S.mean <- apply(pred.MCML$samples, 2, mean);
acf(S.mean, main = '');
plot(S.mean, type = 'l');
plot(ecdf(S.mean[1:5000]), main = '');
lines(ecdf(S.mean[5001:10000]), col = 2, lty = 'dashed');

ind.S <- sample(1:nrow(grid.pred), 2);
acf(pred.MCML$samples[ind.S[1],], main = '');
plot(
  pred.MCML$samples[ind.S[1],],
  ylab = paste('Component n.', ind.S[1]),
  type = 'l'
);
plot(ecdf(pred.MCML$samples[ind.S[1], 1:5000]), main = '');
lines(ecdf(pred.MCML$samples[ind.S[1], 5001:10000]), col = 2, lty = 'dashed');

acf(pred.MCML$samples[ind.S[2],], main = '');
plot(
  pred.MCML$samples[ind.S[2],],
  ylab = paste('Component n.', ind.S[2]),
  type = 'l'
);
plot(ecdf(pred.MCML$samples[ind.S[2], 1:5000]), main = '');
lines(ecdf(pred.MCML$samples[ind.S[2], 5001:10000]), col = 2, lty = 'dashed');
##016

# ----
dev.off();

# {Bayesian Analysis}
cp <- control.prior(
  beta.mean = 0,
  beta.covar = 100^2,
  log.normal.sigma2 = c(1, 5),
  uniform.phi = c(0, 8),
  log.normal.nugget = c(-3, 1)
);
# configure prioir distribution

mcmc.Bayes <- control.mcmc.Bayes(
  n.sim = 6000,
  burnin = 1000,
  thin = 1,
  h.theta1 = 1,
  h.theta2 =  0.7,
  h.theta3 = 0.05,
  L.S.lim = c(5, 50),
  epsilon.S.lim = c(0.03, 0.06),
  start.beta = -2.3,
  start.sigma2 = 2.6,
  start.phi = 0.8,
  start.nugget = 0.05,
  start.S = predict(fit.glm)
);

fit.Bayes <- binomial.logistic.Bayes(
  formula = NO_INF ~ 1,
  units.m = ~ NO_EXAM,
  coords = ~ LONGITUDE + LATITUDE,
  data = loaloa,
  control.prior = cp,
  control.mcmc = mcmc.Bayes,
  kappa = 0.5
);
# Fitting Bayesian Binomial Logistics model
summary(fit.Bayes, hpd.coverage = 0.95);
##017

# {{}}

# {{Autocorrelation plots}}
par(mfrow = c(2, 4));
autocor.plot(fit.Bayes, param = 'beta', component.beta = 1);
autocor.plot(fit.Bayes, param = 'sigma2');
autocor.plot(fit.Bayes, param = 'phi');
autocor.plot(fit.Bayes, param = 'tau2');
i <- sample(1:nrow(loaloa), 4);
autocor.plot(fit.Bayes, param = 'S', component.S = i[1]);
autocor.plot(fit.Bayes, param = 'S', component.S = i[2]);
autocor.plot(fit.Bayes, param = 'S', component.S = i[3]);
autocor.plot(fit.Bayes, param = 'S', component.S = i[4]);
##019

# {{Predictive Distribution}}
pred.Bayes <- spatial.pred.binomial.Bayes(
  fit.Bayes,
  grid.pred,
  type = 'marginal',
  scale.prediction = 'prevalence',
  quantiles = NULL,
  standard.errors = TRUE,
  thresholds = 0.2,
  scale.thresholds = 'prevalence'
);

par(mfrow = c(1, 3));
plot(
  pred.Bayes,
  type = 'prevalence',
  summary = 'predictions',
  zlim = c(0, 0.45),
  main = 'Prevalence - predictions\n(Bayesian analysis)'
);
contour(
  pred.Bayes,
  type = 'prevalence',
  summary = 'predictions',
  levels = c(0.05, 0.1, 0.2, 0.3),
  add = TRUE
);

plot(
  pred.Bayes,
  type = 'prevalence',
  summary = 'standard.errors',
  zlim = c(0, 0.3),
  main = 'Prevalence - standard errors\n(Bayesian analysis)'
);
contour(
  pred.Bayes,
  type = 'prevalence',
  summary = 'standard.errors',
  levels = c(0.05, 0.1, 0.15, 0.2),
  add = TRUE
);

plot(
  pred.Bayes,
  type = 'exceedance.prob',
  zlim = c(0, 1),
  main = 'Prevalence - exceedance probabilities\n(Bayesian analysis)'
);
contour(
  pred.Bayes,
  type = 'exceedance.prob',
  levels = c(0.1, 0.4, 0.5, 0.7),
  add = TRUE
);
##020
