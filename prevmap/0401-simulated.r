# {Preparation}

library('PrevMap');
library('geoR');

data('data_sim');
# Load data

knots1 <- expand.grid(
  seq(-0.2, 1.2, length = 5),
  seq(-0.2, 1.2, length = 5)
);
knots2 <- expand.grid(
  seq(-0.2, 1.2, length = 10),
  seq(-0.2, 1.2, length = 10)
);
knots3 <- expand.grid(
  seq(-0.2, 1.2, length = 15),
  seq(-0.2, 1.2, length = 15)
);

# {MCML}
par0.exact <- c(0, 1, 0.15);
exact.mcmc <- control.mcmc.MCML(
  n.sim = 65000,
  burnin = 5000,
  thin = 12,
  h = 1.65 / (nrow(data_sim) ^ (1 / 6))
);
system.time(fit.MCML.exact <- binomial.logistic.MCML(
  y ~ 1,
  units.m = ~ units.m,
  coords = ~ x1 + x2,
  data = data_sim,
  par0 = par0.exact,
  start.cov.pars = 0.15,
  control.mcmc = exact.mcmc,
  kappa = 2,
  fixed.rel.nugget = 0,
  method = 'nlminb',
  plot.correlogram = FALSE
));
##001

par0.lr <- c(-0.219294, 0.97945, 0.21393);
lr.mcmc <- control.mcmc.MCML(
  n.sim = 65000,
  burnin = 5000,
  thin = 12,
  h = 1.65 / (nrow(knots1) ^ (1 / 6))
);
system.time(fit.MCML.lr1 <- binomial.logistic.MCML(
  y ~ 1,
  units.m = ~ units.m,
  coords = ~ x1 + x2,
  data = data_sim,
  par0 = par0.lr,
  start.cov.pars = par0.lr[3],
  control.mcmc = lr.mcmc,
  low.rank = TRUE,
  knots = knots1,
  kappa = 2,
  method = 'nlminb',
  plot.correlogram = FALSE
));
##002

lr.mcmc$h <- 1.65 / (nrow(knots2) ^ (1 / 6));
par0.lr <- c(-0.017333, 0.16490, 0.16971);
system.time(fit.MCML.lr2 <- binomial.logistic.MCML(
  y ~ 1,
  units.m = ~ units.m,
  coords = ~ x1 + x2,
  data = data_sim,
  par0 = par0.lr,
  start.cov.pars = par0.lr[3],
  control.mcmc = lr.mcmc,
  low.rank = TRUE,
  knots = knots2,
  kappa = 2,
  method = 'nlminb',
  plot.correlogram = FALSE
));
##003

lr.mcmc$h <- 1.65 / (nrow(knots3) ^ (1 / 6));
par0.lr <- c(-0.031759, 0.30572, 0.18854);
system.time(fit.MCML.lr3 <- binomial.logistic.MCML(
  y ~ 1,
  units.m = ~ units.m,
  coords = ~ x1 + x2,
  data = data_sim,
  par0 = par0.lr,
  start.cov.pars = par0.lr[3],
  control.mcmc = lr.mcmc,
  low.rank = TRUE,
  knots = knots3,
  kappa = 2,
  method = 'nlminb',
  plot.correlogram = FALSE
));
##004

par.hat <- coef(fit.MCML.exact);
Sigma.hat <- varcov.spatial(
  coords = data_sim[c('x1', 'x2')],
  cov.pars = par.hat[2:3],
  kappa = 2
)$varcov;
mu.hat = rep(par.hat[1], nrow(data_sim));
system.time(S.cond.sim <- Laplace.sampling(
  mu = mu.hat,
  Sigma = Sigma.hat,
  y = data_sim$y,
  units.m = data_sim$units.m,
  control.mcmc = exact.mcmc,
  plot.correlogram = FALSE
));
##005

prevalence.sim <- exp(S.cond.sim$samples) / (1 + exp(S.cond.sim$samples));
prevalence.exact <- apply(prevalence.sim, 2, mean);

lr.mcmc$h <- 1.65 / (nrow(knots1) ^ (1 / 6));
system.time(pred.MCML.lr1 <- spatial.pred.binomial.MCML(
  fit.MCML.lr1,
  grid.pred = data_sim[c('x1', 'x2')],
  control.mcmc = lr.mcmc,
  type = 'joint',
  scale.predictions = 'prevalence',
  plot.correlogram = FALSE
));
##006

lr.mcmc$h <- 1.65 / (nrow(knots2) ^ (1 / 6));
system.time(pred.MCML.lr2 <- spatial.pred.binomial.MCML(
  fit.MCML.lr2,
  grid.pred = data_sim[c('x1', 'x2')],
  control.mcmc = lr.mcmc,
  type = 'joint',
  scale.predictions = 'prevalence',
  plot.correlogram = FALSE
));
##007

lr.mcmc$h <- 1.65 / nrow(knots3) ^ (1 / 6);
system.time(pred.MCML.lr3 <- spatial.pred.binomial.MCML(
  fit.MCML.lr3,
  grid.pred = data_sim[c('x1', 'x2')],
  control.mcmc = lr.mcmc,
  type = 'joint',
  scale.predictions = 'prevalence',
  plot.correlogram = FALSE
));
##008

par(mfrow = c(2, 2), mar = c(3, 4, 3, 4));

r.exact <- rasterFromXYZ(
  cbind(data_sim[, c('x1', 'x2')], prevalence.exact)
);
plot(r.exact, zlim = c(0, 1), main = 'Exact method');
contour(r.exact, levels = seq(0.1, 0.9, 0.1), add = TRUE);

plot(
  pred.MCML.lr1,
  'prevalence',
  'predictions',
  zlim = c(0, 1),
  main = 'Low-rank: 25 knots'
);
contour(
  pred.MCML.lr1,
  'prevalence',
  'predictions',
  zlim = c(0, 1),
  levels = seq(0.1, 0.9, 0.1),
  add = TRUE
);

plot(
  pred.MCML.lr2,
  'prevalence',
  'predictions',
  zlim = c(0, 1),
  main = 'Low-rank: 100 knots'
);
contour(
  pred.MCML.lr2,
  'prevalence',
  'predictions',
  zlim = c(0, 1),
  levels = seq(0.1, 0.9, 0.1),
  add = TRUE
);

plot(
  pred.MCML.lr3,
  'prevalence',
  'predictions',
  zlim = c(0, 1),
  main = 'Low-rank: 100 knots'
);
contour(
  pred.MCML.lr3,
  'prevalence',
  'predictions',
  zlim = c(0, 1),
  levels = seq(0.1, 0.9, 0.1),
  add = TRUE
);
#009
