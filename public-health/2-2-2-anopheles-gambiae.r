# {Preparation}
library('PrevMap');

# Read data
mosq <- read.csv('data/anopheles.csv');
summary(mosq);
##001

elev <- read.csv('data/anopheles_elevation.csv')

# {Sampled Locations}
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1));
elev <- rasterFromXYZ(elev);
plot(elev);
points(mosq[, c('web_x', 'web_y')], pch=20);
# ↑ Plots of sampled locations
##002

# {Chart of Relationship between Positive Cases and Elevation}
plot(
  log(An.gambiae) ~ elevation,
  data = mosq,
  pch = 20,
  cex = 0.5
);
glm.fit <- glm(
  An.gambiae ~ elevation,
  data = mosq,
  family = poisson
); # ← Create generalised linear model of Anopheles gambiae to elevation data
summary(glm.fit);
##003

abline(glm.fit);
##004

# Diagnostic Plots
check.spat <- spat.corr.diagnostic(
  An.gambiae ~ elevation,
  coords = ~ I(web_x / 1000) + I(web_y / 1000),
  data = mosq,
  likelihood = 'Poisson',
  uvec = seq(0, 10, length = 15),
  n.sim = 10000
);
##005
