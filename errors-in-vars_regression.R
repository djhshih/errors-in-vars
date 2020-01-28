library(rstan)
library(ggplot2)
library(ggsci)
library(io)

# Assess how various regression models cope with measurement error in the 
# predictor variable


# Simulate data with measurement error

set.seed(1)

out.fname <- filename("m-error");
pdf.fname <- insert(out.fname, ext="pdf");

N <- 1000;

params <- list(
	mu_z = 0.5,
	sigma_z = sqrt(0.3),
	sigma_x = sqrt(0.7),
	sigma_y = sqrt(1),
	beta0 = 0,
	beta1 = 2
);

# latent predictor variable
z <- with(params, rnorm(N, mu_z, sigma_z));

# observed predictor variable containing measurement error
x <- with(params, rnorm(N, z, sigma_x));

# observed response variable is based on the latent predictor variable
y <- with(params, rnorm(N, beta0 + beta1*z, sigma_y));

# verify that Var[X] = sigma_z^2 + sigma_x^2
sd(x)

# verify that rho^2 matches sigma_z^2 / (sigma_z^2 + sigma_x^2)
cor(z, x)^2


# Fit a least-squares linear regression

fit.lm <- lm(y ~ x);
summary(fit.lm)

# the slope is biased towards 0

# A Bayesian linear regression that does not account for
# measurement error in x gives the same results (not shown)

####

# Fit a errors-in-variables Bayesian regression model

d <- list(
	N = N,
	x = x,
	y = y
);

fit.eiv <- stan("errors-in-vars_regression.stan", data=d,
	pars=c("mu_z", "sigma_z", "sigma_x", "sigma_y", "beta0", "beta1"),
	iter=20000, control=list(adapt_delta=0.9), thin=2, seed=1
);

check_divergences(fit.eiv)
check_hmc_diagnostics(fit.eiv)

print(fit.eiv, par=c("mu_z", "sigma_z", "sigma_x", "sigma_y", "beta0", "beta1"))
stan_par(fit.eiv, par="beta0")
stan_par(fit.eiv, par="beta1")
pairs(fit.eiv, pars=c("beta0", "beta1"))

qwrite(fit.eiv, "errors-in-vars_regression.stanfit.rds");

####

estimate_y <- function(fit, z.new, alpha=0.05) {
	q <- c(alpha/2, 1 - alpha/2);

	beta0.hat <- extract(fit, "beta0")[[1]];
	beta1.hat <- extract(fit, "beta1")[[1]];

	beta1.hat.mean <- mean(beta1.hat);
	beta1.hat.q <- quantile(beta1.hat, q);

	beta0.hat.mean <- mean(beta0.hat);
	beta0.hat.q <- quantile(beta0.hat, q);

	y.hat <- beta0.hat.mean + z.new * beta1.hat.mean;
	y.tilde <- mapply(function(b0, b1) b0 + z.new * b1, beta0.hat, beta1.hat);

	# y.tild.mean == y.hat
	#y.tilde.mean <- apply(y.tilde, 1, mean);
	y.q <- apply(y.tilde, 1, function(y) quantile(y, q));

	est <- data.frame(
		z = z.new,
		y = y.hat,
		ymin = y.q[1,],
		ymax = y.q[2,]
	);

	list(
		params = c(
			beta0 = beta0.hat.mean, beta0.ci = beta0.hat.q,
			beta1 = beta1.hat.mean, beta1.ci = beta1.hat.q
		),
		estimate = est
	)
}

z.new <- seq(min(x), max(x), 0.1);

pred.eiv <- estimate_y(fit.eiv, z.new)$estimate;

####

# Fit a known-errors-in-variables Bayesian regression model
# where the errors in x are known and given
# This is similar to meta-analysis models

d2 <- d;
d2$sigma_x <- rep(params$sigma_x, N);

fit.keiv <- stan("known-errors-in-vars_regression.stan", data=d2,
	pars=c("mu_z", "sigma_z", "sigma_y", "beta0", "beta1"),
	iter=20000, control=list(adapt_delta=0.9), thin=2, seed=1
);

check_divergences(fit.keiv)
check_hmc_diagnostics(fit.keiv)

print(fit.keiv, par=c("mu_z", "sigma_z", "sigma_y", "beta0", "beta1"))
stan_par(fit.keiv, par="beta0")
stan_par(fit.keiv, par="beta1")
pairs(fit.keiv, pars=c("beta0", "beta1"))

qwrite(fit.keiv, "known-errors-in-vars_regression.stanfit.rds");

pred.keiv <- estimate_y(fit.keiv, z.new)$estimate;

####

# Plot and compare the results

data <- data.frame(
	x = x,
	z = z,
	y = y
);

my.coord <- coord_cartesian(xlim=c(-2, 2.5), ylim=c(-4, 6));
no.grid <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank());

options(plot=list(width = 3, height = 3));

qdraw(
	ggplot(data, aes(x=z, y=y)) +
		geom_point(alpha=0.5) + theme_classic() +
		geom_abline(intercept=params$beta0, slope=params$beta1, colour="grey30", linetype=2) +
		my.coord
	,
	file = insert(pdf.fname, "truth")
)

qdraw(
	ggplot(data, aes(x=x, y=y)) +
		geom_point(aes(x=z, y=y), alpha=0.5, colour="forestgreen") + theme_bw() +
		geom_point(alpha=0.5) + theme_classic() +
		geom_abline(intercept=params$beta0, slope=params$beta1, colour="grey30", linetype=2) +
		my.coord
	,
	file = insert(pdf.fname, "add-noise")
)

qdraw(
	ggplot(data, aes(x=x, y=y)) +
		geom_point(alpha=0.5) + theme_classic() +
		my.coord
	,
	file = insert(pdf.fname, "observed")
);


g <- ggplot(data, aes(x=x, y=y)) +
	geom_point(alpha=0.5) + theme_classic() +
	scale_colour_manual("", breaks=names(group.colours), values=group.colours) +
	scale_fill_manual("", breaks=names(group.colours), values=group.colours) +
	geom_abline(intercept=params$beta0, slope=params$beta1, colour="grey30", linetype=2) +
	my.coord +
	theme(
		legend.justification = c("right", "bottom"),
		legend.position=c(0.99, 0.01),
		legend.background=element_blank()
	)
qdraw(g, file = insert(pdf.fname, "base"));

group.colours <- pal_nejm()(3);
names(group.colours) <- c("regression", "errors-in-variables", "known-eiv");

geom.lm <- geom_smooth(method = "lm", aes(colour="regression", fill="regression"));
geom.eiv <- geom_smooth(data = pred.eiv, aes(x = z, y = y, ymin = ymin, ymax = ymax, colour="errors-in-variables", fill="errors-in-variables"), stat="identity");
geom.keiv <- geom_smooth(data = pred.keiv, aes(x = z, y = y, ymin = ymin, ymax = ymax, colour="known-eiv", fill="known-eiv"), stat="identity");

qdraw(g + geom.lm, file = insert(pdf.fname, "lm"));

qdraw(g + geom.eiv, file = insert(pdf.fname, "eiv"));

qdraw(g + geom.keiv, file = insert(pdf.fname, "keiv"));

qdraw(g + geom.lm + geom.eiv + geom.keiv, file = insert(pdf.fname, "all"));

qdraw(
	g + geom.lm + geom.eiv + geom.keiv +
		theme(legend.position = "none")
	,
	file = insert(pdf.fname, c("all", "no-legend"))
);

