data {
	int<lower=1> N;
	vector[N] x;
	vector[N] y;
}

parameters {
	vector[N] z;

	real mu_z;
	real<lower=0> sigma_z;

	real<lower=0> sigma_x;
	real<lower=0> sigma_y;

	real beta0;
	real beta1;
}

model {
	z ~ normal(mu_z, sigma_z);
	x ~ normal(z, sigma_x);
	y ~ normal(beta0 + beta1 * z, sigma_y);
}

