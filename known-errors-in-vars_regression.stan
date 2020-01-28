data {
	int<lower=1> N;
	vector[N] x;
	vector<lower=0>[N] sigma_x;
	vector[N] y;
}

parameters {
	vector[N] z_tilde;

	real mu_z;
	real<lower=0> sigma_z;
	real<lower=0> sigma_y;

	real beta0;
	real beta1;
}

transformed parameters {
	vector[N] z;
	z = mu_z + sigma_z * z_tilde;
}

model {
	mu_z ~ normal(0, 5);
	sigma_z ~ cauchy(0, 5);
	sigma_y ~ cauchy(0, 5);
	beta0 ~ normal(0, 5);
	beta1 ~ normal(0, 5);

	z_tilde ~ normal(0, 1);

	x ~ normal(z, sigma_x);
	y ~ normal(beta0 + beta1 * z, sigma_y);
}

