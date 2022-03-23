%% Initial Setup
clear variables
% clear all

%% Personal parameters
mu_A    = 0.5;
mu_B    = 0.5;
mu_C    = 0.5;
sigma_a = 0.1;
sigma_w = 0.1;
ro_1    = 0.1;
ro_2    = 0.1;
t_0_ms  = 5;
t_0_s   = t_0_ms * 1e-3;
f_0     = 40;

%% Exercise 1.

% Defining sigma and mu vectors
Sigma_aa = sigma_a^2 * [1, ro_1, ro_1^2; ro_1, 1, ro_1; ro_1^2, ro_1, 1];
mu = [mu_A; mu_B; mu_C];

%% Generating random values for A, B and C

p = (mvnrnd(mu, Sigma_aa, 1))';
A = p(1);
B = p(2);
C = p(3);

% Ex2 defines N
N = 5:5:20*5;

estimateData.p = p;
estimateData.N = N;
estimateData.mu = mu;
estimateData.sigma_w = sigma_w;
estimateData.ro_2 = ro_2;
estimateData.measPercent = 100;
estimateData.Sigma_aa = Sigma_aa;
estimateData.f_0 = f_0;
estimateData.t_0_s = t_0_s;

a_MS_full = MS_estimate(estimateData);

estimateData.measPercent = 10;
a_MS_10pc = MS_estimate(estimateData);










