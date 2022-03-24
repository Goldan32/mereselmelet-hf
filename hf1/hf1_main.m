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

% Ex2 defines N
N = 5:5:20*5;

% noise
[rowIndices, colIndices] = ndgrid(1:N(end), 1:N(end));
C_matrix = sigma_w^2 * (ro_2 .^ (abs(rowIndices - colIndices)));

% Generating noise as a function of time: w(t)
wt = (mvnrnd(zeros(N(end), 1), C_matrix, 1))';

estimateData.p = p;
estimateData.N = N;
estimateData.mu = mu;
estimateData.measPercent = 100;
estimateData.Sigma_aa = Sigma_aa;
estimateData.f_0 = f_0;
estimateData.t_0_s = t_0_s;
estimateData.wt = wt;
estimateData.C_matrix = C_matrix;

a_MS.full.diffro = MS_estimate(estimateData);

estimateData.measPercent = 10;
a_MS.pc10.diffro = MS_estimate(estimateData);

% ro_1 = ro_2 = 0
estimateData.Sigma_aa = sigma_a^2 * [1, 0, 0^2; 0, 1, 0; 0^2, 0, 1];
a_MS.pc10.samero = MS_estimate(estimateData);

estimateData.measPercent = 100;
a_MS.full.samero = MS_estimate(estimateData);

plotEstimate(N, a_MS.full.diffro, p, 'Full cycle sample, ro nonzero');
plotEstimate(N, a_MS.full.samero, p, 'Full cycle sample, ro zero');
plotEstimate(N, a_MS.pc10.diffro, p, '10% cycle sample, ro nonzero');
plotEstimate(N, a_MS.pc10.samero, p, '10% cycle sample, ro zero');




