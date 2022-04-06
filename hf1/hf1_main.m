%% Initial Setup
clear variables
close all
% clear all
global mu_A mu_B mu_C sigma_a sigma_w ro_1 ro_2 t_0_ms t_0_s f_0;
loadPersonalParameters();


%%
[a_MS, p] = measure(1, @MS_estimate, ['est'; 'var'; 'dis']);


% %% Personal parameters
% mu_A    = 0.5;
% mu_B    = 0.5;
% mu_C    = 0.5;
% sigma_a = 0.1;
% sigma_w = 0.1;
% ro_1    = 0.1;
% ro_2    = 0.1;
% t_0_ms  = 5;
% t_0_s   = t_0_ms * 1e-3;
% f_0     = 40;
% 
% 
% % Defining sigma and mu vectors
% Sigma_aa = sigma_a^2 * [1, ro_1, ro_1^2; ro_1, 1, ro_1; ro_1^2, ro_1, 1];
% mu = [mu_A; mu_B; mu_C];

% %% Generating random values for A, B and C
% p = (mvnrnd(mu, Sigma_aa, 1))';
% 
% % Ex2 defines N
% N = 5:5:20*5;

% % noise
% [rowIndices, colIndices] = ndgrid(1:N(end), 1:N(end));
% C_matrix = sigma_w^2 * (ro_2 .^ (abs(rowIndices - colIndices)));
% 
% % Generating noise as a function of time: w(t)
% wt = (mvnrnd(zeros(N(end), 1), C_matrix, 1))';
% 
% estimateData.p = p;
% estimateData.N = N;
% estimateData.mu = mu;
% estimateData.measPercent = 100;
% estimateData.Sigma_aa = Sigma_aa;
% estimateData.f_0 = f_0;
% estimateData.t_0_s = t_0_s;
% estimateData.wt = wt;
% estimateData.C_matrix = C_matrix;
% 
% % full cycle estimation with ro not zero
% a_MS.full.diffro = MS_estimate(estimateData);
% 
% % 10% of cycle with ro not zero
% estimateData.measPercent = 10;
% a_MS.pc10.diffro = MS_estimate(estimateData);
% 
% % ro_1 = ro_2 = 0, 10% of cycle
% estimateData.Sigma_aa = sigma_a^2 * [1, 0, 0^2; 0, 1, 0; 0^2, 0, 1];
% estimateData.C_matrix = sigma_w^2 * (0 .^ (abs(rowIndices - colIndices)));
% a_MS.pc10.samero = MS_estimate(estimateData);
% 
% % ro = 0, full cycle
% estimateData.measPercent = 100;
% a_MS.full.samero = MS_estimate(estimateData);
% 
%% Plot each figure on its own to make them easier to analyze
N = 5:5:20*5;
plotEstimate(N, a_MS.full.diffro, p, 'MS Full cycle sample, ro nonzero');
plotEstimate(N, a_MS.full.samero, p, 'MS Full cycle sample, ro zero');
plotEstimate(N, a_MS.pc10.diffro, p, 'MS 10% cycle sample, ro nonzero');
plotEstimate(N, a_MS.pc10.samero, p, 'MS 10% cycle sample, ro zero');

%% Ex2
% Plot for comparing cases of ro=0 and 10% measurements
plotCompare(N, a_MS, 'est', p, 'MS estimation of ');
plotCompare(N, a_MS, 'var', p, 'MS spread of ');
plotCompare(N, a_MS, 'dis', p, 'MS distortion of ');


%% LS estimation

estimateData.p = p;
estimateData.N = N;
estimateData.mu = mu;
estimateData.measPercent = 100;
estimateData.Sigma_aa = Sigma_aa;
estimateData.f_0 = f_0;
estimateData.t_0_s = t_0_s;
estimateData.wt = wt;
estimateData.C_matrix = C_matrix;

% full cycle estimation with ro not zero
a_LS.full.diffro = LS_estimate(estimateData);

% 10% of cycle, ro not zero
estimateData.measPercent = 10;
a_LS.pc10.diffro = LS_estimate(estimateData);

% ro_1 = ro_2 = 0, 10% cycle
estimateData.Sigma_aa = sigma_a^2 * [1, 0, 0^2; 0, 1, 0; 0^2, 0, 1];
estimateData.C_matrix = sigma_w^2 * (0 .^ (abs(rowIndices - colIndices)));
a_LS.pc10.samero = LS_estimate(estimateData);

% full cycle, ro = 0
estimateData.measPercent = 100;
a_LS.full.samero = LS_estimate(estimateData);

%% Plot each figure on its own to make them easier to analyze
plotEstimate(N, a_LS.full.diffro, p, 'LS Full cycle sample, ro nonzero');
plotEstimate(N, a_LS.full.samero, p, 'LS Full cycle sample, ro zero');
plotEstimate(N, a_LS.pc10.diffro, p, 'LS 10% cycle sample, ro nonzero');
plotEstimate(N, a_LS.pc10.samero, p, 'LS 10% cycle sample, ro zero');

%% Plot for comparing cases of ro=0 and 10% measurements
plotCompare(N, a_LS, 'est', p, 'LS estimation of ');
