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

% Defining C (noise) matrix
M  = size(N,2);
a_est = cell(M, 1);
a_var = cell(M, 1);
a_dis = cell(M, 1);
for i=1:M
    [rowIndices, colIndices] = ndgrid(1:N(i), 1:N(i));
    C_matrix = sigma_w^2 * (ro_2 .^ (abs(rowIndices - colIndices)));

    % Generating noise as a function of time: w(t)
    wt = (mvnrnd(zeros(N(i), 1), C_matrix, 1))';

    % Generating sinus signals
    stepsize = 1/ (N(i)*f_0);
    t = t_0_s: stepsize: t_0_s + stepsize*(N(i)-1);

    U = [(sin(2.*pi.*f_0.*t))' (cos(2.*pi.*f_0.*t))' ones(1, size(t, 2))'];
    
    z = U*p + wt;
    
    a_est_MS = mu + inv(U'*inv(C_matrix)*U + inv(Sigma_aa))*U'*inv(C_matrix)*(z - U*mu);
    var_a_est_MS = inv(U'*inv(C_matrix)*U + inv(Sigma_aa));
    distortion = inv(U'*U)*U'*z - mu;
    
    
end

%% Calculate

