%% Initial Setup

% Clearing workspace
clear variables
close all

% Set verbose level (set to 1 for more figures)
VERBOSE = 0;

% Setting personal parameters as global variables
global sigma_a
loadPersonalParameters();
global mu

% Defining N vector
N = 5:5:20*5;

%% Ex1
[a_MS, a_LS, p] = measureAll(N);

%% Ex2

% Plot each figure on its own if VERBOSE is set (for deeper analysis)
if VERBOSE == 1
    plotEstimate(N, a_MS.full.diffro, p, 'MS Full cycle sample, ro nonzero');
    plotEstimate(N, a_MS.full.samero, p, 'MS Full cycle sample, ro zero');
    plotEstimate(N, a_MS.pc10.diffro, p, 'MS 10% cycle sample, ro nonzero');
    plotEstimate(N, a_MS.pc10.samero, p, 'MS 10% cycle sample, ro zero');
end

% Plot for comparing cases of ro=0 and 10% measurements
plotCompare(N, a_MS, 'est', p, 'MS estimation of ');
plotCompare(N, a_MS, 'var', p, 'MS spread of ');
plotCompare(N, a_MS, 'dis', p, 'MS distortion of ');

%% Ex3

% Plot each figure on its own if VERBOSE is set (for deeper analysis)
if VERBOSE == 1
    plotEstimate(N, a_LS.full.diffro, p, 'LS Full cycle sample, ro nonzero');
    plotEstimate(N, a_LS.full.samero, p, 'LS Full cycle sample, ro zero');
    plotEstimate(N, a_LS.pc10.diffro, p, 'LS 10% cycle sample, ro nonzero');
    plotEstimate(N, a_LS.pc10.samero, p, 'LS 10% cycle sample, ro zero');
end

% Plot for comparing cases of ro=0 and 10% measurements
plotCompare(N, a_LS, 'est', p, 'LS estimation of ');


%% Ex4 and Ex5 preparations
M = [5, 10, 100];

arr_A_MS = zeros(size(N, 2), 4, 100);
arr_B_MS = zeros(size(N, 2), 4, 100);
arr_C_MS = zeros(size(N, 2), 4, 100);

arr_A_LS = zeros(size(N, 2), 4, 100);
arr_B_LS = zeros(size(N, 2), 4, 100);
arr_C_LS = zeros(size(N, 2), 4, 100);
% Get results of 100 measurements in a 3D array
for k = 1:100
    [a_MS_tmp, a_LS_temp] = measureAll(N);

    tmp = extractForCompare(a_MS_tmp, 'est');
    arr_A_MS(:,:,k) = tmp.A;
    arr_B_MS(:,:,k) = tmp.B;
    arr_C_MS(:,:,k) = tmp.C;

    tmp = extractForCompare(a_MS_tmp, 'est');
    arr_A_LS(:,:,k) = tmp.A;
    arr_B_LS(:,:,k) = tmp.B;
    arr_C_LS(:,:,k) = tmp.C;

end

%%Ex4
% Calculate the avarage estimate based on 5,10 and 100 measurements (MS)
avg_A_5_MS   = mean(arr_A_MS(:,:,1:M(1)), 3);
avg_A_10_MS  = mean(arr_A_MS(:,:,1:M(2)), 3);
avg_A_100_MS = mean(arr_A_MS(:,:,1:M(3)), 3);

avg_B_5_MS   = mean(arr_B_MS(:,:,1:M(1)), 3);
avg_B_10_MS  = mean(arr_B_MS(:,:,1:M(2)), 3);
avg_B_100_MS = mean(arr_B_MS(:,:,1:M(3)), 3);

avg_C_5_MS   = mean(arr_C_MS(:,:,1:M(1)), 3);
avg_C_10_MS  = mean(arr_C_MS(:,:,1:M(2)), 3);
avg_C_100_MS = mean(arr_C_MS(:,:,1:M(3)), 3);

figure();
plotAvg(avg_A_5_MS, mu, N, 1, 'MS |Mean(A) - mu_A| based on 5 measurements');
plotAvg(avg_B_5_MS, mu, N, 2, 'MS |Mean(B) - mu_B| based on 5 measurements');
plotAvg(avg_C_5_MS, mu, N, 3, 'MS |Mean(C) - mu_C| based on 5 measurements');

figure();
plotAvg(avg_A_10_MS, mu, N, 1, 'MS |Mean(A) - mu_A| based on 10 measurements');
plotAvg(avg_B_10_MS, mu, N, 2, 'MS |Mean(B) - mu_B| based on 10 measurements');
plotAvg(avg_C_10_MS, mu, N, 3, 'MS |Mean(C) - mu_C| based on 10 measurements');

figure();
plotAvg(avg_A_100_MS, mu, N, 1, 'MS |Mean(A) - mu_A| based on 100 measurements');
plotAvg(avg_B_100_MS, mu, N, 2, 'MS |Mean(B) - mu_B| based on 100 measurements');
plotAvg(avg_C_100_MS, mu, N, 3, 'MS |Mean(C) - mu_C| based on 100 measurements');

% Calculate empirical spread for every case
spr_A_5_MS     = sqrt((1/M(1)).* sum((arr_A_MS(:,:,1:M(1)) - avg_A_5_MS).^2, 3));
spr_B_5_MS     = sqrt((1/M(2)).* sum((arr_A_MS(:,:,1:M(2)) - avg_A_5_MS).^2, 3));
spr_C_5_MS     = sqrt((1/M(3)).* sum((arr_A_MS(:,:,1:M(3)) - avg_A_5_MS).^2, 3));

spr_A_10_MS    = sqrt((1/M(1)).* sum((arr_A_MS(:,:,1:M(1)) - avg_A_10_MS).^2, 3));
spr_B_10_MS    = sqrt((1/M(2)).* sum((arr_A_MS(:,:,1:M(2)) - avg_A_10_MS).^2, 3));
spr_C_10_MS    = sqrt((1/M(3)).* sum((arr_A_MS(:,:,1:M(3)) - avg_A_10_MS).^2, 3));

spr_A_100_MS   = sqrt((1/M(1)).* sum((arr_A_MS(:,:,1:M(1)) - avg_A_100_MS).^2, 3));
spr_B_100_MS   = sqrt((1/M(2)).* sum((arr_A_MS(:,:,1:M(2)) - avg_A_100_MS).^2, 3));
spr_C_100_MS   = sqrt((1/M(3)).* sum((arr_A_MS(:,:,1:M(3)) - avg_A_100_MS).^2, 3));


figure();
plotSpread(spr_A_5_MS, sigma_a, N, 1, 'MS |EmpiricalSpread(A) - sigma_a| based on 5 measurements');
plotSpread(spr_B_5_MS, sigma_a, N, 2, 'MS |EmpiricalSpread(B) - sigma_a| based on 5 measurements');
plotSpread(spr_C_5_MS, sigma_a, N, 3, 'MS |EmpiricalSpread(C) - sigma_a| based on 5 measurements');

figure();
plotSpread(spr_A_10_MS, sigma_a, N, 1, 'MS |EmpiricalSpread(A) - sigma_a| based on 10 measurements');
plotSpread(spr_B_10_MS, sigma_a, N, 2, 'MS |EmpiricalSpread(B) - sigma_a| based on 10 measurements');
plotSpread(spr_C_10_MS, sigma_a, N, 3, 'MS |EmpiricalSpread(C) - sigma_a| based on 10 measurements');

figure();
plotSpread(spr_A_100_MS, sigma_a, N, 1, 'MS |EmpiricalSpread(A) - sigma_a| based on 100 measurements');
plotSpread(spr_B_100_MS, sigma_a, N, 2, 'MS |EmpiricalSpread(B) - sigma_a| based on 100 measurements');
plotSpread(spr_C_100_MS, sigma_a, N, 3, 'MS |EmpiricalSpread(C) - sigma_a| based on 100 measurements');


%% Ex5
% Calculate the avarage estimate based on 5,10 and 100 measurements (LS)
avg_A_5_LS   = mean(arr_A_LS(:,:,1:M(1)), 3);
avg_A_10_LS  = mean(arr_A_LS(:,:,1:M(2)), 3);
avg_A_100_LS = mean(arr_A_LS(:,:,1:M(3)), 3);

avg_B_5_LS   = mean(arr_B_LS(:,:,1:M(1)), 3);
avg_B_10_LS  = mean(arr_B_LS(:,:,1:M(2)), 3);
avg_B_100_LS = mean(arr_B_LS(:,:,1:M(3)), 3);

avg_C_5_LS   = mean(arr_C_LS(:,:,1:M(1)), 3);
avg_C_10_LS  = mean(arr_C_LS(:,:,1:M(2)), 3);
avg_C_100_LS = mean(arr_C_LS(:,:,1:M(3)), 3);

figure();
plotAvg(avg_A_5_LS, mu, N, 1, 'LS |Mean(A) - mu_A| based on 5 measurements');
plotAvg(avg_B_5_LS, mu, N, 2, 'LS |Mean(B) - mu_B| based on 5 measurements');
plotAvg(avg_C_5_LS, mu, N, 3, 'LS |Mean(C) - mu_C| based on 5 measurements');

figure();
plotAvg(avg_A_10_LS, mu, N, 1, 'LS |Mean(A) - mu_A| based on 10 measurements');
plotAvg(avg_B_10_LS, mu, N, 2, 'LS |Mean(B) - mu_B| based on 10 measurements');
plotAvg(avg_C_10_LS, mu, N, 3, 'LS |Mean(C) - mu_C| based on 10 measurements');

figure();
plotAvg(avg_A_100_LS, mu, N, 1, 'LS |Mean(A) - mu_A| based on 100 measurements');
plotAvg(avg_B_100_LS, mu, N, 2, 'LS |Mean(B) - mu_B| based on 100 measurements');
plotAvg(avg_C_100_LS, mu, N, 3, 'LS |Mean(C) - mu_C| based on 100 measurements');

% Calculate empirical spread for every case
spr_A_5_LS     = sqrt((1/M(1)).* sum((arr_A_LS(:,:,1:M(1)) - avg_A_5_LS).^2, 3));
spr_B_5_LS     = sqrt((1/M(2)).* sum((arr_A_LS(:,:,1:M(2)) - avg_A_5_LS).^2, 3));
spr_C_5_LS     = sqrt((1/M(3)).* sum((arr_A_LS(:,:,1:M(3)) - avg_A_5_LS).^2, 3));

spr_A_10_LS    = sqrt((1/M(1)).* sum((arr_A_LS(:,:,1:M(1)) - avg_A_10_LS).^2, 3));
spr_B_10_LS    = sqrt((1/M(2)).* sum((arr_A_LS(:,:,1:M(2)) - avg_A_10_LS).^2, 3));
spr_C_10_LS    = sqrt((1/M(3)).* sum((arr_A_LS(:,:,1:M(3)) - avg_A_10_LS).^2, 3));

spr_A_100_LS   = sqrt((1/M(1)).* sum((arr_A_LS(:,:,1:M(1)) - avg_A_100_LS).^2, 3));
spr_B_100_LS   = sqrt((1/M(2)).* sum((arr_A_LS(:,:,1:M(2)) - avg_A_100_LS).^2, 3));
spr_C_100_LS   = sqrt((1/M(3)).* sum((arr_A_LS(:,:,1:M(3)) - avg_A_100_LS).^2, 3));


figure();
plotSpread(spr_A_5_LS, sigma_a, N, 1, 'LS |EmpiricalSpread(A) - sigma_a| based on 5 measurements');
plotSpread(spr_B_5_LS, sigma_a, N, 2, 'LS |EmpiricalSpread(B) - sigma_a| based on 5 measurements');
plotSpread(spr_C_5_LS, sigma_a, N, 3, 'LS |EmpiricalSpread(C) - sigma_a| based on 5 measurements');

figure();
plotSpread(spr_A_10_LS, sigma_a, N, 1, 'LS |EmpiricalSpread(A) - sigma_a| based on 10 measurements');
plotSpread(spr_B_10_LS, sigma_a, N, 2, 'LS |EmpiricalSpread(B) - sigma_a| based on 10 measurements');
plotSpread(spr_C_10_LS, sigma_a, N, 3, 'LS |EmpiricalSpread(C) - sigma_a| based on 10 measurements');

figure();
plotSpread(spr_A_100_LS, sigma_a, N, 1, 'LS |EmpiricalSpread(A) - sigma_a| based on 100 measurements');
plotSpread(spr_B_100_LS, sigma_a, N, 2, 'LS |EmpiricalSpread(B) - sigma_a| based on 100 measurements');
plotSpread(spr_C_100_LS, sigma_a, N, 3, 'LS |EmpiricalSpread(C) - sigma_a| based on 100 measurements');

