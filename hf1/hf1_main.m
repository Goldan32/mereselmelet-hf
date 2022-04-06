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


%% Ex4 and Ex5
M = [5, 10, 100];

arr_A = zeros(size(N, 2), 4, 100);
arr_B = zeros(size(N, 2), 4, 100);
arr_C = zeros(size(N, 2), 4, 100);
% Get results of 100 measurements in a 3D array
for k = 1:100
    a_MS_tmp = measureAll(N);
    tmp = extractForCompare(a_MS_tmp, 'est');
    arr_A(:,:,k) = tmp.A;
    arr_B(:,:,k) = tmp.B;
    arr_C(:,:,k) = tmp.C;
end

% Calculate the avarage estimate based on 5,10 and 100 measurements
avg_A_5   = mean(arr_A(:,:,1:M(1)), 3);
avg_A_10  = mean(arr_A(:,:,1:M(2)), 3);
avg_A_100 = mean(arr_A(:,:,1:M(3)), 3);

avg_B_5   = mean(arr_B(:,:,1:M(1)), 3);
avg_B_10  = mean(arr_B(:,:,1:M(2)), 3);
avg_B_100 = mean(arr_B(:,:,1:M(3)), 3);

avg_C_5   = mean(arr_C(:,:,1:M(1)), 3);
avg_C_10  = mean(arr_C(:,:,1:M(2)), 3);
avg_C_100 = mean(arr_C(:,:,1:M(3)), 3);

figure();
plotAvg(avg_A_5, mu, N, 1, '|Mean(A) - mu_A| based on 5 measurements');
plotAvg(avg_B_5, mu, N, 2, '|Mean(B) - mu_B| based on 5 measurements');
plotAvg(avg_C_5, mu, N, 3, '|Mean(C) - mu_C| based on 5 measurements');

figure();
plotAvg(avg_A_10, mu, N, 1, '|Mean(A) - mu_A| based on 10 measurements');
plotAvg(avg_B_10, mu, N, 2, '|Mean(B) - mu_B| based on 10 measurements');
plotAvg(avg_C_10, mu, N, 3, '|Mean(C) - mu_C| based on 10 measurements');

figure();
plotAvg(avg_A_100, mu, N, 1, '|Mean(A) - mu_A| based on 100 measurements');
plotAvg(avg_B_100, mu, N, 2, '|Mean(B) - mu_B| based on 100 measurements');
plotAvg(avg_C_100, mu, N, 3, '|Mean(C) - mu_C| based on 100 measurements');

% Calculate empirical spread for every case
spr_A_5     = sqrt((1/M(1)).* sum((arr_A(:,:,1:M(1)) - avg_A_5).^2, 3));
spr_B_5     = sqrt((1/M(2)).* sum((arr_A(:,:,1:M(2)) - avg_A_5).^2, 3));
spr_C_5     = sqrt((1/M(3)).* sum((arr_A(:,:,1:M(3)) - avg_A_5).^2, 3));

spr_A_10    = sqrt((1/M(1)).* sum((arr_A(:,:,1:M(1)) - avg_A_10).^2, 3));
spr_B_10    = sqrt((1/M(2)).* sum((arr_A(:,:,1:M(2)) - avg_A_10).^2, 3));
spr_C_10    = sqrt((1/M(3)).* sum((arr_A(:,:,1:M(3)) - avg_A_10).^2, 3));

spr_A_100   = sqrt((1/M(1)).* sum((arr_A(:,:,1:M(1)) - avg_A_100).^2, 3));
spr_B_100   = sqrt((1/M(2)).* sum((arr_A(:,:,1:M(2)) - avg_A_100).^2, 3));
spr_C_100   = sqrt((1/M(3)).* sum((arr_A(:,:,1:M(3)) - avg_A_100).^2, 3));


figure();
plotSpread(spr_A_5, sigma_a, N, 1, '|EmpiricalSpread(A) - sigma_a| based on 5 measurements');
plotSpread(spr_B_5, sigma_a, N, 2, '|EmpiricalSpread(B) - sigma_a| based on 5 measurements');
plotSpread(spr_C_5, sigma_a, N, 3, '|EmpiricalSpread(C) - sigma_a| based on 5 measurements');

figure();
plotSpread(spr_A_10, sigma_a, N, 1, '|EmpiricalSpread(A) - sigma_a| based on 10 measurements');
plotSpread(spr_B_10, sigma_a, N, 2, '|EmpiricalSpread(B) - sigma_a| based on 10 measurements');
plotSpread(spr_C_10, sigma_a, N, 3, '|EmpiricalSpread(C) - sigma_a| based on 10 measurements');

figure();
plotSpread(spr_A_100, sigma_a, N, 1, '|EmpiricalSpread(A) - sigma_a| based on 100 measurements');
plotSpread(spr_B_100, sigma_a, N, 2, '|EmpiricalSpread(B) - sigma_a| based on 100 measurements');
plotSpread(spr_C_100, sigma_a, N, 3, '|EmpiricalSpread(C) - sigma_a| based on 100 measurements');


