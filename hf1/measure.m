function [a, p] = measure(M, estimateFunction, fields)

    global mu_A mu_B mu_C sigma_a sigma_w ro_1 ro_2 t_0_s f_0;

    % Ex2 defines N
    N = 5:5:20*5;
    a = initCellArray(N, fields);

    for i = 1:M
        
        % Defining sigma and mu vectors
        Sigma_aa = sigma_a^2 * [1, ro_1, ro_1^2; ro_1, 1, ro_1; ro_1^2, ro_1, 1];
        mu = [mu_A; mu_B; mu_C];
        
        p = (mvnrnd(mu, Sigma_aa, 1))';

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
        
        
        
        % full cycle estimation with ro not zero
        temp.full.diffro = estimateFunction(estimateData);
        
        % 10% of cycle with ro not zero
        estimateData.measPercent = 10;
        temp.pc10.diffro = estimateFunction(estimateData);
        
        % ro_1 = ro_2 = 0, 10% of cycle
        estimateData.Sigma_aa = sigma_a^2 * [1, 0, 0^2; 0, 1, 0; 0^2, 0, 1];
        estimateData.C_matrix = sigma_w^2 * (0 .^ (abs(rowIndices - colIndices)));
        temp.pc10.samero = estimateFunction(estimateData);
        
        % ro = 0, full cycle
        estimateData.measPercent = 100;
        temp.full.samero = estimateFunction(estimateData);
        
        a = incrementCellArrayElements(a, temp);

    end
    a = evaluateCellArray(M, a);
    
end