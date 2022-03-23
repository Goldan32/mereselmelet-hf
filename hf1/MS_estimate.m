function [a_MS] = MS_estimate(ed)
    M = size(ed.N, 2);
    
    a_MS.est = cell(M, 1);
    a_MS.var = cell(M, 1);
    a_MS.dis = cell(M, 1);
    
    for i=1:M
        [rowIndices, colIndices] = ndgrid(1:ed.N(i), 1:ed.N(i));
        C_matrix = ed.sigma_w^2 * (ed.ro_2 .^ (abs(rowIndices - colIndices)));

        % Generating noise as a function of time: w(t)
        wt = (mvnrnd(zeros(ed.N(i), 1), C_matrix, 1))';

        % Generating sinus signal for full cycle measurement
        stepsize = (ed.measPercent / 100) / (ed.N(i)*ed.f_0);
        t = ed.t_0_s: stepsize: ed.t_0_s + stepsize*(ed.N(i)-1);
        U = [(sin(2.*pi.*ed.f_0.*t))' (cos(2.*pi.*ed.f_0.*t))' ones(1, size(t, 2))'];
        z = U*ed.p + wt;

        a_MS.est{i, 1} = ed.mu + inv(U'*inv(C_matrix)*U + inv(ed.Sigma_aa))*U'*inv(C_matrix)*(z - U*ed.mu);
        a_MS.var{i, 1} = inv(U'*inv(C_matrix)*U + inv(ed.Sigma_aa));
        a_MS.dis{i, 1} = inv(U'*U)*U'*z - ed.mu;
    end
    
    
end