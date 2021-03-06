function [a_MS] = MS_estimate(ed)
    M = size(ed.N, 2);
    
    for i=1:M

        % Selecting part of noise
        C_matrix_local = ed.C_matrix(1:ed.N(i), 1:ed.N(i));
        wt_local = ed.wt(1:ed.N(i));

        % Generating sinus signal for full cycle measurement
        stepsize = (ed.measPercent / 100) / (ed.N(i)*ed.f_0);
        t = ed.t_0_s: stepsize: ed.t_0_s + stepsize*(ed.N(i)-1);
        U = [(sin(2.*pi.*ed.f_0.*t))' (cos(2.*pi.*ed.f_0.*t))' ones(1, size(t, 2))'];

        z = U*ed.p + wt_local;

        % Estimate
        a_MS.est{i, 1} = ed.mu + inv(U'*inv(C_matrix_local)*U + inv(ed.Sigma_aa))*U'*inv(C_matrix_local)*(z - U*ed.mu);
        
        % Variance
        a_MS.var{i, 1} = diag(inv(U'*inv(C_matrix_local)*U + inv(ed.Sigma_aa)));
        
        % Distortion
        E = ed.mu + inv(U'*inv(C_matrix_local)*U + inv(ed.Sigma_aa)*U'*inv(C_matrix_local)*U*(ed.p - ed.mu));
        a_MS.dis{i, 1} = ed.p - E;
    end
    
    
end