function [a_LS] = LS_estimate(ed)
    M = size(ed.N, 2);
    
    for i = 1:M

        % Selecting part of noise
        C_matrix_local = ed.C_matrix(1:ed.N(i), 1:ed.N(i));
        wt_local = ed.wt(1:ed.N(i));

        % Generating sinus signal measPercent% of a full cycle
        stepsize = (ed.measPercent / 100) / (ed.N(i)*ed.f_0);
        t = ed.t_0_s: stepsize: ed.t_0_s + stepsize*(ed.N(i)-1);
        U = [(sin(2.*pi.*ed.f_0.*t))' (cos(2.*pi.*ed.f_0.*t))' ones(1, size(t, 2))'];
        
        z = U*ed.p + wt_local;
        
        % Calculating LS estimation
        a_LS.est{i, 1} = inv(U'*inv(C_matrix_local)*U)*U'*inv(C_matrix_local)*z;

    end
end