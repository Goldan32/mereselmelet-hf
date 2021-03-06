% This function is used for plotting each function
% in its own graph
function plotEstimate(N, a, p, title_appendum)

    figure();

    % Plot |A_estimate - A| of N
    if isfield(a, 'est')
        data = extractFromCellArray(a.est);

        subplot(3, 3, 1);
        plot(N, abs(data(:,1) - p(1)), ':.', 'markersize', 12);
        title(append('Estimation of A ~ ', title_appendum));
        ylabel('$$|\hat{A} - A|$$', 'interpreter', 'latex');
        xlabel('N');
        grid on;

        subplot(3, 3, 4);
        plot(N, abs(data(:,2) - p(2)), ':.', 'markersize', 12);
        title(append('Estimation of B ~ ', title_appendum));
        ylabel('$$|\hat{B} - B|$$', 'interpreter', 'latex');
        xlabel('N');
        grid on;

        subplot(3, 3, 7);
        plot(N, abs(data(:,3) - p(3)), ':.', 'markersize', 12);
        title(append('Estimation of C ~ ', title_appendum));
        ylabel('$$|\hat{C} - C|$$', 'interpreter', 'latex');
        xlabel('N');
        grid on;
    end
    
    % Plot distortion
    if isfield(a, 'dis')
        data = extractFromCellArray(a.dis);

        subplot(3, 3, 2);
        plot(N, data(:, 1), ':.', 'markersize', 12);
        title(append('Distortion of A ~ ', title_appendum));
        ylabel('$$b(\hat{A})$$', 'interpreter', 'latex');
        xlabel('N');
        grid on;

        subplot(3, 3, 5);
        plot(N, data(:, 2), ':.', 'markersize', 12);
        title(append('Distortion of B ~ ', title_appendum));
        ylabel('$$b(\hat{B})$$', 'interpreter', 'latex');
        xlabel('N');
        grid on;

        subplot(3, 3, 8);
        plot(N, data(:, 3), ':.', 'markersize', 12);
        title(append('Distortion of C ~ ', title_appendum));
        ylabel('$$b(\hat{C})$$', 'interpreter', 'latex');
        xlabel('N');
        grid on;
    end
    
    % Plot spread
    if isfield(a, 'var')
        data = extractFromCellArray(a.var);

        subplot(3, 3, 3);
        plot(N, sqrt(data(:, 1)), ':.', 'markersize', 12);
        title(append('Spread of A ~ ', title_appendum));
        ylabel('$$\sigma_{\hat{A}}$$', 'interpreter', 'latex');
        xlabel('N');
        grid on;

        subplot(3, 3, 6);
        plot(N, sqrt(data(:, 2)), ':.', 'markersize', 12);
        title(append('Spread of B ~ ', title_appendum));
        ylabel('$$\sigma_{\hat{B}}$$', 'interpreter', 'latex');
        xlabel('N');
        grid on;

        subplot(3, 3, 9);
        plot(N, sqrt(data(:, 3)), ':.', 'markersize', 12);
        title(append('Spread of C ~ ', title_appendum));
        ylabel('$$\sigma_{\hat{C}}$$', 'interpreter', 'latex');
        xlabel('N');
        grid on;
    end

end