function plotEstimate(N, a, p, title_appendum)

    figure();

    % Plot A_estimate - A of N
    data = extractFromCellArray(a.est);
    
    subplot(3, 3, 1);
    plot(N, abs(data(:,1) - p(1)));
    title(append('Estimation of A ~ ', title_appendum));
    ylabel('$$|\hat{A} - A|$$', 'interpreter', 'latex');
    xlabel('N');
    grid on;
    
    subplot(3, 3, 4);
    plot(N, abs(data(:,2) - p(2)));
    title(append('Estimation of B ~ ', title_appendum));
    ylabel('$$|\hat{B} - B|$$', 'interpreter', 'latex');
    xlabel('N');
    grid on;

    subplot(3, 3, 7);
    plot(N, abs(data(:,3) - p(3)));
    title(append('Estimation of C ~ ', title_appendum));
    ylabel('$$|\hat{C} - C|$$', 'interpreter', 'latex');
    xlabel('N');
    grid on;
    
    
    % Plot distortion
    data = extractFromCellArray(a.dis);
    
    subplot(3, 3, 2);
    plot(N, data(:, 1));
    title(append('Distortion of A ~ ', title_appendum));
    ylabel('$$b(\hat{A})$$', 'interpreter', 'latex');
    xlabel('N');
    grid on;
    
    subplot(3, 3, 5);
    plot(N, data(:, 2));
    title(append('Distortion of B ~ ', title_appendum));
    ylabel('$$b(\hat{B})$$', 'interpreter', 'latex');
    xlabel('N');
    grid on;

    subplot(3, 3, 8);
    plot(N, data(:, 3));
    title(append('Distortion of C ~ ', title_appendum));
    ylabel('$$b(\hat{C})$$', 'interpreter', 'latex');
    xlabel('N');
    grid on;
    
    % Plot variance
    data = extractFromCellArray(a.var);

    
    subplot(3, 3, 3);
    plot(N, data(:, 1));
    title(append('Variance of A ~ ', title_appendum));
    ylabel('$$var(\hat{A})$$', 'interpreter', 'latex');
    xlabel('N');
    grid on;
    
    subplot(3, 3, 6);
    plot(N, data(:, 2));
    title(append('Variance of B ~ ', title_appendum));
    ylabel('$$var(\hat{B})$$', 'interpreter', 'latex');
    xlabel('N');
    grid on;

    subplot(3, 3, 9);
    plot(N, data(:, 3));
    title(append('Variance of C ~ ', title_appendum));
    ylabel('$$var(\hat{C})$$', 'interpreter', 'latex');
    xlabel('N');
    grid on;
    
    
end