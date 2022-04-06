% Plot spread results of Ex4 and Ex5.
% Calling this function creates a subplot.
% eg. plotting the 4 ways A parameter was measured.
function plotSpread(spr, sig, N, k, t)
    subplot(1, 3, k);
    hold on
    plot(N, abs(spr(:,1) - sig), '--.', 'markersize', 20);
    plot(N, abs(spr(:,2) - sig), '--.', 'markersize', 20);
    plot(N, abs(spr(:,3) - sig), '--.', 'markersize', 20);
    plot(N, abs(spr(:,4) - sig), '--.', 'markersize', 20);
    hold off
    legend('ro1,ro2 != 0 & full cycle', 'ro1,ro2 = 0 & full cycle',...
           'ro1,ro2 != 0 & 10% cycle', 'ro1,ro2 = 0 & 10% cycle');
    title(t);
    xlabel('N')
    grid on
end