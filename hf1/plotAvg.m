function plotAvg(avg, mu, N, k, t)
    
    subplot(1, 3, k);
    hold on
    plot(N, abs(avg(:,1) - mu(k)), '--.', 'markersize', 20);
    plot(N, abs(avg(:,2) - mu(k)), '--.', 'markersize', 20);
    plot(N, abs(avg(:,3) - mu(k)), '--.', 'markersize', 20);
    plot(N, abs(avg(:,4) - mu(k)), '--.', 'markersize', 20);
    hold off
    legend('ro1,ro2 != 0 & full cycle', 'ro1,ro2 = 0 & full cycle',...
           'ro1,ro2 != 0 & 10% cycle', 'ro1,ro2 = 0 & 10% cycle');
    title(t);
    xlabel('N')
    grid on
end