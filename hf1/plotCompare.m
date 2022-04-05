% a has fields A, B and C. Each field is a matrix where rows are ascending
% by N and cols are diffrofull, samerofull, diffro10pc, samero10pc in order

function plotCompare(N, a, field, p, appendum)
    data = extractForCompare(a, field);
    ftable = ['A', 'B', 'C'];
    
    if strcmp(field, 'est')
        data.A = abs(data.A - p(1));
        data.B = abs(data.B - p(2));
        data.C = abs(data.C - p(3));
    elseif strcmp(field, 'var')
        data.A = sqrt(data.A);
        data.B = sqrt(data.B);
        data.C = sqrt(data.C);
    end
    
    
    figure();
    for i = 1:length(ftable)
        curData = data.(ftable(i));
        subplot(1, 3, i);
        hold on
        plot(N, curData(:,1), '--.', 'markersize', 20);
        plot(N, curData(:,2), '--.', 'markersize', 20);
        plot(N, curData(:,3), '--.', 'markersize', 20);
        plot(N, curData(:,4), '--.', 'markersize', 20);
        hold off
        title(append(appendum, ftable(i)));
        legend('ro1,ro2 != 0 & full cycle', 'ro1,ro2 = 0 & full cycle',...
               'ro1,ro2 != 0 & 10% cycle', 'ro1,ro2 = 0 & 10% cycle');
        grid on
    end
    
    
end