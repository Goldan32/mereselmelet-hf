function [matrix] = extractFromCellArray(ca)
    rows = size(ca, 1);
    cols = size(ca{1}, 1);

    matrix = zeros(rows, cols);

    for i = 1:rows
        for j = 1:cols
            matrix(i,j) = ca{i}(j);
        end
    end

end