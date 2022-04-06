function [ca] = initCellArray(N, fields)

    for i = 1:size(fields, 1)
        for k = 1:size(N, 2)
            ca.full.diffro.(fields(i, :)){k, 1} = zeros(3, 1);
            ca.full.samero.(fields(i, :)){k, 1} = zeros(3, 1);
            ca.pc10.diffro.(fields(i, :)){k, 1} = zeros(3, 1);
            ca.pc10.samero.(fields(i, :)){k, 1} = zeros(3, 1);
        end
    end
end