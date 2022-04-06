function [y] = evaluateCellArray(M, ca)
    ftable = ['est'; 'var'; 'dis'];

    for i = 1:size(ftable, 1)
        if isfield(ca.full.diffro, ftable(i, :))
            for k = 1:size(ca.full.diffro.(ftable(i, :)), 1)
                y.full.diffro.(ftable(i, :)){k, 1} = ca.full.diffro.(ftable(i, :)){k, 1} ./ M;
                y.full.samero.(ftable(i, :)){k, 1} = ca.full.samero.(ftable(i, :)){k, 1} ./ M;
                y.pc10.diffro.(ftable(i, :)){k, 1} = ca.pc10.diffro.(ftable(i, :)){k, 1} ./ M;
                y.pc10.samero.(ftable(i, :)){k, 1} = ca.pc10.samero.(ftable(i, :)){k, 1} ./ M;
            end
        end
     end
end