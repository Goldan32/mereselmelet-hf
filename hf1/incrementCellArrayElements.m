function [y] = incrementCellArrayElements(ca, data)

    ftable = ['est'; 'var'; 'dis'];
    
    ca.full.diffro

    for i = 1:size(ftable)
        ftable(i, :)
        if isfield(ca.full.diffro, ftable(i, :))
            disp('entered here');
            for k = 1:size(ca.full.diffro.(ftable(i, :)), 1)
                y.full.diffro.(ftable(i, :)){k, 1} = ca.full.diffro.(ftable(i, :)){k, 1} + data.full.diffro.(ftable(i, :)){k, 1};
                y.full.samero.(ftable(i, :)){k, 1} = ca.full.samero.(ftable(i, :)){k, 1} + data.full.samero.(ftable(i, :)){k, 1};
                y.pc10.diffro.(ftable(i, :)){k, 1} = ca.pc10.diffro.(ftable(i, :)){k, 1} + data.pc10.diffro.(ftable(i, :)){k, 1};
                y.pc10.samero.(ftable(i, :)){k, 1} = ca.pc10.samero.(ftable(i, :)){k, 1} + data.pc10.samero.(ftable(i, :)){k, 1};
            end
        else
            disp('not good');
        end
    end
end