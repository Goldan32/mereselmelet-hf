function [r] = extractForCompare(a, field)
    d1 = extractFromCellArray(a.full.diffro.(field));
    d2 = extractFromCellArray(a.full.samero.(field));
    d3 = extractFromCellArray(a.pc10.diffro.(field));
    d4 = extractFromCellArray(a.pc10.samero.(field));
    
    r.A = [d1(:,1), d2(:,1), d3(:,1), d4(:,1)];
    r.B = [d1(:,2), d2(:,2), d3(:,2), d4(:,2)];
    r.C = [d1(:,3), d2(:,3), d3(:,3), d4(:,3)];
end