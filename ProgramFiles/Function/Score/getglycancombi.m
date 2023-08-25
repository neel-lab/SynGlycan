function [glycancombi,glycancombimass] = getglycancombi(PTM,maxglycannum)
% NOTE: each cell in output "glycancombi" must contain the exact number of
%     glycans, e.g. for cell{2}, all glycan combinations has 2 glycans
PTMisglycanind = find(PTM.ptmisglycan);
PTMrealmass = PTM.ptmmass + PTM.ptmmassfix;
glycancombi = cell(1,maxglycannum);
glycancombimass = cell(1,maxglycannum);
for ii = 1:length(glycancombi)
    tempglycancombi = nmultichoosek([0,PTMisglycanind],ii);
    tempglycancombimass = zeros(size(tempglycancombi,1),1);
    for jj = 1:size(tempglycancombi,1)
        tempglycancombiind = tempglycancombi(jj,:);
        tempglycancombiind = tempglycancombiind(tempglycancombiind > 0);
        if ~isempty(tempglycancombiind)
            tempglycancombimass(jj) = sum(PTMrealmass(tempglycancombiind));
        end
    end
    glycancombi{ii} = tempglycancombi;
    glycancombimass{ii} = tempglycancombimass;
end
end