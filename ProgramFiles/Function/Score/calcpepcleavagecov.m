function [pf,gf,pepcov] = calcpepcleavagecov(theores,theofrag,fragAAind,peplen,glypos)
% CALCPEPCLEAVAGECOV: calculate percentage of peptide bond cleavages
%
% Syntax:
% [pf,gf,pepcov] = calcpepcleavagecov(theores,theofrag,fragAAind,peplen,glypos)
%
% Input: 
% theores: structure, the matching result. Must contain field
%     "ionmatchindex".
% theofrag: structure, theoretical fragments.
% fragAAind: 1 x n cell array, for each theoretical fragment, which amino
%     acid is present.
% peplen: length of peptide.
% glypos: which amino acids are glycans attached to.
% 
% Output:
% (Below "percentage" means [#bond broken]/[#peptide bonds], #peptide
% bonds = #amino acids - 1)
% pf: 1 x n numerical array, for all cleavages that produce peptide frag
% without glycan, is the bond fragmented or not. 1 means fragmented. 0 means
% not fragmented
% gf: similar to pf, for all cleavages that produce peptide frag with glycan.
% pepcov: regardless of glycosylation site, the percentage of peptide
% fragments that are matched.
%
% Note:
% N/A
%
% Example:
% N/A
%
% Children function: 
% N/A
% 

purepep = false(size(theofrag));
for ii = 1:length(theofrag)
    if ~any(ismember(theofrag(ii).unitindex{1},glypos))
        purepep(ii) = true;
    end
end
purepep = purepep & [theofrag.npFrag] > 0;
hasgly = ~purepep;
if iscell(theores.ionmatchindex)
    matchedpepfragind = find(~cellfun(@isempty,theores.ionmatchindex) & logical(purepep));
    matchedglyfragind = find(~cellfun(@isempty,theores.ionmatchindex) & logical(hasgly));
else
    matchedpepfragind = find(theores.ionmatchindex > 0 & logical(purepep));
    matchedglyfragind = find(theores.ionmatchindex > 0 & logical(hasgly));
end
pf = false(peplen - 1,1);
gf = false(peplen - 1,1);
for ii = 1:length(matchedpepfragind)
    aaind = fragAAind{matchedpepfragind(ii)};
    writeind = [min(aaind),max(aaind)] + [-1 0];  
    % when cleavage occurs at mth bond, it's the bond between "m-1" and
    % "m".this event should be written at the (m-1)th place.
    writeind = writeind(writeind > 0 & writeind < peplen);
    if ~any(ismember(aaind,glypos))
        pf(writeind) = 1;
    end
end
for ii = 1:length(matchedglyfragind)
    aaind = fragAAind{matchedglyfragind(ii)};
    if ~isempty(aaind)
        writeind = [min(aaind),max(aaind)] + [-1 0];
        writeind = writeind(writeind > 0 & writeind < peplen);
        if any(ismember(aaind,glypos))
            gf(writeind) = 1;
        end
    end
end
pepcov = sum(pf|gf) / (peplen - 1);
end