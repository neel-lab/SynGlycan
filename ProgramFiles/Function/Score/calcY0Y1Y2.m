function Y0Y1Y2 = calcY0Y1Y2(spectrum,theores,fragglyind,fragAAind,peplen)
% CALCY0Y1Y2: calculate if Y0/1/2 ions are matched, if so, what's the
% ratio of intensity compared to highest peak.
%
% Syntax:
% Y0Y1Y2 = calcY0Y1Y2(spectrum,theores,fragglyind,fragAAind,peplen)
%
% Input:
% spectrum: experimental spectrum, sorted by descending intensity.
% theores: structure, the matching result. Must contain field "ionmatchindex"
% fragglyind: 1 x n cell array, for each theoretical fragment, which
%     monosaccharide is present.
% fragAAind: 1 x n cell array, for each theoretical fragment, which amino
%     acid is present.
% peplen: length of peptide.
% 
% Output:
% Y0Y1Y2: [Y0 intensity, Y1 intensity, Y2 intensity]
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

fragpeplen = cellfun(@length,fragAAind);
fragmsnum = cellfun(@length,fragglyind);
posY0 = fragpeplen == peplen & fragmsnum == 0;
posY1 = fragpeplen == peplen & fragmsnum == 1;
posY2 = fragpeplen == peplen & fragmsnum == 2;
peakmatchY0 = [];
peakmatchY1 = [];
peakmatchY2 = [];
if any(posY0)
    peakmatchY0 = theores.ionmatchindex{posY0};
end
if any(posY1)
    peakmatchY1 = theores.ionmatchindex{posY1};
end
if any(posY2)
    peakmatchY2 = theores.ionmatchindex{posY2};
end
Y0intfrac = 0;
Y1intfrac = 0;
Y2intfrac = 0;
if ~isempty(peakmatchY0)
    Y0intfrac = max(spectrum(peakmatchY0,2))/spectrum(1,2);
end
if ~isempty(peakmatchY1)
    Y1intfrac = max(spectrum(peakmatchY1,2))/spectrum(1,2);
end
if ~isempty(peakmatchY2)
    Y2intfrac = max(spectrum(peakmatchY2,2))/spectrum(1,2);
end
Y0Y1Y2 = [Y0intfrac,Y1intfrac,Y2intfrac];
end