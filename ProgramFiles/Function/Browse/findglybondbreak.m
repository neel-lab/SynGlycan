function [glybondbreakpos,matchedfrag] = findglybondbreak(PGM,theofrag,theores)
% FINDGLYBONDBREAK: find the fragmented glycosidic bonds associated with
%     each fragment ion
%
% Syntax:
% [glybondbreakpos,matchedfrag] = findglybondbreak(PGM,theofrag,theores)
%
% Input:
% PGM: [PGM{1},PGM{2},PGM{3}] = breakGlyPep(SGP)
% theofrag: Theoretical fragments of candidate (glyco)peptide.
% theores: result returned by CALCITHSCORE, must contain field
%     "peakmatchindex".
%
% Output:
% glybondbreakpos: identified glycosidic bond cleavage position, for
%     multiple glycans, the position is counted from the very first
%     monosaccharide in the glycopeptide(starting from N-terminus),
%     so these numbers could be big.
% matchedfrag: part of "theofrag" that are matched.
%
% Note:
% this is a modified version of CALCGLYCLEAVAGECOV
%
% Example:
% N/A
%
% Children function:
% N/A
%
% See also:
% CALCGLYCLEAVAGECOV
%

glybondbreakpos  = {};
if iscell(theores.ionmatchindex)
    matchedfrag = theofrag(~cellfun(@isempty,theores.ionmatchindex));
else
    matchedfrag = theofrag(theores.ionmatchindex);
end

if ~isempty(PGM{2})
    glycanlengths = [PGM{2}.len];
    glybmends = glycanlengths * triu(ones(length(glycanlengths)));
    glybmstarts = glybmends - glycanlengths + 1;
    glycanbondmaps = zeros(max(glybmends));
    for i = 1:length(PGM{2})
        glycanbondmaps(glybmstarts(i):glybmends(i),glybmstarts(i):glybmends(i)) = ...
            getglycanbondmap(PGM{2}(i).struct);
    end
    fragglyind = cell(size(theofrag));
    for i = 1:length(theofrag)
        fragglyind{i} = theofrag(i).unitindex{2}(:,1);
    end
    matchedfragglyind = fragglyind(~cellfun(@isempty,theores.ionmatchindex));
    msconnecttedby = cell(1,size(glycanbondmaps,1));
    msconnectto = cell(1,size(glycanbondmaps,1));
    glybondbreakpos = cell(length(matchedfragglyind),2);
    for i = 1:size(glycanbondmaps,1)
        tempmsconnecttedby = find(glycanbondmaps(i,:));
        tempmsconnectto = find(glycanbondmaps(:,i));
        if ~isempty(tempmsconnecttedby)
            msconnecttedby{i} = tempmsconnecttedby;
        end
        if ~isempty(tempmsconnectto)
            msconnectto{i} = tempmsconnectto;
        end
    end
    for i = 1:length(matchedfragglyind)
        thisfragglyind = matchedfragglyind{i};
        if ~isempty(thisfragglyind)
            NRside = [msconnecttedby{thisfragglyind}];
            Rside = [msconnectto{thisfragglyind}];
            glybondbreakpos{i,1} = setdiff(NRside,thisfragglyind);
            glybondbreakpos{i,2} = setdiff(Rside,thisfragglyind);
        end
    end
end
end