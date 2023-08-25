function glycov = calcglycleavagecov(PGM,fragglyind,theores)
% CALCGLYCLEAVAGECOV: calculate ratio of glycosidic bond cleavages
%
% Syntax:
% glycov = calcglycleavagecov(PGM,fragglyind,theores)
%
% Input:
% PGM: [PGM{1},PGM{2},PGM{3}] = breakGlyPep(SGP)
% fragglyind: 1 x n cell array, for each theoretical fragment, which
%     monosaccharide is present.
% theores: structure, the matching result. Must contain field "ionmatchindex"
%
% Output:
% glycov: ratio of glycosidic bond cleavages. When there are multiple
%     glycans, results for each glycan will be added up to return 1 value.
%
% Note:
% N/A
%
% Example:
% (A. every glycosidic bond was fragmented)
% testglycan = '{n{n{h{h{n{h}}}{h{n{h}}}}}}';
% [PGM{1},PGM{2},PGM{3}] = breakGlyPep(testglycan);
% theofrag = multiSGPFrag(testglycan,0,1,0,1,'allion');
% fragunitind = reshape([theofrag.unitindex],3,[]);
% fragglyind = cellfun(@(x) x(:,1),fragunitind(2,:),'UniformOutPut',false);
% theores.ionmatchindex = true(size(theofrag));
% glycov = calcglycleavagecov(PGM,fragglyind,theores)
% Answer:
% glycov =
% 
%          1
% 
% (B. continued from last example, now no fragmentation at all)
% theores.ionmatchindex = false(size(theofrag));
% glycov = calcglycleavagecov(PGM,fragglyind,theores)
% 
% glycov =
% 
%      0
% 
% Children function: 
% GETGLYCANBONDMAP
%

glycanlengths = [PGM{2}.len];
glybmends = glycanlengths * triu(ones(length(glycanlengths)));  % get number of monosacharides for each glycan
glybmstarts = glybmends - glycanlengths + 1;
glycanbondmaps = zeros(max(glybmends));
for i = 1:length(PGM{2})
    glycanbondmaps(glybmstarts(i):glybmends(i),glybmstarts(i):glybmends(i)) = ...
        getglycanbondmap(PGM{2}(i).struct);  
    % stack bondmap, it's easier when you want to consider every glycan in gp
end
if iscell(theores.ionmatchindex)
    matchedfragglyind = fragglyind(~cellfun(@isempty,theores.ionmatchindex));
else
    matchedfragglyind = fragglyind(theores.ionmatchindex);
end
% theores has a good property: monosac in each frag has unique index
% number, the following "unique" operation won't delete useful info
unitinds = cellfun(@(x) num2str(x'),matchedfragglyind,'UniformOutput',false);  % transform to string before "unique"
[~,keepind]=unique(unitinds);
matchedfragglyind = matchedfragglyind(keepind);
msconnectedby = cell(1,size(glycanbondmaps,1));
msconnectto = cell(1,size(glycanbondmaps,1));
% connected by/connect to: for each monosac, find out upstream/downstream
% connections
for i = 1:size(glycanbondmaps,1)
    tempmsconnecttedby = find(glycanbondmaps(i,:));
    tempmsconnectto = find(glycanbondmaps(:,i));
    if ~isempty(tempmsconnecttedby)
        msconnectedby{i} = tempmsconnecttedby;
    end
    if ~isempty(tempmsconnectto)
        msconnectto{i} = tempmsconnectto;
    end
end
% CALCULATION START
for i = 1:length(matchedfragglyind)
    thisfragglyind = matchedfragglyind{i};
    if ~isempty(thisfragglyind)
        terminalside = msconnectedby{thisfragglyind};  % find out which monosac this one connects to
        rootside = msconnectto{thisfragglyind};
        del_col = setdiff(terminalside,thisfragglyind);
        for j = 1:length(thisfragglyind)  % cleave bond: frag - terminal side
            for k = 1:length(del_col)
                glycanbondmaps(thisfragglyind(j),del_col(k)) = 0;
                % monosac appeared in fragment, this means all bonds
                % connecting this one is cleaved. set bondmap to 0 to show
                % cleavage
            end
        end
        del_row = setdiff(rootside,thisfragglyind);
        for j = 1:length(thisfragglyind)  % cleave bond: frag - root side
            for k = 1:length(del_row)
                glycanbondmaps(del_row(k),thisfragglyind(j)) = 0;
            end
        end
    end
end
glycov = 1 - (sum(sum(glycanbondmaps)))/(sum(glycanlengths)-length(glycanlengths));
% about sum(x) - length(x): because number of bonds = number of monosac - 1
% (peptide glycosylation bond),
end